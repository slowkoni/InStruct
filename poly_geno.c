/*
 *  poly_geno.c
 *  
 *
 *  Created by Hong Gao on 3/7/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "nrutil.h"
#include "random.h"
#include "data_interface.h"
#include "initial.h"
#include "check_converg.h"
#include "mcmc.h"
#include "result_analysis.h"
#include "poly_geno.h"

static void gen_allele_poly(SEQDATA data, POLY *polyp);
static void gen_polyinfo(POLY *polyploid, SEQDATA data);
static void free_poly(POLY *polyp, SEQDATA data);
static void initial_geno(UPMCMC *ptr, SEQDATA data,POLY **polyploid);
static void initial_chn(double ***qqnum,SEQDATA data,UPMCMC **ptr,CHAIN *chain,INIT initial,int chn);
static void update_P_auto(UPMCMC **ptr,SEQDATA data);
static void update_P_allo(UPMCMC **ptr,SEQDATA data);
static void update_geno(UPMCMC *ptr, SEQDATA data, POLY *polypld);
static void update_S_POP(SEQDATA data,UPMCMC **ptr,POLY *polypld);
static double cal_lkd_props(int id,float **genofreq_tmp,UPMCMC *ptr,SEQDATA data,POLY *polypld);
static double cal_lkd(UPMCMC *ptr,SEQDATA data,POLY *polypld);
static void move_genofreq(int id, float **genofreq_tmp,POLY *polypld,SEQDATA data);
static void update_ZQ(UPMCMC **ptr,SEQDATA data,int init_flag,double ***qqnum,POLY *polypld);
static int choose_unif(int temp);
static int choose_two_auto(int indv_id,int loci_id,int temp,UPMCMC *ptr,SEQDATA data, POLY *polypld);
static int choose_tri_auto(int indv_id,int loci_id,int temp,UPMCMC *ptr,SEQDATA data, POLY *polypld);
static int choose_two_allo(int indv_id,int loci_id,int temp,UPMCMC *ptr,SEQDATA data, POLY *polypld);
static int choose_tri_allo(int indv_id,int loci_id,int temp,UPMCMC *ptr,SEQDATA data, POLY *polypld);
static int choose_tetra_allo(int indv_id,int loci_id,int temp,UPMCMC *ptr,SEQDATA data, POLY *polypld);
static void calc_exfreq_auto(UPMCMC *ptr,SEQDATA data,POLY *polypld);
static void calc_exfreq_allo(UPMCMC *ptr,SEQDATA data,POLY *polypld);
static void calc_self_genofreq(double self_rate,float **genofreq,SEQDATA data,POLY *polypld,int k);
static double calc_genofq(int loci_id,int indv_id,int *z, SEQDATA data, UPMCMC *ptr,POLY *polypld);
static int get_index_auto(int locus_id, int ploid,int *geno,POLY *polypld,int allelenum,int *cat_id);
static int get_index_allo(int locus_id, int ploid,int *geno,POLY *polypld,int allelenum,int *cat_id);
static int get_cat_auto(int ploid,int *geno);
static int get_cat_allo(int ploid,int *geno);
static int check_rule_auto(int ploid, int *geno,int cat_id);
static int check_rule_allo(int ploid, int *geno,int cat_id);
static void change_geno_auto(int ploid, int **geno,int cat_id);
static void change_geno_allo(int ploid, int **geno,int cat_id);
static int copy_num(int val, int *vec, int leng);
static void auto_geno_num(POLY *polyp);
static void auto_geno_list(int ploid,POLY *);
static void auto_genfreq(float self, int pop_id,int loci_id,SEQDATA data, POLY *polyp, float *freq);
static void allo_geno_num(POLY *polyp);
static void allo_geno_list(int ploid,POLY *);
static void allo_genfreq(float self, int pop_id,int loci_id,SEQDATA data, POLY *polyp, float *freq);
static int calc_val(int *, int num,int i);
static int calc_val2(int* ,int ,int , int,int );
static int find_id(int num, int *array,int);
static void gaussj(float **a, int n, float **b, int m);
static void two_allele_auto(int num,  int i,int locus_id, UPMCMC *ptr,SEQDATA data);
static void two_allele_allo(int num,  int i,int locus_id, UPMCMC *ptr,SEQDATA data);
static void tri_allele_auto(int num,  int i,int locus_id, UPMCMC *ptr,SEQDATA data);
static void tri_allele_allo(int num,  int i,int locus_id, UPMCMC *ptr,SEQDATA data);
static void tetra_allele_allo(int num,  int i,int locus_id, UPMCMC *ptr,SEQDATA data);


CHAIN mcmc_POP_tetra_selfing(SEQDATA data, INIT initial,int chn, CONVG *cvg)
{
	int i;
	long cnt_step=0,step;
	double **qqnum;	
	CHAIN chain;	
	UPMCMC *ptr=NULL;	
	POLY *tetra;
	
	if((tetra=(POLY *)malloc(sizeof(POLY)))==NULL)
	{	nrerror("Allocation failure in tetra");}
	gen_polyinfo(tetra,data);
	initial_chn(&qqnum,data,&ptr,&chain,initial,chn);
	initial_geno(ptr,data,&tetra);	
	
	for(i=0;i<data.popnum;i++)			//initialize selfing rates
	{
		ptr->self_rates[i]=initial.initd[chn][i];
		if(data.back_refl==0) 
		{	ptr->state[i]=dt_stat(ptr->self_rates[i]);}	//determine whether ptr->self_rates[i] is 0.0 (0)state,(0,1) (1)state or 1.0 (2)state
	}
	update_ZQ(&ptr,data,1,&qqnum,tetra);
	
	for(step=0;step<initial.update;step++)
	{
		if(data.autopoly==1) 
		{	
			update_P_auto(&ptr,data);
			calc_exfreq_auto(ptr,data,tetra);
		}
		else{	
			update_P_allo(&ptr,data);	
			calc_exfreq_allo(ptr,data,tetra);
		}
		update_S_POP(data,&ptr,tetra);
		
		update_ZQ(&ptr,data,0,&qqnum,tetra);
		//printf("yymm=%d\n",tetra->num_allele);
		update_geno(ptr,data,tetra);
		//printf("mm=%d\n",tetra->num_allele);
		ptr->totallkh=cal_lkd(ptr,data,tetra);
		//printf("xxmm=%d\n",tetra->num_allele);
		if(data.print_iter==1) 
		{	print_info(ptr,data,step,initial.update);}
		if(step==initial.burnin-1)
		{	allocate_chn(&chain,data);  }
		if(step>=initial.burnin&&(step+1-initial.burnin)%initial.thinning==0)		//put the node content into the chain at certain thinning interval,o.w. discard it
		{	
			store_chn(&chain,ptr,data);
			if(cnt_step<cvg->ckrep)
			{	cvg->convg_ld[chn*cvg->ckrep+cnt_step]=ptr->totallkh;}
			cnt_step++;			 
		}
		if(cnt_step==data.nstep_check_empty_cluster)
		{
			if((chain.flag_empty_cluster=check_empty_cluster(ptr,data))==1)
			{					
				fprintf(stdout,"Chain %d has an empty cluster, thus discarded!\n",chn+1);
				break;
			}
		}
	}
	free_node(ptr,data);
	free_poly(tetra,data);
	return(chain);
}


void gen_polyinfo(POLY *polyploid, SEQDATA data)
{
	int i,j;
	polyploid->num_allogeno[0]=1;
	polyploid->num_allogeno[1]=7;
	polyploid->num_allogeno[2]=12;
	polyploid->num_allogeno[3]=6;
	polyploid->num_autogeno[0]=1;
	polyploid->num_autogeno[1]=3;
	polyploid->num_autogeno[2]=3;
	polyploid->num_autogeno[3]=1;
	//printf("rr=%d\n",polyploid->num_autogeno[3]);
	gen_allele_poly(data,polyploid);
	if(data.autopoly==1)
	{
		auto_geno_num(polyploid);
		//printf("auto=%d\n",polyploid->genonum[0][0]);
		auto_geno_list(data.ploid,polyploid);
	}
	else{
		allo_geno_num(polyploid);
		allo_geno_list(data.ploid,polyploid);
	}
	if((polyploid->exfreq=(float ***)malloc(data.popnum*sizeof(float **)))==NULL)
	{	nrerror("Allocation error in Polyploid->exfreq");}
	if((polyploid->genofreq=(float ***)malloc(data.popnum*sizeof(float **)))==NULL)
	{	nrerror("Allocation error in Polyploid->exfreq");}
	for(i=0;i<data.popnum;i++)
	{
		if((polyploid->exfreq[i]=(float **)malloc(data.locinum*sizeof(float *)))==NULL)
		{	nrerror("Allocation error in Polyploid->exfreq");}
		if((polyploid->genofreq[i]=(float **)malloc(data.locinum*sizeof(float *)))==NULL)
		{	nrerror("Allocation error in Polyploid->exfreq");}
		for(j=0;j<data.locinum;j++)
		{
			if((polyploid->exfreq[i][j]=(float *)malloc(polyploid->genonum[find_id(data.allelenum[j],polyploid->allele_poly,polyploid->num_allele)][0]*sizeof(float)))==NULL)
			{	nrerror("Allocation error in Polyploid->exfreq");}
			if((polyploid->genofreq[i][j]=(float *)malloc(polyploid->genonum[find_id(data.allelenum[j],polyploid->allele_poly,polyploid->num_allele)][0]*sizeof(float)))==NULL)
			{	nrerror("Allocation error in Polyploid->exfreq");}
		}
	}
}

void free_poly(POLY *polyp, SEQDATA data)
{
	int i,j;
	
	for(i=0;i<data.popnum;i++)
	{
		for(j=0;j<data.locinum;j++)
		{
			free(polyp->exfreq[i][j]);
			free(polyp->genofreq[i][j]);
		}
		free(polyp->exfreq[i]);
		free(polyp->genofreq[i]);
	}
	free(polyp->exfreq);
	free(polyp->genofreq);
	for(i=0;i<polyp->num_allele;i++)
		free(polyp->genolist[i]);
	free(polyp->genolist);
	if(data.autopoly==1)
		free_imatrix(polyp->genonum,0,polyp->num_allele-1,0,5);
	else{	free_imatrix(polyp->genonum,0,polyp->num_allele-1,0,4);}
	free_ivector(polyp->allele_poly,0,polyp->num_allele-1);
}

/*void allocate_node(UPMCMC **ptr,SEQDATA data)

{
	if(((*ptr)=(UPMCMC *)malloc(sizeof(UPMCMC)))==NULL)
	{	nrerror("Allocation failure in ptr");}
	
	(*ptr)->self_rates=dvector(0,data.popnum-1);
	if(data.back_refl==0)
	{	(*ptr)->state=ivector(0,data.popnum-1);}	
		
	(*ptr)->indvlkh=dvector(0,data.totalsize-1);	
	(*ptr)->freq=d3tensor(0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);		
		
	(*ptr)->z=i3tensor(0,data.totalsize-1,0,data.locinum-1,0,data.ploid-1);
	(*ptr)->qq=dmatrix(0,data.totalsize-1,0,data.popnum-1);	//store the percentage of each pop. for each individual
	
}


void free_node(UPMCMC *ptr,SEQDATA data)
{
	free_dvector(ptr->self_rates,0,data.popnum-1);
	if(data.back_refl==0)
	{	free_ivector(ptr->state,0,data.popnum-1);}
		
	free_dvector(ptr->indvlkh,0,data.totalsize-1);
	free_d3tensor(ptr->freq,0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);
	
	free_i3tensor(ptr->z,0,data.totalsize-1,0,data.locinum-1,0,data.ploid-1);
	free_dmatrix(ptr->qq,0,data.totalsize-1,0,data.popnum-1);
	
	free(ptr);
}		

void allocate_chn(CHAIN *chain, SEQDATA data)

{
	int j;
	if((chain->node=(MCMC *)malloc(chain->steps*sizeof(MCMC)))==NULL)
	{	nrerror("Memory allocation for variable \'chain->node\' in function allocate_chn()!");}
	
	for(j=0;j<chain->steps;j++)
	{							
		chain->node[j].indvlkh=dvector(0,data.totalsize-1);
		chain->node[j].qq=dmatrix(0,data.totalsize-1,0,data.popnum-1);
		
		chain->node[j].self_rates=dvector(0,data.popnum-1);
		if(data.back_refl==0)
		{	chain->node[j].state=ivector(0,data.popnum-1);}
		if(data.print_freq==1)
		{	chain->node[j].freq=d3tensor(0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);	}
	}
}


void free_tetra_chain(CHAIN *chain, SEQDATA data)

{
	int j;
	free_cvector(chain->chn_name,0,chain->name_len-1);
	for(j=0;j<chain->steps;j++)
	{
		free_dvector(chain->node[j].indvlkh,0,data.totalsize-1);
		free_dmatrix(chain->node[j].qq,0,data.totalsize-1,0,data.popnum-1);
		free_dvector(chain->node[j].self_rates,0,data.popnum-1);
		if(data.back_refl==0)
		{	free_ivector(chain->node[j].state,0,data.popnum-1);}
		if(data.print_freq==1)
		{
			free_d3tensor(chain->node[j].freq,0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);
		}
	}
	free(chain->node);
}

void store_chn(CHAIN *mchain,UPMCMC *ptr,long cnt_step,SEQDATA data)

{
	int i,j,k;
	mchain->node[cnt_step].totallkh=ptr->totallkh;
	for(i=0;i<data.popnum;i++)
		mchain->node[cnt_step].self_rates[i]=ptr->self_rates[i];
	for(i=0;i<data.totalsize;i++)
	{	
		mchain->node[cnt_step].indvlkh[i]=ptr->indvlkh[i];
		for(k=0;k<data.popnum;k++)
		{
			mchain->node[cnt_step].qq[i][k]=ptr->qq[i][k];
		}
	}
	if(data.print_freq==1)
	{
		for(k=0;k<data.popnum;k++)
		{
			for(j=0;j<data.locinum;j++)
			{
				for(i=0;i<data.allelenum[j];i++)
				{
					mchain->node[cnt_step].freq[k][j][i]=ptr->freq[k][j][i];
				}
			}
		}
	}
}*/

void initial_geno(UPMCMC *ptr, SEQDATA data,POLY **polyploid)
{
	int i,j,k;
	//	printf("b%d %d\t",(*polyploid)->num_allele,(*polyploid)->num_autogeno[1]);
	if(data.autopoly==0) //allopolyploid
	{
		for(i=0;i<data.totalsize;i++)
		{
			for(j=0;j<data.locinum;j++)
			{
				if(data.missindx[i][j]!=1)
				{
					switch(data.alleleid[i][j])
					{
						case 1:	for(k=0;k<data.ploid;k++)
									ptr->geno[i][j][k]=data.seqdata[i][j][0];
								break;
						case 2: two_allele_allo(choose_unif((*polyploid)->num_allogeno[1]),i,j,ptr,data);break;
						case 3: tri_allele_allo(choose_unif((*polyploid)->num_allogeno[2]),i,j,ptr,data);break;
						case 4: tetra_allele_allo(choose_unif((*polyploid)->num_allogeno[3]),i,j,ptr,data);break;
					}
				}
			}
		}
	}
	else{
		if(data.autopoly==1) //autopolyploid
		{
			for(i=0;i<data.totalsize;i++)
			{
				for(j=0;j<data.locinum;j++)
				{
					if(data.missindx[i][j]!=1)
					{
						switch(data.alleleid[i][j])
						{
							case 1:	for(k=0;k<data.ploid;k++)
										ptr->geno[i][j][k]=data.seqdata[i][j][0];
									break;
							case 2: two_allele_auto(choose_unif((*polyploid)->num_autogeno[1]),i,j,ptr,data);break;
							case 3: tri_allele_auto(choose_unif((*polyploid)->num_autogeno[2]),i,j,ptr,data);break;
							case 4: for(k=0;k<data.ploid;k++)
										ptr->geno[i][j][k]=data.seqdata[i][j][k];
									break;
						}
						//for(k=0;k<data.ploid;k++)
							//printf("%d ",ptr->geno[i][j][k]);
					}//printf("\t");
				}//printf("\n");
			}
		}
	}
	//printf("%d\n",ptr->geno[0][0][0]);
}


void initial_chn(double ***qqnum,SEQDATA data,UPMCMC **ptr,CHAIN *chain,INIT initial,int chn)
{
	int j;
	allocate_node(ptr,data);       /*allocate space for the updated information of each step*/

	if(data.mode!=0) 
	{
		(*qqnum)=dmatrix(0,data.totalsize-1,0,data.popnum-1);//store the number of allele copies in each individual for each pop.		
		(*ptr)->alpha=ran1()*10;				/*initial sample of alpha (uniformly distributed upon [0,10])*/	
	}
	chain->name_len=initial.name_len[chn];
	chain->chn_name=cvector(0,chain->name_len-1);
	for(j=0;j<chain->name_len;j++)
		chain->chn_name[j]=initial.chn_name[chn][j];
	chain->steps=(int)((initial.update-initial.burnin)/initial.thinning);
	fprintf(stdout,"\n\n%s Starts:\n",chain->chn_name);     	
}

void update_P_auto(UPMCMC **ptr,SEQDATA data)
/*
 *updating P (allele frequency per locus per subpop) with Dirichlet distribution
 */
{
	int ***seqpop=NULL,i,j,k,l,m,n;
	double *tmp_allele=NULL,lambda=1.0;
		
	seqpop=i3tensor(0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);	//store the # of each allele type at each locus in each population
	tmp_allele=dvector(0,data.allelenum_max-1);
	
	for(n=0;n<data.popnum;n++)
		for(m=0;m<data.locinum;m++)		/*initialize "seqpop" with "0"*/
			for(l=0;l<data.allelenum[m];l++)
				seqpop[n][m][l]=0;				
													        	                        
	for(m=0;m<data.locinum;m++)		/*count the # of each allele type at each locus in each population*/
	{	for(l=0;l<data.allelenum[m];l++)
		{	for(n=0;n<data.totalsize;n++)
			{	if(data.missindx[n][m]!=1)
				{	
					for(j=0;j<data.ploid;j++)
					{
						for(i=0;i<data.popnum;i++)
						{
							if((*ptr)->z[n][m][j]==i)
							{	
								if((*ptr)->geno[n][m][j]==l) 	
									seqpop[i][m][l]++;
							}  	
						}       
					}
				}      
			}	
		}
	}
	for(i=0;i<data.popnum;i++)		/*sample frequencies of alleles at each locus in each population from dirichlet distribution*/
	{							
		for(j=0;j<data.locinum;j++)
		{							
			for(k=0;k<data.allelenum[j];k++)
			{	tmp_allele[k]=(double)seqpop[i][j][k];}
			rdirich(tmp_allele,data.allelenum[j],&((*ptr)->freq[i][j]),lambda);		
		}
	}
	
	free_i3tensor(seqpop,0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);
	free_dvector(tmp_allele,0,data.allelenum_max-1);
}


void update_P_allo(UPMCMC **ptr,SEQDATA data)
{
	int ***seqpop=NULL,***seqpop2=NULL,i,j,k,l,m,n;
	double *tmp_allele=NULL,lambda=1.0;
		
	seqpop=i3tensor(0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);	//store the # of each allele type at each locus in each population
	seqpop2=i3tensor(0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);
	tmp_allele=dvector(0,data.allelenum_max-1);
	
	for(n=0;n<data.popnum;n++)
	{	for(m=0;m<data.locinum;m++)		/*initialize "seqpop" with "0"*/
		{	for(l=0;l<data.allelenum[m];l++)
			{
				seqpop[n][m][l]=0;	
				seqpop2[n][m][l]=0;		
			}
		}
	}		
													        	                        
	for(m=0;m<data.locinum;m++)		/*count the # of each allele type at each locus in each population*/
	{	for(l=0;l<data.allelenum[m];l++)
		{	for(n=0;n<data.totalsize;n++)
			{	if(data.missindx[n][m]!=1)
				{	
					for(j=0;j<data.ploid/2;j++)
					{
						for(i=0;i<data.popnum;i++)
						{
							if((*ptr)->z[n][m][j]==i)
							{	
								if((*ptr)->geno[n][m][j]==l) 	
									seqpop[i][m][l]++;
							}  	
						}													//two systems
					}
					for(j=data.ploid/2;j<data.ploid;j++)
					{
						for(i=0;i<data.popnum;i++)
						{
							if((*ptr)->z[n][m][j]==i)
							{	
								if((*ptr)->geno[n][m][j]==l) 	
									seqpop2[i][m][l]++;
							}  	
						}       
					}
				}      
			}	
		}
	}
	for(i=0;i<data.popnum;i++)		/*sample frequencies of alleles at each locus in each population from dirichlet distribution*/
	{							
		for(j=0;j<data.locinum;j++)
		{							
			for(k=0;k<data.allelenum[j];k++)
			{	tmp_allele[k]=(double)seqpop[i][j][k];}
			rdirich(tmp_allele,data.allelenum[j],&((*ptr)->freq[i][j]),lambda);
			for(k=0;k<data.allelenum[j];k++)
			{	tmp_allele[k]=(double)seqpop2[i][j][k];}
			rdirich(tmp_allele,data.allelenum[j],&((*ptr)->freq2[i][j]),lambda);		
		}
	}
	/*for(i=0;i<data.popnum;i++)	
	{							
		for(j=0;j<data.locinum;j++)
		{							
			for(k=0;k<data.allelenum[j];k++)
			{
				printf("freq%d_%d=%f_%f\n",i,j,(*ptr)->freq[i][j][k],(*ptr)->freq2[i][j][k]);
			}
		}
	}*/
	
	free_i3tensor(seqpop,0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);
	free_i3tensor(seqpop2,0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);
	free_dvector(tmp_allele,0,data.allelenum_max-1);
}


void update_geno(UPMCMC *ptr, SEQDATA data, POLY *polypld)
{
	int i,j,k,cat_id=0;
	
	if(data.autopoly==0) //allopolyploid
	{
		for(i=0;i<data.totalsize;i++)
		{
			for(j=0;j<data.locinum;j++)
			{
				if(data.missindx[i][j]!=1)
				{//printf("d=%d\n",data.alleleid[i][j]);
					switch(data.alleleid[i][j])
					{
						case 1:	for(k=0;k<data.ploid;k++)
									ptr->geno[i][j][k]=data.seqdata[i][j][0];
								break;
						case 2: two_allele_allo(choose_two_allo(i,j,polypld->num_allogeno[1],ptr,data,polypld),i,j,ptr,data);break;
						case 3: tri_allele_allo(choose_tri_allo(i,j,polypld->num_allogeno[2],ptr,data,polypld),i,j,ptr,data);break;
						case 4: tetra_allele_allo(choose_tetra_allo(i,j,polypld->num_allogeno[3],ptr,data,polypld),i,j,ptr,data);break;
					}
					cat_id=get_cat_allo(data.ploid,ptr->geno[i][j]);
					if(check_rule_allo(data.ploid,ptr->geno[i][j],cat_id)==0)
					{
						change_geno_allo(data.ploid,&(ptr->geno[i][j]),cat_id);
					}
				}
			}
		}
	}
	else{
		if(data.autopoly==1) //autopolyploid
		{
			for(i=0;i<data.totalsize;i++)
			{
				for(j=0;j<data.locinum;j++)
				{
					if(data.missindx[i][j]!=1)
					{
						switch(data.alleleid[i][j])
						{
							case 1:	for(k=0;k<data.ploid;k++)
										ptr->geno[i][j][k]=data.seqdata[i][j][0];
									break;
							case 2: two_allele_auto(choose_two_auto(i,j,polypld->num_autogeno[1],ptr,data,polypld),i,j,ptr,data);break;
							case 3: tri_allele_auto(choose_tri_auto(i,j,polypld->num_autogeno[2],ptr,data,polypld),i,j,ptr,data);break;
							case 4: for(k=0;k<data.ploid;k++)
										ptr->geno[i][j][k]=data.seqdata[i][j][k];
									break;
						}
						cat_id=get_cat_auto(data.ploid,ptr->geno[i][j]);
						if(check_rule_auto(data.ploid,ptr->geno[i][j],cat_id)==0)
						{
							change_geno_auto(data.ploid,&(ptr->geno[i][j]),cat_id);
						}
					}
				}
			}
		}
	}
}



void update_S_POP(SEQDATA data,UPMCMC **ptr,POLY *polypld)
{
	int i,j,*tem_stat=NULL;
	double delta0=0.05,mhratio=1.0,tmp=0;
	float **genofreq_tmp;
	if((genofreq_tmp=(float **)malloc(data.locinum*sizeof(float *)))==NULL)
	{	nrerror("Allocation error in 'genofreq_tmp'");}
	for(i=0;i<data.locinum;i++)	
	{
		if((genofreq_tmp[i]=(float *)malloc(polypld->genonum[find_id(data.allelenum[i],polypld->allele_poly,polypld->num_allele)][0]*sizeof(float)))==NULL)
		{	nrerror("Allocation error in 'genofreq_tmp'");}
	}
	if(data.back_refl==0)
	{
		tem_stat=ivector(0,data.popnum-1);
	}	
	
	for(j=0;j<data.popnum;j++)	
		calc_self_genofreq((*ptr)->self_rates[j],polypld->genofreq[j],data,polypld,j);												
	for(j=0;j<data.popnum;j++)			/*updating selfing rates for each population*/
	{	
		if(data.back_refl==1)
		{
			tmp=ran1()*2*delta0-delta0;	/*propose the selfing rates*/
			tmp+=(*ptr)->self_rates[j];							
			if(tmp<=0.000)		
			{	tmp=0.000-tmp;}
			else{if(tmp>=1.000)			/*reflective bound*/
			{	tmp=1.000-(tmp-1.000);}	}			
		}
		if(data.back_refl==0)
		{
			for(i=0;i<data.popnum;i++)
			{	
				tem_stat[i]=(*ptr)->state[i];
			}
			tmp=adpt_indp(tem_stat+j,(*ptr)->state[j]);
			//fprintf(stdout,"%d %d %f\n",tem_stat[j],(*ptr)->state[j],tmp[j]);
		}		
		calc_self_genofreq(tmp,genofreq_tmp,data,polypld,j);
		
		mhratio=cal_lkd_props(j,genofreq_tmp,*ptr,data,polypld)-cal_lkd(*ptr,data,polypld);
		if(data.back_refl==0)
		{	mhratio*=hastings_stat(tem_stat,(*ptr)->state,data.popnum);	}
		//printf("propself=%f %f mhratio=%f\n",tmp,(*ptr)->self_rates[j],mhratio);
		if(ran1()<exp(MIN2(0,mhratio)))
		{
			(*ptr)->self_rates[j]=tmp;
			if(data.back_refl==0)
			{	(*ptr)->state[j]=tem_stat[j];   } 	
			move_genofreq(j,genofreq_tmp,polypld,data);	
		}
	}	
			
	if(data.back_refl==0)
	{	free_ivector(tem_stat,0,data.popnum-1);}
	for(i=0;i<data.locinum;i++)	
		free(genofreq_tmp[i]);
	free(genofreq_tmp);
}

double cal_lkd_props(int id,float **genofreq_tmp,UPMCMC *ptr,SEQDATA data,POLY *polypld)
{
	int i,j,m,cat_id=0,geno_id; //cat_id means "category"
	double ld=0;
	
	for(i=0;i<data.totalsize;i++)
	{
		for(j=0;j<data.locinum;j++)
		{
			if(data.missindx[i][j]!=1)
			{	
				if(data.autopoly==1)
				{	geno_id=get_index_auto(j,data.ploid,ptr->geno[i][j],polypld,data.allelenum[j],&cat_id);}
				else{	geno_id=get_index_allo(j,data.ploid,ptr->geno[i][j],polypld,data.allelenum[j],&cat_id);}
					
				if(chcksame(ptr->z[i][j],data.ploid)==0)
				{
					
					if(id==ptr->z[i][j][0])
					{	ld+=(double)genofreq_tmp[j][geno_id];
					//for(m=0;m<data.ploid;m++)
					//printf("%d ",ptr->geno[i][j][m]);
					//printf(" geno_%d_%d_%d=%d,%f %f\n",i,j,cat_id,geno_id,genofreq_tmp[j][geno_id],polypld->genofreq[ptr->z[i][j][0]][j][geno_id]);
					}
					else{	ld+=(double)polypld->genofreq[ptr->z[i][j][0]][j][geno_id];}
				}
				else{
					if(data.autopoly==1)
					{
						for(m=0;m<data.ploid;m++)
						{
							ld+=log(ptr->freq[ptr->z[i][j][m]][j][ptr->geno[i][j][m]]);
						}		
						switch(cat_id) //heterozygote needs to time 2
						{	
							case 1: ld+=log(4);break;
							case 2: ld+=log(6);break;
							case 3: ld+=log(12);break;
							case 4: ld+=log(24);break;
						}
					}
					else{
						if(data.autopoly==0)
						{
							for(m=0;m<(data.ploid/2);m++)
							{
								ld+=log(ptr->freq[ptr->z[i][j][m]][j][ptr->geno[i][j][m]]);
							}
							for(m=(data.ploid/2);m<data.ploid;m++)
							{
								ld+=log(ptr->freq2[ptr->z[i][j][m]][j][ptr->geno[i][j][m]]);
							}
							switch(cat_id) //heterozygote needs to time 2
							{	
								case 1: ld+=log(2);break;
								case 2: ld+=log(2);break;
								case 3: ld+=log(4);break;
							}
						}
					}
				}
			}
		} 
	}
	//printf("ld=%f\n",ld);
	return(ld);
}



double cal_lkd(UPMCMC *ptr,SEQDATA data,POLY *polypld)
{
	int i,j; //cat_id means "category"
	double ld=0,sum=0;
	
	for(i=0;i<data.totalsize;i++)
	{
		ld=0;
		for(j=0;j<data.locinum;j++)
		{
			if(data.missindx[i][j]!=1)
			{	
				ld+=calc_genofq(j,i,ptr->z[i][j],data,ptr,polypld);
			}
		}
		ptr->indvlkh[i]=ld;
		sum+=ld; 
	}
	//printf("sum=%f\n",sum);
	return(sum);
}

void move_genofreq(int id, float **genofreq_tmp,POLY *polypld,SEQDATA data)
{
	int i,j;
	for(i=0;i<data.locinum;i++)
	{
		for(j=0;j<polypld->genonum[find_id(data.allelenum[i],polypld->allele_poly,polypld->num_allele)][0];j++)
		{
			polypld->genofreq[id][i][j]=genofreq_tmp[i][j];
		}
	}
}


void update_ZQ(UPMCMC **ptr,SEQDATA data,int init_flag,double ***qqnum,POLY *polypld)
/*
 * Update Z and Q according to the Metropolis-Hastings ratio
 */
{	
	int i,j,k,m,l,temp,*tmp1;
	double *tmp=NULL,mhratio=0;
	tmp=dvector(0,data.popnum-1);
	tmp1=ivector(0,data.ploid-1);
	
	for(i=0;i<data.totalsize;i++)
	{
		for(j=0;j<data.locinum;j++)		/*updating Z*/
		{
			if(data.missindx[i][j]!=1)
			{
				for(k=0;k<data.ploid;k++)
				{
					for(m=0;m<data.popnum;m++)
					{	
						if(init_flag==1)
						{	tmp[m]=(double)(m+1)/data.popnum;}
						else{							
							tmp[m]=(*ptr)->qq[i][m]*(*ptr)->freq[m][j][(*ptr)->geno[i][j][k]];
							if(m>=1) {tmp[m]+=tmp[m-1];}
						}
					}
					temp=disc_unif(tmp,data.popnum);
					//if(init_flag==1)
					//{	
					(*ptr)->z[i][j][k]=temp;
					/*}
					else{
						for(m=0;m<data.ploid;m++)
						{
							tmp1[m]=(*ptr)->z[i][j][m];
						}
						tmp1[k]=temp;
						mhratio=log(calc_genofq(j,i,tmp1,data,(*ptr),polypld))-log(calc_genofq(j,i,(*ptr)->z[i][j],data,(*ptr),polypld));
						if(ran1()<exp(MIN2(0,mhratio)))
						{
							(*ptr)->z[i][j][k]=temp;
						}
					}
					for(l=0;l<data.ploid;l++)
							{
								tmp1[l]=(*ptr)->z[i][j][l];
							}
							tmp1[k]=m;
						tmpp=tmp[0];
						for(m=0;m<data.popnum;m++)
						{
							tmp[m]=exp(tmp[m]-tmpp);
							if(m>=1) {tmp[m]+=tmp[m-1];}
						}*/
				}
			}
		}

		for(m=0;m<data.popnum;m++)/*updating Q*/
			(*qqnum)[i][m]=0.0;
		for(l=0;l<data.locinum;l++)		
		{
			if(data.missindx[i][l]!=1)
			{
				for(k=0;k<data.ploid;k++)
				{
					for(m=0;m<data.popnum;m++)
					{
						if((*ptr)->z[i][l][k]==m)
						{
							(*qqnum)[i][m]+=1.0;
							break;
						}
					}		
				}
			}
		}

		for(k=0;k<data.popnum;k++)
		{	tmp[k]=(*qqnum)[i][k];}
		rdirich(tmp,data.popnum,&((*ptr)->qq[i]),(*ptr)->alpha);
	}
	free_dvector(tmp,0,data.popnum-1);
	free_ivector(tmp1,0,data.ploid-1);
	
}



int choose_unif(int temp)
//temp=the number of combination of diploid genotypes.  
{
	int j,num;
	double *tmp;
	//printf("temp=%d\n",temp);
	tmp=dvector(0,temp-1);
	for(j=0;j<temp;j++)	
	{	tmp[j]=(double)(j+1)/(double)temp;}
	num=disc_unif(tmp,temp)+1;
	free_dvector(tmp,0,temp-1);
	return(num);
}

int choose_two_auto(int indv_id,int loci_id,int temp,UPMCMC *ptr,SEQDATA data, POLY *polypld)
{
	int i,j,*num,index,n=0,n_type=2,*id,nid;
	double *tmp,*freq,tm;
	tmp=dvector(0,temp-1);
	freq=dvector(0,n_type-1);
	id=ivector(0,n_type-1);
	num=ivector(0,temp-1);
	
	n=data.allelenum[loci_id];
	nid=find_id(n,polypld->allele_poly,polypld->num_allele);
	for(i=0;i<n_type;i++)
	{	id[i]=data.seqdata[indv_id][loci_id][i];}
	
	if(chcksame(ptr->z[indv_id][loci_id],data.ploid)==0)
	{		
		num[0]=id[0]*n*(n*n+n+1)+id[1];	//A1A1A1A2
				
		num[1]=id[1]*n*(n*n+n+1)+id[0];	//A2A2A2A1	
			
		num[2]=(id[0]*n*n+id[1])*(n+1);	//A1A1A2A2
		
		for(i=0;i<temp;i++)
			tmp[i]=(double)polypld->genofreq[ptr->z[indv_id][loci_id][0]][loci_id][find_id(num[i],polypld->genolist[nid],polypld->genonum[nid][0])];
	}
	else{
		for(i=0;i<n_type;i++)
		{
			freq[i]=0;
			for(j=0;j<data.popnum;j++)
				freq[i]+=ptr->qq[indv_id][j]*ptr->freq[j][loci_id][id[i]];
		}
		
		tmp[0]=log(4)+3*log(freq[0])+log(freq[1]);		//A1A1A1A2
		
		tmp[1]=log(4)+3*log(freq[1])+log(freq[0]);		//A2A2A2A1
		
		tmp[2]=log(6)+2*log(freq[0])+2*log(freq[1]);	//A1A1A2A2
	}
	tm=tmp[0];
	for(i=0;i<temp;i++)
		tmp[i]=exp(tmp[i]-tm);
	for(i=1;i<temp;i++)
		tmp[i]+=tmp[i-1];
		
	index=disc_unif(tmp,temp)+1;
	free_dvector(tmp,0,temp-1);
	free_dvector(freq,0,n_type-1);
	free_ivector(id,0,n_type-1);
	free_ivector(num,0,temp-1);
	return(index);
}

int choose_tri_auto(int indv_id,int loci_id,int temp,UPMCMC *ptr,SEQDATA data, POLY *polypld)
{
	int i,j,*num,index,n=0,n_type=3,*id,nid;
	double *tmp,*freq,tm;
	tmp=dvector(0,temp-1);
	freq=dvector(0,n_type-1);
	id=ivector(0,n_type-1);
	num=ivector(0,temp-1);
	
	n=data.allelenum[loci_id];
	nid=find_id(n,polypld->allele_poly,polypld->num_allele);
	for(i=0;i<n_type;i++)
	{	id[i]=data.seqdata[indv_id][loci_id][i];}
	
	if(chcksame(ptr->z[indv_id][loci_id],data.ploid)==0)
	{		
		num[0]=id[0]*n*n*(n+1)+id[1]*n+id[2];	//A1A1A2A3
		
		num[1]=id[1]*n*n*(n+1)+id[0]*n+id[2];	//A2A2A1A3
		
		num[2]=id[2]*n*n*(n+1)+id[0]*n+id[1];	//A3A3A1A2
		
		for(i=0;i<temp;i++)
			tmp[i]=(double)polypld->genofreq[ptr->z[indv_id][loci_id][0]][loci_id][find_id(num[i],polypld->genolist[nid],polypld->genonum[nid][0])];

	}
	else{
		for(i=0;i<n_type;i++)
		{
			freq[i]=0;
			for(j=0;j<data.popnum;j++)
				freq[i]+=ptr->qq[indv_id][j]*ptr->freq[j][loci_id][id[i]];
		}
		
		tmp[0]=2*log(freq[0])+log(freq[1])+log(freq[2]);	//A1A1A2A3
		
		tmp[1]=2*log(freq[1])+log(freq[0])+log(freq[2]);	//A2A2A1A3
		
		tmp[2]=2*log(freq[2])+log(freq[1])+log(freq[0]);	//A3A3A1A2
	}
	
	tm=tmp[0];
	for(i=0;i<temp;i++)
		tmp[i]=exp(tmp[i]-tm);
	for(i=1;i<temp;i++)
		tmp[i]+=tmp[i-1];
		
	index=disc_unif(tmp,temp)+1;
	free_dvector(tmp,0,temp-1);
	free_dvector(freq,0,n_type-1);
	free_ivector(id,0,n_type-1);
	free_ivector(num,0,temp-1);
	return(index);
}

int choose_two_allo(int indv_id,int loci_id,int temp,UPMCMC *ptr,SEQDATA data, POLY *polypld)
//When {A1,A2}, the possible combinations are 
//{1=A1A1+A1A2, 2=A1A2+A1A1, 3=A1A1+A2A2, 4=A2A2+A1A1, 5=A1A2+A2A2, 6=A2A2+A1A2, 7=A1A2+A1A2}
{
	int i,j,*num,index,n=0,n_type=2,*id,nid;
	double *tmp,*freq,*freq2,tm;
	tmp=dvector(0,temp-1);
	freq=dvector(0,n_type-1);
	freq2=dvector(0,n_type-1);
	id=ivector(0,n_type-1);
	num=ivector(0,temp-1);
		
	n=data.allelenum[loci_id];
	nid=find_id(n,polypld->allele_poly,polypld->num_allele);
	for(i=0;i<n_type;i++)
	{	id[i]=data.seqdata[indv_id][loci_id][i];}
	
	if(chcksame(ptr->z[indv_id][loci_id],data.ploid)==0)	//if this locus come from the same subpop.
	{		
		num[0]=id[0]*n*(n*n+n+1)+id[1];		//A1A1+A1A2
		
		num[1]=id[0]*(n*n*n+n+1)+id[1]*n*n;	//A1A2+A1A1
		
		num[2]=id[0]*n*n*(n+1)+id[1]*(n+1);	//A1A1+A2A2
		
		num[3]=id[0]*(n+1)+id[1]*n*n*(n+1);	//A2A2+A1A1
		
		num[4]=id[0]*n*n*n+id[1]*(n*n+n+1);	//A1A2+A2A2
		
		num[5]=id[0]*n+id[1]*(n*n*n+n*n+1);	//A2A2+A1A2
		
		num[6]=id[0]*n*(n*n+1)+id[1]*(n*n+1);	//A1A2+A1A2
		
		for(i=0;i<temp;i++)
			tmp[i]=(double)polypld->genofreq[ptr->z[indv_id][loci_id][0]][loci_id][find_id(num[i],polypld->genolist[nid],polypld->genonum[nid][0])];	
	}
	else{
		for(i=0;i<n_type;i++)
		{
			freq[i]=0;
			freq2[i]=0;
			for(j=0;j<data.popnum;j++)
			{
				freq[i]+=ptr->qq[indv_id][j]*ptr->freq[j][loci_id][id[i]];
				freq2[i]+=ptr->qq[indv_id][j]*ptr->freq2[j][loci_id][id[i]];
			}
		}
		
		tmp[0]=2*log(freq[0])+log(freq2[0])+log(freq2[1]);	//A1A1+A1A2
		
		tmp[1]=log(freq[0])+log(freq[1])+2*log(freq2[0]);	//A1A2+A1A1
		
		tmp[2]=2*log(freq[0])+2*log(freq2[1]);				//A1A1+A2A2
		
		tmp[3]=2*log(freq[1])+2*log(freq2[0]);				//A2A2+A1A1
		
		tmp[4]=log(freq[0])+log(freq[1])+2*log(freq2[1]);	//A1A2+A2A2
		
		tmp[5]=2*log(freq[1])+log(freq2[0])+log(freq2[1]);	//A2A2+A1A2
		
		tmp[6]=log(2)+log(freq[0])+log(freq[1])+log(freq2[0])+log(freq2[1]);	//A1A2+A1A2
	}
	
	tm=tmp[0];
	for(i=0;i<temp;i++)
	{	tmp[i]=exp(tmp[i]-tm);
		//printf("tmp=%f\t",tmp[i]);
	}//printf("\n");
	for(i=1;i<temp;i++)
		tmp[i]+=tmp[i-1];
		
	index=disc_unif(tmp,temp)+1;
	//printf("index=%d\n",index);
	free_dvector(tmp,0,temp-1);
	free_dvector(freq,0,n_type-1);
	free_dvector(freq2,0,n_type-1);
	free_ivector(id,0,n_type-1);
	free_ivector(num,0,temp-1);
	return(index);
}

int choose_tri_allo(int indv_id,int loci_id,int temp,UPMCMC *ptr,SEQDATA data, POLY *polypld)
//{1=A1A1+A2A3, 2=A2A3+A1A1, 3=A2A2+A1A3, 4=A1A3+A2A2, 5=A3A3+A1A2, 6=A1A2+A3A3,
//7=A1A2+A2A3, 8=A2A3+A1A2, 9=A2A3+A1A3, 10=A1A3+A2A3, 11=A1A3+A1A2, 12=A1A2+A1A3}
{
	int i,j,*num,index,n=0,n_type=3,*id,nid;
	double *tmp,*freq,*freq2,tm;
	tmp=dvector(0,temp-1);
	freq=dvector(0,n_type-1);
	freq2=dvector(0,n_type-1);
	id=ivector(0,n_type-1);
	num=ivector(0,temp-1);
		
	n=data.allelenum[loci_id];
	nid=find_id(n,polypld->allele_poly,polypld->num_allele);
	for(i=0;i<n_type;i++)
	{	id[i]=data.seqdata[indv_id][loci_id][i];}
	
	if(chcksame(ptr->z[indv_id][loci_id],data.ploid)==0)	//if this locus come from the same subpop.
	{		
		num[0]=id[0]*n*n*(n+1)+id[1]*n+id[2];		//A1A1+A2A3
		
		num[1]=id[0]*(n+1)+id[1]*n*n*n+id[2]*n*n;	//A2A3+A1A1
		
		num[2]=id[0]*n+id[1]*n*n*(n+1)+id[2];		//A2A2+A1A3
		
		num[3]=id[0]*n*n*n+id[1]*(n+1)+id[2]*n*n;	//A1A3+A2A2
		
		num[4]=id[0]*n+id[1]+id[2]*n*n*(n+1);		//A3A3+A1A2
		
		num[5]=id[0]*n*n*n+id[1]*n*n+id[2]*(n+1);	//A1A2+A3A3
		
		num[6]=id[0]*n*n*n+id[1]*n*(n+1)+id[2];		//A1A2+A2A3
		
		num[7]=id[0]*n+id[1]*(n*n*n+1)+id[2]*n*n;	//A2A3+A1A2
		
		num[8]=id[0]*n+id[1]*n*n*n+id[2]*(n*n+1);	//A2A3+A1A3
		
		num[9]=id[0]*n*n*n+id[1]*n+id[2]*(n*n+1);	//A1A3+A2A3
		
		num[10]=id[0]*n*(n*n+1)+id[1]+id[2]*n*n;	//A1A3+A1A2
		
		num[11]=id[0]*n*(n*n+1)+id[1]*n*n+id[2];	//A1A2+A1A3
		
		for(i=0;i<temp;i++)
			tmp[i]=(double)polypld->genofreq[ptr->z[indv_id][loci_id][0]][loci_id][find_id(num[i],polypld->genolist[nid],polypld->genonum[nid][0])];	
	}
	else{
		for(i=0;i<n_type;i++)
		{
			freq[i]=0;
			freq2[i]=0;
			for(j=0;j<data.popnum;j++)
			{
				freq[i]+=ptr->qq[indv_id][j]*ptr->freq[j][loci_id][id[i]];
				freq2[i]+=ptr->qq[indv_id][j]*ptr->freq2[j][loci_id][id[i]];
			}
		}
		
		tmp[0]=2*log(freq[0])+log(freq2[1])+log(freq2[2]);				//A1A1+A2A3
		
		tmp[1]=log(freq[1])+log(freq[2])+2*log(freq2[0]);				//A2A3+A1A1
		
		tmp[2]=2*log(freq[1])+log(freq2[0])+log(freq2[2]);				//A2A2+A1A3
		
		tmp[3]=log(freq[0])+log(freq[2])+2*log(freq2[1]);				//A1A3+A2A2
		
		tmp[4]=2*log(freq[2])+log(freq2[0])+log(freq2[1]);				//A3A3+A1A2
		
		tmp[5]=log(freq[0])+log(freq[1])+2*log(freq2[2]);				//A1A2+A3A3
		
		tmp[6]=log(2)+log(freq[0])+log(freq[1])+log(freq2[1])+log(freq2[2]);	//A1A2+A2A3
		
		tmp[7]=log(2)+log(freq[1])+log(freq[2])+log(freq2[0])+log(freq2[1]);	//A2A3+A1A2
		
		tmp[8]=log(2)+log(freq[1])+log(freq[2])+log(freq2[0])+log(freq2[2]);	//A2A3+A1A3
		
		tmp[9]=log(2)+log(freq[0])+log(freq[2])+log(freq2[1])+log(freq2[2]);	//A1A3+A2A3
		
		tmp[10]=log(2)+log(freq[0])+log(freq[2])+log(freq2[0])+log(freq2[1]);	//A1A3+A1A2
		
		tmp[11]=log(2)+log(freq[0])+log(freq[1])+log(freq2[0])+log(freq2[2]);	//A1A2+A1A3
	}
	
	tm=tmp[0];
	for(i=0;i<temp;i++)
	{	tmp[i]=exp(tmp[i]-tm);
		//printf("tmp=%f\t",tmp[i]);
	}//printf("\n");
	for(i=1;i<temp;i++)
		tmp[i]+=tmp[i-1];
		
	index=disc_unif(tmp,temp)+1;
	//printf("index=%d\n",index);
	free_dvector(tmp,0,temp-1);
	free_dvector(freq,0,n_type-1);
	free_dvector(freq2,0,n_type-1);
	free_ivector(id,0,n_type-1);
	free_ivector(num,0,temp-1);
	return(index);
}

int choose_tetra_allo(int indv_id,int loci_id,int temp,UPMCMC *ptr,SEQDATA data, POLY *polypld)
//{1=A1A2+A3A4, 2=A3A4+A1A2, 3=A1A3+A2A4, 4=A2A4+A1A3, 5=A1A4+A2A3, 6=A2A3+A1A4}
{
	int i,j,*num,index,n=0,n_type=4,*id,nid;
	double *tmp,*freq,*freq2,tm;
	tmp=dvector(0,temp-1);
	freq=dvector(0,n_type-1);
	freq2=dvector(0,n_type-1);
	id=ivector(0,n_type-1);
	num=ivector(0,temp-1);
		
	n=data.allelenum[loci_id];
	nid=find_id(n,polypld->allele_poly,polypld->num_allele);
	for(i=0;i<n_type;i++)
	{	id[i]=data.seqdata[indv_id][loci_id][i];}
	
	if(chcksame(ptr->z[indv_id][loci_id],data.ploid)==0)	//if this locus come from the same subpop.
	{		
		num[0]=id[0]*n*n*n+id[1]*n*n+id[2]*n+id[3];	//A1A2+A3A4
		
		num[1]=id[2]*n*n*n+id[3]*n*n+id[0]*n+id[1];	//A3A4+A1A2
		
		num[2]=id[0]*n*n*n+id[2]*n*n+id[1]*n+id[3];	//A1A3+A2A4
		
		num[3]=id[1]*n*n*n+id[3]*n*n+id[0]*n+id[2];	//A2A4+A1A3
		
		num[4]=id[0]*n*n*n+id[3]*n*n+id[1]*n+id[2];	//A1A4+A2A3
		
		num[5]=id[1]*n*n*n+id[2]*n*n+id[0]*n+id[3];	//A2A3+A1A4
		
		for(i=0;i<temp;i++)
			tmp[i]=(double)polypld->genofreq[ptr->z[indv_id][loci_id][0]][loci_id][find_id(num[i],polypld->genolist[nid],polypld->genonum[nid][0])];	
	}
	else{
		for(i=0;i<n_type;i++)
		{
			freq[i]=0;
			freq2[i]=0;
			for(j=0;j<data.popnum;j++)
			{
				freq[i]+=ptr->qq[indv_id][j]*ptr->freq[j][loci_id][id[i]];
				freq2[i]+=ptr->qq[indv_id][j]*ptr->freq2[j][loci_id][id[i]];
			}
		}
		
		tmp[0]=log(freq[0])+log(freq[1])+log(freq2[2])+log(freq2[3]);	//A1A2+A3A4
		
		tmp[1]=log(freq[2])+log(freq[3])+log(freq2[0])+log(freq2[1]);	//A3A4+A1A2
		
		tmp[2]=log(freq[0])+log(freq[2])+log(freq2[1])+log(freq2[3]);	//A1A3+A2A4
		
		tmp[3]=log(freq[1])+log(freq[3])+log(freq2[0])+log(freq2[2]);	//A2A4+A1A3
		
		tmp[4]=log(freq[0])+log(freq[3])+log(freq2[1])+log(freq2[2]);	//A1A4+A2A3
		
		tmp[5]=log(freq[1])+log(freq[2])+log(freq2[0])+log(freq2[3]);	//A2A3+A1A4
	}
	
	tm=tmp[0];
	for(i=0;i<temp;i++)
		tmp[i]=exp(tmp[i]-tm);
	for(i=1;i<temp;i++)
		tmp[i]+=tmp[i-1];
		
	index=disc_unif(tmp,temp)+1;
	free_dvector(tmp,0,temp-1);
	free_dvector(freq,0,n_type-1);
	free_dvector(freq2,0,n_type-1);
	free_ivector(id,0,n_type-1);
	free_ivector(num,0,temp-1);
	return(index);
}



void calc_self_genofreq(double self_rate,float **genofreq,SEQDATA data,POLY *polypld,int k)
//For each subpopulation for each locus calculate the genotype frequencies
{
	int i;
	
	for(i=0;i<data.locinum;i++)
	{
		if(data.autopoly==1)
		{	auto_genfreq((float)self_rate, k,i,data,polypld, genofreq[i]);}
		else{
			if(data.autopoly==0)
			{	allo_genfreq((float)self_rate, k,i,data,polypld, genofreq[i]);}
		}
	}
}

double calc_genofq(int loci_id,int indv_id,int *z, SEQDATA data, UPMCMC *ptr,POLY *polypld)
//calculate the genotype frequencies given the assignment of each allele
{
	int geno_id,cat_id,m;
	double ld=0;
	
	if(data.missindx[indv_id][loci_id]==1) return(0);
	if(data.autopoly==1)
	{	geno_id=get_index_auto(loci_id,data.ploid,ptr->geno[indv_id][loci_id],polypld,data.allelenum[loci_id],&cat_id);}
	else{	geno_id=get_index_allo(loci_id,data.ploid,ptr->geno[indv_id][loci_id],polypld,data.allelenum[loci_id],&cat_id);}
	if(chcksame(z,data.ploid)==0)
	{
		ld=(double)polypld->genofreq[z[0]][loci_id][geno_id];
	}
	else{
		if(data.autopoly==1)
		{
			for(m=0;m<data.ploid;m++)
			{
				ld+=log(ptr->freq[z[m]][loci_id][ptr->geno[indv_id][loci_id][m]]);
			}		
			switch(cat_id) //heterozygote needs to time 2
			{	
				case 1: ld+=log(4);break;
				case 2: ld+=log(6);break;
				case 3: ld+=log(12);break;
				case 4: ld+=log(24);break;
			}
		}
		else{
			if(data.autopoly==0)
			{
				for(m=0;m<(data.ploid/2);m++)
				{
					ld+=log(ptr->freq[z[m]][loci_id][ptr->geno[indv_id][loci_id][m]]);
				}
				for(m=(data.ploid/2);m<data.ploid;m++)
				{
					ld+=log(ptr->freq2[z[m]][loci_id][ptr->geno[indv_id][loci_id][m]]);
				}
			
				switch(cat_id) //heterozygote needs to time 2
				{	
					case 1: ld+=log(2);break;
					case 2: ld+=log(2);break;
					case 3: ld+=log(4);break;
				}
			}
		}
	}
	return(ld);
}


int get_index_auto(int locus_id, int ploid,int *geno,POLY *polypld,int allelenum,int *cat_id)
{
	int i,res=0,temp=0,num;
	
	*cat_id=get_cat_auto(ploid, geno);
	/*for(i=0;i<ploid;i++)
		printf("%d ",geno[i]);
	printf("\tcat_id=%d %d\n",*cat_id,allelenum);*/
	if(check_rule_auto(ploid,geno,*cat_id)==1) //following the rules of writing genotypes
	{
		temp=geno[0];
		for(i=1;i<ploid;i++)
		{
			temp*=allelenum;
			temp+=geno[i];
		}
	
		num=find_id(allelenum,polypld->allele_poly,polypld->num_allele);
		res=find_id(temp,polypld->genolist[num],polypld->genonum[num][0]);	
		return(res);
	}
	else{	return(-1);nrerror("genoypes are not following the rules!");}
}

int get_cat_auto(int ploid,int *geno)
//return the index of category that genotype in autoployploid belongs to
{
	int i,cnt=0,cat_id=0,*tmp;
	tmp=ivector(0,ploid-1);
	
	tmp[0]=geno[0];
	cnt++;
	for(i=1;i<ploid;i++)
	{
		if(exists(geno[i],tmp,cnt)==0) 
		{	tmp[cnt++]=geno[i];}
	}
	switch(cnt)
	{
		case 1: cat_id=0;break;
		case 2:	if(copy_num(tmp[0],geno,ploid)==2) 
				{	cat_id=2;}
				else{	cat_id=1;}
				break;
		case 3: cat_id=3;break;
		case 4: cat_id=4;break;
	}
		//printf("%d %d\n",cnt,cat_id);
	free_ivector(tmp,0,ploid-1);
	return(cat_id);
}

int get_cat_allo(int ploid,int *geno)
//return the index of category that genotype in alloployploid belongs to
{
	int i,cnt1=0,cnt2=0,cat_id=0,*tmp;
	tmp=ivector(0,ploid-1);
	
	tmp[0]=geno[0];
	cnt1++;	//count the number of allele types at the first half loci, the first two alleles in a tetraploid
	for(i=1;i<ploid/2;i++)
	{
		if(exists(geno[i],tmp,cnt1)==0) 
		{	tmp[cnt1++]=geno[i];}
	}
	tmp[0]=geno[ploid/2];
	cnt2++;	//count the number of allele types at the first half loci, the first two alleles in a tetraploid
	for(i=ploid/2+1;i<ploid;i++)
	{
		if(exists(geno[i],tmp,cnt2)==0) 
		{	tmp[cnt2++]=geno[i];}
	}
	switch(cnt1+cnt2)
	{
		case 2: cat_id=0;break;
		case 3:	if(cnt1==1&&cnt2==2) 
				{	cat_id=1;}
				else{	cat_id=2;}
				break;
		case 4: cat_id=3;break;
	}
	free_ivector(tmp,0,ploid-1);
	return(cat_id);
}

int get_index_allo(int locus_id, int ploid,int *geno,POLY *polypld,int allelenum,int *cat_id)
{
	int i,res=0,temp=0,num;
	
	*cat_id=get_cat_allo(ploid, geno);
	
	if(check_rule_allo(ploid,geno,*cat_id)==1) //following the rules of writing genotypes
	{
		temp=geno[0];
		for(i=1;i<ploid;i++)
		{
			temp*=allelenum;
			temp+=geno[i];
		}
	
		num=find_id(allelenum,polypld->allele_poly,polypld->num_allele);
		res=find_id(temp,polypld->genolist[num],polypld->genonum[num][0]);
		return(res);
	}
	else{	return(-1);nrerror("genoypes are not following the rules!");}
}

int check_rule_auto(int ploid, int *geno,int cat_id)
//check whether current genotypes follow the rules of writing genotypes
{
	int flag=0;
	
	switch(cat_id)
	{	
		case 0: if(chcksame(geno,ploid)==0) 
				{	flag=1;}
				break;
		case 1:	if(chcksame(geno,ploid-1)==0&&geno[ploid-1]!=geno[0])
				{	flag=1;}
				break;
		case 2: if(chcksame(geno,ploid/2)==0&&chcksame(geno+ploid/2,ploid/2)==0&&geno[0]<geno[ploid-1])
				{	flag=1;}
				break;
		case 3: if(chcksame(geno,ploid/2)==0&&geno[ploid-2]<geno[ploid-1])
				{	flag=1;}
				break;
		case 4: if(geno[0]<geno[1]&&geno[1]<geno[2]&&geno[2]<geno[3])
				{	flag=1;}
				break;
	}
	
	return(flag);
}

void change_geno_auto(int ploid, int **geno,int cat_id)
//check whether current genotypes follow the rules of writing genotypes
{
	int temp;
	switch(cat_id)
	{	
		case 0: if(chcksame(*geno,ploid)==1) 
				{	nrerror("Bad genotypes are generated");}
				break;
		case 1:	if((*geno)[ploid-1]==(*geno)[0])
				{	nrerror("Bad genotypes are generated");}
				break;
		case 2: if((*geno)[0]>(*geno)[ploid-1])
				{	
					SWAP((*geno)[0],(*geno)[ploid-1])
					(*geno)[1]=(*geno)[0];
					(*geno)[ploid-2]=(*geno)[ploid-1];
				}
				break;
		case 3: if((*geno)[ploid-2]>(*geno)[ploid-1])
				{	SWAP((*geno)[ploid-2],(*geno)[ploid-1])}
				break;
		case 4: sort(ploid,*geno);
				break;
	}
}


int check_rule_allo(int ploid, int *geno,int cat_id)
//check whether current genotypes follow the rules of writing genotypes
{
	int flag=0;
	
	switch(cat_id)
	{	
		case 0: if(geno[0]==geno[1]&&geno[ploid-2]==geno[ploid-1])//chcksame(geno,ploid/2)==0&&chcksame(geno+ploid/2,ploid/2)==0)
				{	flag=1;}
				break;
		case 1: if(geno[0]==geno[1]&&geno[ploid-2]<geno[ploid-1])
				{	flag=1;}
				break;
		case 2: if(geno[0]<geno[1]&&geno[ploid-2]==geno[ploid-1])
				{	flag=1;}
				break;
		case 3: if(geno[0]<geno[1]&&geno[ploid-2]<geno[ploid-1])
				{	flag=1;}
				break;
	}
	
	return(flag);
}

void change_geno_allo(int ploid, int **geno,int cat_id)
//check whether current genotypes follow the rules of writing genotypes
{	
	int temp;
	switch(cat_id)
	{	
		case 0: if(chcksame(*geno,ploid/2)==1||chcksame(*geno+ploid/2,ploid/2)==1) 
				{	nrerror("Bad genotypes are generated");}
				break;
		case 1:	if(chcksame(*geno,ploid/2)==1)
				{	nrerror("Bad genotypes are generated");}
				if((*geno)[ploid-2]>(*geno)[ploid-1])
				{	SWAP((*geno)[ploid-2],(*geno)[ploid-1])}
				break;
		case 2: if(chcksame(*geno+ploid/2,ploid/2)==1)
				{	nrerror("Bad genotypes are generated");}
				if((*geno)[0]>(*geno)[1])
				{	SWAP((*geno)[0],(*geno)[1])}
				break;
		case 3: if((*geno)[0]>(*geno)[1])
				{	SWAP((*geno)[0],(*geno)[1])}
				if((*geno)[ploid-2]>(*geno)[ploid-1])
				{	SWAP((*geno)[ploid-2],(*geno)[ploid-1])}
				break;
	}
}

int copy_num(int val, int *vec, int leng)
//count the times that "val" appears in the vector "vec"
{
	int i,num=0;
	for(i=0;i<leng;i++)
	{
		if(val==vec[i])
		{	num++;}
	}
	return(num);
}


void calc_exfreq_auto(UPMCMC *ptr,SEQDATA data,POLY *polypld)
//calculate the expected genotype frequency
{
	int i,j,k,m,n,*digit,n_allele=0,temp=0,tmp=0;
	digit=ivector(0,data.ploid-1);
	
	for(k=0;k<data.popnum;k++)
	{
		for(i=0;i<data.locinum;i++)
		{
			n=data.allelenum[i];
			n_allele=find_id(n,polypld->allele_poly,polypld->num_allele);
			for(j=0;j<polypld->genonum[n_allele][1];j++) //Monoallele
			{
				tmp=polypld->genolist[n_allele][j];
				digit[0]=tmp%n;
				polypld->exfreq[k][i][j]=(float)log(ptr->freq[k][i][digit[0]])*(float)data.ploid;
			}
			temp=polypld->genonum[n_allele][1];
			for(j=temp;j<temp+polypld->genonum[n_allele][2];j++) //Biallele, simplex
			{
				tmp=polypld->genolist[n_allele][j];
				digit[0]=tmp%n;
				tmp/=n;
				digit[1]=tmp%n;
				polypld->exfreq[k][i][j]=(float)(log(4.0)+log(ptr->freq[k][i][digit[1]])*(float)(data.ploid-1)+log(ptr->freq[k][i][digit[0]]));
			}
			temp+=polypld->genonum[n_allele][2];
			for(j=temp;j<temp+polypld->genonum[n_allele][3];j++) //Biallele, duplex
			{
				tmp=polypld->genolist[n_allele][j];
				digit[0]=tmp%n;
				tmp/=(n*n);
				digit[1]=tmp%n;
				polypld->exfreq[k][i][j]=(float)(log(6.0)+(log(ptr->freq[k][i][digit[1]])+log(ptr->freq[k][i][digit[0]]))*(data.ploid/2));			
			}
			temp+=polypld->genonum[n_allele][3];
			for(j=temp;j<temp+polypld->genonum[n_allele][4];j++) //Triallele
			{
				tmp=polypld->genolist[n_allele][j];
				for(m=0;m<data.ploid-1;m++)
				{
					digit[m]=tmp%n;
					tmp/=n;
				}				
				polypld->exfreq[k][i][j]=(float)(log(12.0)+log(ptr->freq[k][i][digit[2]])*(data.ploid/2)+log(ptr->freq[k][i][digit[0]])+log(ptr->freq[k][i][digit[1]]));			
			}
			temp+=polypld->genonum[n_allele][4];
			for(j=temp;j<temp+polypld->genonum[n_allele][5];j++)	//Quadriallele
			{
				tmp=polypld->genolist[n_allele][j];
				for(m=0;m<data.ploid;m++)
				{
					digit[m]=tmp%n;
					tmp/=n;
				}
				polypld->exfreq[k][i][j]=(float)log(24.0);
				for(m=0;m<data.ploid;m++)
				{	polypld->exfreq[k][i][j]+=(float)log(ptr->freq[k][i][digit[m]]);}
			}
		}
	}
	/*for(k=0;k<data.popnum;k++)
	{
		for(i=0;i<data.locinum;i++)
		{
			printf("k=%d, i=%d ",k,i);
			for(j=0;j<polypld->genonum[n_allele][0];j++)	
			{
				printf("%f ",polypld->exfreq[k][i][j]);
			}
			printf("\n");
		}
	}*/
	free_ivector(digit,0,data.ploid-1);
}

void calc_exfreq_allo(UPMCMC *ptr,SEQDATA data,POLY *polypld)
//calculate the expected genotype frequency
{
	int i,j,k,m,n,*digit,n_allele=0,temp=0,tmp=0;
	digit=ivector(0,data.ploid-1);
	
	for(k=0;k<data.popnum;k++)
	{
		for(i=0;i<data.locinum;i++)
		{
			n=data.allelenum[i];
			n_allele=find_id(n,polypld->allele_poly,polypld->num_allele);
			for(j=0;j<polypld->genonum[n_allele][1];j++) //iikk
			{
				tmp=polypld->genolist[n_allele][j];
				digit[0]=tmp%n;
				tmp/=(n*n);
				digit[1]=tmp%n;
				//printf("%f_%f\n ",ptr->freq[k][i][digit[1]],ptr->freq2[k][i][digit[0]]);
				polypld->exfreq[k][i][j]=(float)((log(ptr->freq[k][i][digit[1]])+log(ptr->freq2[k][i][digit[0]]))*(data.ploid/2));			
			//printf("f%f ",polypld->exfreq[k][i][j]);
			}
			temp=polypld->genonum[n_allele][1];
			for(j=temp;j<temp+polypld->genonum[n_allele][2];j++) //iikl
			{
				tmp=polypld->genolist[n_allele][j];
				for(m=0;m<data.ploid-1;m++)
				{
					digit[m]=tmp%n;
					tmp/=n;
				}				
				polypld->exfreq[k][i][j]=(float)(log(2.0)+log(ptr->freq[k][i][digit[2]])*(data.ploid/2)+log(ptr->freq2[k][i][digit[0]])+log(ptr->freq2[k][i][digit[1]]));			
			//printf("f2%f ",polypld->exfreq[k][i][j]);
			}
			temp+=polypld->genonum[n_allele][2];
			for(j=temp;j<temp+polypld->genonum[n_allele][3];j++) //ijkk
			{
				tmp=polypld->genolist[n_allele][j];
				tmp/=n;
				for(m=0;m<data.ploid-1;m++)
				{
					digit[m]=tmp%n;
					tmp/=n;
				}				
				polypld->exfreq[k][i][j]=(float)(log(2.0)+log(ptr->freq2[k][i][digit[0]])*(data.ploid/2)+log(ptr->freq[k][i][digit[2]])+log(ptr->freq[k][i][digit[1]]));			
			//printf("f3%f ",polypld->exfreq[k][i][j]);
			}
			temp+=polypld->genonum[n_allele][3];
			for(j=temp;j<temp+polypld->genonum[n_allele][4];j++)	//ijkl
			{
				tmp=polypld->genolist[n_allele][j];
				for(m=0;m<data.ploid;m++)
				{
					digit[m]=tmp%n;
					tmp/=n;
				}
				polypld->exfreq[k][i][j]=(float)log(4.0);
				for(m=0;m<data.ploid/2;m++)
				{	polypld->exfreq[k][i][j]+=(float)log(ptr->freq2[k][i][digit[m]]);}
				for(m=data.ploid/2;m<data.ploid;m++)
				{	polypld->exfreq[k][i][j]+=(float)log(ptr->freq[k][i][digit[m]]);}
			//printf("f4%f ",polypld->exfreq[k][i][j]);
			}
		}
	}
	/*for(k=0;k<data.popnum;k++)
	{
		for(i=0;i<data.locinum;i++)
		{
			printf("k=%d, i=%d ",k,i);
			for(j=0;j<polypld->genonum[n_allele][0];j++)	
			{
				printf("%f ",polypld->exfreq[k][i][j]);
			}
			printf("\n");
		}
	}*/
	free_ivector(digit,0,data.ploid-1);
}


void gen_allele_poly(SEQDATA data, POLY *polyp)	
//get the distinct number of alleles across loci
{
	int i,*tmp,cnt=0;
	tmp=ivector(0,data.allelenum_max-1);
	
	tmp[0]=data.allelenum[0];
	cnt=1;
	for(i=1;i<data.locinum;i++)
	{
		if(exists(data.allelenum[i],tmp,cnt)==0)
		{
			tmp[cnt]=data.allelenum[i];
			cnt++;
		}
	}
	polyp->num_allele=cnt;
	sort(cnt,tmp);
	polyp->allele_poly=ivector(0,polyp->num_allele-1);
	for(i=0;i<cnt;i++)
		polyp->allele_poly[i]=tmp[i];
	free_ivector(tmp,0,data.allelenum_max-1);
}


void auto_geno_num(POLY *polyp)
//Calculate the number of genotypes Generate the genotype list for each of the values of "allele_num" from 2 to "max_allele_num", the largest number of alleles at one locus
{
	int i,j;
	polyp->genonum=imatrix(0,polyp->num_allele,0,5);
	
	for(j=0;j<polyp->num_allele;j++)
	{
		i=polyp->allele_poly[j];
		//number of genotypes 
		polyp->genonum[j][1]=i; //0=Monoallele, iiii
		polyp->genonum[j][2]=i*(i-1); //1=Biallele (simplex), iiij
		polyp->genonum[j][3]=i*(i-1)/2; //2=Biallele (duplex), iijj
		polyp->genonum[j][4]=i*(i-1)*(i-2)/2; //3=Triallele, iijk
		polyp->genonum[j][5]=i*(i-1)*(i-2)*(i-3)/24; //4=Quadriallele, ijkl
		polyp->genonum[j][0]=i+i*(i-1)*3/2+i*(i-1)*(i-2)/2+i*(i-1)*(i-2)*(i-3)/24; //total number of genotypes
	}
}


void auto_geno_list(int ploid,POLY *polyp)
{
	int l,i,j,k,m,n,cnt,tmp;	
	
	if((polyp->genolist=(int **)malloc(polyp->num_allele*sizeof(int *)))==NULL)
	{	nrerror("Allocation error in 'genolist'");}
	for(i=0;i<polyp->num_allele;i++)	//5 categories of genotypes
	{
		if((polyp->genolist[i]=(int *)malloc(polyp->genonum[i][0]*sizeof(int)))==NULL)
		{	nrerror("Allocation error in 'genolist'");}
	}
	
	for(l=0;l<polyp->num_allele;l++)
	{
		i=polyp->allele_poly[l];
		//0=Monoallele, iiii
		for(j=0;j<polyp->genonum[l][1];j++)
		{
			polyp->genolist[l][j]=j*(i*i*i+i*i+i+1);
			//printf("%d",genolist[i][0][j]);
		}
		tmp=polyp->genonum[l][1];
		//1=Biallele (simplex), iiij
		cnt=0;
		for(j=0;j<i-1;j++)
		{
			for(k=j+1;k<i;k++)
			{
				polyp->genolist[l][tmp+2*cnt]=j*(i*i*i+i*i+i)+k;
				polyp->genolist[l][tmp+2*cnt+1]=i*(i*i+i+1)*k+j;
				//printf("%d",genolist[i][1][2*cnt]);
				cnt++;
			}
		}
		tmp+=polyp->genonum[l][2];
		//2=Biallele (duplex), iijj
		cnt=0;
		for(j=0;j<i-1;j++)
		{
			for(k=j+1;k<i;k++)
			{
				polyp->genolist[l][tmp+cnt]=j*(i*i*i+i*i)+k*(i+1);
				//printf("%d",genolist[i][2][cnt]);
				cnt++;
			}
		}
		tmp+=polyp->genonum[l][3];
		//3=Triallele, iijk
		cnt=0;
		for(j=0;j<i-2;j++)
		{
			for(k=j+1;k<i-1;k++)
			{
				for(m=k+1;m<i;m++)
				{
					polyp->genolist[l][tmp+3*cnt]=j*(i*i*i+i*i)+k*i+m;
					polyp->genolist[l][tmp+3*cnt+1]=k*(i*i*i+i*i)+j*i+m;
					polyp->genolist[l][tmp+3*cnt+2]=m*(i*i*i+i*i)+j*i+k;
					//printf("%d",genolist[i][3][3*cnt]);
					cnt++;
				}
			}
		}
		tmp+=polyp->genonum[l][4];
		//4=Quadriallele, ijkl
		cnt=0;
		for(j=0;j<i-3;j++)
		{
			for(k=j+1;k<i-2;k++)
			{
				for(m=k+1;m<i-1;m++)
				{
					for(n=m+1;n<i;n++)
					{
						polyp->genolist[l][tmp+cnt]=j*i*i*i+k*i*i+i*m+n;
						//printf("%d",genolist[i][4][cnt]);
						cnt++;
					}
				}
			}
		}
	}
}


void auto_genfreq(float self, int pop_id,int loci_id,SEQDATA data, POLY *polyp, float *freq)
//solve the linear equation (I-sA)P=(1-s)R
//s= "self", selfing rate
//R is the expected genotype frequencies
//freq is the outputed genotype frequencies taking into account selfing
{
	int i,j,k,l,n,id,tmp,*digit,num,tri=3;
	float temp,**matr,**vec;
	matr=matrix(1,tri,1,tri);
	vec=matrix(1,tri,1,1);
	digit=ivector(0,tri-1);
	n=data.allelenum[loci_id];
	id=find_id(n,polyp->allele_poly,polyp->num_allele);
	
	tmp=polyp->genonum[id][0];
	//Quadriallele,ijkl
	if(n>=4)
	{
		for(i=tmp-polyp->genonum[id][5];i<tmp;i++)
		{
			freq[i]=log(1-self)+polyp->exfreq[pop_id][loci_id][i]-log(1-self/6);
			if(freq[i]>0)
			{	
				printf("freq4=%f %d\n",freq[i],loci_id);
				nrerror("Genotype frequencies can not be greater than 1!");
			}
		}
	}
	
	//Triallele,iijk
	if(n>=3)
	{	
		tmp-=polyp->genonum[id][5];
	//printf("genonum=%d\n",polyp->genonum[id][4]);
		for(i=0;i<polyp->genonum[id][4]/tri;i++)
		{
			num=polyp->genolist[id][tmp-polyp->genonum[id][4]+i*3];
			for(j=data.ploid-2;j>=0;j--)
			{					//currently we can know that values in "digit" are ascendly ordered.
				digit[j]=num%n;
				num/=n;
			}
			temp=0;	//calculate \sum_{l} P_{ijkl}
			if(n>=4)
			{
				for(l=0;l<n;l++)
				{
					if(exists(l,digit,tri)==0)	
					{	
						num=find_id(calc_val(digit,l,n),polyp->genolist[id],polyp->genonum[id][0]);
						temp+=exp(freq[num]);
					}
				}//printf("tempm=%f %d\n",temp,loci_id);
				if(temp>1)
				{	
					nrerror("Genotype frequencies can not be greater than 1!");
				}
			}
			for(j=1;j<=tri;j++)
			{
				for(k=1;k<=tri;k++)
				{
					if(j==k)	matr[j][k]=1-self*10.0/36.0;
					else{	matr[j][k]=-self/9.0;}
				}
				vec[j][1]=self/18.0*temp+(1.0-self)*exp(polyp->exfreq[pop_id][loci_id][tmp-polyp->genonum[id][4]+i*3+j-1]);
				//printf("%f %f %f\n",vec[j][1],temp,polyp->exfreq[pop_id][loci_id][tmp-polyp->genonum[id][4]+i*3+j-1]);
			}
			temp=vec[1][1];
			for(j=1;j<=tri;j++)
			{	vec[j][1]/=temp;}
			gaussj(matr,tri,vec,1);
			for(j=0;j<tri;j++)
			{	//printf("vec2=%f %f\n",log(vec[j+1][1]),log(temp));
				freq[tmp-polyp->genonum[id][4]+i*3+j]=log(vec[j+1][1])+log(temp);
				//printf("freq3=%f %d\n",freq[tmp-polyp->genonum[id][4]+i*3+j],loci_id);
				if(freq[tmp-polyp->genonum[id][4]+i*3+j]>0)
				{	
					nrerror("Genotype frequencies can not be greater than 1!");
				}
			}
		}	
	}
	
	//Biallele,duplex, iijj
	tmp-=polyp->genonum[id][4];
	for(i=tmp-polyp->genonum[id][3];i<tmp;i++)
	{
		num=polyp->genolist[id][i];
		digit[0]=num%n; //we know digit[0]>digit[1]
		num/=(n*n);
		digit[1]=num%n;
		
		temp=0;
		if(n>=3)
		{
			for(j=0;j<n;j++)
			{
				if(exists(j,digit,2)==0)
				{
					if(digit[0]<j)	//P_iijk
					{	num=find_id(digit[1]*n*n*(n+1)+digit[0]*n+j,polyp->genolist[id],polyp->genonum[id][0]);}
					else{	if(digit[0]>j)	num=find_id(digit[1]*n*n*(n+1)+j*n+digit[0],polyp->genolist[id],polyp->genonum[id][0]);}
					temp+=exp(freq[num])/9.0*self;
					//printf("%f temp1=%f\n",freq[num],temp);
					if(digit[1]<j)	//P_jjik
					{	num=find_id(digit[0]*n*n*(n+1)+digit[1]*n+j,polyp->genolist[id],polyp->genonum[id][0]);}
					else{	if(digit[1]>j)	num=find_id(digit[0]*n*n*(n+1)+j*n+digit[1],polyp->genolist[id],polyp->genonum[id][0]);}
					temp+=exp(freq[num])/9.0*self;
					//printf("%f temp2=%f\n",freq[num],temp);
					num=find_id(j*n*n*(n+1)+digit[1]*n+digit[0],polyp->genolist[id],polyp->genonum[id][0]);
					temp+=exp(freq[num])/36.0*self;	//P_kkij
					//printf("%f temp3=%f\n",freq[num],temp);
					if(n>=4)
					{
						for(k=j+1;k<n;k++)
						{
							if(exists(k,digit,2)==0)
							{
								//printf("%d %d %d %d %d %d %d\n",j,k,n,digit[0],digit[1],num,calc_val2(digit,2,j,k,n));
								num=find_id(calc_val2(digit,2,j,k,n),polyp->genolist[id],polyp->genonum[id][0]);
								temp+=exp(freq[num])/36.0*self;	//P_ijkl
							//printf("%f temp4=%f\n",freq[num],temp);
							}
						}
					}
				}
			}
		}		
		freq[i]=log((1-self)*exp(polyp->exfreq[pop_id][loci_id][i])+temp)-log(1-self/2.0);
		
		if(freq[i]>0)
		{	
			nrerror("Genotype frequencies can not be greater than 1!");
		}
	}
	
	//Biallele,simplex, iiij
	tmp-=polyp->genonum[id][3];
	for(i=tmp-polyp->genonum[id][2];i<tmp;i++)
	{
		num=polyp->genolist[id][i];
		digit[0]=num%n;
		num/=n;
		digit[1]=num%n;
		if(digit[0]<digit[1])
		{	num=find_id((digit[0]*n*n+digit[1])*(n+1),polyp->genolist[id],polyp->genonum[id][0]);}
		else{
			if(digit[0]>digit[1]) num=find_id((digit[1]*n*n+digit[0])*(n+1),polyp->genolist[id],polyp->genonum[id][0]);
		}
		temp=8.0/36.0*exp(freq[num])*self;	//s8P_iijj/36
		
		if(n>=3)
		{
			for(j=0;j<n;j++)	//s4\sum_k P_iijk/36
			{
				if(exists(j,digit,2)==0)
				{
					if(digit[0]<j)
					{	num=find_id(digit[1]*n*n*(n+1)+digit[0]*n+j,polyp->genolist[id],polyp->genonum[id][0]);}
					else{	
						if(digit[0]>j)	
							num=find_id(digit[1]*n*n*(n+1)+j*n+digit[0],polyp->genolist[id],polyp->genonum[id][0]);
					}
					temp+=exp(freq[num])/9.0*self;
				}
			}
		}
		freq[i]=log((1-self)*exp(polyp->exfreq[pop_id][loci_id][i])+temp)-log(1-self/2.0);
		if(freq[i]>0)
		{	printf("freq1=%f %d\n",freq[i],loci_id);
			nrerror("Genotype frequencies can not be greater than 1!");
		}
	}
	
	//Monoallele, iiii
	tmp-=polyp->genonum[id][2];
	for(i=tmp-polyp->genonum[id][1];i<tmp;i++)
	{
		num=polyp->genolist[id][i];
		digit[0]=num%n;
		temp=0;
		for(j=0;j<n;j++)
		{
			if(j!=digit[0])	
			{
				num=find_id(digit[0]*n*(n*n+n+1)+j,polyp->genolist[id],polyp->genonum[id][0]);
				temp+=exp(freq[num])/4.0*self;	//sP_iiij/4.0
				
				if(digit[0]<j)
				{	num=find_id(digit[0]*n*n*(n+1)+j*(n+1),polyp->genolist[id],polyp->genonum[id][0]);}
				else{
					if(digit[0]<j)
					{	num=find_id(j*n*n*(n+1)+digit[0]*(n+1),polyp->genolist[id],polyp->genonum[id][0]);}
				}
				temp+=exp(freq[num])/36.0*self;	//sP_iijj/36.0
				
				if(n>=3)
				{
					for(k=j+1;k<n;k++)
					{
						if(k!=digit[0])
						{
							num=find_id(digit[0]*n*n*(n+1)+j*n+k,polyp->genolist[id],polyp->genonum[id][0]);
							temp+=exp(freq[num])/36.0*self;	//sP_iijk/36.0
						}
					}
				}
			}		
		}
		freq[i]=log((1-self)*exp(polyp->exfreq[pop_id][loci_id][i])+temp)-log(1-self);
		if(freq[i]>0)
		{	printf("freq0=%f %d\n",freq[i],loci_id);
			nrerror("Genotype frequencies can not be greater than 1!");
		}
	}
	for(i=0;i<polyp->genonum[id][0];i++)
	{
		if(freq[i]>0)
		{	printf("freq=%f %d\n",freq[i],loci_id);
		nrerror("Genotype frequencies can not be greater than 1!");}
	}
	free_ivector(digit,0,data.ploid-2);
	free_matrix(matr,1,tri,1,tri);
	free_matrix(vec,1,tri,1,1);
}


void allo_geno_num(POLY *polyp)
//Calculate the number of genotypes Generate the genotype list for each of the values of "allele_num" from 2 to "max_allele_num", the largest number of alleles at one locus
{
	int i,j;
	polyp->genonum=imatrix(0,polyp->num_allele,0,4);
	
	for(j=0;j<polyp->num_allele;j++)
	{
		i=polyp->allele_poly[j];
		//number of genotypes 
		polyp->genonum[j][1]=i*i; //0=Biallele, iikk
		polyp->genonum[j][2]=i*(i-1)/2*i; //1=Triallele1, iikl
		polyp->genonum[j][3]=i*(i-1)/2*i; //2=Triallele2, ijkk
		polyp->genonum[j][4]=i*(i-1)*i*(i-1)/4; //3=Tetraallele, ijkl
		polyp->genonum[j][0]=i*i+i*(i-1)*i+i*(i-1)*i*(i-1)/4;//sum
	}
}


void allo_geno_list(int ploid,POLY *polyp)
{
	int l,i,j,k,m,n,cnt,tmp;	
	
	if((polyp->genolist=(int **)malloc(polyp->num_allele*sizeof(int *)))==NULL)
	{	nrerror("Allocation error in 'genolist'");}
	for(i=0;i<polyp->num_allele;i++)	//4 categories of genotypes
	{
		if((polyp->genolist[i]=(int *)malloc(polyp->genonum[i][0]*sizeof(int)))==NULL)
		{	nrerror("Allocation error in 'genolist'");}
	}
	
	for(l=0;l<polyp->num_allele;l++)
	{
		i=polyp->allele_poly[l];
		//0=Biallele, iikk
		for(j=0;j<i;j++)
		{
			for(k=0;k<i;k++)
			{
				polyp->genolist[l][j*i+k]=j*i*i*(i+1)+k*(i+1);
			}
		}
		tmp=polyp->genonum[l][1];
		//1=Triallele1, iikl
		cnt=0;
		for(j=0;j<i;j++)
		{
			for(k=0;k<i-1;k++)
			{
				for(m=k+1;m<i;m++)
				{
					polyp->genolist[l][tmp+cnt]=j*i*i*(i+1)+i*k+m;
					cnt++;
				}
			}
		}
		tmp+=polyp->genonum[l][2];
		//2=Biallele (duplex), ijkk
		cnt=0;
		for(j=0;j<i-1;j++)
		{
			for(k=j+1;k<i;k++)
			{
				for(m=0;m<i;m++)
				{
					polyp->genolist[l][tmp+cnt]=(j*i+k)*i*i+m*(i+1);
					cnt++;
				}
			}
		}
		tmp+=polyp->genonum[l][3];
		//3=Tetrallele, ijkl
		cnt=0;
		for(j=0;j<i-1;j++)
		{
			for(k=j+1;k<i;k++)
			{
				for(m=0;m<i-1;m++)
				{
					for(n=m+1;n<i;n++)
					{
						polyp->genolist[l][tmp+cnt]=j*i*i*i+k*i*i+m*i+n;
						cnt++;
					}
				}
			}
		}
	}
}


void allo_genfreq(float self, int pop_id,int loci_id,SEQDATA data, POLY *polyp, float *freq)
//solve the linear equation (I-sA)P=(1-s)R
//s= "self", selfing rate
//R is the expected genotype frequencies
//freq is the outputed genotype frequencies taking into account selfing
{
	int i,j,k,n,id,tmp,*digit,num,tri=3;
	float temp;

	digit=ivector(0,tri-1);
	n=data.allelenum[loci_id];
	id=find_id(n,polyp->allele_poly,polyp->num_allele);
	
	//Quadriallele,ijkl
	tmp=polyp->genonum[id][0];
	for(i=tmp-polyp->genonum[id][4];i<tmp;i++)
		{
			freq[i]=log(1-self)+polyp->exfreq[pop_id][loci_id][i]-log(1-self/4);
			
			if(freq[i]>0)
			{	printf("freq4=%f %d\n",freq[i],loci_id);
				nrerror("Genotype frequencies can not be greater than 1!");
			}
			//printf("freq4=%f\n",freq[i]);
		}
	
	
	//Triallele 2, ijkk
	tmp-=polyp->genonum[id][4];
		for(i=tmp-polyp->genonum[id][3];i<tmp;i++)
		{
			num=polyp->genolist[id][i];
			digit[0]=num%n; 
			num/=(n*n);
			digit[1]=num%n;
			num/=n;
			digit[2]=num%n;
		
			temp=0;
			for(j=0;j<n;j++)
				{
					if(j!=digit[0])
					{
						if(digit[0]<j)	
						{	num=find_id(digit[2]*n*n*n+digit[1]*n*n+digit[0]*n+j,polyp->genolist[id],polyp->genonum[id][0]);}
						else{
							if(digit[0]>j)	
							{	num=find_id(digit[2]*n*n*n+digit[1]*n*n+j*n+digit[0],polyp->genolist[id],polyp->genonum[id][0]);}
						}
						temp+=exp(freq[num])*self/8.0;
					}
				}
			freq[i]=log((1-self)*exp(polyp->exfreq[pop_id][loci_id][i])+temp)-log(1-self/2.0);
			
			if(freq[i]>0)
			{	printf("freq3=%f %d\n",freq[i],loci_id);
				nrerror("Genotype frequencies can not be greater than 1!");
			}

		}
	
	
	//Triallele 1, iikl
	tmp-=polyp->genonum[id][3];
		for(i=tmp-polyp->genonum[id][2];i<tmp;i++)
		{
			num=polyp->genolist[id][i];
			digit[0]=num%n; 
			num/=n;
			digit[1]=num%n;
			num/=n;
			digit[2]=num%n;
		
			temp=0;
			for(j=0;j<n;j++)
				{
					if(j!=digit[2])
					{
						if(digit[2]<j)	
						{	num=find_id(digit[2]*n*n*n+j*n*n+digit[1]*n+digit[0],polyp->genolist[id],polyp->genonum[id][0]);}
						else{
							if(digit[2]>j)	
							{	num=find_id(j*n*n*n+digit[2]*n*n+digit[1]*n+digit[0],polyp->genolist[id],polyp->genonum[id][0]);}
						}
						temp+=exp(freq[num])*self/8.0;
					}
				}
			
			freq[i]=log((1-self)*exp(polyp->exfreq[pop_id][loci_id][i])+temp)-log(1-self/2.0);
			
			if(freq[i]>0)
			{	printf("freq2=%f %d\n",freq[i],loci_id);
				nrerror("Genotype frequencies can not be greater than 1!");
			}

		}
	
	
	//Biallele, iikk
	tmp-=polyp->genonum[id][2];
	for(i=tmp-polyp->genonum[id][1];i<tmp;i++)
	{
		num=polyp->genolist[id][i];
		digit[0]=num%n; 
		num/=(n*n);
		digit[1]=num%n;
		
		temp=0;			//P_iikl
			for(j=0;j<n;j++)
			{
				if(j!=digit[0])
				{
					if(digit[0]<j)	
					{	num=find_id(digit[1]*n*n*(n+1)+digit[0]*n+j,polyp->genolist[id],polyp->genonum[id][0]);}
					else{
						if(digit[0]>j)	
						{	num=find_id(digit[1]*n*n*(n+1)+j*n+digit[0],polyp->genolist[id],polyp->genonum[id][0]);}
					}
					temp+=exp(freq[num])*self/4.0;
				}
			}		
			for(j=0;j<n;j++)	//P_ijkk
			{
				if(j!=digit[1])
				{
					if(digit[1]<j)	
					{	num=find_id(digit[1]*n*n*n+j*n*n+digit[0]*(n+1),polyp->genolist[id],polyp->genonum[id][0]);}
					else{
						if(digit[1]>j)	
						{	num=find_id(j*n*n*n+digit[1]*n*n+digit[0]*(n+1),polyp->genolist[id],polyp->genonum[id][0]);}
					}
					temp+=exp(freq[num])*self/4.0;
				}
			}
		for(j=0;j<n;j++)	//P_ijkl
			{
				for(k=0;k<n;k++)
				{
					if(j!=digit[1]&&k!=digit[0])
					{
						if(digit[1]<j)	
						{	
							if(digit[0]<k) 
								num=find_id(digit[1]*n*n*n+j*n*n+digit[0]*n+k,polyp->genolist[id],polyp->genonum[id][0]);
							else{
								if(digit[0]>k) 
									num=find_id(digit[1]*n*n*n+j*n*n+k*n+digit[0],polyp->genolist[id],polyp->genonum[id][0]);
							}
						}
						else{
							if(digit[1]>j)	
							{	
								if(digit[0]<k) 
									num=find_id(j*n*n*n+digit[1]*n*n+digit[0]*n+k,polyp->genolist[id],polyp->genonum[id][0]);
								else{
									if(digit[0]>k) 
										num=find_id(j*n*n*n+digit[1]*n*n+k*n+digit[0],polyp->genolist[id],polyp->genonum[id][0]);
								}
							}
						}
						temp+=exp(freq[num])*self/16.0;
					}
				}
			}
		
		freq[i]=log((1-self)*exp(polyp->exfreq[pop_id][loci_id][i])+temp)-log(1-self);
		
		if(freq[i]>0)
		{	
			printf("freq1=%f %d\n",freq[i],loci_id);
			nrerror("Genotype frequencies can not be greater than 1!");
		}

	}
	for(i=0;i<polyp->genonum[id][0];i++)
	{
		if(freq[i]>0)
		{	printf("freq=%f %d\n",freq[i],loci_id);
			nrerror("mmmmGenotype frequencies can not be greater than 1!");
		}
	}
	free_ivector(digit,0,data.ploid-2);
}


int calc_val(int *num, int val,int i)
//we know num is ascendly ordered with three elements.
{
	int temp=0;
	
	if(val<num[0])
	{	temp=val*i*i*i+num[0]*i*i+num[1]*i+num[2];}
	else{
		if(val>num[0]&&val<num[1])
		{
			temp=num[0]*i*i*i+val*i*i+num[1]*i+num[2];
		}
		else{
			if(val>num[1]&&val<num[2])
			{
				temp=num[0]*i*i*i+num[1]*i*i+val*i+num[2];
			}
			else{
				if(val>num[2])
				{	temp=num[0]*i*i*i+num[1]*i*i+num[2]*i+val;}
			}
		}
	}
	return(temp);
}

int calc_val2(int *num, int leng,int val1,int val2,int i)
//we know num is descendly ordered with two elements and val1<val2
{
	int temp=0;
	
	if(val2<num[1])//val1<val2<num[1]<num[0]
	{	temp=val1*i*i*i+val2*i*i+num[1]*i+num[0];}
	else{
		if(val2>num[1]&&val2<num[0]&&val1<num[1])//val1<num[1]<val2<num[0]	
		{
			temp=val1*i*i*i+num[1]*i*i+val2*i+num[0];
		}
		else{
			if(val1>num[1]&&val2<num[0])	//num[1]<val1<val2<num[0]
			{
				temp=num[1]*i*i*i+val1*i*i+val2*i+num[0];
			}
			else{				
				if(val1>num[1]&&val1<num[0]&&val2>num[0])	//num[1]<val1<num[0]<val2
				{	temp=num[1]*i*i*i+val1*i*i+num[0]*i+val2;}
				else{
					if(val1>num[0])	//num[1]<num[0]<val1<val2
					{	temp=num[1]*i*i*i+num[0]*i*i+val1*i+val2;}
					else{
						if(val1<num[1]&&val2>num[0])	//val1<num[1]<num[0]<val2
						{	temp=val1*i*i*i+num[1]*i*i+i*num[0]+val2;}
					}
				}
			}
		}
	}
	return(temp);
}

int find_id(int num, int *array,int array_len)
{
	int i,flag=0;

	for(i=0;i<array_len;i++)
	{
		if(array[i]==num)
		{
			flag=1;
			break;
		}
	}
	if(flag==0) nrerror("Cannot find the index in the genotype list");
	return(i);
}


void gaussj(float **a, int n, float **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	float big,dum,pivinv,temp;

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}




void two_allele_auto(int num, int i,int locus_id, UPMCMC *ptr,SEQDATA data)
//When {A1,A2}, the possible combinations are 
//{1=A1A1A1A2, 2=A2A2A2A1, 3=A1A1A2A2
{
	switch(num)
	{
		case 1: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][1];
				break;
		case 2: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][0];
				break;
		case 3: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][1];
				break;
	}
}


void two_allele_allo(int num,  int i,int locus_id, UPMCMC *ptr,SEQDATA data)
//When {A1,A2}, the possible combinations are 
//{1=A1A1+A1A2, 2=A1A2+A1A1, 3=A1A1+A2A2, 4=A2A2+A1A1, 5=A1A2+A2A2, 6=A2A2+A1A2, 7=A1A2+A1A2}
{
	switch(num)
	{
		case 1: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][1];
				break;
		case 2: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][0];
				break;
		case 3: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][1];
				break;
		case 4: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][0];
				break;
		case 5: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][1];
				break;
		case 6: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][1];
				break;
		case 7: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][1];
				break;
	}
}

void tri_allele_auto(int num,  int i,int locus_id, UPMCMC *ptr,SEQDATA data)
//{1=A1A1A2A3, 2=A2A2A1A3, 3=A3A3A1A2		
{
	switch(num)
	{
		case 1: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][2];
				break;
		case 2: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][2];
				break;				
		case 3: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][2];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][2];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][1];
				break;
	}
}


void tri_allele_allo(int num,  int i,int locus_id, UPMCMC *ptr,SEQDATA data)
//{1=A1A1+A2A3, 2=A2A3+A1A1, 3=A2A2+A1A3, 4=A1A3+A2A2, 5=A3A3+A1A2, 6=A1A2+A3A3,
//7=A1A2+A2A3, 8=A2A3+A1A2, 9=A2A3+A1A3, 10=A1A3+A2A3, 11=A1A3+A1A2, 12=A1A2+A1A3}			
{
	switch(num)
	{
		case 1: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][2];
				break;
		case 2: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][2];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][0];
				break;				
		case 3: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][2];
				break;
		case 4: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][2];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][1];
				break;
		case 5: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][2];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][2];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][1];
				break;
		case 6: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][2];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][2];
				break;
		case 7: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][2];
				break;
		case 8: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][2];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][1];
				break;
		case 9: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][2];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][2];
				break;
		case 10:ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][2];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][2];
				break;
		case 11:ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][2];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][1];
				break;
		case 12:ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][2];
				break;
	}
}

void tetra_allele_allo(int num,  int i, int locus_id, UPMCMC *ptr,SEQDATA data)
//{1=A1A2+A3A4, 2=A3A4+A1A2, 3=A1A3+A2A4, 4=A2A4+A1A3, 5=A1A4+A2A3, 6=A2A3+A1A4}
{
	switch(num)
	{
		case 1: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][2];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][3];
				break;
		case 2: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][2];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][3];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][1];
				break;
		case 3: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][2];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][3];
				break;
		case 4: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][3];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][2];
				break;
		case 5: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][3];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][2];
				break;
		case 6: ptr->geno[i][locus_id][0]=data.seqdata[i][locus_id][1];
				ptr->geno[i][locus_id][1]=data.seqdata[i][locus_id][2];
				ptr->geno[i][locus_id][2]=data.seqdata[i][locus_id][0];
				ptr->geno[i][locus_id][3]=data.seqdata[i][locus_id][3];
				break;
	}
}

/*void convmatr(float self, int allele_num,int geno_num, float **matr, float *genofreq,float **Bfreq)
//Solve the matrix linear equation: Ax=b, (I-sA)P=(1-s)R
//geno_num is the dimension of matrix "matr" and vector "genofreq"
//"R" is the expected genotype frequencies
//The result of "P" is returned in "Bfreq".
{
	int i, j;
	float **Amat;Bfreq;
	Amat=matrix(0,geno_num-1,0,geno_num-1);
	
	for(i=0;i<geno_num;i++)
	{
		for(j=0;j<geno_num;j++)
		{
			if(i==j)
			{	Amat[i][j]=1-self*matr[i][j];}
			else{	Amat[i][j]=-self*matr[i][j];}
		}
	}
	
	for(i=0;i<geno_num;i++)
	{
		Bfreq[i][0]=(1-self)*genofreq[i];
	}
	
	gaussj(Amat, geno_num, Bfreq, 1);
	
	free_matrix(Amat,0,geno_num-1,0,geno_num-1);
	//return(Bfreq);
}
*/
	/*if((polyp->selfmat=(float ***)malloc(polyp->num_allele*sizeof(float **)))==NULL)
	{	nrerror("Allocation error in 'selfmat'");}	
	for(i=0;i<polyp->num_allele;i++)
	{
		if((selfmat[i]=(float **)malloc(polyp->genonum[i][0]*sizeof(float *)))==NULL)
		{	nrerror("Allocation error in 'selfmat'");}
		for(j=0;j<genonum[i-2][5];j++)
		{
			if((selfmat[i-2][j]=(float *)malloc(polyp->genonum[i][0]*sizeof(float)))==NULL)
			{	nrerror("Allocation error in 'selfmat'");}
		}
		
		for(j=0;j<genonum[i-2][5];j++)
		{
			for(k=0;k<genonum[i-2][5];k++)
			{
				selfmat[i-2][j][k]=0.0;
			}
		}
		
		//Coefficients for Monoallele
			for(j=0;j<genonum[i-2][0];j++)
			{
				selfmat[i-2][j][j]=1;	//Monoallele
				for(k=0;k<i;k++)
				{
					if(k!=j)//Biallele,simplex
					{
						temp=find_id((j*(i*i*i+i*i+i)+k), genolist[i-2],genonum[i-2][5]);
						selfmat[i-2][j][temp]=1.0/4.0;
					}
					if(k<j)
					{	
						temp2=find_id((k*(i*i*i+i*i)+j*(i+1)), genolist[i-2],genonum[i-2][5]);//Biallele,duplex
						selfmat[i-2][j][temp2]=1.0/36.0;
					}
					else{
						if(k>j)	
						{	
							temp2=find_id((j*(i*i*i+i*i)+k*(i+1)), genolist[i-2],genonum[i-2][5]);//Biallele,duplex
							selfmat[i-2][j][temp2]=1.0/36.0;
						}
					}
				}
				for(k=0;k<i;k++)	//Triallele
				{
					for(m=k+1;m<i;m++)
					{
						if(k!=j&&m!=j)
						{
							temp=find_id((j*(i*i*i+i*i)+k*i+m), genolist[i-2],genonum[i-2][5]);
							selfmat[i-2][j][temp]=1.0/36.0;
						}
					}
				}
			}
			
			//Coefficients for Biallele, simplex
			tmp=genonum[i-2][0];
			for(j=tmp;j<tmp+genonum[i-2][1];j++)
			{
				selfmat[i-2][j][j]=0.5;
				num=genolist[i-2][j];
				digit1=num%i;
				num/=i;
				digit2=num%i;
				if(digit2<digit1) temp=find_id((digit2*(i*i*i+i*i)+digit1*(i+1)), genolist[i-2],genonum[i-2][5]);
				else{
					if(digit2>digit1) temp=find_id((digit1*(i*i*i+i*i)+digit2*(i+1)), genolist[i-2],genonum[i-2][5]);
				}
				selfmat[i-2][j][temp]=2.0/9.0;
				for(k=0;k<i;k++)
				{
					if(k!=digit2)
					{
						if(k>digit1)
						{
							temp=find_id((digit2*(i*i*i+i*i)+digit1*i+k), genolist[i-2],genonum[i-2][5]);
							selfmat[i-2][j][temp]=1.0/9.0;
						}
						else{
							if(k<digit1)
							{
								temp=find_id((digit2*(i*i*i+i*i)+digit1+k*i), genolist[i-2],genonum[i-2][5]);
								selfmat[i-2][j][temp]=1.0/9.0;
							}
						}
					}
				}
			}
			
			//Coefficients for Biallele, duplex
			tmp+=genonum[i-2][1];
			for(j=tmp;j<tmp+genonum[i-2][2];j++)
			{
				selfmat[i-2][j][j]=0.5;
				num=genolist[i-2][j];
				digit1=num%i;
				num/=i;
				num/=i;
				digit2=num%i;
				for(k=0;k<i;k++)
				{
					if(k!=digit2)
					{
						if(k>digit1)
						{
							temp=find_id((digit2*(i*i*i+i*i)+digit1*i+k), genolist[i-2],genonum[i-2][5]);
							selfmat[i-2][j][temp]=1.0/9.0;
						}
						else{
							if(k<digit1)
							{
								temp=find_id((digit2*(i*i*i+i*i)+digit1+k*i), genolist[i-2],genonum[i-2][5]);
								selfmat[i-2][j][temp]=1.0/9.0;
							}
						}
					}
					if(k!=digit1)
					{
						if(k>digit2)
						{
							temp=find_id((digit1*(i*i*i+i*i)+digit2*i+k), genolist[i-2],genonum[i-2][5]);
							selfmat[i-2][j][temp]=1.0/9.0;
						}
						else{
							if(k<digit2)
							{
								temp=find_id((digit1*(i*i*i+i*i)+digit2+k*i), genolist[i-2],genonum[i-2][5]);
								selfmat[i-2][j][temp]=1.0/9.0;
							}
						}
					}
					if(k!=digit1&&k!=digit2)
					{
						if(digit1<digit2)
						{
							temp=find_id((k*(i*i*i+i*i)+digit1*i+digit2), genolist[i-2],genonum[i-2][5]);
							selfmat[i-2][j][temp]=1.0/36.0;
						}
						else{
							if(digit1>digit2)
							{
								temp=find_id((k*(i*i*i+i*i)+digit2*i+digit1), genolist[i-2],genonum[i-2][5]);
								selfmat[i-2][j][temp]=1.0/36.0;
							}
						}
					}
				}
				for(k=digit1+1;k<i-1;k++)
				{
					for(m=k+1;m<i;m++)
					{
						temp=find_id((digit2*i*i*i+i*i*digit1+k*i+m), genolist[i-2],genonum[i-2][5]);
						selfmat[i-2][j][temp]=1.0/36.0;
					}
				}
			}
			
			//Coefficients for Triallele
			tmp+=genonum[i-2][2];
			for(j=tmp;j<tmp+genonum[i-2][3];j++)
			{
				selfmat[i-2][j][j]=10.0/36.0;
				num=genolist[i-2][j];
				digit1=num%i;
				num/=i;
				digit2=num%i;
				num/=i;
				digit3=num%i;
				
				if(digit3<digit2)	//as we know digit3<digit2<digit1
				{
					temp=find_id((digit2*(i*i*i+i*i)+digit3*i+digit1), genolist[i-2],genonum[i-2][5]);
					selfmat[i-2][j][temp]=1.0/9.0;
					temp=find_id((digit1*(i*i*i+i*i)+digit3*i+digit2), genolist[i-2],genonum[i-2][5]);
					selfmat[i-2][j][temp]=1.0/9.0;
					for(k=0;k<i;k++)
					{
						if(k!=digit1&&k!=digit2&&k!=digit3)
						{
							temp=find_id(calc_val(digit3,digit2,digit1,k,i), genolist[i-2],genonum[i-2][5]);
							selfmat[i-2][j][temp]=1.0/18.0;
						}
					}
				}
				if(digit3>digit2&&digit3<digit1) //digit2<digit3<digit1
				{
					temp=find_id((digit2*(i*i*i+i*i)+digit3*i+digit1), genolist[i-2],genonum[i-2][5]);
					selfmat[i-2][j][temp]=1.0/9.0;
					temp=find_id((digit1*(i*i*i+i*i)+digit2*i+digit3), genolist[i-2],genonum[i-2][5]);
					selfmat[i-2][j][temp]=1.0/9.0;
					for(k=0;k<i;k++)
					{
						if(k!=digit1&&k!=digit2&&k!=digit3)
						{
							temp=find_id(calc_val(digit2,digit3,digit1,k,i), genolist[i-2],genonum[i-2][5]);
							selfmat[i-2][j][temp]=1.0/18.0;
						}
					}
				}
				if(digit3>digit1) //digit2<digit1<digit3
				{
					temp=find_id((digit2*(i*i*i+i*i)+digit1*i+digit3), genolist[i-2],genonum[i-2][5]);
					selfmat[i-2][j][temp]=1.0/9.0;
					temp=find_id((digit1*(i*i*i+i*i)+digit2*i+digit3), genolist[i-2],genonum[i-2][5]);
					selfmat[i-2][j][temp]=1.0/9.0;
					for(k=0;k<i;k++)
					{
						if(k!=digit1&&k!=digit2&&k!=digit3)
						{
							temp=find_id(calc_val(digit2,digit1,digit3,k,i), genolist[i-2],genonum[i-2][5]);
							selfmat[i-2][j][temp]=1.0/18.0;
						}
					}
				}
			}
			
			//Coefficients for Quadriallele
			tmp+=genonum[i-2][3];
			for(j=tmp;j<tmp+genonum[i-2][4];j++)
			{
				selfmat[i-2][j][j]=1.0/6.0;
			}
		
		
		printf("\n\nPrint the selfing matrix with allele_num=%d\n",i);
		for(k=0;k<genonum[i-2][5];k++)
		{
			printf("%d ",genolist[i-2][k]);
		}
		
		printf("\n");
		for(j=0;j<genonum[i-2][5];j++)
		{
			
			for(k=0;k<genonum[i-2][5];k++)
			{
				printf("%.2f ",selfmat[i-2][j][k]);
			}
			printf("\n");
		}
	}
		
	return(selfmat);
	
	
	
	//Coefficients for Biallele, iikk
			for(j=0;j<genonum[i-2][0];j++)
			{
				selfmat[i-2][j][j]=1;	
				num=genolist[i-2][j];
				digit1=num%i;
				num/=i;
				digit2=num%i;
				for(k=0;k<i;k++)
				{
					if(k<digit1)
					{
						temp=find_id((digit2*(i*i*i+i*i)+i*k+digit1), genolist[i-2],genonum[i-2][4]);
						selfmat[i-2][j][temp]=1.0/4.0;
					}
					if(k>digit1)
					{
						temp=find_id((digit2*(i*i*i+i*i)+i*digit1+k), genolist[i-2],genonum[i-2][4]);
						selfmat[i-2][j][temp]=1.0/4.0;
					}
				}
				for(k=0;k<i;k++)
				{
					if(k<digit2)
					{
						temp=find_id((k*i*i*i+i*i*digit2+(i+1)*digit1), genolist[i-2],genonum[i-2][4]);
						selfmat[i-2][j][temp]=1.0/4.0;
					}
					if(k>digit2)
					{
						temp=find_id((digit2*i*i*i+i*i*k+(i+1)*digit1), genolist[i-2],genonum[i-2][4]);
						selfmat[i-2][j][temp]=1.0/4.0;
					}
				}
				for(k=0;k<i;k++)	
				{
					for(m=0;m<i;m++)
					{
						if(k<digit2&&m<digit1)
						{
							temp=find_id((k*i*i*i+i*i*digit2+m*i+digit1), genolist[i-2],genonum[i-2][4]);
							selfmat[i-2][j][temp]=1.0/16.0;
						}
						if(k<digit2&&m>digit1)
						{
							temp=find_id((k*i*i*i+i*i*digit2+i*digit1+m), genolist[i-2],genonum[i-2][4]);
							selfmat[i-2][j][temp]=1.0/16.0;
						}
						if(k>digit2&&m<digit1)
						{
							temp=find_id((digit2*i*i*i+i*i*k+m*i+digit1), genolist[i-2],genonum[i-2][4]);
							selfmat[i-2][j][temp]=1.0/16.0;
						}
						if(k>digit2&&m>digit1)
						{
							temp=find_id((digit2*i*i*i+i*i*k+i*digit1+m), genolist[i-2],genonum[i-2][4]);
							selfmat[i-2][j][temp]=1.0/16.0;
						}
					}
				}
			}
			
			//Coefficients for Triallele1, iikl
			tmp=genonum[i-2][0];
			for(j=tmp;j<tmp+genonum[i-2][1];j++)
			{
				selfmat[i-2][j][j]=0.5;
				num=genolist[i-2][j];
				digit1=num%i;
				num/=i;
				digit2=num%i;
				num/=i;
				digit3=num%i;
				for(k=0;k<i;k++)
				{
					if(k>digit3)
					{
						temp=find_id((digit3*i*i*i+i*i*k+digit2*i+digit1), genolist[i-2],genonum[i-2][4]);
						selfmat[i-2][j][temp]=1.0/8.0;
					}
					if(k<digit3)
					{
						temp=find_id((k*i*i*i+i*i*digit3+digit2*i+digit1), genolist[i-2],genonum[i-2][4]);
						selfmat[i-2][j][temp]=1.0/8.0;
					}
				}
			}
			
			//Coefficients for Triallele2, ijkk
			tmp=genonum[i-2][0];
			for(j=tmp;j<tmp+genonum[i-2][1];j++)
			{
				selfmat[i-2][j][j]=0.5;
				num=genolist[i-2][j];
				digit1=num%i;
				for(k=0;k<i;k++)
				{
					if(k>digit1)
					{
						temp=find_id((num-digit1+k), genolist[i-2],genonum[i-2][4]);
						selfmat[i-2][j][temp]=1.0/8.0;
					}
					if(k<digit1)
					{
						temp=find_id((num+(k-digit1)*i), genolist[i-2],genonum[i-2][4]);
						selfmat[i-2][j][temp]=1.0/8.0;
					}
				}
			}
			
			
			
			//Coefficients for Quadriallele,ijlk
			tmp+=genonum[i-2][2];
			for(j=tmp;j<tmp+genonum[i-2][3];j++)
			{
				selfmat[i-2][j][j]=1.0/4.0;
			}
		
		
		printf("\n\nPrint the selfing matrix with allele_num=%d\n",i);
		for(k=0;k<genonum[i-2][4];k++)
		{
			printf("%d ",genolist[i-2][k]);
		}
		
		printf("\n");
		for(j=0;j<genonum[i-2][4];j++)
		{
			
			for(k=0;k<genonum[i-2][4];k++)
			{
				printf("%.2f ",selfmat[i-2][j][k]);
			}
			printf("\n");
		}
	}
		
	return(selfmat);*/
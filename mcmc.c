/*
 *  mcmc.c
 *  
 *
 *  Created by Hong Gao on 7/21/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
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
#include "DPMM.h"
#include "result_analysis.h"
#include "poly_geno.h"


#define PI 3.141592654


static CHAIN mcmc_POP_no_admixture(SEQDATA data, INIT initial,int chn,CONVG *);
static CHAIN mcmc_POP_admixture(SEQDATA data, INIT initial,int chn,CONVG *);
static CHAIN mcmc_POP_selfing(SEQDATA data, INIT initial,int, CONVG *);
static CHAIN mcmc_INDV_selfing(SEQDATA data, INIT initial,int chn, CONVG *);
static CHAIN mcmc_POP_inbreedcoff(SEQDATA data, INIT initial,int chn, CONVG *);
static CHAIN mcmc_INDV_inbreedcoff(SEQDATA data, INIT initial,int, CONVG *);
static void initial_chn(double ***qqnum,SEQDATA data,UPMCMC **ptr,CHAIN *chain,INIT initial,int chn);
static void free_space(int cnt_step,CHAIN *mchain, int chn, double **qqnum,UPMCMC *ptr,SEQDATA data);
static void update_P(UPMCMC **ptr,SEQDATA data);
static void update_S_IND(int totalsize,UPMCMC **ptr);
static void update_S_POP(SEQDATA data,UPMCMC **ptr);
static void update_inbreedcoff_POP(SEQDATA data,UPMCMC **ptr);
static void update_G(SEQDATA data,UPMCMC **ptr);
static void update_ZQ(UPMCMC **ptr,SEQDATA data,int init_flag,double ***qqnum);
static void update_alpha(UPMCMC **ptr,SEQDATA data,double **qqnum);
static double q(int a,int b);
static double proposal(double *,int* ,double **,int,int );
static double genofreq(int *,double *,int,int);
static double log_ld_indv(int gen,UPMCMC *ptr,int index,SEQDATA data);
static double log_ld_F_pop(double *inbreed,UPMCMC *ptr,int i,SEQDATA data);
static double log_ld_F_indv(double inbreed,UPMCMC *ptr,int i,SEQDATA data);
static double log_ld_F_total(double *inbreed,UPMCMC *ptr,SEQDATA data);
static double log_ld_noselfing_indv(int i,UPMCMC *ptr,SEQDATA data);	
static void cal_lkh(UPMCMC **ptr,SEQDATA data);
static void update_F_IND(int totalsize,UPMCMC **ptr,SEQDATA);
static void update_Z(UPMCMC **ptr,SEQDATA data,int init_flag);
static double log_ld_indv_K(UPMCMC *ptr,SEQDATA data,int i,int K);
static double genofreq_z(SEQDATA data,UPMCMC *ptr, int *z,int i,int j);
static void initialize_chn(CHAIN *chain, SEQDATA data);





CHAIN mcmc_updating(SEQDATA data, INIT initial,int chn, CONVG *cvg)
/*
 * Using different MCMC updating functions according to Mode
 */
{
	CHAIN chain;
	
	if(data.ploid==2)
	{
		switch(data.mode)
		{
			case 0: chain=mcmc_POP_no_admixture(data,initial,chn,cvg);break;
			case 1: chain=mcmc_POP_admixture(data,initial,chn,cvg);break;
			case 2: chain=mcmc_POP_selfing(data,initial,chn,cvg);break;
			case 3: chain=mcmc_INDV_selfing(data,initial,chn,cvg);break;
			case 4: chain=mcmc_POP_inbreedcoff(data,initial,chn,cvg);break;
			case 5: chain=mcmc_INDV_inbreedcoff(data,initial,chn,cvg);break;
		}
	}
	else{
		if(data.ploid==4)
		{	chain=mcmc_POP_tetra_selfing(data,initial,chn,cvg);}
	}
	return(chain);
}


CHAIN mcmc_POP_no_admixture(SEQDATA data, INIT initial,int chn, CONVG *cvg)
/*
 * Function mcmc_updating_POP implements MCMC updating for admixture models and estimates the selfing rates of each subpop.
 */
{
	int j;
	long cnt_step=0,step;
	CHAIN mchain;	
	UPMCMC *ptr;	
	
	mchain.name_len=initial.name_len[chn];
	mchain.chn_name=cvector(0,mchain.name_len-1);
	for(j=0;j<mchain.name_len;j++)
		mchain.chn_name[j]=initial.chn_name[chn][j];
	mchain.steps=(int)((initial.update-initial.burnin)/initial.thinning);
	fprintf(stdout,"\n\n%s Starts:\n",mchain.chn_name);    
	allocate_node(&ptr,data);
	
	/* MCMC loops begins*/     
	update_Z(&ptr,data,1);	//use dirichlet distribution to sample admixture proportions qq		/*initially assign each allele copy of each individual*/ 			
			
	for(step=0;step<initial.update;step++)/*loop for each  MCMC chain*/
	{	      
		update_P(&ptr,data);				//updating P
		update_Z(&ptr,data,0);	/* Z (individual assignments) for each individual*/ 				            
		cal_lkh(&ptr,data);
		if(data.print_iter==1) 
		{	print_info(ptr,data,step,initial.update);}
		if(step==initial.burnin-1)
		{	allocate_chn(&mchain,data);  }
		if(step>=initial.burnin&&(step+1-initial.burnin)%initial.thinning==0)		/*put the node content into the chain at certain thinning interval,o.w. discard it*/
		{										
			store_chn(&mchain,ptr,data);
			if(cnt_step<cvg->ckrep)
			{	cvg->convg_ld[chn*cvg->ckrep+cnt_step]=ptr->totallkh;}
			cnt_step++;			 
		}		
	}
	
	free_node(ptr,data);
	return(mchain);
}



CHAIN mcmc_POP_admixture(SEQDATA data, INIT initial,int chn, CONVG *cvg)
/*
 * Function mcmc_updating_POP implements MCMC updating for admixture models and estimates the selfing rates of each subpop.
 */
{
	long cnt_step=0,step;
	double **qqnum;	
	CHAIN mchain;	
	UPMCMC *ptr;	
	
	initial_chn(&qqnum,data,&ptr,&mchain,initial,chn);
	
	/* MCMC loops begins*/     
	update_ZQ(&ptr,data,1,&qqnum);	//use dirichlet distribution to sample admixture proportions qq		/*initially assign each allele copy of each individual*/ 			
			
	for(step=0;step<initial.update;step++)/*loop for each  MCMC chain*/
	{	      
		update_P(&ptr,data);				//updating P
		update_ZQ(&ptr,data,0,&qqnum);	/*Q (genome percentage) and Z (allele assignments) for each individual*/ 				            
		update_alpha(&ptr,data,qqnum);
		cal_lkh(&ptr,data);
		if(data.print_iter==1) 
		{	print_info(ptr,data,step,initial.update);}
		if(step==initial.burnin-1)
		{	allocate_chn(&mchain,data);  }
		if(step>=initial.burnin&&(step+1-initial.burnin)%initial.thinning==0)		/*put the node content into the chain at certain thinning interval,o.w. discard it*/
		{										
			store_chn(&mchain,ptr,data);
			if(cnt_step<cvg->ckrep)
			{	cvg->convg_ld[chn*cvg->ckrep+cnt_step]=ptr->totallkh;}
			cnt_step++;			 
		}
		if(cnt_step==data.nstep_check_empty_cluster)
		{
			if((mchain.flag_empty_cluster=check_empty_cluster(ptr,data))==1)
			{					
				fprintf(stdout,"Chain %d has an empty cluster, thus discarded!\n",chn+1);
				break;
			}
		}
	}
	
	free_space(cnt_step,&mchain,chn,qqnum,ptr,data);
	return(mchain);
}


CHAIN mcmc_POP_selfing(SEQDATA data, INIT initial,int chn, CONVG *cvg)
/*
 * Function mcmc_updating_POP implements MCMC updating for admixture models and estimates the selfing rates of each subpop.
 */
{
	int i;
	long cnt_step=0,step;
	double **qqnum;	
	CHAIN mchain;	
	UPMCMC *ptr;	
			          	
	initial_chn(&qqnum,data,&ptr,&mchain,initial,chn);	 
	
	/* MCMC begins*/     
	for(i=0;i<data.totalsize;i++) {		/*initialize assign Generations to each individual*/
		ptr->generation[i]=rgeom(ran1());	
		if(ptr->generation[i]>50) ptr->generation[i]=50;
	}
	for(i=0;i<data.popnum;i++)			/*initialize selfing rates*/
	{
		ptr->self_rates[i]=initial.initd[chn][i];
		if(data.back_refl==0) 
		{	ptr->state[i]=dt_stat(ptr->self_rates[i]);}	//determine whether ptr->self_rates[i] is 0.0 (0)state,(0,1) (1)state or 1.0 (2)state
	}
	update_ZQ(&ptr,data,1,&qqnum);	//use dirichlet distribution to sample admixture proportions qq		/*initially assign each allele copy of each individual*/ 			

	for(step=0;step<initial.update;step++)/*loop for MCMC chain*/
	{	      
		update_P(&ptr,data);				//updating P
		update_S_POP(data,&ptr);			/*updating selfing rates for each individual and updating G (generations) for each individual*/
		update_G(data,&ptr);	
		update_ZQ(&ptr,data,0,&qqnum);	/*Q (genome percentage) and Z (allele copy origins) for each individual*/ 				            
		update_alpha(&ptr,data,qqnum);
		cal_lkh(&ptr,data);
		if(data.print_iter==1) 
		{	print_info(ptr,data,step,initial.update);}
		if(step==initial.burnin-1)
		{	allocate_chn(&mchain,data);  }
		if(step>=initial.burnin&&(step+1-initial.burnin)%initial.thinning==0)		/*put the node content into the chain at certain thinning interval,o.w. discard it*/
		{	
			store_chn(&mchain,ptr,data);
			if(cnt_step<cvg->ckrep)
			{	cvg->convg_ld[chn*cvg->ckrep+cnt_step]=ptr->totallkh;}
			cnt_step++;			 
		}
		if(cnt_step==data.nstep_check_empty_cluster)
		{
			if((mchain.flag_empty_cluster=check_empty_cluster(ptr,data))==1)
			{					
				fprintf(stdout,"Chain %d has an empty cluster, thus discarded!\n",chn+1);
				break;
			}
		}
	}
	
	free_space(cnt_step,&mchain,chn,qqnum,ptr,data);
	return(mchain);
}


CHAIN mcmc_POP_inbreedcoff(SEQDATA data, INIT initial,int chn, CONVG *cvg)
/*
 * Function mcmc_updating_POP implements MCMC updating for admixture models and estimates the selfing rates of each subpop.
 */
{
	int i;
	long cnt_step=0,step;
	double **qqnum;	
	CHAIN mchain;	
	UPMCMC *ptr;	
	
	initial_chn(&qqnum,data,&ptr,&mchain,initial,chn);
		
	/* MCMC begins*/     
	for(i=0;i<data.popnum;i++)			/*the initial inbreeding coefficients are*/
	{
		ptr->inbreed[i]=initial.initd[chn][i];
		if(data.back_refl==0) ptr->state[i]=dt_stat(ptr->inbreed[i]);	//determine whether ptr->inbreed[i] is 0.0 (0)state,(0,1) (1)state or 1.0 (2)state
	}
	update_ZQ(&ptr,data,1,&qqnum);	//use dirichlet distribution to sample admixture proportions qq		/*initially assign each allele copy of each individual*/ 			
		
	for(step=0;step<initial.update;step++)/*loop for each  MCMC chain*/
	{	      
		update_P(&ptr,data);				//updating P
		update_inbreedcoff_POP(data,&ptr);			/*updating selfing rates for each individual and updating G (generations) for each individual*/
		update_ZQ(&ptr,data,0,&qqnum);		/*Q (genome percentage) and Z (allele copy origins) for each individual*/ 				            
		update_alpha(&ptr,data,qqnum);
		cal_lkh(&ptr,data);
		if(data.print_iter==1) 
		{	print_info(ptr,data,step,initial.update);}
		if(step==initial.burnin-1)
		{	allocate_chn(&mchain,data);  }
		if(step>=initial.burnin&&(step+1-initial.burnin)%initial.thinning==0)		/*put the node content into the chain at certain thinning interval,o.w. discard it*/
		{										
			store_chn(&mchain,ptr,data);
			if(cnt_step<cvg->ckrep)
			{	cvg->convg_ld[chn*cvg->ckrep+cnt_step]=ptr->totallkh;}
			cnt_step++;			 
		}
		if(cnt_step==data.nstep_check_empty_cluster)
		{
			if((mchain.flag_empty_cluster=check_empty_cluster(ptr,data))==1)
			{					
				fprintf(stdout,"Chain %d has an empty cluster, thus discarded!\n",chn+1);
				break;
			}
		}
	}
		
	free_space(cnt_step,&mchain,chn,qqnum,ptr,data);
	return(mchain);
}



CHAIN mcmc_INDV_selfing(SEQDATA data, INIT initial,int chn, CONVG *cvg)
/*
 * Function mcmc_updating_INDV implements MCMC updating for admixture models and 
 * estimates the selfing rates of each individual
 */
{
	int i,j,cnt_node=0;
	long cnt_step=0,step;
	double **qqnum;		
	CHAIN mchain;	
	UPMCMC *ptr;	
	NODE *head=NULL,*qtr1,*qtr2;
	SF *indv_array=NULL;
	if(data.prior_flag==1) 
	{	
		if((indv_array=(SF *)malloc(data.totalsize*sizeof(SF)))==NULL)
		{	nrerror("Allocation failure in indv_array!");}
	}
	initial_chn(&qqnum,data,&ptr,&mchain,initial,chn);
	
	/* MCMCbegins*/ 				
	if(data.prior_flag==1) 
	{
		head=NULL;
		init_DP(&head,data.alpha_dpm,&indv_array,&cnt_node,data.totalsize);
		for(i=0;i<data.totalsize;i++)
			ptr->self_rates[i]=indv_array[i].value;//write the generated selfing rate into "ptr"
	}
	else{
		for(i=0;i<data.totalsize;i++)		/*initial assign Generations to each individual according to Geometric(1-selfing)*/
			ptr->self_rates[i]=ran1();		/*for each individual S is randomly chosen,thus Generation is randomly chosen*/
	}
	for(i=0;i<data.totalsize;i++)
		ptr->generation[i]=rgeom(1-ptr->self_rates[i]);
	if(ptr->generation[i]>50)	ptr->generation[i]=50;
	update_ZQ(&ptr,data,1,&qqnum);	//use dirichlet distribution to sample admixture proportions qq		/*initially assign each allele copy of each individual*/ 			
	
	for(step=0;step<initial.update;step++)	/*loop for MCMC chain*/
	{
		update_P(&ptr,data);
		if(data.prior_flag==1)
		{	
			update_DP(&head,data.alpha_dpm,&indv_array,&cnt_node,data.totalsize,data,ptr);
			for(j=0;j<data.totalsize;j++)
				ptr->self_rates[j]=indv_array[j].value;
		}/*updating selfing rates for each individual*/	
		if(data.prior_flag==0)
		{	update_S_IND(data.totalsize,&ptr);}
		update_G(data,&ptr);	/*updating G (generations) for each individual*/
		update_ZQ(&ptr,data,0,&qqnum);	/*Q (genome percentage) and Z (allele copy origins) for each individual*/ 				            
		update_alpha(&ptr,data,qqnum);																		
		cal_lkh(&ptr,data);		 	
		if(data.print_iter==1) 
		{	print_info(ptr,data,step,initial.update);}
		if(step==initial.burnin-1)
		{	allocate_chn(&mchain,data);  }
		if(step>=initial.burnin&&(step+1-initial.burnin)%initial.thinning==0)		/*put the node content into the chain at certain thinning interval,o.w. discard it*/
		{					
			store_chn(&mchain,ptr,data);
			if(cnt_step<cvg->ckrep)
			{	cvg->convg_ld[chn*cvg->ckrep+cnt_step]=ptr->totallkh;}
			cnt_step++;			 
		}
		if(cnt_step==data.nstep_check_empty_cluster)
		{
			if((mchain.flag_empty_cluster=check_empty_cluster(ptr,data))==1)
			{					
				fprintf(stdout,"Chain %d has an empty cluster, thus discarded!\n",chn+1);
				break;
			}
		}
	}
		
	if(data.prior_flag==1) 
	{
		qtr1=qtr2=head;
		while(qtr1!=NULL)
		{
			qtr2=qtr1;
			qtr1=qtr1->next;
			free(qtr2);
		}
		free(indv_array);
	}	
	free_space(cnt_step,&mchain,chn,qqnum,ptr,data);	
	return(mchain);
}


CHAIN mcmc_INDV_inbreedcoff(SEQDATA data, INIT initial,int chn, CONVG *cvg)
/*
 * Function mcmc_updating_INDV implements MCMC updating for admixture models and 
 * estimates the selfing rates of each individual
 */
{
	int i,j,cnt_node=0;
	long cnt_step=0,step;
	double **qqnum;		
	CHAIN mchain;	
	UPMCMC *ptr;	
	NODE *head=NULL,*qtr1,*qtr2;
	SF *indv_array=NULL;
	if(data.prior_flag==1) 
	{	
		if((indv_array=(SF *)malloc(data.totalsize*sizeof(SF)))==NULL)
		{	nrerror("Allocation failure in indv_array!");}
	}
	initial_chn(&qqnum,data,&ptr,&mchain,initial,chn);
	
	/* MCMCbegins*/ 				
	if(data.prior_flag==1) 
	{
		head=NULL;
		init_DP(&head,data.alpha_dpm,&indv_array,&cnt_node,data.totalsize);
		for(i=0;i<data.totalsize;i++)
			ptr->inbreed[i]=indv_array[i].value;
	}
	else{
		for(i=0;i<data.totalsize;i++)	/*initial assign Generations to each individual according to Geometric(1-selfing)*/
			ptr->inbreed[i]=ran1();		/*for each individual S is randomly chosen,thus Generation is randomly chosen*/
	}												
	update_ZQ(&ptr,data,1,&qqnum);	//use dirichlet distribution to sample admixture proportions qq		/*initially assign each allele copy of each individual*/ 			
	
	for(step=0;step<initial.update;step++)	/*loop for MCMC chain*/
	{
		update_P(&ptr,data);
		if(data.prior_flag==1)
		{	
			update_DP(&head,data.alpha_dpm,&indv_array,&cnt_node,data.totalsize,data,ptr);
			for(j=0;j<data.totalsize;j++)
				ptr->inbreed[j]=indv_array[j].value;
		}/*updating selfing rates for each individual*/	
		if(data.prior_flag==0)
		{	update_F_IND(data.totalsize,&ptr,data);}
			//put the updated values of individuals into MCMC scheme		
		update_ZQ(&ptr,data,0,&qqnum);	/*Q (genome percentage) and Z (allele copy origins) for each individual*/ 				            
		update_alpha(&ptr,data,qqnum);																		
		cal_lkh(&ptr,data);		 	
		if(data.print_iter==1) 
		{	print_info(ptr,data,step,initial.update);}
		if(step==initial.burnin-1)
		{	allocate_chn(&mchain,data);  }
		if(step>=initial.burnin&&(step+1-initial.burnin)%initial.thinning==0)		/*put the node content into the chain at certain thinning interval,o.w. discard it*/
		{					
			store_chn(&mchain,ptr,data);
			if(cnt_step<cvg->ckrep)
			{	cvg->convg_ld[chn*cvg->ckrep+cnt_step]=ptr->totallkh;}
			cnt_step++;			 
		}
		if(cnt_step==data.nstep_check_empty_cluster)
		{
			if((mchain.flag_empty_cluster=check_empty_cluster(ptr,data))==1)
			{					
				fprintf(stdout,"Chain %d has an empty cluster, thus discarded!\n",chn+1);
				break;
			}
		}
	}		
	if(data.prior_flag==1) 
	{
		qtr1=qtr2=head;
		while(qtr1!=NULL)
		{
			qtr2=qtr1;
			qtr1=qtr1->next;
			free(qtr2);
		}
		free(indv_array);
	}	
	free_space(cnt_step,&mchain,chn,qqnum,ptr,data);	
	return(mchain);
}


void initial_chn(double ***qqnum,SEQDATA data,UPMCMC **ptr,CHAIN *chain,INIT initial,int chn)
{
	int j;
	allocate_node(ptr,data);       /*allocate space for the updated information of each step*/

	if(data.ploid==2&&data.mode!=0) 
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


void free_space(int cnt_step,CHAIN *mchain, int chn, double **qqnum,UPMCMC *ptr, SEQDATA data)
{
	if(mchain->flag_empty_cluster==0)
	{
		if(cnt_step!=mchain->steps)
		{
			nrerror("The number of iterations attained is not the same as counted");
		}		
		fprintf(stdout,"\n\nChain %d is finished running.\n",chn+1);
	}
	
	if(data.ploid==2&&data.mode!=0) free_dmatrix(qqnum,0,data.totalsize-1,0,data.popnum-1);
	free_node(ptr,data);
}


void allocate_node(UPMCMC **ptr,SEQDATA data)
/*
 * Allocate the space for UPMCMC structure
 */
{
	if(((*ptr)=(UPMCMC *)malloc(sizeof(UPMCMC)))==NULL)
	{	nrerror("Allocation failure in ptr");}
	
	if(data.ploid==2)
	{
		switch(data.mode)
		{
			case 2: (*ptr)->self_rates=dvector(0,data.popnum-1);break;
			case 3: (*ptr)->self_rates=dvector(0,data.totalsize-1);break;
			case 4: (*ptr)->inbreed=dvector(0,data.popnum-1);break;
			case 5: (*ptr)->inbreed=dvector(0,data.totalsize-1);break;
		}
	}
	else{
		if(data.ploid==4)	(*ptr)->self_rates=dvector(0,data.popnum-1);
	}
	if(data.back_refl==0&&((data.ploid==2&&(data.mode==2||data.mode==4))||data.ploid==4))
	{	(*ptr)->state=ivector(0,data.popnum-1);}	
	if(data.ploid==2&&(data.mode==2||data.mode==3))
	{	(*ptr)->generation=ivector(0,data.totalsize-1);}
		
	(*ptr)->indvlkh=dvector(0,data.totalsize-1);	
	(*ptr)->freq=d3tensor(0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);			
	if((data.ploid==2&&data.mode!=0)||data.ploid==4)
	{
		(*ptr)->z=i3tensor(0,data.totalsize-1,0,data.locinum-1,0,data.ploid-1);
		(*ptr)->qq=dmatrix(0,data.totalsize-1,0,data.popnum-1);	//store the percentage of each pop. for each individual
	}
	else{	if(data.ploid==2&&data.mode==0)	(*ptr)->zz=ivector(0,data.totalsize-1);}
	if(data.ploid==4)
	{
		if(data.autopoly==0) 
		{	(*ptr)->freq2=d3tensor(0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);}
		(*ptr)->geno=i3tensor(0,data.totalsize-1,0,data.locinum-1,0,data.ploid-1);
	}
}


void free_node(UPMCMC *ptr,SEQDATA data)
{
	if(data.ploid==2)
	{
		switch(data.mode)
		{	
			case 2:free_dvector(ptr->self_rates,0,data.popnum-1);break;
			case 3:free_dvector(ptr->self_rates,0,data.totalsize-1);break;
			case 4:free_dvector(ptr->inbreed,0,data.popnum-1);break;
			case 5:free_dvector(ptr->inbreed,0,data.totalsize-1);break;
		}
	}
	else{
		if(data.ploid==4)	free_dvector(ptr->self_rates,0,data.popnum-1);
	}				
	if(data.back_refl==0&&((data.ploid==2&&(data.mode==2||data.mode==4))||data.ploid==4))
	{	free_ivector(ptr->state,0,data.popnum-1);}
	
	if(data.ploid==2&&(data.mode==2||data.mode==3))
	{	free_ivector(ptr->generation,0,data.totalsize-1);}
	
	free_dvector(ptr->indvlkh,0,data.totalsize-1);
	free_d3tensor(ptr->freq,0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);
	if((data.ploid==2&&data.mode!=0)||data.ploid==4)
	{
		free_i3tensor(ptr->z,0,data.totalsize-1,0,data.locinum-1,0,data.ploid-1);
		free_dmatrix(ptr->qq,0,data.totalsize-1,0,data.popnum-1);
	}
	else{	if(data.ploid==2&&data.mode==0) free_ivector(ptr->zz,0,data.totalsize-1);}
	if(data.ploid==4)
	{
		if(data.autopoly==0) 
		{	free_d3tensor(ptr->freq2,0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);}
		free_i3tensor(ptr->geno,0,data.totalsize-1,0,data.locinum-1,0,data.ploid-1);
	}
	free(ptr);
}		
		
		
void allocate_chn(CHAIN *chain, SEQDATA data)
/*
 * Allocate the space for CHAIN structure
 */
{
	chain->indvlkh=dvector(0,data.totalsize-1);

	if((data.ploid==2&&data.mode!=0)||data.ploid==4)
	{	
		chain->qq=dmatrix(0,data.totalsize-1,0,data.popnum-1);
		chain->qq2=dmatrix(0,data.totalsize-1,0,data.popnum-1);
	}
	else{	
		if(data.ploid==2&&data.mode==0)	
		{
			chain->z=lmatrix(0,data.totalsize-1,0,data.popnum-1);
		}
	}
	if(data.ploid==2)
	{
		switch(data.mode)
		{
			case 2: chain->self_rates=dvector(0,data.popnum-1);
					chain->self_rates2=dvector(0,data.popnum-1);
					break;
			case 3: chain->self_rates=dvector(0,data.totalsize-1);
					chain->self_rates2=dvector(0,data.totalsize-1);
					break;
			case 4: chain->inbreed=dvector(0,data.popnum-1);
					chain->inbreed2=dvector(0,data.popnum-1);
					break;
			case 5: chain->inbreed=dvector(0,data.totalsize-1);
					chain->inbreed2=dvector(0,data.totalsize-1);
					break;
		}
	}
	else{
		if(data.ploid==4)	
		{
			chain->self_rates=dvector(0,data.popnum-1);
			chain->self_rates2=dvector(0,data.popnum-1);
		}
	}
	if(data.ploid==2&&(data.mode==2||data.mode==3))
	{	
		chain->gen=dvector(0,data.totalsize-1);
		chain->gen2=dvector(0,data.totalsize-1);
	}
	if(data.print_freq==1&&data.ploid==2)
	{	
		chain->freq=d3tensor(0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);
		chain->freq2=d3tensor(0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);
	}	
	initialize_chn(chain, data);
}

void initialize_chn(CHAIN *chain, SEQDATA data) //let all elements equal zero
{
	int i,j,k;
	chain->step=0;
	chain->totallkh=1;
	chain->totallkh2=1;
	for(i=0;i<data.totalsize;i++)
	{
		chain->indvlkh[i]=1;
	}
	if((data.ploid==2&&data.mode!=0)||data.ploid==4)
	{	
		for(i=0;i<data.totalsize;i++)
		{
			for(j=0;j<data.popnum;j++)
			{
				chain->qq[i][j]=1;
				chain->qq2[i][j]=1;
			}
		}
	}
	else{	
		if(data.ploid==2&&data.mode==0)	
		{
			for(i=0;i<data.totalsize;i++)
			{
				for(j=0;j<data.popnum;j++)
					chain->z[i][j]=0;
			}
		}
	}
	if(data.ploid==2)
	{
		switch(data.mode)
		{
			case 2: for(j=0;j<data.popnum;j++)
					{
						chain->self_rates[j]=1;
						chain->self_rates2[j]=1;
					}
					break;
			case 3: for(j=0;j<data.totalsize;j++)
					{
						chain->self_rates[j]=1;
						chain->self_rates2[j]=1;
					}
					break;
			case 4: for(j=0;j<data.popnum;j++)
					{
						chain->inbreed[j]=1;
						chain->inbreed2[j]=1;
					}
					break;
			case 5: for(j=0;j<data.totalsize;j++)
					{
						chain->inbreed[j]=1;
						chain->inbreed2[j]=1;
					}
					break;
		}
	}
	else{
		if(data.ploid==4)	
		{
			for(j=0;j<data.popnum;j++)
			{
				chain->self_rates[j]=1;
				chain->self_rates2[j]=1;
			}
		}
	}
	if(data.ploid==2&&(data.mode==2||data.mode==3))
	{	
		for(i=0;i<data.totalsize;i++)
		{
			chain->gen[i]=1;
			chain->gen2[i]=1;
		}
	}
	if(data.print_freq==1&&data.ploid==2)
	{	
		for(j=0;j<data.popnum;j++)
		{
			for(i=0;i<data.locinum;i++)
			{
				for(k=0;k<data.allelenum[i];k++)
				{
					chain->freq[j][i][k]=1;	
					chain->freq2[j][i][k]=1;	
				}
			}
		}
	}
	
}

void free_chain(CHAIN *chain, SEQDATA data)
/*
 * free the space for a MCMC chain
 */
{
	free_cvector(chain->chn_name,0,chain->name_len-1);
	
	free_dvector(chain->indvlkh,0,data.totalsize-1);
	if((data.ploid==2&&data.mode!=0)||data.ploid==4)
	{	free_dmatrix(chain->qq,0,data.totalsize-1,0,data.popnum-1);}
	else{	if(data.ploid==2&&data.mode==0)	free_lmatrix(chain->z,0,data.totalsize-1,0,data.popnum-1);}
	if(data.ploid==2)
	{
		switch(data.mode)
		{	
			case 2: free_dvector(chain->self_rates,0,data.popnum-1);break;
			case 3: free_dvector(chain->self_rates,0,data.totalsize-1);break;
			case 4: free_dvector(chain->inbreed,0,data.popnum-1);break;
			case 5: free_dvector(chain->inbreed,0,data.totalsize-1);break;
		}
	}
	else{
		if(data.ploid==4)	free_dvector(chain->self_rates,0,data.popnum-1);
	}
	//if(data.back_refl==0&&((data.ploid==2&&(data.mode==2||data.mode==4))||data.ploid==4))
	//{	free_ivector(chain->node[j].state,0,data.popnum-1);}
	if(data.ploid==2&&(data.mode==2||data.mode==3))
	{	free_dvector(chain->gen,0,data.totalsize-1);}
	if(data.print_freq==1&&data.ploid==2)
	{
		free_d3tensor(chain->freq,0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);
	}
	
	if((data.ploid==2&&data.mode!=0)||data.ploid==4)
	{	free_dmatrix(chain->qq2,0,data.totalsize-1,0,data.popnum-1);}
	if(data.ploid==2)
	{
		switch(data.mode)
		{	
			case 2: free_dvector(chain->self_rates2,0,data.popnum-1);break;
			case 3: free_dvector(chain->self_rates2,0,data.totalsize-1);break;
			case 4: free_dvector(chain->inbreed2,0,data.popnum-1);break;
			case 5: free_dvector(chain->inbreed2,0,data.totalsize-1);break;
		}
	}
	else{
		if(data.ploid==4)	free_dvector(chain->self_rates2,0,data.popnum-1);
	}
	//if(data.back_refl==0&&((data.ploid==2&&(data.mode==2||data.mode==4))||data.ploid==4))
	//{	free_ivector(chain->node[j].state,0,data.popnum-1);}
	if(data.ploid==2&&(data.mode==2||data.mode==3))
	{	free_dvector(chain->gen2,0,data.totalsize-1);}
	if(data.print_freq==1&&data.ploid==2)
	{
		free_d3tensor(chain->freq2,0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);
	}
}


void update_P(UPMCMC **ptr,SEQDATA data)
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
	{	for(n=0;n<data.totalsize;n++)
		{	if(data.missindx[n][m]!=1&&data.allelenum[m]>1)
			{	
				for(l=0;l<data.allelenum[m];l++)
				{
					for(j=0;j<data.ploid;j++)
					{
						for(i=0;i<data.popnum;i++)
						{
							if(data.mode==0)
							{
								if((*ptr)->zz[n]==i)
								{	
									if(data.seqdata[n][m][j]==l) 	
										seqpop[i][m][l]++;
								}
							}
							else{
								if((*ptr)->z[n][m][j]==i)
								{	
									if(data.seqdata[n][m][j]==l) 	
										seqpop[i][m][l]++;
								}
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
			if(data.allelenum[j]>1)
			{
				for(k=0;k<data.allelenum[j];k++)
				{	tmp_allele[k]=(double)seqpop[i][j][k];}
				rdirich(tmp_allele,data.allelenum[j],&((*ptr)->freq[i][j]),lambda);		
			}
		}
	}
	
	free_i3tensor(seqpop,0,data.popnum-1,0,data.locinum-1,0,data.allelenum_max-1);
	free_dvector(tmp_allele,0,data.allelenum_max-1);
}


void update_S_IND(int totalsize,UPMCMC **ptr)
/*
 * Update the selfing rates of individuals using uniform  prior
 */
{
	int j;
	double tmp,delta0=0.05,mhratio=1.0;
	
	for(j=0;j<totalsize;j++)			
	{
		tmp=ran1()*2*delta0-delta0;	/*propose the selfing rates*/
		tmp+=(*ptr)->self_rates[j];	
						
		if(tmp<=0.0)		
		{	tmp=0.0-tmp;}
		if(tmp>=1.0)			/*reflective bound*/
		{	tmp=1.0-(tmp-1);}
		
		mhratio=exp(log(dgeom(tmp,(*ptr)->generation[j]))-log(dgeom((*ptr)->self_rates[j],(*ptr)->generation[j])));
		
		(*ptr)->self_rates[j]=(ran1()<MIN2(1,mhratio))? tmp:(*ptr)->self_rates[j];		
	}
}

void update_F_IND(int totalsize,UPMCMC **ptr,SEQDATA data)
/*
 * Update the selfing rates of individuals using uniform  prior
 */
{
	int j;
	double tmp,delta0=0.05,mhratio=1.0;
	
	for(j=0;j<totalsize;j++)			
	{
		tmp=ran1()*2*delta0-delta0;	/*propose the selfing rates*/
		tmp+=(*ptr)->inbreed[j];	
						
		if(tmp<=0.0)		
		{	tmp=0.0-tmp;}
		if(tmp>=1.0)			/*reflective bound*/
		{	tmp=1.0-(tmp-1);}
		
		mhratio=exp(log_ld_F_indv(tmp,*ptr,j,data)-log_ld_F_indv((*ptr)->inbreed[j],*ptr,j,data));
		
		(*ptr)->inbreed[j]=(ran1()<MIN2(1,mhratio))? tmp:(*ptr)->inbreed[j];		
	}
}


void update_S_POP(SEQDATA data,UPMCMC **ptr)
/*
 * Update the selfing rates of subpopulations and store them in "ptr"
 */
{
	int i,j,*tem_stat=NULL;
	double delta0=0.05,mhratio=1.0,*tmp=NULL;
	tmp=dvector(0,data.popnum-1);
	if(data.back_refl==0)
	{
		tem_stat=ivector(0,data.popnum-1);
	}	
														
	for(j=0;j<data.popnum;j++)			/*updating selfing rates for each population*/
	{	
		if(data.back_refl==1)
		{
			for(i=0;i<j;i++)
			{	
				tmp[i]=(*ptr)->self_rates[i];
			}
			for(i=j+1;i<data.popnum;i++)
			{	
				tmp[i]=(*ptr)->self_rates[i]; 
			}

			tmp[j]=ran1()*2*delta0-delta0;	/*propose the selfing rates*/
			tmp[j]+=(*ptr)->self_rates[j];	
						
			if(tmp[j]<=0.000)		
			{	tmp[j]=0.000-tmp[j];}
			else{if(tmp[j]>=1.000)			/*reflective bound*/
			{	tmp[j]=1.000-(tmp[j]-1.000);}	}
			
		}
		if(data.back_refl==0)
		{
			for(i=0;i<j;i++)
			{	
				tmp[i]=(*ptr)->self_rates[i];
				tem_stat[i]=(*ptr)->state[i];
			}
			for(i=j+1;i<data.popnum;i++)
			{	
				tmp[i]=(*ptr)->self_rates[i]; 
				tem_stat[i]=(*ptr)->state[i];
			}
			tmp[j]=adpt_indp(&(tem_stat[j]),(*ptr)->state[j]);
			//fprintf(stdout,"%d %d %f\n",tem_stat[j],(*ptr)->state[j],tmp[j]);
		}		

		mhratio=exp(proposal(tmp,(*ptr)->generation,(*ptr)->qq,data.totalsize,data.popnum)-proposal((*ptr)->self_rates,(*ptr)->generation,(*ptr)->qq,data.totalsize,data.popnum));
		if(data.back_refl==0)
		{	mhratio*=hastings_stat(tem_stat,(*ptr)->state,data.popnum);	}
		
		if(ran1()<MIN2(1,mhratio))
		{
			(*ptr)->self_rates[j]=tmp[j];
			if(data.back_refl==0)
			{	(*ptr)->state[j]=tem_stat[j];   } 		
		}
		else{
			(*ptr)->self_rates[j]=(*ptr)->self_rates[j];
			if(data.back_refl==0)
			{	(*ptr)->state[j]=(*ptr)->state[j];	}
		}
	}			
	free_dvector(tmp,0,data.popnum-1);	
	if(data.back_refl==0)
	{	free_ivector(tem_stat,0,data.popnum-1);}
}	


void update_inbreedcoff_POP(SEQDATA data,UPMCMC **ptr)
/*
 * Update the inbreeding coefficients of subpopulations and store them in "ptr"
 */
{
	int i,j,*tem_stat=NULL;
	double delta0=0.05,mhratio=1.0,*tmp=NULL;
	tmp=dvector(0,data.popnum-1);
	if(data.back_refl==0)
	{
		tem_stat=ivector(0,data.popnum-1);
	}
														
	for(j=0;j<data.popnum;j++)			/*updating selfing rates for each population*/
	{	
		if(data.back_refl==1)
		{
			for(i=0;i<j;i++)
			{	
				tmp[i]=(*ptr)->inbreed[i];					
			}
			for(i=j+1;i<data.popnum;i++)
			{	
				tmp[i]=(*ptr)->inbreed[i]; 	
			}

			tmp[j]=ran1()*2*delta0-delta0;	/*propose the selfing rates*/
			tmp[j]+=(*ptr)->inbreed[j];	
						
			if(tmp[j]<=0.000)		
			{	tmp[j]=0.000-tmp[j];}
			else{if(tmp[j]>=1.000)			/*reflective bound*/
			{	tmp[j]=1.000-(tmp[j]-1.000);}	}
		}
		if(data.back_refl==0)
		{
			for(i=0;i<j;i++)
			{	
				tmp[i]=(*ptr)->inbreed[i];
				tem_stat[i]=(*ptr)->state[i];
			}
			for(i=j+1;i<data.popnum;i++)
			{	
				tmp[i]=(*ptr)->inbreed[i]; 
				tem_stat[i]=(*ptr)->state[i];
			}
			tmp[j]=adpt_indp(&(tem_stat[j]),(*ptr)->state[j]);
			//fprintf(stdout,"%d %d %f\n",tem_stat[j],(*ptr)->state[j],tmp[j]);
		}		

		mhratio=log_ld_F_total(tmp,*ptr,data)-log_ld_F_total((*ptr)->inbreed,*ptr,data);
		if(data.back_refl==0)
		{	mhratio*=hastings_stat(tem_stat,(*ptr)->state,data.popnum);	}
		
		if(ran1()<exp(MIN2(1,mhratio)))
		{		
			(*ptr)->inbreed[j]=tmp[j];
			if(data.back_refl==0)
			{	(*ptr)->state[j]=tem_stat[j];  }  		
		}		
	}			
	free_dvector(tmp,0,data.popnum-1);	
	if(data.back_refl==0)
	{	free_ivector(tem_stat,0,data.popnum-1);}
}


void update_G(SEQDATA data,UPMCMC **ptr)
/*
 *updating G (generations) for each individual
 */
{	
	int i,j,stat,gen=0;
	double selfing=0,mhratio=1;
	for(i=0;i<data.totalsize;i++)				/*individual-updating loop begins*/
	{	       	    
		if(data.mode==2)
		{
			selfing=0;
			for(j=0;j<data.popnum;j++)				/*update generation for each individual*/
				selfing+=(*ptr)->qq[i][j]*(*ptr)->self_rates[j];								/*grtprpl is the proposed generation*/
		}
		if(data.mode==3)
		{	selfing=(*ptr)->self_rates[i];	}
		
		stat=dt_stat(selfing);
		if(stat==1)
		{	
			gen=rgeom(1-selfing);
			if(gen<1) gen=1;
			if(gen>50) gen=50;    /*have a upbound over the generations until outcross*/
		}
		else{
			if(stat==0)	{	gen=1;}
			else{	
				if(stat==2) { gen=50;}
				else {nrerror("No such state for selfing rate!");}
			}
		}
		mhratio=exp(log_ld_indv(gen,*ptr,i,data)-log_ld_indv((*ptr)->generation[i],*ptr,i,data));
		if(ran1()<MIN2(1,mhratio))
		{
			(*ptr)->generation[i]=gen;
		}
	}
}


void update_Z(UPMCMC **ptr,SEQDATA data,int init_flag)
/*
 * Update Z and Q according to the Metropolis-Hastings ratio
 */
{	
	int i,m;
	double *tmp=NULL,temp=0;
	tmp=dvector(0,data.popnum-1);
	
	for(i=0;i<data.totalsize;i++)
	{
		for(m=0;m<data.popnum;m++)
		{	
			if(init_flag==1)
			{	tmp[m]=(double)(m+1)/data.popnum;}
			else{
				tmp[m]=log_ld_indv_K(*ptr,data,i,m);
				if(m==0)	{temp=tmp[m];}
				tmp[m]=exp(tmp[m]-temp);
				if(m>=1) {tmp[m]+=tmp[m-1];}
			}
		}
		(*ptr)->zz[i]=disc_unif(tmp,data.popnum);
	}
	free_dvector(tmp,0,data.popnum-1);
}


void update_ZQ(UPMCMC **ptr,SEQDATA data,int init_flag,double ***qqnum)
/*
 * Update Z and Q according to the Metropolis-Hastings ratio
 */
{	
	int i,j,k,m,l,temp,*tmp1;
	double *tmp=NULL,mhratio=1,inbrd=0,temp1=0.0,temp2=0.0;
	tmp=dvector(0,data.popnum-1);
	tmp1=ivector(0,data.ploid-1);
	//tmp2=dvector(0,data.ploid-1);
	
	for(i=0;i<data.totalsize;i++)
	{
		for(j=0;j<data.locinum;j++)		/*updating Z*/
		{
			if(data.missindx[i][j]!=1&&data.allelenum[j]>1)
			{
				for(k=0;k<data.ploid;k++)
				{
					for(m=0;m<data.popnum;m++)
					{	
						if(init_flag==1)
						{	tmp[m]=(double)(m+1)/data.popnum;}
						else{
							tmp[m]=(*ptr)->qq[i][m]*(*ptr)->freq[m][j][data.seqdata[i][j][k]];
							if(m>=1) {tmp[m]+=tmp[m-1];}
						}
					}
					temp=disc_unif(tmp,data.popnum);
					//if(init_flag==1||(*ptr)->z[i][j][k]==temp)
					//{	
					(*ptr)->z[i][j][k]=temp;
					/*}
					else{
						//for(m=0;m<data.ploid;m++)
						//	tmp1[m]=(*ptr)->z[i][j][m];
						//tmp1[k]=temp;
						mhratio=exp(genofreq_z(data,*ptr,tmp1,i,j)-genofreq_z(data,*ptr,(*ptr)->z[i][j],i,j))*(*ptr)->qq[i][temp]/(*ptr)->qq[i][(*ptr)->z[i][j][k]];
						//mhratio=(*ptr)->qq[i][temp]*(*ptr)->freq[temp][j][data.seqdata[i][j][k]]/((*ptr)->qq[i][(*ptr)->z[i][j][k]]*(*ptr)->freq[(*ptr)->z[i][j][k]][j][data.seqdata[i][j][k]]);
						if(mhratio<0.05)	
						{
							for(m=0;m<data.ploid;m++)
								printf("%d ",(*ptr)->z[i][j][m]);
							printf("%d tmp=%d %f %f mhratio=%f %f\n",k,tmp1[k],genofreq_z(data,*ptr,tmp1,i,j),genofreq_z(data,*ptr,(*ptr)->z[i][j],i,j),mhratio,(*ptr)->qq[i][temp]/(*ptr)->qq[i][(*ptr)->z[i][j][k]]);
						}
						if(ran1()<MIN2(1,mhratio))
						{	
							(*ptr)->z[i][j][k]=temp;
						}
					}*/
				}				
			}
		}

		for(m=0;m<data.popnum;m++)/*updating Q*/
			(*qqnum)[i][m]=0.0;
		for(l=0;l<data.locinum;l++)		
		{
			if(data.missindx[i][l]!=1&&data.allelenum[l]>1)
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
	//free_dvector(tmp2,0,data.ploid-1);
}


double genofreq_z(SEQDATA data,UPMCMC *ptr, int *z,int i,int j)
{
	int m;
	double temp=0.0,*tmp1,inbrd;
	tmp1=dvector(0,data.ploid-1);
	
	if(chcksame(z,data.ploid)==1)
	{
		temp=0.0;
		for(m=0;m<data.ploid;m++)
		{
			temp+=log(ptr->freq[z[m]][j][data.seqdata[i][j][m]]);
		}
		if(chcksame(data.seqdata[i][j],data.ploid)==1)
		{	temp+=log(2);}
	}
	else{
		for(m=0;m<data.ploid;m++)
		{
			tmp1[m]=ptr->freq[z[m]][j][data.seqdata[i][j][m]];
		}
		switch(data.mode)
		{
			case 1: temp=log(genofreq(data.seqdata[i][j],tmp1,1,data.ploid));break;
			case 2: temp=log(genofreq(data.seqdata[i][j],tmp1,ptr->generation[i],data.ploid));break;
			case 3: temp=log(genofreq(data.seqdata[i][j],tmp1,ptr->generation[i],data.ploid));break; 
			case 4: inbrd=0;
					for(m=0;m<data.popnum;m++)	
					{	inbrd+=ptr->qq[i][m]*ptr->inbreed[m];	}
					temp=log(genofreq_inbreedcoff(data.seqdata[i][j],tmp1,inbrd,data.ploid));break; 
			case 5: temp=log(genofreq_inbreedcoff(data.seqdata[i][j],tmp1,ptr->inbreed[i],data.ploid));break;
		}
	}

	free_dvector(tmp1,0,data.ploid-1);
	return(temp);
}

void update_alpha(UPMCMC **ptr,SEQDATA data,double **qqnum)
/*
 * Update alpha according to the Metropolis-Hastings ratio
 */
{
	double alphasd=1.0,mhratio=1.0,ralpha;					/*updating alphas*/
	int i,m;
	ralpha=rnormal((*ptr)->alpha,alphasd);
	if(ralpha>0)
	{
		for(i=0;i<data.totalsize;i++)
		{
			for(m=0;m<data.popnum;m++)	/*metropolis-hasting ratio= */
			{						/*f(a')/f(a) as q(a'|a)/q(a|a')=1*/
				mhratio*=pow((*ptr)->qq[i][m],ralpha+qqnum[i][m])/pow((*ptr)->qq[i][m],qqnum[i][m]+(*ptr)->alpha);		
			}
		}
		(*ptr)->alpha=(ran1()<MIN2(1,mhratio))?ralpha:(*ptr)->alpha;
	}
}



void print_info(UPMCMC *ptr,SEQDATA data,int step,int maxstep)
/*
 * print the information of each updating step
 */ 
{
	int i,s;
	s=maxstep/100;
	if(step%s!=0) return;
	fprintf(stdout,"\nStep=%d\tlog_likelihood=%f\n",step+1,ptr->totallkh);
	if(data.mode==2||data.ploid==4)
	{
		for(i=0;i<data.popnum;i++)
		{	
			fprintf(stdout,"s_%d=%f",i,ptr->self_rates[i]);
			if(data.back_refl==0)
			{	fprintf(stdout," st_%d=%d",i,ptr->state[i]);}
			if(i<data.popnum-1) fprintf(stdout," ");
		}
		fprintf(stdout,"\n");	
	}
	if(data.mode==4)
	{
		for(i=0;i<data.popnum;i++)
		{	
			fprintf(stdout,"f_%d=%f",i,ptr->inbreed[i]);
			if(data.back_refl==0)
			{	fprintf(stdout," st_%d=%d",i,ptr->state[i]);}
			if(i<data.popnum-1) fprintf(stdout," ");
		}
		fprintf(stdout,"\n");	
	}	 
	if(data.mode==3)
	{
		for(i=0;i<data.totalsize;i++)
			fprintf(stdout,"s_%d=%f ",i,ptr->self_rates[i]);
		fprintf(stdout,"\n");	
	}
	if(data.mode==5)
	{
		for(i=0;i<data.totalsize;i++)
			fprintf(stdout,"f_%d=%f ",i,ptr->inbreed[i]);
		fprintf(stdout,"\n");	
	}
	/*if(data.mode==2||data.mode==3)
	{
		for(i=0;i<data.totalsize;i++) 
			fprintf(stdout,"g_%d=%d ",i,ptr->generation[i]);
		fprintf(stdout,"\n");
	}*/
}



void store_chn(CHAIN *mchain,UPMCMC *ptr,SEQDATA data)
/*
 *store the updated information in CHAINS at each iteration
 */
{
	int i,j,k;
	
	if(mchain->totallkh!=0)
	{	mchain->totallkh=mchain->totallkh*((mchain->step+ptr->totallkh/mchain->totallkh)/(1+mchain->step));}
	else{	mchain->totallkh=ptr->totallkh/(1+mchain->step);}
	if(mchain->totallkh2!=0)
	{	mchain->totallkh2=mchain->totallkh2*((mchain->step+ptr->totallkh*ptr->totallkh/mchain->totallkh2)/(1+mchain->step));}
	else{	mchain->totallkh2=ptr->totallkh*ptr->totallkh/(1+mchain->step);}
	
	for(i=0;i<data.totalsize;i++)
	{
		if(mchain->indvlkh[i]!=0)
		{	mchain->indvlkh[i]=mchain->indvlkh[i]*((mchain->step+ptr->indvlkh[i]/mchain->indvlkh[i])/(1+mchain->step));}
		else{	mchain->indvlkh[i]=ptr->indvlkh[i]/(1+mchain->step);}
	}
	if((data.ploid==2&&data.mode!=0)||data.ploid==4)
	{	
		for(i=0;i<data.totalsize;i++)
		{
			for(j=0;j<data.popnum;j++)
			{
				if(mchain->qq[i][j]!=0)
				{	mchain->qq[i][j]=mchain->qq[i][j]*((mchain->step+ptr->qq[i][j]/mchain->qq[i][j])/(1+mchain->step));}
				else{	mchain->qq[i][j]=ptr->qq[i][j]/(1+mchain->step);	}
				if(mchain->qq2[i][j]!=0)
				{	mchain->qq2[i][j]=mchain->qq2[i][j]*((mchain->step+ptr->qq[i][j]*ptr->qq[i][j]/mchain->qq2[i][j])/(1+mchain->step));}
				else{	mchain->qq2[i][j]=ptr->qq[i][j]*ptr->qq[i][j]/(1+mchain->step);	}
			}
		}
	}
	else{	
		if(data.ploid==2&&data.mode==0)	
		{
			for(i=0;i<data.totalsize;i++)
			{
					mchain->z[i][ptr->zz[i]]+=1;
			}
		}
	}
	if(data.ploid==2)
	{
		switch(data.mode)
		{
			case 2: for(j=0;j<data.popnum;j++)
					{
						if(mchain->self_rates[j]!=0)
						{	mchain->self_rates[j]=mchain->self_rates[j]*((mchain->step+ptr->self_rates[j]/mchain->self_rates[j])/(1+mchain->step));}
						else{	mchain->self_rates[j]=ptr->self_rates[j]/(1+mchain->step);	}
						if(mchain->self_rates2[j]!=0)
						{	mchain->self_rates2[j]=mchain->self_rates2[j]*((mchain->step+ptr->self_rates[j]*ptr->self_rates[j]/mchain->self_rates2[j])/(1+mchain->step));}
						else{	mchain->self_rates2[j]=ptr->self_rates[j]*ptr->self_rates[j]/(1+mchain->step);	}
					}
					break;
			case 3: for(j=0;j<data.totalsize;j++)
					{
						if(mchain->self_rates[j]!=0)
						{	mchain->self_rates[j]=mchain->self_rates[j]*((mchain->step+ptr->self_rates[j]/mchain->self_rates[j])/(1+mchain->step));}
						else{	mchain->self_rates[j]=ptr->self_rates[j]/(1+mchain->step);	}
						if(mchain->self_rates2[j]!=0)
						{	mchain->self_rates2[j]=mchain->self_rates2[j]*((mchain->step+ptr->self_rates[j]*ptr->self_rates[j]/mchain->self_rates2[j])/(1+mchain->step));}
						else{	mchain->self_rates2[j]=ptr->self_rates[j]*ptr->self_rates[j]/(1+mchain->step);	}
					}
					break;
			case 4: for(j=0;j<data.popnum;j++)
					{
						if(mchain->inbreed[j]!=0)
						{	mchain->inbreed[j]=mchain->inbreed[j]*((mchain->step+ptr->inbreed[j]/mchain->inbreed[j])/(1+mchain->step));}
						else{	mchain->inbreed[j]=ptr->self_rates[j]/(1+mchain->step);	}
						if(mchain->inbreed2[j]!=0)
						{	mchain->inbreed2[j]=mchain->inbreed2[j]*((mchain->step+ptr->inbreed[j]*ptr->inbreed[j]/mchain->inbreed2[j])/(1+mchain->step));}
						else{	mchain->inbreed2[j]=ptr->self_rates[j]*ptr->self_rates[j]/(1+mchain->step);	}
					}
					break;
			case 5: for(j=0;j<data.totalsize;j++)
					{
						if(mchain->inbreed[j]!=0)
						{	mchain->inbreed[j]=mchain->inbreed[j]*((mchain->step+ptr->inbreed[j]/mchain->inbreed[j])/(1+mchain->step));}
						else{	mchain->inbreed[j]=ptr->self_rates[j]/(1+mchain->step);	}
						if(mchain->inbreed2[j]!=0)
						{	mchain->inbreed2[j]=mchain->inbreed2[j]*((mchain->step+ptr->inbreed[j]*ptr->inbreed[j]/mchain->inbreed2[j])/(1+mchain->step));}
						else{	mchain->inbreed2[j]=ptr->self_rates[j]*ptr->self_rates[j]/(1+mchain->step);	}
					}
					break;
		}
	}
	else{
		if(data.ploid==4)	
		{
			for(j=0;j<data.popnum;j++)
			{
				if(mchain->self_rates[j]!=0)
				{	mchain->self_rates[j]=mchain->self_rates[j]*((mchain->step+ptr->self_rates[j]/mchain->self_rates[j])/(1+mchain->step));}
				else{	mchain->self_rates[j]=ptr->self_rates[j]/(1+mchain->step);	}
				if(mchain->self_rates2[j]!=0)
				{	mchain->self_rates2[j]=mchain->self_rates2[j]*((mchain->step+ptr->self_rates[j]*ptr->self_rates[j]/mchain->self_rates2[j])/(1+mchain->step));}
				else{	mchain->self_rates2[j]=ptr->self_rates[j]*ptr->self_rates[j]/(1+mchain->step);	}
			}
		}
	}
	if(data.ploid==2&&(data.mode==2||data.mode==3))
	{	
		for(i=0;i<data.totalsize;i++)
		{
			if(mchain->gen[i]!=0)
			{	mchain->gen[i]=mchain->gen[i]*((mchain->step+ptr->generation[i]/mchain->gen[i])/(1+mchain->step));}
			else{	mchain->gen[i]=ptr->generation[i]/(1+mchain->step);	}
			if(mchain->gen2[i]!=0)
			{	mchain->gen2[i]=mchain->gen2[i]*((mchain->step+ptr->generation[i]*ptr->generation[i]/mchain->gen2[i])/(1+mchain->step));}
			else{	mchain->gen2[i]=ptr->generation[i]*ptr->generation[i]/(1+mchain->step);	}
		}
	}
	if(data.print_freq==1&&data.ploid==2)
	{	
		for(j=0;j<data.popnum;j++)
		{
			for(i=0;i<data.locinum;i++)
			{
				for(k=0;k<data.allelenum[i];k++)
				{
					if(mchain->freq[j][i][k]!=0)
					{	mchain->freq[j][i][k]=mchain->freq[j][i][k]*((mchain->step+ptr->freq[j][i][k]/mchain->freq[j][i][k])/(1+mchain->step));	}
					else{	mchain->freq[j][i][k]=ptr->freq[j][i][k]/(1+mchain->step);}
					if(mchain->freq2[j][i][k]!=0)
					{	mchain->freq2[j][i][k]=mchain->freq2[j][i][k]*((mchain->step+ptr->freq[j][i][k]*ptr->freq[j][i][k]/mchain->freq2[j][i][k])/(1+mchain->step));}
					else{	mchain->freq2[j][i][k]=ptr->freq[j][i][k]*ptr->freq[j][i][k]/(1+mchain->step);}	
				}
			}
		}
	}
	
	mchain->step++;
}




double adpt_indp(int *stat_tmp,int stat)
/*
 * Transition probability matrix for Adaptive Independence Sampler for selfing rates or inbreeding coefficients
 */
{
	double tmp=0,tt=0;
	if(stat==0)
	{	
		if(ran1()<0.50)	
		{	
			tmp=0.000; 
			*stat_tmp=0;
		}
		else{	
			tmp=ran1();
			*stat_tmp=1;
		}
	}
	else{
		if(stat==2)
		{	
			if(ran1()<0.5)	
			{
				tmp=1.000;
				*stat_tmp=2;
			}
			else{	
				tmp=ran1();
				*stat_tmp=1;
			}
		}
		else{
			if(stat==1)
			{	
				tt=ran1();
				if(tt<=0.05)	
				{
					tmp=0.0000; 
					*stat_tmp=0;
				}
				else{
					if(tt>=0.95)
					{
						tmp=1.000;
						*stat_tmp=2;
					}
					else{	
						tmp=ran1();
						*stat_tmp=1;
					}
				}
			}
			else{
				fprintf(stdout,"ERROR: State %d is beyond 0, 1 and 2!\n",stat);//there can only three states, 0,1 and 2.
				exit(1);
			}
		}
	}
	return(tmp);
}



int dt_stat(double num)
/*
 * Determine which state a selfing rate falls in, state 0={0}, state 1=(0,1), and state 2={1}.
 */
{
	int stat=-1;
	double epislon=0.001;
	if(num<=0.000+epislon&&num>=0.000-epislon)
	{	stat=0;}
	else{
		if(num>=1.000-epislon&&num<=1.000+epislon)
		{	stat=2;}
		else{
			if(num>=0.0+epislon&&num<1.000-epislon)
			{	stat=1;}
			else{
				fprintf(stdout,"ERROR: The value of selfing rate or inbreeding coefficient %f is beyond [0,1]!\n",num);
				exit(1);
			}
		}
	}
	return(stat);
}



double hastings_stat(int *tmp,int *prev,int num)
/*
 * Calculate the Hastings ratio for the Adaptive Independence Sampler
 * treat selfing rates of subpopulations independent
 */
 {
	int i;
	double temp=1.0;
	for(i=0;i<num;i++)
	{	
		temp*=q(prev[i],tmp[i])/q(tmp[i],prev[i]);
	}
	return(temp);
}


double q(int a,int b)
/*
 * Calculate the numerate or denominator of the Hastings ratio
 * in function hastings_stat()
 */
{
	double temp=0;
	if(a==0)
	{
		if(b==0)	{	temp=0.5;}
		if(b==1)	{	temp=0.5;}
	}
	else{
		if(a==2)
		{
			if(b==2) {	temp=0.5;}
			if(b==1) {	temp=0.5;}
		}
		else{
			if(a==1)
			{
				if(b==0||b==2)	{	temp=0.05;}
				if(b==1)		{	temp=0.90;}
			}
		}
	}
	return(temp);
}


double dgeom(double self,int gen)
/*
 *calculate the density value of geometric distribution
 */
{
	double tmp;
	tmp=pow(self,(double)(gen-1))*(1-self);
	return(tmp);
}


void sample_mu2(double *mudraw, double *sigma2draw, double ave,double pre_mu,double *data, int n, double kappa_0, double nu_0, double mu_0, double sigmasqr_0)
/*
 * draw mu and var from a hierarchical normal prior, mu ~ Normal and var ~ Inverse Gamma
 */
{
  int i;
  double sd, random, kappa_n, mu_n, nu_n, sigmasqr_n,sum=0;
  kappa_n = kappa_0 + n;
  nu_n = nu_0 + n;
  for(i=0;i<n;i++)
  {
	sum+=pow((pre_mu-data[i]),2);
  }
  sigmasqr_n = nu_0*sigmasqr_0 + kappa_0*pow((pre_mu-mu_0),2)+sum;
  random = sigmasqr_n/rgamma(nu_n*0.5, 0.5);
  *sigma2draw = random;
  mu_n = (kappa_0*mu_0 + n*ave)/kappa_n;
  sd = sqrt(random/kappa_n);
  *mudraw = rnormal(mu_n, sd);
}



double proposal(double *inbreed,int* gen,double** qq,int totalsize, int popnum)
/*
 * calculates the log-likelihood of generations of individual given population selfing rates
 */
{
	double ld=0,temp;
	int i,j;
		
	for(i=0;i<totalsize;i++)
	{
		temp=0;
		for(j=0;j<popnum;j++)
		{
			temp+=qq[i][j]*inbreed[j];	//take the expected selfing rate of each individual
		}
		ld+=log(pow(temp,gen[i]-1)*(1-temp));
	}
	return(ld);
}



/*
 * Function chcksame checks whether all the elements in the integer 
 * vector "pop" of length "num" are the same or not.
 * It returns a boolean integer. 0 means all the elements are the same,
 * 1 means some might be different from others
 */
int chcksame(int *pop,int num)
{
	int flag=0;
	int i;
	for(i=1;i<num;i++)
	{
		if(pop[i]!=pop[0]) flag=1;
	}
	return(flag);
}



/*
 * Function genofreq calculates genotype frequency given selfing 
 * generations and current allele frequencies
 * Input argument: 
 *		locus_ata contains the genotype, thus for diploid, 
 *				it contains two elements			
 *		freq contains the allele frequencies for that locus
 *		generation is the number of generation that the individual has selfed
 *		ploid is the number of haplotype in one genome, for diploid, it is 2
 * This function returns a double number that is the frequency of  
 *		the genotype given in "locus_data"
 */
double genofreq(int *locus_data,double *freq,int generation,int ploid)
{
	double result,temp=0;
	int i;
	
	if(chcksame(locus_data,ploid)==0)
	{
		result=pow(freq[0],(double)ploid);
		temp=2*freq[0]*(1-freq[0]);			//homozygote 
		for(i=1;i<generation;i++)
		{
			temp/=2;
			result+=temp/2;
		}					
	}
	else							//heterozygote
	{
		result=2*freq[0]*freq[1]*pow(0.5,(double)(generation-1));			
	}	
	return(result);		
}



double genofreq_inbreedcoff(int *seqdata,double *freq,double inbreed,int ploid)
/*
 * Calculates the genotype frequencies given inbreeding coefficients
 */
{
	double result;

	if(chcksame(seqdata,ploid)==0)		//homozygote
	{
		result=pow(freq[0],(double)ploid)*(1-inbreed)+freq[0]*inbreed;			 
	}
	else								//heterozygote
	{
		result=2*freq[0]*freq[1]*(1-inbreed);	
	}	
	return(result);		
}


double log_ld_indv(int gen,UPMCMC *ptr,int index,SEQDATA data)
/*
 * Calculate the log likelihood of a given generation for individual "index" for mode 1 and mode 2
 */
{
	int j,k,m;
	double temp=0,*tmp;					/*calculate the probability */
	tmp=dvector(0,data.ploid-1);
	
	for(j=0;j<data.locinum;j++)		/*for the proposed generations*/
	{
		if(data.missindx[index][j]!=1&&data.allelenum[j]>1)
		{				
			if(data.type_freq==0)
			{
				for(k=0;k<data.ploid;k++)
				{
					tmp[k]=0;
					for(m=0;m<data.popnum;m++)
						tmp[k]+=ptr->freq[m][j][data.seqdata[index][j][k]]*ptr->qq[index][m];
				}
				temp+=log(genofreq(data.seqdata[index][j],tmp,gen,data.ploid));
				
			}
			if(data.type_freq==1)
			{
				if(chcksame(ptr->z[index][j],data.ploid)==0)	//same assignment
				{
					for(m=0;m<data.ploid;m++)
					{
						tmp[m]=ptr->freq[ptr->z[index][j][m]][j][data.seqdata[index][j][m]];
					}
					temp+=log(genofreq(data.seqdata[index][j],tmp,gen,data.ploid));
				}
				else{
					for(m=0;m<data.ploid;m++)
					{
						temp+=log(ptr->freq[ptr->z[index][j][m]][j][data.seqdata[index][j][m]]);
					}
					if(chcksame(data.seqdata[index][j],data.ploid)==1) //heterozygote needs to time 2
					{	temp+=log(2);}
				}
			}
		}                                     
	}
	free_dvector(tmp,0,data.ploid-1);
	return(temp);
}


double log_ld_F_pop(double *inbreed,UPMCMC *ptr,int i,SEQDATA data)
/*
 * Calculate the log likelihood of an individual "i" given the inbreeding coefficients for all subpop.
 */
{
	int j,m;
	double temp=0,*tmp,ld=0;					/*calculate the probability */
	tmp=dvector(0,data.ploid-1);
	
	for(j=0;j<data.locinum;j++)		
	{
		if(data.missindx[i][j]!=1&&data.allelenum[j]>1)
		{	
			if(chcksame(ptr->z[i][j],data.ploid)==0)
			{
				for(m=0;m<data.ploid;m++)
				{
					tmp[m]=ptr->freq[ptr->z[i][j][m]][j][data.seqdata[i][j][m]];
				}
				temp=genofreq_inbreedcoff(data.seqdata[i][j],tmp,inbreed[ptr->z[i][j][0]],data.ploid);
				ld+=log(temp);
			}
			else{
				for(m=0;m<data.ploid;m++)
				{
					ld+=log(ptr->freq[ptr->z[i][j][m]][j][data.seqdata[i][j][m]]);
				}
				if(chcksame(data.seqdata[i][j],data.ploid)==1) //heterozygote needs to time 2
				{	ld+=log(2);}
			}
		}                                     
	}
	free_dvector(tmp,0,data.ploid-1);
	return(ld);
}


double log_ld_F_indv(double inbreed,UPMCMC *ptr,int i,SEQDATA data)
/*
 * Calculate the log likelihood of an individual "i" given the inbreeding coefficients for that individual
 */
{
	int j,m;
	double temp=0,*tmp,ld=0;					/*calculate the probability */
	tmp=dvector(0,data.ploid-1);
	
	for(j=0;j<data.locinum;j++)		
	{
		if(data.missindx[i][j]!=1&&data.allelenum[j]>1)
		{	
			if(chcksame(ptr->z[i][j],data.ploid)==0)
			{
				for(m=0;m<data.ploid;m++)
				{
					tmp[m]=ptr->freq[ptr->z[i][j][m]][j][data.seqdata[i][j][m]];
				}
				temp=genofreq_inbreedcoff(data.seqdata[i][j],tmp,inbreed,data.ploid);
				ld+=log(temp);
			}
			else{
				for(m=0;m<data.ploid;m++)
				{
					ld+=log(ptr->freq[ptr->z[i][j][m]][j][data.seqdata[i][j][m]]);
				}
				if(chcksame(data.seqdata[i][j],data.ploid)==1) //heterozygote needs to time 2
				{	ld+=log(2);}
			}
		}                                     
	}
	free_dvector(tmp,0,data.ploid-1);
	return(ld);
}



double log_ld_F_total(double *inbreed,UPMCMC *ptr,SEQDATA data)
{
	int i;
	double ld=0;
	
	switch(data.mode)
	{
		case 4: for(i=0;i<data.totalsize;i++)
					ld+=log_ld_F_pop(inbreed,ptr,i,data);
				break;
		case 5: for(i=0;i<data.totalsize;i++)
					ld+=log_ld_F_indv(inbreed[i],ptr,i,data);
				break;
	}
	return(ld);
}


double log_ld_noselfing_indv(int i,UPMCMC *ptr,SEQDATA data)
/*
 * Calculate the log likelihood of an individual
 */
{
	int j,m;
	double ld=0;					
		
	for(j=0;j<data.locinum;j++)		
	{
		if(data.missindx[i][j]!=1&&data.allelenum[j]>1)
		{	
			for(m=0;m<data.ploid;m++)
			{
				ld+=log(ptr->freq[ptr->z[i][j][m]][j][data.seqdata[i][j][m]]);
			}
			if(chcksame(data.seqdata[i][j],data.ploid)==1) //heterozygote needs to time 2
			{	ld+=log(2);}
		}                                     
	}
	return(ld);
}


double log_ld_indv_K(UPMCMC *ptr,SEQDATA data,int i,int K)
/*
 * Calculate the log likelihood of an individual
 */
{
	int j,m;
	double ld=0;					
		
	for(j=0;j<data.locinum;j++)		
	{
		if(data.missindx[i][j]!=1&&data.allelenum[j]>1)
		{	
			for(m=0;m<data.ploid;m++)
			{
				ld+=log(ptr->freq[K][j][data.seqdata[i][j][m]]);
			}
			if(chcksame(data.seqdata[i][j],data.ploid)==1) //heterozygote needs to time 2
			{	ld+=log(2);}
		}                                     
	}
	return(ld);
}

void cal_lkh(UPMCMC **ptr,SEQDATA data)
/*
 * calculate the log-likelihood of individuals and the total log likelihood
 */
{
	int i;
	(*ptr)->totallkh=0;
	for(i=0;i<data.totalsize;i++)			
	{								/*calculate the log likelihood for each individual*/
		switch(data.mode)
		{
			case 0: (*ptr)->indvlkh[i]=log_ld_indv_K((*ptr),data,i,(*ptr)->zz[i]);
					break;
			case 1: (*ptr)->indvlkh[i]=log_ld_noselfing_indv(i,(*ptr),data);
					break;
			case 2: (*ptr)->indvlkh[i]=log_ld_indv((*ptr)->generation[i],*ptr,i,data);
					break;
			case 3: (*ptr)->indvlkh[i]=log_ld_indv((*ptr)->generation[i],*ptr,i,data);
					break;
			case 4: (*ptr)->indvlkh[i]=log_ld_F_pop((*ptr)->inbreed,(*ptr),i,data);
					break;
			case 5: (*ptr)->indvlkh[i]=log_ld_F_indv((*ptr)->inbreed[i],(*ptr),i,data);
					break;
		}
		(*ptr)->totallkh+=(*ptr)->indvlkh[i];
	} 
}

int check_empty_cluster(UPMCMC *ptr,SEQDATA data)
/*
 *check whether there exists any empty cluster, which has no more than one individual. 
 *If empty cluster exists, return 1; otherwise, return 0.
 */
{
	int flag=0,j,k;
	double *sum;
	sum=dvector(0,data.popnum-1);
	
	for(k=0;k<data.popnum;k++)
	{
		sum[k]=0;
		for(j=0;j<data.totalsize;j++)
		{
			sum[k]+=ptr->qq[j][k];
		}		
	}

	for(k=0;k<data.popnum;k++)
	{
		//fprintf(stdout,"pop%d=%f\n",k+1,sum[k]);
		if(sum[k]<0.01)
		{
			flag=1;
			break;
		}
	}
	free_dvector(sum,0,data.popnum-1);
	return(flag);
}


/*
 *  result_analysis.c
 *  
 *
 *  Created by Hong Gao on 7/25/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "nrutil.h"
#include "data_interface.h"
#include "initial.h"
#include "check_converg.h"
#include "mcmc.h"
#include "quantile.h"
#include "result_analysis.h"


static double print_lkh_to_file(FILE* outfile,CHAIN chain,SEQDATA data);
static void print_S_POP_to_file(FILE* outfile,CHAIN chain,SEQDATA data,int **indx);
static void print_S_INDV_to_file(FILE* outfile,CHAIN chain,SEQDATA data);
static void print_F_POP_to_file(FILE* outfile,CHAIN chain,SEQDATA data,int **indx);
static void print_F_INDV_to_file(FILE* outfile,CHAIN chain,SEQDATA data);
static void print_gen_to_file(FILE* outfile,CHAIN chain,SEQDATA data);
//static void fill_in_convg(CHAIN chain,SEQDATA data, CONVG*,int *indx,int chn);
static void print_Q_to_file(FILE* outfile,CHAIN chain,SEQDATA data,int *indx);
static void print_P_to_file(FILE* outfile,CHAIN chain,SEQDATA data,int *indx);
static void print_Z_to_file(FILE* outfile,CHAIN chain,SEQDATA data);


double chain_stat(char *outfilename,CHAIN chain,SEQDATA data,int chn)
/*
 * Function chain_stat summarizes all the results from MCMC chains 
 * and print the final results to output file
 */
{
	int *indx;	
	double dic;														
	FILE *outfile;
	indx=ivector(0,data.popnum+1);
	
	//print file
	if((outfile=fopen(outfilename,"a+"))==NULL)
	{	nrerror("Cannot open output file!");}
	dic=print_lkh_to_file(outfile,chain,data);		
	if(data.ploid==2)
	{	
		switch(data.mode)
		{
			case 2: print_S_POP_to_file(outfile,chain,data,&indx);print_gen_to_file(outfile,chain,data);break;
			case 3: print_S_INDV_to_file(outfile,chain,data);print_gen_to_file(outfile,chain,data);break;
			case 4: print_F_POP_to_file(outfile,chain,data,&indx);break;
			case 5: print_F_INDV_to_file(outfile,chain,data);break;
		}
	}	
	else{	if(data.ploid==4)	{	print_S_POP_to_file(outfile,chain,data,&indx);}}
	if((data.ploid==2&&data.mode!=0)||data.ploid==4)	{	print_Q_to_file(outfile,chain,data,indx);}
	else{	if(data.ploid==2&&data.mode==0)	print_Z_to_file(outfile,chain,data);}
	if(data.print_freq==1&&data.ploid==2)
	{	
		print_P_to_file(outfile,chain,data,indx);
	}
	fclose(outfile);

	free_ivector(indx,0,data.popnum+1);
	return(dic);
}



void print_S_POP_to_file(FILE* outfile,CHAIN chain,SEQDATA data,int **indx)
{
	int j;
	double *F_res;
												//print the summarized information into output file 
	F_res=dvector(0,data.popnum);
	
	fprintf(outfile,"\nThe Posterior distribution of Selfing Rates:\n");
	fprintf(outfile,"\t\tMean\tVar\n");
	for(j=0;j<data.popnum;j++)
	{
		F_res[j+1]=chain.self_rates[j];
		indexx(data.popnum,F_res,*indx);
	}
	for(j=0;j<data.popnum;j++)
	{
		fprintf(outfile,"Cluster %d\t%.3f\t%.3f\n",j+1,chain.self_rates[(*indx)[j+1]-1],chain.self_rates2[(*indx)[j+1]-1]-chain.self_rates[(*indx)[j+1]-1]*chain.self_rates[(*indx)[j+1]-1]);
	}
	
	free_dvector(F_res,0,data.popnum);
}


void print_S_INDV_to_file(FILE* outfile,CHAIN chain,SEQDATA data)
{	
	int j;
	
	fprintf(outfile,"\nThe Posterior distribution of Selfing Rates:\n");
	if(data.label==1) {fprintf(outfile,"\t");}
	fprintf(outfile,"\t\tMean\tVar\n");
	for(j=1;j<=data.totalsize;j++)
	{
		fprintf(outfile,"Indv %d\t\t",j);
		if(data.label==1)
		{	fprintf(outfile, "%s\t",data.indvname[j]);}
		fprintf(outfile,"%.3f\t%.3f\n",chain.self_rates[j],chain.self_rates2[j]-chain.self_rates[j]*chain.self_rates[j]);
	}
}


void print_F_POP_to_file(FILE* outfile,CHAIN chain,SEQDATA data,int **indx)
{
	int j;
	double *F_res;
	F_res=dvector(0,data.popnum);
	
	fprintf(outfile,"\nThe Posterior distribution of Inbreeding Coefficients:\n");
	fprintf(outfile,"\t\tMean\tVar\n");
	for(j=0;j<data.popnum;j++)
	{
		F_res[j+1]=chain.inbreed[j];
		indexx(data.popnum,F_res,*indx);
	}
	for(j=0;j<data.popnum;j++)
	{
		fprintf(outfile,"Cluster %d\t%.3f\t%.3f\n",j+1,chain.inbreed[(*indx)[j+1]-1],chain.inbreed2[(*indx)[j+1]-1]-chain.inbreed[(*indx)[j+1]-1]*chain.inbreed[(*indx)[j+1]-1]);
	}
	
	free_dvector(F_res,0,data.popnum);
}

void print_F_INDV_to_file(FILE* outfile,CHAIN chain,SEQDATA data)
{	
	int j;
		
	fprintf(outfile,"\nThe Posterior distribution of Inbreeding Coefficients:\n");
	fprintf(outfile,"\t\tMean\tVar\n");
	for(j=0;j<data.totalsize;j++)
	{
		fprintf(outfile,"Indv %d\t\t",j);
		if(data.label==1)
		{	fprintf(outfile, "%s\t",data.indvname[j]);}
		fprintf(outfile,"%.3f\t%.3f\n",chain.inbreed[j],chain.inbreed2[j]-chain.inbreed[j]*chain.inbreed[j]);
	}
}



	
void print_Z_to_file(FILE* outfile,CHAIN chain,SEQDATA data)
{
	int i,j;
	double **tmp;
	tmp=dmatrix(0,data.totalsize-1,0,data.popnum-1);
	
	for(i=0;i<data.totalsize;i++)
	{
		for(j=0;j<data.popnum;j++)
			tmp[i][j]=(double)chain.z[i][j]/(double)chain.steps;		
	}

	fprintf(outfile,"\nInferred Classification of individuals:\n");
	fprintf(outfile,"\nIndv\t");
	if(data.label==1)
	{	fprintf(outfile,"Label\t");	}
	fprintf(outfile,"(Miss)\t");
	if(data.popdata==1)
	{	fprintf(outfile,"Pop : ");	}	
	for(j=0;j<data.popnum;j++)
		fprintf(outfile,"Prob in Cluster %d\t",j+1);
	fprintf(outfile,"\n");
	for(j=0;j<data.totalsize;j++)
	{
		fprintf(outfile,"%d\t",j+1);
		if(data.label==1) 
		{	fprintf(outfile,"%s\t",data.indvname[j]);}
		fprintf(outfile,"(%d)\t",(int)data.missvec[j]);
		if(data.popdata==1) 
		{	fprintf(outfile,"%d : ",data.popindx[j]);}
		for(i=0;i<data.popnum;i++)
		{
			fprintf(outfile,"\t%f",tmp[j][i]);
		}
		fprintf(outfile,"\n");
	}
	
	free_dmatrix(tmp,0,data.totalsize-1,0,data.popnum-1);
}
	
	
void print_Q_to_file(FILE* outfile,CHAIN chain,SEQDATA data,int *indx)
{
	int i,j,k,*tm;
	long count;
	double **tmp;
															
	count=chain.steps;																//print the summarized information into output file 
	if(data.popdata==0)
	{	data.pop_count=1;}	
	tmp=dmatrix(0,data.pop_count-1,0,data.popnum-1);
	tm=ivector(0,data.pop_count-1);
	for(i=0;i<data.pop_count;i++)
	{
		for(j=0;j<data.popnum;j++)
			tmp[i][j]=0;
		tm[i]=0;
	}

	fprintf(outfile,"\nInferred ancestry of individuals:\n");
	fprintf(outfile,"\nIndv\t");
	if(data.label==1)
	{	fprintf(outfile,"Label\t");	}
	fprintf(outfile,"(Miss)\tPop : ");	
	if(data.distr_fmt==0)
	{
		for(j=0;j<data.popnum;j++)
			fprintf(outfile,"Cluster %d:Mean\tVar\t\t",j+1);
	}
	else if(data.distr_fmt==1)
	{
		for(j=0;j<data.popnum;j++)
			fprintf(outfile,"\tCluster %d",j+1);
	}
	fprintf(outfile,"\n");
	for(j=0;j<data.totalsize;j++)
	{
		fprintf(outfile,"%d\t",j+1);
		if(data.label==1) 
		{	fprintf(outfile,"%s\t",data.indvname[j]);}
		fprintf(outfile,"(%d)\t",(int)data.missvec[j]);
		if(data.popdata==1) 
		{	fprintf(outfile,"%d : ",data.popindx[j]);}
		else{	fprintf(outfile,"1 : ");}
		for(k=0;k<data.popnum;k++)
		{
			if(data.popdata==1) tmp[data.popindx[j]][k]+=chain.qq[j][k];
			else{	tmp[0][k]+=chain.qq[j][k];}
		}
		if((data.ploid==2&&(data.mode==2||data.mode==4))||data.ploid==4)
		{	
			for(k=0;k<data.popnum;k++)
			{
				if(data.distr_fmt==0)
				{
					fprintf(outfile,"\t%.3f\t%.3f\t",chain.qq[j][k],chain.qq2[j][k]-chain.qq[j][k]*chain.qq[j][k]);
				}
				else if(data.distr_fmt==1)
				{
					fprintf(outfile,"\t%.3f",chain.qq[j][k]);
				}
			}
		}
		if(data.ploid==2&&(data.mode==1||data.mode==3||data.mode==5))
		{	
			for(k=0;k<data.popnum;k++)
			{
				if(data.distr_fmt==0)
				{
					fprintf(outfile,"\t%.3f\t%.3f\t",chain.qq[j][k],chain.qq2[j][k]-chain.qq[j][k]*chain.qq[j][k]);
				}
				else if(data.distr_fmt==1)
				{
					fprintf(outfile,"\t%.3f",chain.qq[j][k]);
				}
			}
		}
		if(data.popdata==1) tm[data.popindx[j]]++;
		else{	tm[0]++;}
		fprintf(outfile,"\n");
	}

	
	fprintf(outfile,"\n\n\nThe index and name of pre-defined populations:\n");
	if(data.popdata==1)
	{
		for(i=0;i<data.pop_count;i++)
		{
			fprintf(outfile,"%d %s\n",i,data.poptype[i]);
		}
	}
	else{	fprintf(outfile,"1\n");}
	fprintf(outfile,"\n\nProportion of membership of each pre-defined population in each of the %d clusters\n",data.popnum);
	fprintf(outfile,"Given Pop\tInferred Clusters\t\tNumber ofIndividuals\n");
	fprintf(outfile,"    \t\t");
	for(i=0;i<data.popnum;i++)
	{
		fprintf(outfile,"%d    ",i+1);
	}
	fprintf(outfile,"\n");
	for(i=0;i<data.pop_count;i++)
	{
		fprintf(outfile,"%d:\t",i);
		for(j=0;j<data.popnum;j++)
		{
			if((data.ploid==2&&(data.mode==2||data.mode==4))||data.ploid==4) {fprintf(outfile,"%.3f ",tmp[i][indx[j+1]-1]/tm[i]);}
			if(data.ploid==2&&(data.mode==1||data.mode==3||data.mode==5)) {fprintf(outfile,"%.3f ",tmp[i][j]/tm[i]);}
		}
		fprintf(outfile,"\t%d\n",tm[i]);
	}
	fprintf(outfile,"\n");
	
	free_dmatrix(tmp,0,data.pop_count-1,0,data.popnum-1);
	free_ivector(tm,0,data.pop_count-1);
}





void print_P_to_file(FILE* outfile,CHAIN chain,SEQDATA data,int *indx)
/*
 * print the summarized allele frequencies at each locus into output file
 */
{
	int i,j,k;
		
	fprintf(outfile,"\n\n\nEstimated allele frequencies:\n");
	fprintf(outfile,"\nLocus_ID\t");
	if(data.markername_flag==1)
	{
		fprintf(outfile,"Marker Name\t");
	}
	fprintf(outfile,"Alleletype\t");
	for(j=0;j<data.popnum;j++)
		fprintf(outfile,"Cluster %d:Mean\tVar\t\t",j+1);
	fprintf(outfile,"\n");
	for(j=0;j<data.locinum;j++)
	{		
		for(i=0;i<data.allelenum[j];i++)
		{
			if(i==0)
			{
				fprintf(outfile,"%d\t",j+1);
				if(data.markername_flag==1)
				{
					fprintf(outfile,"%s\t",data.marker_names[j]);
				}
			}
			else{	
				fprintf(outfile,"\t");
				if(data.markername_flag==1)
				{
					fprintf(outfile,"\t");
				}
			}
			fprintf(outfile,"%s\t",data.alleletype[j][i]);
			
			if(data.mode==2||data.mode==4)
			{	
				for(k=0;k<data.popnum;k++)
				{
					fprintf(outfile,"\t%.3f\t%.3f\t",chain.freq[indx[k+1]-1][j][i],chain.freq2[indx[k+1]-1][j][i]-chain.freq[indx[k+1]-1][j][i]*chain.freq[indx[k+1]-1][j][i]);
				}
			}
			if(data.mode==1||data.mode==3||data.mode==5||data.mode==0)
			{	
				for(k=0;k<data.popnum;k++)
				{
					fprintf(outfile,"\t%.3f\t%.3f\t",chain.freq[k][j][i],chain.freq2[k][j][i]-chain.freq[k][j][i]*chain.freq[k][j][i]);
				}
			}
			fprintf(outfile,"\n");
		}
		fprintf(outfile,"\n");
	}
}



void print_gen_to_file(FILE* outfile,CHAIN chain,SEQDATA data)
{
	int j;

	fprintf(outfile,"\nThe Posterior distribution of Generations:\n");
	fprintf(outfile,"\t\tMean\tVariance\n");
	for(j=0;j<data.totalsize;j++)
	{
		fprintf(outfile,"Indv %d\t\t",j);
		if(data.label==1)
		{	fprintf(outfile, "%s\t",data.indvname[j]);}
		fprintf(outfile,"%.3f\t%.3f\n",chain.gen[j],chain.gen2[j]-chain.gen[j]*chain.gen[j]);
	}
}


double print_lkh_to_file(FILE* outfile,CHAIN chain,SEQDATA data)
{
	int i,j;
	double tempt,DIC=0;
	
	fprintf(outfile,"\n\n\n");
	for(i=0;i<chain.name_len;i++)
		fprintf(outfile,"%c",chain.chn_name[i]);
	fprintf(outfile,":\n");
	fprintf(outfile,"\nThe log Likelihood:\n");
	fprintf(outfile,"    Posterior Mean = %.3f\n",chain.totallkh);
	fprintf(stdout,"    Posterior Mean = %.3f\n",chain.totallkh);
	fprintf(outfile,"    Posterior Variance = %.3f\n",chain.totallkh2-chain.totallkh*chain.totallkh);
		
	//calculate the DIC
	tempt=0.0;
	for(j=0;j<data.totalsize;j++)
	{
		tempt+=chain.indvlkh[j];	// tempt is log(D|E[theta|D]))
	}
	
	DIC=-4*chain.totallkh+2*tempt;
	fprintf(outfile,"\nThe Deviance information criterion of this model is %f.\n", DIC);
	
	return(DIC);
}




double mean(double *vec,int length)
{
	int i;
	double sum=0;
	for(i=0;i<length;i++)
		sum+=vec[i];

	return(sum/(double)length);
}



double variance(double *vec,int length)
{
	int i;
	double sum=0,mu;
	mu=mean(vec,length);
	for(i=0;i<length;i++)
		sum+=pow((vec[i]-mu),2);

	return(sum/(double)(length-1));
}

int find_min(double *vec,int length)
{
	int i,min=0;
	double temp=vec[0];

	for(i=1;i<length;i++)
	{
		if(temp>vec[i])
		{
			min=i;
			temp=vec[i];
		}
	}
	return(min);
}




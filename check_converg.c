/*
 *  check_converg.c
 *  
 *
 *  Created by Hong Gao on 7/25/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include "nrutil.h"
#include "data_interface.h"
#include "initial.h"
#include "quantile.h"
#include "check_converg.h"


static double GelmanRubin(double *vec, int numchains, int totrep);




void allocate_convg(SEQDATA data, CONVG *cvg,int chainnum,int ckrep,char *convgfilename)
/*
 * allocate space for the structure CONVG which stores updates for assessing convergence
 */
{
	cvg->n_chain=chainnum;
	cvg->ckrep=ckrep;
	cvg->convgfilename=convgfilename;
	cvg->convg_ld=dvector(0,(cvg->n_chain*cvg->ckrep)-1);
}		


void free_convg(CONVG *cvg)
/*
 * free space for the array "convg" which stores updates for assessing convergence
 */
{	
	free_dvector(cvg->convg_ld,0,(cvg->n_chain*cvg->ckrep)-1);
}			        

int chain_converg(char *outfilename,CONVG *cvg)
/*
 * Function chain_converg is used to mainly check the convergence of 
 * selfing rates or inbreeding coefficients using the Gelman_Rubin 
 * statistics and print the mixing condition to the output file.
 */
{
	int i,k,flag=0;
	double GR_ld,threshold=1.1;
	FILE *fp;
	
	if((fp=fopen(outfilename,"a+"))==NULL)
	{
		nrerror("ERROR: Cannot open output file!\n");
	}
	
	if(cvg->n_chain==1)
	{
		fprintf(fp,"There is only one MCMC. No need to check the convergence.\n");
	}
	else{
		
		fprintf(stdout,"The Gelman-Rubin statistics of log-likelihood is ");
		GR_ld=GelmanRubin(cvg->convg_ld,cvg->n_chain,cvg->ckrep);
		fprintf(stdout,"%f\n",GR_ld);
		fprintf(fp,"\n\nThe Gelman-Rubin statistics for the convergence of log-likelihood is %f.\n",GR_ld);	
		if(GR_ld>threshold) flag=1;
		fclose(fp);
	}
	
	//output the array "convg" to file "convgfilename"
	if(cvg->convgfilename!=NULL)
	{
		if((fp=fopen(cvg->convgfilename,"w"))==NULL)
		{
			nrerror("ERROR: Cannot open convergence file!\n");
		}
		fprintf(fp,"Values of log-likelihood:\n");
		for(k=0;k<(cvg->ckrep*cvg->n_chain);k++)
		{
			if(k==0) {fprintf(fp,"%f ",cvg->convg_ld[k]);}
			else{fprintf(fp," %f ",cvg->convg_ld[k]);}
		}
		fprintf(fp,"\n");
		fclose(fp);
	}
	return(flag);
}








double GelmanRubin(double *vec, int numchains, int totrep)
/*
 * Function GelmanRubin is used to calculate the the Gelman_Rubin statistics
 * Based on an ANOVA idea for a single variable 
 * Asumme m chains, each of length n 
 * Can estimate the variance of a stationary distribution in two ways 
 *		Ðvariance within a single chain, W 
 *		Ðvariance over all chains, B/n
 * If the chains have converged, both estimates are unbiased, i.e. B=W 
 * If the initial values are overdispersed and have not dispersed, then the Between term is an overestimate 
 * The statistic: R = B/W 
 * If R>1, it have not converged, we estimate R by (m+1)/m*((n-1)/n+B/W)-(n-1)/mn
 */
{
	double *psii, psi, *S, W, B, V;
	int i, j, repperchain;
  
	psii = dvector(0,numchains-1);

	S = dvector(0,numchains-1);
		
	repperchain = totrep/numchains;
	psi=0;
	for (i=0; i<numchains; i++)
	{
		psii[i]=0;
		for (j=0; j<repperchain; j++)
			psii[i]+=vec[i*repperchain+j];
		psii[i]=psii[i]/repperchain;
		psi=psi+psii[i];
	}
	psi=psi/numchains;
	W = 0;
	for (i=0; i<numchains; i++)
	{
		S[i]=0;
		for (j=0; j<repperchain; j++)
			S[i]+=(vec[i*repperchain+j]-psii[i])*(vec[i*repperchain+j]-psii[i]);
		S[i]=S[i]/(repperchain-1);
		W+=S[i];
	}
	W=W/numchains;
	B=0;
	for (i=0; i<numchains; i++)
		B+=(psii[i]-psi)*(psii[i]-psi);
	
	B=(B*repperchain)/(numchains-1);
	V=(W*(repperchain-1))/repperchain+B/repperchain;
   
  //  printf("B: %f. W: %f (%f %f)\n",B,W,(W*(repperchain-1))/repperchain,B/repperchain);
	free_dvector(psii,0,numchains-1);
	free_dvector(S,0,numchains-1);
	return V/W;
}



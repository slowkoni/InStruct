/*
 *  check_converg.h
 *  
 *
 *  Created by Hong Gao on 7/25/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

typedef struct convg
	{
		double *convg_ld;		//check the convergence of log-likelihood
	//	double **convg_param;	//check the convergence of parameters, such as selfing rates or inbreeding coefficients
	//	int n_param;
		int n_chain;
		int ckrep;
		char *convgfilename;	//the name of the file to store values to be checked of convergence
	}CONVG;
	
int chain_converg(char *outfilename,CONVG *cvg);
void free_convg(CONVG *cvg);
void allocate_convg(SEQDATA data, CONVG *cvg,int chainnum,int ckrep,char *convgfilename);

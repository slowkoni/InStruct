/*
 *  poly_geno.h
 *  
 *
 *  Created by Hong Gao on 3/7/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

typedef struct polyploid
	{
		int num_allele;
		int *allele_poly; //store the distinct number of alleles across loci
		int **genonum;
		int **genolist;
		int num_allogeno[4];
		int num_autogeno[4];
		float ***exfreq;
		float ***genofreq;
	}POLY;

CHAIN mcmc_POP_tetra_selfing(SEQDATA data, INIT initial,int chn,CONVG *);
//void free_tetra_chain(CHAIN *chain, SEQDATA data);


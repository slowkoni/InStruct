/*
 *  data_interface.h
 *  
 *
 *  Created by Hong Gao on 7/21/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

typedef struct seqd
	{
		int ploid;
		int popnum;
		int locinum;
		int totalsize;		
		
		int datafmt;
		int ***seqdata;
		int *allelenum;
		int allelenum_max; //record the maximum number of alleles at any locus
		char ***alleletype;
		int **alleleid;
		
		char **indvname;
		char **poptype;
		int *popindx;
		int pop_count;
		int label;
		int popdata;
		int markername_flag;
		char **marker_names;
		int n_extra_col;
		char ***extra_col;
		
		char *missingdata;
		int missingnum;
		int *missvec;
		int **missindx;
		
		double siglevel;
		int back_refl;
		int type_freq;
		int nstep_check_empty_cluster;
		int prior_flag;
		int mode;
		
		double alpha_dpm;
		int print_iter;
		int print_freq;
		int inf_K;
		int distr_fmt;
		int autopoly;
		double max_mem;	//maximum memory allowed
		
		
	}SEQDATA;

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void sort(int leng, int *vec);
int exists(int value, int *vec, int leng);	
SEQDATA read_data(char *infilename,int ploid,int totalsize,int popnum,int nloci,char *missingdata,int label,int popdata, double siglevel,int back_refl,int type_freq,int nstep_check_empty_cluster,int prior_flag,int mode,int n_extra_col,int markername_flag,double alpha_dpm, int print_iter,int print_freq,int,int,int,int,double);

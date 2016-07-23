 /*
 *  mcmc.h
 *  
 *
 *  Created by Hong Gao on 7/21/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */
 
#define MIN2(X, Y) (((X) > (Y)) ? (Y) : (X)) 

typedef struct UPMC
	{							
		int *generation;	/*generation[individual index]*/
		double ***freq;		/*freq[population index][loci number][allelenum]*/
		double ***freq2;	//for allotetraploid
		int ***z;			/*z[individual index][loci number][chromosome index]*/
		int *zz;
		double **qq;		/*qq[individual index][population index]*/
		double alpha;		/*alpha updating uses Metroplis-Hastings algorithm*/
		double *inbreed;	/*inbreed stores the selfing rate for each population*/
		int *state;
		double *self_rates;	/*self_rates stores the selfing rate for each individual*/
		double totallkh;
		double *indvlkh;
		int ***geno;
	}UPMCMC;

typedef struct MC
	{		/*this structure stores updating information*/
		long steps; //store the total iterations
		long step; //record the current iteration
		int name_len;
		char *chn_name;
		int flag_empty_cluster;
		
		double totallkh;
		double *indvlkh;
		double *self_rates;	//record to calculate the mean
		double **qq;
		double *inbreed;
		//int *state;
		double *gen;			
		long **z;
		double ***freq;		/*freq[population index][loci number][allelenum]*/
		
		double totallkh2;	//record the square of the values in order to calculate the variance
		double *self_rates2;
		double **qq2;
		double *inbreed2;
		double *gen2;
		double ***freq2;	
	}CHAIN;


CHAIN mcmc_updating(SEQDATA data, INIT initial,int chn, CONVG *);		
void free_chain(CHAIN *chain, SEQDATA);
int chcksame(int *,int);
double genofreq_inbreedcoff(int *seqdata,double *tem,double inbreed,int ploid);
double dgeom(double ,int gen);
void print_info(UPMCMC *ptr,SEQDATA data,int step,int);
double adpt_indp(int *stat_tmp,int stat);
double hastings_stat(int *tmp,int *prev,int num);
int dt_stat(double num);
void allocate_node(UPMCMC **ptr,SEQDATA);
void free_node(UPMCMC *ptr,SEQDATA);
void allocate_chn(CHAIN *mchain, SEQDATA data);
void store_chn(CHAIN *mchain,UPMCMC *ptr,SEQDATA data);
int check_empty_cluster(UPMCMC *ptr,SEQDATA data);




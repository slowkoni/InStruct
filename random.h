/*
 *  admrand.h
 *  
 *
 *  Created by Hong Gao on 7/18/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

double wichmann();					/*Random number generator for U(0,1) distribution*/

double ran1();						/*Random number generator for U(0,1) distribution*/

void setseeds(int sd1, int sd2,int sd3);/*set the seeds new values by users*/

void printseeds(FILE* fp);			/*print the seeds into a file "fp"*/

float gammln(float xx);				/* returns ln[Gamma[xx]] for xx > 0*/

float factln(int n)	;				/* returns ln(n!)*/

float bico(int n, int k);			/* returns binomial coefficients (nCk)*/

void moment(double data[], int n, float *ave, float *var); /*Calculate the mean and variance of "data"*/

double rexp(double lambda);			/*Generates from an exponential distribution*/

float poidev(float xm);                  // Returns a poisson RV with mean xm 

double rgamma1(double alpha);		/*Generates from a gamma distribution with alpha < 1*/

double rgamma2(double alpha);		/*Generates from a gamma distribution with alpha > 1*/

double rgamma(double alpha, double beta);/*Generates from a general gamma(alpha,beta) distribution*/

double rbeta(double alpha,double beta);	/* Generates from a beta (alpha,beta) distribution*/

void rdirich(double *alpha,int length,double **rand,double add);/* Generates from a Dirichlet distribution*/

double rstd_normal();				/*Generates from a standard normal(0,1) distribution*/

double rnormal(double mean, double sd);	/*Generates from a general normal(mu,sigma) distribution*/

int rgeom(double p);				/*Generates from a geometric distribution*/

int rbinom(double p, int N);		/* Generates from a binomial (N,p) distribution*/

float bnldev(float pp, int n);

int disc_unif(double *vec,int length);	//generate the integer that represent the interval where a uniform random value falls


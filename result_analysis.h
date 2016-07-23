/*
 *  result_analysis.h
 *  
 *
 *  Created by Hong Gao on 7/25/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

		
double chain_stat(char *outfilename,CHAIN chain,SEQDATA data,int chn);

double variance(double *vec,int length);

double mean(double *vec,int length);

int find_min(double *vec,int length);

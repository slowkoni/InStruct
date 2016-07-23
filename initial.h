/*
 *  initial.h
 *  
 *
 *  Created by Hong Gao on 7/25/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */
typedef struct initialdata
		{
			float **initd;
			char **chn_name;
			int *name_len;
			int chainnum;
			long update;
			long burnin;
			int thinning;
			int popnum;
		}INIT;			//INIT only contains the initial values of inbreeding coefficients for each population
		
INIT read_init(char *initialfilename,int chainnum,int popnum,long update,long burnin,int thinning);

void int_split(char *s,int *nums,int num);
void word_split(char *,float *,int );
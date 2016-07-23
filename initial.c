/*
 *  initial.c
 *  
 *
 *  Created by Hong Gao on 7/25/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "nrutil.h"
#include "initial.h"
#include "random.h"

#define MAXLINE 1000
#define MAXLEN 100
static int word_cnt(char *);
static char* int_to_string(int );

/*
 * Function read_init is used to read in the initial values of selfing rates (S) for MCMC updating
 * Input argument:  initialfilename is the directory of the file containing the initial values, 
 *					if it is NULL,then using random number generator to generate the INIT structure
 *					chainnum is the number of MCMC chains, determining the number of sets of initial values
 *					popnum is the number of subpopulation assumed
 *					chainnum and popnum determines the dimensions of the "initd" element of the INIT struture
 * Output argument: this application returns an INIT structure, which contains the "chainnum" sets of initial values for F per subpopulation
 * The data in the file should be ranged as:
 *	>chain_name1							each set of the initial values should begin with ">S"
 *	num1 num2..
 *
 *  >chain_name2
 *	....
 */
INIT read_init(char *initialfilename,int chainnum,int popnum,long update,long burnin,int thinning)
{
	//char sign='#';
	int i,j,temp,cnt_chn=0;
	FILE *initfp;
	char *line;
	INIT initial;
	
	initial.chainnum=chainnum;
	initial.update=update;
	initial.burnin=burnin;
	initial.thinning=thinning;
	initial.popnum=popnum;
	
	initial.initd=matrix(0,chainnum-1,0,popnum-1);
	initial.name_len=ivector(0,chainnum-1);
	initial.chn_name=cmatrix(0,chainnum-1,0,MAXLEN-1);
	
	if(initialfilename==NULL)
	{
		for(i=0;i<chainnum;i++)
		{
			for(j=0;j<popnum;j++)
				initial.initd[i][j]=ran1();
			strcpy(initial.chn_name[i],"Chain#");
			strcat(initial.chn_name[i],int_to_string(i+1));			

			initial.name_len[i]=strlen(initial.chn_name[i])+1;
	
			initial.chn_name[i][initial.name_len[i]-1]='\0';
		}
	}
	else{
		if((initfp=fopen(initialfilename,"r"))==NULL)
		{	nrerror("Cannot open inital file!");}
		line=cvector(0,MAXLINE-1);
		cnt_chn=0;
		for(i=0;i<chainnum&&(!feof(initfp));i++)
		{
			while(!feof(initfp))
			{
				fgets(line,MAXLINE,initfp);
				if(line[0]=='>') break;
			}
			for(j=1;j<strlen(line)&&line[j]!='\n';j++)
			{
				initial.chn_name[i][j-1]=line[j];
			}
			initial.chn_name[j]='\0';
			//printf("%s\n",initial.chn_name[i]);
			initial.name_len[i]=j;
			fgets(line,MAXLINE,initfp);
			temp=word_cnt(line);
			if(temp!=popnum)
			{
				nrerror("The number of initial values for selfing rates is not equal the number of subpopulation assumed!\n");
			}
			word_split(line,initial.initd[cnt_chn],popnum);
			cnt_chn++;
		}
		if(cnt_chn<=chainnum)
		{
			for(i=cnt_chn;i<chainnum;i++)
			{
				for(j=0;j<popnum;j++)
					initial.initd[i][j]=ran1();
				strcpy(initial.chn_name[i],"Chain#");
				strcat(initial.chn_name[i],int_to_string(i+1));
				
				initial.name_len[i]=strlen(initial.chn_name[i])+1;
				initial.chn_name[i][initial.name_len[i]-1]='\0';
	
			}
		}
		if(cnt_chn>chainnum)
		{
			nrerror("The number of chain starting points is greater than the number of chains!\n");
		}
		fclose(initfp);
		free_cvector(line,0,MAXLINE-1);
	}
	/*for(i=0;i<chainnum;i++)
	{
		for(j=0;j<popnum;j++)
			fprintf(stdout,"%f\t",initial.initd[i][j]);
		fprintf(stdout,"\n");
	}*/	
	return(initial);
}




/*
 * Function word_cnt counts the number of words in a string
 */
int word_cnt(char *s)
{
	int cnt=0;
	while(*s !='\0')
	{	while(isspace(*s)) ++s;				//skip white space
		if(*s!='\0')						//found a word
		{	++cnt;
			while(!isspace(*s)&&*s!='\0')	//skip word
			++s;
		}
	}
	return cnt;
}




/*
 * Function word_split splits the string into "num" words and makes each word a floating number
 */
void word_split(char *s,float *nums,int num)
{
	int cnt=0,i;
	char *p;
	p=cvector(0,MAXLINE-1);
	while(*s !='\0'&&cnt<num)
	{	while(isspace(*s)) ++s;				//skip white space
		if(*s!='\0')						//found a word
		{	i=0;
			while(!isspace(*s)&&*s!='\0')	//translate the word into float number
			{	p[i++]=*s++;}
			p[i]='\0';
			nums[cnt++]=atof(p);
		}
	}
	free_cvector(p,0,MAXLINE-1);
}



void int_split(char *s,int *nums,int num)
{
	int cnt=0,i;
	char *p;
	p=cvector(0,MAXLINE-1);
	while(*s !='\0'&&cnt<num)
	{	while(isspace(*s)) ++s;				//skip white space
		if(*s!='\0')						//found a word
		{	i=0;
			while(!isspace(*s)&&*s!='\0')	//translate the word into float number
			{	p[i++]=*s++;}
			p[i]='\0';
			nums[cnt++]=atoi(p);
		}
	}
	free_cvector(p,0,MAXLINE-1);
}



/*
 * Function int_to_string takes in an integer number and change it into a string for print
 */
static char* int_to_string(int num)
{
	static char res[20];
	sprintf(res,"%d",num);
	return(res);
}


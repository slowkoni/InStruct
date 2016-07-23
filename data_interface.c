/*
 *  data_interface.c
 *  
 *
 *  Created by Hong Gao on 7/21/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "nrutil.h"
#include "data_interface.h"

#define MAXLEN 1000000		//maximum length of one line in input file
#define ELMLEN 100
#define INCRE_POP 5		/* Incremental number of original populations for which memory is allocated */


static int word_cnt(char *);
static void word_split(char *,char **,int );
static void get_missing(SEQDATA *);
static int isnew_str(char *str, int len,char **array);
static void read_data_fmt1(char *infilename,SEQDATA *data);
static void transform_data(SEQDATA **data, char ****allele);
static void cnt_lines(FILE *infile,SEQDATA *data);
static void cnt_loci(FILE *infile,SEQDATA *data);
static void read_data_from_file(FILE *infile, SEQDATA *data, char *****allele);
static int find_max(int *vec, int length);
static void read_data_fmt2(char *infilename,SEQDATA *data);
static void get_missing_tetra(SEQDATA *data);


SEQDATA read_data(char *infilename,int ploid,int totalsize,int popnum,
	int nloci,char *missingdata,int label,int popdata, double siglevel,
	int back_refl,int type_freq,int nstep_check_empty_cluster,int prior_flag,
	int mode,int n_extra_col,int markername_flag,double alpha_dpm,int print_iter,
	int print_freq,int inf_K,int distr_fmt,int autopoly,int datafmt,double max_mem)
/*
 * Function read_seqs reads in the sequence data and generate missing data matrix and read other parameters
 */
{	
	SEQDATA data;
	
	data.label=label;
	data.popdata=popdata;
	data.n_extra_col=n_extra_col;
	data.markername_flag=markername_flag;
	data.ploid=ploid;
	data.siglevel=siglevel;
	data.back_refl=back_refl;
	data.type_freq=type_freq;
	data.nstep_check_empty_cluster=nstep_check_empty_cluster;
	data.prior_flag=prior_flag;
	data.mode=mode;
	data.popnum=popnum;
	data.missingdata=missingdata;
	data.totalsize=totalsize;
	data.locinum=nloci;
	data.alpha_dpm=alpha_dpm;
	data.print_iter=print_iter;
	data.print_freq=print_freq;
	data.inf_K=inf_K;
	data.distr_fmt=distr_fmt;
	data.autopoly=autopoly;
	data.datafmt=datafmt;
	data.max_mem=max_mem;
	//data.flag_empty_cluster=flag_check_empty;
	
	if(data.ploid==2)
	{
		if(data.datafmt==0)	read_data_fmt1(infilename,&data);
		else{	if(data.datafmt==1)	read_data_fmt2(infilename,&data);}
		get_missing(&data);
	}
	if(data.ploid==4)
	{
		read_data_fmt2(infilename,&data);
		get_missing_tetra(&data);
	}
	data.allelenum_max=find_max(data.allelenum,data.locinum);
	
	return(data);
}




void read_data_fmt1(char *infilename,SEQDATA *data)
/*
 * Function read_data_fmt1 opens the datafile and reads in the sequence data
 * each individual takes two lines, each haploid takes one line
 */
{
	int i;		
	char ****allele;
	FILE *infile;

	if((infile=fopen(infilename,"r"))==NULL) 
	{	nrerror("Cannot open input file!\n");}
	
	cnt_loci(infile,data);	
	
	cnt_lines(infile,data);

	//allocate the space for "allele" which stores the two alleles for each individual at each locus
	allele=(char ****)malloc(data->totalsize*sizeof(char  ***));		 
	if(!allele) 
	{	nrerror("Memory allocation for variable \'allele\' in function read_seqs()!\n");}
	for(i=0;i<data->totalsize;i++)							//seqdata[individual index][locus index][each allele]
	{
		allele[i]=c3tensor(0,data->locinum-1,0,data->ploid-1,0,ELMLEN-1);
	}
	
	read_data_from_file(infile,data, &allele);
	
	fclose(infile); 
	
	transform_data(&data, allele);
		
	for(i=0;i<data->totalsize;i++)
	{
		free_c3tensor(allele[i],0,data->locinum-1,0,data->ploid-1,0,ELMLEN-1);
	}
	free(allele);
}




void read_data_from_file(FILE *infile, SEQDATA *data, char *****allele)
/*
 * read the data from file and store alleles in the array "allele"
 */
{
	int count=0, pop_cnt=0,i,j,k,cnt_token,indx,max_pop=0;
	char *line,**temp;
	
	cnt_token=data->label+data->popdata+data->locinum+data->n_extra_col;	
	line=cvector(0,MAXLEN-1);	
	temp=cmatrix(0,cnt_token-1,0,ELMLEN-1);
		
	//"popdata" records the original population of each individual
	if(data->popdata==1) 
	{	
		data->popindx=ivector(0,data->totalsize-1);	
		if((data->poptype=(char **)malloc(INCRE_POP*sizeof(char *)))==NULL)
		{	nrerror("Memory allocation for variable \'data->poptype\' in function read_data_from_file()!\n");}
		for(i=0;i<INCRE_POP;i++)
		{
			if((data->poptype[i]=(char *)malloc(ELMLEN*sizeof(char)))==NULL)
			{	nrerror("Memory allocation for variable \'data->poptype[i]\' in function read_data_from_file()!\n");}
		}
		max_pop=INCRE_POP;
	}
	
	//"label" store the name of each individual
	if(data->label==1)	
	{
		data->indvname=cmatrix(0,data->totalsize-1,0,ELMLEN-1);
	}
	
	/*
	 * if there is some extra columns of information about each individual
	 * stores the info in "extral_col"
	 */
	if(data->n_extra_col>0)
	{
		data->extra_col=c3tensor(0,data->totalsize-1,0,data->n_extra_col-1,0,ELMLEN-1);
	}
					
	//read the information from input file
	rewind(infile);
	while(!feof(infile))
	{
		fgets(line,MAXLEN,infile);
		if(strlen(line)!=0) break;
	}
	if(data->markername_flag==1)
	{	fgets(line,MAXLEN,infile);	}
	while(!feof(infile))
	{                
		for(i=0;i<data->ploid;i++)
		{	
			if(word_cnt(line)!=cnt_token)
			{	nrerror("The number of tokens in one line does not match the parameters input from the commandline");}
			word_split(line,temp,cnt_token);            
			if(i==0)
			{
				if(data->label==1) strcpy(data->indvname[count],temp[data->label-1]);
				if(data->popdata==1) 
				{
					indx=isnew_str(temp[data->label+data->popdata-1],pop_cnt,data->poptype);
					if(indx==-1)
					{	
						if (max_pop <= pop_cnt)
						{
							if((data->poptype = realloc(data->poptype, (max_pop+=INCRE_POP)*sizeof(char *)) )==NULL)
							{	nrerror("Memory reallocation for variable \'data->poptype\' in function read_data_from_file()!\n");}
							for(j=max_pop-INCRE_POP;j<max_pop;j++)
							{
								if((data->poptype[j] = malloc(ELMLEN*sizeof(char)) )==NULL)
								{	nrerror("Memory reallocation for variable \'data->poptype[i]\' in function read_data_from_file()!\n");}
							}
						}
						pop_cnt++;//count the number of original populations
						
						strcpy(data->poptype[pop_cnt-1],temp[data->label+data->popdata-1]);//record the original population that each individual is from
						data->popindx[count]=pop_cnt-1;
					}
					else{
						data->popindx[count]=indx;
					}
				}
				if(data->n_extra_col>0)
				{
					for(j=0;j<data->n_extra_col;j++)
						strcpy(data->extra_col[count][j],temp[data->label+data->popdata+j]);
				}
			}
			else{
				if(data->label==1) 
				{
					if(strcmp(data->indvname[count],temp[data->label-1])!=0)
					{	nrerror("Some individuals have different number of haplotypes!\n");}
				}
			}
			for(k=data->label+data->popdata+data->n_extra_col;k<cnt_token;k++)
			{
				strcpy((*allele)[count][k-data->label-data->popdata-data->n_extra_col][i],temp[k]);//store the two alleles for each individual at each locus 
			}
			fgets(line,MAXLEN,infile);
			if(word_cnt(line)!=cnt_token) 
			{	nrerror("The lines of input files do not have the same number of tokens!\n");}
		}
		count++;  		     
	}
	
	data->pop_count=pop_cnt;

	free_cmatrix(temp,0,cnt_token-1,0,ELMLEN-1);
	free_cvector(line,0,MAXLEN-1);
}

void read_data_from_file2(FILE *infile, SEQDATA *data, char *****allele)
/*
 * read the data from file and store alleles in the array "allele"
 */
{
	int count=0, pop_cnt=0,i,j,k,cnt_token,indx,max_pop=0;
	char *line,**temp;
	
	cnt_token=data->label+data->popdata+data->n_extra_col+data->locinum*data->ploid;	
	line=cvector(0,MAXLEN-1);	
	temp=cmatrix(0,cnt_token-1,0,ELMLEN-1);
		
	//"popdata" records the original population of each individual
	if(data->popdata==1) 
	{	
		data->popindx=ivector(0,data->totalsize-1);	
		if((data->poptype=(char **)malloc(INCRE_POP*sizeof(char *)))==NULL)
		{	nrerror("Memory allocation for variable \'data->poptype\' in function read_data_from_file()!\n");}
		for(i=0;i<INCRE_POP;i++)
		{
			if((data->poptype[i]=(char *)malloc(ELMLEN*sizeof(char)))==NULL)
			{	nrerror("Memory allocation for variable \'data->poptype[i]\' in function read_data_from_file()!\n");}
		}
		max_pop=INCRE_POP;
	}
	
	//"label" store the name of each individual
	if(data->label==1)	
	{
		data->indvname=cmatrix(0,data->totalsize-1,0,ELMLEN-1);
	}
	
	/*
	 * if there is some extra columns of information about each individual
	 * stores the info in "extral_col"
	 */
	if(data->n_extra_col>0)
	{
		data->extra_col=c3tensor(0,data->totalsize-1,0,data->n_extra_col-1,0,ELMLEN-1);
	}
					
	//read the information from input file
	rewind(infile);
	while(!feof(infile))
	{
		fgets(line,MAXLEN,infile);
		if(strlen(line)!=0) break;
	}
	if(data->markername_flag==1)
	{	fgets(line,MAXLEN,infile);	}
	while(!feof(infile))
	{                
		word_split(line,temp,cnt_token);            
		if(data->label==1) strcpy(data->indvname[count],temp[data->label-1]);
		if(data->popdata==1) 
		{
			indx=isnew_str(temp[data->label+data->popdata-1],pop_cnt,data->poptype);
			if(indx==-1)
			{	
				if (max_pop <= pop_cnt)
				{
					if((data->poptype = realloc(data->poptype, (max_pop+=INCRE_POP)*sizeof(char *)) )==NULL)
					{	nrerror("Memory reallocation for variable \'data->poptype\' in function read_data_from_file()!\n");}
					for(j=max_pop-INCRE_POP;j<max_pop;j++)
					{
						if((data->poptype[j] = malloc(ELMLEN*sizeof(char)) )==NULL)
						{	nrerror("Memory reallocation for variable \'data->poptype[i]\' in function read_data_from_file()!\n");}
					}
				}
				pop_cnt++;//count the number of original populations
				
				strcpy(data->poptype[pop_cnt-1],temp[data->label+data->popdata-1]);//record the original population that each individual is from
				data->popindx[count]=pop_cnt-1;
			}
			else{
				data->popindx[count]=indx;
			}
		}
		if(data->n_extra_col>0)
		{
			for(j=0;j<data->n_extra_col;j++)
				strcpy(data->extra_col[count][j],temp[data->label+data->popdata+j]);
		}
		for(k=0;k<data->locinum;k++)
		{
			for(i=0;i<data->ploid;i++)
				strcpy((*allele)[count][k][i],temp[k*data->ploid+i+data->label+data->popdata+data->n_extra_col]);//store the two alleles for each individual at each locus 
		}
		fgets(line,MAXLEN,infile);
		if(strlen(line)<cnt_token)
		{	
			break;
			fprintf(stdout,"Too many '\n' at the end of input file!\n");
		}
		if(word_cnt(line)!=cnt_token) 
		{	nrerror("The lines of input files do not have the same number of tokens!\n");}
		count++;  		     
	}
	
	data->pop_count=pop_cnt;

	free_cmatrix(temp,0,cnt_token-1,0,ELMLEN-1);
	free_cvector(line,0,MAXLEN-1);
}





void cnt_loci(FILE *infile,SEQDATA *data)
/*
 * count the number of tokens in each line to determine the number of loci
 */	
{
	int cnt_token=0;
	char *line;
	
	line=cvector(0,MAXLEN-1);
	
	while(!feof(infile))
	{
		fgets(line,MAXLEN,infile);
		if(strlen(line)!=0) break;
	}
	cnt_token=word_cnt(line);
	if(data->markername_flag==1)
	{
		data->locinum=cnt_token;
		data->marker_names=cmatrix(0,cnt_token-1,0,ELMLEN-1);
		word_split(line,data->marker_names,cnt_token);
		fprintf(stdout,"The number of loci is %d now!\n",data->locinum);
	}
	else{
		if((cnt_token-data->label-data->popdata-data->n_extra_col)!=data->locinum)
		{	
			data->locinum=cnt_token-data->label-data->popdata-data->n_extra_col;
			fprintf(stdout,"The Input Number of Loci is wrong!\nThe number of loci is %d now!\n",data->locinum);
		}
	}
	free_cvector(line,0,MAXLEN-1);
}


void cnt_loci2(FILE *infile,SEQDATA *data)
/*
 * count the number of tokens in each line to determine the number of loci
 */	
{
	int cnt_token=0;
	char *line;
	
	line=cvector(0,MAXLEN-1);
	
	while(!feof(infile))
	{
		fgets(line,MAXLEN,infile);
		if(strlen(line)!=0) break;
	}
	cnt_token=word_cnt(line);
	if(data->markername_flag==1)
	{
		data->locinum=cnt_token;
		data->marker_names=cmatrix(0,cnt_token-1,0,ELMLEN-1);
		word_split(line,data->marker_names,cnt_token);
		fprintf(stdout,"The number of loci is %d now!\n",data->locinum);
	}
	else{
		if((cnt_token-data->label-data->popdata-data->n_extra_col)/data->ploid!=data->locinum)
		{	
			data->locinum=(cnt_token-data->label-data->popdata-data->n_extra_col)/data->ploid;
			fprintf(stdout,"The Input Number of Loci is wrong!\nThe number of loci is %d now!\n",data->locinum);
		}
	}
	free_cvector(line,0,MAXLEN-1);
}





void cnt_lines(FILE *infile,SEQDATA *data)
/*
 * count the number of lines in data file to determine the number of individuals
 */
{	
	int cnt_line=0;
	char *line;
	
	line=cvector(0,MAXLEN-1);
	
	rewind(infile);
	while(!feof(infile))
	{
		if(fgets(line,MAXLEN,infile)!=NULL)
		{
			if(isspace(line[0])==0&&strlen(line)<=data->locinum&&feof(infile))	{break;}
			else{cnt_line++;}		
		}
	}
	
	if((cnt_line-data->markername_flag)%data->ploid!=0)
	{	nrerror("Some individuals do not have two copies of haplotype!\n");}
	
	if(data->totalsize!=((cnt_line-data->markername_flag)/data->ploid))
	{
		data->totalsize=(cnt_line-data->markername_flag)/data->ploid;
		fprintf(stdout,"The input population size is incorrect!\nThe population size is %d\n",data->totalsize);
	}
	
	free_cvector(line,0,MAXLEN-1);
}


void cnt_lines2(FILE *infile,SEQDATA *data)
/*
 * count the number of lines in data file to determine the number of individuals
 */
{	
	int cnt_line=0;
	char *line;
	
	line=cvector(0,MAXLEN-1);
	
	rewind(infile);
	while(!feof(infile))
	{
		if(fgets(line,MAXLEN,infile)!=NULL)
		{
			if(isspace(line[0])==0&&strlen(line)<=data->locinum&&feof(infile))	{break;}
			else{cnt_line++;}		
		}
	}
	
	if(data->totalsize!=(cnt_line-data->markername_flag))
	{
		data->totalsize=cnt_line-data->markername_flag;
		fprintf(stdout,"The input population size is incorrect!\ncnt_line is %d.\nThe population size is %d\n",cnt_line,data->totalsize);
	}
	
	free_cvector(line,0,MAXLEN-1);
}

void transform_data(SEQDATA **data, char ****allele)
/*
 * use number to represent each alleletype
 */
{
	int i,j,k,m,missing_num=-9,cnt,allele_cnt=0;
	char **alleletype;
	
	(*data)->seqdata=i3tensor(0,(*data)->totalsize-1,0,(*data)->locinum-1,0,(*data)->ploid-1);		//store the two alleles for each individual at each locus 
	
	alleletype=cmatrix(0,((*data)->totalsize*(*data)->ploid)-1,0,ELMLEN-1);
	(*data)->alleletype=(char ***)malloc((*data)->locinum*sizeof(char **));
	if((*data)->alleletype==NULL)
	{
		nrerror("Memory allocation for variable \'(*data)->alleletype\' in function transform_data()!\n");
	}
	
	//record the number of the types of alleles in each locus
	(*data)->allelenum=ivector(0,(*data)->locinum-1);				

	//store each allele type in array "alleletype"
	for(j=0;j<(*data)->locinum;j++)
	{
		cnt=0;
		for(i=0;i<(*data)->totalsize;i++)
		{
			for(k=0;k<(*data)->ploid;k++)
			{
				if(strcmp(allele[i][j][k],(*data)->missingdata)!=0&&isnew_str(allele[i][j][k],cnt,alleletype)==-1)
				{
					strcpy(alleletype[cnt],allele[i][j][k]);
					cnt++;
				}
			}
		}
		if(cnt>=2)	//at least two types of alleles at one locus
		{
			(*data)->allelenum[allele_cnt]=cnt;
			(*data)->alleletype[allele_cnt]=cmatrix(0,cnt-1,0,ELMLEN-1);
			for(k=0;k<cnt;k++)
			{
				strcpy((*data)->alleletype[allele_cnt][k], alleletype[k]);
			}
		
			for(i=0;i<(*data)->totalsize;i++)
			{
				for(k=0;k<(*data)->ploid;k++)
				{
					if(strcmp(allele[i][j][k],(*data)->missingdata)==0)			
						(*data)->seqdata[i][allele_cnt][k]=missing_num;
					else{
						for(m=0;m<cnt;m++)
							if(strcmp(allele[i][j][k],alleletype[m])==0)
								(*data)->seqdata[i][allele_cnt][k]=m;
					}
				}
			}
			allele_cnt++;	
		}
		else{	fprintf(stdout,"The locus %d is not polymorphic.\n",j+1);}		
	}
	(*data)->missingnum=missing_num;
	(*data)->locinum=allele_cnt;
	fprintf(stdout,"The number of polymorphic loci is %d now.\n",allele_cnt);
	//print alleles to standard output	
	fprintf(stdout,"Print the transformed allele data:\n");
	for(i=0;i<(*data)->totalsize;i++)
	{
		for(k=0;k<(*data)->ploid;k++)
		{
			for(j=0;j<(*data)->locinum;j++)
			{
				fprintf(stdout,"%d ",(*data)->seqdata[i][j][k]);
			}
			fprintf(stdout,"\n");
		}
	}
	fprintf(stdout,"End the printing of the transformed allele data.\n");
	
	free_cmatrix(alleletype,0,((*data)->totalsize*(*data)->ploid)-1,0,ELMLEN-1);
}

void transform_data2(SEQDATA **data, char ****allele)
/*
 * use number to represent each alleletype
 */
{
	int i,j,k,m,missing_num=-9,cnt,flag=0;
	char **alleletype;
	
	(*data)->seqdata=i3tensor(0,(*data)->totalsize-1,0,(*data)->locinum-1,0,(*data)->ploid-1);		//store the two alleles for each individual at each locus 
	(*data)->alleleid=imatrix(0,(*data)->totalsize-1,0,(*data)->locinum-1);
	alleletype=cmatrix(0,((*data)->totalsize*(*data)->ploid)-1,0,ELMLEN-1);
	(*data)->alleletype=(char ***)malloc((*data)->locinum*sizeof(char **));
	if((*data)->alleletype==NULL)
	{
		nrerror("Memory allocation for variable \'(*data)->alleletype\' in function transform_data()!\n");
	}
	
	//record the number of the types of alleles in each locus
	(*data)->allelenum=ivector(0,(*data)->locinum-1);				
	for(j=0;j<(*data)->locinum;j++)
		for(i=0;i<(*data)->totalsize;i++)
			for(k=0;k<(*data)->ploid;k++)
				(*data)->seqdata[i][j][k]=-1;
				
	//store each allele type in array "alleletype"
	for(j=0;j<(*data)->locinum;j++)
	{
		cnt=0;
		for(i=0;i<(*data)->totalsize;i++)
		{
			for(k=0;k<(*data)->ploid;k++)
			{
				if(strcmp(allele[i][j][k],(*data)->missingdata)!=0&&isnew_str(allele[i][j][k],cnt,alleletype)==-1)
				{
					strcpy(alleletype[cnt],allele[i][j][k]);
					cnt++;
				}
			}
		}
		(*data)->allelenum[j]=cnt;
		(*data)->alleletype[j]=cmatrix(0,cnt-1,0,ELMLEN-1);
		for(k=0;k<cnt;k++)
		{
			strcpy((*data)->alleletype[j][k], alleletype[k]);
		}
		
		for(i=0;i<(*data)->totalsize;i++)
		{
			flag=0;
			for(k=0;k<(*data)->ploid;k++)
			{
				if(strcmp(allele[i][j][k],(*data)->missingdata)!=0)			
				{
					for(m=0;m<cnt;m++)
					{
						if(strcmp(allele[i][j][k],alleletype[m])==0)
						{
							if(exists(m,(*data)->seqdata[i][j],flag)==0)
							{	(*data)->seqdata[i][j][flag]=m;
								flag++;
							}
						}
					}
				}
			}
			sort(flag,(*data)->seqdata[i][j]);
			(*data)->alleleid[i][j]=flag;
			if(flag==0)
			{	(*data)->seqdata[i][j][0]=missing_num;}
		}		
	}
	(*data)->missingnum=missing_num;
	
	//print alleles to standard output	
	fprintf(stdout,"Print the number of alleles per individual per locus:\n");
	for(i=0;i<(*data)->totalsize;i++)
	{
		for(j=0;j<(*data)->locinum;j++)
		{
			fprintf(stdout,"%d ",(*data)->alleleid[i][j]);
		}
		fprintf(stdout,"\n");
	}
	fprintf(stdout,"Print the transformed allele data:\n");
	for(i=0;i<(*data)->totalsize;i++)
	{
		for(k=0;k<(*data)->ploid;k++)
		{
			for(j=0;j<(*data)->locinum;j++)
			{
				fprintf(stdout,"%d ",(*data)->seqdata[i][j][k]);
			}
			fprintf(stdout,"\n");
		}
	}
	fprintf(stdout,"End the printing of the transformed allele data.\n");
	
	free_cmatrix(alleletype,0,((*data)->totalsize*(*data)->ploid)-1,0,ELMLEN-1);
}

void read_data_fmt2(char *infilename,SEQDATA *data)
//one individual only takes one line
{
	int i;		
	char ****allele;
	FILE *infile;

	if((infile=fopen(infilename,"r"))==NULL) 
	{	nrerror("Cannot open input file!\n");}
	
	cnt_loci2(infile,data);	
	
	cnt_lines2(infile,data);

	//allocate the space for "allele" which stores the two alleles for each individual at each locus
	allele=(char ****)malloc(data->totalsize*sizeof(char  ***));		 
	if(!allele) 
	{	nrerror("Memory allocation for variable \'allele\' in function read_seqs()!\n");}
	for(i=0;i<data->totalsize;i++)							//seqdata[individual index][locus index][each allele]
	{
		allele[i]=c3tensor(0,data->locinum-1,0,data->ploid-1,0,ELMLEN-1);
	}
	
	read_data_from_file2(infile,data, &allele);
	
	fclose(infile); 
	
	if(data->ploid==2)	transform_data(&data, allele);
	else{	if(data->ploid==4)	transform_data2(&data, allele);}
		
	for(i=0;i<data->totalsize;i++)
	{
		free_c3tensor(allele[i],0,data->locinum-1,0,data->ploid-1,0,ELMLEN-1);
	}
	free(allele);

}

void sort(int leng, int *vec)	//sort the vector ascendly
{
	int i,j,temp;
	for(i=0;i<leng-1;i++)
	{
		for(j=i+1;j<leng;j++)
		{	
			if(vec[i]>vec[j])
			{	SWAP(vec[i],vec[j])}
		}
	}
}

void get_missing_tetra(SEQDATA *data)
{
	int i,j;
	
	data->missvec=ivector(0,data->totalsize-1);
	data->missindx=imatrix(0,data->totalsize-1,0,data->locinum-1);//indicating whether alleles are missing at a locus in one's genome
	
	for(i=0;i<data->totalsize;i++)
	{	
		data->missvec[i]=0;
		for(j=0;j<data->locinum;j++)
		{
			data->missindx[i][j]=0;			//'0' means not missing
			if(data->alleleid[i][j]==0)
				data->missindx[i][j]=1;	//'1'means missing, genotype at that locus is discarded in that individual
			data->missvec[i]+=data->missindx[i][j];//record for one individual the number of missing locus
		}
	}

}


int isnew_str(char *str, int len,char **array)
/*
 * Function isnew checks whether array has the element of the same value as "str"
 * It returns '-1' if none of the elements is the same as "num", otherwise returns index of "str" in the "array".
 */
{
	int i,result=-1;
	if(len==0) result=-1;
	else{
		for(i=0;i<len;i++)
		{
			if(strcmp(str,array[i])==0)
			{	result=i;
				break;
			}
		}
	}			
	return(result);
}


int word_cnt(char *s)
/*
 * Function word_cnt counts the number of words in a string
 */
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




void word_split(char *s,char **words,int num)
/*
 * Function word_split splits the string into "num" words
 */
{
	int cnt=0,i;
	while(*s !='\0')
	{	while(isspace(*s)) ++s;				//skip white space
		if(*s!='\0')						//found a word
		{	i=0;
			while(!isspace(*s)&&*s!='\0')	//translate the word into an integer number
			{	words[cnt][i++]=*s++;}
			words[cnt][i]='\0';
			cnt++;
		}
	}
	if(num!=cnt)
	{
		nrerror("The number of tokens in a line is different from the number given !\n");
	}
}






void get_missing(SEQDATA *data)
/*
 * generate the matrix and vector that stores the missing data
 */
{
	int i,j,k;
	
	data->missvec=ivector(0,data->totalsize-1);
	data->missindx=imatrix(0,data->totalsize-1,0,data->locinum-1);//indicating whether alleles are missing at a locus in one's genome
	
	for(i=0;i<data->totalsize;i++)
	{	
		data->missvec[i]=0;
		for(j=0;j<data->locinum;j++)
		{
			data->missindx[i][j]=0;			//'0' means not missing
			for(k=0;k<data->ploid;k++)
			{
				if(data->seqdata[i][j][k]==data->missingnum)
					data->missindx[i][j]=1;	//'1'means missing, genotype at that locus is discarded in that individual
			}
			data->missvec[i]+=data->missindx[i][j];//record for one individual the number of missing locus
		}
	}
	
	/*print to standard output
	for(i=0;i<data->totalsize;i++)
	{
		for(j=0;j<data->locinum;j++)
		{
			fprintf(stdout,"%d ",data->missindx[i][j]);
		}fprintf(stdout," %d\n",data->missvec[i]);
	}
	*/
}


int find_max(int *vec, int length)
/*
 *	return the largest value in the vector "vec"
 */
{
	int x,i;
	x=vec[0];
	for(i=1;i<length;i++)
	{
		if(x<vec[i])
			x=vec[i];
	}
	return(x);
}


int exists(int value, int *vec, int leng)
//whether "value" exists in the vector "vec" (1) or not (0)
{
	int i,flag=0;
	if(leng==0) return(0);
	else{
		for(i=0;i<leng;i++)
		{
			if(value==vec[i])
			{	flag=1;}
		}
		return(flag);
	}	
}



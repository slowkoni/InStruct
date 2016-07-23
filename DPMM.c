/*
 *  DPMM.c
 *  
 *
 *  Created by Hong Gao on 11/13/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nrutil.h"
#include "random.h"
#include "data_interface.h"
#include "initial.h"
#include "check_converg.h"
#include "mcmc.h"
#include "DPMM.h"

#define EPS 1.0e-6
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5
#define FUNC(x,a,b,c) ((*func)(x,a,b,c))
#define NRANSI

static NODE* find(NODE **head,double self);
static NODE* creat(NODE **head,double self,int index,int *id);
static void insert(NODE *ptr,int index,int *id);
static int delete(NODE **head,SF *pointer,SF** indv_array);
static void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
static double trapzd(double (*func)(double,UPMCMC *,SEQDATA,int), double a, double b, int n,UPMCMC *ptr,SEQDATA data,int id);
static double qromb(double (*func)(double,UPMCMC *,SEQDATA,int), double a, double b,UPMCMC *ptr,SEQDATA data,int id);
static double sample_poster(UPMCMC *ptr, SEQDATA data,int);
static double gen_nonconjg(UPMCMC *ptr,SEQDATA data,int);
static void gen_post_prob(double *prob_interval,NODE **head,int *cnt_node,UPMCMC *ptr,SEQDATA,int id);
static double func(double x,UPMCMC *ptr,SEQDATA data,int i);

void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=dvector(1,n);
	d=dvector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_dvector(d,1,n);
	free_dvector(c,1,n);
}
#undef NRANSI


double trapzd(double (*func)(double,UPMCMC *,SEQDATA,int), double a, double b, int n,UPMCMC *ptr,SEQDATA data,int id)
{
	double x,tnm,sum,del,temp;
	static double s;
	int it,j;	//FUNC is log function, thus needs to exponentiate it.

	if (n == 1) {
		return (s=0.5*(b-a)*(exp(FUNC(a,ptr,data,id))+exp(FUNC(b,ptr,data,id))));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		temp=FUNC(x,ptr,data,id);
		for (sum=0.0,j=1;j<=it;j++,x+=del)
		{ sum += exp(FUNC(x,ptr,data,id)-temp);}
		sum*=exp(temp);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC


double qromb(double (*func)(double,UPMCMC *,SEQDATA,int), double a, double b,UPMCMC *ptr,SEQDATA data,int id)
{
	double ss,dss;
	double s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd(func,a,b,j,ptr,data,id);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K


void init_DP(NODE** head,double alpha,SF **indv_array,int *cnt_node,int num)
/*
 * Initially draw the selfing rate of each individual from the Dirichlet Process Prior
 */
{
	int i,j,k,prob_tmp;
	NODE* ptr_node;
	double *prob_interval;
	prob_interval=dvector(0,num);
	
	*cnt_node=0;									//"cnt_node" is used to record the number of nodes in the list
	for(j=0;j<num;j++)			
	{									
		prob_interval[0]=alpha/(alpha+(double)j);	//The probability that s_j is randomly drawn from a uniform distribution
		ptr_node=*head;
		for(i=1;i<=*cnt_node&&ptr_node!=NULL;i++)
		{											//The probability that s_j is the same as one of existing selfing rates
			prob_interval[i]=prob_interval[i-1]+(double)ptr_node->num/(alpha+(double)j);
			ptr_node=ptr_node->next;
		}													//According to the "prob_tmp" to determine whether to draw a new value for s_j or let it equal to one of previous selfing rates
		prob_tmp=disc_unif(prob_interval,*cnt_node+1);
		if(prob_tmp==0)
		{
			(*indv_array)[j].value=ran1();		//A new node is added the list containing the newly-drawn selfing rate
			(*indv_array)[j].p=creat(head,(*indv_array)[j].value,j,&((*indv_array)[j].index));
			(*cnt_node)++;
		}
		else{										//Insert this individual into one of existing nodes with the same selfing rate
			ptr_node=*head;
			for(k=0;k<prob_tmp-1&&ptr_node!=NULL;k++)
			{	ptr_node=ptr_node->next;}
			(*indv_array)[j].value=ptr_node->value;						
			insert(ptr_node,j,&((*indv_array)[j].index));
			(*indv_array)[j].p=ptr_node;
		}
	}
	free_dvector(prob_interval,0,num);
}



void update_DP(NODE** head,double alpha,SF **indv_array,int *cnt_node,int totalsize,SEQDATA data,UPMCMC *ptr)
/*
 * Update the selfing rates with the Dirichlet Process Prior
 */
{
	int j,k,prob_tmp,int_tmp;
	NODE* ptr_node;
	double *prob_interval;
	prob_interval=dvector(0,totalsize);
	
	for(j=0;j<totalsize;j++)			
	{	
		int_tmp=delete(head,&((*indv_array)[j]),indv_array);	//The individual j is deleted from the list
		*cnt_node-=int_tmp;	//the number of nodes might be decreased by one or not
		
		gen_post_prob(prob_interval,head,cnt_node,ptr,data,j);	//fill the array prob_interval with posterior distribution
				
		prob_tmp=disc_unif(prob_interval,*cnt_node+1);				
		if(prob_tmp==0)					//According to the "prob_tmp" to determine whether to draw a new value for s_j or let it equal to one of previous selfing rates
		{								//A new node is added the list containing the newly-drawn selfing rate
			(*indv_array)[j].value=sample_poster(ptr,data,j);	
			(*indv_array)[j].p=creat(head,(*indv_array)[j].value,j,&((*indv_array)[j].index));
			(*cnt_node)++;
		}
		else{
			ptr_node=*head;				//Insert this individual into one of existing nodes with the same selfing rate
			for(k=0;k<prob_tmp-1;k++)
			{	ptr_node=ptr_node->next;}
			(*indv_array)[j].value=ptr_node->value;
			insert(ptr_node,j,&((*indv_array)[j].index));
			(*indv_array)[j].p=ptr_node;
		}
	}
	free_dvector(prob_interval,0,totalsize);
}		


NODE* find(NODE **head,double self)
/*
 *find the node with element "value" less than or equal "self" and return the found node address.
 */
{
	NODE *ptr=NULL,*qtr=NULL;
	ptr=qtr=*head;
	if(*head==NULL) return(NULL);
	if(self<=(*head)->value)
	{	return(NULL);	}
	else{
		while(ptr!=NULL&&ptr->value<=self)
		{
			qtr=ptr;
			ptr=ptr->next;
		}
		return(qtr);
	}	
}



NODE * creat(NODE **head,double self,int index,int *id)
/*
 * creat a new node and insert it into the list
 */
{
	NODE *ptr=NULL,*pointer=NULL;
	if(*head==NULL)		//If the list is empty
	{	
		if((ptr=(NODE *)malloc(sizeof(NODE)))==NULL)
		{	nrerror("Allocation failure in ptr in creat()!");}
		ptr->value=self;
		ptr->num=1;
		ptr->index=ivector(0,ptr->num-1);
		ptr->index[ptr->num-1]=index;
		
		ptr->next=NULL;
		*head=ptr;
		*id=ptr->num-1;
	}
	else{
		if((ptr=(NODE *)malloc(sizeof(NODE)))==NULL)
		{	nrerror("Allocation failure in ptr in creat()!");}
		ptr->value=self;
		ptr->num=1;
		ptr->index=ivector(0,ptr->num-1);
		ptr->index[ptr->num-1]=index;
		*id=ptr->num-1;
		pointer=find(head,self);
		if(pointer==NULL)
		{
			ptr->next=*head;
			*head=ptr;
		}
		else{			
			ptr->next=pointer->next;
			pointer->next=ptr;			
		}
	}
	return(ptr);
}


void insert(NODE *ptr,int index,int *id)
/*
 *insert an individual into one of existing nodes
 */
{
	if((ptr->index=(int *)realloc(ptr->index,(ptr->num+1)*sizeof(int)))==NULL)
	{	nrerror("Allocation failure in ptr in insert()!");}
	ptr->index[ptr->num]=index;
	*id=ptr->num;
	ptr->num++;
}



int delete(NODE **head,SF *pointer,SF** indv_array)
/*
 * delete the individual and any node that its element "num" is zero
 */
{
	int i,flag=0;
	NODE *ptr,*qtr;
	ptr=pointer->p;
	if(ptr==NULL)
	{	nrerror("This node does not contain the individual to be deleted!");}
	for(i=pointer->index+1;i<ptr->num;i++)
	{
		ptr->index[i-1]=ptr->index[i];
		(*indv_array)[ptr->index[i-1]].index=i-1;
	}	
	ptr->num--;
	
	ptr=qtr=*head;
	if(*head==NULL) return(0);
	if(ptr->num==0)
	{
		*head=ptr->next;
		free(ptr);
		flag=1;
	}
	else{
		while(ptr!=NULL)
		{
			if(ptr->num==0)
			{
				qtr->next=ptr->next;
				free(ptr);
				ptr=qtr->next;
				flag=1;
				if(ptr==NULL) break;
			}
			qtr=ptr;
			ptr=ptr->next;
		}		
	}
	return(flag);
}





double func(double x,UPMCMC *ptr,SEQDATA data,int i)
{
	int j,m;
	double temp=0,*tmp,ld=0;					/*calculate the probability */
	tmp=dvector(0,data.ploid-1);
	
	for(j=0;j<data.locinum;j++)		
	{
		if(data.missindx[i][j]!=1)
		{	
			if(chcksame(ptr->z[i][j],data.ploid)==0)
			{
				for(m=0;m<data.ploid;m++)
				{
					tmp[m]=ptr->freq[ptr->z[i][j][m]][j][data.seqdata[i][j][m]];
				}
				temp=genofreq_inbreedcoff(data.seqdata[i][j],tmp,x,data.ploid);
				ld+=log(temp);
			}
			else{
				for(m=0;m<data.ploid;m++)
				{
					ld+=log(ptr->freq[ptr->z[i][j][m]][j][data.seqdata[i][j][m]]);
				}
				if(chcksame(data.seqdata[i][j],data.ploid)==1) //heterozygote needs to time 2
				{	ld+=log(2);}
			}
		}                                     
	}
	free_dvector(tmp,0,data.ploid-1);
	return(ld);
}


void gen_post_prob(double *prob_interval,NODE **head,int *cnt_node,UPMCMC *ptr,SEQDATA data,int id)
{
	NODE* ptr_node;
	int i;
	double temp,tmp_prob;
	
	if(data.mode==3)	//For selfing rates, prior is uniform and the likelihood is geometric(generation), thus the posterior is beta distribution.
	{	
		*(prob_interval+0)=data.alpha_dpm/(ptr->generation[id]+1)/ptr->generation[id];//The probability that s_j is randomly drawn from a uniform distribution
		ptr_node=*head;
		for(i=1;i<=*cnt_node&&ptr_node!=NULL;i++)
		{		//The probability that s_j is the same as one of existing selfing rates
			tmp_prob=ptr_node->num*dgeom(ptr_node->value,ptr->generation[id]);
			*(prob_interval+i)=*(prob_interval+i-1)+tmp_prob;
			ptr_node=ptr_node->next;
		}
	}	
	if(data.mode==5)
	{
		temp=log(data.alpha_dpm)+log(qromb(func,0,1,ptr,data,id));//The probability that s_j is randomly drawn from a uniform distribution
		*(prob_interval+0)=1;
		ptr_node=*head;
		for(i=1;i<=*cnt_node&&ptr_node!=NULL;i++)
		{	//The probability that s_j is the same as one of existing selfing rates
			tmp_prob=log(ptr_node->num)+func(ptr_node->value,ptr,data,id);
			*(prob_interval+i)=*(prob_interval+i-1)+exp(tmp_prob-temp);
			ptr_node=ptr_node->next;
		}
	}
}

double sample_poster(UPMCMC *ptr, SEQDATA data,int j)
{
	double temp=0;
	if(data.mode==3) temp=rbeta(ptr->generation[j],2);
	if(data.mode==5) temp=gen_nonconjg(ptr,data,j);
	return(temp);
}

//sample from a non-conjugate posterior distribution
double gen_nonconjg(UPMCMC *ptr,SEQDATA data,int id)
{
	double temp=0;
	
	
	return(temp);
}



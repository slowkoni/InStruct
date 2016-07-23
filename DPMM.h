 /*
 *  DPMM.h
 *  
 *
 *  Created by Hong Gao on 11/13/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

typedef struct node		//construct a list of nodes which stores different values
{
	double value;
	int num;
	int *index;
	struct node *next;
} NODE;

typedef struct self_indv
{
	double value;
	int index;
	NODE *p;
} SF;

void init_DP(NODE** head,double alpha,SF **indv_array,int *cnt_node,int);
void update_DP(NODE** head,double alpha,SF **indv_array,int *cnt_node,int totalsize,SEQDATA data,UPMCMC *ptr);

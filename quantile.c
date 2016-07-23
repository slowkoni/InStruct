#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "quantile.h"

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50
#define NR_END 1
#define FREE_ARG char*



/*
 * Function indexx() sorts the array "arr" using quick sort algorithm 
 * and stores the order in "indx" while the array "arr" remains unchanged
 * "n" is the length of the array "arr"
 */
void indexx(int n, double arr[], int indx[])
{
	unsigned long i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0,*istack;
	double a;
	
	istack=ivector(1,NSTACK);
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=l;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[l+1]]) {
				SWAP(indx[l],indx[l+1])
			}
			i=l+1;
			j=ir;
			indxt=indx[l+1];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j])
			}
			indx[l+1]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in indexx.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,NSTACK);
}


/*
 * Function quantile seeks the number with certain quantile in an array "thedata"
 * Input argument:  numrows is the length of the array "thedata"
 *					thedata is an unordered array
 *					indx storess the indeces of elements in "thedata" and when sorting
 *					the array, indx is changed to store the order while "thedata" remains unchanged
 *					percent is the quantile, such as "0.05"
 * It return the number of certain quantile in "thedata".
 */
double quantile(int numrows, double thedata[], int indx[], float percent)
{
  int y;
  y = numrows*percent;
  indexx(numrows,thedata,indx);
  return thedata[indx[y]];
}



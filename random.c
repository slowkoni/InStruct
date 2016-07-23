#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include "random.h"

#define E 2.71828182
#define PI 3.141592654

#define SEED1 13           /* Seed for wichman() */
#define SEED2 4
#define SEED3 1972

static long int seed1=SEED1;  /* Seeds declared externally to avoid passing */
static long int seed2=SEED2;  /* each time wichmann is called               */
static long int seed3=SEED3;


double wichmann()
/*
 * Random number generator for U(0,1) distribution.
 */
{
   double random;
   seed1 = (171 * seed1)%30269;
   seed2 = (172 * seed2)%30307;
   seed3 = (170 * seed3)%30323;
 
   random = fmod(seed1/30269.0 + seed2/30307.0 + seed3/30323.0, 
     1.0);
   return random;
}

double ran1()
/*
 * Random number generator for U(0,1) distribution.
 */
{
   double random;
   seed1 = (171 * seed1)%30269;
   seed2 = (172 * seed2)%30307;
   seed3 = (170 * seed3)%30323;

   random = fmod(seed1/30269.0 + seed2/30307.0 + seed3/30323.0, 
     1.0);
   return random;
}


void setseeds(int sd1, int sd2,int sd3)
/*
 * set the seeds new values by users
 */
{
	seed1=sd1;
	seed2=sd2;
	seed3=sd3;
}

void printseeds(FILE* fp)
{
	fprintf(fp,"%ld %ld %ld\n",seed1,seed2,seed3);
}


float gammln(float xx)					// returns ln[Gamma[xx]] for xx > 0
{
	double x,tmp,ser;
	static double cof[6]={76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.120858003e-2,
		-0.536382e-5};
	int j;
	x = xx-1.0;
	tmp = x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser = 1.0;
	for (j=0; j<=5; j++) {
		x += 1.0;
		ser += cof[j]/x;
	}
	return -tmp+log(2.50662827465*ser);
}

float factln(int n)					// returns ln(n!)
{
	static float a[101];
	if (n <= 1) return 0.0;
	if (n <= 500) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
	else return gammln(n+1.0);
} 

float bico(int n, int k)			// returns binomial coefficients (nCk)
{
	if (k == 1) return n;
	if (k == 0) return 1.0;
	else return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}

void moment(double data[], int n, float *ave, float *var)
/*
 * Calculate the mean and variance of "data"
 */
{
	void nrerror(char error_text[]);
	int j;
	float s;

	if (n <= 1) nrerror("n must be at least 2 in moment");
	s=0.0;
	(*var) = 0.0;
	for (j=0;j<n;j++) 
	  {
	    s += data[j];
	    *var += data[j]*data[j]/n;
	  }
	*ave=s/n;
  	*var=(*var-(*ave)*(*ave))*n/(n-1);
}


/////////////////////////
double rexp(double lambda)
/*
 * Generates from an exponential distribution
 */
{
   double random, uniform;
   uniform = wichmann();
   random = - (1/lambda) * log(uniform);
   return random;
}

float poidev(float xm)                  // Returns a poisson RV with mean xm 
{
        static float sq,alxm,g,oldm=(-1.0);
        float em,t,y;

        if (xm < 12.0) {
                if (xm != oldm) {
                        oldm=xm;
                        g=exp(-xm);
                }
                em = -1;
                t=1.0;
                do {
                        ++em;
                        t *= wichmann();
                } while (t > g);
        } else {
                if (xm != oldm) {
                        oldm=xm;
                        sq=sqrt(2.0*xm);
                        alxm=log(xm);
                        g=xm*alxm-gammln(xm+1.0);
                }
                do {
                        do {
                                y=tan(PI*wichmann());
                                em=sq*y+xm;
                        } while (em < 0.0);
                        em=floor(em);
                        t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
                } while (wichmann() > t);
        }
        return em;
}

double rgamma1(double alpha)
/* 
 * Generates from a gamma distribution with alpha < 1
 */
{
   double uniform0, uniform1;
   double random, x;
   uniform0 = wichmann();
   uniform1 = wichmann();
   if (uniform0 > E/(alpha + E))
   {
      random = -log((alpha + E)*(1-uniform0)/(alpha*E));
      if ( uniform1 > pow(random,alpha - 1))
         return -1;
      else 
         return random;
   }
   else
   {
      x = (alpha + E) * uniform0 / E;
      random = pow(x,1/alpha);
      if ( uniform1 > exp(-random))
         return -1;
      else
         return random;
   } 
}

double rgamma2(double alpha)
/*
 * Generates from a gamma distribution with alpha > 1
 */
{
   double uniform1,uniform2;
   double c1,c2,c3,c4,c5,w;
   double random;
   int done = 1;
   c1 = alpha - 1;
   c2 = (alpha - 1/(6 * alpha))/c1;
   c3 = 2 / c1;
   c4 = c3 + 2;
   c5 = 1 / sqrt(alpha);
   do
   {
      uniform1 = wichmann();
      uniform2 = wichmann();
      if (alpha > 2.5)
      {
          uniform1 = uniform2 + c5 * (1 - 1.86 * uniform1);
      }
   }
   while ((uniform1 >= 1) || (uniform1 <= 0));
   w = c2 * uniform2 / uniform1;
   if ((c3 * uniform1 + w + 1/w) > c4)
   {
      if ((c3 * log(uniform1) - log(w) + w) >= 1)
      {
         done = 0;
      }
   }
   if (done == 0)
      return -1;
   random = c1 * w; 
   return random;
} 

double rgamma(double alpha, double beta)
/*
 * Generates from a general gamma(alpha,beta) distribution
 */   
{
   double random;
   if (alpha < 1)
      do {
      random = rgamma1(alpha)/beta; 
      } while (random < 0 );
   if (alpha == 1)
      random = rexp(1)/beta; 
   if (alpha > 1)
      do {
      random = rgamma2(alpha)/beta; 
      } while (random < 0);
   return random;
}

double rbeta(double alpha,double beta)
/*
 * Generates from a beta (alpha,beta) distribution
 */ 
{
	double tmp=1,random,temp;
	temp=rgamma(alpha,tmp);
	random=temp/(temp+rgamma(beta,tmp));
	return(random);
}


void rdirich(double *alpha,int length,double **rand,double add)
/*
 * Generates from a Dirichlet distribution,the generated random values are stored in "rand",add is psuedocount.
 */ 
{
	double tmp,sum=0,beta=1.0;
	int k;
	for(k=0;k<length;k++)
	{
		tmp=rgamma(alpha[k]+add,beta);
		(*rand)[k]=tmp;
		sum+=tmp;				
	}
	for(k=0;k<length;k++)
		(*rand)[k]/=sum;
			
}


double rstd_normal()
/*
 * Generates from a standard normal(0,1) distribution
 */ 
{
   double uniform1,uniform2;
   double theta,r;
   double random;
   uniform1 = wichmann();
   uniform2 = wichmann();
   theta = 2 * PI * uniform1;
   r = sqrt(2 * ( - log(uniform2)));
   random = r * cos(theta);
   return random;
}
   
double rnormal(double mean, double sd)
/*
 * Generates from a general normal(mu,sigma) distribution
 */
{
   double random;
   random = mean + sd * rstd_normal();
   return random;
}



int rgeom(double p)
/*
 * Generates from a Geometric(p) distribution,p is the probability of success
 */
{
  double u;
  int x;
  u=wichmann();
  x=(int)(log(u)/log(1-p))+1;
  return x;
}

int rbinom(double p, int N)
/*
 * Generates from a binomial (N,p) distribution
 */
{
  int k;
  double alpha, beta, U;
  int sum;
        
  alpha = 1.0/p;
  beta = 1.0/(1-p);
  sum = 0.0;
  k = 0;
  U = wichmann();
  while(k < N)
    {
//           printf("%f %f %i %i %i\n", p, U, k, N, sum);
      k+=1;
      if(U <= p)
        {
          sum += 1;
          U *= alpha;   
        } 
      else
        {
          U = beta*(U - p);
        }   
    }     
  return sum;
}
 
float bnldev(float pp, int n)
{
	int j;
	static int nold=(-1);
	float am,em,g,angle,p,bnl,sq,t,y;
	static float pold=(-1.0),pc,plog,pclog,en,oldg;

	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;
	if (n < 25) {
		bnl=0.0;
		for (j=1;j<=n;j++)
			if (wichmann() < p) ++bnl;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= wichmann();
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg=gammln(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*wichmann();
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=floor(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
				-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (wichmann() > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
}

int disc_unif(double *vec,int length)
/*
 * A interval is divided by values of "vec" into "length" subintervals.
 * Find the subinterval that a random value "x" falls in and return the index of the subinterval.
 */
{
	int i=0,flag=0;
	double x;
	x=ran1();
	for(i=0;i<length;i++)		//make sure vec is within [0,1]
		vec[i]/=vec[length-1];
		
	if(x<0.00||x>vec[length-1])
	{	nrerror("The value x is outside the interval!");}
	if(x<=vec[0]&&x>=0.00)
	{	flag=0;}
	else{
		for(i=1;i<length;i++)
		{
			if(x>vec[i-1]&&x<=vec[i])
			{
				flag=i;
			}
		}
	}
	
	return(flag);
}




















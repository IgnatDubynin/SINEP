#include "sinepvar.h"
#include <stdlib.h>
extern int nd2,nep;
int x0,x1,x2,x3,nd2_n;
long double * divv, * upper;

void ini_Lagrange(void)
{
int i,j;
int * xx;
int N;

/* границы области интерполяции */


nd2=NPOINTS-LAST_EP-1;
x1=FIRST_EP;    x2=LAST_EP;
nep=x2-x1+1;
//x0=x1-nep;
x0=x1-10;
//if (x0>LABEL/2) x0=LABEL/2;
if (x0<20) x0=20;
nd2_n=x1-x0;
x3=x2+nd2_n;
if (x3>NPOINTS-1)
  {
  nd2_n=NPOINTS-1-x2;
  x0=x1-nd2_n;
  x3=x2+nd2_n;
  }
N=nd2_n*2;
divv=(long double*)malloc(N*sizeof(long double));
upper=(long double*)malloc(nep*sizeof(long double));
xx=(int*)malloc(N*sizeof(int));

for (i=0;i<nd2_n;i++)     xx[i]=i;
for (i=nd2_n;i<N;i++) xx[i]=i+nep;


for (i=0;i<nep;i++)
  {
  upper[i]=1.0;
  for (j=0;j<N;j++)
    upper[i]*=(long double)(nd2_n+i-xx[j]);
  }


for (i=0;i<N;i++)
  {
  divv[i]=1.0;
  for (j=0;j<N;j++)
    if (i!=j)
      divv[i]*=(long double)(xx[i]-xx[j]);
  }

free(xx);

}

void Lagrange(double* inp,double* out)
{
int i,j,k;
long double r,rr;

for (i=0;i<x1;i++)         out[i]=inp[i];
for (i=x2+1;i<NPOINTS;i++) out[i]=inp[i];

for (i=x1;i<=x2;i++)
  {
  r=0;
  for (j=0,k=x0;j<nd2_n;j++,k++)
    {
    rr=upper[i-x1]/divv[j]/(i-k);
    r+=upper[i-x1]/divv[j]*inp[k]/(i-k);
    }
  for (j=nd2_n,k=x2+1;j<nd2_n*2;j++,k++)
    {
    rr=upper[i-x1]/divv[j]/(i-k);
    r+=upper[i-x1]/divv[j]*inp[k]/(i-k);
    }
  out[i]=r;
  }

}



void del_Lagrange(void)
{
free(divv);
free(upper);
}


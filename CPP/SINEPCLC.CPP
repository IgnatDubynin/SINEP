#include "sinepvar.h"
#include <math.h>
#include <stdio.h>
#include "snp_drow.h"

void min_max(float* mas, float* minr,float* maxr)
{
int i,j;
long ind=0;
float min,max,r;
min=max=mas[0];
for (i=0;i<NPOINTS;i++)
  for(j=0;j<NCAN;j++,ind++)
    {
    r=mas[ind];
    if (min>r) min=r;
    if (max<r) max=r;
    }
*minr=min;
*maxr=max;
}
/*****************************************************************/
/*********  метод комплексной модуляции-демодуляции  *************/
/*****************************************************************/

double s1,s2,c1,c2,kb;
void ini_oscillator(float betta)
  {
  s1=-sin((double)betta);
  s2=-sin((double)betta*2);
  c1= cos((double)betta);
  c2= cos((double)betta*2);
  kb=c1*2;
  }

float SIN (void)
  {
  double s0=s1*kb-s2;
  s2=s1; s1=s0;
  return s0;
  }
float COS (void)
  {
  double c0=c1*kb-c2;
  c2=c1; c1=c0;
  return c0;
  }

void modulation(float betta,float* mas,float* mas_cos,float* mas_sin)
{
int i,j;
long ind=0;
float sn,cn;

ini_oscillator(betta);
for (i=0;i<NPOINTS;i++)
  {
  sn=SIN();  cn=COS();
  for (j=0;j<NCAN;j++,ind++)
    {
    mas_cos[ind]=mas[ind]*cn;
    mas_sin[ind]=mas[ind]*sn;
    }
  }
}


void demodulation(float betta,float* mas,float* mas_cos,float* mas_sin)
{
int i,j;
long ind=0;
float sn,cn;

ini_oscillator(betta);
for (i=0;i<NPOINTS;i++)
  {
  sn=SIN();  cn=COS();
  for (j=0;j<NCAN;j++,ind++)
    mas[ind]=mas_cos[ind]*cn+mas_sin[ind]*sn;
  }
}

void summing(float* sum,float* mas1,float* mas2)
{
int i,j;
long ind;
for (i=0,ind=0;i<NPOINTS;i++)
  for (j=0;j<NCAN;j++,ind++)
    sum[ind]=mas1[ind]+mas2[ind];
}

void subtracking(float* sum,float* mas1,float* mas2)
{
int i,j;
long ind;
for (i=0,ind=0;i<NPOINTS;i++)
  for (j=0;j<NCAN;j++,ind++)
    sum[ind]=mas1[ind]-mas2[ind];
}



const int N21=21; // максимальный порядок многочлена
float * orto[20];
float * xx;
int M;
int nep,nd2;
const char er0[]="Не хватает памяти";

void ini_Gram(int m)
{
int i,j,k;
float r;
float coeff[N21];
double * rv,rr;
M=m;
/* создадим массив аргументов */
//x1=FIRST_EP;    x2=LAST_EP;
nep=LAST_EP-FIRST_EP+1;
nd2=NPOINTS-LAST_EP-1;
//xx=(float*)farmalloc(NPOINTS*sizeof(float));
xx=new float [NPOINTS];
if (xx==NULL) error((char far *)er0);
rr=(float)2/NPOINTS;
for (i=0,r=-1.0;i<NPOINTS;i++,r+=rr)       xx[i]=r;

for (i=0;i<=M;i++)
//  orto[i]=(float*)farmalloc(NPOINTS*sizeof(float));
    orto[i]=new float [NPOINTS];
if (orto[M-1]==NULL) error((char far *)er0);

/* многочлен нулевой степени */
r=sqrt(1.0/(nd2+FIRST_EP));
for (i=0;i<NPOINTS;i++) orto[0][i]=r;

/* вычисление многочленов старших порядков */
//rv=(double*)farmalloc(NPOINTS*sizeof(double));
rv= new double [NPOINTS];
if (rv==NULL) error((char far *)er0);

for (k=1;k<=M;k++)
  {
  for (i=0;i<NPOINTS;i++) rv[i]=pow(xx[i],k);
  for (j=0;j<k;j++)
    {
    rr=0;
    for (i=0;i<FIRST_EP;i++)         rr+=rv[i]*orto[j][i];
    for (i=LAST_EP+1;i<NPOINTS;i++)  rr+=rv[i]*orto[j][i];
    coeff[j]=rr;
    }

  for (j=0;j<k;j++)
    for (i=0;i<NPOINTS;i++)
      rv[i]-=coeff[j]*orto[j][i];

  rr=0;
  for (i=0;i<FIRST_EP;i++)         rr+=rv[i]*rv[i];
  for (i=LAST_EP+1;i<NPOINTS;i++)  rr+=rv[i]*rv[i];
  rr=sqrt(rr);
  for (i=0;i<NPOINTS;i++) orto[k][i]=rv[i]/rr;

  }

//farfree(rv);
delete rv;
}


void Gram(float* inp,float* out)
{
int i,k;
float coeff[N21],r,rr;

for (k=0;k<=M;k++)
  {
  r=0;
  for (i=0;i<FIRST_EP;i++)        r+=inp[i]*orto[k][i];
  for (i=LAST_EP+1;i<NPOINTS;i++) r+=inp[i]*orto[k][i];
  coeff[k]=r;
  }

for (i=0;i<NPOINTS;i++)
  {
  out[i]=0;
  for (k=0;k<=M;k++)
    out[i]+=coeff[k]*orto[k][i];
  }

rr=1.0/FIRST_EP;
for (i=0,r=1.0;i<FIRST_EP;i++,r-=rr)
  out[i]=inp[i]*r+out[i]*(1.0-r);

k=NPOINTS-LAST_EP-1;
rr=(k<=0) ? 1.0 : 1.0/k;
for (i=LAST_EP+1,r=0;i<NPOINTS;i++,r+=rr)
  out[i]=inp[i]*r+out[i]*(1.0-r);

}




void del_Gram(void)
{
int i;

for (i=M;i>0;i--)
//  farfree(orto[i]);
   delete orto[i];

//farfree(xx);
delete xx;
}

/*****************************************************************/

float     Bz=0.00391,
	  A1=-1.974556,         A3=-0.975182,
	  A2=0.975181;

void filtr(float * mas)
{
int i,j,k;
long ind;
float r;

/* "зануление" концов */
ind=(long)NPOINTS*NCAN-1;

for (i=0,k=0;i<20;i++)
  {
  r=1.0-cos((float)i*M_PI/40.0);
  for (j=0;j<NCAN;j++,k++,ind--)
    {
    mas[k]*=r;
    mas[ind]*=r;
    }
  }

/* фильтрация в прямом направлении */
ind=NCAN*2;
for (i=2;i<NPOINTS;i++)
  for (j=0;j<NCAN;j++,ind++)
    mas[ind]=Bz*mas[ind]-A1*mas[ind-NCAN]-A2*mas[ind-NCAN-NCAN];
ind=NCAN*2;
for (i=2;i<NPOINTS;i++)
  for (j=0;j<NCAN;j++,ind++)
    mas[ind]=Bz*mas[ind]-A3*mas[ind-NCAN];

/* фильтрация в обратном направлении */

ind=(long)NPOINTS*NCAN-1-NCAN*2;
for (i=2;i<NPOINTS;i++)
  for (j=0;j<NCAN;j++,ind--)
    mas[ind]=Bz*mas[ind]-A1*mas[ind+NCAN]-A2*mas[ind+NCAN+NCAN];
ind=(long)NPOINTS*NCAN-1-NCAN*2;
for (i=2;i<NPOINTS;i++)
  for (j=0;j<NCAN;j++,ind--)
    mas[ind]=Bz*mas[ind]-A3*mas[ind+NCAN];


}

/*********************************************************/
/* определения фазы текущего фона по отношению к стимулу */
/* фаза приведена к косинусу                             */

void phase_search(float * Sn, float * Cs, int stim)
{
int j;
float res;

for (j=0;j<NCAN;j++)
  {
  res=atan2(Sn[stim*NCAN+j],Cs[stim*NCAN+j]) -  M_PI_2;
  if (res<-M_PI) res+=M_PI*2.0;
  PHASE[j]=res;
  }
}

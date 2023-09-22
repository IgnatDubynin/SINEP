#include <math.h>
#include <complex>

/* 256-точечное быстрое преобразование цурье */

const int N256=256;
const int NPOW=8;
void fft256(int sign,double T, double_complex * X)
{
int i,j,ii,ij,nn,nw,nw1,mm,layer,loc,ll;
double delta,zz,w;
double_complex cxcs,hold,xa;
int nmax=N256;
int msk[10];

X=X-1;
zz=M_PI*2.0*sign/(double)nmax;
delta=(sign<=0) ? T : 1.0/(T*(double)nmax);

msk[1]=nmax/2;
for (i=2;i<=NPOW;i++) msk[i]=msk[i-1]/2;

nn=nmax;   mm=2;

/* внешний цикл Аля слоев NPOW */
for (layer=1;layer<=NPOW;layer++)
  {
  nn=nn/2;
  nw=0;
  for (i=1;i<=mm;i+=2)
    {
    ii=nn*i;
    w=(double)nw*zz;
    cxcs=double_complex(cos((double)w),sin((double)w));
    /* вычисление элементов Аля обеих половин кажАого блока */
    for (j=1;j<=nn;j++)
      {
      ij=++ii-nn;
      xa=cxcs*X[ii];
      X[ii]=X[ij]-xa;
      X[ij]=X[ij]+xa;
      }
    /* вычисление обратных аАресов */
    loc=2;
    do
      {
      ll=nw-msk[loc];
      if (ll>0) { nw=ll; loc++; }
      }
    while
      ((loc<=NPOW) && (ll>0));

    if (ll!=0)  nw=msk[loc]+nw;
    else nw=msk[loc+1];

    }
  mm*=2;
  }
nw=0;
for (i=1;i<=nmax;i++)
  {
  nw1=nw+1;
  hold=X[nw1];
  if ((nw1-i)>0)  X[nw1]=X[i]*delta;
  if ((nw1-i)>=0) X[i]=hold*delta;
  /* вычисление обратного аАреса */
  for (loc=1;loc<=NPOW;loc++)
    {
    ll=nw-msk[loc];
    if (ll>0) nw=ll;
    if (ll==0) { nw=msk[loc+1];  goto end; }
    if (ll<0)  { nw=msk[loc]+nw; goto end; }
    }
end:
 ;
  }

}


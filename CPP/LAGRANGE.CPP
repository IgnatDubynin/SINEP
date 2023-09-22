
int x0,x1,x2,x3,nep,nd2;
long double * div, * upper;

void ini_Lagrange(void)
{
int i,j;
int * xx;
int N;

/* границы области интерполяции */
x1=FIRST_EP;    x2=LAST_EP;
nep=x2-x1+1;
//x0=x1-nep;
x0=x1-10;
//if (x0>LABEL/2) x0=LABEL/2;
if (x0<20) x0=20;
nd2=x1-x0;
x3=x2+nd2;
if (x3>NPOINTS-1)
  {
  nd2=NPOINTS-1-x2;
  x0=x1-nd2;
  x3=x2+nd2;
  }
N=nd2*2;
div=(long double*)farmalloc(N*sizeof(long double));
upper=(long double*)farmalloc(nep*sizeof(long double));
xx=(int*)farmalloc(N*sizeof(int));

for (i=0;i<nd2;i++)     xx[i]=i;
for (i=nd2;i<N;i++) xx[i]=i+nep;


for (i=0;i<nep;i++)
  {
  upper[i]=1.0;
  for (j=0;j<N;j++)
    upper[i]*=(long double)(nd2+i-xx[j]);
  }


for (i=0;i<N;i++)
  {
  div[i]=1.0;
  for (j=0;j<N;j++)
    if (i!=j)
      div[i]*=(long double)(xx[i]-xx[j]);
  }

farfree(xx);

}

void Lagrange(float* inp,float* out)
{
int i,j,k;
long double r,rr;

for (i=0;i<x1;i++)         out[i]=inp[i];
for (i=x2+1;i<NPOINTS;i++) out[i]=inp[i];

for (i=x1;i<=x2;i++)
  {
  r=0;
  for (j=0,k=x0;j<nd2;j++,k++)
    {
    rr=upper[i-x1]/div[j]/(i-k);
    r+=upper[i-x1]/div[j]*inp[k]/(i-k);
    }
  for (j=nd2,k=x2+1;j<nd2*2;j++,k++)
    {
    rr=upper[i-x1]/div[j]/(i-k);
    r+=upper[i-x1]/div[j]*inp[k]/(i-k);
    }
  out[i]=r;
  }

}



void del_Lagrange(void)
{
farfree(div);
farfree(upper);
}


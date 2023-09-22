#include "sinepvar.h"
//#include <alloc.h>
#include <stdio.h>
#include "get_dat.h"
//#include <complex>
#include <string.h>
#include "fft256.cpp"
#include "sinepclc.h"

FILE * f_dat;

void  open_tmp_file(char * name);
void  mail_c(double * mas,double_complex * buf,int n);
void  mail_f(double * mas,double_complex * buf,int n);
void  spectr(double_complex * buf,double * buf2);
void  aver_spectr(double * buf2, double * smas,int j);
void  aver(double * mas,double * amas);
void  transf_func(int j,double * buf2, double * smas,int cnt_aver);
int   read_tmp_file(char * name,double_complex * bufc);


void  fourie_filtr(char * name_inf, int have_all, int have_pause)
{
unsigned int j,i;
int k;
long ind;
int cnt_aver=0;
char name[65];
double minr,maxr;
double * mas, * amas, * smas, * buf2, * buff;
double_complex * buf, * bufc;

ini_dat(name_inf);
/*
mas=(double*)farmalloc((long)NCAN*NPOINTS*sizeof(double));
amas=(double*)farmalloc((long)NCAN*NPOINTS*sizeof(double));
buff=(double*)farmalloc((long)NPOINTS*sizeof(double));
smas=(double*)farmalloc((long)NCAN*256*sizeof(double));
buf=(double_complex*)farmalloc((long)256*sizeof(double_complex));
buf2=(double*)farmalloc((long)256*sizeof(double));
*/
mas = new double  [NCAN*NPOINTS + NCAN];
amas= new double  [NCAN*NPOINTS + NCAN];
buff= new double  [NPOINTS + 1];
smas= new double  [NCAN*256 + NCAN];
buf = new double_complex[256+1];
buf2= new double  [256 + 1];

for (j=0;j<NCAN*256;j++)     smas[j]=0.0;
for (j=0;j<NCAN*NPOINTS;j++) amas[j]=0.0;

while (next_dat_fltr(mas,name))
  {
  open_tmp_file(name);
  for (j=0;j<NCAN;j++)
    {
    mail_c(mas,buf,j);
    fft256(-1,NMS, buf);
    fwrite(buf,sizeof(double_complex),256,f_dat);

    spectr(buf,buf2);
  
    aver_spectr(buf2,smas,j);
    }

  fclose(f_dat);
  aver(mas,amas);
  cnt_aver++;
  }


/* вычисление переАаточной ф-ции и собственно фильтрация */


if (cnt_aver>0)
  {
  for (i=0,ind=0;i<NPOINTS;i++)
    for (j=0;j<NCAN;j++,ind++)
      amas[ind]/=cnt_aver;
  for (i=0,ind=0;i<256;i++)
    for (j=0;j<NCAN;j++,ind++)
      smas[ind]/=cnt_aver;
  }

/* спектр усреАненного Вх */
for (j=0;j<NCAN;j++)
  {
  mail_c(amas,buf,j);
  fft256(-1,NMS, buf);
  spectr(buf,buf2);
  transf_func(j,buf2,smas,cnt_aver);
  }

/* зануление отрицательных значений */

for(j=0,ind=0;j<NCAN;j++)
  for (i=0;i<256;i++,ind++)
    if (smas[ind]<0.0) smas[ind]=0.0;

//getch();
//farfree(buf2);
delete buf2;

if (have_all==0) goto end;

//bufc=(double_complex*)farmalloc((long)256*NCAN*sizeof(double_complex));
bufc= new double_complex [256*NCAN + NCAN];
if (buf2==NULL) perror("юе хватает памяти");




ini_dat(name_inf);
min_max(amas,&minr,&maxr);

while (next_dat_fltr(mas,name))
  {

  if (read_tmp_file(name,bufc))
    {
    for (j=0;j<NCAN;j++)
      {
      for (i=0,ind=(long)j*256;i<256;i++,ind++) buf[i]=bufc[ind];
      /* умножение на переАаточную функцию */
      for (i=0;i<256;i++) buf[i]*=smas[(long)i*NCAN+j];
      fft256(1,NMS, buf);
      mail_f(mas,buf,j);

      }
    /* запись в выхоАной файл и уАаление временного файла */
    change_ext(name,".sep");
    write_result(name,mas);
    change_ext(name,".tmp");
    remove(name); 
    }

  }

/*
*/
delete bufc;
end:
delete buf;
delete smas;
delete buff;
delete amas;
delete mas;
}


void  open_tmp_file(char * name)
{
char name_dat[128];

//fnsplit(name,drive,dir,file,ext);
//fnmerge(name_dat,drive,dir,file,".tmp");
strcpy(name_dat,name);
change_ext(name_dat,".tmp");
if ((f_dat=fopen(name_dat,"wb"))==NULL) exit(1);
}

void  mail_c(double * mas,double_complex * buf,int n)
{
int first,last,n_zerou,i,cnt=0;
if (NPOINTS>=256)
  {
  i=(NPOINTS-256)/2;
  first=i; last=i+256; n_zerou=0;
  }
else
  {
  first=0; last=NPOINTS; n_zerou=256-NPOINTS;
  }
for (i=first;i<last;i++,cnt++) buf[cnt]=double_complex(mas[i*NCAN+n],0.0);
for (i=0;i<n_zerou;i++,cnt++)  buf[cnt]=double_complex(0.0,0.0);

}

void  mail_f(double * mas,double_complex * buf,int n)
{
int first,i,cnt=0;
if (NPOINTS>=256)
  {
  i=(NPOINTS-256)/2;
  first=i;
  for (cnt=0;cnt<first;cnt++)    mas[NCAN*cnt+n]=0.0;
  for (i=0;i<256;  i++,cnt++)    mas[NCAN*cnt+n]=real(buf[i]);
  for (;cnt<NPOINTS;cnt++)       mas[NCAN*cnt+n]=0.0;
  }
else
  for (i=0;i<NPOINTS;i++)      mas[NCAN*i+n]=real(buf[i]);

}


void  spectr(double_complex * buf,double * buf2)
{
int i;
for (i=0;i<256;i++)
  buf2[i]=(double)(real(buf[i])*real(buf[i])+imag(buf[i])*imag(buf[i]));
}

void aver_spectr(double * buf2, double * smas,int j)
{
unsigned int i;
for (i=0;i<256;i++)  smas[(long)i*NCAN+j]+=buf2[i];
}

void  aver(double * mas,double * amas)
{
int i,j;
long ind=0;
for (i=0;i<NPOINTS;i++)
  for (j=0;j<NCAN;j++,ind++)
    amas[ind]+=mas[ind];
}

void  transf_func(int j,double * buf2, double * smas,int cnt_aver)
{
unsigned int i;
long k;
double r;
if (cnt_aver>1)
  for (i=0;i<256;i++)
    {
    k=(long)i*NCAN+j;
    r=smas[k];
    if (r>1.0e-20)
      smas[k]= (buf2[i]*cnt_aver/(cnt_aver-1)-r/(cnt_aver-1))/r;
    else
      smas[k]=0.0;
    }
else
  for (i=0;i<256;i++)
    {
    k=(long)i*NCAN+j;
    smas[k]=1.0;
    }
}

int   read_tmp_file(char * name,double_complex * bufc)
{
char name_dat[65];
FILE * f_dat;

//fnsplit(name,drive,dir,file,ext);
//fnmerge(name_dat,drive,dir,file,".tmp");
change_ext(name,".tmp");
strcpy(name_dat,name);
if ((f_dat=fopen(name_dat,"rb"))==NULL) return 0;
if (fread(bufc,sizeof(double_complex)*NCAN,256,f_dat)!=256) return 0;
fclose(f_dat);
return 1;

}

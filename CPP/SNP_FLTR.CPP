#include "sinepvar.h"
#include <conio.h>
//#include <alloc.h>
#include <stdio.h>
#include <dir.h>
#include <process.h>
#include "get_dat.h"
#include "snp_drow.h"
#include <complex.h>
#include <string.h>
#include "fft256.cpp"
#include "sinepclc.h"

FILE * f_dat;

void  open_tmp_file(char * name);
void  mail_c(float * mas,complex * buf,int n);
void  mail_f(float * mas,complex * buf,int n);
void  spectr(complex * buf,float * buf2);
void  aver_spectr(float * buf2, float * smas,int j);
void  aver(float * mas,float * amas);
void  transf_func(int j,float * buf2, float * smas,int cnt_aver);
int   read_tmp_file(char * name,complex * bufc);


void  fourie_filtr(char * name_inf, int have_all, int have_pause)
{
unsigned int j,i;
int k;
long ind;
int cnt_aver=0;
char name[65];
float minr,maxr;
float * mas, * amas, * smas, * buf2, * buff;
complex * buf, * bufc;
if (LEVEL_SHOW) drow_frame("FOURIE Filtration pass 1","");
else
  if (have_all)
    {
    printf("%s\n","ВЫДЕЛЕНИЕ ЕДИНИЧНЫХ ВП МЕТОДОМ ФУРЬЕ ФИЛЬТРАЦИИ");
    printf("%s\n","ПРОХОД 1");
    }
if (LEVEL_SHOW) ini_fourie_filtr();
ini_dat(name_inf);
/*
mas=(float*)farmalloc((long)NCAN*NPOINTS*sizeof(float));
amas=(float*)farmalloc((long)NCAN*NPOINTS*sizeof(float));
buff=(float*)farmalloc((long)NPOINTS*sizeof(float));
smas=(float*)farmalloc((long)NCAN*256*sizeof(float));
buf=(complex*)farmalloc((long)256*sizeof(complex));
buf2=(float*)farmalloc((long)256*sizeof(float));
*/
mas = new float  [NCAN*NPOINTS];
amas= new float  [NCAN*NPOINTS];
buff= new float  [NPOINTS];
smas= new float  [NCAN*256];
buf = new complex[256];
buf2= new float  [256];

for (j=0;j<NCAN*256;j++)     smas[j]=0.0;
for (j=0;j<NCAN*NPOINTS;j++) amas[j]=0.0;

while (next_dat_fltr(mas,name))
  {
  if (LEVEL_SHOW) drow_name_file(name);
  else printf("%s\n",name);
  open_tmp_file(name);
  for (j=0;j<NCAN;j++)
    {
    mail_c(mas,buf,j);
    fft256(-1,NMS, buf);
    fwrite(buf,sizeof(complex),256,f_dat);

    spectr(buf,buf2);
    if (LEVEL_SHOW)
      out_spectr(buf2,j,cnt_aver);
    aver_spectr(buf2,smas,j);
    }

  fclose(f_dat);
  aver(mas,amas);
  cnt_aver++;
  }

if ((LEVEL_SHOW) && (have_pause))
  {
  i=getch();
  if (i==0) i=getch();
  if (i==27) exit(0);
  }

/* вычисление передаточной ф-ции и собственно фильтрация */


if (cnt_aver>0)
  {
  for (i=0,ind=0;i<NPOINTS;i++)
    for (j=0;j<NCAN;j++,ind++)
      amas[ind]/=cnt_aver;
  for (i=0,ind=0;i<256;i++)
    for (j=0;j<NCAN;j++,ind++)
      smas[ind]/=cnt_aver;
  }

/* спектр усредненного ВП */
for (j=0;j<NCAN;j++)
  {
  mail_c(amas,buf,j);
  fft256(-1,NMS, buf);
  spectr(buf,buf2);
  transf_func(j,buf2,smas,cnt_aver);
  if (LEVEL_SHOW) out_spectr(buf2,j,-1);
  }
if ((LEVEL_SHOW) && (have_pause))
  {
  i=getch();
  if (i==0) i=getch();
  if (i==27) exit(0);
  }
/* зануление отрицательных значений */

for(j=0,ind=0;j<NCAN;j++)
  for (i=0;i<256;i++,ind++)
    if (smas[ind]<0.0) smas[ind]=0.0;

//getch();
//farfree(buf2);
delete buf2;

if (have_all==0) goto end;

//bufc=(complex*)farmalloc((long)256*NCAN*sizeof(complex));
bufc= new complex [256*NCAN];
if (buf2==NULL) error("Не хватает памяти");


if (LEVEL_SHOW) drow_frame("FOURIE Filtration pass 2","");
else printf("%s\n","ПРОХОД 2");


ini_dat(name_inf);
min_max(amas,&minr,&maxr);

while (next_dat_fltr(mas,name))
  {

  if (read_tmp_file(name,bufc))
    {
    if (LEVEL_SHOW)
      {
      drow_frame("FOURIE Filtration pass 2",name);
      drow_subhead(0,"Averaging EP");
      drow_labels(shift_curve);
      drow_ampl(shift_curve,minr,maxr,0);
      for (j=0;j<NCAN;j++)
	{
	for (i=0,ind=j;i<NPOINTS;i++,ind+=NCAN)
	buff[i]=amas[ind];
	out_curve(buff,5,h_ch*5+shift_curve/2+shift_curve*j,NPOINTS,h_graph,minr,maxr);
	}
      drow_subhead(1,"Single EP");
      }
    else printf("%s\n",name);
    for (j=0;j<NCAN;j++)
      {
      for (i=0,ind=(long)j*256;i<256;i++,ind++) buf[i]=bufc[ind];
      /* умножение на передаточную функцию */
      for (i=0;i<256;i++) buf[i]*=smas[(long)i*NCAN+j];
      fft256(1,NMS, buf);
      mail_f(mas,buf,j);

      if (LEVEL_SHOW)
	{
	for (i=0;i<256;i++)
	buff[i]=real(buf[i])*2.0;
	k=(NPOINTS-256)/2;
	if (k<0) k=0;
	out_curve(buff,335+k,h_ch*5+shift_curve/2+shift_curve*j,256,h_graph,minr,maxr);
	}


      }
    /* запись в выходной файл и удаление временного файла */
    change_ext(name,".sep");
    write_result(name,mas);
    change_ext(name,".tmp");
    remove(name);

   if ((LEVEL_SHOW) && (have_pause))
     {
     i=getch();
     if (i==0) i=getch();
     if (i==27) exit(0);
     }
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
char name_dat[65];
char drive[MAXDRIVE];
char dir[MAXDIR];
char file[MAXFILE];
char ext[MAXEXT];

fnsplit(name,drive,dir,file,ext);
fnmerge(name_dat,drive,dir,file,".tmp");
if ((f_dat=fopen(name_dat,"wb"))==NULL) exit(1);
}

void  mail_c(float * mas,complex * buf,int n)
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
for (i=first;i<last;i++,cnt++) buf[cnt]=complex(mas[i*NCAN+n],0.0);
for (i=0;i<n_zerou;i++,cnt++)  buf[cnt]=complex(0.0,0.0);

}

void  mail_f(float * mas,complex * buf,int n)
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


void  spectr(complex * buf,float * buf2)
{
int i;
for (i=0;i<256;i++)
  buf2[i]=(float)(real(buf[i])*real(buf[i])+imag(buf[i])*imag(buf[i]));
}

void aver_spectr(float * buf2, float * smas,int j)
{
unsigned int i;
for (i=0;i<256;i++)  smas[(long)i*NCAN+j]+=buf2[i];
}

void  aver(float * mas,float * amas)
{
int i,j;
long ind=0;
for (i=0;i<NPOINTS;i++)
  for (j=0;j<NCAN;j++,ind++)
    amas[ind]+=mas[ind];
}

void  transf_func(int j,float * buf2, float * smas,int cnt_aver)
{
unsigned int i;
long k;
float r;
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

int   read_tmp_file(char * name,complex * bufc)
{
char name_dat[65];
char drive[MAXDRIVE];
char dir[MAXDIR];
char file[MAXFILE];
char ext[MAXEXT];
FILE * f_dat;

fnsplit(name,drive,dir,file,ext);
fnmerge(name_dat,drive,dir,file,".tmp");
if ((f_dat=fopen(name_dat,"rb"))==NULL) return 0;
if (fread(bufc,sizeof(complex)*NCAN,256,f_dat)!=256) return 0;
fclose(f_dat);
return 1;

}

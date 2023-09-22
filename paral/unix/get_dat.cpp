#include <stdio.h>
#include "sinepvar.h"
#include <string.h>
#include <dirent.h>


FILE* f_inf=NULL;
char ext[4]="";

int ini_dat(char* name)
{
#define ERROR printf(SER,s2,name),getchar(),exit(1)
const char SER[]="\nOшибка в параметре %s файла %s";
//const char SER2[]="\nВ информационном файле %s не заАан параметр %s";
char str[129];
char * s2;
if (f_inf==NULL)
{
if ((f_inf=fopen(name,"rt"))!=NULL)
  {
  while (feof(f_inf)==0)
    {
    fscanf(f_inf,"%s",str);
    if (strlen(str)>1)
      if (s2=strchr(str,'='),s2!=NULL)
	{
	s2[0]='\x0';
	s2+=1;
	if ((strcmp(str,"NCAN")==0) || (strcmp(str,"nkan")==0))
	  if (sscanf(s2,"%i",&NCAN)!=1)     ERROR;
	if (strcmp(str,"NPOINTS")==0)
	  if (sscanf(s2,"%i",&NPOINTS)!=1)  ERROR;
	if (strcmp(str,"NMS")==0)
	  if (sscanf(s2,"%f",&NMS)!=1)      ERROR;
	if (strcmp(str,"LABEL")==0)
	  if (sscanf(s2,"%i",&LABEL)!=1) ERROR;
	if (strcmp(str,"WORD")==0)
	  if (sscanf(s2,"%s",WORD)!=1)      ERROR;

	}
    } // end of   while (feof(f_inf))

//  if (NCAN==0)    printf(SER2,name,"NCAN"),   getch(),exit(1);
//  if (NPOINTS==0) printf(SER2,name,"NPOINTS"),getch(),exit(1);
//  if (NMS==0)     printf(SER2,name,"NMS"),    getch(),exit(1);

  fseek(f_inf,0,SEEK_SET);
  if ((DEMODULATION) && (FOURIE))
    { ext[0]='s'; ext[1]='e'; ext[2]='p'; ext[3]='\x0'; }

  return 1;
  }
else // if ((f_inf=fopen(name,"rt"))!=NULL)
 return 0;
}
else
 {
 fseek(f_inf,0,SEEK_SET);
  if ((DEMODULATION) || (FOURIE))
//    { ext[0]='s'; ext[1]='e'; ext[2]='p'; ext[3]='\x0'; }
        strcpy(ext,"sep");
 return 0;
 }
}


int next_dat(double* mas, char * name)
{
int i,j;
char str[129];
FILE * f_dat=NULL;
short* buf;
char * bufc;
long ind=0;

str[0]='\x0';

while ((f_dat==NULL) && (feof(f_inf)==0))
  {
  fscanf(f_inf,"%s",str);
  if ((strlen(str)>1)          &&
      (strchr(str,'=')==NULL)  &&
      (strchr(str,'.')!=NULL) )
    {
    sscanf(str,"%s",name);
//    strcpy(name,str);
    if ((f_dat=fopen(str,"rb"))!=NULL)
      {
      if (strcmp(FORMAT,"byte")==0)
	{
	bufc=new char [NPOINTS*NCAN + NCAN];
	fread(bufc,NPOINTS,NCAN,f_dat);
	fclose(f_dat);
	for (i=0;i<NPOINTS;i++)
	  for (j=0;j<NCAN;j++,ind++)
	    mas[ind]=(double)(bufc[ind]);
	delete bufc;
	}
      else
	{
	buf= new short [NPOINTS*NCAN + NCAN];
	fread(buf,NPOINTS,NCAN*2,f_dat);
	fclose(f_dat);
	for (i=0;i<NPOINTS;i++)
	  for (j=0;j<NCAN;j++,ind++)
	    mas[ind]=(double)(buf[ind]-2048);
	delete buf;
	}
      }
    }
  }



if (f_dat==NULL) return 0;
else             return 1;

}

void  write_result(char * name,double * mas)
{
char ss[128];
FILE * f_dat;
int i,j;
long ind;

strcpy(ss,name);
//fnsplit(name,drive,dir,file,ext);
//fnmerge(ss,drive,dir,file,".sep");
change_ext(ss,".sep");
f_dat=fopen(ss,"wt");
for (i=0,ind=0;i<NPOINTS;i++,fprintf(f_dat,"\n"))
  for (j=0;j<NCAN;j++,ind++)
    fprintf(f_dat," %7.2f",mas[ind]);
fclose(f_dat);

}

int next_dat_fltr(double* mas,char * name)
{
int i,j;
FILE * f_dat=NULL;
long ind=0;
char str[129];

if (ext[0]=='s')
  {
  str[0]='\x0';

  while ((f_dat==NULL) && (feof(f_inf)==0))
    {
    fscanf(f_inf,"%s",str);
    if ((strlen(str)>1)          &&
	(strchr(str,'=')==NULL)  &&
	(strchr(str,'.')!=NULL) )
      {

//      fnsplit(str,drive,dir,file,ex);
//      fnmerge(str,drive,dir,file,".sep");
       change_ext(str,".sep");

      if ((f_dat=fopen(str,"rt"))!=NULL)
	{
	for (i=0;i<NPOINTS;i++)
	  for (j=0;j<NCAN;j++,ind++)
	      fscanf(f_dat,"%f",&mas[ind]);
	fclose(f_dat);
	strcpy(name,str);
	}
      }
    }



  if (f_dat==NULL) return 0;
  else             return 1;


  }
else
  return next_dat(mas,name);

}



void create_inf(int h)
{
FILE * f_inf;
f_inf=fopen(out_inf,"wt");

    fprintf(f_inf,"%s=%i\n","NCAN",NCAN);
    fprintf(f_inf,"%s=%i\n","NPOINTS",NPOINTS);
    fprintf(f_inf,"%s=%f\n","NMS",NMS);
    fprintf(f_inf,"%s=%i\n","FIRST_EP",FIRST_EP);
    fprintf(f_inf,"%s=%i\n","LAST_EP",LAST_EP);
    fprintf(f_inf,"%s=%i\n","LABEL",LABEL);
    fprintf(f_inf,"%s=%f\n","MIN_FRQ",MIN_FRQ);
    fprintf(f_inf,"%s=%f\n","MAX_FRQ",MAX_FRQ);
    fprintf(f_inf,"%s=%s\n","FORMAT","ascii");
    if (h)
      fprintf(f_inf,"%s=%f\n","ORD_FRQ",ORD_FRQ);
    if (COEF_MUL!=1.0)
      fprintf(f_inf,"%s=%f\n","coef_mul",COEF_MUL);


fclose(f_inf);
}

void add_to_inf(char *name,int h)
{
int i;
FILE * f_inf;
f_inf=fopen(out_inf,"at");
fprintf(f_inf,"%s",name);
/*if (h)
  for (i=0;i<NCAN;i++)
    fprintf(f_inf," %5.2f",PHASE[i]);*/
fprintf(f_inf,"\n");
fclose(f_inf);
}

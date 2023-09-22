#include <string.h>
//#include <alloc.h>
#include <stdarg.h>
//#include <varargs.h>
#include <stdio.h>
#include "param.h"

extern unsigned _floatconvert;
#pragma extref _floatconvert


abscent_value GAP={0x8000,0x8000,0,0};

PARAM * ptr_PARAM=NULL;



PARAM::PARAM(char * str)
  {
  PARAM * p;
//mod=ind=dim=0;
  mod=ind=dim=0;
  strncpy(name,str,Len_name-1);
  next=NULL;
  p=ptr_PARAM;
  if (p==NULL)
    ptr_PARAM=this;
  else
    {
    while ((p->next)!=NULL) p=p->next;
    p->next=this;
    }
  }
void PARAM::clear(void) {};

//PARAM::~PARAM() {};


PARAM_B::PARAM_B(char * str) : PARAM(str) {  mod=ind=dim=0;};
PARAM_B::~PARAM_B() {};

PARAM_I::PARAM_I(char * str) : PARAM(str) {  mod=ind=dim=0;};

PARAM_I::PARAM_I(char * str,int n) : PARAM(str)
  {
  mod=ind=dim=0;
  dim=n;
  v.ptri= new int [n+1];
  v.ptri[0]=0;
  }
void PARAM_I::clear(void)
{
if (ind)
  {
  if (dim!=0)
    delete v.ptri;
  dim=ind=mod='\x0';
  }
}
PARAM_I::~PARAM_I() { if (ind) if (dim) delete v.ptri;}

PARAM_R::PARAM_R(char * str) : PARAM(str) { mod=ind=dim=0;};

PARAM_R::PARAM_R(char * str,int n) : PARAM(str)
  {
  mod=ind=dim=0;
  dim=n;
  v.ptrr=new float [n+1];
  v.ptrr[0]=0.0;
  }
void PARAM_R::clear(void)
{
if (ind)
  {
  if (dim!=0)
    delete v.ptrr;
  dim=ind=mod='\x0';
  }
}
PARAM_R::~PARAM_R() { if (ind) if (dim) delete v.ptrr;}

PARAM_S::PARAM_S(char * str) : PARAM(str) {mod=ind=dim=0;};

PARAM_S::PARAM_S(char * str,int n) : PARAM(str)
  {
  int i;
  mod=ind=dim=0;
  dim=n;
//  i=sizeof(char*);
  v.ss=/*(char**)farmalloc(n*sizeof(char*));*/ new char * [n];
  for (i=0;i<n;i++)
    v.ss[i]=NULL;
  }

void PARAM_S::clear(void)
  {
  int i;
  if (ind)
    {
    if (dim)
      {
      for (i=0;i<dim;i++)
	if (v.ss[i]!=NULL) delete v.ss[i];
      delete v.ss;
      }
    else
      delete v.s;
    }
    mod=ind=dim=0;
  }

PARAM_S::~PARAM_S()
  {
  int i;
  if (ind)
  {
  if (dim)
    {
    for (i=0;i<dim;i++)
      if (v.ss[i]!=NULL) delete v.ss[i];
    delete v.ss;
    }
  else
    delete v.s;
  }
  mod=ind=dim=0;
  }
void PARAM_S::set(char * ss)
{
char *pp= new char [strlen(ss)+1];
v.s=pp;
strcpy(v.s,ss);
}

void PARAM_S::set(int n,char * ss)
{
char *pp= new char [strlen(ss)+1];
strcpy(pp,ss);
v.ss[n-1]=pp;
}


int PARAM_S::PAR_READ(char *ss)
{
int n,k,i;
char * pp;
while ((strlen(ss)>0) && ((unsigned)ss[0]<=(unsigned)' ')) ss++;
while (((n=strlen(ss))>0) && ((unsigned)ss[n-1]<=(unsigned)' ')) ss[n-1]='\x0';
if (dim==0)
   {
   pp=new char [strlen(ss)+1];
   v.s=pp; strcpy(v.s,ss);
   }
else
  {
  i=sscanf(ss,"%i",&k);
  if (i!=1) return 1;
  if ((k<1) || (k>dim)) return 1;
  while ((strlen(ss)>0) && ((unsigned)ss[0]>(unsigned)' ')) ss++;
  while ((strlen(ss)>0) && ((unsigned)ss[0]<=(unsigned)' ')) ss++;
  pp=new char [strlen(ss)+1];
  strcpy(pp,ss);
  v.ss[k-1]=pp;
  }
ind='\x1';
return 0;
}

void PARAM_I::set(int lc,...)
{
int i;
va_list ap;
va_start(ap,lc);
if (dim==0)
  v.i=va_arg(ap,int);
else
  for (i=0;i<dim;i++)
    if ((v.ptri[i]=va_arg(ap,int))!=GAP.INT)
      v.ptri[0]+=1;
    else
      i=dim;
va_end(ap);
}



void PARAM_R::set(int lc,...)
{
int i;
va_list ap;
va_start(ap,lc);
if (dim==0)
  v.r=va_arg(ap,float);
else
  for (i=0;i<dim;i++)
    if ((v.ptrr[i]=va_arg(ap,float))!=GAP.FLT)
      v.ptrr+=1;
    else
      i=dim;
va_end(ap);
}


int PARAM::PAR_READ(char*str)
  {
  int i;
  while ((strlen(str)>0) && ((unsigned)str[0]<=(unsigned)' ')) str++;
  if ((str[0]=='<') || (str[0]=='>')) mod=*str++;
  if (dim==0)
    { if (sscanf(str,"%i",&v.i)) { ind='\x1'; return 0; } else return 1; }
  else
    for (i=1;i<=dim;i++)
      {
      switch (sscanf(str,"%i",&v.ptri[i]))
	{
	case 1  : v.ptri[0]+=1;
		 break;
	case 0  : return 1;
	case EOF: i=dim+1;
	}
      /* пропуск ввеАенного числа */
      while ((strlen(str)>0) && ((unsigned)str[0]>(unsigned)' ')) str++;
      while ((strlen(str)>0) && ((unsigned)str[0]<=(unsigned)' ')) str++;
      }
  ind='\x1';
  return 0;
  }

int PARAM_R::PAR_READ(char*str)
  {
  int i;
  while ((strlen(str)>0) && ((unsigned)str[0]<=(unsigned)' ')) str++;
  if ((str[0]=='<') || (str[0]=='>')) mod=*str++;
  if (dim==0)
    { if (sscanf(str,"%f",&v.r)) { ind='\x1'; return 0; } else return 1; }
  else
    for (i=1;i<=dim;i++)
      {
      switch (sscanf(str,"%f",&v.ptrr[i]))
	{
	case 1  : v.ptri[0]+=1;
		 break;
	case 0  : return 1;
	case EOF: i=dim+1;
	}
      /* пропуск ввеАенного числа */
      while ((strlen(str)>0) && ((unsigned)str[0]>(unsigned)' ')) str++;
      while ((strlen(str)>0) && ((unsigned)str[0]<=(unsigned)' ')) str++;
      }
  ind='\x1';
  return 0;
  }

int PARAM_B::PAR_READ(char*str)
  {
  int n;
  while ((strlen(str)>0) && ((unsigned)str[0]<=(unsigned)' ')) str++;
  while (((n=strlen(str))>0) && ((unsigned)str[n-1]<=(unsigned)' ')) str[n-1]='\x0';
  if ((strcmp(str,"yes")==0) ||
      (strcmp(str,"y")  ==0) ||
      (strcmp (str,"Аа") ==0) ||
      (strcmp (str,"фа") ==0) ||
      (strcmp (str,"фA") ==0) ||
      (strcmp(str,"on") ==0) )
    v.i=1;
  else
    v.i=0;

  ind='\x1';
  return 0;
  }





char * ER[]={"",
	     "файл не найАен",
	     "параметры не заАаны"
	     };
char ERp[129]="неверное значение параметра\n";

int read_par(char * str);

int read_PARAM(char * name,char** error)
{
const int N257=257;  // максимальная Алина строки
FILE * par;
ERp[28]='\x0';
*error=ER[0];
char str[N257];
int i,k;

if (ptr_PARAM==NULL) { *error=ER[2]; return -1; }
par=fopen(name,"rt");
if (par==NULL)       { *error=ER[1]; return -2; }
do
  {
  fgets(str,N257-1,par);
  k=read_par(str);
  }
while( (feof(par)==0) && (k==0));

fclose(par);

if (k==0) return 0;
else
  {
  int n=strlen(str)+28;
  if (n>128) n=128;
  for (i=28;i<n;i++) ERp[i]=str[i-28];
  ERp[i]='\x0';
  *error=ERp;
  return -3;
  }

}


int read_par(char * str)
{
char * ptr,*p;
PARAM * PAR;
int n;
/* уАаление комментариев */
ptr=strchr(str,';');
if (ptr!=NULL) while (*ptr>'\x0')
		 if (*(ptr++)==';') *(--ptr)='\x0';
/* уАаление начальных пробелов */
ptr=str;
while (((unsigned)*ptr>(unsigned)'\x0') && ((unsigned)*ptr<=(unsigned)' ')) ptr++;
if (strlen(ptr)==0) return 0;
p=strchr(ptr,'=');
if (p==NULL) return 0; // если не параметр
else *p++='\x0';
//strlwr(ptr);
while ((ptr[0]!='\x0') && ((unsigned)ptr[strlen(ptr)-1]<=(unsigned)' '))
  ptr[strlen(ptr)-1]='\x0';  // уАаление конечных пробелов в имени параметра

n=0;
PAR=ptr_PARAM;
/* Аля совместимости со старой версией */
if (strcmp("nkan",ptr)==0)    ptr[1]='c';
if (strcmp("typeval",ptr)==0) strcpy(ptr,"format");

while (PAR!=NULL)
 if (strcmp(PAR->name,ptr))
   PAR=PAR->next;
 else
   {
   PAR->clear();
   n=PAR->PAR_READ(p);
   PAR=NULL;
   }

return n;
}


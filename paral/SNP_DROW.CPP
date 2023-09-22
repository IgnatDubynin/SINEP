#include <graphics.h>
#include "sinepvar.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <process.h>

int curve_color=EGA_WHITE,
    curve_style=SOLID_LINE,
    curve_thick=1,
    curve_paceox=1;

void out_curve(
		 float* buf,
		 int x0,int y0,int n,int h,
		 float minr,float maxr)
{
int i;
int x,y,y_old;
float r=maxr-minr;

setcolor(curve_color);
setlinestyle(curve_style,0,curve_thick);
x=x0;
y_old=y0-(int)((buf[0]-minr)*h/r);
for (i=1;i<n;i++)
  {
  y=y0-(int)((buf[i]-minr)*h/r);
  line(x,y_old,x+curve_paceox,y);
  x+=curve_paceox;
  y_old=y;
  }

}

void drow_name_file(char * name)
{
const char SF[]="File ";
char str[65];
int x;
int y=h_ch>>1;
x=maxx-2-textwidth(SF)-textwidth(name);
setfillstyle(SOLID_FILL,EGA_BLACK);
bar(x-textwidth(" "),y-1,maxx-1,y+h_ch+1);
strcpy(str,SF);
strcat(str,name);
moveto(x,y);
setcolor(EGA_LIGHTGRAY);
if (name[0]) outtext(SF);
setcolor(EGA_WHITE);
outtext(name);
}

void drow_frame(char * part, char * name)
{
const char SN[]="Single Evoked Potentials";
char str[129];
int y=h_ch>>1;
cleardevice();
setcolor(EGA_BLUE);
setlinestyle(SOLID_LINE,0,1);
rectangle(0,0,maxx,maxy);
line(1,h_ch*2+1,maxx-1,h_ch*2+1);
setfillstyle(SOLID_FILL,EGA_DARKGRAY);
bar(1,1,maxx-1,h_ch*2);
setcolor(EGA_LIGHTBLUE);
outtextxy(5,y,SN);
setcolor(EGA_LIGHTGREEN);
moveto((maxx-textwidth(part))/2,y+1);
outtext(part);
drow_name_file(name);
}

void drow_subhead(int n,char * text)
{
int x,y;
setcolor(EGA_CYAN);
y=h_ch*3;
x=(maxx/2-textwidth(text))/2;
x+=maxx/2*n;
setfillstyle(SOLID_FILL,EGA_BLACK);
bar(2+maxx/2*n,y,maxx/2*(n+1)-5,y+h_ch);
outtextxy(x,y,text);

}
void drow_labels(int shifty)
{
const char SN[]="Canal N ";
int i;
int x=5,y;
char s3[3];

int j,k,col;
int xx,y1,y2,x2;
int n_value;
float lr,rr,r,r0,r00,pacer_ox;
char ss[3]="ms";
char str[6],s_format[6]="%3.0f";
char *pss;


y=h_ch*4+5;

for (i=0;i<NCAN;i++,y+=shifty)
  {
  moveto(x,y);
  setcolor(EGA_LIGHTGRAY);
  outtext(SN);
  sprintf(s3,"%i",i+1);
  setcolor(EGA_CYAN);
  outtext(s3);
  }

/* вывод шкалы времени */

y2=maxy-4;
y1=y2-h_ch*2;
x2=x+curve_paceox*NPOINTS;
pacer_ox=NMS*(NPOINTS-1)/(x2-x);

setfillstyle(SOLID_FILL,EGA_DARKGRAY);
bar(5,maxy-h_ch*2-6,maxx/2-2,maxy-2);

setcolor(EGA_LIGHTGRAY);
col=EGA_LIGHTGRAY;

pss=ss+1;
lr=NMS*(NPOINTS-1);
if (lr<400)
  pss--;
else
 if (lr<10000)  s_format[3]='1';

n_value=3;

rr=lr/n_value;

if (rr<10)
  rr=(float)((int)(rr+0.99));
else
if (rr<100)
  rr=(float)((int)(rr/10+0.99)*10);
else
if (rr<1000)
  rr=(float)((int)(rr/100+0.99)*100);
else
if (rr<10000)
  rr=(float)((int)(rr/1000+0.99)*1000);
else
  rr=(float)((int)(rr/10000+0.99)*10000);

  for (j=0;j<n_value;j++)
    {
    r=rr*j;
    if (pss==ss)
      sprintf(str,s_format,r);
    else
      sprintf(str,s_format,r/1000.0);
    x=5+LABEL*curve_paceox+(int)(r/pacer_ox);
    if (x<x2+1)  line(x,y1+1,x,y2-h_ch-2);
    r=rr/10;
    for (i=1;i<10;i++) // вывод меток без чисел на шкале
      {
      xx=x+r*i/pacer_ox;
      if (xx<=x2)
	{
	if (i % 5 == 0)
	  line(xx,y1+2,xx,y2-h_ch-4);
	else
	  putpixel(xx,y1+3,col);
	}
      }
    if (x+1+textwidth(str)+textwidth(pss)<x2+1)
      outtextxy(x+1,y2-h_ch,str);
    }
  outtextxy(x2+1-textwidth(pss),y2-h_ch,pss);

  for (j=-1;j>-n_value+1;j--)
    {
    r=rr*j;
    if (pss==ss)
      sprintf(str,s_format,r);
    else
      sprintf(str,s_format,r/1000.0);
    x=5+LABEL*curve_paceox+(int)(r/pacer_ox);
    if (x>=5)  line(x,y1+1,x,y2-h_ch-2);
    r=rr/10;
    for (i=1;i<10;i++) // вывод меток без чисел на шкале
      {
      xx=x+r*i/pacer_ox;
      if (xx>=5)
	{
	if (i % 5 == 0)
	  line(xx,y1+2,xx,y2-h_ch-4);
	else
	  putpixel(xx,y1+3,col);
	}
      }
    if (x+1+textwidth(str)+textwidth(pss)<x2+1)
      outtextxy(x+1,y2-h_ch,str);
    }

x=5+LABEL*curve_paceox;
r=0;
setlinestyle(DOTTED_LINE,0,1);
setcolor(EGA_DARKGRAY);
for (i=LABEL;i<=LAST_EP+1;i++)
  {
  if (r<pacer_ox) line(x,h_ch*5,x,maxy-h_ch*4);
  x+=curve_paceox;
  r+=pacer_ox;
  if (r>=100.0) r=r-100.0;
  }

}

void drow_ampl(int shifty,float minr,float maxr,int n)
{
float max,min;
char s1[13],s2[13];
int i,x,y0,y1,y2;
max=maxr/2048*5120/CA*1000;
min=minr/2048*5120/CA*1000;
sprintf(s1,"%.0f",max);
sprintf(s2,"%.0f",min);

x=5+curve_paceox*NPOINTS+maxx/2*n+1;

setlinestyle(SOLID_LINE,0,1);
setcolor(EGA_LIGHTGRAY);

for (i=0;i<NCAN;i++)
  {
  y1=h_ch*5+shifty*i+shifty/2-h_graph;
  y2=y1+h_graph;
  line(x,y1,x,y2);
  line(x,y1,x+2,y1);
  line(x,y2,x+2,y2);
  outtextxy(x+3,y1,s1);
  outtextxy(x+3,y2,s2);
  if ((max>0) && (min<0))
    {
    y0=y1+(int)((float)h_graph*maxr/(maxr-minr));
    line(x,y0,x+2,y0);
    if ((y0-y1>=h_ch) && (y2-y0>=h_ch))
      outtextxy(x+4,y0-h_ch/2,"mkV");
    else
      outtextxy(x+10,(y2-h_ch/2-y1)/2,"mkV");
    }
  }

}

void subscript(int blue,int yellow)
{
int x=maxx*5/8;
int y=maxy-h_ch*3/2;

moveto(x,y);
setcolor(blue);
setlinestyle(0,0,3);
linerel(20,0);
setcolor(EGA_LIGHTGRAY);
outtext(" Row Data  ");

setcolor(yellow);
setlinestyle(0,0,1);
linerel(20,0);
setcolor(EGA_LIGHTGRAY);
outtext(" Single EP");

}


/*******************************************************************/
/**************** рисование при Фурье фильтрации *******************/
/*******************************************************************/

int len_ox=128,len_oy,cnt_n=0;
const int pace_ox_sp=4;

int pace_oz=4;
int max_n;
struct POINT {
       int x,y;
	     } coor[16];

void ini_fourie_filtr(void)
{
const char SS[]="PW functions (Hz):";
int i,x,y;
int h;
float freq;
char ss[21];
extern void out_name_can(int n,int x,int y,float tg);

freq=(1.0/(256*NMS/1000.0))*128/pace_ox_sp;
sprintf(ss,"%5.1f",freq);

i=(NCAN+1)/2;
h=(maxy-h_ch*2-2)/i;
len_oy=h/2;
max_n=h/2/pace_oz-2;
y=3;
for (i=0;i<NCAN;i++)
  {
  if (i%2==0) { x=50; y+=h; }
  else          x=370;
  coor[i].x=x;
  coor[i].y=y;
  }
setcolor(EGA_LIGHTMAGENTA);
setlinestyle(SOLID_LINE,0,3);
for (i=0;i<NCAN;i++)
  {
  line(coor[i].x,coor[i].y,coor[i].x+len_ox,coor[i].y);
  line(coor[i].x,coor[i].y,coor[i].x,coor[i].y-len_oy);
  line(coor[i].x,coor[i].y,coor[i].x+len_ox,coor[i].y-max_n*pace_oz);
  out_name_can(i,coor[i].x,coor[i].y-len_oy,(float)max_n*pace_oz/len_ox);
  outtextxy(coor[i].x+len_ox,coor[i].y+2,ss);
  }

setcolor(EGA_LIGHTRED);
outtextxy((maxx-textwidth(SS))/2,h_ch*3,SS);

}

void out_name_can(int n,int x,int y,float tg)
{
char ss[]="Canal 00";
char s2[2]=" ";
int i,j,l;
int xx=x+5,yy;
float yr=y;

l=strlen(ss);
n++;
j=n%10;
if (n>9) { ss[l-2]='1'; ss[l-1]=(char)(j+'0');}
else     { ss[l-2]=(char)(j+'0'); ss[l-1]='\x0'; l--; }

for (i=0;i<l;i++)
  {
  yy=(int)yr+0.5;
  s2[0]=ss[i];
  outtextxy(xx,yy,s2);
  j=textwidth(s2);
  xx+=j-1;
  yr-=(float)j*tg;
  }
}

float coef_y;

void out_spectr(float * buf,int canal,int cnt_aver)
{
int i;
int x,y0,y,base;
int y_min;
int n;
float max;

if ((cnt_aver<1) && (canal==0))
  {
  for (i=0,max=0.0;i<128;i++) if (max<buf[i])  max=buf[i];
  coef_y=max*4.0;
  }
if (cnt_aver>=0)
  {
  setlinestyle(SOLID_LINE,0,1);
  setcolor(7+(cnt_aver/max_n)*2);
  }
else
  {
  setlinestyle(SOLID_LINE,0,3);
  setcolor(15+cnt_aver);
  cnt_aver=-cnt_aver-1;
  }
x=coor[canal].x+128*(cnt_aver%max_n)/max_n+4;
base=coor[canal].y-pace_oz*(cnt_aver%max_n)-2;
y0=base-(int)((exp(buf[0]/coef_y)-1.0)*len_oy);
y_min=base-len_oy;
n=128/pace_ox_sp;
for (i=1;i<n;i++,x+=pace_ox_sp,y0=y)
  {
  if (y0<y_min) y0=y_min;
  y=base-(int)((exp(buf[i]/coef_y)-1.0)*len_oy);
  if (y>=y_min) line(x,y0,x+pace_ox_sp,y);
  }

}

void error(char * str)
{
closegraph();
printf(str);
exit(1);
}

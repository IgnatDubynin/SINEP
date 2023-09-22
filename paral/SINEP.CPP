#include <conio.h>
#include <stdio.h>
#include <graphics.h>
#include "get_dat.h"
#include "sinepvar.h"
#include "sinepclc.h"
#include "snp_drow.h"
#include "snp_fltr.h"
//#include <alloc.h>
#include <process.h>
#include <math.h>
#include "param.h"
/***/
#include <dos.h>

 PARAM_I Ncan("ncan");
 PARAM_I Npoints("npoints");
 PARAM_I Label("label");
 PARAM_R Nms("nms");
 PARAM_I Level("level");
 PARAM_I First_EP("first_ep");
 PARAM_I Last_EP("last_ep");
 PARAM_R Campl("ca");
 PARAM_R Min_FRQ("min_frq");
 PARAM_R Max_FRQ("max_frq");
 PARAM_R Ord_FRQ("ord_frq");
 PARAM_B Demodulation("demodulation");
 PARAM_B Fourie("fourie");
 PARAM_B Pause("pause");
 PARAM_B Spectr("spectr");
 PARAM_I Order("order");
 PARAM_S Format("format");
 PARAM_R Coef_MUL("coef_mul");





float* mas,* buf,* buf1;
float * mas_sin, * mas_cos, * sum;
float * mas_sini,* mas_cosi,* mas1, *mas2;
float minr,maxr;

void revealing(float* mas,int h);



void main(int argc,char * argv[])
{


char name_dat[65];
char name_cfg[65];
int gdriver = DETECT, gmode=0, errorcode;
long ind;
int i,j;
int pass;
char * er;


if (argc>1)
  {
  i=read_PARAM(argv[1],&er);

  if (i!=0)
    { printf("%s",er);getch();exit(0); }
  do name_cfg[i]=argv[0][i];
  while (name_cfg[i++] !='\x0');

  change_ext(name_cfg,".cfg");
  i=read_PARAM(name_cfg,&er);
  if (i!=0)
    { printf("%s",er);getch();exit(0); }
  if (ini_dat(argv[1]))
    {

    NCAN=Ncan.val();
    NPOINTS=Npoints.val();
    LABEL=Label.val();
    NMS=Nms.val();
    LEVEL_SHOW=Level.val();
    FIRST_EP=First_EP.val();
    LAST_EP=Last_EP.val();
    CA=Campl.val();
    MIN_FRQ=Min_FRQ.val();
    MAX_FRQ=Max_FRQ.val();
    ORD_FRQ=Ord_FRQ.val();
    DEMODULATION=Demodulation.val();
    FOURIE=Fourie.val();
    ORDER=Order.val();
    FORMAT=Format.val();
    if (Coef_MUL.ind) COEF_MUL=Coef_MUL.val();

const char ER1[]="¥ § ¤ ­ ¯ à ¬¥âà %s";
    if (Ncan.ind=='\x0')
       { printf(ER1,Ncan.name); getch();exit(0); }
    if (Nms.ind=='\x0')
       { printf(ER1,Nms.name); getch();exit(0); }
    if (Npoints.ind=='\x0')
       { printf(ER1,Npoints.name); getch();exit(0); }


    if (Ord_FRQ.ind) PHASE= new float [NCAN];

    mas     = new float [NPOINTS*NCAN];
    buf     = new float [NPOINTS];
    buf1    = new float [NPOINTS];
    mas_sin = new float [NPOINTS*NCAN];
    mas_cos = new float [NPOINTS*NCAN];
    sum     = new float [NPOINTS*NCAN];
    mas_sini= new float [NPOINTS*NCAN];
    mas_cosi= new float [NPOINTS*NCAN];
    mas1    = new float [NPOINTS*NCAN];
    mas2    = new float [NPOINTS*NCAN];
/*
    mas=(float*)farmalloc((long)NPOINTS*NCAN*sizeof(float));
    buf=(float*)farmalloc((long)NPOINTS*sizeof(float));
    buf1=(float*)farmalloc((long)NPOINTS*sizeof(float));
    mas_sin =(float*)farmalloc((long)NPOINTS*NCAN*sizeof(float));
    mas_cos =(float*)farmalloc((long)NPOINTS*NCAN*sizeof(float));
    sum     =(float*)farmalloc((long)NPOINTS*NCAN*sizeof(float));
    mas_sini=(float*)farmalloc((long)NPOINTS*NCAN*sizeof(float));
    mas_cosi=(float*)farmalloc((long)NPOINTS*NCAN*sizeof(float));
    mas1    =(float*)farmalloc((long)NPOINTS*NCAN*sizeof(float));
    mas2    =(float*)farmalloc((long)NPOINTS*NCAN*sizeof(float));
*/




    i=0;
    do out_inf[i]=argv[1][i];
    while
       (out_inf[i++]!='\x0');

    change_ext(out_inf,".fin");
    printf("ˆ¬ï ¢ëå®¤­®£® ¨­ä®à¬ æ¨®­­®£® ä ©«  [%s] ? ",out_inf);
    gets(name_dat);
    i=0;
    if (name_dat[0]>' ')
      do out_inf[i]=name_dat[i];
      while
	 (out_inf[i++]!='\x0');

    if (LEVEL_SHOW)
      {
      initgraph(&gdriver, &gmode, "");
      errorcode = graphresult();
      if (errorcode != grOk)
      {
	 printf("Graphics error: %s\n", grapherrormsg(errorcode));
	 printf("Press any key to halt:");
	 getch();
	 exit(1);
      }
      maxx=getmaxx();    maxy=getmaxy();
      h_ch=textheight("H");
      shift_curve=(maxy-h_ch*7)/NCAN;
      h_graph=(shift_curve-h_ch*3)/2;
      switch (NCAN) {
	     case 1 : h_graph/=2;
		      break;
	     case 2 : h_graph=h_graph*2/3;
		      break;
		    }
      }
    create_inf(Ord_FRQ.ind);

if (DEMODULATION==0) goto only_fourie;

    ini_Gram(ORDER);
    if (LEVEL_SHOW==0)
      printf("%s\n","‚›„…‹…ˆ… …„ˆˆ—›• ‚ Œ…’Ž„ŽŒ ŠŽŒ‹…Š‘Ž‰ „…ŒŽ„“‹Ÿ–ˆˆ");

    while (next_dat(mas,name_dat))
      {
      min_max(mas,&minr,&maxr);
      if (LEVEL_SHOW==0) printf("%s\n",name_dat);
      if (LEVEL_SHOW>0)
	{
	drow_frame("COMPLEX DEMODULATION",name_dat);
	if (LEVEL_SHOW<3)
	{
	for (j=0;j<NCAN;j++)
	  {
	  for (i=0,ind=j;i<NPOINTS;i++,ind+=NCAN)
	    buf[i]=mas[ind];
	  out_curve(buf,5,h_ch*5+shift_curve/2+shift_curve*j,NPOINTS,h_graph,minr,maxr);
	  }
	drow_subhead(0,"RAW DATA");
	drow_labels(shift_curve);
	drow_ampl(shift_curve,minr,maxr,0);
	}
	else
	{
	drow_subhead(0,"INTERPOLATION");
	drow_labels(shift_curve);
	}
	drow_subhead(1,"PROCESSING ...");
	}

      revealing(mas,Ord_FRQ.ind);
      if (LEVEL_SHOW) drow_ampl(shift_curve,minr,maxr,0);
      write_result(name_dat,sum);
      change_ext(name_dat,".sep");
      add_to_inf(name_dat,Ord_FRQ.ind);
      if ((LEVEL_SHOW) && Pause.val())
	{
	drow_subhead(1,"RESULT");
	i=getch();
	if (i==0) i=getch();
	if (i==27) exit(0);
	}

      } // while (next_dat(mas),name_dat)

    del_Gram();

only_fourie :
/*
    farfree(mas2);
    farfree(mas1);
    farfree(mas_cosi);
    farfree(mas_sini);
    farfree(sum);
    farfree(mas_cos);
    farfree(mas_sin);
    farfree(buf1);
    farfree(buf);
    farfree(mas);
    if (Ord_FRQ.ind) farfree(PHASE);
*/
    delete mas2;
    delete mas1;
    delete mas_cosi;
    delete mas_sini;
    delete sum;
    delete mas_cos;
    delete mas_sin;
    delete buf1;
    delete buf;
    delete mas;
    if (Ord_FRQ.ind) delete PHASE;


    if ((FOURIE) || Spectr.val()) fourie_filtr(argv[1],FOURIE,Pause.val());

    if (LEVEL_SHOW) closegraph();

    }
  else  //  if (ini_dat(argv[1]))
    {
     printf("\nˆ­ä®à¬ æ¨®­­ë© ä ©« %s ­¥ ­ ©¤¥­\n",argv[1]);
     getch();
    }

  }
else // if (argv>1)
  {
  printf("\n‚ë§®¢ ¯à®£à ¬¬ë:\nsinep.exe name.inf\n");
  printf(",£¤¥ name.inf - ¨¬ï ¨­ä®à¬ æ¨®­­®£® ä ©« \n");
  getch();
  }

fcloseall();

 Ncan.clear();
 Npoints.clear();
 Label.clear();
 Nms.clear();
 Level.clear();
 First_EP.clear();
 Last_EP.clear();
 Campl.clear();
 Min_FRQ.clear();
 Max_FRQ.clear();
 Ord_FRQ.clear();
 Demodulation.clear();
 Fourie.clear();
 Pause.clear();
 Spectr.clear();
 Order.clear();
 Format.clear();
 Coef_MUL.clear();


}




void revealing(float* mas,int have_ord)
{
int max_freq,min_freq;
float pace_freq,freq;
int i,j,k;
long ind;
char ss[65];

for (i=0,ind=0;i<NPOINTS;i++)
  for (j=0;j<NCAN;j++,ind++) sum[ind]=0;

pace_freq=1000.0/(NPOINTS*NMS);
min_freq=(int)(MIN_FRQ/pace_freq)+1;
max_freq=(int)(MAX_FRQ/pace_freq+0.9);
if (max_freq>NPOINTS/2) max_freq=NPOINTS/2;

for (k=min_freq,freq=min_freq*pace_freq;k<max_freq;k++,freq+=pace_freq)
  {
  modulation(freq*M_PI*2*NMS/1000.0,mas,mas_cos,mas_sin);
  sprintf(ss,"PROCESSING ...%4.1f Hz",freq);


  if (LEVEL_SHOW) drow_subhead(1,ss);

  /* ä¨«ìâà æ¨ï (âã¤ -®¡à â­®) */
  filtr(mas_cos);
  filtr(mas_sin);

  for (j=0,ind=0;j<NCAN;j++)
    {
    for (i=0,ind=j;i<NPOINTS;i++,ind+=NCAN)
      buf[i]=mas_cos[ind];
//      interpolation(buf,buf);
      Gram(buf,buf1);
    for (i=0,ind=j;i<NPOINTS;i++,ind+=NCAN)
      mas_cosi[ind]=buf1[i];
    }
  for (j=0,ind=0;j<NCAN;j++)
    {
    for (i=0,ind=j;i<NPOINTS;i++,ind+=NCAN)
      buf[i]=mas_sin[ind];
//      interpolation(buf,buf);
      Gram(buf,buf1);
    for (i=0,ind=j;i<NPOINTS;i++,ind+=NCAN)
      mas_sini[ind]=buf1[i];
    }

/****/
  extern void debug(float * m1,float * m2,int shift);
if (LEVEL_SHOW>2) debug(mas_sin,mas_sini,shift_curve);
    if (kbhit()) i=getch();
    if (i==0x1b) exit(1);
/****/
  if ( (have_ord) && (ORD_FRQ<=freq+pace_freq/2) && (ORD_FRQ>freq-pace_freq/2))
    phase_search(mas_sin,mas_cos,LABEL);

  subtracking(mas_cosi,mas_cos,mas_cosi);
  subtracking(mas_sini,mas_sin,mas_sini);

//  demodulation(freq*M_PI*2*NMS/1000.0,mas1,mas_cos,mas_sin);
  demodulation(freq*M_PI*2*NMS/1000.0,mas2,mas_cosi,mas_sini);


//  subtracking(mas2,mas1,mas2);
  summing(sum,sum,mas2);
  summing(sum,sum,mas2);

  if (LEVEL_SHOW>0)
    {
    setfillstyle(SOLID_FILL,EGA_BLACK);
    bar( 335,h_ch*2+5,maxx-1,maxy-h_ch*3-1);
    }
  if (LEVEL_SHOW>2)
  {
  subscript(EGA_LIGHTBLUE,EGA_YELLOW);
  curve_color=EGA_LIGHTBLUE;
  curve_thick=3;
  for (j=0,ind=0;j<4;j++)
    {
    for (i=0,ind=j;i<NPOINTS;i++,ind+=NCAN)
      buf[i]=mas[ind];
    out_curve(buf,335,h_ch*5+shift_curve/2+shift_curve*j,NPOINTS,h_graph,minr,maxr);
    }
  curve_thick=1;
  }
  if (LEVEL_SHOW>0)
    {
    curve_color=EGA_YELLOW;
    for (j=0,ind=0;j<NCAN;j++)
      {
      for (i=0,ind=j;i<NPOINTS;i++,ind+=NCAN)
	buf[i]=sum[ind];
      out_curve(buf,335,h_ch*5+shift_curve/2+shift_curve*j,NPOINTS,h_graph,minr,maxr);
      }
    }
    curve_color=15;

  }

}
/****/
void debug(float * m1,float * m2,int shift)
{
int nc;   // ­®¬¥à ª ­ « 
int i,j;
int x1,x2,y1,y2,y11,zerou,x;
float coef,zerr,max,min,r;
int h=shift-h_ch*2;
y11=h_ch*4+5+h_ch;

for (nc=0;nc<NCAN;nc++)
  {
  x1=0;    x2=x1+NPOINTS;
  y1=y11;
  y2=y1+h;
  zerou=(y2+y1)/2;

  setfillstyle(SOLID_FILL,0);
  bar(x1,y1,x2,y2);
  setlinestyle(SOLID_LINE,0,1);
  setcolor(EGA_BLUE);
  rectangle(x1,y1,x2,y2);

  max=min=m1[nc];
  for (i=1;i<NPOINTS;i++)
    {
    r=m1[i*NCAN+nc];
    if (min>r) min=r;
    if (max<r) max=r;
    }
  for (i=0;i<NPOINTS;i++)
    {
    r=m2[i*NCAN+nc];
    if (min>r) min=r;
    if (max<r) max=r;
    }
  zerr=(max+min)/2;
  coef=(float)(h-5)/(max-min);

  setlinestyle(SOLID_LINE,0,3);
  setcolor(EGA_LIGHTRED);
  for (i=1,x=x1;i<NPOINTS;i++,x++)
    {
    x2=x+1;
    r=m1[(i-1)*NCAN+nc];
    y1=zerou-(r-zerr)*coef;
    r=m1[i*NCAN+nc];
    y2=zerou-(r-zerr)*coef;
    line(x,y1,x2,y2);
    }

  setlinestyle(SOLID_LINE,0,1);
  setcolor(EGA_LIGHTGREEN);
  for (i=1,x=x1;i<NPOINTS;i++,x++)
    {
    x2=x+1;
    r=m2[(i-1)*NCAN+nc];
    y1=zerou-(r-zerr)*coef;
    r=m2[i*NCAN+nc];
    y2=zerou-(r-zerr)*coef;
    line(x,y1,x2,y2);
    }
  y11=y11+shift;
  }
}


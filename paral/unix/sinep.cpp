#include <cstdio>
#include "get_dat.h"
#include "sinepvar.h"
#include "sinepclc.h"
#include "snp_fltr.h"
#include <math.h>
#include <mpi.h>
#include "param.h"
/***/

//extern int fcloseall(void);
 
 PARAM_I Ncan("NCAN");
 PARAM_I Npoints("NPOINTS");
 PARAM_I Label("LABEL");
 PARAM_R Nms("NMS");
 PARAM_I Level("LEVEL");
 PARAM_I First_EP("FIRST_EP");
 PARAM_I Last_EP("LAST_EP");
 PARAM_R Campl("CA");
 PARAM_R Min_FRQ("MIN_FRQ");
 PARAM_R Max_FRQ("MAX_FRQ");
 PARAM_R Ord_FRQ("ORD_FRQ");
 PARAM_B Demodulation("DEMODULATION");
 PARAM_B Fourie("FOURIE");
 PARAM_B Pause("PAUSE");
 PARAM_B Spectr("SPECTR");
 PARAM_I Order("ORDER");
 PARAM_S Format("FORMAT");
 PARAM_R Coef_MUL("COEF_MUL");
 
double* mas,* buf,* buf1;
double * mas_sin, * mas_cos, * sum, *lsum ;
double * mas_sini,* mas_cosi,* mas1, *mas2;
double minr,maxr;
int iam,npr,nsum;

void revealing(int iam, int npr,double* mas,int h);

void main(int argc,char * argv[])
{
//----------
double time;
//---------
char name_dat[65];
char name_cfg[65];
long ind;
int i,j;
int pass;
char * er;

MPI_Init(&argc,&argv);

MPI_Comm_size(MPI_COMM_WORLD,&npr);
MPI_Comm_rank(MPI_COMM_WORLD,&iam);
time=MPI_Wtime();
if (argc>1)
  {
  i=read_PARAM(argv[1],&er);
  if (i!=0)
    { printf("%s",er);getchar();exit(0); }
  do name_cfg[i]=argv[0][i];
  while (name_cfg[i++] !='\x0');

  change_ext(name_cfg,".cfg");
  i=read_PARAM(name_cfg,&er);
  if (i!=0)
    { printf("%s",er);getchar();exit(0); }
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

const char ER1[]="юе заАан параметр %s";
    if (Ncan.ind=='\x0')
       { printf(ER1,Ncan.name); getchar();exit(0); }
    if (Nms.ind=='\x0')
       { printf(ER1,Nms.name); getchar();exit(0); }
    if (Npoints.ind=='\x0')
       { printf(ER1,Npoints.name); getchar();exit(0); }


    if (Ord_FRQ.ind) PHASE= new double [NCAN];

    mas     = new double [NPOINTS*NCAN+NCAN];
    buf     = new double [NPOINTS+1];
    buf1    = new double [NPOINTS+1];
    mas_sin = new double [NPOINTS*NCAN+NCAN];
    mas_cos = new double [NPOINTS*NCAN+NCAN];
    sum     = new double [NPOINTS*NCAN+NCAN];
    lsum    = new double [NPOINTS*NCAN+NCAN];
    mas_sini= new double [NPOINTS*NCAN+NCAN];
    mas_cosi= new double [NPOINTS*NCAN+NCAN];
    mas1    = new double [NPOINTS*NCAN+NCAN];
    mas2    = new double [NPOINTS*NCAN+NCAN];

//    i=0;
//    do out_inf[i]=argv[1][i];
//    while
//       (out_inf[i++]!='\x0');

    strcpy(out_inf,argv[1]);
    change_ext(out_inf,".fin");
//    printf("щмя выхоАного информационного файла [%s] ? ",out_inf);
//    sscanf("%s",name_dat);
//    i=0;
//    if (name_dat[0]>' ')
//      do out_inf[i]=name_dat[i];
//      while
//	 (out_inf[i++]!='\x0');

    create_inf(Ord_FRQ.ind);

if (DEMODULATION==0) goto only_fourie;

    ini_Gram(ORDER);
    if (LEVEL_SHOW==0)
    if(iam==0) printf("%s\n","ВЫДЕЛЕНИЕ ЕДИНИЧНЫХ ВП МЕТОДОМ КОМПЛЕКСНОЙ ДЕМОДУЛЯЦИИ");

    while (next_dat(mas,name_dat))
      {
      min_max(mas,&minr,&maxr);
      if ((iam==0) && (LEVEL_SHOW==0)) printf("%s\n",name_dat);
      revealing(iam,npr,mas,Ord_FRQ.ind);
      if(iam==0) write_result(name_dat,sum);
      if(iam==0) change_ext(name_dat,".sep");
      if(iam==0) add_to_inf(name_dat,Ord_FRQ.ind);
        
      } // while (next_dat(mas),name_dat)

    del_Gram();
        
only_fourie :
    delete mas2;
    delete mas1;
    delete mas_cosi;
    delete mas_sini;
    delete sum;
    delete lsum;
    delete mas_cos;
    delete mas_sin;
    delete buf1;
    delete buf;
    delete mas;
    if (Ord_FRQ.ind) delete PHASE;


    if ((FOURIE) || Spectr.val()) fourie_filtr(argv[1],FOURIE,Pause.val());

/*    if (LEVEL_SHOW) closegraph();*/

    }
  else  //  if (ini_dat(argv[1]))
    {
     printf("\nщнформационный файл %s не найАен\n",argv[1]);
     getchar();
    }

  }
else // if (argv>1)
  {
  printf("\nВызов программы:\nsinep.exe name.inf\n");
  printf(",гАе name.inf - имя информационного файла\n");
  getchar();
  }
 time=MPI_Wtime()-time;
 if(iam==0)
{ 
printf("Time calculation: %10.5f\n",time);
}
//fcloseall();

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

MPI_Finalize();        
 
}




void revealing(int iam,int npr, double* mas,int have_ord)
{ 
int max_freq,min_freq;
double pace_freq,freq;
int i,j,k,nsum;
long ind;
char ss[65];

for (i=0,ind=0;i<NPOINTS;i++)
  for (j=0;j<NCAN;j++,ind++) {lsum[ind]=0; sum[ind]=0.0;}
 nsum=ind;
pace_freq=1000.0/(NPOINTS*NMS);
min_freq=(int)(MIN_FRQ/pace_freq)+1;
max_freq=(int)(MAX_FRQ/pace_freq+0.9);
if (max_freq>NPOINTS/2) max_freq=NPOINTS/2;

for (k=min_freq+iam;k<max_freq;k+=npr)
  {
  freq=k*pace_freq; 
  modulation(freq*M_PI*2*NMS/1000.0,mas,mas_cos,mas_sin);
  printf("%i %i %6.2f \n",iam, k,freq);

  /* фильтрация (туАа-обратно) */
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
  extern void debug(double * m1,double * m2,int shift);
  
  if ( (have_ord) && (ORD_FRQ<=freq+pace_freq/2) && (ORD_FRQ>freq-pace_freq/2))
    phase_search(mas_sin,mas_cos,LABEL);

  subtracking(mas_cosi,mas_cos,mas_cosi);
  subtracking(mas_sini,mas_sin,mas_sini);

//  demodulation(freq*M_PI*2*NMS/1000.0,mas1,mas_cos,mas_sin);
  demodulation(freq*M_PI*2*NMS/1000.0,mas2,mas_cosi,mas_sini);


//  subtracking(mas2,mas1,mas2);
  summing(lsum,lsum,mas2);
  summing(lsum,lsum,mas2);
/* freq+=pace_freq*npr; */
 }
MPI_Allreduce(lsum,sum,nsum,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

}
/****/

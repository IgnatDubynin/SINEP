#include <string.h>
#include "sinepvar.h"

char WORD[17];
char * FORMAT;

double * PHASE;
char out_inf[65];

/* параметры графики */
int  maxx,maxy;
int  shift_curve,h_graph;
int h_ch;

int NPOINTS,NCAN,LABEL;
float NMS;

/* параметры программы */
int LEVEL_SHOW;
int FIRST_EP;     // выАеляемый Вх
int LAST_EP;
float CA;         //коэффициент усиления сигнала
float MIN_FRQ;    // частотный Аиапазон выАеляемого сигнала в чц
float MAX_FRQ;
float ORD_FRQ;
int ORDER;
int FOURIE;         // фильтрация цурье
int DEMODULATION;   // выАеление Вх АемоАуляцией
float COEF_MUL=1.0; // вспомогательный коэф. в .inf Аля EEGIN


void change_ext(char * name, const char * ext)
{
char *ptr_old,*ptr,*str;
ptr_old=name;

while ((ptr=strchr(ptr_old,'.'))>0) ptr_old=ptr+1;
if (ptr_old>name) *(ptr_old-1)='\0';
strcat(name,ext);
}

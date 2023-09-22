#include <dir.h>

char WORD[17];
char * FORMAT;

float * PHASE;
char out_inf[65];

/* параметры графики */
int  maxx,maxy;
int  shift_curve,h_graph;
int h_ch;

int NPOINTS,NCAN,LABEL;
float NMS;

/* параметры программы */
int LEVEL_SHOW;
int FIRST_EP;     // выделяемый ВП
int LAST_EP;
float CA;         //коэффициент усиления сигнала
float MIN_FRQ;    // частотный диапазон выделяемого сигнала в Гц
float MAX_FRQ;
float ORD_FRQ;
int ORDER;
int FOURIE;         // фильтрация Фурье
int DEMODULATION;   // выделение ВП демодуляцией
float COEF_MUL=1.0; // вспомогательный коэф. в .inf для EEGIN


void change_ext(char * name, const char * ex)
{
char s[MAXPATH];
char drive[MAXDRIVE];
char dir[MAXDIR];
char file[MAXFILE];
char ext[MAXEXT];

fnsplit(name,drive,dir,file,ext);
fnmerge(name,drive,dir,file,ex);

}

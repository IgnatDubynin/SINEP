#include <string.h>
#include "sinepvar.h"

char WORD[17];
char * FORMAT;

double * PHASE;
char out_inf[65];

/* ��������� ������� */
int  maxx,maxy;
int  shift_curve,h_graph;
int h_ch;

int NPOINTS,NCAN,LABEL;
float NMS;

/* ��������� ��������� */
int LEVEL_SHOW;
int FIRST_EP;     // ���������� ��
int LAST_EP;
float CA;         //����������� �������� �������
float MIN_FRQ;    // ��������� �������� ����������� ������� � ��
float MAX_FRQ;
float ORD_FRQ;
int ORDER;
int FOURIE;         // ���������� �����
int DEMODULATION;   // ��������� �� ������������
float COEF_MUL=1.0; // ��������������� ����. � .inf ��� EEGIN


void change_ext(char * name, const char * ext)
{
char *ptr_old,*ptr,*str;
ptr_old=name;

while ((ptr=strchr(ptr_old,'.'))>0) ptr_old=ptr+1;
if (ptr_old>name) *(ptr_old-1)='\0';
strcat(name,ext);
}

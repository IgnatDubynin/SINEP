#include <dir.h>

char WORD[17];
char * FORMAT;

float * PHASE;
char out_inf[65];

/* ��ࠬ���� ��䨪� */
int  maxx,maxy;
int  shift_curve,h_graph;
int h_ch;

int NPOINTS,NCAN,LABEL;
float NMS;

/* ��ࠬ���� �ணࠬ�� */
int LEVEL_SHOW;
int FIRST_EP;     // �뤥�塞� ��
int LAST_EP;
float CA;         //�����樥�� �ᨫ���� ᨣ����
float MIN_FRQ;    // ����� �������� �뤥�塞��� ᨣ���� � ��
float MAX_FRQ;
float ORD_FRQ;
int ORDER;
int FOURIE;         // 䨫����� ����
int DEMODULATION;   // �뤥����� �� �������樥�
float COEF_MUL=1.0; // �ᯮ����⥫�� ����. � .inf ��� EEGIN


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

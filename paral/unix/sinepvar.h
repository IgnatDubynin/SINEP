
extern int NPOINTS,NCAN,LABEL;
extern float NMS;
extern char WORD[17];
extern char * FORMAT;
extern double * PHASE; // ������ ��� ������ �������
extern char out_inf[65];

/* ��������� ������� */
extern int  maxx,maxy;
extern int  shift_curve,h_graph;
extern int h_ch;


/* ��������� ��������� */
extern int LEVEL_SHOW; // 0,1,2,3
extern float CA;  //����������� �������� �������
extern int FIRST_EP,LAST_EP;
extern float MIN_FRQ,MAX_FRQ; // ��������� �������� ����������� �������
extern float ORD_FRQ;
extern int ORDER,FOURIE,DEMODULATION;
extern float COEF_MUL;

extern void change_ext(char * name, const char * ex);

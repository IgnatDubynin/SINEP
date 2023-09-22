extern void out_curve(
		 float* buf,
		 int x0,int y0,int n,int h,
		 float minr,float maxr);
extern int curve_color,
	   curve_style,
	   curve_thick,
	   curve_paceox;
extern void drow_frame(char * part, char * name);
extern void drow_subhead(int n,char * text);
extern void drow_labels(int shifty);
extern void drow_ampl(int shifty,float minr,float maxr,int n);
extern void ini_fourie_filtr(void);
extern void drow_name_file(char * name);
extern void out_spectr(float * buf,int canal,int cnt_aver);
extern void error(char * str);
extern void subscript(int,int);


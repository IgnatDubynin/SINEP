extern void min_max(double* mas, double* minr, double* maxr);
extern void modulation(double betta,double* mas,double* mas_cos,double* mas_sin);
extern void demodulation(double betta,double* mas,double* mas_cos,double* mas_sin);
extern void filtr(double* mas);
extern void summing(double* sum,double* mas1,double* mas2);
extern void subtracking(double* sum,double* mas1,double* mas2);
extern void Gram(double* inp,double* out);
extern void ini_Gram(int m); // m  - максимальный поряАок многочлена
extern void del_Gram(void);
extern void phase_search(double * Sn, double * Cs, int stim);



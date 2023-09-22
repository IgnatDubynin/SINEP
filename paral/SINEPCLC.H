extern void min_max(float* mas, float* minr, float* maxr);
extern void modulation(float betta,float* mas,float* mas_cos,float* mas_sin);
extern void demodulation(float betta,float* mas,float* mas_cos,float* mas_sin);
extern void filtr(float* mas);
extern void summing(float* sum,float* mas1,float* mas2);
extern void subtracking(float* sum,float* mas1,float* mas2);
extern void Gram(float* inp,float* out);
extern void ini_Gram(int m); // m  - максимальный порядок многочлена
extern void del_Gram(void);
extern void phase_search(float * Sn, float * Cs, int stim);



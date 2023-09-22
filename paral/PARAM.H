union abscent_value{
      int i[4];
      int INT;
      float FLT;
       };

extern abscent_value GAP;
const Len_name=17;

class PARAM  {
      public :
      int dim;
      char ind,mod;        // ind - признак установки параметра
			     // mod - модификатор ( "<" или ">" ) при dim=0
      char name[Len_name];
      union {
	int i;
	int * ptri;
	float r;
	float * ptrr;
	char * s;
	char ** ss;
	     } v;
      PARAM * next;
      PARAM(char*);
//      ~PARAM();
virtual
  int PAR_READ(char*);
virtual
  void clear(void);
	      };

class PARAM_B : public PARAM{
      public:
      PARAM_B(char * str);
      inline int val(void) { return v.i; };
      void set(int);
      void clear(void) {ind=mod=dim=0;};
      int PAR_READ(char *);
      ~PARAM_B();
	       };

class PARAM_I : public PARAM{
      public:
      PARAM_I(char * str);
      PARAM_I(char * str,int i);
      inline int val(void) { return v.i; };
      inline int val(int n){ return v.ptri[n]; };
      void set(...);
      void clear(void);
      ~PARAM_I();
	       };

class PARAM_R : public PARAM{
      public:
      PARAM_R(char * str);
      PARAM_R(char * str,int i);
      inline float val(void) { return v.r; };
      inline float val(int n){ return v.ptrr[n]; };
      void set(...);
      void clear(void);
      int  PAR_READ(char*);
      ~PARAM_R();
	       };

class PARAM_S : public PARAM{
      public:
      PARAM_S(char * str);
      PARAM_S(char * str,int n);
      inline char* val(void)  { return v.s; };
      inline char* val(int n) { return v.ss[n-1]; };
      void set(char*);
      void set(int n,char*);
      void clear(void);
      int PAR_READ(char *);
      ~PARAM_S();
	       };
extern int read_PARAM(char * name,char** error);

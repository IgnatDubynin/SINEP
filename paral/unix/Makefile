MACHINE = LINUX

PROG = sinep
TESTOBJ = sinep.o get_dat.o lagrange.o sinepclc.o \
	  sinepvar.o snp_fltr.o param.o

CXX_SUN         = mpiCC
CXX_ALPHA       = cxx
CXX_NCUBE       = ncc
CXX_LINUX       = mpiCC

CCFLAGS_SUN      = -fast -O4
CCFLAGS_ALPHA    = -fast -O4 -I/usr/local/include
CCFLAGS_NCUBE    = -O 
CCFLAGS_LINUX    = -O

LDFLAGS_SUN       = -fast -O4
LDFLAGS_ALPHA     = -fast -O4 
LDFLAGS_NCUBE     = -O -Ncomm 3000000 
LDFLAGS_LINUX     = -O

LIB_SUN         =   
LIB_ALPHA       = -lmpi 
LIB_NCUBE       = -laztec -lmpi -lblasn -lf -lm 
LIB_LINUX       = 

CXX            = $(CXX_${MACHINE})
CCFLAGS        = $(CCFLAGS_${MACHINE})
LIBS           = $(LIB_${MACHINE}) 
LDFLAGS        = $(LDFLAGS_${MACHINE})
INCLUDE   = -I/usr/local/include

.SUFFIXES:

.SUFFIXES: .cpp .o


all: exe


exe: $(TESTOBJ)
#	${CXX} ${LDFLAGS} -o $(PROG)_${MACHINE} ${TESTOBJ} ${LIBS}
	${CXX} ${LDFLAGS} -o $(PROG) ${TESTOBJ} ${LIBS}

clean :
	rm -f *.o

.f.o : ; $(F77) -c ${FFLAGS} $*.f

.cpp.o : ; $(CXX) -c $(CCFLAGS) $(CDEFS) $*.cpp


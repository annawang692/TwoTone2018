#
# Unix/Linux GCC Makefile for Levenberg - Marquardt minimization
# Under windows, use Makefile.vc for MSVC
#

CC=gcc
CONFIGFLAGS=-ULINSOLVERS_RETAIN_MEMORY# YOU MIGHT WANT TO COMMENT THIS OUT FOR MULTITHREADED APPLICATIONS - see http://www.ics.forth.gr/~lourakis/levmar/changelog.html#LSMEM
#ARCHFLAGS=-march=pentium4 # YOU MIGHT WANT TO UNCOMMENT THIS FOR P4
FLAGS_64BIT=-fPIC# YOU SHOULD UNCOMMENT THIS FOR 64-bit
CFLAGS=$(FLAGS_64BIT) $(CONFIGFLAGS) $(ARCHFLAGS) -O3 -funroll-loops -Wall #-pg
#LAPACKLIBS_PATH=/usr/lib # WHEN USING LAPACK, CHANGE THIS TO WHERE YOUR COMPILED LIBS ARE!
LAPACKLIBS_PATH=/home/HoldenS/Dropbox/Private/development/twotone3alphaFinal/lapack-3.2.2# WHEN USING LAPACK, CHANGE THIS TO WHERE YOUR COMPILED LIBS ARE!
LDFLAGS=-L$(LAPACKLIBS_PATH) -L.
LIBOBJS=lm.o Axb.o misc.o lmlec.o lmbc.o lmblec.o
LIBSRCS=lm.c Axb.c misc.c lmlec.c lmbc.c lmblec.c
DEMOBJS=lmdemo.o
DEMOSRCS=lmdemo.c
AR=ar
RANLIB=ranlib
#LAPACKLIBS=-llapack -lblas -lgfortran	 # comment this line if you are not using LAPACK.
LAPACKLIBS=$(LAPACKLIBS_PATH)/lapack_LINUX.a $(LAPACKLIBS_PATH)/blas_LINUX.a  -lgfortran	
                                 # On systems with a FORTRAN (not f2c'ed) version of LAPACK, -lf2c is
                                 # not necessary; on others, -lf2c is equivalent to -lF77 -lI77

#LAPACKLIBS=-L/usr/local/atlas/lib -llapack -lcblas -lf77blas -latlas -lf2c # This works with the ATLAS updated lapack and Linux_P4SSE2
                                                                            # from http://www.netlib.org/atlas/archives/linux/

#LAPACKLIBS=-llapack -lgoto -lpthread -lf2c # This works with GotoBLAS
                                            # from http://www.tacc.utexas.edu/resources/software/

#LAPACKLIBS=-L/opt/intel/mkl/8.0.1/lib/32/ -lmkl_lapack -lmkl_ia32 -lguide -lf2c # This works with MKL 8.0.1 from
                                            # http://www.intel.com/cd/software/products/asmo-na/eng/perflib/mkl/index.htm

LIBS=$(LAPACKLIBS)

all: liblevmar.a lmdemo

liblevmar.a: $(LIBOBJS)
	$(AR) crv liblevmar.a $(LIBOBJS)
	$(RANLIB) liblevmar.a

lmdemo: $(DEMOBJS) liblevmar.a
	$(CC) $(LDFLAGS) $(DEMOBJS) -o lmdemo -llevmar $(LIBS) -lm

lm.o: lm.c lm_core.c lm.h misc.h compiler.h
Axb.o: Axb.c Axb_core.c lm.h misc.h
misc.o: misc.c misc_core.c lm.h misc.h
lmlec.o: lmlec.c lmlec_core.c lm.h misc.h
lmbc.o: lmbc.c lmbc_core.c lm.h misc.h  compiler.h
lmblec.o: lmblec.c lmblec_core.c lm.h misc.h

lmdemo.o: lm.h

clean:
	@rm -f $(LIBOBJS) $(DEMOBJS)

cleanall: clean
	@rm -f lmdemo
	@rm -f liblevmar.a

depend:
	makedepend -f Makefile $(LIBSRCS) $(DEMOSRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.


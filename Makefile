.c.o :
	$(CC) -c $(CFLAGS) $*.c
CC = gcc
# set CC to what makes most sense for you
# CC = /opt/local/bin/gcc-mp-4.8
#
CFLAGS = -O2
#
#
# If you have OpenMP support in the compiler, uncomment the following line
#
FOPENMP = -fopenmp
#
#
CFLAGS += $(FOPENMP)
all : pfield 
PFIELDO = pfield.o integrate_path.o phix.o calc_derivs.o loadmodels.o calcfpsi.o \
	bsstep.o odeint.o rk4.o rkqc.o nrutil.o \
	mmid.o rzextr.o
pfield : $(PFIELDO)
	$(CC) $(FOPENMP) -o pfield $(PFIELDO) -lm
clean :
	rm pfield *.o






.c.o :
	$(CC) -c $(CFLAGS) $*.c
CC=gcc
CFLAGS = -O2
#CFLAGS = -g
#
all : pfield redopfield 
PFIELDO = pfield.o integrate_path.o phix.o calc_derivs.o loadmodels.o calcfpsi.o \
	qromb.o bsstep.o odeint.o rk4.o rkqc.o polint.o trapzd.o nrutil.o \
	mmid.o rzextr.o
pfield : $(PFIELDO)
	$(CC) -o pfield $(PFIELDO) -lm
REDOPFIELDO = redopfield.o loadmodels.o nrutil.o 
redopfield : $(REDOPFIELDO)
	$(CC) -o redopfield $(REDOPFIELDO) -lm
clean :
	rm pfield redopfield *.o






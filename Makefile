.c.o :
	$(CC) -c $(CFLAGS) $*.c
CC=gcc  
#
all : pfield 
PFIELDO = pfield.o integrate_path.o phix.o calc_derivs.o loadmodels.o calcfpsi.o \
	qromb.o bsstep.o odeint.o rk4.o rkqc.o polint.o trapzd.o nrutil.o \
	mmid.o rzextr.o
pfield : $(PFIELDO)
	$(CC) -o pfield $(PFIELDO) -lm
clean :
	rm pfield *.o






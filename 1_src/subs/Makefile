FFLAGS = -g
FC = gfortran
SUBS = dfftr.o fft.o locast.o rdlyrs.o \
	npowr2.o yesno.o 		

subs.a: $(SUBS)
	ar cr subs.a $(SUBS)
	ranlib subs.a
clean:
	/bin/rm *.o


dfftr.o:  dfftr.f 
fft.o:  fft.f 
locast.o:  locast.f 
rdlyrs.o:  rdlyrs.f 
yesno.o	: yesno.f
npowr2.o : npowr2.f
FFLAGS = -O3
FC = gfortran
SACLIBRARY = $(SACLIB)
SACLIBRARY = /usr/local/sac/lib/sacio.a
SUBS = dfftr.o fft.o locast.o rdlyrs.o \
	npowr2.o yesno.o 
OBJ = respknt.o kntfun.o rcvrfn.o ifmat.o abm.o

respknt: $(OBJ) $(SUBS) kennet.inc
	$(FC) $(FFLAGS) -o respknt  $(OBJ) $(SUBS) $(SACLIBRARY)

respknt.o: respknt.f kennet.inc
	$(FC) -c respknt.f kennet.inc
rcvrfn.o: rcvrfn.f kennet.inc
	$(FC) -c rcvrfn.f kennet.inc
ifmat.o:  ifmat.f kennet.inc
	$(FC) -c ifmat.f kennet.inc
kntfun.o: kntfun.f 
	$(FC) -c kntfun.f 
abm.o: abm.f
	$(FC) -c abm.f

dfftr.o: dfftr.f 
fft.o: fft.f 
locast.o: locast.f 
rdlyrs.o: rdlyrs.f 
yesno.o	: yesno.f
npowr2.o : npowr2.f

clean:
	/bin/rm *.o

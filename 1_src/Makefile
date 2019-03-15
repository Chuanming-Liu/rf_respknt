FFLAGS = -O3
FC = gfortran
SACLIBRARY = $(SACLIB)
TJOLIB = ./Subs/subs.a
SACLIBRARY = /usr/local/sac/lib/sacio.a

OBJ = respknt.o kntfun.o rcvrfn.o ifmat.o abm.o
respknt: $(OBJ) kennet.inc
	$(FC) $(FFLAGS) -o respknt  $(OBJ) $(TJOLIB) $(SACLIBRARY)
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

clean:
	/bin/rm *.o

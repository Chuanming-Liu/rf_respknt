# f2py -m module_name -h sig_file.pyf list_of_fortran_files
rm respknt.pyf
f2py -m respknt -h respknt.pyf  respknt_interface.f kntfun.f rcvrfn.f ifmat.f abm.f dfftr.f fft.f npowr2.f kennet.inc

# f2py -c --fcompiler=gnu95 sig_file.pyf list_of_fortran_files
f2py -c --f77flags="-ffixed-line-length-none -O3" --f90flags="-O3" respknt.pyf respknt_interface.f kntfun.f rcvrfn.f ifmat.f abm.f dfftr.f fft.f npowr2.f kennet.inc

cc = icc
cflags=/C /DICC /DOS2 /Sd /Gd /Tl5 /Gf /O /Oi64 /Gi /G4 /Fi /Si

fit2.exe:fit2.obj mrqfit2.obj funclib2.obj fitutil2.obj fitcmds2.obj solve_da.obj plot.obj linear2.obj popen.obj
	$(cc) /B/stack:128000 fit2.obj mrqfit2.obj funclib2.obj fitutil2.obj fitcmds2.obj solve_da.obj plot.obj linear2.obj popen.obj

mrqfit2.obj:mrqfit2.c fit.h makefile
	$(cc) $(cflags) mrqfit2.c

fitcmds2.obj:fitcmds2.c fit.h makefile
	$(cc) $(cflags) fitcmds2.c

funclib2.obj:funclib2.c fit.h makefile
	$(cc) $(cflags) funclib2.c

fitutil2.obj:fitutil2.c fit.h makefile
	$(cc) $(cflags) fitutil2.c

solve_da.obj:solve_da.c fit.h makefile
	$(cc) $(cflags) solve_da.c

linear2.obj:linear2.c fit.h makefile
	$(cc) $(cflags) linear2.c

popen.obj:popen.c fit.h makefile
	$(cc) $(cflags) popen.c

plot.obj:plot.c makefile
	$(cc) $(cflags) plot.c

fit2.obj:fit2.c fit.h makefile
	$(cc) $(cflags) fit2.c

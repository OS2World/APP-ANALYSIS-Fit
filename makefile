cc = gcc
cflags=-c -O2 -DOS2 -m486

fit2.exe:fit2.o mrqfit2.o funclib2.o fitutil2.o fitcmds2.o solve_da.o plot.o linear2.o
	$(cc) fit2.o mrqfit2.o funclib2.o fitutil2.o fitcmds2.o solve_da.o plot.o linear2.o -lm -o fit2.exe

.c.o:
	$(cc) $(cflags) $<

mrqfit2.o:mrqfit2.c fit.h
	$(cc) $(cflags) mrqfit2.c

fitcmds2.o:fitcmds2.c fit.h
	$(cc) $(cflags) fitcmds2.c

fitutil2.o:fitutil2.c fit.h
	$(cc) $(cflags) fitutil2.c

solve_da.o:solve_da.c fit.h
	$(cc) $(cflags) solve_da.c

linear2.o:linear2.c fit.h
	$(cc) $(cflags) linear2.c

plot.o:plot.c fit.h
	$(cc) $(cflags) plot.c

fit2.o:fit2.c fit.h
	$(cc) $(cflags) fit2.c


fitplot.exe: fitplot.ex fitplot.res
	rc fitplot.res fitplot.ex
	copy fitplot.ex fitplot.exe

fitplot.res:  fitplot.rc  fitplot.h  fitplot.h
	rc -r fitplot.rc
  
fitplot.ex:  fitplot.c fitplot.h
	icc /B/A:16 /B/NOI /B/PM:PM /Fefitplot.ex /O /G4 /Gf /Gi /Gm /Gd- /Tl30 /Fi /Si fitplot.c

fp.exe: fp.c
	gcc fp.c -lm -O2 -o fp.exe


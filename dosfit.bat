gcc -c -DDOS fit2.c
gcc -c -DDOS fitcmds2.c 
gcc -c -DDOS fitutil2.c 
gcc -c -DDOS funclib2.c 
gcc -c -DDOS linear2.c 
gcc -c -DDOS mrqfit2.c
gcc -c -DDOS plot.c 
gcc -c -DDOS solve_da.c
gcc -lm fit2.o fitcmds2.o fitutil2.o funclib2.o linear2.o mrqfit2.o plot.o solve_da.o -lm -o a.out

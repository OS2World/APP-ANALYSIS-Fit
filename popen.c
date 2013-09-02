#include <stdio.h>
#include <process.h>
#define INCL_DOS
#include <os2.h>
#include <io.h>

/*courtesy of Roger Fearick */
FILE *popen(char *file, char mode[5]){

HFILE    hWritePipe ;
HFILE    hR2 ;
int      hR, save ;
FILE     *gnu ;
int error;
    
/* retain copy of old file descriptor */  
save = _dup( 0 ) ;         

/* create pipe */
DosCreatePipe( &hR2, &hWritePipe, 1024 ) ;

/* make '0' refer to pipe */
hR = 0 ;                   
_dup2( hR2, hR ) ;

/* start gnuplot */
error = _spawnlp( P_NOWAIT, "gnuplot.exe", "gnuplot.exe", NULL ) ;
if(error == -1) printf("error in _spawnl()\n");

/* open a file handle to write to gnuplot */
gnu = _fdopen( hWritePipe, "w" ) ;
_dup2( save, 0 ) ;        /* restore original stdin */

return gnu;
}

void pclose(FILE *stream){

}

void sleep(int num){
unsigned long milliseconds;

milliseconds = 1000 * num;
DosSleep(milliseconds);
}


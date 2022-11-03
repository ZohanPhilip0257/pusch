#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "wnNrExeTimeX86Intel.h"

/*
RDTSC is a current time-stamp counter variable ,which is a 64-bit variable, into registers (edx:eax).
TSC(time stamp counter) is incremented every cpu tick (1/CPU_HZ)
*/

unsigned long long int getCycles()
{
 unsigned int low;
 unsigned int high;
 __asm__ __volatile__("rdtsc" : "=a" (low), "=d" (high));

 return ((unsigned long long int)high << 32) | low;
}



void wnNrExeTimeX86Intel(int Iter)
{

    FILE *estPtr = fopen("estTime.txt","w+");
    FILE *corrPtr = fopen("corrTime.txt","w+");
    for (int iterIdx = 0; iterIdx < PROF_NUM_OF_ITER; iterIdx++)
    {
        fprintf (estPtr, "%lu\n",exeTimeEst[iterIdx]);
        fprintf (corrPtr, "%lu\n",exeTimeCorr[iterIdx]);
    }

    fclose(estPtr);
    fclose(corrPtr);
}

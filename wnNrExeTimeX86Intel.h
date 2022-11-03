
#ifndef WN_NR_X86_INTEL_PROFILE
#define WN_NR_X86_INTEL_PROFILE

#define PROF_NUM_OF_ITER (10000)
#include<time.h>
#include<stdint.h>
#define PRE_PROFILE_FLAG  (0)
#define POST_PROFILE_FLAG (1)

void wnNrExeTimeX86Intel(int);

#ifdef PROFILE_GNB
    static struct timespec estStTime, estEndTime;
    static struct timespec corrStTime, corrEndTime;
#endif

uint64_t exeTimeCorr[PROF_NUM_OF_ITER]; // Time 
uint64_t exeTimeEst[PROF_NUM_OF_ITER]; // Time 

// uint64_t exeTimeEst[PROF_NUM_OF_ITER]; // Time 

#endif
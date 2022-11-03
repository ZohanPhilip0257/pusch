

#ifndef WN_NR_PHY_FFT_CP_REMOVAL_H
#define WN_NR_PHY_FFT_CP_REMOVAL_H

#include "wnNrPhygNBUlChnlConfigInit.h"
#include <fftw3.h>

#define ALIGN_BOUNDARY		(32)

fftwf_plan planFft;
fftwf_complex *fftShift, *fftIn;

__attribute__ ((aligned(ALIGN_BOUNDARY))) cmplx IfftCpInsert[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_IFFT_CP_DATA_IN_SLOT];
__attribute__ ((aligned(ALIGN_BOUNDARY))) cmplx IfftCpIn[MAX_IFFT_CP_DATA_IN_SLOT];

#endif
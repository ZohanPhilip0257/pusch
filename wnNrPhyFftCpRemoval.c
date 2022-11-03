#include <stdio.h>
#include<string.h>
#include<stdlib.h>
#include "wnNrPhyFftCpRemoval.h"
#include "wnNrPhygNBReDeMapper.h"


wnVoid *wnNrPhyFftCpRemoval(commonUlConfigT* commonUlConfig)
{

    wnUInt16 idx, idx2 = 0, antennaPort, halfLengthPrb = (commonUlConfig->bwInPrbs*12)>>1, halfLengthPrbBytes = halfLengthPrb*sizeof(fftwf_complex), nSample, fullLengthBytes = commonUlConfig->fftSize*sizeof(fftwf_complex);
    wnUInt16 nAddCpSamples = (commonUlConfig->fftSize>>7)<<commonUlConfig->numerology;

    // Plan for FFT Operation is created
    planFft = fftwf_plan_dft_1d(commonUlConfig->fftSize, fftIn, fftShift, FFTW_FORWARD, FFTW_ESTIMATE);

    // CP Removal, FFT, FFT Shift
	for (antennaPort = 0; antennaPort < commonUlConfig->nPhysicalAnt; antennaPort++)
	{
        idx2 = 0;
	    for (idx = 0; idx < NUMBER_OF_OFDM_SYMBOLS_SLOT ; idx++)
	    {   
	        /********************** CP Addition *********************/
	        nSample = commonUlConfig->cpSyms + nAddCpSamples * ((NUMBER_OF_OFDM_SYMBOLS_SLOT*commonUlConfig->nSlotNumber + idx) % (7<<commonUlConfig->numerology) == 0);
	        memcpy (&fftIn[0], &IfftCpInsert[antennaPort][idx2 + nSample], fullLengthBytes);
	        idx2 += commonUlConfig->fftSize + nSample; 

            if(idx == 0)
            {
                if(antennaPort == 0)
                {
                    printf("IfftCpInsert[antennaPort][idx2 + nSample].real: %f\n", IfftCpInsert[antennaPort][idx2 + nSample].real);
                    printf("IfftCpInsert[antennaPort][idx2 + nSample].imag: %f\n", IfftCpInsert[antennaPort][idx2 + nSample].imag);
                }
            }

	        /********************** FFT ****************************/
	        fftwf_execute(planFft);  

	        /********************** FFT Shift ************************/
	        memcpy(&reDemapperIn[antennaPort][idx][halfLengthPrb], fftShift, halfLengthPrbBytes); // commonUlConfig->bwInPrbs/2 implies half size and halfLengthPrbBytes = (Length of 2nd half = commonUlConfig->fftSize/2)*(sizeof(fftwf_complex) = 8)
	        memcpy(&reDemapperIn[antennaPort][idx][0], &fftShift[commonUlConfig->fftSize - halfLengthPrb], halfLengthPrbBytes); // commonUlConfig->bwInPrbs/2 implies half size and halfLengthPrbBytes = (Length of 1st half = commonUlConfig->fftSize/2)*(sizeof(fftwf_complex) = 8)
            if(idx == 0)
            {
                if(antennaPort == 0)
                {
                    printf("reDemapperIn[antennaPort][idx][halfLengthPrb].real: %f\n", reDemapperIn[antennaPort][idx][halfLengthPrb].real);
                    printf("reDemapperIn[antennaPort][idx][halfLengthPrb].imag: %f\n", reDemapperIn[antennaPort][idx][halfLengthPrb].imag);
                }
            }
	    }  
	}

	// Termination of FFT Plan
	fftwf_destroy_plan(planFft);
    fftwf_free(fftIn);     
    fftwf_free(fftShift);


    #ifdef DEBUG
    FILE *fid, *fid1;

    fid = fopen("IOs/cpremovalReal.txt","w");
    fid1 = fopen("IOs/cpremovalImag.txt","w");
    
    for (antennaPort = 0; antennaPort < commonUlConfig->nPhysicalAnt; antennaPort++)
    {
        for (idx = 0; idx < NUMBER_OF_OFDM_SYMBOLS_SLOT; idx++)
        {
            for (idx2 = 0; idx2 < commonUlConfig->bwInPrbs*12; idx2++)
            {
                fprintf(fid,"%f\n", reDemapperIn[antennaPort][idx][idx2].real);
                fprintf(fid1,"%f\n", reDemapperIn[antennaPort][idx][idx2].imag);
            }
        }
    }
    fclose (fid);
    fclose (fid1);
    #endif

	return NULL;
}
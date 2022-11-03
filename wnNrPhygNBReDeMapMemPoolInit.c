#include <stdio.h>
#include<string.h>
#include "wnNrPhygNBReDeMapper.h"
#include "wnNrPhyFftCpRemoval.h"


wnVoid* wnNrPhygNBReDeMapMemPoolInit()
{
	// Assign memory elements related to demapper and FFTCP module with 0. 
	memset (&reDemapperIn[0][0][0], 0, MAX_SUMIMO_PUSCH_ANTENNA_PORTS
										*NUMBER_OF_OFDM_SYMBOLS_SLOT
										*MAX_SUB_CARRIERS_OFDM_SYMBOL
										*sizeof (cmplx));

	memset (&IfftCpInsert[0][0], 0, MAX_SUMIMO_PUSCH_ANTENNA_PORTS
									*MAX_IFFT_CP_DATA_IN_SLOT
									*sizeof (cmplx));


	fftShift = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * MAX_ZERO_PADDED_SUB_CARRIERS_OFDM_SYMBOL);

    fftIn = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * MAX_ZERO_PADDED_SUB_CARRIERS_OFDM_SYMBOL);

	wnUInt32 iqSamples;

	// printf("gnbUlChnlSt.commonUlConfig.allocatedBw: %d\n", gnbUlChnlSt.commonUlConfig.allocatedBw);

	// if(gnbUlChnlSt.commonUlConfig.allocatedBw == 100)
	// {
	// 	iqSamples = 61440;
	// }
	// else if(gnbUlChnlSt.commonUlConfig.allocatedBw == 50)
	// {
	// 	iqSamples = 30720;
	// }

	// Read Channel out real and imaginary values.
    FILE *fid1 = fopen ("IOs/channelOutReal.txt","r");
    FILE *fid2 = fopen ("IOs/channelOutImag.txt","r");
	printf("fid1: %d\n", fid1);
	printf("fid2: %d\n", fid2);

    for (wnUInt8 antIdx = 0; antIdx < 1; antIdx++)
    {
	    for (wnUInt32 idx = 0; idx < 30720; idx++)
	    {
			// printf("[%d][%d] \n", antIdx, idx);
	    	fscanf (fid1,"%f\n",&IfftCpInsert[antIdx][idx].real);
	    	fscanf (fid2,"%f\n",&IfftCpInsert[antIdx][idx].imag);
	    }
	}
    fclose(fid1);
    fclose(fid2);

	printf("IfftCpInsert[0][0].real: %f\n", IfftCpInsert[0][0].real );
	printf("IfftCpInsert[0][0].imag: %f\n", IfftCpInsert[0][0].imag );


	return NULL;
}
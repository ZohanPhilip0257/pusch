#include "wnNrPhyPuschHeader.h"
#include "wnNrPhygNBReDeMapper.h"

/** 
* @brief 
* 1: Remove's location offset's in both Time and Frequency(In frequency domain, not within BWPart).
* 2: This function changes received data format from "Real Imag Real Imag ..." to
* "Real Real ... Imag Imag "
**/
wnVoid wnNrPhyPuschdemapping(puschConfigT* nrUlPhyPuschPDUs)
{
    wnInt32 loopIdx, portIdx, reIdx, newreIdx, dmrsSymIdx, reStartIdx = (nrUlPhyPuschPDUs->nBWPStart + nrUlPhyPuschPDUs->nRBStart)*12,\
    reEndIdx = (nrUlPhyPuschPDUs->nBWPStart + nrUlPhyPuschPDUs->nRBStart + nrUlPhyPuschPDUs->nRBSize)*12,\
    numberOfDataInOfdmSymbol = MAX_PUSCH_INFO_PER_SYMBOL, numberOfHalfDataInOfdmSymbol = numberOfDataInOfdmSymbol>>1,\
    numberOf3quartersDataInOfdmSymbol = numberOfHalfDataInOfdmSymbol + numberOfDataInOfdmSymbol;

    __attribute__ ((aligned(ALIGN_BOUNDARY))) wnFlt temp[8];
    __m256 inter;
    __m256i dataShift = _mm256_set_epi32 (7,5,3,1,6,4,2,0);
    __m256i dmrsShift = _mm256_set_epi32 (7,3,6,2,5,1,4,0);
    wnUInt8 symIdx;
    
    // Loop runs for Demapping all the Antenna Ports
    for ( portIdx = 0; portIdx < nrUlPhyPuschPDUs->nAntennaPorts; portIdx++)
    {
        // printf("Before rxFdSamples: %f\n", rxFdSamples[portIdx][0][0]);
        dmrsSymIdx = 0;
        symIdx = 0;
        // Loop runs to demap and adjust data in every symbol such that it is packed appropriately for Channel Estimator
        for ( loopIdx = nrUlPhyPuschPDUs->nStartSymbolIndex; loopIdx < (nrUlPhyPuschPDUs->nStartSymbolIndex + nrUlPhyPuschPDUs->nNrOfSymbols); loopIdx++)
        {
            //Condition check whether the symbol has DMRS or Data. Differentiation is made as DMRS and Data Packing fo Channel Estimation is different
            if ( loopIdx == dmrs_s[dmrsSymIdx])
            {
                for (reIdx = reStartIdx; reIdx < reEndIdx; reIdx+=4)
                {
                    newreIdx = (reIdx - reStartIdx)>>1;
                    inter = _mm256_loadu_ps ((wnFlt *)&reDemapperIn[portIdx][loopIdx][reIdx]);
                    inter = _mm256_permutevar8x32_ps (inter, dmrsShift);
                    _mm256_store_ps (&temp[0], inter);
                    memcpy(&rxFdSamples[portIdx][symIdx][newreIdx], &temp[0], 8);
                    memcpy(&rxFdSamples[portIdx][symIdx][numberOfHalfDataInOfdmSymbol + newreIdx], &temp[2], 8);
                    memcpy(&rxFdSamples[portIdx][symIdx][numberOfDataInOfdmSymbol + newreIdx], &temp[4], 8);
                    memcpy(&rxFdSamples[portIdx][symIdx][numberOf3quartersDataInOfdmSymbol + newreIdx], &temp[6], 8);
                }
                dmrsSymIdx++;
                symIdx++;
            }
            else
            {             
                for (reIdx = reStartIdx; reIdx < reEndIdx; reIdx+=4)
                {
                    newreIdx = reIdx - reStartIdx;
                    inter = _mm256_loadu_ps ((wnFlt *)&reDemapperIn[portIdx][loopIdx][reIdx]);
                    inter = _mm256_permutevar8x32_ps (inter, dataShift);
                    _mm256_store_ps (&temp[0], inter);
                    memcpy(&rxFdSamples[portIdx][symIdx][newreIdx], &temp[0], 16);
                    memcpy(&rxFdSamples[portIdx][symIdx][numberOfDataInOfdmSymbol + newreIdx], &temp[4], 16);
                }
                symIdx++;
            }
        }
        // printf("After rxFdSamples: %f\n", rxFdSamples[portIdx][0][0]);
    }   // End of loop that demaps all the ports


    #ifdef DEBUG
        FILE *fid;
        fid = fopen("IOs/puschDemapperOut.txt","w");

        for ( portIdx = 0; portIdx < nrUlPhyPuschPDUs->nAntennaPorts; portIdx++)
        {
            // printf("rxFdSamples: %f\n", rxFdSamples[portIdx][0][0]);
            for ( loopIdx = nrUlPhyPuschPDUs->nStartSymbolIndex; loopIdx < (nrUlPhyPuschPDUs->nStartSymbolIndex + nrUlPhyPuschPDUs->nNrOfSymbols); loopIdx++)
            {
                for (reIdx = 0; reIdx < numberOfDataInOfdmSymbol*2; reIdx++)//numberOfDataInOfdmSymbol 270*12*
                {
                    fprintf(fid,"%+f\n",rxFdSamples[portIdx][loopIdx][reIdx]);
                }
            }
        }
        fclose (fid);
    #endif

return;
}

/*
Doubts:
    1. BWP Start and RB Start Inclusion
*/

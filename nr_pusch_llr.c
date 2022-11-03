/********************************************************************************
 *
 * WISIG NETWORKS CONFIDENTIAL
 * Copyright (c) 2016-2020, WiSig Networks Pvt Ltd. All rights reserved.
 * www.wisig.com
 *
 * All information contained herein is property of WiSig Networks Pvt Ltd.
 * unless otherwise explicitly mentioned.
 *
 * The intellectual and technical concepts in this file are proprietary
 * to WiSig Networks and may be covered by granted or in process national
 * and international patents and are protect by trade secrets and
 * copyright law.
 *
 * Redistribution and use in source and binary forms of the content in
 * this file, with or without modification are not permitted unless
 * permission is explicitly granted by WiSig Networks.
 * If WiSig Networks permits this source code to be used as a part of open
 * source project, the terms and conditions of CC-By-ND (No Derivative) license
 * (https://creativecommons.org/licenses/by-nd/4.0/) shall apply.
 ********************************************************************************/

 /**
  * @brief This file includes all the function defintions
  *         for NR PUSCH Channel Estimation and Equalization.
  * @file    : nr_pusch.c
  * @ingroup : nr_pusch
  * @author  : MIRZA SAMI BAIG
  **/

// #include "common_def.h"
#include "nr_pusch_llr.h"
// #define SNR24

void nr_pusch_llr_compute_avx2(wnUInt8 segPart, wnUInt8* startSym, wnUInt8* endSym, commonUlConfigT* pNrUlCommonParams,
                                puschConfigT* pNrPuschInParams,
                                P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams )
{
    wnUInt16 nRBSize = pNrPuschInParams->nRBSize;
    wnUInt16 modOrder = pNrPuschInParams->modulationOrder;
    wnUInt8 nNrOfLayers = pNrPuschInParams->nNrOfLayers;
    __m256 SGN = _mm256_set1_ps(-0.0f);
    // wnFlt* MemTemp;
    __m256 absReal, absImag;
    __m256 MEM1, MEM2, MEM3, MEM4, MEM5, MEM6, MEM7, MEM8;
    __m256 MEM11, MEM21, MEM31, MEM41, MEM51, MEM61, MEM71, MEM81;
    // __m256 MEM22, MEM33;
    // __m256 MEM12;
    __m256 Mask1Real, Mask1Imag;
    // __m256 RES1;
    __m256 RES1Real, RES1Imag;
    __m256 RES2Real, RES2Imag;
    // __m256 RES3;
    // __m256d RES0d, RES1d, RES2d, RES3d;
    // __m256d MEM1d, MEM2d, MEM3d;
    // __m256d RES11d, RES22d, RES33d, RES44d;
    // __m256 RES11, RES22, RES33, RES44;

    __m256 *AVXReal  = (__m256*) pNrPuschOutParams->layerDemapperOutReal;
    __m256 *AVXImag  = (__m256*) pNrPuschOutParams->layerDemapperOutImag;

    wnUInt32 idxOffset = 0;
    wnInt32 inpLen;
    wnUInt8 symsComp;
    if(segPart == SEG1)
    {
        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
        {
            symsComp = 0;
            idxOffset = 0;
            inpLen = nRBSize*12*nNrOfLayers*(endSym[0]-startSym[0]);
        }
        else //if(pNrPuschInParams->nNrOfDMRSSymbols == 2)
        {
            symsComp = 0;
            idxOffset = 0;
            inpLen = nRBSize*12*nNrOfLayers*(endSym[0]-startSym[0]-1);
        }
        
    }
    else if(segPart == SEG2)
    {
        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
        {
            symsComp = endSym[0]-startSym[0];
            idxOffset = nRBSize*12*modOrder*nNrOfLayers*symsComp;
            inpLen = nRBSize*12*nNrOfLayers*(endSym[1]-startSym[1]);
        }
        else //if(pNrPuschInParams->nNrOfDMRSSymbols == 2)
        {
            symsComp = endSym[0]-startSym[0]-1;
            idxOffset = nRBSize*12*modOrder*nNrOfLayers*symsComp;
            inpLen = nRBSize*12*nNrOfLayers*(endSym[1]-startSym[1]-1);
        }
    }
    else if(segPart == SEG3)
    {
        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
        {
            symsComp = (endSym[0]-startSym[0])+(endSym[1]-startSym[1]);
            idxOffset = nRBSize*12*modOrder*nNrOfLayers*symsComp;
            inpLen = nRBSize*12*nNrOfLayers*(endSym[2]-startSym[2]);
        }
        else //if(pNrPuschInParams->nNrOfDMRSSymbols == 2)
        {
            symsComp = (endSym[0]-startSym[0]-1)+(endSym[1]-startSym[1]-1);
            idxOffset = nRBSize*12*modOrder*nNrOfLayers*symsComp;
            inpLen = nRBSize*12*nNrOfLayers*(endSym[2]-startSym[2]-1);
        }
        
    }
    else //if(segPart == SEG4)
    {
        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
        {
            symsComp = (endSym[0]-startSym[0])+(endSym[1]-startSym[1])+(endSym[2]-startSym[2]);
            idxOffset = nRBSize*12*modOrder*nNrOfLayers*symsComp;
            inpLen = nRBSize*12*nNrOfLayers*(endSym[3]-startSym[3]);
        }
        else //if(pNrPuschInParams->nNrOfDMRSSymbols == 2)
        {
            symsComp = (endSym[0]-startSym[0]-1)+(endSym[1]-startSym[1]-1)+(endSym[2]-startSym[2]-1);
            idxOffset = nRBSize*12*modOrder*nNrOfLayers*symsComp;
            inpLen = nRBSize*12*nNrOfLayers*(endSym[3]-startSym[3]-1);
        }
    }


    __m256 *FinalOut1 = (__m256*)(&llrOut[idxOffset]);//llrOut[idxOffset]
    // __m256 *FinalOut2 = FinalOut1+1;
    // __m256 *FinalOut3 = FinalOut1+2;
    // __m256 *FinalOut4 = FinalOut1+3;
    // __m256 *FinalOut5 = FinalOut1+4;
    // __m256 *FinalOut6 = FinalOut1+5;
    // __m256 *FinalOut7 = FinalOut1+6;
    // __m256 *FinalOut8 = FinalOut1+7;

    __m256 tempStore;
    wnFlt *tempStore_f = (wnFlt*)&tempStore;

    wnUInt32 i;
    wnUInt8 cnt = 0;
    wnUInt8 rem;
    if(modOrder == 2)
    {
        for(i = 0;i<(inpLen>>3);i++)
        {
            MEM1 = _mm256_unpacklo_ps( *AVXReal, *AVXImag);
            MEM2 = _mm256_unpackhi_ps( *AVXReal, *AVXImag);
            _mm256_storeu_ps(&llrOut[idxOffset], _mm256_permute2f128_ps(MEM1, MEM2, 0x20));
            _mm256_storeu_ps(&llrOut[idxOffset+8], _mm256_permute2f128_ps(MEM1, MEM2, 0x31));
            
            AVXReal += 1;
            AVXImag += 1;
            idxOffset += 16;
        }
        cnt = 0;
        for(i = ((inpLen>>3)<<3); i<inpLen;i++)
        {
            llrOut[idxOffset+cnt] = pNrPuschOutParams->layerDemapperOutReal[i];
            llrOut[idxOffset+cnt+1] = pNrPuschOutParams->layerDemapperOutImag[i];
            cnt += 2;
        }

    }
    else if(modOrder == 4)
    {

        __m256 Scalar1 = _mm256_set1_ps(0.6324555f);//2/sqrt(10)
        for(i = 0;i<(inpLen>>3);i++)
        {
            absReal = _mm256_andnot_ps(SGN, *AVXReal);
            absImag = _mm256_andnot_ps(SGN, *AVXImag);
            RES1Real = _mm256_add_ps( _mm256_or_ps(SGN, absReal), Scalar1);
            RES1Imag = _mm256_add_ps( _mm256_or_ps(SGN, absImag), Scalar1);

            MEM1 = _mm256_unpacklo_ps(*AVXReal, *AVXImag);
            MEM2 = _mm256_unpacklo_ps(RES1Real, RES1Imag);
            MEM3 = _mm256_unpackhi_ps(*AVXReal, *AVXImag);
            MEM4 = _mm256_unpackhi_ps(RES1Real, RES1Imag);

            MEM11 = _mm256_shuffle_ps(MEM1, MEM2, 0b01000100);
            MEM21 = _mm256_shuffle_ps(MEM1, MEM2, 0b11101110);
            MEM31 = _mm256_shuffle_ps(MEM3, MEM4, 0b01000100);
            MEM41 = _mm256_shuffle_ps(MEM3, MEM4, 0b11101110);

            // *FinalOut1 = _mm256_permute2f128_ps(MEM11, MEM21, 0x20);
            // *FinalOut2 = _mm256_permute2f128_ps(MEM31, MEM41, 0x20);
            // *FinalOut3 = _mm256_permute2f128_ps(MEM11, MEM21, 0x31);
            // *FinalOut4 = _mm256_permute2f128_ps(MEM31, MEM41, 0x31);

            _mm256_storeu_ps(&llrOut[idxOffset], _mm256_permute2f128_ps(MEM11, MEM21, 0x20));
            _mm256_storeu_ps(&llrOut[idxOffset+8], _mm256_permute2f128_ps(MEM31, MEM41, 0x20));
            _mm256_storeu_ps(&llrOut[idxOffset+16], _mm256_permute2f128_ps(MEM11, MEM21, 0x31));
            _mm256_storeu_ps(&llrOut[idxOffset+24], _mm256_permute2f128_ps(MEM31, MEM41, 0x31));

            AVXReal = AVXReal + 1;
            AVXImag = AVXImag + 1;
            idxOffset += 32;
            // FinalOut1 = FinalOut1 + 4;
            // FinalOut2 = FinalOut2 + 4;
            // FinalOut3 = FinalOut3 + 4;
            // FinalOut4 = FinalOut4 + 4;
        }

        rem = inpLen&0x0007;//mod(inpLen, 8)

        if(rem != 0)
        {
            absReal = _mm256_andnot_ps(SGN, *AVXReal);
            absImag = _mm256_andnot_ps(SGN, *AVXImag);
            RES1Real = _mm256_add_ps( _mm256_or_ps(SGN, absReal), Scalar1);
            RES1Imag = _mm256_add_ps( _mm256_or_ps(SGN, absImag), Scalar1);

            MEM1 = _mm256_unpacklo_ps(*AVXReal, *AVXImag);
            MEM2 = _mm256_unpacklo_ps(RES1Real, RES1Imag);
            MEM3 = _mm256_unpackhi_ps(*AVXReal, *AVXImag);
            MEM4 = _mm256_unpackhi_ps(RES1Real, RES1Imag);

            MEM11 = _mm256_shuffle_ps(MEM1, MEM2, 0b01000100);
            MEM21 = _mm256_shuffle_ps(MEM1, MEM2, 0b11101110);
            MEM31 = _mm256_shuffle_ps(MEM3, MEM4, 0b01000100);
            MEM41 = _mm256_shuffle_ps(MEM3, MEM4, 0b11101110);

            if(rem == 7)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM31, MEM41, 0x20);
                memcpy(&llrOut[idxOffset+8], tempStore_f, 32); 

                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x31);
                memcpy(&llrOut[idxOffset+16], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM31, MEM41, 0x31);
                memcpy(&llrOut[idxOffset+32], tempStore_f, 16);
            }
            else if(rem == 6)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM31, MEM41, 0x20);
                memcpy(&llrOut[idxOffset+8], tempStore_f, 32); 

                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x31);
                memcpy(&llrOut[idxOffset+16], tempStore_f, 32);
            }
            else if(rem == 5)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM31, MEM41, 0x20);
                memcpy(&llrOut[idxOffset+8], tempStore_f, 32); 

                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x31);
                memcpy(&llrOut[idxOffset+16], tempStore_f, 16);
            }
            else if(rem == 4)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM31, MEM41, 0x20);
                memcpy(&llrOut[idxOffset+8], tempStore_f, 32);
            }
            else if(rem == 3)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM31, MEM41, 0x20);
                memcpy(&llrOut[idxOffset+8], tempStore_f, 16);
            }
            else if(rem == 2)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
            }
            else
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 16);
            }
        }
    }
    else if(modOrder == 6)
    {
        __m256 Scalar2 = _mm256_set1_ps(0.308606);//2/sqrt(42)
        __m256 Scalar3 = _mm256_set1_ps(0.617213);//4/sqrt(42)
        __m256 Scalar4 = _mm256_set1_ps(0.925820);//6/sqrt(42)

        for(i = 0;i<(inpLen>>3);i++)
        {
            absReal = _mm256_andnot_ps(SGN, *AVXReal);
            absImag = _mm256_andnot_ps(SGN, *AVXImag);
            MEM1 = _mm256_or_ps(SGN, absReal);
            MEM2 = _mm256_or_ps(SGN, absImag);
            RES1Real = _mm256_add_ps(MEM1 , Scalar3);
            RES1Imag = _mm256_add_ps(MEM2 , Scalar3);

            Mask1Real = _mm256_cmp_ps(absReal, Scalar3, 1);
            Mask1Imag = _mm256_cmp_ps(absImag, Scalar3, 1);
            MEM3 = _mm256_sub_ps(absReal, Scalar2);
            MEM4 = _mm256_sub_ps(absImag, Scalar2);
            MEM5 = _mm256_add_ps(MEM1, Scalar4);
            MEM6 = _mm256_add_ps(MEM2, Scalar4);
            RES2Real = _mm256_blendv_ps(MEM5, MEM3, Mask1Real);
            RES2Imag = _mm256_blendv_ps(MEM6, MEM4, Mask1Imag);


            MEM1 = _mm256_unpacklo_ps(*AVXReal, *AVXImag);
            MEM2 = _mm256_unpackhi_ps(*AVXReal, *AVXImag);
            MEM3 = _mm256_unpacklo_ps(RES1Real, RES1Imag);
            MEM4 = _mm256_unpackhi_ps(RES1Real, RES1Imag);
            MEM5 = _mm256_unpacklo_ps(RES2Real, RES2Imag);
            MEM6 = _mm256_unpackhi_ps(RES2Real, RES2Imag);

            MEM11 = _mm256_shuffle_ps(MEM1, MEM3, 0b01000100);
            MEM21 = _mm256_shuffle_ps(MEM3, MEM5, 0b11101110);
            MEM31 = _mm256_shuffle_ps(MEM2, MEM4, 0b01000100);
            MEM41 = _mm256_shuffle_ps(MEM4, MEM6, 0b11101110);
            MEM51 = _mm256_shuffle_ps(MEM5, MEM1, 0b11100100);
            MEM61 = _mm256_shuffle_ps(MEM6, MEM2, 0b11100100);

            // *FinalOut1 = _mm256_permute2f128_ps(MEM11, MEM51, 0x20);
            // *FinalOut2 = _mm256_permute2f128_ps(MEM21, MEM31, 0x20);
            // *FinalOut3 = _mm256_permute2f128_ps(MEM61, MEM41, 0x20);
            // *FinalOut4 = _mm256_permute2f128_ps(MEM11, MEM51, 0x31);
            // *FinalOut5 = _mm256_permute2f128_ps(MEM21, MEM31, 0x31);
            // *FinalOut6 = _mm256_permute2f128_ps(MEM61, MEM41, 0x31);

            _mm256_storeu_ps(&llrOut[idxOffset], _mm256_permute2f128_ps(MEM11, MEM51, 0x20) );
            _mm256_storeu_ps(&llrOut[idxOffset+8], _mm256_permute2f128_ps(MEM21, MEM31, 0x20) );
            _mm256_storeu_ps(&llrOut[idxOffset+16], _mm256_permute2f128_ps(MEM61, MEM41, 0x20) );
            _mm256_storeu_ps(&llrOut[idxOffset+24], _mm256_permute2f128_ps(MEM11, MEM51, 0x31) );
            _mm256_storeu_ps(&llrOut[idxOffset+32], _mm256_permute2f128_ps(MEM21, MEM31, 0x31) );
            _mm256_storeu_ps(&llrOut[idxOffset+40], _mm256_permute2f128_ps(MEM61, MEM41, 0x31) );

            AVXReal += 1;
            AVXImag += 1;
            idxOffset += 48;
            // FinalOut1 += 6;
            // FinalOut2 += 6;
            // FinalOut3 += 6;
            // FinalOut4 += 6;
            // FinalOut5 += 6;
            // FinalOut6 += 6;
        }
        rem = inpLen&0x0007;//mod(inpLen, 8);

        if(rem != 0)
        {
            absReal = _mm256_andnot_ps(SGN, *AVXReal);
            absImag = _mm256_andnot_ps(SGN, *AVXImag);
            MEM1 = _mm256_or_ps(SGN, absReal);
            MEM2 = _mm256_or_ps(SGN, absImag);
            RES1Real = _mm256_add_ps(MEM1 , Scalar3);
            RES1Imag = _mm256_add_ps(MEM2 , Scalar3);

            Mask1Real = _mm256_cmp_ps(absReal, Scalar3, 1);
            Mask1Imag = _mm256_cmp_ps(absImag, Scalar3, 1);
            MEM3 = _mm256_sub_ps(absReal, Scalar2);
            MEM4 = _mm256_sub_ps(absImag, Scalar2);
            MEM5 = _mm256_add_ps(MEM1, Scalar4);
            MEM6 = _mm256_add_ps(MEM2, Scalar4);
            RES2Real = _mm256_blendv_ps(MEM5, MEM3, Mask1Real);
            RES2Imag = _mm256_blendv_ps(MEM6, MEM4, Mask1Imag);


            MEM1 = _mm256_unpacklo_ps(*AVXReal, *AVXImag);
            MEM2 = _mm256_unpackhi_ps(*AVXReal, *AVXImag);
            MEM3 = _mm256_unpacklo_ps(RES1Real, RES1Imag);
            MEM4 = _mm256_unpackhi_ps(RES1Real, RES1Imag);
            MEM5 = _mm256_unpacklo_ps(RES2Real, RES2Imag);
            MEM6 = _mm256_unpackhi_ps(RES2Real, RES2Imag);

            MEM11 = _mm256_shuffle_ps(MEM1, MEM3, 0b01000100);
            MEM21 = _mm256_shuffle_ps(MEM3, MEM5, 0b11101110);
            MEM31 = _mm256_shuffle_ps(MEM2, MEM4, 0b01000100);
            MEM41 = _mm256_shuffle_ps(MEM4, MEM6, 0b11101110);
            MEM51 = _mm256_shuffle_ps(MEM5, MEM1, 0b11100100);
            MEM61 = _mm256_shuffle_ps(MEM6, MEM2, 0b11100100);
            
            if(rem == 7)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM51, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM21, MEM31, 0x20);
                memcpy(&llrOut[idxOffset+8], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM61, MEM41, 0x20);
                memcpy(&llrOut[idxOffset+16], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM11, MEM51, 0x31);
                memcpy(&llrOut[idxOffset+24], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM21, MEM31, 0x31);
                memcpy(&llrOut[idxOffset+32], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM61, MEM41, 0x31);
                memcpy(&llrOut[idxOffset+40], tempStore_f, 8);
            }
            else if(rem == 6)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM51, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM21, MEM31, 0x20);
                memcpy(&llrOut[idxOffset+8], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM61, MEM41, 0x20);
                memcpy(&llrOut[idxOffset+16], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM11, MEM51, 0x31);
                memcpy(&llrOut[idxOffset+24], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM21, MEM31, 0x31);
                memcpy(&llrOut[idxOffset+32], tempStore_f, 16);
            }
            else if(rem == 5)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM51, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM21, MEM31, 0x20);
                memcpy(&llrOut[idxOffset+8], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM61, MEM41, 0x20);
                memcpy(&llrOut[idxOffset+16], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM11, MEM51, 0x31);
                memcpy(&llrOut[idxOffset+24], tempStore_f, 24);
            }
            else if(rem == 4)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM51, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM21, MEM31, 0x20);
                memcpy(&llrOut[idxOffset+8], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM61, MEM41, 0x20);
                memcpy(&llrOut[idxOffset+16], tempStore_f, 32);
            }
            else if(rem == 3)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM51, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM21, MEM31, 0x20);
                memcpy(&llrOut[idxOffset+8], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM61, MEM41, 0x20);
                memcpy(&llrOut[idxOffset+16], tempStore_f, 8);
            }
            else if(rem == 2)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM51, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM21, MEM31, 0x20);
                memcpy(&llrOut[idxOffset+8], tempStore_f, 16);
            }
            else
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM51, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 8);
            }
        }

    }
    else if(modOrder == 8)
    {
        __m256 Scalar5 = _mm256_set1_ps(0.6135719f);//8/sqrt(170)
        __m256 Scalar6 = _mm256_set1_ps(0.3067859f);//4/sqrt(170)
        __m256 Scalar7 = _mm256_set1_ps(0.920357986f);//12/sqrt(170)
        __m256 Scalar8 = _mm256_set1_ps(0.153392f);//2/sqrt(170)
        __m256 Scalar9 = _mm256_set1_ps(0.4601789f);//6/sqrt(170)
        __m256 Scalar10 = _mm256_set1_ps(0.7669649f);//10/sqrt(170)
        __m256 Scalar11 = _mm256_set1_ps(1.07375098f);//14/sqrt(170)

        __m256 Mask2Real, Mask2Imag, Mask3Real, Mask3Imag, RES3Real, RES3Imag;


        for(i = 0;i<(inpLen>>3);i++)
        {
            absReal = _mm256_andnot_ps(SGN, *AVXReal);
            absImag = _mm256_andnot_ps(SGN, *AVXImag);
            MEM1 = _mm256_or_ps(SGN, absReal);
            MEM2 = _mm256_or_ps(SGN, absImag);
            RES1Real = _mm256_add_ps(MEM1 , Scalar5);
            RES1Imag = _mm256_add_ps(MEM2 , Scalar5);

            Mask1Real = _mm256_cmp_ps(absReal, Scalar5, 1);
            Mask1Imag = _mm256_cmp_ps(absImag, Scalar5, 1);
            MEM3 = _mm256_sub_ps(absReal, Scalar6);
            MEM4 = _mm256_sub_ps(absImag, Scalar6);
            MEM5 = _mm256_add_ps(MEM1, Scalar7);
            MEM6 = _mm256_add_ps(MEM2, Scalar7);
            RES2Real = _mm256_blendv_ps(MEM5, MEM3, Mask1Real);
            RES2Imag = _mm256_blendv_ps(MEM6, MEM4, Mask1Imag);

            Mask2Real = _mm256_cmp_ps(absReal, Scalar6, 1);
            Mask2Imag = _mm256_cmp_ps(absImag, Scalar6, 1);
            MEM1 = _mm256_sub_ps(absReal, _mm256_and_ps(Mask2Real, Scalar8));
            MEM2 = _mm256_sub_ps(absImag, _mm256_and_ps(Mask2Imag, Scalar8));

            Mask1Real = _mm256_cmp_ps(Scalar6, MEM1, 2);
            Mask1Imag = _mm256_cmp_ps(Scalar6, MEM2, 2);
            Mask2Real = _mm256_cmp_ps(MEM1, Scalar5, 1);
            Mask2Imag = _mm256_cmp_ps(MEM2, Scalar5, 1);
            Mask3Real = _mm256_and_ps(Mask1Real, Mask2Real);
            Mask3Imag = _mm256_and_ps(Mask1Imag, Mask2Imag);

            MEM3 = _mm256_add_ps( _mm256_or_ps(MEM1 , _mm256_and_ps(SGN, Mask3Real)), _mm256_and_ps(Scalar9, Mask3Real));
            MEM4 = _mm256_add_ps( _mm256_or_ps(MEM2 , _mm256_and_ps(SGN, Mask3Imag)), _mm256_and_ps(Scalar9, Mask3Imag));

            Mask1Real = _mm256_cmp_ps(Scalar5, MEM3, 2);
            Mask1Imag = _mm256_cmp_ps(Scalar5, MEM4, 2);
            Mask2Real = _mm256_cmp_ps(MEM3, Scalar7, 1);
            Mask2Imag = _mm256_cmp_ps(MEM4, Scalar7, 1);
            Mask3Real = _mm256_and_ps(Mask1Real, Mask2Real);
            Mask3Imag = _mm256_and_ps(Mask1Imag, Mask2Imag);

            MEM1 = _mm256_sub_ps( MEM3, _mm256_and_ps(Scalar10, Mask3Real));
            MEM2 = _mm256_sub_ps( MEM4, _mm256_and_ps(Scalar10, Mask3Imag));

            Mask1Real = _mm256_cmp_ps(Scalar7, MEM1, 2);
            Mask1Imag = _mm256_cmp_ps(Scalar7, MEM2, 2);
            RES3Real =  _mm256_add_ps( _mm256_or_ps(MEM1 , _mm256_and_ps(SGN, Mask1Real)), _mm256_and_ps(Scalar11, Mask1Real));//d
            RES3Imag =  _mm256_add_ps( _mm256_or_ps(MEM2 , _mm256_and_ps(SGN, Mask1Imag)), _mm256_and_ps(Scalar11, Mask1Imag));


            MEM1 = _mm256_unpacklo_ps(*AVXReal, *AVXImag);
            MEM2 = _mm256_unpackhi_ps(*AVXReal, *AVXImag);
            MEM3 = _mm256_unpacklo_ps(RES1Real, RES1Imag);
            MEM4 = _mm256_unpackhi_ps(RES1Real, RES1Imag);
            MEM5 = _mm256_unpacklo_ps(RES2Real, RES2Imag);
            MEM6 = _mm256_unpackhi_ps(RES2Real, RES2Imag);
            MEM7 = _mm256_unpacklo_ps(RES3Real, RES3Imag);
            MEM8 = _mm256_unpackhi_ps(RES3Real, RES3Imag);

            MEM11 = _mm256_shuffle_ps(MEM1, MEM3, 0b01000100);
            MEM21 = _mm256_shuffle_ps(MEM5, MEM7, 0b01000100);
            MEM31 = _mm256_shuffle_ps(MEM1, MEM3, 0b11101110);
            MEM41 = _mm256_shuffle_ps(MEM5, MEM7, 0b11101110);
            MEM51 = _mm256_shuffle_ps(MEM2, MEM4, 0b01000100);
            MEM61 = _mm256_shuffle_ps(MEM6, MEM8, 0b01000100);
            MEM71 = _mm256_shuffle_ps(MEM2, MEM4, 0b11101110);
            MEM81 = _mm256_shuffle_ps(MEM6, MEM8, 0b11101110);

            // *FinalOut1 = _mm256_permute2f128_ps(MEM11, MEM21, 0x20);
            // *FinalOut2 = _mm256_permute2f128_ps(MEM31, MEM41, 0x20);
            // *FinalOut3 = _mm256_permute2f128_ps(MEM51, MEM61, 0x20);
            // *FinalOut4 = _mm256_permute2f128_ps(MEM71, MEM81, 0x20);
            // *FinalOut5 = _mm256_permute2f128_ps(MEM11, MEM21, 0x31);
            // *FinalOut6 = _mm256_permute2f128_ps(MEM31, MEM41, 0x31);
            // *FinalOut7 = _mm256_permute2f128_ps(MEM51, MEM61, 0x31);
            // *FinalOut8 = _mm256_permute2f128_ps(MEM71, MEM81, 0x31);

            _mm256_storeu_ps(&llrOut[idxOffset], _mm256_permute2f128_ps(MEM11, MEM21, 0x20) );
            _mm256_storeu_ps(&llrOut[idxOffset+8], _mm256_permute2f128_ps(MEM31, MEM41, 0x20) );
            _mm256_storeu_ps(&llrOut[idxOffset+16], _mm256_permute2f128_ps(MEM51, MEM61, 0x20) );
            _mm256_storeu_ps(&llrOut[idxOffset+24], _mm256_permute2f128_ps(MEM71, MEM81, 0x20) );
            _mm256_storeu_ps(&llrOut[idxOffset+32], _mm256_permute2f128_ps(MEM11, MEM21, 0x31) );
            _mm256_storeu_ps(&llrOut[idxOffset+40], _mm256_permute2f128_ps(MEM31, MEM41, 0x31) );
            _mm256_storeu_ps(&llrOut[idxOffset+48], _mm256_permute2f128_ps(MEM51, MEM61, 0x31) );
            _mm256_storeu_ps(&llrOut[idxOffset+56], _mm256_permute2f128_ps(MEM71, MEM81, 0x31) );


            AVXReal += 1;
            AVXImag += 1;
            idxOffset += 64;

            // FinalOut1 += 8;
            // FinalOut2 += 8;
            // FinalOut3 += 8;
            // FinalOut4 += 8;
            // FinalOut5 += 8;
            // FinalOut6 += 8;
            // FinalOut7 += 8;
            // FinalOut8 += 8;
        }
        rem = inpLen&0x0007;//mod(inpLen, 8);

        if(rem != 0)
        {
            absReal = _mm256_andnot_ps(SGN, *AVXReal);
            absImag = _mm256_andnot_ps(SGN, *AVXImag);
            MEM1 = _mm256_or_ps(SGN, absReal);
            MEM2 = _mm256_or_ps(SGN, absImag);
            RES1Real = _mm256_add_ps(MEM1 , Scalar5);
            RES1Imag = _mm256_add_ps(MEM2 , Scalar5);

            Mask1Real = _mm256_cmp_ps(absReal, Scalar5, 1);
            Mask1Imag = _mm256_cmp_ps(absImag, Scalar5, 1);
            MEM3 = _mm256_sub_ps(absReal, Scalar6);
            MEM4 = _mm256_sub_ps(absImag, Scalar6);
            MEM5 = _mm256_add_ps(MEM1, Scalar7);
            MEM6 = _mm256_add_ps(MEM2, Scalar7);
            RES2Real = _mm256_blendv_ps(MEM5, MEM3, Mask1Real);
            RES2Imag = _mm256_blendv_ps(MEM6, MEM4, Mask1Imag);

            Mask2Real = _mm256_cmp_ps(absReal, Scalar6, 1);
            Mask2Imag = _mm256_cmp_ps(absImag, Scalar6, 1);
            MEM1 = _mm256_sub_ps(absReal, _mm256_and_ps(Mask2Real, Scalar8));
            MEM2 = _mm256_sub_ps(absImag, _mm256_and_ps(Mask2Imag, Scalar8));

            Mask1Real = _mm256_cmp_ps(Scalar6, MEM1, 2);
            Mask1Imag = _mm256_cmp_ps(Scalar6, MEM2, 2);
            Mask2Real = _mm256_cmp_ps(MEM1, Scalar5, 1);
            Mask2Imag = _mm256_cmp_ps(MEM2, Scalar5, 1);
            Mask3Real = _mm256_and_ps(Mask1Real, Mask2Real);
            Mask3Imag = _mm256_and_ps(Mask1Imag, Mask2Imag);

            MEM3 = _mm256_add_ps( _mm256_or_ps(MEM1 , _mm256_and_ps(SGN, Mask3Real)), _mm256_and_ps(Scalar9, Mask3Real));
            MEM4 = _mm256_add_ps( _mm256_or_ps(MEM2 , _mm256_and_ps(SGN, Mask3Imag)), _mm256_and_ps(Scalar9, Mask3Imag));

            Mask1Real = _mm256_cmp_ps(Scalar5, MEM3, 2);
            Mask1Imag = _mm256_cmp_ps(Scalar5, MEM4, 2);
            Mask2Real = _mm256_cmp_ps(MEM3, Scalar7, 1);
            Mask2Imag = _mm256_cmp_ps(MEM4, Scalar7, 1);
            Mask3Real = _mm256_and_ps(Mask1Real, Mask2Real);
            Mask3Imag = _mm256_and_ps(Mask1Imag, Mask2Imag);

            MEM1 = _mm256_sub_ps( MEM3, _mm256_and_ps(Scalar10, Mask3Real));
            MEM2 = _mm256_sub_ps( MEM4, _mm256_and_ps(Scalar10, Mask3Imag));

            Mask1Real = _mm256_cmp_ps(Scalar7, MEM1, 2);
            Mask1Imag = _mm256_cmp_ps(Scalar7, MEM2, 2);
            RES3Real =  _mm256_add_ps( _mm256_or_ps(MEM1 , _mm256_and_ps(SGN, Mask1Real)), _mm256_and_ps(Scalar11, Mask1Real));//d
            RES3Imag =  _mm256_add_ps( _mm256_or_ps(MEM2 , _mm256_and_ps(SGN, Mask1Imag)), _mm256_and_ps(Scalar11, Mask1Imag));


            MEM1 = _mm256_unpacklo_ps(*AVXReal, *AVXImag);
            MEM2 = _mm256_unpackhi_ps(*AVXReal, *AVXImag);
            MEM3 = _mm256_unpacklo_ps(RES1Real, RES1Imag);
            MEM4 = _mm256_unpackhi_ps(RES1Real, RES1Imag);
            MEM5 = _mm256_unpacklo_ps(RES2Real, RES2Imag);
            MEM6 = _mm256_unpackhi_ps(RES2Real, RES2Imag);
            MEM7 = _mm256_unpacklo_ps(RES3Real, RES3Imag);
            MEM8 = _mm256_unpackhi_ps(RES3Real, RES3Imag);

            MEM11 = _mm256_shuffle_ps(MEM1, MEM3, 0b01000100);
            MEM21 = _mm256_shuffle_ps(MEM5, MEM7, 0b01000100);
            MEM31 = _mm256_shuffle_ps(MEM1, MEM3, 0b11101110);
            MEM41 = _mm256_shuffle_ps(MEM5, MEM7, 0b11101110);
            MEM51 = _mm256_shuffle_ps(MEM2, MEM4, 0b01000100);
            MEM61 = _mm256_shuffle_ps(MEM6, MEM8, 0b01000100);
            MEM71 = _mm256_shuffle_ps(MEM2, MEM4, 0b11101110);
            MEM81 = _mm256_shuffle_ps(MEM6, MEM8, 0b11101110);

            if(rem == 7)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM31, MEM41, 0x20);
                memcpy(&llrOut[idxOffset+8], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM51, MEM61, 0x20);
                memcpy(&llrOut[idxOffset+16], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM71, MEM81, 0x20);
                memcpy(&llrOut[idxOffset+24], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x31);
                memcpy(&llrOut[idxOffset+32], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM31, MEM41, 0x31);
                memcpy(&llrOut[idxOffset+40], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM51, MEM61, 0x31);
                memcpy(&llrOut[idxOffset+48], tempStore_f, 32);
            }
            else if(rem == 6)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM31, MEM41, 0x20);
                memcpy(&llrOut[idxOffset+8], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM51, MEM61, 0x20);
                memcpy(&llrOut[idxOffset+16], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM71, MEM81, 0x20);
                memcpy(&llrOut[idxOffset+24], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x31);
                memcpy(&llrOut[idxOffset+32], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM31, MEM41, 0x31);
                memcpy(&llrOut[idxOffset+40], tempStore_f, 32);
            }
            else if(rem == 5)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM31, MEM41, 0x20);
                memcpy(&llrOut[idxOffset+8], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM51, MEM61, 0x20);
                memcpy(&llrOut[idxOffset+16], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM71, MEM81, 0x20);
                memcpy(&llrOut[idxOffset+24], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x31);
                memcpy(&llrOut[idxOffset+32], tempStore_f, 32);
            }
            else if(rem == 4)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM31, MEM41, 0x20);
                memcpy(&llrOut[idxOffset+8], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM51, MEM61, 0x20);
                memcpy(&llrOut[idxOffset+16], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM71, MEM81, 0x20);
                memcpy(&llrOut[idxOffset+24], tempStore_f, 32);
            }
            else if(rem == 3)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM31, MEM41, 0x20);
                memcpy(&llrOut[idxOffset+8], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM51, MEM61, 0x20);
                memcpy(&llrOut[idxOffset+16], tempStore_f, 32);
            }
            else if(rem == 2)
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
                tempStore = _mm256_permute2f128_ps(MEM31, MEM41, 0x20);
                memcpy(&llrOut[idxOffset+8], tempStore_f, 32);
            }
            else
            {
                tempStore = _mm256_permute2f128_ps(MEM11, MEM21, 0x20);
                memcpy(&llrOut[idxOffset], tempStore_f, 32);
            }
        }
    }
}



void llrScalingQuant_24(int symPart, unsigned char nsymb, int nPRB, unsigned short startToneIdx, char modOrder,\
                     char nLayers, char bwl,  float mseRealLow[][12*MAX_NUM_PRB], float *inp, int8_t *out)
{
    uint32_t prbIdx = 0;

    __m256* softOut     = (__m256*)inp;
    __m256i* scaleOut    = (__m256i*)out;

    int8_t clip = (1<<(bwl-1))-1;

    __m256 mseVec;
    int8_t* scaleOut8;
    __m256i inter1, inter2, inter3, order = _mm256_setr_epi32 (0, 4, 1, 5, 2, 6, 3, 7);
    __m256i minVal = _mm256_set1_epi8 (-1*clip), maxVal = _mm256_set1_epi8 (clip), zeros = _mm256_setzero_si256 ();;
    // Layers = 1
    if (nLayers == 1)
    {
        float* mseValues = &mseRealLow[0][0];

        // QPSK
        if (modOrder == 2)
        {
            for (uint8_t symIdx = 0; symIdx < nsymb; symIdx++)
            {
                // There are 2 PRBs with same MSE
                // So, we need to work on 2 PRBs*12 Tones*2 Mod Order = 48 LLRs
                // AVX2 works on 32 at a time. => We need 1.5 AVX loops to work on => Generates Edge effects
                // So, we work on 4 PRBs at a time to increase efficiency or restrict edge effects
                for (prbIdx = 0; (prbIdx + 4) <= nPRB; prbIdx+=4)
                {
                    mseVec = _mm256_set1_ps (mseValues[prbIdx>>1]);

                    // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                    inter3 = _mm256_packs_epi32(inter1, inter2);

                    // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                    inter2 = _mm256_packs_epi32(inter1, inter2);
                    
                    // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                    // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                    // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                    inter1 = _mm256_packs_epi16(inter3, inter2);

                    // The above order is interleaved. So, we put them back in the correct format
                    // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                    // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                    inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                    // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                    // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                    inter1 = _mm256_max_epi8(inter1, minVal);
                    inter1 = _mm256_min_epi8(inter1, maxVal);

                    // Storing the Data in the output array
                    _mm256_storeu_si256 (scaleOut, inter1);

                    scaleOut++;


                    // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                    inter3 = _mm256_packs_epi32(inter1, inter2);

                    mseVec = _mm256_set1_ps (mseValues[(prbIdx>>1)+1]);

                    // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                    inter2 = _mm256_packs_epi32(inter1, inter2);
                    
                    // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                    // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                    // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                    inter1 = _mm256_packs_epi16(inter3, inter2);

                    // The above order is interleaved. So, we put them back in the correct format
                    // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                    // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                    inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                    // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                    // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                    inter1 = _mm256_max_epi8(inter1, minVal);
                    inter1 = _mm256_min_epi8(inter1, maxVal);

                    // Storing the Data in the output array
                    _mm256_storeu_si256 (scaleOut, inter1);

                    scaleOut++;

                    // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                    inter3 = _mm256_packs_epi32(inter1, inter2);

                    // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                    inter2 = _mm256_packs_epi32(inter1, inter2);
                    
                    // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                    // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                    // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                    inter1 = _mm256_packs_epi16(inter3, inter2);

                    // The above order is interleaved. So, we put them back in the correct format
                    // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                    // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                    inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                    // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                    // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                    inter1 = _mm256_max_epi8(inter1, minVal);
                    inter1 = _mm256_min_epi8(inter1, maxVal);

                    // Storing the Data in the output array
                    _mm256_storeu_si256 (scaleOut, inter1);

                    scaleOut++;
                }

                // When we work on 4 PRBs at a time, we have edge effects of <4 PRBs remaining.
                uint8_t excessPrbs = nPRB - prbIdx;

                // printf("excessPrbs: %u\n", excessPrbs);
                // In case, we have 1 PRB remaining.
                // 1 PRB*12 Tones*2 Mod Order = 24 LLRs
                // AVX works on 32 LLRs at a time.
                // So, we place last 8 LLRs as 0s in this case
                if (excessPrbs == 1)
                {
                    scaleOut8 = (int8_t*)scaleOut;

                    mseVec = _mm256_set1_ps (mseValues[prbIdx>>1]);

                    // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                    inter3 = _mm256_packs_epi32(inter1, inter2);

                    // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                    inter2 = _mm256_packs_epi32(inter1, zeros);
                    
                    // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                    // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                    // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                    inter1 = _mm256_packs_epi16(inter3, inter2);

                    // The above order is interleaved. So, we put them back in the correct format
                    // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                    // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                    inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                    // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                    // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                    inter1 = _mm256_max_epi8(inter1, minVal);
                    inter1 = _mm256_min_epi8(inter1, maxVal);

                    // Storing the Data in the output array
                    memcpy (scaleOut8, &inter1, 24);

                    scaleOut8 += 24;

                    scaleOut = (__m256i*)scaleOut8;
                }

                // In case, we have 2 PRB remaining.
                // 2 PRB*12 Tones*2 Mod Order = 48 LLRs
                // AVX works on 32 LLRs at a time. So, we have 1.5 Iterations
                else if (excessPrbs == 2)
                {
                    // We divide the 48 LLR Processing into 2 parts.

                    // 1. We deal with 32 LLRs as follows
                    mseVec = _mm256_set1_ps (mseValues[prbIdx>>1]);

                    // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                    inter3 = _mm256_packs_epi32(inter1, inter2);

                    // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                    inter2 = _mm256_packs_epi32(inter1, inter2);
                    
                    // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                    // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                    // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                    inter1 = _mm256_packs_epi16(inter3, inter2);

                    // The above order is interleaved. So, we put them back in the correct format
                    // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                    // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                    inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                    // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                    // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                    inter1 = _mm256_max_epi8(inter1, minVal);
                    inter1 = _mm256_min_epi8(inter1, maxVal);

                    // Storing the Data in the output array
                    _mm256_storeu_si256 (scaleOut, inter1);

                    scaleOut++;

                    // 2. Once the 32 LLR Processing is done, we work on remaining 16 LLRs
                    // As AVX wors on 32 LLRs, we append 16 0s to make it a perfect vector
                    scaleOut8 = (int8_t*)scaleOut;

                    // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                    inter3 = _mm256_packs_epi32(inter1, inter2);

                    inter1 = _mm256_packs_epi16(inter3, zeros);

                    // The above order is interleaved. So, we put them back in the correct format
                    // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                    // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                    inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                    // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                    // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                    inter1 = _mm256_max_epi8(inter1, minVal);
                    inter1 = _mm256_min_epi8(inter1, maxVal);

                    // Storing the Data in the output array
                    memcpy (scaleOut8, &inter1, 16);

                    scaleOut8 += 16;

                    scaleOut = (__m256i *)scaleOut8;
                }

                // In case, we have 3 PRB remaining.
                // 3 PRB*12 Tones*2 Mod Order = 72 LLRs
                // AVX works on 32 LLRs at a time. So, we have 2.25 Iterations
                // But after the 1st 2 PRBs, i.e., 1.5 AVX Iterations, we have to change the MSE for
                // the 3rd remaining PRB.
                // The whole process is done in 3 steps
                else if (excessPrbs == 3)
                {
                    // 1. We process the 1st 2 PRBs
                    // Here 32 LLRs from 48 LLRs of 1st 2 PRBs are processed
                    mseVec = _mm256_set1_ps (mseValues[prbIdx>>1]);

                    // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                    inter3 = _mm256_packs_epi32(inter1, inter2);

                    // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                    inter2 = _mm256_packs_epi32(inter1, inter2);
                    
                    // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                    // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                    // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                    inter1 = _mm256_packs_epi16(inter3, inter2);

                    // The above order is interleaved. So, we put them back in the correct format
                    // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                    // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                    inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                    // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                    // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                    inter1 = _mm256_max_epi8(inter1, minVal);
                    inter1 = _mm256_min_epi8(inter1, maxVal);

                    // Storing the Data in the output array
                    _mm256_storeu_si256 (scaleOut, inter1);

                    scaleOut++;

                    // 2. Now, the 16 remaining LLRs of 1st 2 PRBs and 16 LLRs of 3rd PRB are
                    // managed here.

                    // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                    inter3 = _mm256_packs_epi32(inter1, inter2);

                    mseVec = _mm256_set1_ps (mseValues[(prbIdx>>1)+1]);

                    // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                    inter2 = _mm256_packs_epi32(inter1, inter2);
                    
                    // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                    // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                    // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                    inter1 = _mm256_packs_epi16(inter3, inter2);

                    // The above order is interleaved. So, we put them back in the correct format
                    // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                    // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                    inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                    // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                    // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                    inter1 = _mm256_max_epi8(inter1, minVal);
                    inter1 = _mm256_min_epi8(inter1, maxVal);

                    // Storing the Data in the output array
                    _mm256_storeu_si256 (scaleOut, inter1);

                    scaleOut++;

                    // 3. 3rd PRB has 1 PRB*12 Tones*2 Mod Order = 24 LLRs
                    // Out of 24 LLRs, we have processed the 1st 16 in the previous step(2)
                    // Now we shall process the last 8 LLRs, by appending 24 zros to fit the AVX vector
                    scaleOut8 = (int8_t*)scaleOut;

                    // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                    inter3 = _mm256_packs_epi32(inter1, zeros);
                    
                    // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                    // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                    // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                    inter1 = _mm256_packs_epi16(inter3, zeros);

                    // The above order is interleaved. So, we put them back in the correct format
                    // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                    // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                    inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                    // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                    // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                    inter1 = _mm256_max_epi8(inter1, minVal);
                    inter1 = _mm256_min_epi8(inter1, maxVal);

                    // Storing the Data in the output array
                    memcpy (scaleOut8, &inter1, 8);

                    scaleOut8 += 8;

                    scaleOut = (__m256i *)scaleOut8;
                }
            }
        }
        // 16 QAM
        else if (modOrder == 4)
        {
            // printf("modOrder: %d\n", modOrder);
            for (uint8_t symIdx = 0; symIdx < nsymb; symIdx++)
            {
                // There are 2 PRBs with same MSE
                // So, we need to work on 2 PRBs*12 Tones*4 Mod Order = 96 LLRs
                // AVX2 works on 32 at a time. => We need 3 AVX loops to work on => No Edge effects
                // So, we work on 2 PRBs at a time to increase efficiency or restrict edge effects
                for (prbIdx = 0; (prbIdx + 2) <= nPRB; prbIdx+=2)
                {
                    // Fix the MSE Vector used for the next 3 AVX iterations
                    mseVec = _mm256_set1_ps (mseValues[prbIdx>>1]);

                    // We perform the 3 AVX Iterations here
                    for (uint8_t idx = 0; idx < 3; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }

                // When we work on 2 PRBs at a time, we have edge effects of <2 PRBs remaining.
                uint8_t excessPrbs = nPRB - prbIdx;

                // In case, we have 1 PRB remaining.
                // 1 PRB*12 Tones*4 Mod Order = 48 LLRs
                // AVX operates on 32 LLRs at a time. Thus, we require 1.5 Iterations
                // We shall perform the whole process in 2 steps.
                if (excessPrbs)
                {
                    // 1. We shall perform 1st 32 LLR processing normally
                    mseVec = _mm256_set1_ps (mseValues[prbIdx>>1]);

                    // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                    inter3 = _mm256_packs_epi32(inter1, inter2);

                    // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                    inter2 = _mm256_packs_epi32(inter1, inter2);
                    
                    // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                    // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                    // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                    inter1 = _mm256_packs_epi16(inter3, inter2);

                    // The above order is interleaved. So, we put them back in the correct format
                    // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                    // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                    inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                    // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                    // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                    inter1 = _mm256_max_epi8(inter1, minVal);
                    inter1 = _mm256_min_epi8(inter1, maxVal);

                    // Storing the Data in the output array
                    _mm256_storeu_si256 (scaleOut, inter1);

                    scaleOut++;

                    // 2. We shall perform the same operation as above on last 16 LLRs, by appending
                    // 16 0s in the AVX Vector
                    scaleOut8 = (int8_t*)scaleOut;

                    // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                    inter3 = _mm256_packs_epi32(inter1, inter2);

                    inter1 = _mm256_packs_epi16(inter3, zeros);

                    // The above order is interleaved. So, we put them back in the correct format
                    // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                    // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                    inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                    // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                    // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                    inter1 = _mm256_max_epi8(inter1, minVal);
                    inter1 = _mm256_min_epi8(inter1, maxVal);

                    // Storing the Data in the output array
                    memcpy (scaleOut8, &inter1, 16);

                    scaleOut8 += 16;

                    scaleOut = (__m256i *)scaleOut8;
                }
            }
        }
        // 64 QAM
        else if (modOrder == 6)
        {
            for (uint8_t symIdx = 0; symIdx < nsymb; symIdx++)
            {
                // There are 2 PRBs with same MSE
                // So, we need to work on 2 PRBs*12 Tones*6 Mod Order = 144 LLRs
                // AVX2 works on 32 at a time. => We need 4.5 AVX loops to work on => Generates Edge effects
                // So, we work on 4 PRBs at a time to increase efficiency or restrict edge effects
                // This is done in 3 steps
                for (prbIdx = 0; (prbIdx + 4) <= nPRB; prbIdx+=4)
                {
                    // 1. We shall perform AVX operation on 1st 128 LLRs which
                    // are the part of 1st 2 PRBs, using 4 AVX Iterations as follows
                    mseVec = _mm256_set1_ps (mseValues[prbIdx>>1]);

                    for (uint8_t idx = 0; idx < 4; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }

                    // 2. We perform the 2nd step by placing remaining 16 LLRs of 1st 2 PRBs and
                    // appending them with 1st 16 LLRs of next 2 PRBsmaking it a perfect vector as follows
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        mseVec = _mm256_set1_ps (mseValues[(prbIdx>>1)+1]);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }

                    // 3. We shall perform AVX operation on remaining 128 LLRs which
                    // are the part of 2nd set of 2 PRBs, using 4 AVX Iterations as follows
                    // and complete the operation
                    for (uint8_t idx = 0; idx < 4; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }

                // When we work on 4 PRBs at a time, we have edge effects of <4 PRBs remaining.
                uint8_t excessPrbs = nPRB - prbIdx;


                // In case, we have 1 PRB remaining.
                // 1 PRB*12 Tones*6 Mod Order = 72 LLRs
                // AVX works on 32 LLRs at a time. Hence we need to perform 2.25 AVX Iterations
                // So, we do this operation in 2 steps
                if (excessPrbs == 1)
                {
                    // 1. We shall work on 1st 64 LLRs in normal way which we are following in above steps
                    mseVec = _mm256_set1_ps (mseValues[prbIdx>>1]);

                    for (uint8_t idx = 0; idx < 2; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }

                    // 2. We work on remaining 8 LLRs, by padding 24 0s to make it a perfect AVX Vector
                    {
                        scaleOut8 = (int8_t*)scaleOut;

                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, zeros);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, zeros);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        memcpy (scaleOut8, &inter1, 8);

                        scaleOut8 += 8;

                        scaleOut = (__m256i *)scaleOut8;
                    }
                }

                // In case, we have 2 PRB remaining.
                // 2 PRB*12 Tones*6 Mod Order = 144 LLRs
                // AVX works on 32 LLRs at a time. Hence we need to perform 4.5 AVX Iterations
                // So, we do this operation in 2 steps, as no change in MSE
                else if (excessPrbs == 2)
                {
                    // 1. We work on 1st 128 LLRs in normal way
                    mseVec = _mm256_set1_ps (mseValues[prbIdx>>1]);

                    for (uint8_t idx = 0; idx < 4; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }

                    // 2. We shall work on remaining 16 LLRs by apending 0s at the end
                    {
                        scaleOut8 = (int8_t*)scaleOut;

                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        inter1 = _mm256_packs_epi16(inter3, zeros);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        memcpy (scaleOut8, &inter1, 16);

                        scaleOut8 += 16;

                        scaleOut = (__m256i *)scaleOut8;
                    }
                }

                // In case, we have 3 PRB remaining.
                // 3 PRB*12 Tones*6 Mod Order = 216 LLRs
                // AVX works on 32 LLRs at a time. Hence we need to perform 6.75 AVX Iterations
                // So, we do this operation in 4 steps, as there is a change in MSE after 2 PRBs
                else if (excessPrbs == 3)
                {
                    // Has 3 PRBs*12 Tones*6 Mod Order = 216 LLRs
                    // Each Iteration we manage 32 
                    mseVec = _mm256_set1_ps (mseValues[prbIdx>>1]);

                    // 1. 1st 2 PRBs*12 Tones*6 Mod Order = 144 LLRs have same MSE
                    // Among them we manage 1st 128 LLRs by the below loop,
                    // where we manage 32 LLRs per iteration and with 4 iterations,
                    // we complete 128 LLRs
                    for (uint8_t idx = 0; idx < 4; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }

                    // 2. This part covers the remaining 16 LLRs of 1st 2 PRBs and 1st 12 LLRs of 3rd PRB
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        mseVec = _mm256_set1_ps (mseValues[(prbIdx>>1)+1]);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }

                    // 3. 3rd PRB has 12 Tones*6 Mod Order = 72 LLRs
                    // Among which 6 were done in the previous loop
                    // Among the remaining 72-16 = 56 LLRs, we manage 32 here in below loop
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }

                    // 4. 3rd PRB has 56-32=24 LLRs remaing which are managed here
                    {
                        scaleOut8 = (int8_t*)scaleOut;

                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, zeros);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        memcpy (scaleOut8, &inter1, 24);

                        scaleOut8 += 24;

                        scaleOut = (__m256i *)scaleOut8;
                    }
                }       
            }
        }
        // 256 QAM
        else if (modOrder == 8)
        {
            for (uint8_t symIdx = 0; symIdx < nsymb; symIdx++)
            {
                // There are 2 PRBs with same MSE
                // So, we need to work on 2 PRBs*12 Tones*8 Mod Order = 192 LLRs
                // AVX2 works on 32 at a time. => We need 6 AVX loops to work on => No Edge effects
                for (prbIdx = 0; (prbIdx + 2) <= nPRB; prbIdx+=2)
                {
                    mseVec = _mm256_set1_ps (mseValues[prbIdx>>1]);

                    for (uint8_t idx = 0; idx < 6; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }

                // When we work on 2 PRBs, we have an edge case of 1 PRB remaining,
                // which is checked here
                uint8_t excessPrbs = nPRB - prbIdx;

                // 1 PRB has 1 PRB*12 Tones*8 Mod Order = 96 LLRs
                // AVX2 performs 32 LLR Processing, so we need 3 AVX Iterations
                if (excessPrbs)
                {
                    mseVec = _mm256_set1_ps (mseValues[prbIdx>>1]);

                    for (uint8_t idx = 0; idx < 3; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }
            }
        }
    }
    // Layers = 2
    else if (nLayers == 2)
    {
        // printf("nLayers: %d\n", nLayers);
        // QPSK
        if (modOrder == 2)
        {
            // printf("mseRealLow[1][0]: %f\n", mseRealLow[1][0]);
            for (uint8_t symIdx = 0; symIdx < nsymb; symIdx++)
            {
                for (prbIdx = 0; (prbIdx + 2) <= nPRB; prbIdx+=2)
                {
                    mseVec = _mm256_setr_ps (mseRealLow[0][prbIdx>>1], mseRealLow[0][prbIdx>>1],
                                             mseRealLow[1][prbIdx>>1], mseRealLow[1][prbIdx>>1],
                                             mseRealLow[0][prbIdx>>1], mseRealLow[0][prbIdx>>1],
                                             mseRealLow[1][prbIdx>>1], mseRealLow[1][prbIdx>>1]);

                    for (uint8_t idx = 0; idx < 3; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }

                uint8_t excessPrbs = nPRB - prbIdx;

                if (excessPrbs)
                {
                    mseVec = _mm256_setr_ps (mseRealLow[0][prbIdx>>1], mseRealLow[0][prbIdx>>1],
                                             mseRealLow[1][prbIdx>>1], mseRealLow[1][prbIdx>>1],
                                             mseRealLow[0][prbIdx>>1], mseRealLow[0][prbIdx>>1],
                                             mseRealLow[1][prbIdx>>1], mseRealLow[1][prbIdx>>1]);

                    // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                    inter3 = _mm256_packs_epi32(inter1, inter2);

                    // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                    inter2 = _mm256_packs_epi32(inter1, inter2);
                    
                    // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                    // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                    // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                    inter1 = _mm256_packs_epi16(inter3, inter2);

                    // The above order is interleaved. So, we put them back in the correct format
                    // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                    // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                    inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                    // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                    // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                    inter1 = _mm256_max_epi8(inter1, minVal);
                    inter1 = _mm256_min_epi8(inter1, maxVal);

                    // Storing the Data in the output array
                    _mm256_storeu_si256 (scaleOut, inter1);

                    scaleOut++;

                    scaleOut8 = (int8_t*)scaleOut;

                    // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                    inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                    inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                    softOut++;

                    // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                    // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                    inter3 = _mm256_packs_epi32(inter1, inter2);

                    inter1 = _mm256_packs_epi16(inter3, zeros);

                    // The above order is interleaved. So, we put them back in the correct format
                    // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                    // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                    inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                    // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                    // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                    inter1 = _mm256_max_epi8(inter1, minVal);
                    inter1 = _mm256_min_epi8(inter1, maxVal);

                    // Storing the Data in the output array
                    memcpy (scaleOut8, &inter1, 16);

                    scaleOut8 += 16;

                    scaleOut = (__m256i *)scaleOut8;
                }
            }
        }
        // 16 QAM
        else if (modOrder == 4)
        {
            for (uint8_t symIdx = 0; symIdx < nsymb; symIdx++)
            {
                for (prbIdx = 0; (prbIdx + 2) <= nPRB; prbIdx+=2)
                {
                    mseVec = _mm256_setr_ps (mseRealLow[0][prbIdx>>1], mseRealLow[0][prbIdx>>1],
                                             mseRealLow[0][prbIdx>>1], mseRealLow[0][prbIdx>>1],
                                             mseRealLow[1][prbIdx>>1], mseRealLow[1][prbIdx>>1],
                                             mseRealLow[1][prbIdx>>1], mseRealLow[1][prbIdx>>1]);

                    for (uint8_t idx = 0; idx < 6; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }

                uint8_t excessPrbs = nPRB - prbIdx;

                if (excessPrbs)
                {
                    mseVec = _mm256_setr_ps (mseRealLow[0][prbIdx>>1], mseRealLow[0][prbIdx>>1],
                                             mseRealLow[0][prbIdx>>1], mseRealLow[0][prbIdx>>1],
                                             mseRealLow[1][prbIdx>>1], mseRealLow[1][prbIdx>>1],
                                             mseRealLow[1][prbIdx>>1], mseRealLow[1][prbIdx>>1]);

                    for (uint8_t idx = 0; idx < 3; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }
            }
        }
        // 64 QAM
        else if (modOrder == 6)
        {
            __m256 mseScaleVec[3];
            for (uint8_t symIdx = 0; symIdx < nsymb; symIdx++)
            {
                for (prbIdx = 0; (prbIdx + 2) <= nPRB; prbIdx+=2)
                {
                    float a = mseRealLow[0][prbIdx>>1], b = mseRealLow[1][prbIdx>>1];

                    mseScaleVec[0] = _mm256_setr_ps (a, a, a, a,
                                                     a, a, b, b);

                    mseScaleVec[1] = _mm256_setr_ps (b, b, b, b,
                                                     a, a, a, a);

                    mseScaleVec[2] = _mm256_setr_ps (a, a, b, b,
                                                     b, b, b, b);

                    uint8_t inc = 0;

                    for (uint8_t idx = 0; idx < 9; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }

                uint8_t excessPrbs = nPRB - prbIdx;

                if (excessPrbs)
                {
                    float a = mseRealLow[0][prbIdx>>1], b = mseRealLow[1][prbIdx>>1];

                    mseScaleVec[0] = _mm256_setr_ps (a, a, a, a,
                                                     a, a, b, b);

                    mseScaleVec[1] = _mm256_setr_ps (b, b, b, b,
                                                     a, a, a, a);

                    mseScaleVec[2] = _mm256_setr_ps (a, a, b, b,
                                                     b, b, b, b);

                    uint8_t inc = 0;

                    for (uint8_t idx = 0; idx < 4; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }

                    
                        scaleOut8 = (int8_t*)scaleOut;

                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        inter1 = _mm256_packs_epi16(inter3, zeros);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        memcpy (scaleOut8, &inter1, 16);

                        scaleOut8 += 16;

                        scaleOut = (__m256i *)scaleOut8;
                    
                }
            }
        }
        // 64 QAM
        else if (modOrder == 8)
        {
            __m256 mseScaleVec[2];
            for (uint8_t symIdx = 0; symIdx < nsymb; symIdx++)
            {
                for (prbIdx = 0; (prbIdx + 2) <= nPRB; prbIdx+=2)
                {
                    float a = mseRealLow[0][prbIdx>>1], b = mseRealLow[1][prbIdx>>1];

                    mseScaleVec[0] = _mm256_set1_ps (a);

                    mseScaleVec[1] = _mm256_set1_ps (b);

                    for (uint8_t idx = 0; idx < 12; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[0]));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[1]));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[0]));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[1]));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }

                uint8_t excessPrbs = nPRB - prbIdx;

                if (excessPrbs)
                {
                    float a = mseRealLow[0][prbIdx>>1], b = mseRealLow[1][prbIdx>>1];

                    mseScaleVec[0] = _mm256_set1_ps (a);

                    mseScaleVec[1] = _mm256_set1_ps (b);

                    for (uint8_t idx = 0; idx < 6; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[0]));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[1]));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[0]));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[1]));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }
            }
        }
    }
    // Layers = 4
    else if (nLayers == 4)
    {
        // QPSK
        if (modOrder == 2)
        {
            for (uint8_t symIdx = 0; symIdx < nsymb; symIdx++)
            {
                for (prbIdx = 0; (prbIdx + 2) <= nPRB; prbIdx+=2)
                {
                    mseVec = _mm256_setr_ps (mseRealLow[0][prbIdx>>1], mseRealLow[0][prbIdx>>1],
                                             mseRealLow[1][prbIdx>>1], mseRealLow[1][prbIdx>>1],
                                             mseRealLow[2][prbIdx>>1], mseRealLow[2][prbIdx>>1],
                                             mseRealLow[3][prbIdx>>1], mseRealLow[3][prbIdx>>1]);

                    for (uint8_t idx = 0; idx < 6; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }

                uint8_t excessPrbs = nPRB - prbIdx;

                if (excessPrbs)
                {
                    mseVec = _mm256_setr_ps (mseRealLow[0][prbIdx>>1], mseRealLow[0][prbIdx>>1],
                                             mseRealLow[1][prbIdx>>1], mseRealLow[1][prbIdx>>1],
                                             mseRealLow[2][prbIdx>>1], mseRealLow[2][prbIdx>>1],
                                             mseRealLow[3][prbIdx>>1], mseRealLow[3][prbIdx>>1]);

                    for (uint8_t idx = 0; idx < 3; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseVec));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }
            }
        }

        // 16 QAM
        else if (modOrder == 4)
        {
            __m256 mseScaleVec[2];
            for (uint8_t symIdx = 0; symIdx < nsymb; symIdx++)
            {
                for (prbIdx = 0; (prbIdx + 2) <= nPRB; prbIdx+=2)
                {
                    float a = mseRealLow[0][prbIdx>>1], b = mseRealLow[1][prbIdx>>1], c = mseRealLow[2][prbIdx>>1], d = mseRealLow[3][prbIdx>>1];
                    mseScaleVec[0] = _mm256_setr_ps (a, a, a, a,
                                                     b, b, b, b);

                    mseScaleVec[1] = _mm256_setr_ps (c, c, c, c,
                                                     d, d, d, d);

                    for (uint8_t idx = 0; idx < 12; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[0]));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[1]));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[0]));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[1]));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }

                uint8_t excessPrbs = nPRB - prbIdx;

                if (excessPrbs)
                {
                    float a = mseRealLow[0][prbIdx>>1], b = mseRealLow[1][prbIdx>>1], c = mseRealLow[2][prbIdx>>1], d = mseRealLow[3][prbIdx>>1];
                    mseScaleVec[0] = _mm256_setr_ps (a, a, a, a,
                                                     b, b, b, b);

                    mseScaleVec[1] = _mm256_setr_ps (c, c, c, c,
                                                     d, d, d, d);

                    for (uint8_t idx = 0; idx < 6; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[0]));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[1]));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[0]));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[1]));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }
            }
        }


        // 64 QAM
        else if (modOrder == 6)
        {
            __m256 mseScaleVec[3];
            for (uint8_t symIdx = 0; symIdx < nsymb; symIdx++)
            {
                for (prbIdx = 0; (prbIdx + 2) <= nPRB; prbIdx+=2)
                {
                    float a = mseRealLow[0][prbIdx>>1], b = mseRealLow[1][prbIdx>>1], c = mseRealLow[2][prbIdx>>1], d = mseRealLow[3][prbIdx>>1];

                    mseScaleVec[0] = _mm256_setr_ps (a, a, a, a,
                                                     a, a, b, b);

                    mseScaleVec[1] = _mm256_setr_ps (b, b, b, b,
                                                     c, c, c, c);

                    mseScaleVec[2] = _mm256_setr_ps (c, c, d, d,
                                                     d, d, d, d);

                    uint8_t inc = 0;

                    for (uint8_t idx = 0; idx < 18; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }

                uint8_t excessPrbs = nPRB - prbIdx;

                if (excessPrbs)
                {
                    float a = mseRealLow[0][prbIdx>>1], b = mseRealLow[1][prbIdx>>1], c = mseRealLow[2][prbIdx>>1], d = mseRealLow[3][prbIdx>>1];

                    mseScaleVec[0] = _mm256_setr_ps (a, a, a, a,
                                                     a, a, b, b);

                    mseScaleVec[1] = _mm256_setr_ps (b, b, b, b,
                                                     c, c, c, c);

                    mseScaleVec[2] = _mm256_setr_ps (c, c, d, d,
                                                     d, d, d, d);

                    uint8_t inc = 0;

                    for (uint8_t idx = 0; idx < 9; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[inc%3]));
                        softOut++;
                        inc++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }
            }
        }

        // 256 QAM
        else if (modOrder == 8)
        {
            __m256 mseScaleVec[4]; 
            for (uint8_t symIdx = 0; symIdx < nsymb; symIdx++)
            {
                for (prbIdx = 0; (prbIdx + 2) <= nPRB; prbIdx+=2)
                {
                    mseScaleVec[0] = _mm256_set1_ps (mseRealLow[0][prbIdx>>1]);
                    mseScaleVec[1] = _mm256_set1_ps (mseRealLow[1][prbIdx>>1]);
                    mseScaleVec[2] = _mm256_set1_ps (mseRealLow[2][prbIdx>>1]);
                    mseScaleVec[3] = _mm256_set1_ps (mseRealLow[3][prbIdx>>1]);

                    for (uint8_t idx = 0; idx < 24; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[0]));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[1]));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[2]));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[3]));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }

                uint8_t excessPrbs = nPRB - prbIdx;

                if (excessPrbs)
                {
                    mseScaleVec[0] = _mm256_set1_ps (mseRealLow[0][prbIdx>>1]);
                    mseScaleVec[1] = _mm256_set1_ps (mseRealLow[1][prbIdx>>1]);
                    mseScaleVec[2] = _mm256_set1_ps (mseRealLow[2][prbIdx>>1]);
                    mseScaleVec[3] = _mm256_set1_ps (mseRealLow[3][prbIdx>>1]);

                    for (uint8_t idx = 0; idx < 12; idx++)
                    {
                        // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[0]));
                        softOut++;

                        // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[1]));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                        inter3 = _mm256_packs_epi32(inter1, inter2);

                        // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                        inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[2]));
                        softOut++;

                        // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                        inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*softOut, mseScaleVec[3]));
                        softOut++;

                        // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                        // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                        inter2 = _mm256_packs_epi32(inter1, inter2);
                        
                        // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                        // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                        // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                        inter1 = _mm256_packs_epi16(inter3, inter2);

                        // The above order is interleaved. So, we put them back in the correct format
                        // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                        // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                        inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                        // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                        // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                        inter1 = _mm256_max_epi8(inter1, minVal);
                        inter1 = _mm256_min_epi8(inter1, maxVal);

                        // Storing the Data in the output array
                        _mm256_storeu_si256 (scaleOut, inter1);

                        scaleOut++;
                    }
                }
            }
        }
    }
}


wnVoid mseSet(wnUInt8 segPart,
            wnUInt8 nsymb,
            wnUInt16 nPRB,
            wnUInt16 startToneIdx,
            wnUInt8 modOrder,
            wnUInt8 nLayers,
            wnUInt8 bwl,
            wnFlt mseRealLow[][12*MAX_NUM_PRB],
            wnFlt *mseReal,
            wnFlt *oneBymse_avg_re)
{

    __m256* mseVec = (__m256 *)mseReal;

    // Layers = 1
    if (nLayers == 1)
    {
        __m256 mseScaleVec = _mm256_set1_ps (oneBymse_avg_re[0]);
        wnFlt* mseValues = &mseRealLow[0][0];

        // QPSK
        if (modOrder == 2)
        {
            for (uint32_t scIdx = startToneIdx; scIdx < (startToneIdx+ 12*nPRB); scIdx+=4)
            {
                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseValues[scIdx], mseValues[scIdx],
                                                          mseValues[scIdx+1], mseValues[scIdx+1],
                                                          mseValues[scIdx+2], mseValues[scIdx+2],
                                                          mseValues[scIdx+3], mseValues[scIdx+3]),
                                          mseScaleVec);
                mseVec++;
            }
        }
        // 16 QAM
        else if (modOrder == 4)
        {
            for (uint32_t scIdx = startToneIdx; scIdx < (startToneIdx+ 12*nPRB); scIdx+=4)
            {
                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseValues[scIdx], mseValues[scIdx],
                                                          mseValues[scIdx], mseValues[scIdx],
                                                          mseValues[scIdx+1], mseValues[scIdx+1],
                                                          mseValues[scIdx+1], mseValues[scIdx+1]),
                                        mseScaleVec);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseValues[scIdx+2], mseValues[scIdx+2],
                                                          mseValues[scIdx+2], mseValues[scIdx+2],
                                                          mseValues[scIdx+3], mseValues[scIdx+3],
                                                          mseValues[scIdx+3], mseValues[scIdx+3]),
                                          mseScaleVec);
                mseVec++;
            }
        }
        // 64 QAM
        else if (modOrder == 6)
        {
            for (uint32_t scIdx = startToneIdx; scIdx < (startToneIdx+ 12*nPRB); scIdx+=4)
            {
                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseValues[scIdx], mseValues[scIdx],
                                                          mseValues[scIdx], mseValues[scIdx],
                                                          mseValues[scIdx], mseValues[scIdx],
                                                          mseValues[scIdx+1], mseValues[scIdx+1]),
                                          mseScaleVec);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseValues[scIdx+1], mseValues[scIdx+1],
                                                          mseValues[scIdx+1], mseValues[scIdx+1],
                                                          mseValues[scIdx+2], mseValues[scIdx+2],
                                                          mseValues[scIdx+2], mseValues[scIdx+2]),
                                          mseScaleVec);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseValues[scIdx+2], mseValues[scIdx+2],
                                                          mseValues[scIdx+3], mseValues[scIdx+3],
                                                          mseValues[scIdx+3], mseValues[scIdx+3],
                                                          mseValues[scIdx+3], mseValues[scIdx+3]),
                                          mseScaleVec);
                mseVec++;
            }
        }
        // 256 QAM
        else if (modOrder == 8)
        {
            for (uint32_t scIdx = startToneIdx; scIdx < (startToneIdx+ 12*nPRB); scIdx+=4)
            {
                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseValues[scIdx]), mseScaleVec);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseValues[scIdx+1]), mseScaleVec);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseValues[scIdx+2]), mseScaleVec);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseValues[scIdx+3]), mseScaleVec);
                mseVec++;
            }
        }
    }
    // Layers = 2
    else if (nLayers == 2)
    {
        // QPSK
        if (modOrder == 2)
        {
            __m256 mseScaleVec = _mm256_setr_ps (oneBymse_avg_re[0], oneBymse_avg_re[0],
                                                 oneBymse_avg_re[1], oneBymse_avg_re[1],
                                                 oneBymse_avg_re[0], oneBymse_avg_re[0],
                                                 oneBymse_avg_re[1], oneBymse_avg_re[1]);

            for (uint32_t scIdx = startToneIdx; scIdx < (startToneIdx+ 12*nPRB); scIdx+=4)
            {
                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[0][scIdx], mseRealLow[0][scIdx],
                                                          mseRealLow[1][scIdx], mseRealLow[1][scIdx],
                                                          mseRealLow[0][scIdx+1], mseRealLow[0][scIdx+1],
                                                          mseRealLow[1][scIdx+1], mseRealLow[1][scIdx+1]),
                                          mseScaleVec);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[0][scIdx+2], mseRealLow[0][scIdx+2],
                                                          mseRealLow[1][scIdx+2], mseRealLow[1][scIdx+2],
                                                          mseRealLow[0][scIdx+3], mseRealLow[0][scIdx+3],
                                                          mseRealLow[1][scIdx+3], mseRealLow[1][scIdx+3]),
                                          mseScaleVec);
                mseVec++;
            }
        }
        // 16 QAM
        else if (modOrder == 4)
        {
            __m256 mseScaleVec = _mm256_setr_ps (oneBymse_avg_re[0], oneBymse_avg_re[0],
                                                 oneBymse_avg_re[0], oneBymse_avg_re[0],
                                                 oneBymse_avg_re[1], oneBymse_avg_re[1],
                                                 oneBymse_avg_re[1], oneBymse_avg_re[1]);

            for (uint32_t scIdx = startToneIdx; scIdx < (startToneIdx + 12*nPRB); scIdx+=4)
            {
                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[0][scIdx], mseRealLow[0][scIdx],
                                                          mseRealLow[0][scIdx], mseRealLow[0][scIdx],
                                                          mseRealLow[1][scIdx], mseRealLow[1][scIdx],
                                                          mseRealLow[1][scIdx], mseRealLow[1][scIdx]),
                                          mseScaleVec);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[0][scIdx+1], mseRealLow[0][scIdx+1],
                                                          mseRealLow[0][scIdx+1], mseRealLow[0][scIdx+1],
                                                          mseRealLow[1][scIdx+1], mseRealLow[1][scIdx+1],
                                                          mseRealLow[1][scIdx+1], mseRealLow[1][scIdx+1]),
                                          mseScaleVec);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[0][scIdx+2], mseRealLow[0][scIdx+2],
                                                          mseRealLow[0][scIdx+2], mseRealLow[0][scIdx+2],
                                                          mseRealLow[1][scIdx+2], mseRealLow[1][scIdx+2],
                                                          mseRealLow[1][scIdx+2], mseRealLow[1][scIdx+2]),
                                          mseScaleVec);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[0][scIdx+3], mseRealLow[0][scIdx+3],
                                                          mseRealLow[0][scIdx+3], mseRealLow[0][scIdx+3],
                                                          mseRealLow[1][scIdx+3], mseRealLow[1][scIdx+3],
                                                          mseRealLow[1][scIdx+3], mseRealLow[1][scIdx+3]),
                                          mseScaleVec);
                mseVec++;
            }
        }
        // 64 QAM
        else if (modOrder == 6)
        {
            __m256 mseScaleVec1 = _mm256_setr_ps (oneBymse_avg_re[0], oneBymse_avg_re[0],
                                                  oneBymse_avg_re[0], oneBymse_avg_re[0],
                                                  oneBymse_avg_re[0], oneBymse_avg_re[0],
                                                  oneBymse_avg_re[1], oneBymse_avg_re[1]);

            __m256 mseScaleVec2 = _mm256_setr_ps (oneBymse_avg_re[1], oneBymse_avg_re[1],
                                                  oneBymse_avg_re[1], oneBymse_avg_re[1],
                                                  oneBymse_avg_re[0], oneBymse_avg_re[0],
                                                  oneBymse_avg_re[0], oneBymse_avg_re[0]);

            __m256 mseScaleVec3 = _mm256_setr_ps (oneBymse_avg_re[0], oneBymse_avg_re[0],
                                                  oneBymse_avg_re[1], oneBymse_avg_re[1],
                                                  oneBymse_avg_re[1], oneBymse_avg_re[1],
                                                  oneBymse_avg_re[1], oneBymse_avg_re[1]);

            for (uint32_t scIdx = startToneIdx; scIdx < (startToneIdx + 12*nPRB); scIdx+=2)
            {
                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[0][scIdx], mseRealLow[0][scIdx],
                                                          mseRealLow[0][scIdx], mseRealLow[0][scIdx],
                                                          mseRealLow[0][scIdx], mseRealLow[0][scIdx],
                                                          mseRealLow[1][scIdx], mseRealLow[1][scIdx]),
                                          mseScaleVec1);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[1][scIdx], mseRealLow[1][scIdx],
                                                          mseRealLow[1][scIdx], mseRealLow[1][scIdx],
                                                          mseRealLow[0][scIdx+1], mseRealLow[0][scIdx+1],
                                                          mseRealLow[0][scIdx+1], mseRealLow[0][scIdx+1]),
                                          mseScaleVec2);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[0][scIdx+1], mseRealLow[0][scIdx+1],
                                                          mseRealLow[1][scIdx+1], mseRealLow[1][scIdx+1],
                                                          mseRealLow[1][scIdx+1], mseRealLow[1][scIdx+1],
                                                          mseRealLow[1][scIdx+1], mseRealLow[1][scIdx+1]),
                                          mseScaleVec3);
                mseVec++;
            }
        }
        // 64 QAM
        else if (modOrder == 8)
        {
            __m256 mseScaleVec1 = _mm256_set1_ps (oneBymse_avg_re[0]);
            __m256 mseScaleVec2 = _mm256_set1_ps (oneBymse_avg_re[1]);

            for (uint32_t scIdx = startToneIdx; scIdx < (startToneIdx + 12*nPRB); scIdx+=4)
            {
                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseRealLow[0][scIdx]), mseScaleVec1);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseRealLow[1][scIdx]), mseScaleVec2);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseRealLow[0][scIdx+1]), mseScaleVec1);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseRealLow[1][scIdx+1]), mseScaleVec2);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseRealLow[0][scIdx+2]), mseScaleVec1);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseRealLow[1][scIdx+2]), mseScaleVec2);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseRealLow[0][scIdx+3]), mseScaleVec1);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseRealLow[1][scIdx+3]), mseScaleVec2);
                mseVec++;
            }
        }
    }
    // Layers = 4
    else if (nLayers == 4)
    {
        // QPSK
        if (modOrder == 2)
        {
            __m256 mseScaleVec = _mm256_setr_ps (oneBymse_avg_re[0], oneBymse_avg_re[0],
                                                 oneBymse_avg_re[1], oneBymse_avg_re[1],
                                                 oneBymse_avg_re[2], oneBymse_avg_re[2],
                                                 oneBymse_avg_re[3], oneBymse_avg_re[3]);

            for (uint32_t scIdx = startToneIdx; scIdx < (startToneIdx + 12*nPRB); scIdx+=4)
            {
                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[0][scIdx], mseRealLow[0][scIdx],
                                                          mseRealLow[1][scIdx], mseRealLow[1][scIdx],
                                                          mseRealLow[2][scIdx], mseRealLow[2][scIdx],
                                                          mseRealLow[3][scIdx], mseRealLow[3][scIdx]),
                                          mseScaleVec);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[0][scIdx+1], mseRealLow[0][scIdx+1],
                                                          mseRealLow[1][scIdx+1], mseRealLow[1][scIdx+1],
                                                          mseRealLow[2][scIdx+1], mseRealLow[2][scIdx+1],
                                                          mseRealLow[3][scIdx+1], mseRealLow[3][scIdx+1]),
                                          mseScaleVec);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[0][scIdx+2], mseRealLow[0][scIdx+2],
                                                          mseRealLow[1][scIdx+2], mseRealLow[1][scIdx+2],
                                                          mseRealLow[2][scIdx+2], mseRealLow[2][scIdx+2],
                                                          mseRealLow[3][scIdx+2], mseRealLow[3][scIdx+2]),
                                          mseScaleVec);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[0][scIdx+3], mseRealLow[0][scIdx+3],
                                                          mseRealLow[1][scIdx+3], mseRealLow[1][scIdx+3],
                                                          mseRealLow[2][scIdx+3], mseRealLow[2][scIdx+3],
                                                          mseRealLow[3][scIdx+3], mseRealLow[3][scIdx+3]),
                                          mseScaleVec);
                mseVec++;
            }
        }

        // 16 QAM
        else if (modOrder == 4)
        {
            __m256 mseScaleVec1 = _mm256_setr_ps (oneBymse_avg_re[0], oneBymse_avg_re[0],
                                                  oneBymse_avg_re[0], oneBymse_avg_re[0],
                                                  oneBymse_avg_re[1], oneBymse_avg_re[1],
                                                  oneBymse_avg_re[1], oneBymse_avg_re[1]);

            __m256 mseScaleVec2 = _mm256_setr_ps (oneBymse_avg_re[2], oneBymse_avg_re[2],
                                                  oneBymse_avg_re[2], oneBymse_avg_re[2],
                                                  oneBymse_avg_re[3], oneBymse_avg_re[3],
                                                  oneBymse_avg_re[3], oneBymse_avg_re[3]);

            for (uint32_t scIdx = startToneIdx; scIdx < (startToneIdx + 12*nPRB); scIdx+=2)
            {
                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[0][scIdx], mseRealLow[0][scIdx],
                                                          mseRealLow[0][scIdx], mseRealLow[0][scIdx],
                                                          mseRealLow[1][scIdx], mseRealLow[1][scIdx],
                                                          mseRealLow[1][scIdx], mseRealLow[1][scIdx]),
                                          mseScaleVec1);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[2][scIdx], mseRealLow[2][scIdx],
                                                          mseRealLow[2][scIdx], mseRealLow[2][scIdx],
                                                          mseRealLow[3][scIdx], mseRealLow[3][scIdx],
                                                          mseRealLow[3][scIdx], mseRealLow[3][scIdx]),
                                          mseScaleVec2);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[0][scIdx+1], mseRealLow[0][scIdx+1],
                                                          mseRealLow[0][scIdx+1], mseRealLow[0][scIdx+1],
                                                          mseRealLow[1][scIdx+1], mseRealLow[1][scIdx+1],
                                                          mseRealLow[1][scIdx+1], mseRealLow[1][scIdx+1]),
                                          mseScaleVec1);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[2][scIdx+1], mseRealLow[2][scIdx+1],
                                                          mseRealLow[2][scIdx+1], mseRealLow[2][scIdx+1],
                                                          mseRealLow[3][scIdx+1], mseRealLow[3][scIdx+1],
                                                          mseRealLow[3][scIdx+1], mseRealLow[3][scIdx+1]),
                                          mseScaleVec2);
                mseVec++;
            }
        }


        // 64 QAM
        else if (modOrder == 6)
        {
            __m256 mseScaleVec1 = _mm256_setr_ps (oneBymse_avg_re[0], oneBymse_avg_re[0],
                                                  oneBymse_avg_re[0], oneBymse_avg_re[0],
                                                  oneBymse_avg_re[0], oneBymse_avg_re[0],
                                                  oneBymse_avg_re[1], oneBymse_avg_re[1]);

            __m256 mseScaleVec2 = _mm256_setr_ps (oneBymse_avg_re[1], oneBymse_avg_re[1],
                                                  oneBymse_avg_re[1], oneBymse_avg_re[1],
                                                  oneBymse_avg_re[2], oneBymse_avg_re[2],
                                                  oneBymse_avg_re[2], oneBymse_avg_re[2]);


            __m256 mseScaleVec3 = _mm256_setr_ps (oneBymse_avg_re[2], oneBymse_avg_re[2],
                                                  oneBymse_avg_re[3], oneBymse_avg_re[3],
                                                  oneBymse_avg_re[3], oneBymse_avg_re[3],
                                                  oneBymse_avg_re[3], oneBymse_avg_re[3]);

            for (uint32_t scIdx = startToneIdx; scIdx < (startToneIdx + 12*nPRB); scIdx++)
            {
                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[0][scIdx], mseRealLow[0][scIdx],
                                                          mseRealLow[0][scIdx], mseRealLow[0][scIdx],
                                                          mseRealLow[0][scIdx], mseRealLow[0][scIdx],
                                                          mseRealLow[1][scIdx], mseRealLow[1][scIdx]),
                                          mseScaleVec1);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[1][scIdx], mseRealLow[1][scIdx],
                                                          mseRealLow[1][scIdx], mseRealLow[1][scIdx],
                                                          mseRealLow[2][scIdx], mseRealLow[2][scIdx],
                                                          mseRealLow[2][scIdx], mseRealLow[2][scIdx]),
                                          mseScaleVec2);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_setr_ps (mseRealLow[2][scIdx], mseRealLow[2][scIdx],
                                                          mseRealLow[3][scIdx], mseRealLow[3][scIdx],
                                                          mseRealLow[3][scIdx], mseRealLow[3][scIdx],
                                                          mseRealLow[3][scIdx], mseRealLow[3][scIdx]),
                                          mseScaleVec3);
                mseVec++;
            }
        }

        // 256 QAM
        else if (modOrder == 8)
        {
            __m256 mseScaleVec1 = _mm256_set1_ps (oneBymse_avg_re[0]);
            __m256 mseScaleVec2 = _mm256_set1_ps (oneBymse_avg_re[1]);
            __m256 mseScaleVec3 = _mm256_set1_ps (oneBymse_avg_re[2]);
            __m256 mseScaleVec4 = _mm256_set1_ps (oneBymse_avg_re[3]);

            for (uint32_t scIdx = startToneIdx; scIdx < (startToneIdx + 12*nPRB); scIdx+=2)
            {
                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseRealLow[0][scIdx]), mseScaleVec1);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseRealLow[1][scIdx]), mseScaleVec2);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseRealLow[2][scIdx]), mseScaleVec3);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseRealLow[3][scIdx]), mseScaleVec4);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseRealLow[0][scIdx+1]), mseScaleVec1);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseRealLow[1][scIdx+1]), mseScaleVec2);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseRealLow[2][scIdx+1]), mseScaleVec3);
                mseVec++;

                *mseVec = _mm256_mul_ps ( _mm256_set1_ps (mseRealLow[3][scIdx+1]), mseScaleVec4);
                mseVec++;
            }
        }
    }
}



void llrScalingQuant (wnUInt8 nsymb, wnUInt16 nPRB, wnUInt8 modOrder, wnUInt8 nLayers, wnUInt8 bwl, \
                      wnFlt *mseReal, wnFlt *inp, wnInt8 *out, puschConfigT* pNrPuschInParams)
{
    
    wnInt8 clip, *scaleOut8;
    wnUInt8 symIdx;
    wnUInt32 totalLlrs, llrsSet32, posIdx;
    
    // AVX Variables used
    __m256 *softOut, *mseVec;
    __m256i inter1, inter2, inter3, order, minVal, maxVal, zeros, *scaleOut;
    __m256 scaleVal;

    totalLlrs = nPRB*12*nLayers*modOrder;
    llrsSet32 = totalLlrs&0xFFFFFFE0;
    clip = (1<<(bwl-1))-1;

    softOut = (__m256*)inp;
    scaleOut = (__m256i*)out;
    order = _mm256_setr_epi32(0, 4, 1, 5, 2, 6, 3, 7);
    minVal = _mm256_set1_epi8(-1*clip);
    maxVal = _mm256_set1_epi8(clip);
    zeros = _mm256_setzero_si256();

    // LLR Scaling and Quantisation
    if(pNrPuschInParams->transformPrecode == 0)//If transform precoding is disabled.
    {
        for (symIdx = 0; symIdx < nsymb; symIdx++)
        {
            mseVec = (__m256*)mseReal;

            // All Tones
            for (posIdx = 0; posIdx < llrsSet32; posIdx+=32)
            {
                // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*mseVec, *softOut));
                softOut++;
                mseVec++;

                // Does Scaling and converts to fixed point for 2nd set of 32 LLRs
                inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*mseVec, *softOut));
                softOut++;
                mseVec++;

                // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                inter3 = _mm256_packs_epi32(inter1, inter2);

                // Does Scaling and converts to fixed point for 3rd set of 32 LLRs
                inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*mseVec, *softOut));
                softOut++;
                mseVec++;

                // Does Scaling and converts to fixed point for 4th set of 32 LLRs
                inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps (*mseVec, *softOut));
                softOut++;
                mseVec++;

                // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                inter2 = _mm256_packs_epi32(inter1, inter2);
                
                // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                inter1 = _mm256_packs_epi16(inter3, inter2);

                // The above order is interleaved. So, we put them back in the correct format
                // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                inter1 = _mm256_max_epi8(inter1, minVal);
                inter1 = _mm256_min_epi8(inter1, maxVal);

                // Storing the Data in the output array
                _mm256_storeu_si256 (scaleOut, inter1);

                scaleOut++;
            } // All Tones packed in 256 Words

            // Remaining Tones
            // Works on the edge cases 
            scaleOut8 = (wnInt8*)scaleOut;
            for (posIdx = llrsSet32; posIdx < totalLlrs; posIdx+=8)
            {
                // Does Scaling and converts to fixed point for 1st set of 32 LLRs
                inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps (*mseVec, *softOut));
                softOut++;
                mseVec++;

                inter2 = _mm256_packs_epi32(inter1, zeros);
                inter1 = _mm256_packs_epi16(inter2, zeros);

                inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                inter1 = _mm256_max_epi8(inter1, minVal);
                inter1 = _mm256_min_epi8(inter1, maxVal);

                memcpy (scaleOut8, &inter1, 8);

                scaleOut8+=8;
            } // Tones less than 256
            scaleOut = (__m256i *)scaleOut8;
        }   
    }
    else if(pNrPuschInParams->transformPrecode == 1)//If transform precoding is enabled.
    {
        scaleVal = _mm256_set1_ps(mseReal[0]);
        for (symIdx = 0; symIdx < nsymb; symIdx++)
        {
            // All Tones
            for (posIdx = 0; posIdx < llrsSet32; posIdx+=32)
            {
                // Converts to fixed point for 1st set of 32 LLRs
                inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps(*softOut, scaleVal));
                softOut++;

                // Converts to fixed point for 2nd set of 32 LLRs
                inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps(*softOut, scaleVal));
                softOut++;

                // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                // a1 a2 a3 a4 b1 b2 b3 b4 a5 a6 a7 a8 b5 b6 b7 b8
                inter3 = _mm256_packs_epi32(inter1, inter2);

                // Converts to fixed point for 3rd set of 32 LLRs
                inter1 = _mm256_cvtps_epi32 (_mm256_mul_ps(*softOut, scaleVal));
                softOut++;

                // Converts to fixed point for 4th set of 32 LLRs
                inter2 = _mm256_cvtps_epi32 (_mm256_mul_ps(*softOut, scaleVal));
                softOut++;

                // Converts 32 bit LLRs packed in above 2 sets into 16 bit LLRs
                // c1 c2 c3 c4 d1 d2 d3 d4 c5 c6 c7 c8 d5 d6 d7 d8
                inter2 = _mm256_packs_epi32(inter1, inter2);
                
                // Converts 16 bit LLRs packed in above 2 sets into 8 bit LLRs
                // (a1 a2 a3 a4) (b1 b2 b3 b4) (c1 c2 c3 c4) (d1 d2 d3 d4) ...
                // (a5 a6 a7 a8) (b5 b6 b7 b8) (c5 c6 c7 c8) (d5 d6 d7 d8)
                inter1 = _mm256_packs_epi16(inter3, inter2);

                // The above order is interleaved. So, we put them back in the correct format
                // (a1 a2 a3 a4) (a5 a6 a7 a8) (b1 b2 b3 b4) (b5 b6 b7 b8) ...
                // (c1 c2 c3 c4) (c5 c6 c7 c8) (d1 d2 d3 d4) (d5 d6 d7 d8)
                inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                inter1 = _mm256_max_epi8(inter1, minVal);
                inter1 = _mm256_min_epi8(inter1, maxVal);

                // Storing the Data in the output array
                _mm256_storeu_si256 (scaleOut, inter1);

                scaleOut++;
            } // All Tones packed in 256 Words

            // Remaining Tones
            // Works on the edge cases 
            scaleOut8 = (wnInt8*)scaleOut;
            for (posIdx = llrsSet32; posIdx < totalLlrs; posIdx+=8)
            {
                // Converts to fixed point for 1st set of 32 LLRs
                inter1 = _mm256_cvtps_epi32 (*softOut);
                softOut++;

                inter2 = _mm256_packs_epi32(inter1, zeros);
                inter1 = _mm256_packs_epi16(inter2, zeros);

                inter1 = _mm256_permutevar8x32_epi32 (inter1, order);

                // Clips them in the range we require [-B, B]: B = pow2(bwl)-1 where -1 is for techincal reasons
                // i.e., if we have 128, then -1*128 still gives 128 as we use 8 bit LLRs(we restrict this problem)
                inter1 = _mm256_max_epi8(inter1, minVal);
                inter1 = _mm256_min_epi8(inter1, maxVal);

                memcpy (scaleOut8, &inter1, 8);

                scaleOut8+=8;
            } // Tones less than 256
            scaleOut = (__m256i *)scaleOut8;
        }
    }
    
    //return 0;
}


void nr_pusch_llr_scaling_avx2(wnUInt8 segPart,
                               wnUInt8* startSym,
                               wnUInt8* endSym,
                               commonUlConfigT* pNrUlCommonParams,
                               puschConfigT* pNrPuschInParams,
                               P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams,
                               wnUInt8 bwl,\
                               wnFlt Lmax)
{
    #ifdef SNR24
        uint8_t scaleVal = (1<<(bwl-1));
        float scalingValTemp = (float)scaleVal/Lmax;
        float scalingVal;

        __m256 scalingVec;
        __m256 inter;
        uint16_t prbIdx = 0;
        wnUInt8 nLayers = pNrPuschInParams->nNrOfLayers;
        wnUInt8 nsymb = 6;
        wnInt32 nPRB = pNrPuschInParams->nRBSize;
        wnUInt16 startToneIdx = pNrPuschInParams->nBWPStart + pNrPuschInParams->nRBStart;
        wnInt8 modOrder = pNrPuschInParams->modulationOrder;


        for (uint8_t lyrIdx = 0; lyrIdx < nLayers; lyrIdx++)
        {
            scalingVal = scalingValTemp*pNrPuschOutParams->oneBymse_avg_re[lyrIdx];// /round((float)nPRB/2)

            scalingVec = _mm256_set1_ps (scalingVal);
            
            for (prbIdx = 0; (prbIdx+16) <= nPRB; prbIdx+=16)
            {
                inter = _mm256_loadu_ps (&pNrPuschOutParams->oneBymse_re[lyrIdx][prbIdx>>1]);
                inter = _mm256_mul_ps (inter, scalingVec);
                _mm256_storeu_ps (&pNrPuschOutParams->oneBymse_re[lyrIdx][prbIdx>>1], inter);
                
            }

            for (uint16_t nPrbIdx = prbIdx; nPrbIdx < nPRB; nPrbIdx+=2)
            {
                pNrPuschOutParams->oneBymse_re[lyrIdx][nPrbIdx>>1] *= scalingVal;
            }
        }

        wnUInt32 idxOffset = 0;
        if(symPart == LOW_PUSCH_SYM)
        {
            idxOffset = 0;
        }
        else
        {
            idxOffset = nPRB*12*6*modOrder*nLayers;
        }

        // Does MSE Scaling and places them in appropriate format for multiplication
        llrScalingQuant_24(symPart, nsymb, nPRB, startToneIdx, modOrder, nLayers, bwl, &pNrPuschOutParams->oneBymse_re[0][0],\
                        llrOut+idxOffset, llrFxd+idxOffset);
    #else
        __attribute__ ((aligned(32))) wnFlt mseReal[105600] = {0};//275*12*4*8
        wnUInt8 scaleVal, nLayers, nsymb, modOrder, symsComp, lyrIdx;
        wnUInt16 nPRB, startToneIdx;
        wnUInt32 idxOffset;
        wnFlt scalingVal;

        scaleVal = (1<<(bwl-1));
        scalingVal = (wnFlt)scaleVal/Lmax;
        nLayers = pNrPuschInParams->nNrOfLayers;
        nPRB = pNrPuschInParams->nRBSize;
        startToneIdx = pNrPuschInParams->nBWPStart + pNrPuschInParams->nRBStart;
        modOrder = pNrPuschInParams->modulationOrder;

        if(segPart == SEG1)
        {
            if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
            {
                nsymb = endSym[0] - startSym[0];
                symsComp = 0;
                idxOffset = 0;
            }
            else //if(pNrPuschInParams->nNrOfDMRSSymbols == 2)
            {
                nsymb = endSym[0] - startSym[0] - 1;
                symsComp = 0;
                idxOffset = 0;
            }
        }
        else if(segPart == SEG2)
        {
            if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
            {
                nsymb = endSym[1] - startSym[1];
                symsComp = endSym[0] - startSym[0];
                idxOffset = nPRB*12*modOrder*nLayers*symsComp;
            }
            else //if(pNrPuschInParams->nNrOfDMRSSymbols == 2)
            {
                nsymb = endSym[1] - startSym[1] - 1;
                symsComp = endSym[0] - startSym[0] - 1;
                idxOffset = nPRB*12*modOrder*nLayers*symsComp;
            }
        }
        else if(segPart == SEG3)
        {
            if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
            {
                nsymb = endSym[2] - startSym[2];
                symsComp = (endSym[0]-startSym[0])+(endSym[1]-startSym[1]);
                idxOffset = nPRB*12*modOrder*nLayers*symsComp;
            }
            else //if(pNrPuschInParams->nNrOfDMRSSymbols == 2)
            {
                nsymb = endSym[2] - startSym[2] - 1;
                symsComp = (endSym[0] - startSym[0] - 1)+(endSym[1] - startSym[1] - 1);
                idxOffset = nPRB*12*modOrder*nLayers*symsComp;
            }
            
        }
        else //if(segPart == SEG4)
        {
            if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
            {
                nsymb = endSym[3] - startSym[3];
                symsComp = (endSym[0] - startSym[0])+(endSym[1] - startSym[1])+(endSym[2] - startSym[2]);
                idxOffset = nPRB*12*modOrder*nLayers*symsComp;
            }
            else //if(pNrPuschInParams->nNrOfDMRSSymbols == 2)
            {
                nsymb = endSym[3] - startSym[3] - 1;
                symsComp = (endSym[0] - startSym[0] - 1)+(endSym[1] - startSym[1] - 1)+(endSym[2] - startSym[2] - 1);
                idxOffset = nPRB*12*modOrder*nLayers*symsComp;
            }
        }

        //If transform precoding is disabled.
        if(pNrPuschInParams->transformPrecode == 0)
        {
            // for (lyrIdx = 0; lyrIdx < nLayers; lyrIdx++)
            // pNrPuschOutParams->oneBymse_avg_re[lyrIdx] = scalingVal/pNrPuschOutParams->oneBymse_avg_re[lyrIdx];

            for (lyrIdx = 0; lyrIdx < nLayers; lyrIdx++)
                pNrPuschOutParams->oneBymse_avg_re[lyrIdx] = scalingVal*pNrPuschOutParams->oneBymse_avg_re[lyrIdx];

            // Does MSE Scaling and places them in appropriate format for multiplication
            mseSet(segPart, nsymb, nPRB, startToneIdx, modOrder, nLayers,\
                   bwl, &pNrPuschOutParams->oneBymse_re[0][0],\
                   mseReal, &pNrPuschOutParams->oneBymse_avg_re[0]);

            // LLR Scaling.
            llrScalingQuant(nsymb, nPRB, modOrder, nLayers, bwl,\
                            mseReal, llrOut+idxOffset, llrFxd+idxOffset, pNrPuschInParams);

        }
        else //if(pNrPuschInParams->transformPrecode == 1)
        {
            mseReal[0] = scalingVal;
            // LLR Quantisation.
            llrScalingQuant(nsymb, nPRB, modOrder, nLayers, bwl,\
                            mseReal, llrOut+idxOffset, llrFxd+idxOffset, pNrPuschInParams);
        }
        
        

    #endif



    #ifdef LLRScaleIRC
    int idxAnt, idxLyr, idxAlloc, idxSym, idxPrb, idxSc, idxRow, idxCol, itr, idx;
    int startPrbIdx, endPrbIdx, startScIdx, endScIdx;
    int ctPrbIdx, ctScIdx, avxPrbIdx, ctPrgIdx;

    int startSym, endSym, dmrsSym;
    int nRxAnt, nLyr, nSbAlloc, modOrder, nPRB;
    int nPrb_alloc, nPrb_total;

    unsigned short nRBSize = pNrPuschInParams->nRBSize;

    nRxAnt   = pNrPuschInParams->nAntennaPorts;
    nLyr     = pNrPuschInParams->nNrOfLayers;
    nSbAlloc = pNrPuschInParams->nSbAlloc;
    modOrder = pNrPuschInParams->modOrder;
    nPRB = pNrPuschInParams->nRBSize;

    float L1_llrScale, L2_llrScale, L3_llrScale, L4_llrScale;
    __m256 llr_vec, llrScale_v1, llrScale_v2, llrScale_v3, llrScale_v4;

    wnUInt32 idxOffset = 0;
    if(symPart == LOW_PUSCH_SYM)
    {
        idxOffset = 0;
    }
    else
    {
        idxOffset = nPRB*12*6*modOrder*nLyr;
    }

    nPrb_total = 0;
    for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
    {
        nPrb_alloc = pNrPuschInParams->nPrb[idxAlloc];

        nPrb_total += nPrb_alloc;

        startPrbIdx = pNrPuschInParams->rbStart[idxAlloc];
        endPrbIdx   = startPrbIdx + (pNrPuschInParams->nPrb[idxAlloc]) - 1;

        for(idxSym = 0; idxSym <= 6; idxSym++)
        {
            float *llrSymVecIn = &llrOut[idxOffset+nRBSize*12*modOrder*nLyr*idxSym];
            float *llrSymVecOut = &llrOutScaledFlt[idxOffset+nRBSize*12*modOrder*nLyr*idxSym];

            switch(modOrder)
            {
                case 2:
                    if(nLyr == 1)
                    {
                        for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb++)
                        {
                            llrScale_v1 = _mm256_set1_ps(pNrPuschOutParams->msePerPrb[0][idxPrb]);

                            for(idx = 0; idx < 3; idx++)
                            {
                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v1);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;
                            }
                        }
                    }
                    else if(nLyr == 2)
                    {
                        for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb++)
                        {
                            L1_llrScale = pNrPuschOutParams->msePerPrb[0][idxPrb];
                            L2_llrScale = pNrPuschOutParams->msePerPrb[1][idxPrb];

                            llrScale_v1 = _mm256_setr_ps( L1_llrScale, L1_llrScale,
                                                           L2_llrScale, L2_llrScale,
                                                           L1_llrScale, L1_llrScale,
                                                           L2_llrScale, L2_llrScale );

                            for(idx = 0; idx < 6; idx++)
                            {
                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v1);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;
                            }
                        }
                    }
                    else if(nLyr == 4)
                    {
                        for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb++)
                        {
                            L1_llrScale = pNrPuschOutParams->msePerPrb[0][idxPrb];
                            L2_llrScale = pNrPuschOutParams->msePerPrb[1][idxPrb];
                            L3_llrScale = pNrPuschOutParams->msePerPrb[2][idxPrb];
                            L4_llrScale = pNrPuschOutParams->msePerPrb[3][idxPrb];

                            llrScale_v1 = _mm256_setr_ps( L1_llrScale, L1_llrScale,
                                                           L2_llrScale, L2_llrScale,
                                                           L3_llrScale, L3_llrScale,
                                                           L4_llrScale, L4_llrScale );

                            for(idx = 0; idx < 12; idx++)
                            {
                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v1);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;
                            }
                        }
                    }

                break;

                case 4:
                    if(nLyr == 1)
                    {
                        for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb++)
                        {
                            llrScale_v1 = _mm256_set1_ps(pNrPuschOutParams->msePerPrb[0][idxPrb]);

                            for(idx = 0; idx < 6; idx++)
                            {
                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v1);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;
                            }
                        }
                    }
                    else if(nLyr == 2)
                    {
                        for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb++)
                        {
                            L1_llrScale = pNrPuschOutParams->msePerPrb[0][idxPrb];
                            L2_llrScale = pNrPuschOutParams->msePerPrb[1][idxPrb];

                            llrScale_v1 = _mm256_setr_ps( L1_llrScale, L1_llrScale,
                                                          L1_llrScale, L1_llrScale,
                                                          L2_llrScale, L2_llrScale,
                                                          L2_llrScale, L2_llrScale );

                            for(idx = 0; idx < 12; idx++)
                            {
                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v1);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;
                            }
                        }
                    }
                    else if(nLyr == 4)
                    {
                        for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb++)
                        {
                            L1_llrScale = pNrPuschOutParams->msePerPrb[0][idxPrb];
                            L2_llrScale = pNrPuschOutParams->msePerPrb[1][idxPrb];
                            L3_llrScale = pNrPuschOutParams->msePerPrb[2][idxPrb];
                            L4_llrScale = pNrPuschOutParams->msePerPrb[3][idxPrb];

                            llrScale_v1 = _mm256_setr_ps( L1_llrScale, L1_llrScale,
                                                          L1_llrScale, L1_llrScale,
                                                          L2_llrScale, L2_llrScale,
                                                          L2_llrScale, L2_llrScale );

                            llrScale_v2 = _mm256_setr_ps( L3_llrScale, L3_llrScale,
                                                          L3_llrScale, L3_llrScale,
                                                          L4_llrScale, L4_llrScale,
                                                          L4_llrScale, L4_llrScale );

                            for(idx = 0; idx < 12; idx++)
                            {
                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v1);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;

                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v2);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;
                            }
                        }
                    }

                break;

                case 6:
                    if(nLyr == 1)
                    {
                        for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb++)
                        {
                            llrScale_v1 = _mm256_set1_ps(pNrPuschOutParams->msePerPrb[0][idxPrb]);

                            for(idx = 0; idx < 9; idx++)
                            {
                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v1);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;
                            }
                        }
                    }
                    else if(nLyr == 2)
                    {
                        for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb++)
                        {
                            L1_llrScale = pNrPuschOutParams->msePerPrb[0][idxPrb];
                            L2_llrScale = pNrPuschOutParams->msePerPrb[1][idxPrb];

                            llrScale_v1 = _mm256_setr_ps( L1_llrScale, L1_llrScale,
                                                          L1_llrScale, L1_llrScale,
                                                          L1_llrScale, L1_llrScale,
                                                          L2_llrScale, L2_llrScale );

                            llrScale_v2 = _mm256_setr_ps( L2_llrScale, L2_llrScale,
                                                          L2_llrScale, L2_llrScale,
                                                          L1_llrScale, L1_llrScale,
                                                          L1_llrScale, L1_llrScale );

                            llrScale_v3 = _mm256_setr_ps( L1_llrScale, L1_llrScale,
                                                          L2_llrScale, L2_llrScale,
                                                          L2_llrScale, L2_llrScale,
                                                          L2_llrScale, L2_llrScale );

                            for(idx = 0; idx < 6; idx++)
                            {
                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v1);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;

                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v2);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;

                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v3);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;
                            }
                        }
                    }
                    else if(nLyr == 4)
                    {
                        for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb++)
                        {
                            L1_llrScale = pNrPuschOutParams->msePerPrb[0][idxPrb];
                            L2_llrScale = pNrPuschOutParams->msePerPrb[1][idxPrb];
                            L3_llrScale = pNrPuschOutParams->msePerPrb[2][idxPrb];
                            L4_llrScale = pNrPuschOutParams->msePerPrb[3][idxPrb];

                            llrScale_v1 = _mm256_setr_ps( L1_llrScale, L1_llrScale,
                                                          L1_llrScale, L1_llrScale,
                                                          L1_llrScale, L1_llrScale,
                                                          L2_llrScale, L2_llrScale );

                            llrScale_v2 = _mm256_setr_ps( L2_llrScale, L2_llrScale,
                                                          L2_llrScale, L2_llrScale,
                                                          L3_llrScale, L3_llrScale,
                                                          L3_llrScale, L3_llrScale );

                            llrScale_v3 = _mm256_setr_ps( L3_llrScale, L3_llrScale,
                                                          L4_llrScale, L4_llrScale,
                                                          L4_llrScale, L4_llrScale,
                                                          L4_llrScale, L4_llrScale );

                            for(idx = 0; idx < 12; idx++)
                            {
                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v1);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;

                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v2);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;

                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v3);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;
                            }
                        }
                    }

                break;

                case 8:
                    if(nLyr == 1)
                    {
                        for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb++)
                        {
                            llrScale_v1 = _mm256_set1_ps(pNrPuschOutParams->msePerPrb[0][idxPrb]);

                            for(idx = 0; idx < 12; idx++)
                            {
                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v1);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;
                            }
                        }
                    }
                    else if(nLyr == 2)
                    {
                        for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb++)
                        {
                            L1_llrScale = pNrPuschOutParams->msePerPrb[0][idxPrb];
                            L2_llrScale = pNrPuschOutParams->msePerPrb[1][idxPrb];

                            llrScale_v1 = _mm256_setr_ps( L1_llrScale, L1_llrScale,
                                                          L1_llrScale, L1_llrScale,
                                                          L1_llrScale, L1_llrScale,
                                                          L1_llrScale, L1_llrScale );

                            llrScale_v2 = _mm256_setr_ps( L2_llrScale, L2_llrScale,
                                                          L2_llrScale, L2_llrScale,
                                                          L2_llrScale, L2_llrScale,
                                                          L2_llrScale, L2_llrScale );

                            for(idx = 0; idx < 12; idx++)
                            {
                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v1);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;

                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v2);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;
                            }
                        }
                    }
                    else if(nLyr == 4)
                    {
                        for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb++)
                        {
                            L1_llrScale = pNrPuschOutParams->msePerPrb[0][idxPrb];
                            L2_llrScale = pNrPuschOutParams->msePerPrb[1][idxPrb];
                            L3_llrScale = pNrPuschOutParams->msePerPrb[2][idxPrb];
                            L4_llrScale = pNrPuschOutParams->msePerPrb[3][idxPrb];

                            llrScale_v1 = _mm256_setr_ps( L1_llrScale, L1_llrScale,
                                                          L1_llrScale, L1_llrScale,
                                                          L1_llrScale, L1_llrScale,
                                                          L1_llrScale, L1_llrScale );

                            llrScale_v2 = _mm256_setr_ps( L2_llrScale, L2_llrScale,
                                                          L2_llrScale, L2_llrScale,
                                                          L2_llrScale, L2_llrScale,
                                                          L2_llrScale, L2_llrScale );

                            llrScale_v3 = _mm256_setr_ps( L3_llrScale, L3_llrScale,
                                                          L3_llrScale, L3_llrScale,
                                                          L3_llrScale, L3_llrScale,
                                                          L3_llrScale, L3_llrScale );

                            llrScale_v4 = _mm256_setr_ps( L4_llrScale, L4_llrScale,
                                                          L4_llrScale, L4_llrScale,
                                                          L4_llrScale, L4_llrScale,
                                                          L4_llrScale, L4_llrScale );

                            for(idx = 0; idx < 12; idx++)
                            {
                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v1);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;

                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v2);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;

                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v3);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;

                                llr_vec = _mm256_loadu_ps(llrSymVecIn);

                                llr_vec = _mm256_mul_ps(llr_vec, llrScale_v4);

                                _mm256_storeu_ps(llrSymVecOut, llr_vec);

                                llrSymVecIn += 8;
                                llrSymVecOut += 8;
                            }
                        }
                    }

                break;
            }
        }
    }


    wnUInt8 scaleVal = (1<<(bwl-1));
    __m256 scaleValAVX = _mm256_set1_ps(scaleVal/(float)Lmax);

    __m256i minVal = _mm256_set1_epi8(-(scaleVal - 1));
    __m256i maxVal = _mm256_set1_epi8((scaleVal - 1));

    unsigned int  floorUpperLimit = (nPRB*12*modOrder*6*nLyr) / (8*4);

    __m256 *llrFloatAVX = (__m256*) &llrOutScaledFlt[idxOffset];
    __m256i *llrFxdAVX = (__m256*)&llrFxd[idxOffset];
    __m256i tempFxd0, tempFxd1, tempFxd2, tempFxd3, tempFxd4, temp1, temp2;

   for(wnUInt32 i = 0;i< floorUpperLimit ;i++)
    {
        tempFxd0 = _mm256_cvtps_epi32(_mm256_mul_ps(*(llrFloatAVX),  scaleValAVX));
        tempFxd1 = _mm256_cvtps_epi32(_mm256_mul_ps(*(llrFloatAVX+1),  scaleValAVX));
        tempFxd2 = _mm256_packs_epi32(tempFxd0, tempFxd1);
        tempFxd0 = _mm256_cvtps_epi32(_mm256_mul_ps(*(llrFloatAVX+2),  scaleValAVX));
        tempFxd1 = _mm256_cvtps_epi32(_mm256_mul_ps(*(llrFloatAVX+3),  scaleValAVX));

        tempFxd3 = _mm256_packs_epi32(tempFxd0, tempFxd1);
        tempFxd4 = _mm256_packs_epi16(tempFxd2, tempFxd3);

        temp1 = _mm256_permutevar8x32_epi32(tempFxd4, _mm256_setr_epi32(0,4, 1,5, 2,6, 3,7));

        temp2 = _mm256_max_epi8(temp1, minVal);
        *llrFxdAVX = _mm256_min_epi8(temp2, maxVal);
        // temp2 = _mm256_min_epi8(temp2, maxVal);
        // _mm256_storeu_si256(llrFxdAVX, temp2);

        llrFloatAVX += 4;
        llrFxdAVX += 1;
    }


    float tempFlt;
    short tempFxd;
    for(wnUInt32 i = floorUpperLimit*8*4;i<(nPRB*12*modOrder*6);i++)//Convert remainder llrs. (max loop iteration = 31)
    {
        tempFlt = llrOutScaledFlt[idxOffset+i]*scaleVal;
        tempFlt = tempFlt>0?(tempFlt+0.5):(tempFlt-0.5);//Rounding
        tempFxd = tempFlt;

        if(tempFxd > (scaleVal - 1))
        tempFxd = (scaleVal - 1);
        else if(tempFxd < (scaleVal - 1))
        tempFxd = -(scaleVal - 1);

        llrFxd[idxOffset+i] = tempFxd;
    }

    #endif
}//*/

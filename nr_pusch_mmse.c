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
  *         for NR PUSCH Channel Estimation and Equalization using MMSE Algorithm.
  * @file    : nr_pusch_mmse.c
  * @ingroup : nr_pusch
  * @author  : MIRZA SAMI BAIG
  **/

#include "nr_pusch_mmse.h"
// #define SNR24
#include "wnNrExeTimeX86Intel.h"
#include <fftw3.h>

#define MAX_RX_ANT MAX_SUMIMO_PUSCH_ANTENNA_PORTS 
void mat_inv(int n, float complex mat[][MAX_RX_ANT], float complex inv[][MAX_RX_ANT])
{
    int i,j,k;
    float complex ratio,a;

    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            if(i==j)
            {
                 inv[i][j] = 1.0;
            }
            else
            {
                inv[i][j] = 0.0;
            }
        }
    }

    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            if(i!=j)
            {
                ratio = mat[j][i] / mat[i][i];
                for(k=0; k<n; k++)
                {
                    mat[j][k] = mat[j][k] - (ratio * mat[i][k]);
                    inv[j][k] = inv[j][k] - (ratio * inv[i][k]);
                }
            }
        }
    }

    for(i=0; i<n; i++)
    {
        a = mat[i][i];
        for(j=0; j<n; j++)
        {
            mat[i][j] = mat[i][j] / a;
            inv[i][j] = inv[i][j] / a;
        }
    }

    return;
}

void mat_mul(int m1, int n1, float complex mat1[][MAX_RX_ANT], int m2, int n2, float complex mat2[][MAX_RX_ANT], float complex res[][MAX_RX_ANT])
{
    int r,c,i;

    for(r=0; r<m1; r++)
    {
        for(c=0; c<n2; c++)
        {
            res[r][c] = 0;
            for(i=0; i<n1; i++)
            {
                res[r][c] += (mat1[r][i] * mat2[i][c]);
            }
        }
    }
}



// { perPRB Start
/*void nr_pusch_mmse_perPrb_est_avx2( int symPart,
                                    commonUlConfigT* pNrUlCommonParams,
                                    puschConfigT* pNrPuschInParams,
                                    P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams)
{
    int idxAnt, idxLyr, idxAlloc, idxPrb, idxSc, idxRow, idxCol, itr, idx;
    int startPrbIdx, endPrbIdx, startScIdx, endScIdx;
    int ctPrbIdx, ctScIdx, avxPrbIdx, ctPrgIdx;
    int printPrbIdx = 0, avxPrintIdx = 0;

    int prgScStartIdx, prbScStartIdx[8];

    int dmrsSym, nRxAnt, nLyr, nSbAlloc, nPrb_alloc, nPrb_total, remPrb, remSc;
    int nEvenTonesLyr=0, nOddTonesLyr=0;

    nRxAnt   = pNrPuschInParams->nAntennaPorts;
    nLyr     = pNrPuschInParams->nNrOfLayers;
    nSbAlloc = pNrPuschInParams->nSbAlloc;

    __m256 genDmrs_re[6], genDmrs_im[6];
    __m256 genDmrs2_re[6], genDmrs2_im[6];
    __m256 recDmrs_evenRe[nRxAnt][6], recDmrs_evenIm[nRxAnt][6];
    __m256 recDmrs_oddRe[nRxAnt][6], recDmrs_oddIm[nRxAnt][6];
    __m256 lsEst_evenRe[nRxAnt][6], lsEst_evenIm[nRxAnt][6];
    __m256 lsEst_oddRe[nRxAnt][6], lsEst_oddIm[nRxAnt][6];

    __m256 bd_re, ad_im;
    __m256 mid_re, mid_im;
    __m256 mid2_re, mid2_im;
    __m256 mid3_re, mid3_im;
    __m256 out_re, out_im;

    __m256 scale_lsEst = _mm256_set1_ps(0.5012);
    __m256 div5Vec = _mm256_set1_ps(0.2);
    __m256 div4Vec = _mm256_set1_ps(0.25);
    __m256 div6Vec = _mm256_set1_ps(0.16667);

    float complex toe_vec[6];
    __m256 toe_reVec[6], toe_imVec[6];
    __m256 toe_perAnt_re, toe_perAnt_im;
    __m256 toe_perPrb_re, toe_perPrb_im;
    __m256 toe_perPrb_abs;

    complex toe_evenTones, toe_oddTones, toe_overAll;

    __m256 h1Avg_re[nRxAnt], h1Avg_im[nRxAnt];
    __m256 h2Avg_re[nRxAnt], h2Avg_im[nRxAnt];

    __m256 Inf_re[nRxAnt], Inf_im[nRxAnt];

    __m256 Rnn_re[nRxAnt][nRxAnt], Rnn_im[nRxAnt][nRxAnt];
    __m256 RnnTemp_re[nRxAnt][nRxAnt], RnnTemp_im[nRxAnt][nRxAnt];

    complex h1Avg[nRxAnt][8], h2Avg[nRxAnt][8], hAvg;
    complex Rnn[nRxAnt][nRxAnt], Inf[nRxAnt];

    for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
    {
        if((pNrPuschInParams->nPortIndex[idxLyr] == 0) || (pNrPuschInParams->nPortIndex[idxLyr] == 1))
        {
            nEvenTonesLyr += 1;
        }
        else if((pNrPuschInParams->nPortIndex[idxLyr] == 2) || (pNrPuschInParams->nPortIndex[idxLyr] == 3))
        {
            nOddTonesLyr += 1;
        }
    }
    pNrPuschOutParams->nEvenTonesLyr = nEvenTonesLyr;
    pNrPuschOutParams->nOddTonesLyr = nOddTonesLyr;

    if(symPart == LOW_PUSCH_SYM)
    {
        dmrsSym = pNrPuschInParams->dmrs1_sym;
    }
    else
    {
        dmrsSym = pNrPuschInParams->dmrs2_sym;
    }

    for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
    {
        nPrb_alloc = pNrPuschInParams->nPrb[idxAlloc];

        startPrbIdx = pNrPuschInParams->rbStart[idxAlloc];
        endPrbIdx = startPrbIdx + nPrb_alloc - 1;

        startScIdx = startPrbIdx * 6;
        endScIdx   = endPrbIdx * 6;

        ctPrgIdx = 0;
        for(idxPrb = startPrbIdx; idxPrb <= endPrbIdx; idxPrb+=8, ctPrgIdx++)
        {
            prgScStartIdx = idxPrb * 6;
            ctScIdx = prgScStartIdx;
            for(idx = 0; idx < 8; idx++)
            {
                prbScStartIdx[idx] = ctScIdx;

                ctScIdx += 6;
            }

            //Loading Generated DMRS in local buffer
            float *genDmrs_rePtr = &pNrPuschOutParams->genDMRS[startScIdx];
            float *genDmrs_imPtr = &pNrPuschOutParams->genDMRS[HALF_MAX_NUM_SC + startScIdx];

            float *genDmrs2_rePtr = &pNrPuschOutParams->genDMRS2[startScIdx];
            float *genDmrs2_imPtr = &pNrPuschOutParams->genDMRS2[HALF_MAX_NUM_SC + startScIdx];

            for(idxSc = 0; idxSc < 6; idxSc++)
            {
                genDmrs_re[idxSc] = _mm256_set_ps( genDmrs_rePtr[prbScStartIdx[7] + idxSc],
                                                   genDmrs_rePtr[prbScStartIdx[6] + idxSc],
                                                   genDmrs_rePtr[prbScStartIdx[5] + idxSc],
                                                   genDmrs_rePtr[prbScStartIdx[4] + idxSc],
                                                   genDmrs_rePtr[prbScStartIdx[3] + idxSc],
                                                   genDmrs_rePtr[prbScStartIdx[2] + idxSc],
                                                   genDmrs_rePtr[prbScStartIdx[1] + idxSc],
                                                   genDmrs_rePtr[prbScStartIdx[0] + idxSc] );

                genDmrs_im[idxSc] = _mm256_set_ps( genDmrs_imPtr[prbScStartIdx[7] + idxSc],
                                                   genDmrs_imPtr[prbScStartIdx[6] + idxSc],
                                                   genDmrs_imPtr[prbScStartIdx[5] + idxSc],
                                                   genDmrs_imPtr[prbScStartIdx[4] + idxSc],
                                                   genDmrs_imPtr[prbScStartIdx[3] + idxSc],
                                                   genDmrs_imPtr[prbScStartIdx[2] + idxSc],
                                                   genDmrs_imPtr[prbScStartIdx[1] + idxSc],
                                                   genDmrs_imPtr[prbScStartIdx[0] + idxSc] );

                if(nEvenTonesLyr == 2)
                {
                    genDmrs2_re[idxSc] = _mm256_set_ps( genDmrs2_rePtr[prbScStartIdx[7] + idxSc],
                                                        genDmrs2_rePtr[prbScStartIdx[6] + idxSc],
                                                        genDmrs2_rePtr[prbScStartIdx[5] + idxSc],
                                                        genDmrs2_rePtr[prbScStartIdx[4] + idxSc],
                                                        genDmrs2_rePtr[prbScStartIdx[3] + idxSc],
                                                        genDmrs2_rePtr[prbScStartIdx[2] + idxSc],
                                                        genDmrs2_rePtr[prbScStartIdx[1] + idxSc],
                                                        genDmrs2_rePtr[prbScStartIdx[0] + idxSc] );

                    genDmrs2_im[idxSc] = _mm256_set_ps( genDmrs2_imPtr[prbScStartIdx[7] + idxSc],
                                                        genDmrs2_imPtr[prbScStartIdx[6] + idxSc],
                                                        genDmrs2_imPtr[prbScStartIdx[5] + idxSc],
                                                        genDmrs2_imPtr[prbScStartIdx[4] + idxSc],
                                                        genDmrs2_imPtr[prbScStartIdx[3] + idxSc],
                                                        genDmrs2_imPtr[prbScStartIdx[2] + idxSc],
                                                        genDmrs2_imPtr[prbScStartIdx[1] + idxSc],
                                                        genDmrs2_imPtr[prbScStartIdx[0] + idxSc] );
                }
            }

            if(nEvenTonesLyr > 0)
            {
                //Loading Received DMRS in local buffer and Performing LS Estimation
                for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                {
                    float *recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                    float *recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC/2) + startScIdx];

                    for(idxSc = 0; idxSc < 6; idxSc++)
                    {
                        recDmrs_evenRe[idxAnt][idxSc] = _mm256_set_ps( recDmrs_rePtr[prbScStartIdx[7] + idxSc],
                                                                       recDmrs_rePtr[prbScStartIdx[6] + idxSc],
                                                                       recDmrs_rePtr[prbScStartIdx[5] + idxSc],
                                                                       recDmrs_rePtr[prbScStartIdx[4] + idxSc],
                                                                       recDmrs_rePtr[prbScStartIdx[3] + idxSc],
                                                                       recDmrs_rePtr[prbScStartIdx[2] + idxSc],
                                                                       recDmrs_rePtr[prbScStartIdx[1] + idxSc],
                                                                       recDmrs_rePtr[prbScStartIdx[0] + idxSc] );

                        recDmrs_evenIm[idxAnt][idxSc] = _mm256_set_ps( recDmrs_imPtr[prbScStartIdx[7] + idxSc],
                                                                       recDmrs_imPtr[prbScStartIdx[6] + idxSc],
                                                                       recDmrs_imPtr[prbScStartIdx[5] + idxSc],
                                                                       recDmrs_imPtr[prbScStartIdx[4] + idxSc],
                                                                       recDmrs_imPtr[prbScStartIdx[3] + idxSc],
                                                                       recDmrs_imPtr[prbScStartIdx[2] + idxSc],
                                                                       recDmrs_imPtr[prbScStartIdx[1] + idxSc],
                                                                       recDmrs_imPtr[prbScStartIdx[0] + idxSc] );

                        bd_re = _mm256_mul_ps(recDmrs_evenIm[idxAnt][idxSc], genDmrs_im[idxSc]);
                        ad_im = _mm256_mul_ps(recDmrs_evenRe[idxAnt][idxSc], genDmrs_im[idxSc]);

                        out_re = _mm256_fmadd_ps(recDmrs_evenRe[idxAnt][idxSc], genDmrs_re[idxSc], bd_re);
                        out_im = _mm256_fmsub_ps(recDmrs_evenIm[idxAnt][idxSc], genDmrs_re[idxSc], ad_im);

                        lsEst_evenRe[idxAnt][idxSc] = _mm256_mul_ps(out_re, scale_lsEst);
                        lsEst_evenIm[idxAnt][idxSc] = _mm256_mul_ps(out_im, scale_lsEst);
                    }
                }

                //Time offset estimation (TOE)
                toe_evenTones = 0;

                toe_perAnt_re = _mm256_setzero_ps();
                toe_perAnt_im = _mm256_setzero_ps();

                for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                {
                    toe_perPrb_re = _mm256_setzero_ps();
                    toe_perPrb_im = _mm256_setzero_ps();

                    if(nEvenTonesLyr == 1)
                    {
                        for(idxSc = 0; idxSc < 5; idxSc++)
                        {
                            bd_re = _mm256_mul_ps(lsEst_evenIm[idxAnt][idxSc+1], lsEst_evenIm[idxAnt][idxSc]);
                            ad_im = _mm256_mul_ps(lsEst_evenRe[idxAnt][idxSc+1], lsEst_evenIm[idxAnt][idxSc]);

                            out_re = _mm256_fmadd_ps(lsEst_evenRe[idxAnt][idxSc+1], lsEst_evenRe[idxAnt][idxSc], bd_re);
                            out_im = _mm256_fmsub_ps(lsEst_evenIm[idxAnt][idxSc+1], lsEst_evenRe[idxAnt][idxSc], ad_im);

                            toe_perPrb_re = _mm256_add_ps(toe_perPrb_re, out_re);
                            toe_perPrb_im = _mm256_add_ps(toe_perPrb_im, out_im);
                        }
                        toe_perPrb_re = _mm256_mul_ps(toe_perPrb_re, div5Vec);
                        toe_perPrb_im = _mm256_mul_ps(toe_perPrb_im, div5Vec);
                    }
                    else if(nEvenTonesLyr == 2)
                    {
                        for(idxSc = 0; idxSc < 4; idxSc++)
                        {
                            bd_re = _mm256_mul_ps(lsEst_evenIm[idxAnt][idxSc+1], lsEst_evenIm[idxAnt][idxSc]);
                            ad_im = _mm256_mul_ps(lsEst_evenRe[idxAnt][idxSc+1], lsEst_evenIm[idxAnt][idxSc]);

                            out_re = _mm256_fmadd_ps(lsEst_evenRe[idxAnt][idxSc+1], lsEst_evenRe[idxAnt][idxSc], bd_re);
                            out_im = _mm256_fmsub_ps(lsEst_evenIm[idxAnt][idxSc+1], lsEst_evenRe[idxAnt][idxSc], ad_im);

                            toe_perPrb_re = _mm256_add_ps(toe_perPrb_re, out_re);
                            toe_perPrb_im = _mm256_add_ps(toe_perPrb_im, out_im);
                        }
                        toe_perPrb_re = _mm256_mul_ps(toe_perPrb_re, div4Vec);
                        toe_perPrb_im = _mm256_mul_ps(toe_perPrb_im, div4Vec);
                    }

                    toe_perAnt_re = _mm256_add_ps(toe_perAnt_re, toe_perPrb_re);
                    toe_perAnt_im = _mm256_add_ps(toe_perAnt_im, toe_perPrb_im);
                }

                mid_re = _mm256_permute2f128_ps (toe_perAnt_re, toe_perAnt_re, 1);
                toe_perAnt_re = _mm256_add_ps(toe_perAnt_re, mid_re);
                toe_perAnt_re = _mm256_hadd_ps(toe_perAnt_re, toe_perAnt_re);
                toe_perAnt_re = _mm256_hadd_ps(toe_perAnt_re, toe_perAnt_re);

                mid_im = _mm256_permute2f128_ps (toe_perAnt_im, toe_perAnt_im, 1);
                toe_perAnt_im = _mm256_add_ps(toe_perAnt_im, mid_im);
                toe_perAnt_im = _mm256_hadd_ps(toe_perAnt_im, toe_perAnt_im);
                toe_perAnt_im = _mm256_hadd_ps(toe_perAnt_im, toe_perAnt_im);

                toe_evenTones = (toe_perAnt_re[0] + I *  toe_perAnt_im[0]) / (float)(nRxAnt*8);

                if(nEvenTonesLyr == 2)
                {
                    toe_evenTones = csqrt(toe_evenTones);
                }
            }

            if(nOddTonesLyr > 0)
            {
                //Loading Received DMRS in local buffer and Performing LS Estimation
                for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                {
                    float *recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC + startScIdx];
                    float *recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][ONE_HALF_MAX_NUM_SC + startScIdx];

                    for(idxSc = 0; idxSc < 6; idxSc++)
                    {
                        recDmrs_oddRe[idxAnt][idxSc] = _mm256_set_ps( recDmrs_rePtr[prbScStartIdx[7] + idxSc],
                                                                      recDmrs_rePtr[prbScStartIdx[6] + idxSc],
                                                                      recDmrs_rePtr[prbScStartIdx[5] + idxSc],
                                                                      recDmrs_rePtr[prbScStartIdx[4] + idxSc],
                                                                      recDmrs_rePtr[prbScStartIdx[3] + idxSc],
                                                                      recDmrs_rePtr[prbScStartIdx[2] + idxSc],
                                                                      recDmrs_rePtr[prbScStartIdx[1] + idxSc],
                                                                      recDmrs_rePtr[prbScStartIdx[0] + idxSc] );

                        recDmrs_oddIm[idxAnt][idxSc] = _mm256_set_ps( recDmrs_imPtr[prbScStartIdx[7] + idxSc],
                                                                      recDmrs_imPtr[prbScStartIdx[6] + idxSc],
                                                                      recDmrs_imPtr[prbScStartIdx[5] + idxSc],
                                                                      recDmrs_imPtr[prbScStartIdx[4] + idxSc],
                                                                      recDmrs_imPtr[prbScStartIdx[3] + idxSc],
                                                                      recDmrs_imPtr[prbScStartIdx[2] + idxSc],
                                                                      recDmrs_imPtr[prbScStartIdx[1] + idxSc],
                                                                      recDmrs_imPtr[prbScStartIdx[0] + idxSc] );

                        bd_re = _mm256_mul_ps(recDmrs_oddIm[idxAnt][idxSc], genDmrs_im[idxSc]);
                        ad_im = _mm256_mul_ps(recDmrs_oddRe[idxAnt][idxSc], genDmrs_im[idxSc]);

                        out_re = _mm256_fmadd_ps(recDmrs_oddRe[idxAnt][idxSc], genDmrs_re[idxSc], bd_re);
                        out_im = _mm256_fmsub_ps(recDmrs_oddIm[idxAnt][idxSc], genDmrs_re[idxSc], ad_im);

                        lsEst_oddRe[idxAnt][idxSc] = _mm256_mul_ps(out_re, scale_lsEst);
                        lsEst_oddIm[idxAnt][idxSc] = _mm256_mul_ps(out_im, scale_lsEst);
                    }
                }

                //Time offset estimation (TOE)
                toe_oddTones = 0;

                toe_perAnt_re = _mm256_setzero_ps();
                toe_perAnt_im = _mm256_setzero_ps();

                for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                {
                    toe_perPrb_re = _mm256_setzero_ps();
                    toe_perPrb_im = _mm256_setzero_ps();

                    if(nOddTonesLyr == 1)
                    {
                        for(idxSc = 0; idxSc < 5; idxSc++)
                        {
                            bd_re = _mm256_mul_ps(lsEst_oddIm[idxAnt][idxSc+1], lsEst_oddIm[idxAnt][idxSc]);
                            ad_im = _mm256_mul_ps(lsEst_oddRe[idxAnt][idxSc+1], lsEst_oddIm[idxAnt][idxSc]);

                            out_re = _mm256_fmadd_ps(lsEst_oddRe[idxAnt][idxSc+1], lsEst_oddRe[idxAnt][idxSc], bd_re);
                            out_im = _mm256_fmsub_ps(lsEst_oddIm[idxAnt][idxSc+1], lsEst_oddRe[idxAnt][idxSc], ad_im);

                            toe_perPrb_re = _mm256_add_ps(toe_perPrb_re, out_re);
                            toe_perPrb_im = _mm256_add_ps(toe_perPrb_im, out_im);
                        }
                        toe_perPrb_re = _mm256_mul_ps(toe_perPrb_re, div5Vec);
                        toe_perPrb_im = _mm256_mul_ps(toe_perPrb_im, div5Vec);
                    }
                    else if(nOddTonesLyr == 2)
                    {
                        for(idxSc = 0; idxSc < 4; idxSc++)
                        {
                            bd_re = _mm256_mul_ps(lsEst_oddIm[idxAnt][idxSc+1], lsEst_oddIm[idxAnt][idxSc]);
                            ad_im = _mm256_mul_ps(lsEst_oddRe[idxAnt][idxSc+1], lsEst_oddIm[idxAnt][idxSc]);

                            out_re = _mm256_fmadd_ps(lsEst_oddRe[idxAnt][idxSc+1], lsEst_oddRe[idxAnt][idxSc], bd_re);
                            out_im = _mm256_fmsub_ps(lsEst_oddIm[idxAnt][idxSc+1], lsEst_oddRe[idxAnt][idxSc], ad_im);

                            toe_perPrb_re = _mm256_add_ps(toe_perPrb_re, out_re);
                            toe_perPrb_im = _mm256_add_ps(toe_perPrb_im, out_im);
                        }
                        toe_perPrb_re = _mm256_mul_ps(toe_perPrb_re, div4Vec);
                        toe_perPrb_im = _mm256_mul_ps(toe_perPrb_im, div4Vec);
                    }

                    toe_perAnt_re = _mm256_add_ps(toe_perAnt_re, toe_perPrb_re);
                    toe_perAnt_im = _mm256_add_ps(toe_perAnt_im, toe_perPrb_im);
                }

                mid_re = _mm256_permute2f128_ps (toe_perAnt_re, toe_perAnt_re, 1);
                toe_perAnt_re = _mm256_add_ps(toe_perAnt_re, mid_re);
                toe_perAnt_re = _mm256_hadd_ps(toe_perAnt_re, toe_perAnt_re);
                toe_perAnt_re = _mm256_hadd_ps(toe_perAnt_re, toe_perAnt_re);

                mid_im = _mm256_permute2f128_ps (toe_perAnt_im, toe_perAnt_im, 1);
                toe_perAnt_im = _mm256_add_ps(toe_perAnt_im, mid_im);
                toe_perAnt_im = _mm256_hadd_ps(toe_perAnt_im, toe_perAnt_im);
                toe_perAnt_im = _mm256_hadd_ps(toe_perAnt_im, toe_perAnt_im);

                toe_oddTones = (toe_perAnt_re[0] + I *  toe_perAnt_im[0]) / (float)(nRxAnt*8);

                if(nOddTonesLyr == 2)
                {
                    toe_oddTones = csqrt(toe_oddTones);
                }
            }

            if((nEvenTonesLyr > 0) && (nOddTonesLyr > 0))
            {
                toe_overAll = ((toe_evenTones) + (toe_oddTones)) * 0.5;
            }
            else if(nEvenTonesLyr > 0)
            {
                toe_overAll = toe_evenTones;
            }
            else if(nOddTonesLyr > 0)
            {
                toe_overAll = toe_oddTones;
            }

            toe_overAll = (conj(toe_overAll)/cabs(toe_overAll));

            pNrPuschOutParams->toc_overAll_re[ctPrgIdx] = creal(toe_overAll);
            pNrPuschOutParams->toc_overAll_im[ctPrgIdx] = cimag(toe_overAll);

            if(nEvenTonesLyr > 0)
            {
                //Time Offset Correction
                toe_vec[0] = 1;
                toe_vec[1] = (toe_overAll);
                toe_vec[2] = (toe_vec[1] * toe_vec[1]);
                toe_vec[3] = (toe_vec[2] * toe_vec[1]);
                toe_vec[4] = (toe_vec[3] * toe_vec[1]);
                toe_vec[5] = (toe_vec[4] * toe_vec[1]);

                toe_reVec[0] = _mm256_set1_ps(1);
                toe_reVec[1] = _mm256_set1_ps(creal(toe_vec[1]));
                toe_reVec[2] = _mm256_set1_ps(creal(toe_vec[2]));
                toe_reVec[3] = _mm256_set1_ps(creal(toe_vec[3]));
                toe_reVec[4] = _mm256_set1_ps(creal(toe_vec[4]));
                toe_reVec[5] = _mm256_set1_ps(creal(toe_vec[5]));

                toe_imVec[0] = _mm256_set1_ps(0);
                toe_imVec[1] = _mm256_set1_ps(cimag(toe_vec[1]));
                toe_imVec[2] = _mm256_set1_ps(cimag(toe_vec[2]));
                toe_imVec[3] = _mm256_set1_ps(cimag(toe_vec[3]));
                toe_imVec[4] = _mm256_set1_ps(cimag(toe_vec[4]));
                toe_imVec[5] = _mm256_set1_ps(cimag(toe_vec[5]));

                for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                {
                    for(idxSc = 0; idxSc < 6; idxSc++)
                    {
                        bd_re = _mm256_mul_ps(lsEst_evenIm[idxAnt][idxSc], toe_imVec[idxSc]);
                        ad_im = _mm256_mul_ps(lsEst_evenRe[idxAnt][idxSc], toe_imVec[idxSc]);

                        lsEst_evenRe[idxAnt][idxSc] = _mm256_fmsub_ps(lsEst_evenRe[idxAnt][idxSc], toe_reVec[idxSc], bd_re);
                        lsEst_evenIm[idxAnt][idxSc] = _mm256_fmadd_ps(lsEst_evenIm[idxAnt][idxSc], toe_reVec[idxSc], ad_im);
                    }
                }

                //Channel Estimation
                for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                {
                    float *chEst_eLyr1_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][0][idxPrb];
                    float *chEst_eLyr1_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][0][idxPrb];

                    float *chEst_eLyr2_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][1][idxPrb];
                    float *chEst_eLyr2_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][1][idxPrb];

                    h1Avg_re[idxAnt] = _mm256_setzero_ps();
                    h1Avg_im[idxAnt] = _mm256_setzero_ps();

                    for(idxSc = 0; idxSc < 6; idxSc++)
                    {
                        h1Avg_re[idxAnt] = _mm256_add_ps(h1Avg_re[idxAnt], lsEst_evenRe[idxAnt][idxSc]);
                        h1Avg_im[idxAnt] = _mm256_add_ps(h1Avg_im[idxAnt], lsEst_evenIm[idxAnt][idxSc]);
                    }

                    h1Avg_re[idxAnt] = _mm256_mul_ps(h1Avg_re[idxAnt], div6Vec);
                    h1Avg_im[idxAnt] = _mm256_mul_ps(h1Avg_im[idxAnt], div6Vec);

                    _mm256_storeu_ps(chEst_eLyr1_rePtr, h1Avg_re[idxAnt]);
                    _mm256_storeu_ps(chEst_eLyr1_imPtr, h1Avg_im[idxAnt]);

                    if(nEvenTonesLyr == 2)
                    {
                        h2Avg_re[idxAnt] = _mm256_setzero_ps();
                        h2Avg_im[idxAnt] = _mm256_setzero_ps();

                        for(idxSc = 0; idxSc < 6; idxSc++)
                        {
                            if((idxSc==0) || (idxSc==2) || (idxSc==4))
                            {
                                h2Avg_re[idxAnt] = _mm256_add_ps(h2Avg_re[idxAnt], lsEst_evenRe[idxAnt][idxSc]);
                                h2Avg_im[idxAnt] = _mm256_add_ps(h2Avg_im[idxAnt], lsEst_evenIm[idxAnt][idxSc]);
                            }
                            else
                            {
                                h2Avg_re[idxAnt] = _mm256_sub_ps(h2Avg_re[idxAnt], lsEst_evenRe[idxAnt][idxSc]);
                                h2Avg_im[idxAnt] = _mm256_sub_ps(h2Avg_im[idxAnt], lsEst_evenIm[idxAnt][idxSc]);
                            }
                        }

                        h2Avg_re[idxAnt] = _mm256_mul_ps(h2Avg_re[idxAnt], div6Vec);
                        h2Avg_im[idxAnt] = _mm256_mul_ps(h2Avg_im[idxAnt], div6Vec);

                        _mm256_storeu_ps(chEst_eLyr2_rePtr, h2Avg_re[idxAnt]);
                        _mm256_storeu_ps(chEst_eLyr2_imPtr, h2Avg_im[idxAnt]);
                    }
                }

                //Covariance Matrix (Rnn) computation
                for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                {
                    for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                    {
                        Rnn_re[idxRow][idxCol] = _mm256_setzero_ps();
                        Rnn_im[idxRow][idxCol] = _mm256_setzero_ps();
                    }
                }

                for(idxSc = 0; idxSc < 6; idxSc++)
                {
                    bd_re = _mm256_mul_ps(genDmrs_im[idxSc], toe_imVec[idxSc]);
                    ad_im = _mm256_mul_ps(genDmrs_re[idxSc], toe_imVec[idxSc]);

                    mid_re = _mm256_fmadd_ps(genDmrs_re[idxSc], toe_reVec[idxSc], bd_re);
                    mid_im = _mm256_fmsub_ps(genDmrs_im[idxSc], toe_reVec[idxSc], ad_im);

                    bd_re = _mm256_mul_ps(genDmrs2_im[idxSc], toe_imVec[idxSc]);
                    ad_im = _mm256_mul_ps(genDmrs2_re[idxSc], toe_imVec[idxSc]);

                    mid2_re = _mm256_fmadd_ps(genDmrs2_re[idxSc], toe_reVec[idxSc], bd_re);
                    mid2_im = _mm256_fmsub_ps(genDmrs2_im[idxSc], toe_reVec[idxSc], ad_im);

                    for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                    {
                        bd_re = _mm256_mul_ps(h1Avg_im[idxAnt], mid_im);
                        ad_im = _mm256_mul_ps(h1Avg_re[idxAnt], mid_im);

                        mid3_re = _mm256_fmsub_ps(h1Avg_re[idxAnt], mid_re, bd_re);
                        mid3_im = _mm256_fmadd_ps(h1Avg_im[idxAnt], mid_re, ad_im);

                        Inf_re[idxAnt] = _mm256_sub_ps(recDmrs_evenRe[idxAnt][idxSc], mid3_re);
                        Inf_im[idxAnt] = _mm256_sub_ps(recDmrs_evenIm[idxAnt][idxSc], mid3_im);

                        if(nEvenTonesLyr == 2)
                        {
                            bd_re = _mm256_mul_ps(h2Avg_im[idxAnt], mid2_im);
                            ad_im = _mm256_mul_ps(h2Avg_re[idxAnt], mid2_im);

                            mid3_re = _mm256_fmsub_ps(h2Avg_re[idxAnt], mid2_re, bd_re);
                            mid3_im = _mm256_fmadd_ps(h2Avg_im[idxAnt], mid2_re, ad_im);

                            Inf_re[idxAnt] = _mm256_sub_ps(Inf_re[idxAnt], mid3_re);
                            Inf_im[idxAnt] = _mm256_sub_ps(Inf_im[idxAnt], mid3_im);
                        }
                    }

                    for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                    {
                        for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                        {
                            bd_re = _mm256_mul_ps(Inf_im[idxRow], Inf_im[idxCol]);
                            ad_im = _mm256_mul_ps(Inf_re[idxRow], Inf_im[idxCol]);

                            mid_re = _mm256_fmadd_ps(Inf_re[idxRow], Inf_re[idxCol], bd_re);
                            mid_im = _mm256_fmsub_ps(Inf_im[idxRow], Inf_re[idxCol], ad_im);

                            Rnn_re[idxRow][idxCol] = _mm256_add_ps(Rnn_re[idxRow][idxCol], mid_re);
                            Rnn_im[idxRow][idxCol] = _mm256_add_ps(Rnn_im[idxRow][idxCol], mid_im);

                            RnnTemp_re[idxRow][idxCol] = _mm256_fmadd_ps(Inf_re[idxRow], Inf_re[idxCol], bd_re);
                            RnnTemp_im[idxRow][idxCol] = _mm256_fmsub_ps(Inf_im[idxRow], Inf_re[idxCol], ad_im);
                        }
                    }
                }

                for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                {
                    for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                    {
                        if(pNrPuschInParams->irc_enable)
                        {
                            Rnn_re[idxRow][idxCol] = _mm256_mul_ps(Rnn_re[idxRow][idxCol], div6Vec);
                            Rnn_im[idxRow][idxCol] = _mm256_mul_ps(Rnn_im[idxRow][idxCol], div6Vec);
                        }

                        for(avxPrbIdx = 0; avxPrbIdx < 8; avxPrbIdx++)
                        {
                            pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+avxPrbIdx] = Rnn_re[idxRow][idxCol][avxPrbIdx];
                            pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+avxPrbIdx] = Rnn_im[idxRow][idxCol][avxPrbIdx];
                        }
                    }
                }

                //Covariance Matrix (Rnn) computation of Null Tones
                if((pNrPuschInParams->nullTone_enable == 1) && (nOddTonesLyr == 0))
                {
                    for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                    {
                        for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                        {
                            Rnn_re[idxRow][idxCol] = _mm256_setzero_ps();
                            Rnn_im[idxRow][idxCol] = _mm256_setzero_ps();
                        }
                    }

                    for(idxSc = 0; idxSc < 6; idxSc++)
                    {
                        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                        {
                            float *recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC + startScIdx];
                            float *recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][ONE_HALF_MAX_NUM_SC + startScIdx];

                            Inf_re[idxAnt] = _mm256_set_ps( recDmrs_rePtr[prbScStartIdx[7] + idxSc],
                                                            recDmrs_rePtr[prbScStartIdx[6] + idxSc],
                                                            recDmrs_rePtr[prbScStartIdx[5] + idxSc],
                                                            recDmrs_rePtr[prbScStartIdx[4] + idxSc],
                                                            recDmrs_rePtr[prbScStartIdx[3] + idxSc],
                                                            recDmrs_rePtr[prbScStartIdx[2] + idxSc],
                                                            recDmrs_rePtr[prbScStartIdx[1] + idxSc],
                                                            recDmrs_rePtr[prbScStartIdx[0] + idxSc] );

                            Inf_im[idxAnt] = _mm256_set_ps( recDmrs_imPtr[prbScStartIdx[7] + idxSc],
                                                            recDmrs_imPtr[prbScStartIdx[6] + idxSc],
                                                            recDmrs_imPtr[prbScStartIdx[5] + idxSc],
                                                            recDmrs_imPtr[prbScStartIdx[4] + idxSc],
                                                            recDmrs_imPtr[prbScStartIdx[3] + idxSc],
                                                            recDmrs_imPtr[prbScStartIdx[2] + idxSc],
                                                            recDmrs_imPtr[prbScStartIdx[1] + idxSc],
                                                            recDmrs_imPtr[prbScStartIdx[0] + idxSc] );

                        }


                        for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                        {
                            for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                            {
                                bd_re = _mm256_mul_ps(Inf_im[idxRow], Inf_im[idxCol]);
                                ad_im = _mm256_mul_ps(Inf_re[idxRow], Inf_im[idxCol]);

                                mid_re = _mm256_fmadd_ps(Inf_re[idxRow], Inf_re[idxCol], bd_re);
                                mid_im = _mm256_fmsub_ps(Inf_im[idxRow], Inf_re[idxCol], ad_im);


                                Rnn_re[idxRow][idxCol] = _mm256_add_ps(Rnn_re[idxRow][idxCol], mid_re);
                                Rnn_im[idxRow][idxCol] = _mm256_add_ps(Rnn_im[idxRow][idxCol], mid_im);

                                RnnTemp_re[idxRow][idxCol] = _mm256_fmadd_ps(Inf_re[idxRow], Inf_re[idxCol], bd_re);
                                RnnTemp_im[idxRow][idxCol] = _mm256_fmsub_ps(Inf_im[idxRow], Inf_re[idxCol], ad_im);
                            }
                        }
                    }

                    for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                    {
                        for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                        {
                            if(pNrPuschInParams->irc_enable)
                            {
                                Rnn_re[idxRow][idxCol] = _mm256_mul_ps(Rnn_re[idxRow][idxCol], div6Vec);
                                Rnn_im[idxRow][idxCol] = _mm256_mul_ps(Rnn_im[idxRow][idxCol], div6Vec);
                            }

                            for(avxPrbIdx = 0; avxPrbIdx < 8; avxPrbIdx++)
                            {
                                pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+avxPrbIdx] = Rnn_re[idxRow][idxCol][avxPrbIdx];
                                pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+avxPrbIdx] = Rnn_im[idxRow][idxCol][avxPrbIdx];
                            }
                        }
                    }

                }

            }

            if(nOddTonesLyr > 0)
            {
                //Time Offset Correction
                toe_overAll = csqrt(toe_overAll);
                toe_vec[0] = (toe_overAll);
                toe_vec[1] = (toe_vec[0] * toe_overAll * toe_overAll);
                toe_vec[2] = (toe_vec[1] * toe_overAll * toe_overAll);
                toe_vec[3] = (toe_vec[2] * toe_overAll * toe_overAll);
                toe_vec[4] = (toe_vec[3] * toe_overAll * toe_overAll);
                toe_vec[5] = (toe_vec[4] * toe_overAll * toe_overAll);

                toe_reVec[0] = _mm256_set1_ps(creal(toe_vec[0]));
                toe_reVec[1] = _mm256_set1_ps(creal(toe_vec[1]));
                toe_reVec[2] = _mm256_set1_ps(creal(toe_vec[2]));
                toe_reVec[3] = _mm256_set1_ps(creal(toe_vec[3]));
                toe_reVec[4] = _mm256_set1_ps(creal(toe_vec[4]));
                toe_reVec[5] = _mm256_set1_ps(creal(toe_vec[5]));

                toe_imVec[0] = _mm256_set1_ps(cimag(toe_vec[0]));
                toe_imVec[1] = _mm256_set1_ps(cimag(toe_vec[1]));
                toe_imVec[2] = _mm256_set1_ps(cimag(toe_vec[2]));
                toe_imVec[3] = _mm256_set1_ps(cimag(toe_vec[3]));
                toe_imVec[4] = _mm256_set1_ps(cimag(toe_vec[4]));
                toe_imVec[5] = _mm256_set1_ps(cimag(toe_vec[5]));

                for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                {
                    for(idxSc = 0; idxSc < 6; idxSc++)
                    {
                        bd_re = _mm256_mul_ps(lsEst_oddIm[idxAnt][idxSc], toe_imVec[idxSc]);
                        ad_im = _mm256_mul_ps(lsEst_oddRe[idxAnt][idxSc], toe_imVec[idxSc]);

                        lsEst_oddRe[idxAnt][idxSc] = _mm256_fmsub_ps(lsEst_oddRe[idxAnt][idxSc], toe_reVec[idxSc], bd_re);
                        lsEst_oddIm[idxAnt][idxSc] = _mm256_fmadd_ps(lsEst_oddIm[idxAnt][idxSc], toe_reVec[idxSc], ad_im);
                    }
                }

                //Channel Estimation
                for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                {
                    float *chEst_eLyr1_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][2][idxPrb];
                    float *chEst_eLyr1_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][2][idxPrb];

                    float *chEst_eLyr2_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][3][idxPrb];
                    float *chEst_eLyr2_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][3][idxPrb];

                    h1Avg_re[idxAnt] = _mm256_setzero_ps();
                    h1Avg_im[idxAnt] = _mm256_setzero_ps();

                    for(idxSc = 0; idxSc < 6; idxSc++)
                    {
                        h1Avg_re[idxAnt] = _mm256_add_ps(h1Avg_re[idxAnt], lsEst_oddRe[idxAnt][idxSc]);
                        h1Avg_im[idxAnt] = _mm256_add_ps(h1Avg_im[idxAnt], lsEst_oddIm[idxAnt][idxSc]);
                    }

                    h1Avg_re[idxAnt] = _mm256_mul_ps(h1Avg_re[idxAnt], div6Vec);
                    h1Avg_im[idxAnt] = _mm256_mul_ps(h1Avg_im[idxAnt], div6Vec);

                    _mm256_storeu_ps(chEst_eLyr1_rePtr, h1Avg_re[idxAnt]);
                    _mm256_storeu_ps(chEst_eLyr1_imPtr, h1Avg_im[idxAnt]);

                    if(nOddTonesLyr == 2)
                    {
                        h2Avg_re[idxAnt] = _mm256_setzero_ps();
                        h2Avg_im[idxAnt] = _mm256_setzero_ps();

                        for(idxSc = 0; idxSc < 6; idxSc++)
                        {
                            if((idxSc==0) || (idxSc==2) || (idxSc==4))
                            {
                                h2Avg_re[idxAnt] = _mm256_add_ps(h2Avg_re[idxAnt], lsEst_oddRe[idxAnt][idxSc]);
                                h2Avg_im[idxAnt] = _mm256_add_ps(h2Avg_im[idxAnt], lsEst_oddIm[idxAnt][idxSc]);
                            }
                            else
                            {
                                h2Avg_re[idxAnt] = _mm256_sub_ps(h2Avg_re[idxAnt], lsEst_oddRe[idxAnt][idxSc]);
                                h2Avg_im[idxAnt] = _mm256_sub_ps(h2Avg_im[idxAnt], lsEst_oddIm[idxAnt][idxSc]);
                            }
                        }

                        h2Avg_re[idxAnt] = _mm256_mul_ps(h2Avg_re[idxAnt], div6Vec);
                        h2Avg_im[idxAnt] = _mm256_mul_ps(h2Avg_im[idxAnt], div6Vec);

                        _mm256_storeu_ps(chEst_eLyr2_rePtr, h2Avg_re[idxAnt]);
                        _mm256_storeu_ps(chEst_eLyr2_imPtr, h2Avg_im[idxAnt]);
                    }
                }

                //Covariance Matrix (Rnn) computation
                for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                {
                    for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                    {
                        Rnn_re[idxRow][idxCol] = _mm256_setzero_ps();
                        Rnn_im[idxRow][idxCol] = _mm256_setzero_ps();
                    }
                }

                for(idxSc = 0; idxSc < 6; idxSc++)
                {
                    bd_re = _mm256_mul_ps(genDmrs_im[idxSc], toe_imVec[idxSc]);
                    ad_im = _mm256_mul_ps(genDmrs_re[idxSc], toe_imVec[idxSc]);

                    mid_re = _mm256_fmadd_ps(genDmrs_re[idxSc], toe_reVec[idxSc], bd_re);
                    mid_im = _mm256_fmsub_ps(genDmrs_im[idxSc], toe_reVec[idxSc], ad_im);

                    bd_re = _mm256_mul_ps(genDmrs2_im[idxSc], toe_imVec[idxSc]);
                    ad_im = _mm256_mul_ps(genDmrs2_re[idxSc], toe_imVec[idxSc]);

                    mid2_re = _mm256_fmadd_ps(genDmrs2_re[idxSc], toe_reVec[idxSc], bd_re);
                    mid2_im = _mm256_fmsub_ps(genDmrs2_im[idxSc], toe_reVec[idxSc], ad_im);

                    for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                    {
                        bd_re = _mm256_mul_ps(h1Avg_im[idxAnt], mid_im);
                        ad_im = _mm256_mul_ps(h1Avg_re[idxAnt], mid_im);

                        mid3_re = _mm256_fmsub_ps(h1Avg_re[idxAnt], mid_re, bd_re);
                        mid3_im = _mm256_fmadd_ps(h1Avg_im[idxAnt], mid_re, ad_im);

                        Inf_re[idxAnt] = _mm256_sub_ps(recDmrs_oddRe[idxAnt][idxSc], mid3_re);
                        Inf_im[idxAnt] = _mm256_sub_ps(recDmrs_oddIm[idxAnt][idxSc], mid3_im);

                        if(nEvenTonesLyr == 2)
                        {
                            bd_re = _mm256_mul_ps(h2Avg_im[idxAnt], mid2_im);
                            ad_im = _mm256_mul_ps(h2Avg_re[idxAnt], mid2_im);

                            mid3_re = _mm256_fmsub_ps(h2Avg_re[idxAnt], mid2_re, bd_re);
                            mid3_im = _mm256_fmadd_ps(h2Avg_im[idxAnt], mid2_re, ad_im);

                            Inf_re[idxAnt] = _mm256_sub_ps(Inf_re[idxAnt], mid3_re);
                            Inf_im[idxAnt] = _mm256_sub_ps(Inf_im[idxAnt], mid3_im);
                        }
                    }

                    for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                    {
                        for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                        {
                            bd_re = _mm256_mul_ps(Inf_im[idxRow], Inf_im[idxCol]);
                            ad_im = _mm256_mul_ps(Inf_re[idxRow], Inf_im[idxCol]);

                            mid_re = _mm256_fmadd_ps(Inf_re[idxRow], Inf_re[idxCol], bd_re);
                            mid_im = _mm256_fmsub_ps(Inf_im[idxRow], Inf_re[idxCol], ad_im);


                            Rnn_re[idxRow][idxCol] = _mm256_add_ps(Rnn_re[idxRow][idxCol], mid_re);
                            Rnn_im[idxRow][idxCol] = _mm256_add_ps(Rnn_im[idxRow][idxCol], mid_im);
                        }
                    }
                }

                for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                {
                    for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                    {
                        if(pNrPuschInParams->irc_enable)
                        {
                            Rnn_re[idxRow][idxCol] = _mm256_mul_ps(Rnn_re[idxRow][idxCol], div6Vec);
                            Rnn_im[idxRow][idxCol] = _mm256_mul_ps(Rnn_im[idxRow][idxCol], div6Vec);
                        }


                        for(avxPrbIdx = 0; avxPrbIdx < 8; avxPrbIdx++)
                        {
                            pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+avxPrbIdx] = Rnn_re[idxRow][idxCol][avxPrbIdx];
                            pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+avxPrbIdx] = Rnn_im[idxRow][idxCol][avxPrbIdx];
                        }
                    }
                }

                //Covariance Matrix (Rnn) computation of Null Tones
                if((pNrPuschInParams->nullTone_enable == 1) && (nEvenTonesLyr == 0))
                {
                    for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                    {
                        for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                        {
                            Rnn_re[idxRow][idxCol] = _mm256_setzero_ps();
                            Rnn_im[idxRow][idxCol] = _mm256_setzero_ps();
                        }
                    }

                    for(idxSc = 0; idxSc < 6; idxSc++)
                    {
                        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                        {
                            float *recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                            float *recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][HALF_MAX_NUM_SC + startScIdx];

                            Inf_re[idxAnt] = _mm256_set_ps( recDmrs_rePtr[prbScStartIdx[7] + idxSc],
                                                            recDmrs_rePtr[prbScStartIdx[6] + idxSc],
                                                            recDmrs_rePtr[prbScStartIdx[5] + idxSc],
                                                            recDmrs_rePtr[prbScStartIdx[4] + idxSc],
                                                            recDmrs_rePtr[prbScStartIdx[3] + idxSc],
                                                            recDmrs_rePtr[prbScStartIdx[2] + idxSc],
                                                            recDmrs_rePtr[prbScStartIdx[1] + idxSc],
                                                            recDmrs_rePtr[prbScStartIdx[0] + idxSc] );

                            Inf_im[idxAnt] = _mm256_set_ps( recDmrs_imPtr[prbScStartIdx[7] + idxSc],
                                                            recDmrs_imPtr[prbScStartIdx[6] + idxSc],
                                                            recDmrs_imPtr[prbScStartIdx[5] + idxSc],
                                                            recDmrs_imPtr[prbScStartIdx[4] + idxSc],
                                                            recDmrs_imPtr[prbScStartIdx[3] + idxSc],
                                                            recDmrs_imPtr[prbScStartIdx[2] + idxSc],
                                                            recDmrs_imPtr[prbScStartIdx[1] + idxSc],
                                                            recDmrs_imPtr[prbScStartIdx[0] + idxSc] );

                        }


                        for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                        {
                            for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                            {
                                bd_re = _mm256_mul_ps(Inf_im[idxRow], Inf_im[idxCol]);
                                ad_im = _mm256_mul_ps(Inf_re[idxRow], Inf_im[idxCol]);

                                mid_re = _mm256_fmadd_ps(Inf_re[idxRow], Inf_re[idxCol], bd_re);
                                mid_im = _mm256_fmsub_ps(Inf_im[idxRow], Inf_re[idxCol], ad_im);

                                Rnn_re[idxRow][idxCol] = _mm256_add_ps(Rnn_re[idxRow][idxCol], mid_re);
                                Rnn_im[idxRow][idxCol] = _mm256_add_ps(Rnn_im[idxRow][idxCol], mid_im);

                                RnnTemp_re[idxRow][idxCol] = _mm256_fmadd_ps(Inf_re[idxRow], Inf_re[idxCol], bd_re);
                                RnnTemp_im[idxRow][idxCol] = _mm256_fmsub_ps(Inf_im[idxRow], Inf_re[idxCol], ad_im);
                            }
                        }

                    }

                    for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                    {
                        for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                        {
                            if(pNrPuschInParams->irc_enable)
                            {
                                Rnn_re[idxRow][idxCol] = _mm256_mul_ps(Rnn_re[idxRow][idxCol], div6Vec);
                                Rnn_im[idxRow][idxCol] = _mm256_mul_ps(Rnn_im[idxRow][idxCol], div6Vec);
                            }

                            for(avxPrbIdx = 0; avxPrbIdx < 8; avxPrbIdx++)
                            {
                                pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+avxPrbIdx] = Rnn_re[idxRow][idxCol][avxPrbIdx];
                                pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+avxPrbIdx] = Rnn_im[idxRow][idxCol][avxPrbIdx];
                            }
                        }
                    }

                }


            }

        }

        if((pNrPuschInParams->nullTone_enable == 1) || ((pNrPuschOutParams->nEvenTonesLyr > 0) && (pNrPuschOutParams->nOddTonesLyr > 0)))
        {
            for(idxPrb = startPrbIdx; idxPrb <= endPrbIdx; idxPrb++)
            {
                for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                {
                    for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                    {
                        pNrPuschOutParams->Rnn_prb_re[idxRow][idxCol][idxPrb] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb] +
                                                                                  pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb] ) * 0.5;
                        pNrPuschOutParams->Rnn_prb_im[idxRow][idxCol][idxPrb] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb] +
                                                                                  pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb] ) * 0.5;
                    }
                }
            }
        }

    }
}

void nr_pusch__perPrb_nVarCalc_avx2( int symPart,
                                     commonUlConfigT* pNrUlCommonParams,
                                     puschConfigT* pNrPuschInParams,
                                     P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams )
{
    int idxAnt, idxAlloc, idxPrb, ctPrgIdx;
    int startPrbIdx, endPrbIdx;

    int nRxAnt, nLyr, nSbAlloc, nPrb_alloc;

    __m256 real, imag;
    __m256 nVar, nVar_out, mid_val;
    float nVar_overAll;

    nRxAnt   = pNrPuschInParams->nAntennaPorts;
    nLyr     = pNrPuschInParams->nNrOfLayers;
    nSbAlloc = pNrPuschInParams->nSbAlloc;

    int nEvenTonesLyr = pNrPuschOutParams->nEvenTonesLyr;
    int nOddTonesLyr = pNrPuschOutParams->nOddTonesLyr;

    for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
    {
        nPrb_alloc = pNrPuschInParams->nPrb[idxAlloc];

        startPrbIdx = pNrPuschInParams->rbStart[idxAlloc];
        endPrbIdx = startPrbIdx + nPrb_alloc;

        if((nEvenTonesLyr > 0) && (nOddTonesLyr > 0))
        {
            nVar = _mm256_setzero_ps();
            ctPrgIdx = 0;
            for(idxPrb = startPrbIdx; idxPrb < endPrbIdx; idxPrb += 8, ctPrgIdx++)
            {
                for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                {
                    real = _mm256_loadu_ps(&pNrPuschOutParams->Rnn_prb_re[idxAnt][idxAnt][idxPrb]);
                    imag = _mm256_loadu_ps(&pNrPuschOutParams->Rnn_prb_im[idxAnt][idxAnt][idxPrb]);

                    nVar = _mm256_add_ps(nVar, real);
                }

                mid_val = _mm256_permute2f128_ps (nVar, nVar, 1);
                nVar_out = _mm256_add_ps(nVar, mid_val);
                nVar_out = _mm256_hadd_ps(nVar_out, nVar_out);
                nVar_out = _mm256_hadd_ps(nVar_out, nVar_out);

                nVar_overAll = nVar_out[0] / (1632 * nRxAnt);

                pNrPuschOutParams->nVar_overAll[ctPrgIdx] = nVar_overAll;
            }
        }
        else if(nEvenTonesLyr > 0)
        {
            nVar = _mm256_setzero_ps();
            ctPrgIdx = 0;
            for(idxPrb = startPrbIdx; idxPrb < endPrbIdx; idxPrb += 8, ctPrgIdx++)
            {

                for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                {
                    real = _mm256_loadu_ps(&pNrPuschOutParams->RnnEven_prb_re[idxAnt][idxAnt][idxPrb]);

                    nVar = _mm256_add_ps(nVar, real);
                }

                mid_val = _mm256_permute2f128_ps (nVar, nVar, 1);
                nVar_out = _mm256_add_ps(nVar, mid_val);
                nVar_out = _mm256_hadd_ps(nVar_out, nVar_out);
                nVar_out = _mm256_hadd_ps(nVar_out, nVar_out);

                nVar_overAll = nVar_out[0] / (1632 * nRxAnt);

                pNrPuschOutParams->nVar_overAll[ctPrgIdx] = nVar_overAll;
            }
        }
        else if(nOddTonesLyr > 0)
        {
            for(idxPrb = startPrbIdx; idxPrb < endPrbIdx; idxPrb++)
            {
                nVar = _mm256_setzero_ps();
                for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                {
                    real = _mm256_loadu_ps(&pNrPuschOutParams->RnnOdd_prb_re[idxAnt][idxAnt][idxPrb]);
                    imag = _mm256_loadu_ps(&pNrPuschOutParams->RnnOdd_prb_im[idxAnt][idxAnt][idxPrb]);

                    real = _mm256_mul_ps(real, real);
                    imag = _mm256_mul_ps(imag, imag);

                    nVar = _mm256_add_ps(nVar, real);
                    nVar = _mm256_add_ps(nVar, imag);
                }

                mid_val = _mm256_permute2f128_ps (nVar, nVar, 1);
                nVar = _mm256_add_ps(nVar, mid_val);
                nVar = _mm256_hadd_ps(nVar, nVar);
                nVar = _mm256_hadd_ps(nVar, nVar);

                nVar_overAll = nVar[0] / (nRxAnt);

                pNrPuschOutParams->nVar_overAll[ctPrgIdx] = nVar_overAll;
            }
        }
    }
}

void nr_pusch_mmse_perPrb_equ_avx2( int symPart,
                                    commonUlConfigT* pNrUlCommonParams,
                                    puschConfigT* pNrPuschInParams,
                                    P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams)
{
    nr_pusch__perPrb_nVarCalc_avx2(symPart, pNrUlCommonParams, pNrPuschInParams, pNrPuschOutParams);

    int idxAnt, idxLyr, idxAlloc, idxPrb, idxSym, idxSc, rIndx, cIndx, i, idx;
    int startScIdx, endScIdx, remSc;
    int startPrbIdx, endPrbIdx, remPrb, CtPrbScStartIdx, avxPrbIdx, ctPrgIdx, lyrPort, lyrType;

    int startSym, endSym, dmrsSym;

    float zeroBuf[6] = {0, 0, 0, 0, 0, 0};

    int prbScStartIdx_re[8], prbScStartIdx_im[8];

    int nRxAnt, nLyr, nSbAlloc, nPrb_alloc, nPrb_total;

    char prbTemp;

    g_negate = VECT_SET1(-0.0);

    nRxAnt   = pNrPuschInParams->nAntennaPorts;
    nLyr     = pNrPuschInParams->nNrOfLayers;
    nSbAlloc = pNrPuschInParams->nSbAlloc;

    VECT_T H_re[nLyr][nRxAnt];
	VECT_T H_im[nLyr][nRxAnt];
	VECT_T HHh_re[nLyr][nLyr];
	VECT_T HHh_im[nLyr][nLyr];
	VECT_T invHHh_re[nLyr][nLyr];
	VECT_T invHHh_im[nLyr][nLyr];
	VECT_T weights_re[nRxAnt][nLyr];
	VECT_T weights_im[nRxAnt][nLyr];
	VECT_T Y_vect_re[1][nRxAnt];
	VECT_T Y_vect_im[1][nRxAnt];
	VECT_T Y_tilda_re[1][nLyr];
	VECT_T Y_tilda_im[1][nLyr];
	VECT_T Y_cap_re;
	VECT_T Y_cap_im;

	VECT_T nVar_vec;
	VECT_T sq_toc_re[2][12];
	VECT_T sq_toc_im[2][12];

    float SNRPerPrb;
    wnUInt16 mseCnt = 0;
    float sum = 0;

    for(rIndx = 0; rIndx < nLyr; rIndx++)
    {
        pNrPuschOutParams->oneBymse_avg_re[rIndx] = 0;
    }

    int total_Prb = 0;
    float one_By_total_Prb = 0;

    for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
    {
        total_Prb = total_Prb + pNrPuschInParams->nPrb[idxAlloc];
    }
    one_By_total_Prb = 1.0/total_Prb;

	float complex toc_overAll, sqrt_toc[12];

    if(symPart == LOW_PUSCH_SYM)
    {
        dmrsSym = pNrPuschInParams->dmrs1_sym;
        startSym = 0;
        endSym = 6;
    }
    else
    {
        dmrsSym = pNrPuschInParams->dmrs2_sym;
        startSym = 7;
        endSym = 13;
    }

    nPrb_total = 0;

    for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
    {
        nPrb_alloc = pNrPuschInParams->nPrb[idxAlloc];

        nPrb_total += nPrb_alloc;

        startPrbIdx = pNrPuschInParams->rbStart[idxAlloc];
        endPrbIdx   = startPrbIdx + (pNrPuschInParams->nPrb[idxAlloc]) - 1;
        CtPrbScStartIdx = startPrbIdx * 12;

        ctPrgIdx = 0;
        for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb+=8, CtPrbScStartIdx += 12, ctPrgIdx++)
        {
            prbScStartIdx_re[0] = idxPrb * 12;
            prbScStartIdx_im[0] = MAX_NUM_SC + prbScStartIdx_re[0];

            for(idx = 1; idx < 8; idx++)
            {
                prbScStartIdx_re[idx] = prbScStartIdx_re[idx-1] + 12;
                prbScStartIdx_im[idx] = prbScStartIdx_im[idx-1] + 12;
            }

           for(rIndx = 0; rIndx < nLyr; rIndx++)
            {
                for(cIndx = 0; cIndx < nRxAnt; cIndx++)
                {
                    lyrPort = pNrPuschInParams->nPortIndex[rIndx];
                    H_re[rIndx][cIndx] = _mm256_loadu_ps(&pNrPuschOutParams->chEst_re[cIndx][lyrPort][idxPrb]);
                    H_im[rIndx][cIndx] = _mm256_loadu_ps(&pNrPuschOutParams->chEst_im[cIndx][lyrPort][idxPrb]);
                }
            }

            if((idxPrb + 8) >= endPrbIdx)
            {
                remPrb = endPrbIdx - idxPrb;

                for(rIndx = 0; rIndx < nLyr; rIndx++)
                {
                    for(cIndx = 0; cIndx < nRxAnt; cIndx++)
                    {
                        for(idx = remPrb; idx < 8; idx++)
                        {
                            H_re[rIndx][cIndx][idx] = 1;
                            H_im[rIndx][cIndx][idx] = 1;
                        }
                    }
                }
            }
            else
            {
                remPrb = 8;
            }

            // HH'
            MATRIXMULTTH(H_re, H_im, HHh_re, HHh_im, nLyr, nRxAnt, nLyr);

            nVar_vec = VECT_SET1(pNrPuschOutParams->nVar_overAll[ctPrgIdx]);

            for(idx = 0; idx < nLyr; idx++)
            {
                HHh_re[idx][idx] = VECT_ADD(nVar_vec, HHh_re[idx][idx]);
            }

            // Inverse of (HH' + (epsilon)I)
            if(nLyr == 1)
            {
                InverseMatrix1x1HTranspose(HHh_re, HHh_im, invHHh_re, invHHh_im);
            }
            else if(nLyr == 2)
            {
                InverseMatrix2x2H(HHh_re, HHh_im, invHHh_re, invHHh_im);
            }
            else if(nLyr == 4)
            {
                InverseMatrix4x4HTranspose(HHh_re, HHh_im, invHHh_re, invHHh_im);
            }

            // for(idx = 0; idx < remPrb; idx++)
            // {
            //     for(rIndx = 0; rIndx < nLyr; rIndx++)
            //     {
            //         pNrPuschOutParams->msePerPrb[rIndx][idxPrb + idx] = invHHh_re[rIndx][rIndx][idx];
            //     }
            // }
            
            for(idx = 0; idx < remPrb; idx++)
            {
                for(rIndx = 0; rIndx < nLyr; rIndx++)
                {
                    SNRPerPrb = 1/invHHh_re[rIndx][rIndx][idx];
                    for(wnUInt8 toneIdx = 0;toneIdx<12;toneIdx++)
                    {
                        pNrPuschOutParams->oneBymse_re[rIndx][mseCnt+toneIdx] = SNRPerPrb;
                    }
                    sum = pNrPuschOutParams->oneBymse_avg_re[rIndx];
                    sum = sum + invHHh_re[rIndx][rIndx][idx];
                    pNrPuschOutParams->oneBymse_avg_re[rIndx] = sum;
                }
                mseCnt += 12;
            }

            // Multiplication of H' and inverse of (HH' + (epsilon)I)
            MATRIXMULTHTTRANSPOSE2(H_re, H_im, invHHh_re, invHHh_im, weights_re, weights_im, nRxAnt, nLyr, nLyr);

            if(pNrPuschOutParams->nEvenTonesLyr > 0)
            {
                toc_overAll = pNrPuschOutParams->toc_overAll_re[ctPrgIdx] + I * pNrPuschOutParams->toc_overAll_im[ctPrgIdx];

                sqrt_toc[0] = 1;
                sqrt_toc[1] = csqrtf( toc_overAll );

                sq_toc_re[EVEN_PUSCH_LYR][0] = VECT_SET1(1);
                sq_toc_im[EVEN_PUSCH_LYR][0] = VECT_SET1(0);

                sq_toc_re[EVEN_PUSCH_LYR][1] = VECT_SET1(creal(sqrt_toc[1]));
                sq_toc_im[EVEN_PUSCH_LYR][1] = VECT_SET1(cimag(sqrt_toc[1]));

                for(idx = 2; idx < 12; idx++)
                {
                    sqrt_toc[idx] = sqrt_toc[idx -1] * sqrt_toc[1];

                    sq_toc_re[EVEN_PUSCH_LYR][idx] = VECT_SET1(creal(sqrt_toc[idx]));
                    sq_toc_im[EVEN_PUSCH_LYR][idx] = VECT_SET1(cimag(sqrt_toc[idx]));
                }
            }

            if(pNrPuschOutParams->nOddTonesLyr > 0)
            {
                toc_overAll = pNrPuschOutParams->toc_overAll_re[ctPrgIdx] + I * pNrPuschOutParams->toc_overAll_im[ctPrgIdx];

                sqrt_toc[0] = 1;
                sqrt_toc[1] = csqrt( toc_overAll );

                sq_toc_re[ODD_PUSCH_LYR][0] = VECT_SET1(1);
                sq_toc_im[ODD_PUSCH_LYR][0] = VECT_SET1(0);

                sq_toc_re[ODD_PUSCH_LYR][1] = VECT_SET1(creal(sqrt_toc[1]));
                sq_toc_im[ODD_PUSCH_LYR][1] = VECT_SET1(cimag(sqrt_toc[1]));

                for(idx = 2; idx < 12; idx++)
                {
                    sqrt_toc[idx] = sqrt_toc[idx -1] * sqrt_toc[1];

                    sq_toc_re[ODD_PUSCH_LYR][idx] = VECT_SET1(creal(sqrt_toc[idx]));
                    sq_toc_im[ODD_PUSCH_LYR][idx] = VECT_SET1(cimag(sqrt_toc[idx]));
                }
            }

            int equaSym = 0;
            for(idxSym = startSym; idxSym <= endSym; idxSym++)
            {
                if(idxSym != dmrsSym)
                {
                    for(idxSc = 0; idxSc < 12; idxSc++)
                    {
                        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                        {
                            Y_vect_re[0][idxAnt] = VECT_SET(rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[7] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[6] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[5] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[4] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[3] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[2] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[1] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[0] + idxSc] );

                            Y_vect_im[0][idxAnt] = VECT_SET(rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[7] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[6] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[5] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[4] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[3] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[2] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[1] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[0] + idxSc] );
                        }

                        MATRIXMULNN(Y_vect_re, Y_vect_im, weights_re, weights_im, Y_tilda_re, Y_tilda_im, 1, nRxAnt, nLyr);

                        for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
                        {
                            lyrPort = pNrPuschInParams->nPortIndex[idxLyr];
                            if(lyrPort < 2)
                            {
                                lyrType = EVEN_PUSCH_LYR;
                            }
                            else
                            {
                                lyrType = ODD_PUSCH_LYR;
                            }

                            C_VECT_MUL_NA_NB(Y_tilda_re[0][idxLyr], Y_tilda_im[0][idxLyr], sq_toc_re[lyrType][idxSc], sq_toc_im[lyrType][idxSc], Y_cap_re, Y_cap_im);

                            for(idx = 0; idx < remPrb; idx++)
                            {
                                pNrPuschOutParams->equaOutSamples_re[idxLyr][equaSym][prbScStartIdx_re[idx] + idxSc] = Y_cap_re[idx];
                                pNrPuschOutParams->equaOutSamples_im[idxLyr][equaSym][prbScStartIdx_re[idx] + idxSc] = Y_cap_im[idx];
                            }
                        }
                    }
                    equaSym = equaSym + 1;
                }
            }
        }
    }

    for(rIndx = 0; rIndx < nLyr; rIndx++)
    {
        pNrPuschOutParams->oneBymse_avg_re[rIndx] = \
            pNrPuschOutParams->oneBymse_avg_re[rIndx]*one_By_total_Prb;
    }
}//*/

// } perPRB End



/*//{ Per HPrb Start
void nr_pusch_mmse_perHPrb_modRemoval_avx2( int symPart,
                               commonUlConfigT* pNrUlCommonParams,
                               puschConfigT* pNrPuschInParams,
                               P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams)
{
    int idxAnt, idxLyr, idxAlloc, idxSc, i;
    int startScIdx, endScIdx, remSc;

    int dmrsSym, nRxAnt, nLyr, nSbAlloc, nPrb_alloc, nPrb_total;

    __m256 rDmrs_re, rDmrs_im;
    __m256 gDmrs_re, gDmrs_im;
    __m256 bd_re, ad_im;
    __m256 out_re, out_im;

    __m256 alternateSigns = _mm256_setr_ps(1, -1, 1, -1, 1, -1, 1, -1);
    float* alternateSigns_f = (float *)&alternateSigns;

    __m256 zeroVec = _mm256_setzero_ps();
    __m256 div2Vec = _mm256_set1_ps(0.5);

    float* rDmrs_re_f = (float*)&rDmrs_re;
    float* rDmrs_im_f = (float*)&rDmrs_im;

    float* gDmrs_re_f = (float*)&gDmrs_re;
    float* gDmrs_im_f = (float*)&gDmrs_im;

    if(symPart == LOW_PUSCH_SYM)
    {
        dmrsSym = pNrPuschInParams->dmrs1_sym;
    }
    else
    {
        dmrsSym = pNrPuschInParams->dmrs2_sym;
    }

    nRxAnt   = pNrPuschInParams->nAntennaPorts;
    nLyr     = pNrPuschInParams->nNrOfLayers;
    nSbAlloc = pNrPuschInParams->nSbAlloc;

    nPrb_total = 0;
    for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
    {
        nPrb_total += pNrPuschInParams->nPrb[idxAlloc];
    }

    pNrPuschOutParams->total_nPRB = nPrb_total;

    for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
    {
        for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
        {
            for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
            {
                startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                if(pNrPuschInParams->nPortIndex[idxLyr] == 0)
                {
                    float *recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                    float *recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][(MAX_NUM_SC/2) + startScIdx];

                    float *genDmrs_rePtr = &pNrPuschOutParams->genDMRS[startScIdx];
                    float *genDmrs_imPtr = &pNrPuschOutParams->genDMRS[(MAX_NUM_SC/2) + startScIdx];
// printf ("%f%+fi\t%f%+fi\t%d\n",recDmrs_rePtr[0],recDmrs_imPtr[0],genDmrs_rePtr[0],genDmrs_imPtr[0],dmrsSym);
                    float *lsEst_rePtr = &pNrPuschOutParams->lsEst[idxAnt][idxLyr][startScIdx];
                    float *lsEst_imPtr =  &pNrPuschOutParams->lsEst[idxAnt][idxLyr][(MAX_NUM_SC/2) + startScIdx];

                    for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                    {
                        if((idxSc + 8) >= (endScIdx))
                        {
                            remSc = (endScIdx) - idxSc;

                            rDmrs_re = _mm256_setzero_ps();
                            rDmrs_im = _mm256_setzero_ps();

                            gDmrs_re = _mm256_setzero_ps();
                            gDmrs_im = _mm256_setzero_ps();

                            for(i = 0; i < remSc; i++)
                            {
                                rDmrs_re_f[i] = recDmrs_rePtr[i];
                                rDmrs_im_f[i] = recDmrs_imPtr[i];

                                gDmrs_re_f[i] = genDmrs_rePtr[i];
                                gDmrs_im_f[i] = genDmrs_imPtr[i];
                            }
                        }
                        else
                        {
                            rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                            rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);

                            gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                            gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                        }

                        bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                        ad_im = _mm256_mul_ps(rDmrs_re, gDmrs_im);

                        out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                        out_im = _mm256_fmsub_ps(rDmrs_im, gDmrs_re, ad_im);

                        out_re = _mm256_mul_ps(out_re, div2Vec);
                        out_im = _mm256_mul_ps(out_im, div2Vec);

                        if((idxSc + 8) >= (endScIdx))
                        {
                            remSc = (endScIdx) - idxSc;

                            for(i = 0; i < remSc; i++)
                            {
                                lsEst_rePtr[i] = out_re[i];
                                lsEst_imPtr[i] = out_im[i];
                            }

                        }
                        else
                        {
                            _mm256_storeu_ps(lsEst_rePtr, out_re);
                            _mm256_storeu_ps(lsEst_imPtr, out_im);
                        }

                        recDmrs_rePtr += 8;
                        recDmrs_imPtr += 8;

                        genDmrs_rePtr += 8;
                        genDmrs_imPtr += 8;

                        lsEst_rePtr += 8;
                        lsEst_imPtr += 8;

                    }   //for idxSc
                }
                else if(pNrPuschInParams->nPortIndex[idxLyr] == 1)
                {
                    float *recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                    float *recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][(MAX_NUM_SC/2) + startScIdx];

                    float *genDmrs_rePtr = &pNrPuschOutParams->genDMRS[startScIdx];
                    float *genDmrs_imPtr = &pNrPuschOutParams->genDMRS[(MAX_NUM_SC/2) + startScIdx];

                    float *lsEst_rePtr = &pNrPuschOutParams->lsEst[idxAnt][idxLyr][startScIdx];
                    float *lsEst_imPtr = &pNrPuschOutParams->lsEst[idxAnt][idxLyr][(MAX_NUM_SC/2) + startScIdx];

                    for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                    {
                        if((idxSc + 8) >= endScIdx)
                        {
                            remSc = endScIdx - idxSc;

                            rDmrs_re = _mm256_setzero_ps();
                            rDmrs_im = _mm256_setzero_ps();

                            gDmrs_re = _mm256_setzero_ps();
                            gDmrs_im = _mm256_setzero_ps();

                            for(i = 0; i < remSc; i++)
                            {
                                rDmrs_re_f[i] = recDmrs_rePtr[i];
                                rDmrs_im_f[i] = recDmrs_imPtr[i];

                                gDmrs_re_f[i] = alternateSigns_f[i]*genDmrs_rePtr[i];
                                gDmrs_im_f[i] = alternateSigns_f[i]*genDmrs_imPtr[i];
                            }
                        }
                        else
                        {
                            rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                            rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);

                            // gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                            // gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);

                            gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                            gDmrs_re = _mm256_mul_ps(gDmrs_re, alternateSigns);
                            gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                            gDmrs_im = _mm256_mul_ps(gDmrs_im, alternateSigns);
                        }

                        // bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                        // ad_im = _mm256_mul_ps(rDmrs_im, gDmrs_re);

                        // out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                        // out_im = _mm256_fmsub_ps(rDmrs_re, gDmrs_im, ad_im);

                        bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                        ad_im = _mm256_mul_ps(rDmrs_re, gDmrs_im);

                        out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                        out_im = _mm256_fmsub_ps(rDmrs_im, gDmrs_re, ad_im);

                        out_re = _mm256_mul_ps(out_re, div2Vec);
                        out_im = _mm256_mul_ps(out_im, div2Vec);

                        if((idxSc + 8) >= (endScIdx))
                        {
                            remSc = (endScIdx) - idxSc;

                            for(i = 0; i < remSc; i++)
                            {
                                lsEst_rePtr[i] = out_re[i];
                                lsEst_imPtr[i] = out_im[i];
                            }

                        }
                        else
                        {
                            _mm256_storeu_ps(lsEst_rePtr, out_re);
                            _mm256_storeu_ps(lsEst_imPtr, out_im);
                        }

                        recDmrs_rePtr += 8;
                        recDmrs_imPtr += 8;

                        genDmrs_rePtr += 8;
                        genDmrs_imPtr += 8;

                        lsEst_rePtr += 8;
                        lsEst_imPtr += 8;

                    }   //for idxSc
                }
                else if(pNrPuschInParams->nPortIndex[idxLyr] == 2)
                {
                    float *recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC + startScIdx];
                    float *recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];

                    float *genDmrs_rePtr = &pNrPuschOutParams->genDMRS[startScIdx];
                    float *genDmrs_imPtr = &pNrPuschOutParams->genDMRS[(MAX_NUM_SC/2) + startScIdx];

                    float *lsEst_rePtr = &pNrPuschOutParams->lsEst[idxAnt][idxLyr][startScIdx];
                    float *lsEst_imPtr = &pNrPuschOutParams->lsEst[idxAnt][idxLyr][(MAX_NUM_SC/2) + startScIdx];

                    for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                    {
                        if((idxSc + 8) >= endScIdx)
                        {
                            remSc = endScIdx - idxSc;

                            rDmrs_re = _mm256_setzero_ps();
                            rDmrs_im = _mm256_setzero_ps();

                            gDmrs_re = _mm256_setzero_ps();
                            gDmrs_im = _mm256_setzero_ps();

                            for(i = 0; i < remSc; i++)
                            {
                                rDmrs_re_f[i] = recDmrs_rePtr[i];
                                rDmrs_im_f[i] = recDmrs_imPtr[i];

                                gDmrs_re_f[i] = genDmrs_rePtr[i];
                                gDmrs_im_f[i] = genDmrs_imPtr[i];
                            }
                        }
                        else
                        {
                            rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                            rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);

                            gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                            gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                        }

                        bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                        ad_im = _mm256_mul_ps(rDmrs_re, gDmrs_im);

                        out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                        out_im = _mm256_fmsub_ps(rDmrs_im, gDmrs_re, ad_im);

                        out_re = _mm256_mul_ps(out_re, div2Vec);
                        out_im = _mm256_mul_ps(out_im, div2Vec);

                        if((idxSc + 8) >= (endScIdx))
                        {
                            remSc = (endScIdx) - idxSc;

                            for(i = 0; i < remSc; i++)
                            {
                                lsEst_rePtr[i] = out_re[i];
                                lsEst_imPtr[i] = out_im[i];
                            }

                        }
                        else
                        {
                            _mm256_storeu_ps(lsEst_rePtr, out_re);
                            _mm256_storeu_ps(lsEst_imPtr, out_im);
                        }
                        recDmrs_rePtr += 8;
                        recDmrs_imPtr += 8;

                        genDmrs_rePtr += 8;
                        genDmrs_imPtr += 8;

                        lsEst_rePtr += 8;
                        lsEst_imPtr += 8;

                    }   //for idxSc
                }
                else if(pNrPuschInParams->nPortIndex[idxLyr] == 3)
                {
                    float *recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC + startScIdx];
                    float *recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];

                    float *genDmrs_rePtr = &pNrPuschOutParams->genDMRS[startScIdx];
                    float *genDmrs_imPtr = &pNrPuschOutParams->genDMRS[(MAX_NUM_SC/2) + startScIdx];

                    float *lsEst_rePtr = &pNrPuschOutParams->lsEst[idxAnt][idxLyr][startScIdx];
                    float *lsEst_imPtr = &pNrPuschOutParams->lsEst[idxAnt][idxLyr][(MAX_NUM_SC/2) + startScIdx];

                    for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                    {
                        if((idxSc + 8) >= endScIdx)
                        {
                            remSc = endScIdx - idxSc;

                            rDmrs_re = _mm256_setzero_ps();
                            rDmrs_im = _mm256_setzero_ps();

                            gDmrs_re = _mm256_setzero_ps();
                            gDmrs_im = _mm256_setzero_ps();

                            for(i = 0; i < remSc; i++)
                            {
                                rDmrs_re_f[i] = recDmrs_rePtr[i];
                                rDmrs_im_f[i] = recDmrs_imPtr[i];

                                gDmrs_re_f[i] = alternateSigns_f[i]*genDmrs_rePtr[i];
                                gDmrs_im_f[i] = alternateSigns_f[i]*genDmrs_imPtr[i];
                            }
                        }
                        else
                        {
                            rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                            rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);

                            // gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                            // gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);

                            gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                            gDmrs_re = _mm256_mul_ps(gDmrs_re, alternateSigns);
                            gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                            gDmrs_im = _mm256_mul_ps(gDmrs_im, alternateSigns);
                        }

                        // bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                        // ad_im = _mm256_mul_ps(rDmrs_im, gDmrs_re);

                        // out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                        // out_im = _mm256_fmsub_ps(rDmrs_re, gDmrs_im, ad_im);

                        bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                        ad_im = _mm256_mul_ps(rDmrs_re, gDmrs_im);

                        out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                        out_im = _mm256_fmsub_ps(rDmrs_im, gDmrs_re, ad_im);

                        out_re = _mm256_mul_ps(out_re, div2Vec);
                        out_im = _mm256_mul_ps(out_im, div2Vec);

                        if((idxSc + 8) >= (endScIdx))
                        {
                            remSc = (endScIdx) - idxSc;

                            for(i = 0; i < remSc; i++)
                            {
                                lsEst_rePtr[i] = out_re[i];
                                lsEst_imPtr[i] = out_im[i];
                            }

                        }
                        else
                        {
                            _mm256_storeu_ps(lsEst_rePtr, out_re);
                            _mm256_storeu_ps(lsEst_imPtr, out_im);
                        }

                        recDmrs_rePtr += 8;
                        recDmrs_imPtr += 8;

                        genDmrs_rePtr += 8;
                        genDmrs_imPtr += 8;

                        lsEst_rePtr += 8;
                        lsEst_imPtr += 8;

                    }   //for idxSc
                }
            }
        }
    }

    #ifdef DEBUG
    // int idxAnt;
    nRxAnt   = pNrPuschInParams->nAntennaPorts;
    // int idxLyr;
    nLyr     = pNrPuschInParams->nNrOfLayers;
    FILE* fp;
    if(symPart == 0)
    {
        // printf("Real: %f\n", pNrPuschOutParams->lsEst[0][2][0]);
        fp = fopen("IOs/modulationRemovalLow_mmse_HPrb.txt", "w");
        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
            {
                for(int i11 = 0;i11<275*6*2;i11++)
                {
                    fprintf(fp, "%f\n", pNrPuschOutParams->lsEst[idxAnt][idxLyr][i11]);
                }
            }
        }
        fclose(fp);
    }
    else if(symPart == 1)
    {
        fp = fopen("IOs/modulationRemovalHigh_mmse_HPrb.txt", "w");
        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
            {
                for(int i11 = 0;i11<275*6*2;i11++)
                {
                    fprintf(fp, "%f\n", pNrPuschOutParams->lsEst[idxAnt][idxLyr][i11]);
                }
            }
        }
        fclose(fp);
    }
    #endif
}

void nr_pusch_mmse_perHPrb_toe_avx2( int symPart,
                        commonUlConfigT* pNrUlCommonParams,
                        puschConfigT* pNrPuschInParams,
                        P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams)
{
    int idxAnt, idxLyr, idxAlloc, idxSc, idxHPrb, i;
    int startScIdx, endScIdx, remSc;
    int startHPrbIdx, endHPrbIdx, remHPrb;

    int dmrsSym, nRxAnt, nLyr, nSbAlloc, nHPrb_alloc, nHPrb_total;

    __m256 h_t1_re, h_t1_im;
    __m256 h_t2_re, h_t2_im;
    __m256 h_t3_re, h_t3_im;

    __m256 B1_re, B1_im;
    __m256 B2_re, B2_im;
    __m256 B3_re, B3_im;

    __m256 C_re, C_im;
    __m256 toe_re, toe_im;
    __m256 div_2;

    __m256 mid_re, mid_im;

    float* h_t1_re_f = (float*)&h_t1_re;
    float* h_t1_im_f = (float*)&h_t1_im;

    float* h_t2_re_f = (float*)&h_t2_re;
    float* h_t2_im_f = (float*)&h_t2_im;

    float* h_t3_re_f = (float*)&h_t3_re;
    float* h_t3_im_f = (float*)&h_t3_im;

    float* toe_re_f = (float*)&toe_re;
    float* toe_im_f = (float*)&toe_im;

    nRxAnt   = pNrPuschInParams->nAntennaPorts;
    nLyr     = pNrPuschInParams->nNrOfLayers;
    nSbAlloc = pNrPuschInParams->nSbAlloc;
    nHPrb_total = pNrPuschOutParams->total_nPRB * 2;

    toe_re = _mm256_setzero_ps();
    toe_im = _mm256_setzero_ps();

    div_2 = _mm256_set1_ps(0.5);

    int testCnt  = 0;

    for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
    {
        for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
        {
            for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
            {
                nHPrb_alloc = pNrPuschInParams->nPrb[idxAlloc] * 2;

                startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                startHPrbIdx = pNrPuschInParams->rbStart[idxAlloc];
                endHPrbIdx   = startHPrbIdx + (pNrPuschInParams->nPrb[idxAlloc] * 2);

                float *lsEst_rePtr = &pNrPuschOutParams->lsEst[idxAnt][idxLyr][startScIdx];
                float *lsEst_imPtr = &pNrPuschOutParams->lsEst[idxAnt][idxLyr][(MAX_NUM_SC/2) + startScIdx];

                for(idxHPrb = 0; idxHPrb < (nHPrb_alloc); idxHPrb += 8)
                {
                    if((idxHPrb + 8) >= nHPrb_alloc)
                    {
                        remHPrb = endHPrbIdx - idxHPrb;

                        h_t1_re = _mm256_setzero_ps();
                        h_t1_im = _mm256_setzero_ps();

                        h_t2_re = _mm256_setzero_ps();
                        h_t2_im = _mm256_setzero_ps();

                        h_t3_re = _mm256_setzero_ps();
                        h_t3_im = _mm256_setzero_ps();
                    }
                    else
                    {
                        remHPrb = 8;
                    }

                    for(i = 0; i < remHPrb; i++)
                    {
                        h_t1_re_f[i] = lsEst_rePtr[3*i];
                        h_t1_im_f[i] = lsEst_imPtr[3*i];

                        h_t2_re_f[i] = lsEst_rePtr[3*i+1];
                        h_t2_im_f[i] = lsEst_imPtr[3*i+1];

                        h_t3_re_f[i] = lsEst_rePtr[3*i+2];
                        h_t3_im_f[i] = lsEst_imPtr[3*i+2];
                    }

                    mid_re = _mm256_mul_ps(h_t1_im, h_t2_im);
                    mid_im = _mm256_mul_ps(h_t1_im, h_t2_re);


                    B1_re = _mm256_fmadd_ps(h_t1_re, h_t2_re, mid_re);
                    B1_im = _mm256_fmsub_ps(h_t1_re, h_t2_im, mid_im);

                    mid_re = _mm256_mul_ps(h_t2_im, h_t3_im);
                    mid_im = _mm256_mul_ps(h_t2_im, h_t3_re);

                    B2_re = _mm256_fmadd_ps(h_t2_re, h_t3_re, mid_re);
                    B2_im = _mm256_fmsub_ps(h_t2_re, h_t3_im, mid_im);

                    B3_re = _mm256_add_ps(B1_re, B2_re);
                    B3_im = _mm256_add_ps(B1_im, B2_im);

                    B3_re = _mm256_mul_ps(B3_re, div_2);
                    B3_im = _mm256_mul_ps(B3_im, div_2);


                    mid_im = _mm256_mul_ps(B3_im, B3_im);
                    mid_re = _mm256_fmadd_ps(B3_re, B3_re, mid_im);

                    mid_re = _mm256_sqrt_ps(mid_re);


                    if(remHPrb < 8)
                    {
                        for(i = remHPrb; i < 8; i++)
                        {
                            mid_re[i] = 1;
                        }

                    }

                    C_re = _mm256_div_ps(B3_re, mid_re);
                    C_im = _mm256_div_ps(B3_im, mid_re);

                    toe_re = _mm256_add_ps(toe_re, C_re);
                    toe_im = _mm256_add_ps(toe_im, C_im);

                    lsEst_rePtr += 24;
                    lsEst_imPtr += 24;

                }
            }
        }
    }

    mid_re = _mm256_permute2f128_ps (toe_re,toe_re,1);
    toe_re = _mm256_add_ps(toe_re, mid_re);
    toe_re = _mm256_hadd_ps(toe_re,toe_re);
    toe_re = _mm256_hadd_ps(toe_re,toe_re);

    mid_im = _mm256_permute2f128_ps (toe_im,toe_im,1);
    toe_im = _mm256_add_ps(toe_im, mid_im);
    toe_im = _mm256_hadd_ps(toe_im,toe_im);
    toe_im = _mm256_hadd_ps(toe_im,toe_im);
    
    pNrPuschOutParams->toc_overAll[0] = toe_re_f[0]  /  (nHPrb_total * nRxAnt * nLyr);
    pNrPuschOutParams->toc_overAll[1] = toe_im_f[0]  /  (nHPrb_total * nRxAnt * nLyr);

    // printf("TO : %f\n", pNrPuschOutParams->toc_overAll[0]);
    // printf("TO : %f\n", pNrPuschOutParams->toc_overAll[1]);
}


void nr_pusch_mmse_perHPrb_est_avx2( int symPart,
                                commonUlConfigT* pNrUlCommonParams,
                                puschConfigT* pNrPuschInParams,
                                P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams)
{
    nr_pusch_mmse_perHPrb_modRemoval_avx2(symPart, pNrUlCommonParams, pNrPuschInParams, pNrPuschOutParams);
    nr_pusch_mmse_perHPrb_toe_avx2(symPart, pNrUlCommonParams, pNrPuschInParams, pNrPuschOutParams);

    int idxAnt, idxLyr, idxAlloc, idxSc, idxHPrb, i;
    int startScIdx, endScIdx, remSc;
    int startHPrbIdx, endHPrbIdx, remHPrb;

    int dmrsSym, nRxAnt, nLyr, nSbAlloc, nHPrb_alloc, nHPrb_total;

    float nVar_overAll = 0;

    __m256 mid_re, mid_im;

    __m256 h_t1_re, h_t1_im;
    __m256 h_t2_re, h_t2_im;
    __m256 h_t3_re, h_t3_im;
    __m256 h_avg_re, h_avg_im;

    __m256 B1_re, B1_im;
    __m256 B2_re, B2_im;
    __m256 B3_re, B3_im;

    __m256 toc_re, toc_im;
    __m256 toc_sq_re, toc_sq_im;
    __m256 sq_con_toc_re, sq_con_toc_im;
    __m256 div_3, scale_vec;

    __m256 n_t1_re, n_t1_im;
    __m256 n_t2_re, n_t2_im;
    __m256 n_t3_re, n_t3_im;
    __m256 n_avg;
    __m256 nVar_alloc_vec;

    float* nVar_alloc_vec_f = (float*)&nVar_alloc_vec;

    float* h_t1_re_f = (float*)&h_t1_re;
    float* h_t1_im_f = (float*)&h_t1_im;

    float* h_t2_re_f = (float*)&h_t2_re;
    float* h_t2_im_f = (float*)&h_t2_im;

    float* h_t3_re_f = (float*)&h_t3_re;
    float* h_t3_im_f = (float*)&h_t3_im;

    nRxAnt   = pNrPuschInParams->nAntennaPorts;
    nLyr     = pNrPuschInParams->nNrOfLayers;
    nSbAlloc = pNrPuschInParams->nSbAlloc;

    toc_re = _mm256_set1_ps(pNrPuschOutParams->toc_overAll[0]);
    toc_im = _mm256_set1_ps(pNrPuschOutParams->toc_overAll[1]);

    mid_re = _mm256_mul_ps(toc_im,toc_im);
    mid_im = _mm256_mul_ps(toc_im,toc_re);

    toc_sq_re = _mm256_fmsub_ps(toc_re,toc_re,mid_re);
    toc_sq_im = _mm256_fmadd_ps(toc_re,toc_im,mid_im);

    float complex con_toc = pNrPuschOutParams->toc_overAll[0] - I * pNrPuschOutParams->toc_overAll[1];

    float complex sq_con_toc = csqrtf(con_toc);

    sq_con_toc_re = _mm256_set1_ps(crealf(sq_con_toc));
    sq_con_toc_im = _mm256_set1_ps(cimagf(sq_con_toc));

    div_3 = _mm256_set1_ps(0.333333);

    for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
    {
        for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
        {
            nVar_alloc_vec = _mm256_setzero_ps();

            for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
            {
                nHPrb_alloc = pNrPuschInParams->nPrb[idxAlloc] * 2;

//                scale_vec = _mm256_set1_ps(1.0 / (3.0 * nHPrb_alloc));
                scale_vec = _mm256_set1_ps(0.333333);

                startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                startHPrbIdx = pNrPuschInParams->rbStart[idxAlloc];
                endHPrbIdx   = startHPrbIdx + (pNrPuschInParams->nPrb[idxAlloc] * 2);

                float *lsEst_rePtr = &pNrPuschOutParams->lsEst[idxAnt][idxLyr][startScIdx];
                float *lsEst_imPtr = &pNrPuschOutParams->lsEst[idxAnt][idxLyr][(MAX_NUM_SC/2) + startScIdx];

                float *chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][2*startHPrbIdx];
                float *chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][2*startHPrbIdx];

                for(idxHPrb = 0; idxHPrb < nHPrb_alloc; idxHPrb += 8)
                {
                    if((idxHPrb + 8) >= nHPrb_alloc)
                    {
                        remHPrb = endHPrbIdx - idxHPrb;

                        h_t1_re = _mm256_setzero_ps();
                        h_t1_im = _mm256_setzero_ps();

                        h_t2_re = _mm256_setzero_ps();
                        h_t2_im = _mm256_setzero_ps();

                        h_t3_re = _mm256_setzero_ps();
                        h_t3_im = _mm256_setzero_ps();
                    }
                    else
                    {
                        remHPrb = 8;
                    }

                    for(i = 0; i < remHPrb; i++)
                    {
                        h_t1_re_f[i] = lsEst_rePtr[3*i];
                        h_t1_im_f[i] = lsEst_imPtr[3*i];

                        h_t2_re_f[i] = lsEst_rePtr[3*i+1];
                        h_t2_im_f[i] = lsEst_imPtr[3*i+1];

                        h_t3_re_f[i] = lsEst_rePtr[3*i+2];
                        h_t3_im_f[i] = lsEst_imPtr[3*i+2];
                    }

                    mid_re = _mm256_mul_ps(h_t2_im, toc_im);
                    mid_im = _mm256_mul_ps(h_t2_re, toc_im);

                    h_t2_re = _mm256_fmadd_ps(h_t2_re, toc_re, mid_re);
                    h_t2_im = _mm256_fmsub_ps(h_t2_im, toc_re, mid_im);

                    mid_re = _mm256_mul_ps(h_t3_im, toc_sq_im);
                    mid_im = _mm256_mul_ps(h_t3_re, toc_sq_im);

                    h_t3_re = _mm256_fmadd_ps(h_t3_re, toc_sq_re, mid_re);
                    h_t3_im = _mm256_fmsub_ps(h_t3_im, toc_sq_re, mid_im);

                    h_avg_re = _mm256_add_ps(h_t1_re, h_t2_re);
                    h_avg_im = _mm256_add_ps(h_t1_im, h_t2_im);

                    h_avg_re = _mm256_add_ps(h_avg_re, h_t3_re);
                    h_avg_im = _mm256_add_ps(h_avg_im, h_t3_im);

                    if((pNrPuschInParams->nPortIndex[idxLyr] == 2) || (pNrPuschInParams->nPortIndex[idxLyr] == 3))
                    {
                        mid_re = _mm256_mul_ps(h_avg_im,sq_con_toc_im);
                        mid_im = _mm256_mul_ps(h_avg_re,sq_con_toc_im);

                        h_avg_re = _mm256_fmsub_ps(h_avg_re,sq_con_toc_re,mid_re);
                        h_avg_im = _mm256_fmadd_ps(h_avg_im,sq_con_toc_re,mid_im);
                    }

                    h_avg_re = _mm256_mul_ps(h_avg_re, div_3);
                    h_avg_im = _mm256_mul_ps(h_avg_im, div_3);

                    _mm256_storeu_ps(chEst_rePtr, h_avg_re);
                    _mm256_storeu_ps(chEst_imPtr, h_avg_im);

                    n_t1_re = _mm256_sub_ps(h_t1_re, h_avg_re);
                    n_t1_im = _mm256_sub_ps(h_t1_im, h_avg_im);

                    n_t2_re = _mm256_sub_ps(h_t2_re, h_avg_re);
                    n_t2_im = _mm256_sub_ps(h_t2_im, h_avg_im);

                    n_t3_re = _mm256_sub_ps(h_t3_re, h_avg_re);
                    n_t3_im = _mm256_sub_ps(h_t3_im, h_avg_im);

                    n_avg = _mm256_setzero_ps();
                    n_avg = _mm256_fmadd_ps(n_t1_re, n_t1_re, n_avg);
                    n_avg = _mm256_fmadd_ps(n_t1_im, n_t1_im, n_avg);
                    n_avg = _mm256_fmadd_ps(n_t2_re, n_t2_re, n_avg);
                    n_avg = _mm256_fmadd_ps(n_t2_im, n_t2_im, n_avg);
                    n_avg = _mm256_fmadd_ps(n_t3_re, n_t3_re, n_avg);
                    n_avg = _mm256_fmadd_ps(n_t3_im, n_t3_im, n_avg);
                    n_avg = _mm256_mul_ps(n_avg, scale_vec);

                    nVar_alloc_vec = _mm256_add_ps(nVar_alloc_vec, n_avg);

                    lsEst_rePtr += 24;
                    lsEst_imPtr += 24;

                    chEst_rePtr += 8;
                    chEst_imPtr += 8;
                }

                if(symPart == 0)
                {
                    // printf("nVar_alloc_vec_f: %f %f %f %f %f %f %f %f\n", nVar_alloc_vec_f[0], nVar_alloc_vec_f[1],\
                    //         nVar_alloc_vec_f[2], nVar_alloc_vec_f[3],\
                    //         nVar_alloc_vec_f[4], nVar_alloc_vec_f[5],\
                    //         nVar_alloc_vec_f[6], nVar_alloc_vec_f[7]);
                    // printf("\n");
                }
                
                mid_re = _mm256_permute2f128_ps (nVar_alloc_vec, nVar_alloc_vec, 1);
                nVar_alloc_vec = _mm256_add_ps(nVar_alloc_vec, mid_re);
                nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);
                nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);

                pNrPuschOutParams->nVar[idxAnt][idxLyr][idxAlloc] = nVar_alloc_vec_f[0] ;

                nVar_overAll += nVar_alloc_vec_f[0];
            }
        }
    }

    
    pNrPuschOutParams->nVar_overAll[0] = nVar_overAll / (nLyr * nHPrb_alloc * nRxAnt * nSbAlloc);
    // pNrPuschOutParams->nVar_overAll[0] = 1.3979;
    // printf("nVar_overAll: %f\n", pNrPuschOutParams->nVar_overAll[0]);
}


void nr_pusch_mmse_perHPrb_equ_avx2( int symPart,
                                    commonUlConfigT* pNrUlCommonParams,
                                    puschConfigT* pNrPuschInParams,
                                    P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams)
{

    int idxAnt, idxLyr, idxAlloc, idxHPrb, idxSym, idxSc, rIndx, cIndx, i, idx;
    int startScIdx, endScIdx, remSc;
    int startHPrbIdx, endHPrbIdx, remHPrb, rbStartIdx;

    int startSym, endSym, dmrsSym;

    float zeroBuf[6] = {0, 0, 0, 0, 0, 0};

    int hPrbScStartIdx_re[8], hPrbScStartIdx_im[8];

    int nRxAnt, nLyr, nSbAlloc, nHPrb_alloc, nHPrb_total, completedHPrb;

    char prbTemp;

    g_negate = VECT_SET1(-0.0);

    nRxAnt   = pNrPuschInParams->nAntennaPorts;
    nLyr     = pNrPuschInParams->nNrOfLayers;
    nSbAlloc = pNrPuschInParams->nSbAlloc;

    VECT_T H_re[nLyr][nRxAnt];
	VECT_T H_im[nLyr][nRxAnt];
	VECT_T HHh_re[nLyr][nLyr];
	VECT_T HHh_im[nLyr][nLyr];
	VECT_T invHHh_re[nLyr][nLyr];
	VECT_T invHHh_im[nLyr][nLyr];
	VECT_T weights_re[nRxAnt][nLyr];
	VECT_T weights_im[nRxAnt][nLyr];
	VECT_T Y_vect_re[1][nRxAnt];
	VECT_T Y_vect_im[1][nRxAnt];
	VECT_T Y_tilda_re[1][nLyr];
	VECT_T Y_tilda_im[1][nLyr];
	VECT_T Y_cap_re;
	VECT_T Y_cap_im;
	VECT_T conj_sq_toc_re[6];
	VECT_T conj_sq_toc_im[6];

	VECT_T nVar_vec = VECT_SET1(pNrPuschOutParams->nVar_overAll[0]);
//printf ("[Noise Variance : %f]\n",pNrPuschOutParams->nVar_overAll);
    if(symPart == LOW_PUSCH_SYM)
    {
        dmrsSym = pNrPuschInParams->dmrs1_sym;
        startSym = 0;
        endSym = 6;
    }
    else
    {
        dmrsSym = pNrPuschInParams->dmrs2_sym;
        startSym = 7;
        endSym = 13;
    }

    float nVar_overAll        = pNrPuschOutParams->nVar_overAll[0];
    float complex toc_overAll = pNrPuschOutParams->toc_overAll[0] + I * pNrPuschOutParams->toc_overAll[1];

    float complex sqrt_toc_1, sqrt_toc_2, sqrt_toc_3, sqrt_toc_4, sqrt_toc_5;

    sqrt_toc_1 = csqrtf( toc_overAll );
    sqrt_toc_2 = sqrt_toc_1 * sqrt_toc_1;
    sqrt_toc_3 = sqrt_toc_2 * sqrt_toc_1;
    sqrt_toc_4 = sqrt_toc_3 * sqrt_toc_1;
    sqrt_toc_5 = sqrt_toc_4 * sqrt_toc_1;

    conj_sq_toc_re[0] = VECT_SET1(1);
    conj_sq_toc_im[0] = VECT_SET1(0);

    conj_sq_toc_re[1] = VECT_SET1(creal(conj(sqrt_toc_1)));
    conj_sq_toc_im[1] = VECT_SET1(cimag(conj(sqrt_toc_1)));

    conj_sq_toc_re[2] = VECT_SET1(creal(conj(sqrt_toc_2)));
    conj_sq_toc_im[2] = VECT_SET1(cimag(conj(sqrt_toc_2)));

    conj_sq_toc_re[3] = VECT_SET1(creal(conj(sqrt_toc_3)));
    conj_sq_toc_im[3] = VECT_SET1(cimag(conj(sqrt_toc_3)));

    conj_sq_toc_re[4] = VECT_SET1(creal(conj(sqrt_toc_4)));
    conj_sq_toc_im[4] = VECT_SET1(cimag(conj(sqrt_toc_4)));

    conj_sq_toc_re[5] = VECT_SET1(creal(conj(sqrt_toc_5)));
    conj_sq_toc_im[5] = VECT_SET1(cimag(conj(sqrt_toc_5)));

    nHPrb_total = 0;

   float inv_nVar_overAll;
   int testCnt = 0;
   float sum = 0;
   int total_HPrb = 0;
   float one_By_total_HPrb = 0;

   for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
   {
      total_HPrb = total_HPrb + pNrPuschInParams->nPrb[idxAlloc] * 2;
   }
   one_By_total_HPrb = 1.0/total_HPrb;

   for(rIndx = 0; rIndx < nLyr; rIndx++)
   {
      pNrPuschOutParams->oneBymse_avg_re[rIndx] = 0;
   }

    wnFlt SNRPerHPrb;
    wnUInt16 mseCnt = 0;

    for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
    {
        nHPrb_alloc = pNrPuschInParams->nPrb[idxAlloc] * 2;
        completedHPrb = nHPrb_total;
        nHPrb_total += nHPrb_alloc;

        rbStartIdx = pNrPuschInParams->rbStart[idxAlloc];
        startHPrbIdx = rbStartIdx * 2;
        endHPrbIdx   = startHPrbIdx + (pNrPuschInParams->nPrb[idxAlloc] * 2);

        float *chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][startHPrbIdx];
        float *chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][startHPrbIdx];

        for(idxHPrb = startHPrbIdx; idxHPrb < (endHPrbIdx); idxHPrb += 8)
        {
            for(idx = 0; idx < 8; idx++)
            {
                hPrbScStartIdx_re[idx] = (idxHPrb + idx) * 6;
                hPrbScStartIdx_im[idx] = MAX_NUM_SC + (idxHPrb + idx) * 6;
            }

           for(rIndx = 0; rIndx < nLyr; rIndx++)
            {
                for(cIndx = 0; cIndx < nRxAnt; cIndx++)
                {
                    H_re[rIndx][cIndx] = _mm256_loadu_ps(&pNrPuschOutParams->chEst_re[cIndx][rIndx][idxHPrb]);
                    H_im[rIndx][cIndx] = _mm256_loadu_ps(&pNrPuschOutParams->chEst_im[cIndx][rIndx][idxHPrb]);

                }
            }
            

            if((idxHPrb + 8) >= endHPrbIdx)
            {
                remHPrb = endHPrbIdx - idxHPrb;

                for(rIndx = 0; rIndx < nLyr; rIndx++)
                {
                    for(cIndx = 0; cIndx < nRxAnt; cIndx++)
                    {
                        for(idx = remHPrb; idx < 8; idx++)
                        {
                            H_re[rIndx][cIndx][idx] = 1;
                            H_im[rIndx][cIndx][idx] = 1;
                        }
                    }
                }
            }
            else
            {
                remHPrb = 8;
            }

            // HH'
            MATRIXMULTTH(H_re, H_im, HHh_re, HHh_im, nLyr, nRxAnt, nLyr);

            // inv_nVar_overAll = 1.0f/nVar_overAll;

            // for(rIndx = 0; rIndx < nLyr; rIndx++)
            // {
            //     for(idx = 0; idx < 8; idx++)
            //     {
            //         // pNrPuschOutParams->mse_re[rIndx][idxHPrb + idx] = HHh_re[rIndx][rIndx][idx]*inv_nVar_overAll + 1;
            //     }
            // }

            for(idx = 0; idx < nLyr; idx++)
            {
                HHh_re[idx][idx] = VECT_ADD(nVar_vec, HHh_re[idx][idx]);
            }

            // Inverse of (HH' + (epsilon)I)
            if(nLyr == 1)
            {
                InverseMatrix1x1HTranspose(HHh_re, HHh_im, invHHh_re, invHHh_im);
            }
            else if(nLyr == 2)
            {
                InverseMatrix2x2H(HHh_re, HHh_im, invHHh_re, invHHh_im);
            }
            else if(nLyr == 4)
            {
                InverseMatrix4x4H(HHh_re, HHh_im, invHHh_re, invHHh_im);
            }

            
            for(idx = 0; idx < remHPrb; idx++)
            {
                for(rIndx = 0; rIndx < nLyr; rIndx++)
                {
                    SNRPerHPrb = 1/invHHh_re[rIndx][rIndx][idx];
                    for(wnUInt8 toneIdx = 0;toneIdx<6;toneIdx++)
                    {
                        pNrPuschOutParams->oneBymse_re[rIndx][mseCnt+toneIdx] = SNRPerHPrb;
                    }
                    sum = pNrPuschOutParams->oneBymse_avg_re[rIndx];
                    sum = sum + invHHh_re[rIndx][rIndx][idx];
                    pNrPuschOutParams->oneBymse_avg_re[rIndx] = sum;
                }
                mseCnt += 6;
            }
                     

            // for(rIndx = 0; rIndx < nLyr; rIndx++)
            // {
            //     sum = pNrPuschOutParams->oneBymse_avg_re[rIndx];

            //     for(idx = 0; idx < 8; idx++)
            //     {
            //         sum = sum + invHHh_re[rIndx][rIndx][idx]*nVar_overAll*one_By_total_HPrb;
            //     }
            //     pNrPuschOutParams->oneBymse_avg_re[rIndx] = sum;
            // }

            // Multiplication of H' and inverse of (HH' + (epsilon)I)
            MATRIXMULTHTTRANSPOSE2(H_re, H_im, invHHh_re, invHHh_im, weights_re, weights_im, nRxAnt, nLyr, nLyr);

            int equaSym = 0;
            for(idxSym = startSym; idxSym <= endSym; idxSym++)
            {
                if(idxSym != dmrsSym)
                {
                    int layerMapIdx,temp1 = (6*(startSym>0)+equaSym)*pNrPuschOutParams->total_nPRB*12*nLyr + completedHPrb*6*nLyr;
                    for(idxSc = 0; idxSc < 6; idxSc++)
                    {
                        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                        {
                            Y_vect_re[0][idxAnt] = VECT_SET(rxFdSamples[idxAnt][idxSym][hPrbScStartIdx_re[7] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][hPrbScStartIdx_re[6] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][hPrbScStartIdx_re[5] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][hPrbScStartIdx_re[4] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][hPrbScStartIdx_re[3] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][hPrbScStartIdx_re[2] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][hPrbScStartIdx_re[1] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][hPrbScStartIdx_re[0] + idxSc] );

                            Y_vect_im[0][idxAnt] = VECT_SET(rxFdSamples[idxAnt][idxSym][hPrbScStartIdx_im[7] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][hPrbScStartIdx_im[6] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][hPrbScStartIdx_im[5] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][hPrbScStartIdx_im[4] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][hPrbScStartIdx_im[3] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][hPrbScStartIdx_im[2] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][hPrbScStartIdx_im[1] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][hPrbScStartIdx_im[0] + idxSc] );
                        }

                        MATRIXMULNN(Y_vect_re, Y_vect_im, weights_re, weights_im, Y_tilda_re, Y_tilda_im, 1, nRxAnt, nLyr);

                        for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
                        {
                            C_VECT_MUL_NA_NB(Y_tilda_re[0][idxLyr], Y_tilda_im[0][idxLyr], conj_sq_toc_re[idxSc],\
                             conj_sq_toc_im[idxSc], Y_cap_re, Y_cap_im);

                            for(idx = 0; idx < remHPrb; idx++)
                            {
                                pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][hPrbScStartIdx_re[idx] + idxSc] = Y_cap_re[idx];
                                pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][hPrbScStartIdx_im[idx] + idxSc] = Y_cap_im[idx];
                                // Directly stores data in their respective Layer Demapped Positions
                                layerMapIdx = temp1 + (hPrbScStartIdx_re[idx]-startHPrbIdx*6+idxSc)*nLyr + idxLyr;
                                pNrPuschOutParams->layerDemapperOutReal[layerMapIdx] = Y_cap_re[idx];
                                pNrPuschOutParams->layerDemapperOutImag[layerMapIdx] = Y_cap_im[idx];
                            }

                        }
                    }
                    equaSym = equaSym + 1;
                }
            }

        }
    }
    totalLength = 12*pNrPuschOutParams->total_nPRB*12*nLyr;

    for(rIndx = 0; rIndx < nLyr; rIndx++)
    {
        pNrPuschOutParams->oneBymse_avg_re[rIndx] = \
        pNrPuschOutParams->oneBymse_avg_re[rIndx]*one_By_total_HPrb;
    }

    #ifdef DEBUG
        FILE *fp, *pf;
        if(symPart == 0)
        {
            fp = fopen("IOs/EqualiserLowReal_HPrb.txt", "w");
            pf = fopen("IOs/EqualiserLowImag_HPrb.txt", "w");
        }
        else
        {
            fp = fopen("IOs/EqualiserHighReal_HPrb.txt", "w");
            pf = fopen("IOs/EqualiserHighImag_HPrb.txt", "w");
        }

        for(int idx = 0;idx<6*pNrPuschOutParams->total_nPRB*12*nLyr;idx++)
        {
            fprintf(fp, "%f\n", pNrPuschOutParams->layerDemapperOutReal[idx]);
            fprintf(pf, "%f\n", pNrPuschOutParams->layerDemapperOutImag[idx]);
        }

        fclose(fp);
        fclose(pf);
    #endif
}

//} Per HPrb End */


// { //Per Tone Start
void nr_pusch_mmse_perTone_est_avx2(wnUInt8 segPart,
                               commonUlConfigT* pNrUlCommonParams,
                               puschConfigT* pNrPuschInParams,
                               P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams)
{
    wnUInt16 idxAnt, idxLyr, idxAlloc, idxSc, i;
    wnUInt16 startScIdx, endScIdx, remSc;
    wnUInt16 toneStartIdx, toneEndIdx;
    wnUInt8 nRxAnt, nLyr, nSbAlloc, dmrsSym;

    __m256 rDmrs_re, rDmrs_im;
    __m256 gDmrs_re, gDmrs_im;
    __m256 bd_re, ad_im;
    __m256 out_re, out_im;
    wnFlt *out_re_f = (wnFlt*)&out_re;
    wnFlt *out_im_f = (wnFlt*)&out_im;
    __m256 temp_re1, temp_im1;
    __m256 temp_re2, temp_im2;

    __m256 zeroVec = _mm256_setzero_ps();
    __m256 div2Vec = _mm256_set1_ps(0.5);
    __m256* AVX_chEst_rePtr, *AVX_chEst_imPtr;
    __m256 AVX_chEst_re, AVX_chEst_im;

    wnFlt* rDmrs_re_f = (wnFlt*)&rDmrs_re;
    wnFlt* rDmrs_im_f = (wnFlt*)&rDmrs_im;
    wnFlt* gDmrs_re_f = (wnFlt*)&gDmrs_re;
    wnFlt* gDmrs_im_f = (wnFlt*)&gDmrs_im;
    wnFlt *recDmrs_rePtr;
    wnFlt *recDmrs_imPtr;
    wnFlt *genDmrs_rePtr;
    wnFlt *genDmrs_imPtr;
    // wnFlt *lsEst_rePtr;
    // wnFlt *lsEst_imPtr;
    wnFlt *chEst_rePtr, *chEst_imPtr;
    wnUInt8 CDMFlag;
    wnUInt8 cdmTotal[2] = {0, 0};

    __m256 alternateSigns = _mm256_setr_ps(1, -1, 1, -1, 1, -1, 1, -1);
    wnFlt* alternateSigns_f = (wnFlt *)&alternateSigns;

    //Calculate number of layers in  a given CDM group.
    // cdmTotal[0] => number of layers in CDM group 0.
    // cdmTotal[1] => number of layers in CDM group 1.
    for(idxLyr = 0;idxLyr<pNrPuschInParams->nNrOfLayers;idxLyr++)
    {
        if(pNrPuschInParams->nPortIndex[idxLyr] == 0 || pNrPuschInParams->nPortIndex[idxLyr] == 1)
        {
            cdmTotal[0] = cdmTotal[0]+1;
        }
        else
        {
            cdmTotal[1] = cdmTotal[1]+1;
        }
    }

    __m256i const1 = _mm256_set_epi32(1, 1, 6, 1, 4, 1, 2, 1);
    __m256i const2 = _mm256_set_epi32(0, 6, 5, 4, 3, 2, 1, 2);

    // Get DM-RS Index
    if(segPart == SEG1)
    {
        dmrsSym = pNrPuschInParams->dmrs1_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG2)
    {
        dmrsSym = pNrPuschInParams->dmrs2_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG3)
    {
        dmrsSym = pNrPuschInParams->dmrs3_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG4)
    {
        dmrsSym = pNrPuschInParams->dmrs4_sym - pNrPuschInParams->nStartSymbolIndex;
    }

    wnFlt cdm0EstRe[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_NUM_ALLOC*6];
    wnFlt cdm0EstIm[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_NUM_ALLOC*6];
    wnFlt *cdm0Est_rePtr, *cdm0Est_imPtr;
    __m256 cdm0EstRe1, cdm0EstRe2, cdm0EstIm1, cdm0EstIm2;

    wnFlt cdm1EstRe[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_NUM_ALLOC*6];
    wnFlt cdm1EstIm[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_NUM_ALLOC*6];
    wnFlt *cdm1Est_rePtr, *cdm1Est_imPtr;
    __m256 cdm1EstRe1, cdm1EstRe2, cdm1EstIm1, cdm1EstIm2;


    nRxAnt   = pNrPuschInParams->nAntennaPorts;
    nLyr     = pNrPuschInParams->nNrOfLayers;
    nSbAlloc = pNrPuschInParams->nSbAlloc;

    // If there exist CDM based transmission using CDM Group 0.
    if(cdmTotal[0] == 2)
    {
        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
            {
                startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
                toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);

                // Received DM-RS Signal                          
                recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][(MAX_NUM_SC/2) + startScIdx];
                
                // Generated DM-RS Signal
                genDmrs_rePtr = &pNrPuschOutParams->genDMRS[0][startScIdx];
                genDmrs_imPtr = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];

                cdm0Est_rePtr = &cdm0EstRe[idxAnt][startScIdx];
                cdm0Est_imPtr = &cdm0EstIm[idxAnt][startScIdx];

                // Modulation Removal
                for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                {
                    remSc = endScIdx - idxSc;
                    if( remSc < 8)
                    {
                        rDmrs_re = _mm256_setzero_ps();
                        rDmrs_im = _mm256_setzero_ps();

                        gDmrs_re = _mm256_setzero_ps();
                        gDmrs_im = _mm256_setzero_ps();

                        for(i = 0; i < remSc; i++)
                        {
                            rDmrs_re_f[i] = recDmrs_rePtr[i];
                            rDmrs_im_f[i] = recDmrs_imPtr[i];

                            gDmrs_re_f[i] = genDmrs_rePtr[i];
                            gDmrs_im_f[i] = genDmrs_imPtr[i];
                        }
                        bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                        ad_im = _mm256_mul_ps(rDmrs_re, gDmrs_im);

                        out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                        out_im = _mm256_fmsub_ps(rDmrs_im, gDmrs_re, ad_im);

                        // out = ( Y*conj(Xdmrs) )/Beta_dmrs;
                        out_re = _mm256_mul_ps(out_re, div2Vec);
                        out_im = _mm256_mul_ps(out_im, div2Vec);

                        for(i = 0;i<remSc;i++)
                        {
                            cdm0Est_rePtr[i] = out_re_f[i];
                            cdm0Est_imPtr[i] = out_im_f[i];
                        }
                    }
                    else
                    {
                        rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                        rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);

                        gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                        gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);

                        bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                        ad_im = _mm256_mul_ps(rDmrs_re, gDmrs_im);

                        out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                        out_im = _mm256_fmsub_ps(rDmrs_im, gDmrs_re, ad_im);

                        // out = ( Y*conj(Xdmrs) )/Beta_dmrs;
                        out_re = _mm256_mul_ps(out_re, div2Vec);
                        out_im = _mm256_mul_ps(out_im, div2Vec);

                        _mm256_storeu_ps(cdm0Est_rePtr, out_re);
                        _mm256_storeu_ps(cdm0Est_imPtr, out_im);
                    }

                    recDmrs_rePtr += 8;
                    recDmrs_imPtr += 8;

                    genDmrs_rePtr += 8;
                    genDmrs_imPtr += 8;

                    cdm0Est_rePtr += 8;
                    cdm0Est_imPtr += 8;
                }//for idxSc.
            }
        }

        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
            {
                for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                {
                    startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                    endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                    toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
                    toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);

                    cdm0Est_rePtr = &cdm0EstRe[idxAnt][startScIdx];
                    cdm0Est_imPtr = &cdm0EstIm[idxAnt][startScIdx];

                    // Channel Estimates after Linear Interpolation
                    chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                    chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];

                    // port 0 processing
                    if(pNrPuschInParams->nPortIndex[idxLyr] == 0)
                    {
                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc <= 8)
                            {
                                for(i = 0; i < remSc-1; i++)
                                {
                                    chEst_rePtr[2*i] = (cdm0Est_rePtr[i] + cdm0Est_rePtr[i+1])*0.5;
                                    chEst_rePtr[2*i+1] = 0;
                                    chEst_imPtr[2*i] = (cdm0Est_imPtr[i] + cdm0Est_imPtr[i+1])*0.5;
                                    chEst_imPtr[2*i+1] = 0;
                                }
                                chEst_rePtr[2*(remSc-1)] = chEst_rePtr[2*(remSc-2)];
                                chEst_rePtr[2*(remSc-1)+1] = 0;
                                chEst_imPtr[2*(remSc-1)] = chEst_imPtr[2*(remSc-2)];
                                chEst_imPtr[2*(remSc-1)+1] = 0;
                            }
                            else
                            {
                                cdm0EstRe1 = _mm256_loadu_ps(cdm0Est_rePtr);
                                cdm0EstIm1 = _mm256_loadu_ps(cdm0Est_imPtr);
                                cdm0EstRe2 = _mm256_loadu_ps(cdm0Est_rePtr+1);
                                cdm0EstIm2 = _mm256_loadu_ps(cdm0Est_imPtr+1);

                                out_re = _mm256_add_ps(cdm0EstRe1, cdm0EstRe2);
                                out_im = _mm256_add_ps(cdm0EstIm1, cdm0EstIm2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);

                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);

                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);

                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);

                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );

                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );

                            }
                            cdm0Est_rePtr += 8;
                            cdm0Est_imPtr += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }

                        //Frequency Interpolation 
                        //Calculate estimates for null tones using Left and Right estimates.
                        //Inp : [ a 0 b 0 c ...]
                        //Out : [a (a+b)/2 b (b+c)/2 c...]
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];
                        AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                        AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                        for(idxSc = toneStartIdx; idxSc < toneEndIdx; idxSc += 8)
                        {
                            if(toneEndIdx - idxSc<=8)
                            {
                                remSc = toneEndIdx - idxSc;

                                for(i = 0; i < remSc; i += 2)
                                {
                                    if(remSc - i == 2)
                                    {
                                        chEst_rePtr[i+1] = (chEst_rePtr[i]+chEst_rePtr[i-1])*0.5;
                                        chEst_imPtr[i+1] = (chEst_imPtr[i]+chEst_imPtr[i-1])*0.5;
                                    }
                                    else
                                    {
                                        chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                        chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                    }
                                }
                            }
                            else
                            {
                                AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                chEst_rePtr += 8;
                                chEst_imPtr += 8;
                            }
                        }
                    }
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 1)
                    {
                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc <= 8)
                            {
                                for(i = 0; i < remSc-1; i++)
                                {
                                    chEst_rePtr[2*i] = alternateSigns_f[i]*(cdm0Est_rePtr[i] - cdm0Est_rePtr[i+1])*0.5;
                                    chEst_rePtr[2*i+1] = 0;
                                    chEst_imPtr[2*i] = alternateSigns_f[i]*(cdm0Est_imPtr[i] - cdm0Est_imPtr[i+1])*0.5;
                                    chEst_imPtr[2*i+1] = 0;
                                }
                                chEst_rePtr[2*(remSc-1)] = chEst_rePtr[2*(remSc-2)];
                                chEst_rePtr[2*(remSc-1)+1] = 0;
                                chEst_imPtr[2*(remSc-1)] = chEst_imPtr[2*(remSc-2)];
                                chEst_imPtr[2*(remSc-1)+1] = 0;
                            }
                            else
                            {
                                cdm0EstRe1 = _mm256_loadu_ps(cdm0Est_rePtr);
                                cdm0EstIm1 = _mm256_loadu_ps(cdm0Est_imPtr);
                                cdm0EstRe2 = _mm256_loadu_ps(cdm0Est_rePtr+1);
                                cdm0EstIm2 = _mm256_loadu_ps(cdm0Est_imPtr+1);

                                out_re = _mm256_sub_ps(cdm0EstRe1, cdm0EstRe2);
                                out_re = _mm256_mul_ps(out_re, alternateSigns);
                                out_im = _mm256_sub_ps(cdm0EstIm1, cdm0EstIm2);
                                out_im = _mm256_mul_ps(out_im, alternateSigns);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);

                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);

                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);

                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);

                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );

                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                
                            }
                            cdm0Est_rePtr += 8;
                            cdm0Est_imPtr += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }

                        //Frequency Interpolation 
                        //Calculate estimates for null tones using Left and Right estimates.
                        //Inp : [ a 0 b 0 c ...]
                        //Out : [a (a+b)/2 b (b+c)/2 c...]
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];
                        AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                        AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                        for(idxSc = toneStartIdx; idxSc < toneEndIdx; idxSc += 8)
                        {
                            remSc = toneEndIdx - idxSc;
                            if( remSc<=8 )
                            {
                                for(i = 0; i < remSc; i += 2)
                                {
                                    if(remSc - i == 2)
                                    {
                                        chEst_rePtr[i+1] = (chEst_rePtr[i]+chEst_rePtr[i-1])*0.5;
                                        chEst_imPtr[i+1] = (chEst_imPtr[i]+chEst_imPtr[i-1])*0.5;
                                    }
                                    else
                                    {
                                        chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                        chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                    }
                                }
                            }
                            else
                            {
                                AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                chEst_rePtr += 8;
                                chEst_imPtr += 8;
                            }
                        }
                    }
                }
            }
        }
    }

    // If there exist CDM based transmission using CDM Group 1.
    if(cdmTotal[1] == 2)
    {
        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
            {
                startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
                toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);

                recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC + startScIdx];
                recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];

                genDmrs_rePtr = &pNrPuschOutParams->genDMRS[0][startScIdx];
                genDmrs_imPtr = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];

                cdm1Est_rePtr = &cdm1EstRe[idxAnt][startScIdx];
                cdm1Est_imPtr = &cdm1EstIm[idxAnt][startScIdx];

                for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                {
                    remSc = endScIdx - idxSc;
                    if( remSc < 8)
                    {
                        rDmrs_re = _mm256_setzero_ps();
                        rDmrs_im = _mm256_setzero_ps();

                        gDmrs_re = _mm256_setzero_ps();
                        gDmrs_im = _mm256_setzero_ps();

                        for(i = 0; i < remSc; i++)
                        {
                            rDmrs_re_f[i] = recDmrs_rePtr[i];
                            rDmrs_im_f[i] = recDmrs_imPtr[i];

                            gDmrs_re_f[i] = genDmrs_rePtr[i];
                            gDmrs_im_f[i] = genDmrs_imPtr[i];
                        }

                        bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                        ad_im = _mm256_mul_ps(rDmrs_re, gDmrs_im);

                        out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                        out_im = _mm256_fmsub_ps(rDmrs_im, gDmrs_re, ad_im);

                        // out = ( Y*conj(Xdmrs) )/beta_DMRS;
                        out_re = _mm256_mul_ps(out_re, div2Vec);
                        out_im = _mm256_mul_ps(out_im, div2Vec);

                        for(i = 0; i < remSc; i++)
                        {
                            cdm1Est_rePtr[i] = out_re_f[i];
                            cdm1Est_imPtr[i] = out_im_f[i];
                        }
                    }
                    else
                    {
                        rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                        rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);

                        gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                        gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);

                        bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                        ad_im = _mm256_mul_ps(rDmrs_re, gDmrs_im);

                        out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                        out_im = _mm256_fmsub_ps(rDmrs_im, gDmrs_re, ad_im);

                        // out = ( Y*conj(Xdmrs) )/beta_DMRS;
                        out_re = _mm256_mul_ps(out_re, div2Vec);
                        out_im = _mm256_mul_ps(out_im, div2Vec);

                        _mm256_storeu_ps(cdm1Est_rePtr, out_re);
                        _mm256_storeu_ps(cdm1Est_imPtr, out_im);
                    }

                    recDmrs_rePtr += 8;
                    recDmrs_imPtr += 8;

                    genDmrs_rePtr += 8;
                    genDmrs_imPtr += 8;

                    cdm1Est_rePtr += 8;
                    cdm1Est_imPtr += 8;
                }
            }
        }


        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
            {
                for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                {
                    startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                    endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                    toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
                    toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);

                    cdm1Est_rePtr = &cdm1EstRe[idxAnt][startScIdx];
                    cdm1Est_imPtr = &cdm1EstIm[idxAnt][startScIdx];

                    // Channel Estimates after Linear Interpolation
                    chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                    chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                    if(pNrPuschInParams->nPortIndex[idxLyr] == 2)
                    {
                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;

                            if( remSc <= 8)
                            {
                                for(i = 0; i < remSc-1; i++)
                                {
                                    chEst_rePtr[2*i] = (cdm1Est_rePtr[i]+cdm1Est_rePtr[i+1])*0.5;
                                    chEst_rePtr[2*i+1] = 0;
                                    chEst_imPtr[2*i] = (cdm1Est_imPtr[i]+cdm1Est_imPtr[i+1])*0.5;
                                    chEst_imPtr[2*i+1] = 0;
                                }

                                chEst_rePtr[2*(remSc-1)] = chEst_rePtr[2*(remSc-2)];
                                chEst_imPtr[2*(remSc-1)] = chEst_imPtr[2*(remSc-2)];
                            }
                            else
                            {   
                                cdm1EstRe1 = _mm256_loadu_ps(cdm1Est_rePtr);
                                cdm1EstIm1 = _mm256_loadu_ps(cdm1Est_imPtr);
                                cdm1EstRe2 = _mm256_loadu_ps(cdm1Est_rePtr+1);
                                cdm1EstIm2 = _mm256_loadu_ps(cdm1Est_imPtr+1);

                                out_re = _mm256_add_ps(cdm1EstRe1, cdm1EstRe2);
                                out_im = _mm256_add_ps(cdm1EstIm1, cdm1EstIm2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);

                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);

                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);

                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);

                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );

                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                
                            }
                            cdm1Est_rePtr += 8;
                            cdm1Est_imPtr += 8;
                            
                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }
                        

                        //Frequency Interpolation.
                        //Calculate estimates for null tones using Left and Right estimates.
                        //Inp : [0 a 0 b 0 c ...]
                        //Out : [(a+(a+b)/2)/2 a (a+b)/2 b (b+c)/2 c ...]
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];
                        AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                        AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                        for(idxSc = toneStartIdx+1; idxSc < toneEndIdx; idxSc += 8)
                        {
                            remSc = toneEndIdx - idxSc;
                            if(remSc<8)
                            {
                                for(i = 0; i < remSc; i += 2)
                                {
                                    chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                    chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                }

                                pNrPuschOutParams->chEst_re[idxAnt][idxLyr][0] = \
                                (pNrPuschOutParams->chEst_re[idxAnt][idxLyr][1]+\
                                pNrPuschOutParams->chEst_re[idxAnt][idxLyr][2])*0.5;

                                pNrPuschOutParams->chEst_im[idxAnt][idxLyr][0] = \
                                (pNrPuschOutParams->chEst_im[idxAnt][idxLyr][1]+\
                                pNrPuschOutParams->chEst_im[idxAnt][idxLyr][2])*0.5;
                            }
                            else
                            {
                                AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                chEst_rePtr += 8;
                                chEst_imPtr += 8;
                            }
                        }
                    }
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 3)
                    {
                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc <= 8)
                            {
                                for(i = 0; i < remSc-1; i++)
                                {
                                    chEst_rePtr[2*i] = alternateSigns_f[i]*(cdm1Est_rePtr[i]-cdm1Est_rePtr[i+1])*0.5;
                                    chEst_rePtr[2*i+1] = 0;
                                    chEst_imPtr[2*i] = alternateSigns_f[i]*(cdm1Est_imPtr[i]-cdm1Est_imPtr[i+1])*0.5;
                                    chEst_imPtr[2*i+1] = 0;
                                }
                                chEst_rePtr[2*(remSc-1)] = chEst_rePtr[2*(remSc-2)];
                                chEst_imPtr[2*(remSc-1)] = chEst_imPtr[2*(remSc-2)];
                            }
                            else
                            {   
                                cdm1EstRe1 = _mm256_loadu_ps(cdm1Est_rePtr);
                                cdm1EstIm1 = _mm256_loadu_ps(cdm1Est_imPtr);
                                cdm1EstRe2 = _mm256_loadu_ps(cdm1Est_rePtr+1);
                                cdm1EstIm2 = _mm256_loadu_ps(cdm1Est_imPtr+1);

                                out_re = _mm256_sub_ps(cdm1EstRe1, cdm1EstRe2);
                                out_re = _mm256_mul_ps(out_re, alternateSigns);
                                out_im = _mm256_sub_ps(cdm1EstIm1, cdm1EstIm2);
                                out_im = _mm256_mul_ps(out_im, alternateSigns);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);

                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);

                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);

                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);

                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );

                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                
                            }
                            cdm1Est_rePtr += 8;
                            cdm1Est_imPtr += 8;
                            
                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }

                        //Frequency Interpolation.
                        //Calculate estimates for null tones using Left and Right estimates.
                        //Inp : [0 a 0 b 0 c ...]
                        //Out : [(a+(a+b)/2)/2 a (a+b)/2 b (b+c)/2 c ...]
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];
                        AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                        AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                        for(idxSc = toneStartIdx+1; idxSc < toneEndIdx; idxSc += 8)
                        {
                            remSc = toneEndIdx - idxSc;
                            if(remSc<8)
                            {
                                for(i = 0; i < remSc; i += 2)
                                {
                                    chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                    chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                }

                                pNrPuschOutParams->chEst_re[idxAnt][idxLyr][0] = \
                                (pNrPuschOutParams->chEst_re[idxAnt][idxLyr][1]+\
                                pNrPuschOutParams->chEst_re[idxAnt][idxLyr][2])*0.5;

                                pNrPuschOutParams->chEst_im[idxAnt][idxLyr][0] = \
                                (pNrPuschOutParams->chEst_im[idxAnt][idxLyr][1]+\
                                pNrPuschOutParams->chEst_im[idxAnt][idxLyr][2])*0.5;
                            }
                            else
                            {
                                AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                chEst_rePtr += 8;
                                chEst_imPtr += 8;
                            }
                        }
                    }
                }
            }
        }

    }

    //If there exit single layer operation or two layer FDM tranmission.
    if(nLyr == 1 || (nLyr == 2 && cdmTotal[0] == 1))
    {
        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
            {
                for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                {
                    startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                    endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                    toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
                    toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);

                    // port 0 processing
                    if(pNrPuschInParams->nPortIndex[idxLyr] == 0)
                    {          
                        // Received DM-RS Signal                          
                        recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                        recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][(MAX_NUM_SC/2) + startScIdx];
                        
                        // Generated DM-RS Signal
                        genDmrs_rePtr = &pNrPuschOutParams->genDMRS[0][startScIdx];
                        genDmrs_imPtr = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                        
                        // Channel Estimates after Linear Interpolation
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];

                        // Modulation Removal
                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc < 8)
                            {
                                rDmrs_re = _mm256_setzero_ps();
                                rDmrs_im = _mm256_setzero_ps();

                                gDmrs_re = _mm256_setzero_ps();
                                gDmrs_im = _mm256_setzero_ps();

                                for(i = 0; i < remSc; i++)
                                {
                                    rDmrs_re_f[i] = recDmrs_rePtr[i];
                                    rDmrs_im_f[i] = recDmrs_imPtr[i];

                                    gDmrs_re_f[i] = genDmrs_rePtr[i];
                                    gDmrs_im_f[i] = genDmrs_imPtr[i];
                                }

                                bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                                ad_im = _mm256_mul_ps(rDmrs_re, gDmrs_im);

                                out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                                out_im = _mm256_fmsub_ps(rDmrs_im, gDmrs_re, ad_im);

                                // out = ( Y*conj(Xdmrs) )/Beta_dmrs;
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);
                                
                                // Pad zeros between estimates [h0 0 h1 0 h2 ...].
                                for(i = 0; i < remSc; i++)
                                {
                                    chEst_rePtr[2*i] = out_re_f[i];
                                    chEst_rePtr[2*i+1] = 0;

                                    chEst_imPtr[2*i] = out_im_f[i];
                                    chEst_imPtr[2*i+1] = 0;
                                }
                            }
                            else
                            {
                                rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                                rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);

                                gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);

                                bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                                ad_im = _mm256_mul_ps(rDmrs_re, gDmrs_im);

                                out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                                out_im = _mm256_fmsub_ps(rDmrs_im, gDmrs_re, ad_im);

                                // out = ( Y*conj(Xdmrs) )/Beta_dmrs;
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                // Pad zeros between estimates [h0 0 h1 0 h2 ...].
                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                            }

                            recDmrs_rePtr += 8;
                            recDmrs_imPtr += 8;

                            genDmrs_rePtr += 8;
                            genDmrs_imPtr += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }//for idxSc


                        //Frequency Interpolation 
                        //Calculate estimates for null tones using Left and Right estimates.
                        //Inp : [ a 0 b 0 c ...]
                        //Out : [a (a+b)/2 b (b+c)/2 c...]
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];
                        AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                        AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                        for(idxSc = toneStartIdx; idxSc < toneEndIdx; idxSc += 8)
                        {
                            remSc = toneEndIdx - idxSc;
                            if(remSc<=8)
                            {
                                for(i = 0; i < remSc; i += 2)
                                {
                                    if(remSc - i == 2)
                                    {
                                        chEst_rePtr[i+1] = (chEst_rePtr[i]+chEst_rePtr[i-1])*0.5;
                                        chEst_imPtr[i+1] = (chEst_imPtr[i]+chEst_imPtr[i-1])*0.5;
                                    }
                                    else
                                    {
                                        chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                        chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                    }
                                }
                            }
                            else
                            {
                                AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                chEst_rePtr += 8;
                                chEst_imPtr += 8;
                            }
                        }                    
                    }
                    // port 1 processing
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 1)
                    {
                        recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                        recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][(MAX_NUM_SC/2) + startScIdx];

                        genDmrs_rePtr = &pNrPuschOutParams->genDMRS[0][startScIdx];
                        genDmrs_imPtr = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];

                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];

                        
                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if(remSc<8)
                            {
                                rDmrs_re = _mm256_setzero_ps();
                                rDmrs_im = _mm256_setzero_ps();

                                gDmrs_re = _mm256_setzero_ps();
                                gDmrs_im = _mm256_setzero_ps();

                                for(i = 0; i < remSc; i++)
                                {
                                    rDmrs_re_f[i] = recDmrs_rePtr[i];
                                    rDmrs_im_f[i] = recDmrs_imPtr[i];

                                    gDmrs_re_f[i] = alternateSigns_f[i]*genDmrs_rePtr[i];
                                    gDmrs_im_f[i] = alternateSigns_f[i]*genDmrs_imPtr[i];
                                }

                                bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                                ad_im = _mm256_mul_ps(rDmrs_re, gDmrs_im);

                                out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                                out_im = _mm256_fmsub_ps(rDmrs_im, gDmrs_re, ad_im);

                                // out = ( Y*conj(Xdmrs) )/beta_dmrs;
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                // Pad zeros between estimates [h0 0 h1 0 h2 ...].
                                for(i = 0; i < remSc; i++)
                                {
                                    chEst_rePtr[2*i] = out_re_f[i];
                                    chEst_rePtr[2*i+1] = 0;

                                    chEst_imPtr[2*i] = out_im_f[i];
                                    chEst_imPtr[2*i+1] = 0;
                                }
                            }
                            else
                            {
                                rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                                rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);

                                gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                gDmrs_re = _mm256_mul_ps(gDmrs_re, alternateSigns);
                                gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                gDmrs_im = _mm256_mul_ps(gDmrs_im, alternateSigns);

                                bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                                ad_im = _mm256_mul_ps(rDmrs_re, gDmrs_im);

                                out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                                out_im = _mm256_fmsub_ps(rDmrs_im, gDmrs_re, ad_im);

                                // out = ( Y*conj(Xdmrs) )/beta_dmrs;
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                // Pad zeros between estimates [h0 0 h1 0 h2 ...].
                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                            }

                            recDmrs_rePtr += 8;
                            recDmrs_imPtr += 8;

                            genDmrs_rePtr += 8;
                            genDmrs_imPtr += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }

                        //Frequency Interpolation
                        //Calculate estimates for null tones using Left and Right estimates.
                        //For dead end null tone, its two adjacent estimates will be used.
                        //Inp : [ a 0 b 0 c ...]
                        //Out : [a (a+b)/2 b (b+c)/2 c ...]
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];

                        AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                        AVX_chEst_imPtr = (__m256*)chEst_imPtr;


                        for(idxSc = toneStartIdx; idxSc < toneEndIdx; idxSc += 8)
                        {
                            remSc = toneEndIdx - idxSc;
                            if(remSc<=8)
                            {
                                for(i = 0; i < remSc; i += 2)
                                {
                                    if(remSc - i == 2)
                                    {
                                        chEst_rePtr[i+1] = (chEst_rePtr[i]+chEst_rePtr[i-1])*0.5;
                                        chEst_imPtr[i+1] = (chEst_imPtr[i]+chEst_imPtr[i-1])*0.5;
                                    }
                                    else
                                    {
                                        chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                        chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                    }
                                }
                            }
                            else
                            {
                                AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                chEst_rePtr += 8;
                                chEst_imPtr += 8;
                            }
                        }
                    }
                    // port 2 processing
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 2)
                    {
                        recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC + startScIdx];
                        recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];

                        genDmrs_rePtr = &pNrPuschOutParams->genDMRS[0][startScIdx];
                        genDmrs_imPtr = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];

                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc < 8)
                            {
                                rDmrs_re = _mm256_setzero_ps();
                                rDmrs_im = _mm256_setzero_ps();

                                gDmrs_re = _mm256_setzero_ps();
                                gDmrs_im = _mm256_setzero_ps();

                                for(i = 0; i < remSc; i++)
                                {
                                    rDmrs_re_f[i] = recDmrs_rePtr[i];
                                    rDmrs_im_f[i] = recDmrs_imPtr[i];

                                    gDmrs_re_f[i] = genDmrs_rePtr[i];
                                    gDmrs_im_f[i] = genDmrs_imPtr[i];
                                }

                                bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                                ad_im = _mm256_mul_ps(rDmrs_re, gDmrs_im);

                                out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                                out_im = _mm256_fmsub_ps(rDmrs_im, gDmrs_re, ad_im);

                                // out = ( Y*conj(Xdmrs) )/beta_DMRS;
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                // Pad zeros between estimates [0 h0 0 h1 0 h2 ...].
                                for(i = 0; i < remSc; i++)
                                {
                                    chEst_rePtr[2*i] = out_re_f[i];
                                    chEst_rePtr[2*i+1] = 0;

                                    chEst_imPtr[2*i] = out_im_f[i];
                                    chEst_imPtr[2*i+1] = 0;
                                }
                            }
                            else
                            {
                                rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                                rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);

                                gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);

                                bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                                ad_im = _mm256_mul_ps(rDmrs_re, gDmrs_im);

                                out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                                out_im = _mm256_fmsub_ps(rDmrs_im, gDmrs_re, ad_im);

                                // out = ( Y*conj(Xdmrs) )/beta_DMRS;
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                // Pad zeros between estimates [0 h0 0 h1 0 h2 ...].
                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                            }

                            recDmrs_rePtr += 8;
                            recDmrs_imPtr += 8;

                            genDmrs_rePtr += 8;
                            genDmrs_imPtr += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }

                        //Frequency Interpolation
                        //Calculate estimates for null tones using Left and Right estimates.
                        //For dead end null tone, its two adjacent estimates will be used.
                        //Inp : [0 a 0 b 0 c ...]
                        //Out : [(a+((a+b)/2))/2 a (a+b)/2 b (b+c)/2 c ...]

                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                        AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                        AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                        for(idxSc = toneStartIdx+1; idxSc < toneEndIdx; idxSc += 8)
                        {
                            remSc = toneEndIdx - idxSc;
                            if(remSc<8)
                            {
                                for(i = 0; i < remSc; i += 2)
                                {
                                    chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                    chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                }

                                pNrPuschOutParams->chEst_re[idxAnt][idxLyr][0] = \
                                (pNrPuschOutParams->chEst_re[idxAnt][idxLyr][1]+\
                                pNrPuschOutParams->chEst_re[idxAnt][idxLyr][2])*0.5;

                                pNrPuschOutParams->chEst_im[idxAnt][idxLyr][0] = \
                                (pNrPuschOutParams->chEst_im[idxAnt][idxLyr][1]+\
                                pNrPuschOutParams->chEst_im[idxAnt][idxLyr][2])*0.5;
                            }
                            else
                            {
                                AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);


                                AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                chEst_rePtr += 8;
                                chEst_imPtr += 8;
                            }
                        }
                    }//*/
                    //port 3 processing
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 3)
                    {
                        recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC + startScIdx];
                        recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];

                        genDmrs_rePtr = &pNrPuschOutParams->genDMRS[0][startScIdx];
                        genDmrs_imPtr = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];

                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc < 8)
                            {
                                rDmrs_re = _mm256_setzero_ps();
                                rDmrs_im = _mm256_setzero_ps();

                                gDmrs_re = _mm256_setzero_ps();
                                gDmrs_im = _mm256_setzero_ps();

                                for(i = 0; i < remSc; i++)
                                {
                                    rDmrs_re_f[i] = recDmrs_rePtr[i];
                                    rDmrs_im_f[i] = recDmrs_imPtr[i];

                                    gDmrs_re_f[i] = alternateSigns_f[i]*genDmrs_rePtr[i];
                                    gDmrs_im_f[i] = alternateSigns_f[i]*genDmrs_imPtr[i];
                                }

                                bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                                ad_im = _mm256_mul_ps(rDmrs_re, gDmrs_im);

                                out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                                out_im = _mm256_fmsub_ps(rDmrs_im, gDmrs_re, ad_im);

                                // out = ( Y*conj(Xdmrs) )/beta_DMRS;
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                // Pad zeros between estimates [0 h0 0 h1 0 h2 ...].
                                for(i = 0; i < remSc; i++)
                                {
                                    chEst_rePtr[2*i] = out_re_f[i];
                                    chEst_rePtr[2*i+1] = 0;

                                    chEst_imPtr[2*i] = out_im_f[i];
                                    chEst_imPtr[2*i+1] = 0;
                                }
                            }
                            else
                            {
                                rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                                rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);

                                gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                gDmrs_re = _mm256_mul_ps(gDmrs_re, alternateSigns);
                                gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                gDmrs_im = _mm256_mul_ps(gDmrs_im, alternateSigns);

                                bd_re = _mm256_mul_ps(rDmrs_im, gDmrs_im);
                                ad_im = _mm256_mul_ps(rDmrs_re, gDmrs_im);

                                out_re = _mm256_fmadd_ps(rDmrs_re, gDmrs_re, bd_re);
                                out_im = _mm256_fmsub_ps(rDmrs_im, gDmrs_re, ad_im);

                                // out = ( Y*conj(Xdmrs) )/beta_DMRS;
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                // Pad zeros between estimates [0 h0 0 h1 0 h2 ...].
                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                            }
                            recDmrs_rePtr += 8;
                            recDmrs_imPtr += 8;

                            genDmrs_rePtr += 8;
                            genDmrs_imPtr += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }


                        //Frequency Interpolation
                        //Calculate estimates for null tones using Left and Right estimates.
                        //For dead end null tone, its two adjacent estimates will be used.
                        //Inp : [0 a 0 b 0 c ...]
                        //Out : [(a+((a+b)/2))/2 a (a+b)/2 b (b+c)/2 c ...]

                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                        AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                        AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                        for(idxSc = toneStartIdx+1; idxSc < toneEndIdx; idxSc += 8)
                        {
                            remSc = toneEndIdx - idxSc;
                            if(remSc<8)
                            {
                                for(i = 0; i < remSc; i += 2)
                                {
                                    chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                    chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                }

                                pNrPuschOutParams->chEst_re[idxAnt][idxLyr][0] = \
                                (pNrPuschOutParams->chEst_re[idxAnt][idxLyr][1]+\
                                pNrPuschOutParams->chEst_re[idxAnt][idxLyr][2])*0.5;

                                pNrPuschOutParams->chEst_im[idxAnt][idxLyr][0] = \
                                (pNrPuschOutParams->chEst_im[idxAnt][idxLyr][1]+\
                                pNrPuschOutParams->chEst_im[idxAnt][idxLyr][2])*0.5;
                            }
                            else
                            {
                                AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                chEst_rePtr += 8;
                                chEst_imPtr += 8;
                            }
                        }
                    }
                }
            }
        }
    }
    

    #ifdef DEBUG
        FILE *fp1;
        if( segPart == 1)
        {
            fp1 = fopen("IOs/modulationRemovalSeg1_mmse_PerTone.txt", "w");
        }
        else if( segPart == 2)
        {
            fp1 = fopen("IOs/modulationRemovalSeg2_mmse_PerTone.txt", "w");
        }
        else if( segPart == 3)
        {
            fp1 = fopen("IOs/modulationRemovalSeg3_mmse_PerTone.txt", "w");
        }
        else //if( segPart == 4)
        {
            fp1 = fopen("IOs/modulationRemovalSeg4_mmse_PerTone.txt", "w");
        }

        for(int idxAnt = 0 ;idxAnt<pNrPuschInParams->nAntennaPorts;idxAnt++)
        {
            for(int idxLyr = 0;idxLyr<pNrPuschInParams->nNrOfLayers;idxLyr++)
            {
                for(int idxScTest = 0; idxScTest < 3300; idxScTest++)
                {
                    fprintf(fp1 , "%f\n",  pNrPuschOutParams->lsEst[idxAnt][idxLyr][idxScTest]);
                }
            }
        }
        fclose(fp1);
    
        FILE *fp2, *pf2;
        if( segPart == 1)
        {
            fp2 = fopen("IOs/estimatesRealInterpolatedSeg1_mmse_PerTone.txt", "w");
            pf2 = fopen("IOs/estimatesImagInterpolatedSeg1_mmse_PerTone.txt", "w");
        }
        else if( segPart == 2)
        {
            fp2 = fopen("IOs/estimatesRealInterpolatedSeg2_mmse_PerTone.txt", "w");
            pf2 = fopen("IOs/estimatesImagInterpolatedSeg2_mmse_PerTone.txt", "w");
        }
        else if( segPart == 3)
        {
            fp2 = fopen("IOs/estimatesRealInterpolatedSeg3_mmse_PerTone.txt", "w");
            pf2 = fopen("IOs/estimatesImagInterpolatedSeg3_mmse_PerTone.txt", "w");
        }
        else //if( segPart == 4)
        {
            fp2 = fopen("IOs/estimatesRealInterpolatedSeg4_mmse_PerTone.txt", "w");
            pf2 = fopen("IOs/estimatesImagInterpolatedSeg4_mmse_PerTone.txt", "w");
        }

        for(int idxAntTest = 0; idxAntTest < pNrPuschInParams->nAntennaPorts; idxAntTest++)
        {
            for(int idxLyrTest = 0; idxLyrTest < pNrPuschInParams->nNrOfLayers; idxLyrTest++)
            {
                for(int idxScTest = 0; idxScTest < 3300; idxScTest++)
                {
                    fprintf(fp2, "%f\n",\
                    pNrPuschOutParams->chEst_re[idxAntTest][idxLyrTest][idxScTest]);
                    fprintf(pf2, "%f\n",\
                    pNrPuschOutParams->chEst_im[idxAntTest][idxLyrTest][idxScTest]);
                }
            }
        }
        fclose(fp2);
        fclose(pf2);
    #endif
}

void nr_pusch_mmse_perTone_FO_est(puschConfigT* pNrPuschInParams,\
                                 P_NR_PUSCH_OUT_PARAMS pNrPuschOutParamsSeg1,\
                                 P_NR_PUSCH_OUT_PARAMS pNrPuschOutParamsSeg2)
{
    wnUInt16 idxAnt, idxLyr, idxAlloc, idxSc;//, i, PRGIdx
    wnUInt16 remSc;//startScIdx, endScIdx,
    wnUInt8 nRxAnt, nLyr, nSbAlloc;//dmrsSym,
    wnUInt16 toneStartIdx, toneEndIdx;

    nRxAnt   = pNrPuschInParams->nAntennaPorts;
    nLyr     = pNrPuschInParams->nNrOfLayers;
    nSbAlloc = pNrPuschInParams->nSbAlloc;

    // __m256* AVX_chEst_rePtrSeg1, *AVX_chEst_imPtrSeg1;
    // __m256* AVX_chEst_rePtrSeg2, *AVX_chEst_imPtrSeg2;
    wnFlt *chEst_rePtrSeg1, *chEst_imPtrSeg1;
    wnFlt *chEst_rePtrSeg2, *chEst_imPtrSeg2;
    __m256 AVX_chEst_reSeg1, AVX_chEst_imSeg1;
    // wnFlt *AVX_chEst_reSeg1_f = (wnFlt *)&AVX_chEst_reSeg1;
    // wnFlt *AVX_chEst_imSeg1_f = (wnFlt *)&AVX_chEst_imSeg1;

    __m256 AVX_chEst_reSeg2, AVX_chEst_imSeg2;
    // wnFlt *AVX_chEst_reSeg2_f = (wnFlt *)&AVX_chEst_reSeg2;
    // wnFlt *AVX_chEst_imSeg2_f = (wnFlt *)&AVX_chEst_imSeg2;

    __m256 bd_re, ad_im, out_re, out_im;

    __m256 FOEstReal = _mm256_set1_ps(0.0);
    __m256 FOEstImag = _mm256_set1_ps(0.0);
    wnFlt* FOEstReal_f = (wnFlt*)&FOEstReal;
    wnFlt* FOEstImag_f = (wnFlt*)&FOEstImag;

    wnFlt chEst_re_Seg1_f, chEst_im_Seg1_f;
    wnFlt chEst_re_Seg2_f, chEst_im_Seg2_f;
    wnUInt16 AVXitrCnt;
    wnUInt32 idxSCLowerLimit;
    // wnUInt8 remTones;

    wnFlt OneBydmrsPosDiff = 1/(wnFlt)(pNrPuschInParams->dmrs2_sym - pNrPuschInParams->dmrs1_sym);

    if(pNrPuschInParams->FOFullAlloc == 1)
    {
        //FO estimate for entire allocation.
        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
            {
                for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                {
                    toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
                    toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);

                    chEst_rePtrSeg1 = &pNrPuschOutParamsSeg1->chEst_re[idxAnt][idxLyr][toneStartIdx];
                    chEst_imPtrSeg1 = &pNrPuschOutParamsSeg1->chEst_im[idxAnt][idxLyr][toneStartIdx];

                    chEst_rePtrSeg2 = &pNrPuschOutParamsSeg2->chEst_re[idxAnt][idxLyr][toneStartIdx];
                    chEst_imPtrSeg2 = &pNrPuschOutParamsSeg2->chEst_im[idxAnt][idxLyr][toneStartIdx];

                    AVXitrCnt = ((toneEndIdx - toneStartIdx)>>3);

                    // if (toneEndIdx - toneStartIdx) = 8*X+Y, following for-loop will complete 8*X part.
                    for(idxSc = 0; idxSc < AVXitrCnt; idxSc++)
                    {
                        AVX_chEst_reSeg1 = _mm256_loadu_ps(chEst_rePtrSeg1);
                        AVX_chEst_imSeg1 = _mm256_loadu_ps(chEst_imPtrSeg1);

                        AVX_chEst_reSeg2 = _mm256_loadu_ps(chEst_rePtrSeg2);
                        AVX_chEst_imSeg2 = _mm256_loadu_ps(chEst_imPtrSeg2);

                        bd_re = _mm256_mul_ps(AVX_chEst_imSeg2, AVX_chEst_imSeg1);
                        ad_im = _mm256_mul_ps(AVX_chEst_reSeg2, AVX_chEst_imSeg1);

                        out_re = _mm256_fmadd_ps(AVX_chEst_reSeg2, AVX_chEst_reSeg1, bd_re);
                        out_im = _mm256_fmsub_ps(AVX_chEst_imSeg2, AVX_chEst_reSeg1, ad_im);

                        FOEstReal = _mm256_add_ps(FOEstReal, out_re);
                        FOEstImag = _mm256_add_ps(FOEstImag, out_im);
                        
                        chEst_rePtrSeg1 += 8;
                        chEst_imPtrSeg1 += 8;
                        chEst_rePtrSeg2 += 8;
                        chEst_imPtrSeg2 += 8;
                    }

                    #if(0)
                        remTones = (toneEndIdx - (toneStartIdx+(AVXitrCnt<<3)));
                        // if (toneEndIdx - toneStartIdx) = 8*X+Y, following if-loop will complete Y part.
                        if( remTones > 0)
                        {
                            AVX_chEst_reSeg1 = _mm256_set1_ps(0.0);
                            AVX_chEst_imSeg1 = _mm256_set1_ps(0.0);

                            AVX_chEst_reSeg2 = _mm256_set1_ps(0.0);
                            AVX_chEst_imSeg2 = _mm256_set1_ps(0.0);

                            for(idxSc = 0;idxSc<remTones;idxSc++)
                            {
                                AVX_chEst_reSeg1_f[idxSc] = pNrPuschOutParamsSeg1->chEst_re[idxAnt][idxLyr][idxSc];
                                AVX_chEst_imSeg1_f[idxSc] = pNrPuschOutParamsSeg1->chEst_im[idxAnt][idxLyr][idxSc];
                                AVX_chEst_reSeg2_f[idxSc] = pNrPuschOutParamsSeg2->chEst_re[idxAnt][idxLyr][idxSc];
                                AVX_chEst_imSeg2_f[idxSc] = pNrPuschOutParamsSeg2->chEst_im[idxAnt][idxLyr][idxSc];
                            }

                            bd_re = _mm256_mul_ps(AVX_chEst_imSeg2, AVX_chEst_imSeg1);
                            ad_im = _mm256_mul_ps(AVX_chEst_reSeg2, AVX_chEst_imSeg1);

                            out_re = _mm256_fmadd_ps(AVX_chEst_reSeg2, AVX_chEst_reSeg1, bd_re);
                            out_im = _mm256_fmsub_ps(AVX_chEst_imSeg2, AVX_chEst_reSeg1, ad_im);

                            FOEstReal = _mm256_add_ps(FOEstReal, out_re);
                            FOEstImag = _mm256_add_ps(FOEstImag, out_im);
                        }
                    #else
                    // ELSE part is faster than IF part.
                        for(idxSc = toneStartIdx+(AVXitrCnt<<3);idxSc<toneEndIdx;idxSc++)
                        {
                            chEst_re_Seg1_f = pNrPuschOutParamsSeg1->chEst_re[idxAnt][idxLyr][idxSc];
                            chEst_im_Seg1_f = pNrPuschOutParamsSeg1->chEst_im[idxAnt][idxLyr][idxSc];
                            chEst_re_Seg2_f = pNrPuschOutParamsSeg2->chEst_re[idxAnt][idxLyr][idxSc];
                            chEst_im_Seg2_f = pNrPuschOutParamsSeg2->chEst_im[idxAnt][idxLyr][idxSc];

                            *FOEstReal_f = *FOEstReal_f + chEst_re_Seg2_f*chEst_re_Seg1_f + chEst_im_Seg2_f*chEst_im_Seg1_f;
                            *FOEstImag_f = *FOEstImag_f + chEst_im_Seg2_f*chEst_re_Seg1_f - chEst_re_Seg2_f*chEst_im_Seg1_f;
                        }
                    #endif
                }
            }
        }

        for(idxSc = 1;idxSc<8;idxSc++)
        {
            FOEstReal_f[0] = FOEstReal_f[0]+FOEstReal_f[idxSc];
            FOEstImag_f[0] = FOEstImag_f[0]+FOEstImag_f[idxSc];
        }

        // FO Estimate for entire allocation.
        FOEst = atan2(FOEstImag_f[0], FOEstReal_f[0])*OneBydmrsPosDiff;
    }
    else
    {
        wnUInt16 numofPRGs;
        wnUInt32 PRGIdx, idxSCUpperLimit;
        wnUInt16 cmpldPRGs;
        
        cmpldPRGs = 0;
        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
        {
            toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
            toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);

            // numofPRGs = ceil(nPrb[idxAlloc]/FOPRG);
            numofPRGs = (wnUInt16)(pNrPuschInParams->nPrb[idxAlloc]/pNrPuschInParams->FOPRG);
            if((pNrPuschInParams->nPrb[idxAlloc] - (wnUInt16)(numofPRGs*pNrPuschInParams->FOPRG))>0)
            {
                numofPRGs = numofPRGs+1;
            }
             
            for(PRGIdx = 0;PRGIdx<numofPRGs;PRGIdx++)
            {
                FOEstReal = _mm256_set1_ps(0.0);
                FOEstImag = _mm256_set1_ps(0.0);

                idxSCLowerLimit = toneStartIdx+(wnUInt16)(PRGIdx*pNrPuschInParams->FOPRG*12);
                idxSCUpperLimit = toneStartIdx+(wnUInt16)((PRGIdx+1)*pNrPuschInParams->FOPRG*12);
                if(PRGIdx == numofPRGs-1)
                {
                    if(toneEndIdx<idxSCUpperLimit)
                    {
                        idxSCUpperLimit = toneEndIdx;
                    }
                }
                AVXitrCnt = ((idxSCUpperLimit - idxSCLowerLimit)>>3);

                for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                {
                    for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
                    {
                        chEst_rePtrSeg1 = &pNrPuschOutParamsSeg1->chEst_re[idxAnt][idxLyr][(wnUInt16)(PRGIdx*pNrPuschInParams->FOPRG*12+toneStartIdx)];
                        chEst_imPtrSeg1 = &pNrPuschOutParamsSeg1->chEst_im[idxAnt][idxLyr][(wnUInt16)(PRGIdx*pNrPuschInParams->FOPRG*12+toneStartIdx)];

                        chEst_rePtrSeg2 = &pNrPuschOutParamsSeg2->chEst_re[idxAnt][idxLyr][(wnUInt16)(PRGIdx*pNrPuschInParams->FOPRG*12+toneStartIdx)];
                        chEst_imPtrSeg2 = &pNrPuschOutParamsSeg2->chEst_im[idxAnt][idxLyr][(wnUInt16)(PRGIdx*pNrPuschInParams->FOPRG*12+toneStartIdx)];

                        // if (idxSCUpperLimit - idxSCLowerLimit) = 8*X+Y, following for-loop will complete 8*X part.
                        for(idxSc = 0;idxSc<AVXitrCnt;idxSc++)
                        {
                            AVX_chEst_reSeg1 = _mm256_loadu_ps(chEst_rePtrSeg1);
                            AVX_chEst_imSeg1 = _mm256_loadu_ps(chEst_imPtrSeg1);

                            AVX_chEst_reSeg2 = _mm256_loadu_ps(chEst_rePtrSeg2);
                            AVX_chEst_imSeg2 = _mm256_loadu_ps(chEst_imPtrSeg2);

                            bd_re = _mm256_mul_ps(AVX_chEst_imSeg2, AVX_chEst_imSeg1);
                            ad_im = _mm256_mul_ps(AVX_chEst_reSeg2, AVX_chEst_imSeg1);

                            out_re = _mm256_fmadd_ps(AVX_chEst_reSeg2, AVX_chEst_reSeg1, bd_re);
                            out_im = _mm256_fmsub_ps(AVX_chEst_imSeg2, AVX_chEst_reSeg1, ad_im);

                            FOEstReal = _mm256_add_ps(FOEstReal, out_re);
                            FOEstImag = _mm256_add_ps(FOEstImag, out_im);
                            
                            chEst_rePtrSeg1 += 8;
                            chEst_imPtrSeg1 += 8;
                            chEst_rePtrSeg2 += 8;
                            chEst_imPtrSeg2 += 8;
                        }

                        // if (idxSCUpperLimit - idxSCLowerLimit) = 8*X+Y, following for-loop will complete Y part.
                        #if(0)
                            remTones = (idxSCUpperLimit - (idxSCLowerLimit+(AVXitrCnt*8)));
                            if( remTones > 0)
                            {
                                AVX_chEst_reSeg1 = _mm256_set1_ps(0.0);
                                AVX_chEst_imSeg1 = _mm256_set1_ps(0.0);

                                AVX_chEst_reSeg2 = _mm256_set1_ps(0.0);
                                AVX_chEst_imSeg2 = _mm256_set1_ps(0.0);

                                for(idxSc = 0;idxSc<remTones;idxSc++)
                                {
                                    AVX_chEst_reSeg1_f[idxSc] = pNrPuschOutParamsSeg1->chEst_re[idxAnt][idxLyr][idxSc];
                                    AVX_chEst_imSeg1_f[idxSc] = pNrPuschOutParamsSeg1->chEst_im[idxAnt][idxLyr][idxSc];
                                    AVX_chEst_reSeg2_f[idxSc] = pNrPuschOutParamsSeg2->chEst_re[idxAnt][idxLyr][idxSc];
                                    AVX_chEst_imSeg2_f[idxSc] = pNrPuschOutParamsSeg2->chEst_im[idxAnt][idxLyr][idxSc];
                                }

                                bd_re = _mm256_mul_ps(AVX_chEst_imSeg2, AVX_chEst_imSeg1);
                                ad_im = _mm256_mul_ps(AVX_chEst_reSeg2, AVX_chEst_imSeg1);

                                out_re = _mm256_fmadd_ps(AVX_chEst_reSeg2, AVX_chEst_reSeg1, bd_re);
                                out_im = _mm256_fmsub_ps(AVX_chEst_imSeg2, AVX_chEst_reSeg1, ad_im);

                                FOEstReal = _mm256_add_ps(FOEstReal, out_re);
                                FOEstImag = _mm256_add_ps(FOEstImag, out_im);
                            }

                        #else
                        // ELSE part is executing faster than IF part.
                        for(idxSc = idxSCLowerLimit+(AVXitrCnt*8);idxSc<idxSCUpperLimit;idxSc++)
                        {
                            chEst_re_Seg1_f = pNrPuschOutParamsSeg1->chEst_re[idxAnt][idxLyr][idxSc];
                            chEst_im_Seg1_f = pNrPuschOutParamsSeg1->chEst_im[idxAnt][idxLyr][idxSc];
                            chEst_re_Seg2_f = pNrPuschOutParamsSeg2->chEst_re[idxAnt][idxLyr][idxSc];
                            chEst_im_Seg2_f = pNrPuschOutParamsSeg2->chEst_im[idxAnt][idxLyr][idxSc];

                            *FOEstReal_f = *FOEstReal_f + chEst_re_Seg2_f*chEst_re_Seg1_f + chEst_im_Seg2_f*chEst_im_Seg1_f;
                            *FOEstImag_f = *FOEstImag_f + chEst_im_Seg2_f*chEst_re_Seg1_f - chEst_re_Seg2_f*chEst_im_Seg1_f;
                        }
                        #endif
                    }
                }
                
                for(idxSc = 1;idxSc<8;idxSc++)
                {
                    FOEstReal_f[0] = FOEstReal_f[0]+FOEstReal_f[idxSc];
                    FOEstImag_f[0] = FOEstImag_f[0]+FOEstImag_f[idxSc];
                }

                // FO Estimate on PRG basis.
                FOEst_PRG[cmpldPRGs+PRGIdx] = atan2(FOEstImag_f[0], FOEstReal_f[0])*OneBydmrsPosDiff;
            }
            cmpldPRGs = cmpldPRGs + numofPRGs;
        }

        #ifdef DEBUG
            FILE *fpFO_PRG = fopen("IOs/FOPRGPerTone.txt", "w");
            for(wnUInt32 PRGIdx = 0;PRGIdx<cmpldPRGs;PRGIdx++)
            {
                fprintf(fpFO_PRG, "%f\n", FOEst_PRG[PRGIdx]);
            }
            fclose(fpFO_PRG);
        #endif
    }
}


/*void nr_pusch_mmse_perTone_nVarCalc_avx2(wnUInt8 segPart,
                               commonUlConfigT* pNrUlCommonParams,
                               puschConfigT* pNrPuschInParams,
                               P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams)
{
    wnUInt8 cdmTotal[2] = {0, 0};
    wnUInt8 cdm0LyrIdx[2], cdm1LyrIdx[2], idxLyr, cnt0, cnt1, dmrsSym, nRxAnt, nLyr, nSbAlloc;
    wnUInt16 nPrb_total, idxAlloc, idxAnt;
    wnUInt16 nullToneIdx, nullToneStart, nullToneEnd, idxTemp, toneStartIdx, toneEndIdx;
    wnUInt8 remSc;

    cnt0 = 0;
    cnt1 = 0;
    for(idxLyr = 0;idxLyr<pNrPuschInParams->nNrOfLayers;idxLyr++)
    {
        if(pNrPuschInParams->nPortIndex[idxLyr] == 0 || pNrPuschInParams->nPortIndex[idxLyr] == 1)
        {
            cdmTotal[0] = cdmTotal[0]+1;
            if(cnt0 == 0)
            {
                cdm0LyrIdx[0] = idxLyr;
                cdm0LyrIdx[1] = idxLyr+1;
                cnt0++;
            }
            else
            {
                cdm0LyrIdx[1] = idxLyr+1;
            }
        }
        else
        {
            cdmTotal[1] = cdmTotal[1]+1;
            if(cnt1 == 0)
            {
                cdm1LyrIdx[0] = idxLyr;
                cdm1LyrIdx[1] = idxLyr+1;
                cnt1++;
            }
            else
            {
                cdm1LyrIdx[1] = idxLyr+1;
            }
        }
    }

    nRxAnt   = pNrPuschInParams->nAntennaPorts;
    nLyr     = pNrPuschInParams->nNrOfLayers;
    nSbAlloc = pNrPuschInParams->nSbAlloc;

    // Get DM-RS Symbol Index.
    if(segPart == SEG1)
    {
        dmrsSym = pNrPuschInParams->dmrs1_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG2)
    {
        dmrsSym = pNrPuschInParams->dmrs2_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG3)
    {
        dmrsSym = pNrPuschInParams->dmrs3_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG4)
    {
        dmrsSym = pNrPuschInParams->dmrs4_sym - pNrPuschInParams->nStartSymbolIndex;
    }

    nPrb_total = 0;
    for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
    {
        nPrb_total += pNrPuschInParams->nPrb[idxAlloc];
    }
    pNrPuschOutParams->total_nPRB = nPrb_total;



    wnUInt16 toneStartIdx2, PRBStart, PRBEnd;
    wnFlt *nVar_alloc_vec_f, *n_avg_f, *chEst_rePtr, *chEst_imPtr, *temp12_f,\
          *temp22_f, *chEstAvg_rePtr, *chEstAvg_imPtr, *temp11_128_f, *alternateSigns_f,\
          finalValReal, finalValImag;
    
    __m256 nVar_alloc_vec, n_avg, n1_re, n1_im, NVTemp, inpReal1, inpReal2,\
           inpReal3, inpImag1, inpImag2, inpImag3, temp11, temp12, temp21,\
           temp22, oneBy12;
    
    wnUInt16 startScIdx, endScIdx, idxSc, i;
    wnFlt *recDmrs_rePtr, *recDmrs_imPtr, *genDmrs_rePtr, *genDmrs_imPtr, *realOut_f, *imagOut_f,\
          *noise_re_f, *noise_im_f, *rDmrs_re_f, *rDmrs_im_f, *gDmrs_re_f, *gDmrs_im_f,\
          *hAvg_re_f, *hAvg_im_f, *HX_AVX_re_f, *HX_AVX_im_f, *HX_rePtr, *HX_imPtr, HX_re[6 * MAX_NUM_PRB], HX_im[6 * MAX_NUM_PRB];
    __m256 alternateSigns, rDmrs_re, rDmrs_im, gDmrs_re, gDmrs_im, hAvg_re,\
           hAvg_im, bd_re, ad_im, realOut, imagOut, noise_re, noise_im, navg, HX_AVX_re, HX_AVX_im;
    

    nVar_alloc_vec = _mm256_set1_ps(0);
    nVar_alloc_vec_f = (wnFlt *)&nVar_alloc_vec;
    n_avg_f = (wnFlt*)&n_avg;
    temp12_f = (wnFlt*)&temp12;
    temp22_f = (wnFlt*)&temp22;
    oneBy12 = _mm256_set1_ps(0.083333333333333);

    __m128 inpImag1_128, inpImag2_128, inpImag3_128, inpImag4_128, inpImag5_128, inpImag6_128;
    __m128 inpReal1_128, inpReal2_128, inpReal3_128, inpReal4_128, inpReal5_128, inpReal6_128;
    __m128 temp11_128, temp12_128, temp13_128;
    __m128 oneBy12_128 = _mm_set1_ps(0.083333333333333);
    temp11_128_f = (wnFlt*)&temp11_128;

    alternateSigns = _mm256_setr_ps(1, -1, 1, -1, 1, -1, 1, -1);
    alternateSigns_f = (wnFlt *)&alternateSigns;

    realOut_f = (wnFlt*)&realOut;
    imagOut_f = (wnFlt*)&imagOut;

    noise_re_f = (wnFlt*)&noise_re;
    noise_im_f = (wnFlt*)&noise_im;
    navg = _mm256_setzero_ps();
    
    rDmrs_re_f = &rDmrs_re;
    rDmrs_im_f = &rDmrs_im;
    gDmrs_re_f = &gDmrs_re;
    gDmrs_im_f = &gDmrs_im;
    hAvg_re_f = &hAvg_re;
    hAvg_im_f = &hAvg_im;
    
    HX_AVX_re_f = (wnFlt*)&HX_AVX_re;
    HX_AVX_im_f = (wnFlt*)&HX_AVX_im;
    

    if(cdmTotal[1] == 0)// if CDM Group 2 is Empty.
    {
        nVar_alloc_vec = _mm256_setzero_ps();
        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            n_avg = _mm256_setzero_ps();
            for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
            {
                nullToneStart = MAX_NUM_SC;
                nullToneEnd = nullToneStart+ (pNrPuschInParams->nPrb[idxAlloc] * 6);

                for(nullToneIdx = nullToneStart;nullToneIdx<nullToneEnd;nullToneIdx += 8)
                {
                    if(nullToneEnd - nullToneIdx < 8)
                    {
                        remSc = nullToneEnd - nullToneIdx;

                        for(idxTemp = 0;idxTemp<remSc;idxTemp++)
                        {
                            n_avg_f[idxTemp] = n_avg_f[idxTemp]+ (rxFdSamples[idxAnt][dmrsSym][nullToneIdx+idxTemp]*\
                                                                rxFdSamples[idxAnt][dmrsSym][nullToneIdx+idxTemp]);
                        }
                    }
                    else
                    {
                        n1_re = _mm256_loadu_ps(&rxFdSamples[idxAnt][dmrsSym][nullToneIdx]);
                        n_avg = _mm256_fmadd_ps(n1_re, n1_re, n_avg);
                    }
                }

                nullToneStart = (MAX_NUM_SC)+(MAX_NUM_SC/2);
                nullToneEnd = nullToneStart+ (pNrPuschInParams->nPrb[idxAlloc] * 6);
                for(nullToneIdx = nullToneStart;nullToneIdx<nullToneEnd;nullToneIdx += 8)
                {
                    if(nullToneEnd - nullToneIdx <= 8)
                    {
                        remSc = nullToneEnd - nullToneIdx;

                        for(idxTemp = 0;idxTemp<remSc;idxTemp++)
                        {
                            n_avg_f[idxTemp] = n_avg_f[idxTemp]+ (rxFdSamples[idxAnt][dmrsSym][nullToneIdx+idxTemp]*\
                                                                rxFdSamples[idxAnt][dmrsSym][nullToneIdx+idxTemp]);
                        }
                    }
                    else
                    {
                        n1_im = _mm256_loadu_ps(&rxFdSamples[idxAnt][dmrsSym][nullToneIdx]);
                        n_avg = _mm256_fmadd_ps(n1_im, n1_im, n_avg);
                    }
                }
                nVar_alloc_vec = _mm256_add_ps(nVar_alloc_vec, n_avg);
            }
        }

        NVTemp = _mm256_permute2f128_ps (nVar_alloc_vec, nVar_alloc_vec, 1);
        nVar_alloc_vec = _mm256_add_ps(nVar_alloc_vec, NVTemp);
        nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);
        nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);

        pNrPuschOutParams->nVar_overAll = nVar_alloc_vec_f[0] / (nPrb_total * 6 * nRxAnt);
    }
    else if(cdmTotal[0] == 0)// if CDM Group 1 is Empty.
    {
        nVar_alloc_vec = _mm256_setzero_ps();
        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
            {
                nullToneStart = 0;
                n_avg = _mm256_setzero_ps();
                nullToneEnd = nullToneStart+ (pNrPuschInParams->nPrb[idxAlloc] * 6);
                for(nullToneIdx = nullToneStart;nullToneIdx<nullToneEnd;nullToneIdx += 8)
                {
                    if(nullToneEnd - nullToneIdx <= 8)
                    {
                        remSc = nullToneEnd - nullToneIdx;

                        for(idxTemp = 0;idxTemp<remSc;idxTemp++)
                        {
                            n_avg_f[idxTemp] = n_avg_f[idxTemp]+ (rxFdSamples[idxAnt][dmrsSym][nullToneIdx+idxTemp]*\
                                                                rxFdSamples[idxAnt][dmrsSym][nullToneIdx+idxTemp]);
                        }
                    }
                    else
                    {
                        n1_re = _mm256_loadu_ps(&rxFdSamples[idxAnt][dmrsSym][nullToneIdx]);
                        n_avg = _mm256_fmadd_ps(n1_re, n1_re, n_avg);
                    }
                }

                nullToneStart = (MAX_NUM_SC/2);
                nullToneEnd = nullToneStart+ (pNrPuschInParams->nPrb[idxAlloc] * 6);
                for(nullToneIdx = nullToneStart;nullToneIdx<nullToneEnd;nullToneIdx += 8)
                {
                    if(nullToneEnd - nullToneIdx <= 8)
                    {
                        remSc = nullToneEnd - nullToneIdx;

                        for(idxTemp = 0;idxTemp<remSc;idxTemp++)
                        {
                            n_avg_f[idxTemp] = n_avg_f[idxTemp]+ (rxFdSamples[idxAnt][dmrsSym][nullToneIdx+idxTemp]*\
                                                                rxFdSamples[idxAnt][dmrsSym][nullToneIdx+idxTemp]);
                        }
                    }
                    else
                    {
                        n1_im = _mm256_loadu_ps(&rxFdSamples[idxAnt][dmrsSym][nullToneIdx]);
                        n_avg = _mm256_fmadd_ps(n1_im, n1_im, n_avg);
                    }
                }
                nVar_alloc_vec = _mm256_add_ps(nVar_alloc_vec, n_avg);
            }
        }
        NVTemp = _mm256_permute2f128_ps (nVar_alloc_vec, nVar_alloc_vec, 1);
        nVar_alloc_vec = _mm256_add_ps(nVar_alloc_vec, NVTemp);
        nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);
        nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);

        pNrPuschOutParams->nVar_overAll = nVar_alloc_vec_f[0] / (nPrb_total * 6 * nRxAnt);
    }
    else
    {
        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        { 
            for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
            {
                for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                {
                    toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
                    toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);
                    chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                    chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];
                    PRBStart = pNrPuschInParams->rbStart[idxAlloc];
                    PRBEnd = pNrPuschInParams->nPrb[idxAlloc];

                    toneStartIdx2 = pNrPuschInParams->rbStart[idxAlloc] * 6;
                    chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][toneStartIdx2];
                    chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][toneStartIdx2];

                    for(int PRBIdx = PRBStart;PRBIdx<PRBEnd;PRBIdx+=2)
                    {
                        if(toneEndIdx - PRBIdx == 1)
                        {
                            #ifdef PROCESSOne
                                inpReal1_128 = _mm_loadu_ps(chEst_rePtr);
                                inpReal2_128 = _mm_loadu_ps(chEst_rePtr+4);
                                inpReal3_128 = _mm_loadu_ps(chEst_rePtr+8);

                                temp11_128 = _mm_add_ps(inpReal1_128, inpReal2_128);
                                temp11_128 = _mm_add_ps(temp11_128, inpReal3_128);
                                temp11_128 = _mm_hadd_ps(temp11_128, temp11_128);
                                temp11_128 = _mm_hadd_ps(temp11_128, temp11_128);
                                temp11_128 = _mm_mul_ps(temp11_128, oneBy12_128);
                                temp12_128 = _mm_set1_ps(temp11_128_f[0]);
                                _mm_storeu_ps(chEstAvg_rePtr, temp12_128);

                                chEstAvg_rePtr[4] = temp11_128_f[0];
                                chEstAvg_rePtr[5] = temp11_128_f[0];

                                inpImag1_128 = _mm_loadu_ps(chEst_imPtr);
                                inpImag2_128 = _mm_loadu_ps(chEst_imPtr+4);
                                inpImag3_128 = _mm_loadu_ps(chEst_imPtr+8);

                                temp11_128 = _mm_add_ps(inpImag1_128, inpImag2_128);
                                temp11_128 = _mm_add_ps(temp11_128, inpImag3_128);
                                temp11_128 = _mm_hadd_ps(temp11_128, temp11_128);
                                temp11_128 = _mm_hadd_ps(temp11_128, temp11_128);
                                temp11_128 = _mm_mul_ps(temp11_128, oneBy12_128);
                                temp12_128 = _mm_set1_ps(temp11_128_f[0]);
                                _mm_storeu_ps(chEstAvg_imPtr, temp12_128);

                                chEstAvg_imPtr[4] = temp11_128_f[0];
                                chEstAvg_imPtr[5] = temp11_128_f[0];

                            #else
                            // Performance wise 256 Byte vector processing is lagging behind.
                            // Use else part. 
                                inpReal1 = _mm256_loadu_ps(chEst_rePtr);
                                temp12 = _mm256_hadd_ps(inpReal1, inpReal1);
                                temp12 = _mm256_permutevar8x32_ps(temp12, _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7));

                                inpImag1 = _mm256_loadu_ps(chEst_imPtr);
                                temp22 = _mm256_hadd_ps(inpImag1, inpImag1);
                                temp22 = _mm256_permutevar8x32_ps(temp22, _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7));

                                finalValReal = 0;
                                finalValImag = 0;
                                for(int idx = 0;idx<4;idx++)
                                {
                                    finalValReal = finalValReal+temp12_f[idx]+chEst_rePtr[idx+8];
                                    finalValImag = finalValImag+temp22_f[idx]+chEst_imPtr[idx+8];
                                }
                                chEstAvg_rePtr[0] = finalValReal*0.083333333333333;//1/12 = 0.083333333333333;
                                chEstAvg_imPtr[0] = finalValImag*0.083333333333333;
                                for(int idx = 1;idx<6;idx++)
                                {
                                    chEstAvg_rePtr[idx] = chEstAvg_rePtr[0];
                                    chEstAvg_imPtr[idx] = chEstAvg_imPtr[0];
                                }
                            #endif
                        }
                        else
                        {

                            #ifdef PROCESSOne
                                inpReal1_128 = _mm_loadu_ps(chEst_rePtr);
                                inpReal2_128 = _mm_loadu_ps(chEst_rePtr+4);
                                inpReal3_128 = _mm_loadu_ps(chEst_rePtr+8);
                                inpReal4_128 = _mm_loadu_ps(chEst_rePtr+12);
                                inpReal5_128 = _mm_loadu_ps(chEst_rePtr+16);
                                inpReal6_128 = _mm_loadu_ps(chEst_rePtr+20);

                                
                                temp11_128 = _mm_add_ps(inpReal1_128, inpReal2_128);
                                temp11_128 = _mm_add_ps(temp11_128, inpReal3_128);
                                temp11_128 = _mm_hadd_ps(temp11_128, temp11_128);
                                temp11_128 = _mm_hadd_ps(temp11_128, temp11_128);
                                temp11_128 = _mm_mul_ps(temp11_128, oneBy12_128);
                                temp12_128 = _mm_set1_ps(temp11_128_f[0]);
                                _mm_storeu_ps(chEstAvg_rePtr, temp12_128);

                                temp11_128 = _mm_add_ps(inpReal4_128, inpReal5_128);
                                temp11_128 = _mm_add_ps(temp11_128, inpReal6_128);
                                temp11_128 = _mm_hadd_ps(temp11_128, temp11_128);
                                temp11_128 = _mm_hadd_ps(temp11_128, temp11_128);
                                temp11_128 = _mm_mul_ps(temp11_128, oneBy12_128);
                                temp13_128 = _mm_set1_ps(temp11_128_f[0]);
                                temp11_128 = _mm_shuffle_ps(temp12_128, temp13_128, 0b01000100);
                                _mm_storeu_ps(chEstAvg_rePtr+4, temp11_128);
                                _mm_storeu_ps(chEstAvg_rePtr+8, temp13_128);


                                inpImag1_128 = _mm_loadu_ps(chEst_imPtr);
                                inpImag2_128 = _mm_loadu_ps(chEst_imPtr+4);
                                inpImag3_128 = _mm_loadu_ps(chEst_imPtr+8);
                                inpImag4_128 = _mm_loadu_ps(chEst_imPtr+12);
                                inpImag5_128 = _mm_loadu_ps(chEst_imPtr+16);
                                inpImag6_128 = _mm_loadu_ps(chEst_imPtr+20);

                                temp11_128 = _mm_add_ps(inpImag1_128, inpImag2_128);
                                temp11_128 = _mm_add_ps(temp11_128, inpImag3_128);
                                temp11_128 = _mm_hadd_ps(temp11_128, temp11_128);
                                temp11_128 = _mm_hadd_ps(temp11_128, temp11_128);
                                temp11_128 = _mm_mul_ps(temp11_128, oneBy12_128);
                                temp12_128 = _mm_set1_ps(temp11_128_f[0]);
                                _mm_storeu_ps(chEstAvg_imPtr, temp12_128);

                                temp11_128 = _mm_add_ps(inpImag4_128, inpImag5_128);
                                temp11_128 = _mm_add_ps(temp11_128, inpImag6_128);
                                temp11_128 = _mm_hadd_ps(temp11_128, temp11_128);
                                temp11_128 = _mm_hadd_ps(temp11_128, temp11_128);
                                temp11_128 = _mm_mul_ps(temp11_128, oneBy12_128);
                                temp13_128 = _mm_set1_ps(temp11_128_f[0]);
                                temp11_128 = _mm_shuffle_ps(temp12_128, temp13_128, 0b01000100);
                                _mm_storeu_ps(chEstAvg_imPtr+4, temp11_128);
                                _mm_storeu_ps(chEstAvg_imPtr+8, temp13_128);
                            #else
                            // Performance wise 256 Byte vector processing is lagging behind.
                            // Use else part. 
                            inpReal1 = _mm256_loadu_ps(chEst_rePtr);
                            inpReal2 = _mm256_loadu_ps(chEst_rePtr+8);
                            inpReal3 = _mm256_loadu_ps(chEst_rePtr+16);

                            temp11 = _mm256_hadd_ps(inpReal1, inpReal1);
                            temp11 = _mm256_permutevar8x32_ps(temp11, _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7));
                            temp12 = _mm256_hadd_ps(temp11, inpReal2);
                            temp12 = _mm256_hadd_ps(temp12, temp12);
                            temp12 = _mm256_hadd_ps(temp12, temp12);
                            temp12 = _mm256_mul_ps(temp12, oneBy12);
                            _mm256_storeu_ps(chEstAvg_rePtr, _mm256_set1_ps(temp12_f[0]));

                            temp21 = _mm256_hadd_ps(inpReal3, inpReal3);
                            temp21 = _mm256_permutevar8x32_ps(temp21, _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7));
                            temp22 = _mm256_hadd_ps(temp21, inpReal2);
                            temp22 = _mm256_hadd_ps(temp22, temp22);
                            temp22 = _mm256_hadd_ps(temp22, temp22);
                            temp22 = _mm256_mul_ps(temp22, oneBy12);
                            
                            _mm256_storeu_ps(chEstAvg_rePtr+6, _mm256_set1_ps(temp22_f[4]));

                            inpImag1 = _mm256_loadu_ps(chEst_imPtr);
                            inpImag2 = _mm256_loadu_ps(chEst_imPtr+8);
                            inpImag3 = _mm256_loadu_ps(chEst_imPtr+16);

                            temp11 = _mm256_hadd_ps(inpImag1, inpImag1);
                            temp11 = _mm256_permutevar8x32_ps(temp11, _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7));
                            temp12 = _mm256_hadd_ps(temp11, inpImag2);
                            temp12 = _mm256_hadd_ps(temp12, temp12);
                            temp12 = _mm256_hadd_ps(temp12, temp12);
                            temp12 = _mm256_mul_ps(temp12, oneBy12);
                            _mm256_storeu_ps(chEstAvg_imPtr, _mm256_set1_ps(temp12_f[0]));

                            temp21 = _mm256_hadd_ps(inpImag3, inpImag3);
                            temp21 = _mm256_permutevar8x32_ps(temp21, _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7));
                            temp22 = _mm256_hadd_ps(temp21, inpImag2);
                            temp22 = _mm256_hadd_ps(temp22, temp22);
                            temp22 = _mm256_hadd_ps(temp22, temp22);
                            temp22 = _mm256_mul_ps(temp22, oneBy12);
                            _mm256_storeu_ps(chEstAvg_imPtr+6, _mm256_set1_ps(temp22_f[4]));

                            #endif

                            chEstAvg_rePtr += 12;
                            chEstAvg_imPtr += 12;
                            chEst_rePtr += 24;
                            chEst_imPtr += 24;
                        }
                    }
                }

            }
        }

        #ifdef DEBUG
            FILE *fpReal, *fpImag;
            if(segPart == SEG1)
            {
                fpReal = fopen("IOs/HAvgRealSeg1_mmse_PerTone.txt", "w");
                fpImag = fopen("IOs/HAvgImagSeg1_mmse_PerTone.txt", "w");
            }
            else if(segPart == SEG2)
            {
                fpReal = fopen("IOs/HAvgRealSeg2_mmse_PerTone.txt", "w");
                fpImag = fopen("IOs/HAvgImagSeg2_mmse_PerTone.txt", "w");
            }
            else if(segPart == SEG3)
            {
                fpReal = fopen("IOs/HAvgRealSeg3_mmse_PerTone.txt", "w");
                fpImag = fopen("IOs/HAvgImagSeg3_mmse_PerTone.txt", "w");
            }
            else if(segPart == SEG4)
            {
                fpReal = fopen("IOs/HAvgRealSeg4_mmse_PerTone.txt", "w");
                fpImag = fopen("IOs/HAvgImagSeg4_mmse_PerTone.txt", "w");
            }

            for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
            {
                for(int idxLyr = 0; idxLyr < nLyr; idxLyr++)
                {
                    for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                    {
                        toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                        toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                        for(int toneIdx = toneStartIdx;toneIdx<toneEndIdx;toneIdx++)
                        {
                            fprintf(fpReal, "%f\n",\
                            pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][toneIdx]);
                            fprintf(fpImag, "%f\n",\
                            pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][toneIdx]);
                        }
                    }
                }
            }

            fclose(fpReal);
            fclose(fpImag);
        #endif

        if(pNrPuschInParams->nNrOfLayers == 2)
        {
            n_avg = _mm256_setzero_ps();
            nVar_alloc_vec = _mm256_setzero_ps();
            if( (cdm0LyrIdx[1] - cdm0LyrIdx[0]) == 1)// FDM with one layer from CDM Group 0
            {
                // do 
                // 1) NoiseRe + j*NoiseIm = Y-Havg*X; Where, Havg = mean(H(1:12));
                // 2) NV = (NoiseRe)^2 + (NoiseIm)^2;

                idxLyr = cdm0LyrIdx[0];

                if(pNrPuschInParams->nPortIndex[idxLyr] == 0)
                {
                    for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                    {
                        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                        {
                            startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                            endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);
                            recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                            recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][(MAX_NUM_SC/2) + startScIdx];
                            chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][startScIdx];
                            chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][startScIdx];
                            genDmrs_rePtr = &pNrPuschOutParams->genDMRS[0][startScIdx];
                            genDmrs_imPtr = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];

                            
                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                if( endScIdx - idxSc < 8)
                                {
                                    remSc = (endScIdx) - idxSc;

                                    rDmrs_re = _mm256_setzero_ps();
                                    rDmrs_im = _mm256_setzero_ps();
                                    gDmrs_re = _mm256_setzero_ps();
                                    gDmrs_im = _mm256_setzero_ps();
                                    hAvg_re = _mm256_setzero_ps();
                                    hAvg_im = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        rDmrs_re_f[i] = recDmrs_rePtr[i];
                                        rDmrs_im_f[i] = recDmrs_imPtr[i];
                                        gDmrs_re_f[i] = genDmrs_rePtr[i];
                                        gDmrs_im_f[i] = genDmrs_imPtr[i];
                                        hAvg_re_f[i] = chEstAvg_rePtr[i];
                                        hAvg_im_f[i] = chEstAvg_imPtr[i];
                                    }
                                }
                                else
                                {
                                    rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                                    rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);
                                    gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                    gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                    hAvg_re  = _mm256_loadu_ps(chEstAvg_rePtr);
                                    hAvg_im  = _mm256_loadu_ps(chEstAvg_imPtr);
                                }

                                bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);
                                noise_re = _mm256_sub_ps(rDmrs_re, realOut);
                                noise_im = _mm256_sub_ps(rDmrs_im, imagOut);

                                n_avg = _mm256_add_ps(n_avg, _mm256_mul_ps(noise_re, noise_re));
                                n_avg = _mm256_add_ps(n_avg, _mm256_mul_ps(noise_im, noise_im));

                                recDmrs_rePtr += 8;
                                recDmrs_imPtr += 8;
                                genDmrs_rePtr += 8;
                                genDmrs_imPtr += 8;
                                chEstAvg_rePtr += 8;
                                chEstAvg_imPtr += 8;
                            }
                        }
                    }
                }
                else if(pNrPuschInParams->nPortIndex[idxLyr] == 1)
                {
                    for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                    {
                        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                        {
                            startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                            endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);
                            recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                            recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][(MAX_NUM_SC/2) + startScIdx];
                            chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][startScIdx];
                            chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][startScIdx];
                            genDmrs_rePtr = &pNrPuschOutParams->genDMRS[0][startScIdx];
                            genDmrs_imPtr = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];

                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                if( endScIdx - idxSc < 8)
                                {
                                    remSc = (endScIdx) - idxSc;

                                    rDmrs_re = _mm256_setzero_ps();
                                    rDmrs_im = _mm256_setzero_ps();
                                    gDmrs_re = _mm256_setzero_ps();
                                    gDmrs_im = _mm256_setzero_ps();
                                    hAvg_re = _mm256_setzero_ps();
                                    hAvg_im = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        rDmrs_re_f[i] = recDmrs_rePtr[i];
                                        rDmrs_im_f[i] = recDmrs_imPtr[i];
                                        gDmrs_re_f[i] = alternateSigns_f[i]*genDmrs_rePtr[i];
                                        gDmrs_im_f[i] = alternateSigns_f[i]*genDmrs_imPtr[i];
                                        hAvg_re_f[i] = chEstAvg_rePtr[i];
                                        hAvg_im_f[i] = chEstAvg_imPtr[i];
                                    }
                                }
                                else
                                {
                                    rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                                    rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);
                                    gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                    gDmrs_re = _mm256_mul_ps(gDmrs_re, alternateSigns);

                                    gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                    gDmrs_im = _mm256_mul_ps(gDmrs_im, alternateSigns);

                                    hAvg_re  = _mm256_loadu_ps(chEstAvg_rePtr);
                                    hAvg_im  = _mm256_loadu_ps(chEstAvg_imPtr);
                                }

                                bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);
                                noise_re = _mm256_sub_ps(rDmrs_re, realOut);
                                noise_im = _mm256_sub_ps(rDmrs_im, imagOut);

                                n_avg = _mm256_add_ps(n_avg, _mm256_mul_ps(noise_re, noise_re));
                                n_avg = _mm256_add_ps(n_avg, _mm256_mul_ps(noise_im, noise_im));

                                recDmrs_rePtr += 8;
                                recDmrs_imPtr += 8;
                                genDmrs_rePtr += 8;
                                genDmrs_imPtr += 8;
                                chEstAvg_rePtr += 8;
                                chEstAvg_imPtr += 8;
                            }
                        }
                    }
                }
            }

            if( (cdm1LyrIdx[1] - cdm1LyrIdx[0]) == 1)// FDM with one layer from CDM Group 1
            {
                // do 
                // 1) NoiseRe + j*NoiseIm = Y-Havg*X; Where, Havg = mean(H(1:12));
                // 2) NV = (NoiseRe)^2 + (NoiseIm)^2;
                idxLyr = cdm1LyrIdx[0];
                if(pNrPuschInParams->nPortIndex[idxLyr] == 2)
                {
                    for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                    {
                        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                        {
                            startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                            endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);
                            recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC +startScIdx];
                            recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                            chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][startScIdx];
                            chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][startScIdx];
                            genDmrs_rePtr = &pNrPuschOutParams->genDMRS[0][startScIdx];
                            genDmrs_imPtr = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];

                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                if( endScIdx - idxSc < 8)
                                {
                                    remSc = (endScIdx) - idxSc;

                                    rDmrs_re = _mm256_setzero_ps();
                                    rDmrs_im = _mm256_setzero_ps();
                                    gDmrs_re = _mm256_setzero_ps();
                                    gDmrs_im = _mm256_setzero_ps();
                                    hAvg_re = _mm256_setzero_ps();
                                    hAvg_im = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        rDmrs_re_f[i] = recDmrs_rePtr[i];
                                        rDmrs_im_f[i] = recDmrs_imPtr[i];
                                        gDmrs_re_f[i] = genDmrs_rePtr[i];
                                        gDmrs_im_f[i] = genDmrs_imPtr[i];
                                        hAvg_re_f[i] = chEstAvg_rePtr[i];
                                        hAvg_im_f[i] = chEstAvg_imPtr[i];
                                    }
                                }
                                else
                                {
                                    rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                                    rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);
                                    gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                    gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                    hAvg_re  = _mm256_loadu_ps(chEstAvg_rePtr);
                                    hAvg_im  = _mm256_loadu_ps(chEstAvg_imPtr);
                                }


                                bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);
                                noise_re = _mm256_sub_ps(rDmrs_re, realOut);
                                noise_im = _mm256_sub_ps(rDmrs_im, imagOut);                                

                                n_avg = _mm256_add_ps(n_avg, _mm256_mul_ps(noise_re, noise_re));
                                n_avg = _mm256_add_ps(n_avg, _mm256_mul_ps(noise_im, noise_im));

                                recDmrs_rePtr += 8;
                                recDmrs_imPtr += 8;
                                genDmrs_rePtr += 8;
                                genDmrs_imPtr += 8;
                                chEstAvg_rePtr += 8;
                                chEstAvg_imPtr += 8;
                            }
                        }
                    }
                }
                else if(pNrPuschInParams->nPortIndex[idxLyr] == 3)
                {
                    for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                    {
                        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                        {
                            startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                            endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);
                            recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC +startScIdx];
                            recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                            chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][startScIdx];
                            chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][startScIdx];
                            genDmrs_rePtr = &pNrPuschOutParams->genDMRS[0][startScIdx];
                            genDmrs_imPtr = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];

                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                if( endScIdx - idxSc < 8)
                                {
                                    remSc = (endScIdx) - idxSc;

                                    rDmrs_re = _mm256_setzero_ps();
                                    rDmrs_im = _mm256_setzero_ps();
                                    gDmrs_re = _mm256_setzero_ps();
                                    gDmrs_im = _mm256_setzero_ps();
                                    hAvg_re = _mm256_setzero_ps();
                                    hAvg_im = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        rDmrs_re_f[i] = recDmrs_rePtr[i];
                                        rDmrs_im_f[i] = recDmrs_imPtr[i];
                                        gDmrs_re_f[i] = alternateSigns_f[i]*genDmrs_rePtr[i];
                                        gDmrs_im_f[i] = alternateSigns_f[i]*genDmrs_imPtr[i];
                                        hAvg_re_f[i] = chEstAvg_rePtr[i];
                                        hAvg_im_f[i] = chEstAvg_imPtr[i];
                                    }
                                }
                                else
                                {
                                    rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                                    rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);
                                    gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                    gDmrs_re = _mm256_mul_ps(gDmrs_re, alternateSigns);
                                    gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                    gDmrs_im = _mm256_mul_ps(gDmrs_im, alternateSigns);
                                    hAvg_re  = _mm256_loadu_ps(chEstAvg_rePtr);
                                    hAvg_im  = _mm256_loadu_ps(chEstAvg_imPtr);
                                }

                                bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);
                                noise_re = _mm256_sub_ps(rDmrs_re, realOut);
                                noise_im = _mm256_sub_ps(rDmrs_im, imagOut);

                                n_avg = _mm256_add_ps(n_avg, _mm256_mul_ps(noise_re, noise_re));
                                n_avg = _mm256_add_ps(n_avg, _mm256_mul_ps(noise_im, noise_im));

                                recDmrs_rePtr += 8;
                                recDmrs_imPtr += 8;
                                genDmrs_rePtr += 8;
                                genDmrs_imPtr += 8;
                                chEstAvg_rePtr += 8;
                                chEstAvg_imPtr += 8;
                            }
                        }
                    }
                }
            }
            
            NVTemp = _mm256_permute2f128_ps (n_avg, n_avg, 1);
            nVar_alloc_vec = _mm256_add_ps(n_avg, NVTemp);
            nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);
            nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);

            pNrPuschOutParams->nVar_overAll = nVar_alloc_vec_f[0] / (nPrb_total * 6 * nRxAnt*2);
            // printf("FDM nVar_overAll: %f\n", pNrPuschOutParams->nVar_overAll);
        }
        
        if(pNrPuschInParams->nNrOfLayers == 4)
        {
            // do 
            // 1) NoiseRe + j*NoiseIm = Y- ( (Havg0*X0) + (Havg1*X1) ); 
            //                                          Where, Havgi = mean(H(1:12)) of i^th CDM;
            //                                                    Xi = DMRS  from i^th CDM;
            // 2) NV = (NoiseRe)^2 + (NoiseIm)^2;
            
            n_avg = _mm256_setzero_ps();
            nVar_alloc_vec = _mm256_setzero_ps();
            for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
            {
                for(idxLyr = cdm0LyrIdx[0];idxLyr<cdm0LyrIdx[1];idxLyr++)
                {
                    if(pNrPuschInParams->nPortIndex[idxLyr] == 0)
                    {
                        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                        {
                            startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                            endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                            chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][startScIdx];
                            chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][startScIdx];
                            genDmrs_rePtr = &pNrPuschOutParams->genDMRS[0][startScIdx];
                            genDmrs_imPtr = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                            HX_rePtr = &HX_re[startScIdx];
                            HX_imPtr = &HX_im[startScIdx];

                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                if(endScIdx - startScIdx<8)
                                {
                                    remSc = (endScIdx) - idxSc;

                                    gDmrs_re = _mm256_setzero_ps();
                                    gDmrs_im = _mm256_setzero_ps();
                                    hAvg_re = _mm256_setzero_ps();
                                    hAvg_im = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        gDmrs_re_f[i] = genDmrs_rePtr[i];
                                        gDmrs_im_f[i] = genDmrs_imPtr[i];
                                        hAvg_re_f[i] = chEstAvg_rePtr[i];
                                        hAvg_im_f[i] = chEstAvg_imPtr[i];
                                    }
                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    for(i = 0; i < remSc; i++)
                                    {
                                        HX_rePtr[i] = realOut_f[i];
                                        HX_imPtr[i] = imagOut_f[i];
                                    }
                                }
                                else
                                {
                                    gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                    gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                    hAvg_re  = _mm256_loadu_ps(chEstAvg_rePtr);
                                    hAvg_im  = _mm256_loadu_ps(chEstAvg_imPtr);

                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    _mm256_storeu_ps(HX_rePtr, realOut);
                                    _mm256_storeu_ps(HX_imPtr, imagOut);

                                    genDmrs_rePtr += 8;
                                    genDmrs_imPtr += 8;
                                    chEstAvg_rePtr +=8;
                                    chEstAvg_imPtr +=8;
                                    HX_rePtr += 8;
                                    HX_imPtr += 8;
                                }
                            }
                        }
                    }
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 1)
                    {
                        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                        {
                            startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                            endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                            chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][startScIdx];
                            chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][startScIdx];
                            genDmrs_rePtr = &pNrPuschOutParams->genDMRS[0][startScIdx];
                            genDmrs_imPtr = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                            HX_rePtr = &HX_re[startScIdx];
                            HX_imPtr = &HX_im[startScIdx];


                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                if(endScIdx - startScIdx<8)
                                {
                                    remSc = (endScIdx) - idxSc;

                                    gDmrs_re = _mm256_setzero_ps();
                                    gDmrs_im = _mm256_setzero_ps();
                                    hAvg_re = _mm256_setzero_ps();
                                    hAvg_im = _mm256_setzero_ps();
                                    HX_AVX_re = _mm256_setzero_ps();
                                    HX_AVX_im = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        gDmrs_re_f[i] = alternateSigns_f[i]*genDmrs_rePtr[i];
                                        gDmrs_im_f[i] = alternateSigns_f[i]*genDmrs_imPtr[i];
                                        hAvg_re_f[i] = chEstAvg_rePtr[i];
                                        hAvg_im_f[i] = chEstAvg_imPtr[i];
                                    }
                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    for(i = 0; i < remSc; i++)
                                    {
                                        HX_rePtr[i] = HX_rePtr[i]+realOut_f[i];
                                        HX_imPtr[i] = HX_imPtr[i]+imagOut_f[i];
                                    }
                                }
                                else
                                {
                                    gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                    gDmrs_re = _mm256_mul_ps(gDmrs_re, alternateSigns);
                                    gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                    gDmrs_im = _mm256_mul_ps(gDmrs_im, alternateSigns);
                                    hAvg_re  = _mm256_loadu_ps(chEstAvg_rePtr);
                                    hAvg_im  = _mm256_loadu_ps(chEstAvg_imPtr);
                                    HX_AVX_re = _mm256_loadu_ps(HX_rePtr);
                                    HX_AVX_im = _mm256_loadu_ps(HX_imPtr);

                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);
                                    _mm256_storeu_ps(HX_rePtr, _mm256_add_ps(HX_AVX_re,realOut));
                                    _mm256_storeu_ps(HX_imPtr, _mm256_add_ps(HX_AVX_im,imagOut));

                                    genDmrs_rePtr += 8;
                                    genDmrs_imPtr += 8;
                                    chEstAvg_rePtr +=8;
                                    chEstAvg_imPtr +=8;
                                    HX_rePtr += 8;
                                    HX_imPtr += 8;
                                }
                            }
                        }
                    }
                }

                for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                {
                    
                    startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                    endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);
                    recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                    recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC/2) + startScIdx];
                    HX_rePtr = &HX_re[startScIdx];
                    HX_imPtr = &HX_im[startScIdx];

                    for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                    {
                        if(endScIdx - startScIdx<8)
                        {
                            remSc = endScIdx - idxSc;

                            rDmrs_re = _mm256_setzero_ps();
                            rDmrs_im = _mm256_setzero_ps();

                            HX_AVX_re = _mm256_setzero_ps();
                            HX_AVX_im = _mm256_setzero_ps();

                            for(i = 0; i < remSc; i++)
                            {
                                rDmrs_re_f[i] = recDmrs_rePtr[i];
                                rDmrs_im_f[i] = recDmrs_imPtr[i];
                                HX_AVX_re_f[i] = HX_rePtr[i];
                                HX_AVX_im_f[i] = HX_imPtr[i];
                            }

                            noise_re = _mm256_sub_ps(rDmrs_re, HX_AVX_re);
                            noise_im = _mm256_sub_ps(rDmrs_im, HX_AVX_im);

                            n_avg = _mm256_fmadd_ps(noise_re, noise_re, n_avg);
                            n_avg = _mm256_fmadd_ps(noise_im, noise_im, n_avg);

                        }
                        else
                        {
                            rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                            rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);
                            HX_AVX_re = _mm256_loadu_ps(HX_rePtr);
                            HX_AVX_im = _mm256_loadu_ps(HX_imPtr);

                            noise_re = _mm256_sub_ps(rDmrs_re, HX_AVX_re);
                            noise_im = _mm256_sub_ps(rDmrs_im, HX_AVX_im);

                            n_avg = _mm256_fmadd_ps(noise_re, noise_re, n_avg);
                            n_avg = _mm256_fmadd_ps(noise_im, noise_im, n_avg);

                            recDmrs_rePtr += 8;
                            recDmrs_imPtr += 8;
                            HX_rePtr += 8;
                            HX_imPtr += 8;
                        }

                    }
                }

                // CDM 2 CODE here
                for(idxLyr = cdm1LyrIdx[0];idxLyr<cdm1LyrIdx[1];idxLyr++)
                {
                    if(pNrPuschInParams->nPortIndex[idxLyr] == 2)
                    {
                        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                        {
                            startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                            endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                            chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][startScIdx];
                            chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][startScIdx];
                            genDmrs_rePtr = &pNrPuschOutParams->genDMRS[0][startScIdx];
                            genDmrs_imPtr = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                            HX_rePtr = &HX_re[startScIdx];
                            HX_imPtr = &HX_im[startScIdx];


                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                if(endScIdx - startScIdx<8)
                                {
                                    remSc = (endScIdx) - idxSc;

                                    gDmrs_re = _mm256_setzero_ps();
                                    gDmrs_im = _mm256_setzero_ps();
                                    hAvg_re = _mm256_setzero_ps();
                                    hAvg_im = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        gDmrs_re_f[i] = genDmrs_rePtr[i];
                                        gDmrs_im_f[i] = genDmrs_imPtr[i];
                                        hAvg_re_f[i] = chEstAvg_rePtr[i];
                                        hAvg_im_f[i] = chEstAvg_imPtr[i];
                                    }
                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    for(i = 0; i < remSc; i++)
                                    {
                                        HX_rePtr[i] = realOut_f[i];
                                        HX_imPtr[i] = imagOut_f[i];
                                    }
                                }
                                else
                                {
                                    gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                    gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                    hAvg_re  = _mm256_loadu_ps(chEstAvg_rePtr);
                                    hAvg_im  = _mm256_loadu_ps(chEstAvg_imPtr);



                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    _mm256_storeu_ps(HX_rePtr, realOut);
                                    _mm256_storeu_ps(HX_imPtr, imagOut);

                                    genDmrs_rePtr += 8;
                                    genDmrs_imPtr += 8;
                                    chEstAvg_rePtr +=8;
                                    chEstAvg_imPtr +=8;
                                    HX_rePtr += 8;
                                    HX_imPtr += 8;
                                }
                            }
                        }
                    }
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 3)
                    {
                        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                        {
                            startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                            endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                            chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][startScIdx];
                            chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][startScIdx];
                            genDmrs_rePtr = &pNrPuschOutParams->genDMRS[0][startScIdx];
                            genDmrs_imPtr = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                            HX_rePtr = &HX_re[startScIdx];
                            HX_imPtr = &HX_im[startScIdx];


                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                if(endScIdx - startScIdx<8)
                                {
                                    remSc = (endScIdx) - idxSc;

                                    gDmrs_re = _mm256_setzero_ps();
                                    gDmrs_im = _mm256_setzero_ps();
                                    hAvg_re = _mm256_setzero_ps();
                                    hAvg_im = _mm256_setzero_ps();
                                    HX_AVX_re = _mm256_setzero_ps();
                                    HX_AVX_im = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        gDmrs_re_f[i] = alternateSigns_f[i]*genDmrs_rePtr[i];
                                        gDmrs_im_f[i] = alternateSigns_f[i]*genDmrs_imPtr[i];
                                        hAvg_re_f[i] = chEstAvg_rePtr[i];
                                        hAvg_im_f[i] = chEstAvg_imPtr[i];
                                    }
                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    for(i = 0; i < remSc; i++)
                                    {
                                        HX_rePtr[i] = HX_rePtr[i]+realOut_f[i];
                                        HX_imPtr[i] = HX_imPtr[i]+imagOut_f[i];
                                    }
                                }
                                else
                                {
                                    gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                    gDmrs_re = _mm256_mul_ps(gDmrs_re, alternateSigns);
                                    gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                    gDmrs_im = _mm256_mul_ps(gDmrs_im, alternateSigns);
                                    hAvg_re  = _mm256_loadu_ps(chEstAvg_rePtr);
                                    hAvg_im  = _mm256_loadu_ps(chEstAvg_imPtr);
                                    HX_AVX_re = _mm256_loadu_ps(HX_rePtr);
                                    HX_AVX_im = _mm256_loadu_ps(HX_imPtr);

                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);
                                    _mm256_storeu_ps(HX_rePtr, _mm256_add_ps(HX_AVX_re,realOut));
                                    _mm256_storeu_ps(HX_imPtr, _mm256_add_ps(HX_AVX_im,imagOut));

                                    genDmrs_rePtr += 8;
                                    genDmrs_imPtr += 8;
                                    chEstAvg_rePtr +=8;
                                    chEstAvg_imPtr +=8;
                                    HX_rePtr += 8;
                                    HX_imPtr += 8;
                                }
                            }
                        }
                    }
                }

                for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                {
                    startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                    endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);
                    recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC +startScIdx];
                    recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                    HX_rePtr = &HX_re[startScIdx];
                    HX_imPtr = &HX_im[startScIdx];

                    for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                    {
                        if(endScIdx - startScIdx<8)
                        {
                            remSc = endScIdx - idxSc;

                            rDmrs_re = _mm256_setzero_ps();
                            rDmrs_im = _mm256_setzero_ps();

                            HX_AVX_re = _mm256_setzero_ps();
                            HX_AVX_im = _mm256_setzero_ps();

                            for(i = 0; i < remSc; i++)
                            {
                                rDmrs_re_f[i] = recDmrs_rePtr[i];
                                rDmrs_im_f[i] = recDmrs_imPtr[i];
                                HX_AVX_re_f[i] = HX_rePtr[i];
                                HX_AVX_im_f[i] = HX_imPtr[i];
                            }

                            noise_re = _mm256_sub_ps(rDmrs_re, HX_AVX_re);
                            noise_im = _mm256_sub_ps(rDmrs_im, HX_AVX_im);


                            n_avg = _mm256_fmadd_ps(noise_re, noise_re, n_avg);
                            n_avg = _mm256_fmadd_ps(noise_im, noise_im, n_avg);

                        }
                        else
                        {
                            rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                            rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);
                            HX_AVX_re = _mm256_loadu_ps(HX_rePtr);
                            HX_AVX_im = _mm256_loadu_ps(HX_imPtr);

                            noise_re = _mm256_sub_ps(rDmrs_re, HX_AVX_re);
                            noise_im = _mm256_sub_ps(rDmrs_im, HX_AVX_im);

                            n_avg = _mm256_fmadd_ps(noise_re, noise_re, n_avg);
                            n_avg = _mm256_fmadd_ps(noise_im, noise_im, n_avg);

                            recDmrs_rePtr += 8;
                            recDmrs_imPtr += 8;
                            HX_rePtr += 8;
                            HX_imPtr += 8;
                        }

                    }
                }


            }

            NVTemp = _mm256_permute2f128_ps (n_avg, n_avg, 1);
            nVar_alloc_vec = _mm256_add_ps(n_avg, NVTemp);
            nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);
            nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);

            pNrPuschOutParams->nVar_overAll = nVar_alloc_vec_f[0] / (nPrb_total * 6 * nRxAnt*2);
            // printf("nVar_overAll: %f\n", pNrPuschOutParams->nVar_overAll);
        }
    }
}
//*/


wnVoid nr_pusch_mmse_perTone_Havg_avx2(P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams,
                                       puschConfigT* pNrPuschInParams)
{
    wnUInt8 nRxAnt   = pNrPuschInParams->nAntennaPorts;
    wnUInt8 nLyr     = pNrPuschInParams->nNrOfLayers;
    wnUInt8 nSbAlloc = pNrPuschInParams->nSbAlloc;
    wnUInt8 idxLyr;
    wnUInt16 idxAlloc, idxAnt, idx;
    wnUInt16 toneStartIdx, toneStartIdx2;//toneEndIdx
    wnUInt16 PRBIdx, PRBStart, PRBEnd;

    __m256 inpReal1, inpReal2, inpReal3;
    __m256 inpImag1, inpImag2, inpImag3;
    __m256 temp11, temp12, temp21, temp22;

    wnFlt *chEst_rePtr, *chEst_imPtr;
    wnFlt *chEstAvg_rePtr, *chEstAvg_imPtr;
    wnFlt *temp12_f = (wnFlt*)&temp12;
    wnFlt *temp22_f = (wnFlt*)&temp22;
    wnFlt finalValReal, finalValImag;

    __m256 oneBy12 = _mm256_set1_ps(0.083333333333333);
    //Find mean of twelve estimates and store it in six consecutive locations of output buffer.
    for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
    { 
        for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
        {
            for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
            {
                toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
                // toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);
                chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];
                PRBStart = pNrPuschInParams->rbStart[idxAlloc];
                PRBEnd = pNrPuschInParams->nPrb[idxAlloc];

                toneStartIdx2 = pNrPuschInParams->rbStart[idxAlloc] * 6;
                chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][toneStartIdx2];
                chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][toneStartIdx2];

                for(PRBIdx = PRBStart;PRBIdx<PRBEnd;PRBIdx+=2)
                {
                    if(PRBEnd - PRBIdx == 1)
                    { 
                        inpReal1 = _mm256_loadu_ps(chEst_rePtr);
                        temp12 = _mm256_hadd_ps(inpReal1, inpReal1);
                        temp12 = _mm256_permutevar8x32_ps(temp12, _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7));

                        inpImag1 = _mm256_loadu_ps(chEst_imPtr);
                        temp22 = _mm256_hadd_ps(inpImag1, inpImag1);
                        temp22 = _mm256_permutevar8x32_ps(temp22, _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7));

                        finalValReal = 0;
                        finalValImag = 0;
                        for(idx = 0;idx<4;idx++)
                        {
                            finalValReal = finalValReal+temp12_f[idx]+chEst_rePtr[idx+8];
                            finalValImag = finalValImag+temp22_f[idx]+chEst_imPtr[idx+8];
                        }
                        chEstAvg_rePtr[0] = finalValReal*0.083333333333333;//1/12 = 0.083333333333333;
                        chEstAvg_imPtr[0] = finalValImag*0.083333333333333;
                        for(idx = 1;idx<6;idx++)
                        {
                            chEstAvg_rePtr[idx] = chEstAvg_rePtr[0];
                            chEstAvg_imPtr[idx] = chEstAvg_imPtr[0];
                        }
                    }
                    else
                    { 
                        inpReal1 = _mm256_loadu_ps(chEst_rePtr);
                        inpReal2 = _mm256_loadu_ps(chEst_rePtr+8);
                        inpReal3 = _mm256_loadu_ps(chEst_rePtr+16);

                        temp11 = _mm256_hadd_ps(inpReal1, inpReal1);
                        temp11 = _mm256_permutevar8x32_ps(temp11, _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7));
                        temp12 = _mm256_hadd_ps(temp11, inpReal2);
                        temp12 = _mm256_hadd_ps(temp12, temp12);
                        temp12 = _mm256_hadd_ps(temp12, temp12);
                        temp12 = _mm256_mul_ps(temp12, oneBy12);
                        _mm256_storeu_ps(chEstAvg_rePtr, _mm256_set1_ps(temp12_f[0]));

                        temp21 = _mm256_hadd_ps(inpReal3, inpReal3);
                        temp21 = _mm256_permutevar8x32_ps(temp21, _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7));
                        temp22 = _mm256_hadd_ps(temp21, inpReal2);
                        temp22 = _mm256_hadd_ps(temp22, temp22);
                        temp22 = _mm256_hadd_ps(temp22, temp22);
                        temp22 = _mm256_mul_ps(temp22, oneBy12);
                        
                        _mm256_storeu_ps(chEstAvg_rePtr+6, _mm256_set1_ps(temp22_f[4]));

                        inpImag1 = _mm256_loadu_ps(chEst_imPtr);
                        inpImag2 = _mm256_loadu_ps(chEst_imPtr+8);
                        inpImag3 = _mm256_loadu_ps(chEst_imPtr+16);

                        temp11 = _mm256_hadd_ps(inpImag1, inpImag1);
                        temp11 = _mm256_permutevar8x32_ps(temp11, _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7));
                        temp12 = _mm256_hadd_ps(temp11, inpImag2);
                        temp12 = _mm256_hadd_ps(temp12, temp12);
                        temp12 = _mm256_hadd_ps(temp12, temp12);
                        temp12 = _mm256_mul_ps(temp12, oneBy12);
                        _mm256_storeu_ps(chEstAvg_imPtr, _mm256_set1_ps(temp12_f[0]));

                        temp21 = _mm256_hadd_ps(inpImag3, inpImag3);
                        temp21 = _mm256_permutevar8x32_ps(temp21, _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7));
                        temp22 = _mm256_hadd_ps(temp21, inpImag2);
                        temp22 = _mm256_hadd_ps(temp22, temp22);
                        temp22 = _mm256_hadd_ps(temp22, temp22);
                        temp22 = _mm256_mul_ps(temp22, oneBy12);
                        _mm256_storeu_ps(chEstAvg_imPtr+6, _mm256_set1_ps(temp22_f[4]));

                        chEstAvg_rePtr += 12;
                        chEstAvg_imPtr += 12;
                        chEst_rePtr += 24;
                        chEst_imPtr += 24;
                    }
                }
            }
        }
    }
}


wnVoid nr_pusch_mmse_perTone_nVarCalc_avx2(wnUInt8 segPart,
                               commonUlConfigT* pNrUlCommonParams,
                               puschConfigT* pNrPuschInParams,
                               P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams)
{
    wnUInt8 cdmTotal[2] = {0, 0};
    wnUInt8 idxLyr, idxSym;
    // wnUInt8 nullToneFlag = 1;
    wnUInt8 dmrsSym;
    wnUInt8 nRxAnt   = pNrPuschInParams->nAntennaPorts;
    wnUInt8 nLyr     = pNrPuschInParams->nNrOfLayers;
    wnUInt8 nSbAlloc = pNrPuschInParams->nSbAlloc;
    wnUInt8 numOfDMRSSymbs = pNrPuschInParams->nNrOfDMRSSymbols;
    
    wnUInt16 idxAlloc, idxAnt;
    wnUInt16 nullToneStart, nullToneEnd;
    wnUInt16 nullToneIdx, remSc;
    wnUInt16 idxTemp;
    wnUInt16 startScIdx, endScIdx, idxSc, i;

    wnUInt16 nPrb_total = pNrPuschOutParams->total_nPRB;

    __m256 n1_re, n1_im, NVTemp;
    __m256 nVar_alloc_vec = _mm256_set1_ps(0);
    __m256 n_avg = _mm256_set1_ps(0);
    __m256 alternateSigns = _mm256_setr_ps(1, -1, 1, -1, 1, -1, 1, -1);
    __m256 rDmrs_re, rDmrs_im, gDmrs_re, gDmrs_im;
    __m256 hAvg_re, hAvg_im;
    __m256 bd_re, ad_im;
    __m256 realOut, imagOut;
    __m256 noise_re, noise_im;
    __m256 HX_AVX_re, HX_AVX_im;
    __m256 neg = _mm256_set1_ps(-0.0);
    
    
    wnFlt *nVar_alloc_vec_f = (wnFlt *)&nVar_alloc_vec;
    wnFlt *n_avg_f = (wnFlt*)&n_avg;
    wnFlt *alternateSigns_f = (wnFlt *)&alternateSigns;
    wnFlt *chEstAvg_rePtr, *chEstAvg_imPtr;
    wnFlt *recDmrs_rePtr, *recDmrs_imPtr, *genDmrs_rePtr, *genDmrs_imPtr;
    // wnFlt *noise_re_f = (wnFlt*)&noise_re;
    // wnFlt *noise_im_f = (wnFlt*)&noise_im;

    wnFlt *realOut_f = (wnFlt*)&realOut;
    wnFlt *imagOut_f = (wnFlt*)&imagOut;
    wnFlt *rDmrs_re_f = (wnFlt*)&rDmrs_re;
    wnFlt *rDmrs_im_f = (wnFlt*)&rDmrs_im;
    wnFlt *gDmrs_re_f = (wnFlt*)&gDmrs_re;
    wnFlt *gDmrs_im_f = (wnFlt*)&gDmrs_im;
    wnFlt *hAvg_re_f = (wnFlt*)&hAvg_re;
    wnFlt *hAvg_im_f = (wnFlt*)&hAvg_im;
    wnFlt HX_re[6 * MAX_NUM_PRB];
    wnFlt HX_im[6 * MAX_NUM_PRB];
    wnFlt *HX_AVX_re_f = (wnFlt*)&HX_AVX_re;
    wnFlt *HX_AVX_im_f = (wnFlt*)&HX_AVX_im;
    wnFlt *HX_rePtr;
    wnFlt *HX_imPtr;

    //cdmTotal[0] => number of antenna ports used from CDM Group 0.
    //cdmTotal[1] => number of antenna ports used from CDM Group 1.
    for(idxLyr = 0;idxLyr<pNrPuschInParams->nNrOfLayers;idxLyr++)
    {
        if(pNrPuschInParams->nPortIndex[idxLyr] == 0 || pNrPuschInParams->nPortIndex[idxLyr] == 1 ||
           pNrPuschInParams->nPortIndex[idxLyr] == 4 || pNrPuschInParams->nPortIndex[idxLyr] == 5)
        {
            cdmTotal[0] = cdmTotal[0]+1;
        }
        else
        {
            cdmTotal[1] = cdmTotal[1]+1;
        }
    }

    // Get DMRS Symbol Index.
    if(segPart == SEG1)
    {
        dmrsSym = pNrPuschInParams->dmrs1_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG2)
    {
        dmrsSym = pNrPuschInParams->dmrs2_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG3)
    {
        dmrsSym = pNrPuschInParams->dmrs3_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else //if(segPart == SEG4)
    {
        dmrsSym = pNrPuschInParams->dmrs4_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    

    //Get ready Havg for Noise Variance(NV) calculations.
    //Havg = mean of estimates from one PRB.
    nr_pusch_mmse_perTone_Havg_avx2(pNrPuschOutParams, pNrPuschInParams);

    #ifdef DEBUG
        FILE *fpReal, *fpImag;
        if(segPart == SEG1)
        {
            fpReal = fopen("IOs/HAvgRealSeg1_mmse_PerTone.txt", "w");
            fpImag = fopen("IOs/HAvgImagSeg1_mmse_PerTone.txt", "w");
        }
        else if(segPart == SEG2)
        {
            fpReal = fopen("IOs/HAvgRealSeg2_mmse_PerTone.txt", "w");
            fpImag = fopen("IOs/HAvgImagSeg2_mmse_PerTone.txt", "w");
        }
        else if(segPart == SEG3)
        {
            fpReal = fopen("IOs/HAvgRealSeg3_mmse_PerTone.txt", "w");
            fpImag = fopen("IOs/HAvgImagSeg3_mmse_PerTone.txt", "w");
        }
        else //if(segPart == SEG4)
        {
            fpReal = fopen("IOs/HAvgRealSeg4_mmse_PerTone.txt", "w");
            fpImag = fopen("IOs/HAvgImagSeg4_mmse_PerTone.txt", "w");
        }

        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            for(int idxLyr = 0; idxLyr < nLyr; idxLyr++)
            {
                for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                {
                    for(int toneIdx = 0;toneIdx<1650;toneIdx++)
                    {
                        fprintf(fpReal, "%f\n",\
                        pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][toneIdx]);
                        fprintf(fpImag, "%f\n",\
                        pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][toneIdx]);
                    }
                }
            }
        }

        fclose(fpReal);
        fclose(fpImag);
    #endif

    if(cdmTotal[1] == 0)// if CDM Group 1 is Empty.
    {
        //Received signal at null tones are considered as noise samples
        // and NV is calculated using them.
        //Do NV = (NoiseReal).^2 + (NoiseImag).^2;
        nVar_alloc_vec = _mm256_setzero_ps();
        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            for(idxSym = 0;idxSym<numOfDMRSSymbs;idxSym++)
            {
                for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                {
                    n_avg = _mm256_setzero_ps();
                    nullToneStart = MAX_NUM_SC;
                    nullToneEnd = nullToneStart+ (pNrPuschInParams->nPrb[idxAlloc] * 6);

                    for(nullToneIdx = nullToneStart;nullToneIdx<nullToneEnd;nullToneIdx += 8)
                    {
                        remSc = nullToneEnd - nullToneIdx;
                        if(remSc < 8)
                        {
                            for(idxTemp = 0;idxTemp<remSc;idxTemp++)
                            {
                                n_avg_f[idxTemp] = n_avg_f[idxTemp]+ (rxFdSamples[idxAnt][dmrsSym+idxSym][nullToneIdx+idxTemp]*\
                                                                      rxFdSamples[idxAnt][dmrsSym+idxSym][nullToneIdx+idxTemp]);
                            }
                        }
                        else
                        {
                            n1_re = _mm256_loadu_ps(&rxFdSamples[idxAnt][dmrsSym+idxSym][nullToneIdx]);
                            n_avg = _mm256_fmadd_ps(n1_re, n1_re, n_avg);
                        }
                    }

                    nullToneStart = (MAX_NUM_SC)+(MAX_NUM_SC/2);
                    nullToneEnd = nullToneStart+ (pNrPuschInParams->nPrb[idxAlloc] * 6);
                    for(nullToneIdx = nullToneStart;nullToneIdx<nullToneEnd;nullToneIdx += 8)
                    {
                        remSc = nullToneEnd - nullToneIdx;
                        if(remSc < 8)
                        {
                            for(idxTemp = 0;idxTemp<remSc;idxTemp++)
                            {
                                n_avg_f[idxTemp] = n_avg_f[idxTemp]+ (rxFdSamples[idxAnt][dmrsSym+idxSym][nullToneIdx+idxTemp]*\
                                                                      rxFdSamples[idxAnt][dmrsSym+idxSym][nullToneIdx+idxTemp]);
                            }
                        }
                        else
                        {
                            n1_im = _mm256_loadu_ps(&rxFdSamples[idxAnt][dmrsSym+idxSym][nullToneIdx]);
                            n_avg = _mm256_fmadd_ps(n1_im, n1_im, n_avg);
                        }
                    }
                    nVar_alloc_vec = _mm256_add_ps(nVar_alloc_vec, n_avg);
                }
            }
        }

        //Below four lines of code, will do sum of all elements in AVX vector "nVar_alloc_vec" and 
        //its final value is positioned in zeroth position of "nVar_alloc_vec" vector.
        NVTemp = _mm256_permute2f128_ps (nVar_alloc_vec, nVar_alloc_vec, 1);
        nVar_alloc_vec = _mm256_add_ps(nVar_alloc_vec, NVTemp);
        nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);
        nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);

        pNrPuschOutParams->nVar_overAll = nVar_alloc_vec_f[0] / (nPrb_total * 6 * nRxAnt * numOfDMRSSymbs);
    }
    else if(cdmTotal[0] == 0)// if CDM Group 0 is Empty.
    {
        //Received signal at null tones are considered as noise samples and
        // NV is calculated using them.
        //Do NV = (NoiseReal).^2 + (NoiseImag).^2;
        nVar_alloc_vec = _mm256_setzero_ps();
        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            for(idxSym = 0;idxSym<numOfDMRSSymbs;idxSym++)
            {
                for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                {
                    nullToneStart = 0;
                    n_avg = _mm256_setzero_ps();
                    nullToneEnd = nullToneStart+ (pNrPuschInParams->nPrb[idxAlloc] * 6);

                    for(nullToneIdx = nullToneStart;nullToneIdx<nullToneEnd;nullToneIdx += 8)
                    {
                        remSc = nullToneEnd - nullToneIdx;
                        if(remSc < 8)
                        {
                            for(idxTemp = 0;idxTemp<remSc;idxTemp++)
                            {
                                n_avg_f[idxTemp] = n_avg_f[idxTemp]+ (rxFdSamples[idxAnt][dmrsSym+idxSym][nullToneIdx+idxTemp]*\
                                                                      rxFdSamples[idxAnt][dmrsSym+idxSym][nullToneIdx+idxTemp]);
                            }
                        }
                        else
                        {
                            n1_re = _mm256_loadu_ps(&rxFdSamples[idxAnt][dmrsSym+idxSym][nullToneIdx]);
                            n_avg = _mm256_fmadd_ps(n1_re, n1_re, n_avg);
                        }
                    }
                    

                    nullToneStart = (MAX_NUM_SC/2);
                    nullToneEnd = nullToneStart+ (pNrPuschInParams->nPrb[idxAlloc] * 6);
                    for(nullToneIdx = nullToneStart;nullToneIdx<nullToneEnd;nullToneIdx += 8)
                    {
                        remSc = nullToneEnd - nullToneIdx;
                        if(remSc < 8)
                        {
                            for(idxTemp = 0;idxTemp<remSc;idxTemp++)
                            {
                                n_avg_f[idxTemp] = n_avg_f[idxTemp]+ (rxFdSamples[idxAnt][dmrsSym+idxSym][nullToneIdx+idxTemp]*\
                                                                      rxFdSamples[idxAnt][dmrsSym+idxSym][nullToneIdx+idxTemp]);
                            }
                        }
                        else
                        {
                            n1_im = _mm256_loadu_ps(&rxFdSamples[idxAnt][dmrsSym+idxSym][nullToneIdx]);
                            n_avg = _mm256_fmadd_ps(n1_im, n1_im, n_avg);
                        }
                    }
                    nVar_alloc_vec = _mm256_add_ps(nVar_alloc_vec, n_avg);
                }
            }
        }

        //Below four lines of code, will do sum of all elements in AVX vector "nVar_alloc_vec" and 
        //its final value is positioned in zeroth position of "nVar_alloc_vec" vector.
        NVTemp = _mm256_permute2f128_ps (nVar_alloc_vec, nVar_alloc_vec, 1);
        nVar_alloc_vec = _mm256_add_ps(nVar_alloc_vec, NVTemp);
        nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);
        nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);

        pNrPuschOutParams->nVar_overAll =  nVar_alloc_vec_f[0] / (nPrb_total * 6 * nRxAnt * numOfDMRSSymbs);
    }
    else
    {
        // Noise Variance estimates will be obtained from two cdm groups and
        // two estimates are averaged.

        //Do   i) NoiseReal_0+NoiseImag_0*1i = Y - (Havg(0)*Xdmrs(0)+ ... + Havg(n)*Xdmrs(n));n = Number of ports from CDM Group 0.
        //    ii) n_avg_0 = (NoiseReal_0).^2 + (NoiseImag_0).^2;//noise variance estimate from CDM Group 0.
        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        { 
            for(idxSym = 0;idxSym<numOfDMRSSymbs;idxSym++)
            {
                memset(HX_re, 0, 4*6*MAX_NUM_PRB);//set 4bytes*nLocations to zero.
                memset(HX_im, 0, 4*6*MAX_NUM_PRB);//set 4bytes*nLocations to zero.

                for(idxLyr = 0;idxLyr<nLyr;idxLyr++)
                {
                    if(pNrPuschInParams->nPortIndex[idxLyr] == 0)
                    {
                        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                        {
                            startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                            endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                            chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][startScIdx];
                            chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][startScIdx];
                            genDmrs_rePtr = &pNrPuschOutParams->genDMRS[idxSym][startScIdx];
                            genDmrs_imPtr = &pNrPuschOutParams->genDMRS[idxSym][(MAX_NUM_SC/2) + startScIdx];
                            HX_rePtr = &HX_re[startScIdx];
                            HX_imPtr = &HX_im[startScIdx];

                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                if( (endScIdx - startScIdx) <8)
                                {
                                    remSc = endScIdx - idxSc;

                                    gDmrs_re = _mm256_setzero_ps();
                                    gDmrs_im = _mm256_setzero_ps();
                                    hAvg_re = _mm256_setzero_ps();
                                    hAvg_im = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        gDmrs_re_f[i] = genDmrs_rePtr[i];
                                        gDmrs_im_f[i] = genDmrs_imPtr[i];
                                        hAvg_re_f[i] = chEstAvg_rePtr[i];
                                        hAvg_im_f[i] = chEstAvg_imPtr[i];
                                    }
                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    for(i = 0; i < remSc; i++)
                                    {
                                        HX_rePtr[i] = realOut_f[i];
                                        HX_imPtr[i] = imagOut_f[i];
                                    }
                                }
                                else
                                {
                                    gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                    gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                    hAvg_re  = _mm256_loadu_ps(chEstAvg_rePtr);
                                    hAvg_im  = _mm256_loadu_ps(chEstAvg_imPtr);

                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    _mm256_storeu_ps(HX_rePtr, realOut);
                                    _mm256_storeu_ps(HX_imPtr, imagOut);

                                    genDmrs_rePtr += 8;
                                    genDmrs_imPtr += 8;
                                    chEstAvg_rePtr +=8;
                                    chEstAvg_imPtr +=8;
                                    HX_rePtr += 8;
                                    HX_imPtr += 8;
                                }
                            }
                        }
                    }
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 1 )
                    {
                        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                        {
                            startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                            endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                            chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][startScIdx];
                            chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][startScIdx];
                            genDmrs_rePtr = &pNrPuschOutParams->genDMRS[idxSym][startScIdx];
                            genDmrs_imPtr = &pNrPuschOutParams->genDMRS[idxSym][(MAX_NUM_SC/2) + startScIdx];
                            HX_rePtr = &HX_re[startScIdx];
                            HX_imPtr = &HX_im[startScIdx];


                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                if(endScIdx - startScIdx<8)
                                {
                                    remSc = (endScIdx) - idxSc;

                                    gDmrs_re = _mm256_setzero_ps();
                                    gDmrs_im = _mm256_setzero_ps();
                                    hAvg_re = _mm256_setzero_ps();
                                    hAvg_im = _mm256_setzero_ps();
                                    HX_AVX_re = _mm256_setzero_ps();
                                    HX_AVX_im = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        gDmrs_re_f[i] = alternateSigns_f[i]*genDmrs_rePtr[i];
                                        gDmrs_im_f[i] = alternateSigns_f[i]*genDmrs_imPtr[i];
                                        hAvg_re_f[i] = chEstAvg_rePtr[i];
                                        hAvg_im_f[i] = chEstAvg_imPtr[i];
                                    }
                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    for(i = 0; i < remSc; i++)
                                    {
                                        HX_rePtr[i] = HX_rePtr[i]+realOut_f[i];
                                        HX_imPtr[i] = HX_imPtr[i]+imagOut_f[i];
                                    }
                                }
                                else
                                {
                                    gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                    gDmrs_re = _mm256_mul_ps(gDmrs_re, alternateSigns);
                                    gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                    gDmrs_im = _mm256_mul_ps(gDmrs_im, alternateSigns);
                                    hAvg_re  = _mm256_loadu_ps(chEstAvg_rePtr);
                                    hAvg_im  = _mm256_loadu_ps(chEstAvg_imPtr);
                                    HX_AVX_re = _mm256_loadu_ps(HX_rePtr);
                                    HX_AVX_im = _mm256_loadu_ps(HX_imPtr);

                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);
                                    _mm256_storeu_ps(HX_rePtr, _mm256_add_ps(HX_AVX_re,realOut));
                                    _mm256_storeu_ps(HX_imPtr, _mm256_add_ps(HX_AVX_im,imagOut));

                                    genDmrs_rePtr += 8;
                                    genDmrs_imPtr += 8;
                                    chEstAvg_rePtr +=8;
                                    chEstAvg_imPtr +=8;
                                    HX_rePtr += 8;
                                    HX_imPtr += 8;
                                }
                            }
                        }
                    }
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 4)
                    {
                        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                        {
                            startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                            endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                            chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][startScIdx];
                            chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][startScIdx];
                            genDmrs_rePtr = &pNrPuschOutParams->genDMRS[idxSym][startScIdx];
                            genDmrs_imPtr = &pNrPuschOutParams->genDMRS[idxSym][(MAX_NUM_SC/2) + startScIdx];
                            HX_rePtr = &HX_re[startScIdx];
                            HX_imPtr = &HX_im[startScIdx];

                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                if(endScIdx - startScIdx<8)
                                {
                                    remSc = (endScIdx) - idxSc;

                                    gDmrs_re = _mm256_setzero_ps();
                                    gDmrs_im = _mm256_setzero_ps();
                                    hAvg_re = _mm256_setzero_ps();
                                    hAvg_im = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        gDmrs_re_f[i] = genDmrs_rePtr[i];
                                        gDmrs_im_f[i] = genDmrs_imPtr[i];
                                        hAvg_re_f[i] = chEstAvg_rePtr[i];
                                        hAvg_im_f[i] = chEstAvg_imPtr[i];
                                    }
                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);
                                    
                                    if(idxSym == 1)
                                    {
                                        //realOut = -1*realOut;
                                        realOut = _mm256_xor_ps(realOut, neg);
                                        //imagOut = -1*imagOut;
                                        imagOut = _mm256_xor_ps(imagOut, neg);
                                    }

                                    for(i = 0; i < remSc; i++)
                                    {
                                        HX_rePtr[i] = HX_rePtr[i] + realOut_f[i];
                                        HX_imPtr[i] = HX_imPtr[i] + imagOut_f[i];
                                    }
                                }
                                else
                                {
                                    gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                    gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                    hAvg_re  = _mm256_loadu_ps(chEstAvg_rePtr);
                                    hAvg_im  = _mm256_loadu_ps(chEstAvg_imPtr);
                                    HX_AVX_re = _mm256_loadu_ps(HX_rePtr);
                                    HX_AVX_im = _mm256_loadu_ps(HX_imPtr);

                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    if(idxSym == 1)
                                    {
                                        //realOut = -1*realOut;
                                        realOut = _mm256_xor_ps(realOut, neg);
                                        //imagOut = -1*imagOut;
                                        imagOut = _mm256_xor_ps(imagOut, neg);
                                    }

                                    _mm256_storeu_ps(HX_rePtr, _mm256_add_ps(HX_AVX_re,realOut));
                                    _mm256_storeu_ps(HX_imPtr, _mm256_add_ps(HX_AVX_im,imagOut));

                                    genDmrs_rePtr += 8;
                                    genDmrs_imPtr += 8;
                                    chEstAvg_rePtr +=8;
                                    chEstAvg_imPtr +=8;
                                    HX_rePtr += 8;
                                    HX_imPtr += 8;
                                }
                            }
                        }
                    
                    }
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 5)
                    {
                        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                        {
                            startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                            endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                            chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][startScIdx];
                            chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][startScIdx];
                            genDmrs_rePtr = &pNrPuschOutParams->genDMRS[idxSym][startScIdx];
                            genDmrs_imPtr = &pNrPuschOutParams->genDMRS[idxSym][(MAX_NUM_SC/2) + startScIdx];
                            HX_rePtr = &HX_re[startScIdx];
                            HX_imPtr = &HX_im[startScIdx];


                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                if(endScIdx - startScIdx<8)
                                {
                                    remSc = (endScIdx) - idxSc;

                                    gDmrs_re = _mm256_setzero_ps();
                                    gDmrs_im = _mm256_setzero_ps();
                                    hAvg_re = _mm256_setzero_ps();
                                    hAvg_im = _mm256_setzero_ps();
                                    HX_AVX_re = _mm256_setzero_ps();
                                    HX_AVX_im = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        gDmrs_re_f[i] = alternateSigns_f[i]*genDmrs_rePtr[i];
                                        gDmrs_im_f[i] = alternateSigns_f[i]*genDmrs_imPtr[i];
                                        hAvg_re_f[i] = chEstAvg_rePtr[i];
                                        hAvg_im_f[i] = chEstAvg_imPtr[i];
                                    }
                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    if(idxSym == 1)
                                    {
                                        //realOut = -1*realOut;
                                        realOut = _mm256_xor_ps(realOut, neg);
                                        //imagOut = -1*imagOut;
                                        imagOut = _mm256_xor_ps(imagOut, neg);
                                    }

                                    for(i = 0; i < remSc; i++)
                                    {
                                        HX_rePtr[i] = HX_rePtr[i]+realOut_f[i];
                                        HX_imPtr[i] = HX_imPtr[i]+imagOut_f[i];
                                    }
                                }
                                else
                                {
                                    gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                    gDmrs_re = _mm256_mul_ps(gDmrs_re, alternateSigns);
                                    gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                    gDmrs_im = _mm256_mul_ps(gDmrs_im, alternateSigns);
                                    hAvg_re  = _mm256_loadu_ps(chEstAvg_rePtr);
                                    hAvg_im  = _mm256_loadu_ps(chEstAvg_imPtr);
                                    HX_AVX_re = _mm256_loadu_ps(HX_rePtr);
                                    HX_AVX_im = _mm256_loadu_ps(HX_imPtr);

                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    if(idxSym == 1)
                                    {
                                        //realOut = -1*realOut;
                                        realOut = _mm256_xor_ps(realOut, neg);
                                        //imagOut = -1*imagOut;
                                        imagOut = _mm256_xor_ps(imagOut, neg);
                                    }

                                    _mm256_storeu_ps(HX_rePtr, _mm256_add_ps(HX_AVX_re,realOut));
                                    _mm256_storeu_ps(HX_imPtr, _mm256_add_ps(HX_AVX_im,imagOut));

                                    genDmrs_rePtr += 8;
                                    genDmrs_imPtr += 8;
                                    chEstAvg_rePtr +=8;
                                    chEstAvg_imPtr +=8;
                                    HX_rePtr += 8;
                                    HX_imPtr += 8;
                                }
                            }
                        }
                    
                    }
                }

                for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                {
                    startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                    endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);
                    recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym+idxSym][startScIdx];
                    recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym+idxSym][(int)(MAX_NUM_SC/2) + startScIdx];
                    HX_rePtr = &HX_re[startScIdx];
                    HX_imPtr = &HX_im[startScIdx];

                    for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                    {
                        if((endScIdx - startScIdx)<8)
                        {
                            remSc = endScIdx - idxSc;

                            rDmrs_re = _mm256_setzero_ps();
                            rDmrs_im = _mm256_setzero_ps();

                            HX_AVX_re = _mm256_setzero_ps();
                            HX_AVX_im = _mm256_setzero_ps();

                            for(i = 0; i < remSc; i++)
                            {
                                rDmrs_re_f[i] = recDmrs_rePtr[i];
                                rDmrs_im_f[i] = recDmrs_imPtr[i];
                                HX_AVX_re_f[i] = HX_rePtr[i];
                                HX_AVX_im_f[i] = HX_imPtr[i];
                            }

                            noise_re = _mm256_sub_ps(rDmrs_re, HX_AVX_re);
                            noise_im = _mm256_sub_ps(rDmrs_im, HX_AVX_im);

                            n_avg = _mm256_fmadd_ps(noise_re, noise_re, n_avg);
                            n_avg = _mm256_fmadd_ps(noise_im, noise_im, n_avg);
                        }
                        else
                        {
                            rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                            rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);
                            HX_AVX_re = _mm256_loadu_ps(HX_rePtr);
                            HX_AVX_im = _mm256_loadu_ps(HX_imPtr);

                            noise_re = _mm256_sub_ps(rDmrs_re, HX_AVX_re);
                            noise_im = _mm256_sub_ps(rDmrs_im, HX_AVX_im);

                            n_avg = _mm256_fmadd_ps(noise_re, noise_re, n_avg);
                            n_avg = _mm256_fmadd_ps(noise_im, noise_im, n_avg);

                            recDmrs_rePtr += 8;
                            recDmrs_imPtr += 8;
                            HX_rePtr += 8;
                            HX_imPtr += 8;
                        }
                    }
                }
            }
        }

        NVTemp = _mm256_permute2f128_ps (n_avg, n_avg, 1);
        nVar_alloc_vec = _mm256_add_ps(n_avg, NVTemp);
        nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);
        nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);
        pNrPuschOutParams->nVar_overAll = nVar_alloc_vec_f[0] / (nPrb_total * 6 * nRxAnt * numOfDMRSSymbs);
        n_avg = _mm256_set1_ps(0);
        
        //Do   i) NoiseReal_1+NoiseImag_1*1i = Y - (Havg(0)*Xdmrs(0)+ ... + Havg(n)*Xdmrs(n));n = Number of ports from CDM Group 1.
        //    ii) n_avg_1 = (NoiseReal_1).^2 + (NoiseImag_1).^2;//noise variance estimate from CDM Group 1.
        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        { 
            for(idxSym = 0;idxSym<numOfDMRSSymbs;idxSym++)
            {
                memset(HX_re, 0, 4*6*MAX_NUM_PRB);//set 4bytes*nLocations to zero.
                memset(HX_im, 0, 4*6*MAX_NUM_PRB);//set 4bytes*nLocations to zero.
                for(idxLyr = 0;idxLyr<nLyr;idxLyr++)
                {
                    if(pNrPuschInParams->nPortIndex[idxLyr] == 2)
                    {
                        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                        {
                            startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                            endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                            chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][startScIdx];
                            chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][startScIdx];
                            genDmrs_rePtr = &pNrPuschOutParams->genDMRS[idxSym][startScIdx];
                            genDmrs_imPtr = &pNrPuschOutParams->genDMRS[idxSym][(MAX_NUM_SC/2) + startScIdx];
                            HX_rePtr = &HX_re[startScIdx];
                            HX_imPtr = &HX_im[startScIdx];


                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                remSc = endScIdx - idxSc;
                                if(remSc<8)
                                {
                                    gDmrs_re = _mm256_setzero_ps();
                                    gDmrs_im = _mm256_setzero_ps();
                                    hAvg_re = _mm256_setzero_ps();
                                    hAvg_im = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        gDmrs_re_f[i] = genDmrs_rePtr[i];
                                        gDmrs_im_f[i] = genDmrs_imPtr[i];
                                        hAvg_re_f[i] = chEstAvg_rePtr[i];
                                        hAvg_im_f[i] = chEstAvg_imPtr[i];
                                    }
                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    for(i = 0; i < remSc; i++)
                                    {
                                        HX_rePtr[i] = realOut_f[i];
                                        HX_imPtr[i] = imagOut_f[i];
                                    }
                                }
                                else
                                {
                                    gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                    gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                    hAvg_re  = _mm256_loadu_ps(chEstAvg_rePtr);
                                    hAvg_im  = _mm256_loadu_ps(chEstAvg_imPtr);

                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    _mm256_storeu_ps(HX_rePtr, realOut);
                                    _mm256_storeu_ps(HX_imPtr, imagOut);

                                    genDmrs_rePtr += 8;
                                    genDmrs_imPtr += 8;

                                    chEstAvg_rePtr +=8;
                                    chEstAvg_imPtr +=8;

                                    HX_rePtr += 8;
                                    HX_imPtr += 8;
                                }
                            }
                        }
                    }
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 3)
                    {
                        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                        {
                            startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                            endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                            chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][startScIdx];
                            chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][startScIdx];
                            genDmrs_rePtr = &pNrPuschOutParams->genDMRS[idxSym][startScIdx];
                            genDmrs_imPtr = &pNrPuschOutParams->genDMRS[idxSym][(MAX_NUM_SC/2) + startScIdx];
                            HX_rePtr = &HX_re[startScIdx];
                            HX_imPtr = &HX_im[startScIdx];


                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                if(endScIdx - startScIdx<8)
                                {
                                    remSc = (endScIdx) - idxSc;

                                    gDmrs_re = _mm256_setzero_ps();
                                    gDmrs_im = _mm256_setzero_ps();
                                    hAvg_re = _mm256_setzero_ps();
                                    hAvg_im = _mm256_setzero_ps();
                                    HX_AVX_re = _mm256_setzero_ps();
                                    HX_AVX_im = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        gDmrs_re_f[i] = alternateSigns_f[i]*genDmrs_rePtr[i];
                                        gDmrs_im_f[i] = alternateSigns_f[i]*genDmrs_imPtr[i];
                                        hAvg_re_f[i] = chEstAvg_rePtr[i];
                                        hAvg_im_f[i] = chEstAvg_imPtr[i];
                                    }
                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    for(i = 0; i < remSc; i++)
                                    {
                                        HX_rePtr[i] = HX_rePtr[i]+realOut_f[i];
                                        HX_imPtr[i] = HX_imPtr[i]+imagOut_f[i];
                                    }
                                }
                                else
                                {
                                    gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                    gDmrs_re = _mm256_mul_ps(gDmrs_re, alternateSigns);
                                    gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                    gDmrs_im = _mm256_mul_ps(gDmrs_im, alternateSigns);
                                    hAvg_re  = _mm256_loadu_ps(chEstAvg_rePtr);
                                    hAvg_im  = _mm256_loadu_ps(chEstAvg_imPtr);
                                    HX_AVX_re = _mm256_loadu_ps(HX_rePtr);
                                    HX_AVX_im = _mm256_loadu_ps(HX_imPtr);

                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);
                                    _mm256_storeu_ps(HX_rePtr, _mm256_add_ps(HX_AVX_re,realOut));
                                    _mm256_storeu_ps(HX_imPtr, _mm256_add_ps(HX_AVX_im,imagOut));

                                    genDmrs_rePtr += 8;
                                    genDmrs_imPtr += 8;
                                    chEstAvg_rePtr +=8;
                                    chEstAvg_imPtr +=8;
                                    HX_rePtr += 8;
                                    HX_imPtr += 8;
                                }
                            }
                        }
                    }
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 6)
                    {
                        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                        {
                            startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                            endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                            chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][startScIdx];
                            chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][startScIdx];
                            genDmrs_rePtr = &pNrPuschOutParams->genDMRS[idxSym][startScIdx];
                            genDmrs_imPtr = &pNrPuschOutParams->genDMRS[idxSym][(MAX_NUM_SC/2) + startScIdx];
                            HX_rePtr = &HX_re[startScIdx];
                            HX_imPtr = &HX_im[startScIdx];


                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                if(endScIdx - startScIdx<8)
                                {
                                    remSc = (endScIdx) - idxSc;

                                    gDmrs_re = _mm256_setzero_ps();
                                    gDmrs_im = _mm256_setzero_ps();
                                    hAvg_re = _mm256_setzero_ps();
                                    hAvg_im = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        gDmrs_re_f[i] = genDmrs_rePtr[i];
                                        gDmrs_im_f[i] = genDmrs_imPtr[i];
                                        hAvg_re_f[i] = chEstAvg_rePtr[i];
                                        hAvg_im_f[i] = chEstAvg_imPtr[i];
                                    }
                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);
                                    
                                    if(idxSym == 1)
                                    {
                                        //realOut = -1*realOut;
                                        realOut = _mm256_xor_ps(realOut, neg);
                                        //imagOut = -1*imagOut;
                                        imagOut = _mm256_xor_ps(imagOut, neg);
                                    }

                                    for(i = 0; i < remSc; i++)
                                    {
                                        HX_rePtr[i] = HX_rePtr[i]+realOut_f[i];
                                        HX_imPtr[i] = HX_imPtr[i]+imagOut_f[i];
                                    }
                                }
                                else
                                {
                                    gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                    gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                    hAvg_re  = _mm256_loadu_ps(chEstAvg_rePtr);
                                    hAvg_im  = _mm256_loadu_ps(chEstAvg_imPtr);
                                    HX_AVX_re = _mm256_loadu_ps(HX_rePtr);
                                    HX_AVX_im = _mm256_loadu_ps(HX_imPtr);



                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    if(idxSym == 1)
                                    {
                                        //realOut = -1*realOut;
                                        realOut = _mm256_xor_ps(realOut, neg);
                                        //imagOut = -1*imagOut;
                                        imagOut = _mm256_xor_ps(imagOut, neg);
                                    }

                                    _mm256_storeu_ps(HX_rePtr, _mm256_add_ps(HX_AVX_re,realOut));
                                    _mm256_storeu_ps(HX_imPtr, _mm256_add_ps(HX_AVX_im,imagOut));

                                    genDmrs_rePtr += 8;
                                    genDmrs_imPtr += 8;
                                    chEstAvg_rePtr +=8;
                                    chEstAvg_imPtr +=8;
                                    HX_rePtr += 8;
                                    HX_imPtr += 8;
                                }
                            }
                        }
                    }
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 7)
                    {
                        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                        {
                            startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                            endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                            chEstAvg_rePtr = &pNrPuschOutParams->chEstAvg_re[idxAnt][idxLyr][startScIdx];
                            chEstAvg_imPtr = &pNrPuschOutParams->chEstAvg_im[idxAnt][idxLyr][startScIdx];
                            genDmrs_rePtr = &pNrPuschOutParams->genDMRS[idxSym][startScIdx];
                            genDmrs_imPtr = &pNrPuschOutParams->genDMRS[idxSym][(MAX_NUM_SC/2) + startScIdx];
                            HX_rePtr = &HX_re[startScIdx];
                            HX_imPtr = &HX_im[startScIdx];


                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                if(endScIdx - startScIdx<8)
                                {
                                    remSc = (endScIdx) - idxSc;

                                    gDmrs_re = _mm256_setzero_ps();
                                    gDmrs_im = _mm256_setzero_ps();
                                    hAvg_re = _mm256_setzero_ps();
                                    hAvg_im = _mm256_setzero_ps();
                                    HX_AVX_re = _mm256_setzero_ps();
                                    HX_AVX_im = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        gDmrs_re_f[i] = alternateSigns_f[i]*genDmrs_rePtr[i];
                                        gDmrs_im_f[i] = alternateSigns_f[i]*genDmrs_imPtr[i];
                                        hAvg_re_f[i] = chEstAvg_rePtr[i];
                                        hAvg_im_f[i] = chEstAvg_imPtr[i];
                                    }
                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    if(idxSym == 1)
                                    {
                                        //realOut = -1*realOut;
                                        realOut = _mm256_xor_ps(realOut, neg);
                                        //imagOut = -1*imagOut;
                                        imagOut = _mm256_xor_ps(imagOut, neg);
                                    }

                                    for(i = 0; i < remSc; i++)
                                    {
                                        HX_rePtr[i] = HX_rePtr[i]+realOut_f[i];
                                        HX_imPtr[i] = HX_imPtr[i]+imagOut_f[i];
                                    }
                                }
                                else
                                {
                                    gDmrs_re = _mm256_loadu_ps(genDmrs_rePtr);
                                    gDmrs_re = _mm256_mul_ps(gDmrs_re, alternateSigns);
                                    gDmrs_im = _mm256_loadu_ps(genDmrs_imPtr);
                                    gDmrs_im = _mm256_mul_ps(gDmrs_im, alternateSigns);
                                    hAvg_re  = _mm256_loadu_ps(chEstAvg_rePtr);
                                    hAvg_im  = _mm256_loadu_ps(chEstAvg_imPtr);
                                    HX_AVX_re = _mm256_loadu_ps(HX_rePtr);
                                    HX_AVX_im = _mm256_loadu_ps(HX_imPtr);

                                    bd_re = _mm256_mul_ps(hAvg_im, gDmrs_im);
                                    ad_im = _mm256_mul_ps(hAvg_re, gDmrs_im);
                                    realOut = _mm256_fmsub_ps(hAvg_re, gDmrs_re, bd_re);
                                    imagOut = _mm256_fmadd_ps(hAvg_im, gDmrs_re, ad_im);

                                    if(idxSym == 1)
                                    {
                                        //realOut = -1*realOut;
                                        realOut = _mm256_xor_ps(realOut, neg);
                                        //imagOut = -1*imagOut;
                                        imagOut = _mm256_xor_ps(imagOut, neg);
                                    }

                                    _mm256_storeu_ps(HX_rePtr, _mm256_add_ps(HX_AVX_re,realOut));
                                    _mm256_storeu_ps(HX_imPtr, _mm256_add_ps(HX_AVX_im,imagOut));

                                    genDmrs_rePtr += 8;
                                    genDmrs_imPtr += 8;
                                    chEstAvg_rePtr +=8;
                                    chEstAvg_imPtr +=8;
                                    HX_rePtr += 8;
                                    HX_imPtr += 8;
                                }
                            }
                        }
                    }
                }

                for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                {
                    startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                    endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);
                    recDmrs_rePtr = &rxFdSamples[idxAnt][dmrsSym+idxSym][MAX_NUM_SC +startScIdx];
                    recDmrs_imPtr = &rxFdSamples[idxAnt][dmrsSym+idxSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                    HX_rePtr = &HX_re[startScIdx];
                    HX_imPtr = &HX_im[startScIdx];

                    for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                    {
                        if(endScIdx - startScIdx<8)
                        {
                            remSc = endScIdx - idxSc;

                            rDmrs_re = _mm256_setzero_ps();
                            rDmrs_im = _mm256_setzero_ps();

                            HX_AVX_re = _mm256_setzero_ps();
                            HX_AVX_im = _mm256_setzero_ps();

                            for(i = 0; i < remSc; i++)
                            {
                                rDmrs_re_f[i] = recDmrs_rePtr[i];
                                rDmrs_im_f[i] = recDmrs_imPtr[i];
                                HX_AVX_re_f[i] = HX_rePtr[i];
                                HX_AVX_im_f[i] = HX_imPtr[i];
                            }

                            noise_re = _mm256_sub_ps(rDmrs_re, HX_AVX_re);
                            noise_im = _mm256_sub_ps(rDmrs_im, HX_AVX_im);


                            n_avg = _mm256_fmadd_ps(noise_re, noise_re, n_avg);
                            n_avg = _mm256_fmadd_ps(noise_im, noise_im, n_avg);

                        }
                        else
                        {
                            rDmrs_re = _mm256_loadu_ps(recDmrs_rePtr);
                            rDmrs_im = _mm256_loadu_ps(recDmrs_imPtr);
                            HX_AVX_re = _mm256_loadu_ps(HX_rePtr);
                            HX_AVX_im = _mm256_loadu_ps(HX_imPtr);

                            noise_re = _mm256_sub_ps(rDmrs_re, HX_AVX_re);
                            noise_im = _mm256_sub_ps(rDmrs_im, HX_AVX_im);

                            n_avg = _mm256_fmadd_ps(noise_re, noise_re, n_avg);
                            n_avg = _mm256_fmadd_ps(noise_im, noise_im, n_avg);

                            recDmrs_rePtr += 8;
                            recDmrs_imPtr += 8;
                            HX_rePtr += 8;
                            HX_imPtr += 8;
                        }
                    }
                }
            }
        }

        NVTemp = _mm256_permute2f128_ps (n_avg, n_avg, 1);
        nVar_alloc_vec = _mm256_add_ps(n_avg, NVTemp);
        nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);
        nVar_alloc_vec = _mm256_hadd_ps(nVar_alloc_vec, nVar_alloc_vec);
        
        //Find n_avg_1 and average it with n_avg_0.
        pNrPuschOutParams->nVar_overAll = (pNrPuschOutParams->nVar_overAll + ( nVar_alloc_vec_f[0] / (nPrb_total * 6 * nRxAnt * numOfDMRSSymbs) ))*0.5;
    }
}

wnVoid nr_pusch_mmse_perTone_equ_avx2(wnUInt8 nSeg,
                                    wnUInt8 segPart,
                                    wnUInt8* startSym,
                                    wnUInt8* endSym,
                                    commonUlConfigT* pNrUlCommonParams,
                                    puschConfigT* pNrPuschInParams,
                                    P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams)
{
    wnUInt8 idxAlloc, idxSym, idx;
    wnUInt8 idxAnt = 0, idxLyr = 0;
    wnUInt8 rIndx, cIndx;
    wnUInt16 idxTone;

    wnUInt8 remTones;
    wnInt32 toneStartIdx, toneEndIdx;

    wnUInt8 dmrsSym;
    // wnFlt zeroBuf[6] = {0, 0, 0, 0, 0, 0};
    wnUInt32 ScStartIdx_re[8], ScStartIdx_im[8];

    wnUInt32 ntones_alloc, ntones_total, completedTones;
    wnUInt8 nRxAnt, nLyr, nSbAlloc;
    wnUInt8 numOfDMRSSymbs = pNrPuschInParams->nNrOfDMRSSymbols;

    // char prbTemp;

    g_negate = VECT_SET1(-0.0);

    nRxAnt   = pNrPuschInParams->nAntennaPorts;
    nLyr     = pNrPuschInParams->nNrOfLayers;
    nSbAlloc = pNrPuschInParams->nSbAlloc;

    VECT_T H_re[nLyr][nRxAnt];
	VECT_T H_im[nLyr][nRxAnt];
	VECT_T HHh_re[nLyr][nLyr];
	VECT_T HHh_im[nLyr][nLyr];
	VECT_T invHHh_re[nLyr][nLyr];
	VECT_T invHHh_im[nLyr][nLyr];
	VECT_T weights_re[nRxAnt][nLyr];
	VECT_T weights_im[nRxAnt][nLyr];
	VECT_T Y_vect_re[1][nRxAnt];
	VECT_T Y_vect_im[1][nRxAnt];
	VECT_T Y_tilda_re[1][nLyr];
	VECT_T Y_tilda_im[1][nLyr];
	// VECT_T Y_cap_re;
	// VECT_T Y_cap_im;
	// VECT_T conj_sq_toc_re[6];
	// VECT_T conj_sq_toc_im[6];

    // wnFlt* invHHh_re_f = (wnFlt*)&invHHh_re[0][0];
	
    
    //Get DM-RS Symbol Index.
    if(segPart == SEG1)
    {
        dmrsSym = pNrPuschInParams->dmrs1_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG2)
    {
        dmrsSym = pNrPuschInParams->dmrs2_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG3)
    {
        dmrsSym = pNrPuschInParams->dmrs3_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else //if(segPart == SEG4)
    {
        dmrsSym = pNrPuschInParams->dmrs4_sym - pNrPuschInParams->nStartSymbolIndex;
    }


    int total_Prb = 0;
    wnFlt one_By_total_Prb = 0;  
    //Find total number of PRBs allocated.
    for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
    {
        total_Prb += pNrPuschInParams->nPrb[idxAlloc];
    }
    pNrPuschOutParams->total_nPRB = total_Prb;
    one_By_total_Prb = 1.0/total_Prb;

    //Calculate noise variance.
    if(pNrPuschInParams->ZF == 0)
    {
        nr_pusch_mmse_perTone_nVarCalc_avx2(segPart, pNrUlCommonParams, pNrPuschInParams, pNrPuschOutParams);
    }
    else
    {
        pNrPuschOutParams->nVar_overAll = 0;
    }
    VECT_T nVar_vec = VECT_SET1(pNrPuschOutParams->nVar_overAll);
    
    // mexPrintf("nVar_overAll: %f\n", pNrPuschOutParams->nVar_overAll);
    printf("nVar_overAll: %f\n", pNrPuschOutParams->nVar_overAll);

    for(rIndx = 0; rIndx < nLyr; rIndx++)
    {
        pNrPuschOutParams->oneBymse_avg_re[rIndx] = 0;
    }

    ntones_total = 0;
    wnFlt sum;
    wnUInt8 counterOf2;
    wnFlt firstToneSNR;
    wnUInt16 mseCnt = 0;
    wnUInt8 equaSym;
    wnUInt32 layerMapIdx, temp1;

    for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
    {
        completedTones = ntones_total;
        ntones_alloc = pNrPuschInParams->nPrb[idxAlloc]*12;
        ntones_total += ntones_alloc;

        toneStartIdx = pNrPuschInParams->rbStart[idxAlloc]*12;
        toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);

        // wnFlt *chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
        // wnFlt *chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];


        counterOf2 = 1;
        for(idxTone = toneStartIdx; idxTone < toneEndIdx; idxTone += 8)
        {
            for(idx = 0; idx < 8; idx++)
            {
                ScStartIdx_re[idx] = (idxTone + idx) * 1;
                ScStartIdx_im[idx] = MAX_NUM_SC + (idxTone + idx) * 1;
            }

            for(rIndx = 0; rIndx < nLyr; rIndx++)
            {
                for(cIndx = 0; cIndx < nRxAnt; cIndx++)
                {
                    H_re[rIndx][cIndx] = \
                    _mm256_loadu_ps(&pNrPuschOutParams->chEst_re[cIndx][rIndx][idxTone]);
                    
                    H_im[rIndx][cIndx] = \
                    _mm256_loadu_ps(&pNrPuschOutParams->chEst_im[cIndx][rIndx][idxTone]);
                }
            }


            if((idxTone + 8) >= toneEndIdx)
            {
                remTones = toneEndIdx - idxTone;

                for(rIndx = 0; rIndx < nLyr; rIndx++)
                {
                    for(cIndx = 0; cIndx < nRxAnt; cIndx++)
                    {
                        for(idx = remTones; idx < 8; idx++)
                        {
                            H_re[rIndx][cIndx][idx] = 1;
                            H_im[rIndx][cIndx][idx] = 1;
                        }
                    }
                }
            }
            else
            {
                remTones = 8;
            }

            // HH'
            MATRIXMULTTH(H_re, H_im, HHh_re, HHh_im, nLyr, nRxAnt, nLyr);

            for(idx = 0; idx < nLyr; idx++)
            {
                HHh_re[idx][idx] = VECT_ADD(nVar_vec, HHh_re[idx][idx]);
            }


            // Inverse of (HH' + (epsilon)I)
            if(nLyr == 1)
            {
                InverseMatrix1x1HTranspose(HHh_re, HHh_im, invHHh_re, invHHh_im);
            }
            else if(nLyr == 2)
            {
                InverseMatrix2x2H(HHh_re, HHh_im, invHHh_re, invHHh_im);
            }
            else if(nLyr == 4)
            {
                ;
                //InverseMatrix4x4H(HHh_re, HHh_im, invHHh_re, invHHh_im);
            }

            if(pNrPuschInParams->transformPrecode == 0)//If transform precoding is disabled.
            {
                #ifdef SNR24 
                // Calculate one SNR per two PRBs.
                if(counterOf2 == 1)
                {
                    for(rIndx = 0; rIndx < nLyr; rIndx++)
                    {
                        firstToneSNR = 1/invHHh_re[rIndx][rIndx][0];
                        pNrPuschOutParams->oneBymse_re[rIndx][mseCnt] = firstToneSNR;

                        sum = pNrPuschOutParams->oneBymse_avg_re[rIndx];
                        sum = sum + firstToneSNR;
                        pNrPuschOutParams->oneBymse_avg_re[rIndx] = sum;
                    }
                    mseCnt += 1;
                    counterOf2 += 1;
                }
                else
                {
                    if(counterOf2 == 3)
                    {
                        counterOf2 = 1;
                    }
                    else
                    {
                        counterOf2 += 1;
                    }
                }
                #else
                    // Calculate SNR per Tone.
                    for(idx = 0; idx < remTones; idx++)
                    {
                        for(rIndx = 0; rIndx < nLyr; rIndx++)
                        {
                            firstToneSNR = 1/invHHh_re[rIndx][rIndx][idx];
                            pNrPuschOutParams->oneBymse_re[rIndx][mseCnt] = firstToneSNR;
                            sum = pNrPuschOutParams->oneBymse_avg_re[rIndx];
                            sum = sum + invHHh_re[rIndx][rIndx][idx];//firstToneSNR;
                            pNrPuschOutParams->oneBymse_avg_re[rIndx] = sum;
                        }
                        mseCnt += 1;
                    }
                #endif
            }

            // Multiplication of H' and inverse of (HH' + (epsilon)I)
            MATRIXMULTHTTRANSPOSE2(H_re, H_im, invHHh_re, invHHh_im, weights_re, weights_im, nRxAnt, nLyr, nLyr);
            
            equaSym = 0;
            for(idxSym = startSym[segPart - 1]; idxSym <= endSym[segPart - 1]; idxSym++)
            {
                if(idxSym != dmrsSym)
                {
                    temp1 = (equaSym)*pNrPuschOutParams->total_nPRB*12*nLyr + completedTones*nLyr;//(6*(startSym>0)+equaSym)*pNrPuschOutParams->total_nPRB*12*nLyr + 
                        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                        {
                            if( (idxTone + 8) >= toneEndIdx )
                            {
                                Y_vect_re[0][idxAnt] = _mm256_set1_ps(0);
                                Y_vect_im[0][idxAnt] = _mm256_set1_ps(0);
                                remTones = toneEndIdx - idxTone;

                                for(idx = 0; idx < remTones; idx++)
                                {
                                    Y_vect_re[0][idxAnt][idx] = rxFdSamples[idxAnt][idxSym][ScStartIdx_re[idx]];
                                    Y_vect_im[0][idxAnt][idx] = rxFdSamples[idxAnt][idxSym][ScStartIdx_im[idx]];

                                }
                            }
                            else
                            {
                                Y_vect_re[0][idxAnt] = VECT_SET(rxFdSamples[idxAnt][idxSym][ScStartIdx_re[7]],
                                                            rxFdSamples[idxAnt][idxSym][ScStartIdx_re[6]],
                                                            rxFdSamples[idxAnt][idxSym][ScStartIdx_re[5]],
                                                            rxFdSamples[idxAnt][idxSym][ScStartIdx_re[4]],
                                                            rxFdSamples[idxAnt][idxSym][ScStartIdx_re[3]],
                                                            rxFdSamples[idxAnt][idxSym][ScStartIdx_re[2]],
                                                            rxFdSamples[idxAnt][idxSym][ScStartIdx_re[1]],
                                                            rxFdSamples[idxAnt][idxSym][ScStartIdx_re[0]] );


                                Y_vect_im[0][idxAnt] = VECT_SET(rxFdSamples[idxAnt][idxSym][ScStartIdx_im[7]],
                                                                rxFdSamples[idxAnt][idxSym][ScStartIdx_im[6]],
                                                                rxFdSamples[idxAnt][idxSym][ScStartIdx_im[5]],
                                                                rxFdSamples[idxAnt][idxSym][ScStartIdx_im[4]],
                                                                rxFdSamples[idxAnt][idxSym][ScStartIdx_im[3]],
                                                                rxFdSamples[idxAnt][idxSym][ScStartIdx_im[2]],
                                                                rxFdSamples[idxAnt][idxSym][ScStartIdx_im[1]],
                                                                rxFdSamples[idxAnt][idxSym][ScStartIdx_im[0]] );
                                
                            }
                        }

                        MATRIXMULNN(Y_vect_re, Y_vect_im, weights_re, weights_im, Y_tilda_re, Y_tilda_im, 1, nRxAnt, nLyr);


                        for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
                        {
                            for(idx = 0; idx < remTones; idx++)
                            {
                                if(pNrPuschInParams->FOEnable == 1)
                                {
                                    pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][ScStartIdx_re[idx]] = Y_tilda_re[0][idxLyr][idx];
                                    pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][ScStartIdx_im[idx]] = Y_tilda_im[0][idxLyr][idx];
                                }
                                else
                                {
                                    // Directly stores data in their respective Layer Demapped Positions
                                    layerMapIdx = temp1 + (ScStartIdx_re[idx]-toneStartIdx)*nLyr + idxLyr;
                                    pNrPuschOutParams->layerDemapperOutReal[layerMapIdx] = Y_tilda_re[0][idxLyr][idx];
                                    pNrPuschOutParams->layerDemapperOutImag[layerMapIdx] = Y_tilda_im[0][idxLyr][idx];
                                }
                            }
                        }
                    equaSym = equaSym + 1;
                }
                else if(numOfDMRSSymbs == 2)
                {
                    idxSym = idxSym+1;
                }
            }
        }
    }
    totalLength = 12*pNrPuschOutParams->total_nPRB*12*nLyr;

    // Calculate mean of all SNR values.
    #ifdef SNR24
        for(rIndx = 0; rIndx < nLyr; rIndx++)
        {
            pNrPuschOutParams->oneBymse_avg_re[rIndx] = \
                pNrPuschOutParams->oneBymse_avg_re[rIndx]*one_By_total_Prb*0.5;
        }
    #else
        wnFlt oneByNumOfTones = one_By_total_Prb/12;
        for(rIndx = 0; rIndx < nLyr; rIndx++)
        {
            pNrPuschOutParams->oneBymse_avg_re[rIndx] = \
            pNrPuschOutParams->oneBymse_avg_re[rIndx]*oneByNumOfTones;
        }
    #endif



    #ifdef DEBUG
        FILE *fpMSE;
        // FILE *fpMSEAvg;
        // FILE *fpLD, *pfLD;
        wnUInt8 nSymbs;

        if( segPart == SEG1)
        {
            fpMSE = fopen("IOs/normalised_SNR_perTone_Seg1.txt", "w");
        }
        else if(segPart == SEG2)
        {
            fpMSE = fopen("IOs/normalised_SNR_perTone_Seg2.txt", "w");
        }
        else if(segPart == SEG3)
        {
            fpMSE = fopen("IOs/normalised_SNR_perTone_Seg3.txt", "w");
        }
        else if(segPart == SEG4)
        {
            fpMSE = fopen("IOs/normalised_SNR_perTone_Seg4.txt", "w");
        }

        
        #ifdef SNR24
            for(int lyrIdx = 0;lyrIdx<pNrPuschInParams->nNrOfLayers;lyrIdx++)
            {
                for(int toneIdx = 0;toneIdx<pNrPuschInParams->nRBSize;toneIdx++)
                {
                    fprintf(fpMSE, "%f\n ", pNrPuschOutParams->oneBymse_re[lyrIdx][toneIdx]);
                }
            }
        #else
            for(int lyrIdx = 0;lyrIdx<pNrPuschInParams->nNrOfLayers;lyrIdx++)
            {
                for(int toneIdx = 0;toneIdx<pNrPuschInParams->nRBSize*12;toneIdx++)
                {
                    fprintf(fpMSE, "%f\n ", pNrPuschOutParams->oneBymse_re[lyrIdx][toneIdx]/pNrPuschOutParams->oneBymse_avg_re[lyrIdx]);
                }
            }
        #endif

        fclose(fpMSE);

            /*nSymbs = endSym[segPart-1] - startSym[segPart-1];
            if( segPart == SEG1)
            {
                // fpLD = fopen("IOs/layerDemapperReal_mmse_perToneSeg1.txt", "w");
                // pfLD = fopen("IOs/layerDemapperImag_mmse_perToneSeg1.txt", "w");
            }
            else if(segPart == SEG2)
            {
                // fpLD = fopen("IOs/layerDemapperReal_mmse_perToneSeg2.txt", "w");
                // pfLD = fopen("IOs/layerDemapperImag_mmse_perToneSeg2.txt", "w");
            }
            else if(segPart == SEG3)
            {
                fpLD = fopen("IOs/layerDemapperReal_mmse_perToneSeg3.txt", "w");
                pfLD = fopen("IOs/layerDemapperImag_mmse_perToneSeg3.txt", "w");
            }
            else if(segPart == SEG4)
            {
                fpLD = fopen("IOs/layerDemapperReal_mmse_perToneSeg4.txt", "w");
                pfLD = fopen("IOs/layerDemapperImag_mmse_perToneSeg4.txt", "w");
            }

            for(int I1 = 0;I1< pNrPuschOutParams->total_nPRB*12*nLyr*nSymbs;I1++)
            {
                // fprintf(fpLD, "%f\n", pNrPuschOutParams->layerDemapperOutReal[I1]);
                // fprintf(pfLD, "%f\n", pNrPuschOutParams->layerDemapperOutImag[I1]);
            }
            if( segPart == SEG1)
            {
                // fclose(fpLD);
                // fclose(pfLD);
            }
            else
            {
                fclose(fpLD);
                fclose(pfLD);
            }//*/

        if(pNrPuschInParams->FOEnable == 1)
        {
            // printf("equaOutSamples[0][0][0]: %f\n", pNrPuschOutParams->equaOutSamples[0][0][0]);
            FILE *eqRePtr, *eqImPtr;
            nSymbs = endSym[segPart-1] - startSym[segPart-1];
            if(segPart == SEG1)
            {
                eqRePtr = fopen("IOs/equalisedRealSeg1.txt", "w");
                eqImPtr = fopen("IOs/equalisedImagSeg1.txt", "w");
            }
            else if(segPart == SEG2)
            {
                eqRePtr = fopen("IOs/equalisedRealSeg2.txt", "w");
                eqImPtr = fopen("IOs/equalisedImagSeg2.txt", "w");
            }
            else if(segPart == SEG3)
            {
                eqRePtr = fopen("IOs/equalisedRealSeg3.txt", "w");
                eqImPtr = fopen("IOs/equalisedImagSeg3.txt", "w");
            }
            else if(segPart == SEG4)
            {
                eqRePtr = fopen("IOs/equalisedRealSeg4.txt", "w");
                eqImPtr = fopen("IOs/equalisedImagSeg4.txt", "w");
            }

            for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
            {
                toneStartIdx = pNrPuschInParams->rbStart[idxAlloc]*12;
                toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);
                for(idxSym = 0;idxSym<nSymbs;idxSym++)
                {
                    for(idxTone = toneStartIdx; idxTone < toneEndIdx; idxTone++)
                    {
                        for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
                        {
                            fprintf(eqRePtr, "%f\n", pNrPuschOutParams->equaOutSamples[idxLyr][idxSym][idxTone]);
                            fprintf(eqImPtr, "%f\n", pNrPuschOutParams->equaOutSamples[idxLyr][idxSym][MAX_NUM_SC+idxTone]);
                        }
                    }
                }
            }
            fclose(eqRePtr);
            fclose(eqImPtr);
        }
    #endif
}


wnVoid nr_pusch_mmse_perTone_FO_Corr(wnUInt8   segPart,
                                wnUInt8* startSym,
                                wnUInt8* endSym,
                                commonUlConfigT* pNrUlCommonParams,
                                puschConfigT* pNrPuschInParams,
                                P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams)
{
    wnUInt8 dmrsSym;
    wnUInt8 symIdx;
    // wnUInt8 equDataIdx[6];
    wnUInt8 idxCnt;
    wnUInt8 idxAlloc, idxSym, idx;
    wnInt32 toneStartIdx, toneEndIdx;
    // wnUInt16 idxTone;    
    // wnUInt32 ScStartIdx_re[8], ScStartIdx_im[8];  
    wnUInt8 idxLyr = 0;//idxAnt = 0,
    // wnUInt8 remTones;
    wnUInt32 ntones_alloc, ntones_total, completedTones;
    wnUInt8 nRxAnt, nLyr, nSbAlloc;
    wnUInt8 equaSym;
    wnUInt32 layerMapIdx, temp1;

    wnUInt16 numofPRGs;
    wnUInt16 PRGIdx, idxSCUpperLimit;
    wnUInt16 cmpldPRGs;
    wnUInt16 idxSc;

    wnFlt cosVal[6];
    wnFlt sinVal[6];
    __m256 cosVec[6];
    __m256 sinVec[6];
    __m256 eqRe, eqIm;
    // wnFlt* eqRe_f = (wnFlt*)&eqRe;
    // wnFlt* eqIm_f = (wnFlt*)&eqIm;
    __m256 eqReCor, eqImCor;

    wnFlt eqReVal, eqImVal;
    wnUInt32 idxSCLowerLimit;
    wnUInt16 AVXitrCnt;

    nRxAnt   = pNrPuschInParams->nAntennaPorts;
    nLyr     = pNrPuschInParams->nNrOfLayers;
    nSbAlloc = pNrPuschInParams->nSbAlloc;

    if(segPart == SEG1)
    {
        dmrsSym = pNrPuschInParams->dmrs1_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG2)
    {
        dmrsSym = pNrPuschInParams->dmrs2_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG3)
    {
        dmrsSym = pNrPuschInParams->dmrs3_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else //if(segPart == SEG4)
    {
        dmrsSym = pNrPuschInParams->dmrs4_sym - pNrPuschInParams->nStartSymbolIndex;
    }

    if(pNrPuschInParams->FOFullAlloc == 1)
    {
        // The following if-else loops will calculate "cosine" and "sine" values required for FO Correction.
        // if(segPart == SEG1)
        // {
            idxCnt = 0;
            
            for(symIdx = startSym[segPart];symIdx<=endSym[segPart];symIdx++)
            {
                if(symIdx != dmrsSym)
                {                    
                    cosVal[idxCnt] = wnCosRad((double)-1*(symIdx-dmrsSym)*FOEst);
                    sinVal[idxCnt] = wnSinRad((double)-1*(symIdx-dmrsSym)*FOEst);
                    cosVec[idxCnt] = _mm256_set1_ps(cosVal[idxCnt]);
                    sinVec[idxCnt] = _mm256_set1_ps(sinVal[idxCnt]);
                    idxCnt += 1;
                }
            }
        // }
        // else
        // {
        //     idxCnt = 0;
        //     for(symIdx = startSym[segPart];symIdx<=endSym[segPart];symIdx++)
        //     {
        //         if(symIdx != dmrsSym)
        //         {
        //             cosVal[idxCnt] = wnCosRad(-1*(symIdx-dmrsSym)*FOEst);
        //             sinVal[idxCnt] = wnSinRad(-1*(symIdx-dmrsSym)*FOEst);
        //             cosVec[idxCnt] = _mm256_set1_ps(cosVal[idxCnt]);
        //             sinVec[idxCnt] = _mm256_set1_ps(sinVal[idxCnt]);
        //             idxCnt += 1;
        //         }
        //     }
        // }

        ntones_total = 0;
        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
        {
            completedTones = ntones_total;
            ntones_alloc = pNrPuschInParams->nPrb[idxAlloc]*12;
            ntones_total += ntones_alloc;

            toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
            toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);

            idxSCLowerLimit = toneStartIdx;
            idxSCUpperLimit = toneEndIdx;

            for(idxLyr = 0;idxLyr<nLyr;idxLyr++)
            {
                equaSym = 0;
                for(idxSym = startSym[segPart]; idxSym <= endSym[segPart]; idxSym++)
                {
                    if(idxSym != dmrsSym)
                    {
                        temp1 = (equaSym)*pNrPuschOutParams->total_nPRB*12*nLyr + completedTones*nLyr;                          

                        AVXitrCnt = (idxSCUpperLimit - idxSCLowerLimit)>>3;
                        
                        // if (idxSCUpperLimit - idxSCLowerLimit) = 8*X+Y, following for-loop will complete 8*X part.
                        for(idxSc = 0;idxSc<AVXitrCnt;idxSc++)
                        {
                            // FO Correction using AVX.
                            eqRe = _mm256_loadu_ps((wnFlt const*) &(pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][idxSCLowerLimit+idxSc*8]));
                            eqIm = _mm256_loadu_ps((wnFlt const*) &(pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][MAX_NUM_SC+idxSCLowerLimit+idxSc*8]));
                            eqReCor = _mm256_sub_ps(_mm256_mul_ps(eqRe, cosVec[equaSym]), _mm256_mul_ps(eqIm, sinVec[equaSym]));
                            eqImCor = _mm256_add_ps(_mm256_mul_ps(eqRe, sinVec[equaSym]), _mm256_mul_ps(eqIm, cosVec[equaSym]));
                            _mm256_storeu_ps((wnFlt*) &(pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][idxSCLowerLimit+idxSc*8]), eqReCor);
                            _mm256_storeu_ps((wnFlt*) &(pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][MAX_NUM_SC+idxSCLowerLimit+idxSc*8]), eqImCor);

                            for(idx = 0; idx < 8; idx++)
                            {
                                // Directly stores data in their respective Layer Demapped Positions
                                layerMapIdx = temp1 + (idxSCLowerLimit+idxSc*8+idx-toneStartIdx)*nLyr + idxLyr;
                                pNrPuschOutParams->layerDemapperOutReal[layerMapIdx] = pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][idxSCLowerLimit+idxSc*8+idx];
                                pNrPuschOutParams->layerDemapperOutImag[layerMapIdx] = pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][MAX_NUM_SC+idxSCLowerLimit+idxSc*8+idx];
                            }
                        }

                        // if (idxSCUpperLimit - idxSCLowerLimit) = 8*X+Y, following for-loop will complete Y part.
                        for(idxSc = idxSCLowerLimit+(AVXitrCnt<<3);idxSc<idxSCUpperLimit;idxSc++)
                        {
                            // FO Correction using normal C code.
                            eqReVal = pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][idxSc]*cosVal[equaSym] -\
                                                                                        pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][MAX_NUM_SC+idxSc]*sinVal[equaSym];
                            eqImVal = pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][idxSc]*sinVal[equaSym]+\
                                                                                                pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][MAX_NUM_SC+idxSc]*cosVal[equaSym];

                            pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][idxSc] = eqReVal;    
                            pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][MAX_NUM_SC+idxSc] = eqImVal;

                            // Directly stores data in their respective Layer Demapped Positions
                            layerMapIdx = temp1 + (idxSc-toneStartIdx)*nLyr + idxLyr;
                            pNrPuschOutParams->layerDemapperOutReal[layerMapIdx] = pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][idxSc];
                            pNrPuschOutParams->layerDemapperOutImag[layerMapIdx] = pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][MAX_NUM_SC+idxSc];
                        }
                        equaSym += 1;
                    }
                }
	        }
        }
    }
    else
    {
        cmpldPRGs = 0;
        ntones_total = 0;
        for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
        {
            completedTones = ntones_total;
            ntones_alloc = pNrPuschInParams->nPrb[idxAlloc]*12;
            ntones_total += ntones_alloc;

            toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
            toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);//

            // numofPRGs = ceil(nPrb[idxAlloc]/FOPRG);
            numofPRGs = (wnUInt16)(pNrPuschInParams->nPrb[idxAlloc]/pNrPuschInParams->FOPRG);
            if((pNrPuschInParams->nPrb[idxAlloc] - (wnUInt16)(numofPRGs*pNrPuschInParams->FOPRG))>0)
            {
                numofPRGs = numofPRGs+1;
            }

            for(idxLyr = 0;idxLyr<nLyr;idxLyr++)
            {
                equaSym = 0;
                for(idxSym = startSym[segPart]; idxSym <= endSym[segPart]; idxSym++)
                {
                    if(idxSym != dmrsSym)
                    {
                        temp1 = (equaSym)*pNrPuschOutParams->total_nPRB*12*nLyr + completedTones*nLyr;

                        for(PRGIdx = 0;PRGIdx<numofPRGs;PRGIdx++)
                        {   
                            cosVal[equaSym] = wnCosRad(-1*(idxSym-dmrsSym)*FOEst_PRG[cmpldPRGs+PRGIdx]);
                            sinVal[equaSym] = wnSinRad(-1*(idxSym-dmrsSym)*FOEst_PRG[cmpldPRGs+PRGIdx]);

                            cosVec[equaSym] = _mm256_set1_ps(cosVal[equaSym]);
                            sinVec[equaSym] = _mm256_set1_ps(sinVal[equaSym]);                            
                            
                            idxSCLowerLimit = toneStartIdx+(wnUInt16)(PRGIdx*pNrPuschInParams->FOPRG*12);
                            idxSCUpperLimit = toneStartIdx+(wnUInt16)((PRGIdx+1)*pNrPuschInParams->FOPRG*12);

                            // If number of PRBS allocated are multiples of FOPRG, then take care of upper limit.
                            if(PRGIdx == numofPRGs-1)
                            {
                                if(toneEndIdx<idxSCUpperLimit)
                                {
                                    idxSCUpperLimit = toneEndIdx;
                                }
                            }

                            AVXitrCnt = (idxSCUpperLimit - idxSCLowerLimit)>>3;
                            
                            // if (idxSCUpperLimit - idxSCLowerLimit) = 8*X+Y, following for-loop will complete 8*X part.
                            for(idxSc = 0;idxSc<AVXitrCnt;idxSc++)
                            {
                                // FO Correction using AVX.
                                eqRe = _mm256_loadu_ps((float const*) &(pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][idxSCLowerLimit+idxSc*8]));
                                eqIm = _mm256_loadu_ps((float const*) &(pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][MAX_NUM_SC+idxSCLowerLimit+idxSc*8]));
                                eqReCor = _mm256_sub_ps(_mm256_mul_ps(eqRe, cosVec[equaSym]), _mm256_mul_ps(eqIm, sinVec[equaSym]));
                                eqImCor = _mm256_add_ps(_mm256_mul_ps(eqRe, sinVec[equaSym]), _mm256_mul_ps(eqIm, cosVec[equaSym]));
                                _mm256_storeu_ps((float*) &(pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][idxSCLowerLimit+idxSc*8]), eqReCor);
                                _mm256_storeu_ps((float*) &(pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][MAX_NUM_SC+idxSCLowerLimit+idxSc*8]), eqImCor);


                                for(idx = 0; idx < 8; idx++)
                                {
                                    // Directly stores data in their respective Layer Demapped Positions
                                    layerMapIdx = temp1 + (idxSCLowerLimit+idxSc*8+idx-toneStartIdx)*nLyr + idxLyr;
                                    pNrPuschOutParams->layerDemapperOutReal[layerMapIdx] = pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][idxSCLowerLimit+idxSc*8+idx];
                                    pNrPuschOutParams->layerDemapperOutImag[layerMapIdx] = pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][MAX_NUM_SC+idxSCLowerLimit+idxSc*8+idx];
                                }
                            }

                            // if (idxSCUpperLimit - idxSCLowerLimit) = 8*X+Y, following for-loop will complete Y part.
                            for(idxSc = idxSCLowerLimit+(AVXitrCnt<<3);idxSc<idxSCUpperLimit;idxSc++)
                            {
                                // FO Correction using normal C code.
                                eqReVal = pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][idxSc]*cosVal[equaSym] -\
                                                                                            pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][MAX_NUM_SC+idxSc]*sinVal[equaSym];
                                eqImVal = pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][idxSc]*sinVal[equaSym]+\
                                                                                                    pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][MAX_NUM_SC+idxSc]*cosVal[equaSym];

                                pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][idxSc] = eqReVal;    
                                pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][MAX_NUM_SC+idxSc] = eqImVal;

                                // Directly stores data in their respective Layer Demapped Positions
                                layerMapIdx = temp1 + (idxSc-toneStartIdx)*nLyr + idxLyr;
                                pNrPuschOutParams->layerDemapperOutReal[layerMapIdx] = pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][idxSc];
                                pNrPuschOutParams->layerDemapperOutImag[layerMapIdx] = pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][MAX_NUM_SC+idxSc];
                            }
                            
                        }
                        equaSym += 1;
                    }
                }
	        }
            cmpldPRGs = cmpldPRGs+numofPRGs;
        }
    }

    #ifdef DEBUG
        FILE *fpReal, *fpImag;
        wnUInt32 dataIdx;
        if(segPart == SEG1)
        {
            fpReal = fopen("IOs/LayerDemapperRealSeg1PerTone.txt", "w");
            fpImag = fopen("IOs/LayerDemapperImagSeg1PerTone.txt", "w");
        }
        else if(segPart == SEG2)
        {
            fpReal = fopen("IOs/LayerDemapperRealSeg2PerTone.txt", "w");
            fpImag = fopen("IOs/LayerDemapperImagSeg2PerTone.txt", "w");
        }
        else if(segPart == SEG3)
        {
            fpReal = fopen("IOs/LayerDemapperRealSeg3PerTone.txt", "w");
            fpImag = fopen("IOs/LayerDemapperImagSeg3PerTone.txt", "w");
        }
        else if(segPart == SEG4)
        {
            fpReal = fopen("IOs/LayerDemapperRealSeg4PerTone.txt", "w");
            fpImag = fopen("IOs/LayerDemapperImagSeg4PerTone.txt", "w");
        }

        for(dataIdx = 0;dataIdx<pNrPuschInParams->nRBSize*12*pNrPuschInParams->nNrOfLayers*6;dataIdx++)
        {
            fprintf(fpReal, "%f\n", pNrPuschOutParams->layerDemapperOutReal[dataIdx]);
            fprintf(fpImag, "%f\n", pNrPuschOutParams->layerDemapperOutImag[dataIdx]);
        }
        fclose(fpReal);
        fclose(fpImag);
    #endif


}
// }//Per Tone End

void nr_pusch_mmse_perTone_DS_est_avx2(wnUInt8 segPart,
                               commonUlConfigT* pNrUlCommonParams,
                               puschConfigT* pNrPuschInParams,
                               P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams)
{
    wnUInt16 idxAnt, idxLyr, idxAlloc, idxSc, i;
    wnUInt16 startScIdx, endScIdx, remSc;
    wnUInt16 toneStartIdx, toneEndIdx;
    wnUInt8 dmrsSym, nRxAnt, nLyr, nSbAlloc;

    __m256 bd_re, ad_im;
    __m256 out_re, out_im;
    wnFlt *out_re_f = (wnFlt*)&out_re;
    wnFlt *out_im_f = (wnFlt*)&out_im;

    __m256 rDmrs_re1, rDmrs_im1;
    __m256 gDmrs_re1, gDmrs_im1;
    __m256 rDmrs_re2, rDmrs_im2;
    __m256 gDmrs_re2, gDmrs_im2;
    __m256 out_re1, out_im1;
    __m256 out_re2, out_im2;
    __m256 temp_re1, temp_im1;
    __m256 temp_re2, temp_im2;

    wnFlt *out_re1_f = (wnFlt*)&out_re1;
    wnFlt *out_im1_f = (wnFlt*)&out_im1;
    wnFlt *out_re2_f = (wnFlt*)&out_re2;
    wnFlt *out_im2_f = (wnFlt*)&out_im2;

    __m256 zeroVec = _mm256_setzero_ps();
    __m256 div2Vec = _mm256_set1_ps(0.5);
    __m256 dmrsScale = _mm256_set1_ps(0.501187);
    // __m256 ls_chEst_re, ls_chEst_im;
    // __m256 ls_chEst_re2, ls_chEst_im2;
    __m256* AVX_chEst_rePtr, *AVX_chEst_imPtr;
    __m256 AVX_chEst_re, AVX_chEst_im;

    wnFlt* rDmrs_re_f1 = (wnFlt*)&rDmrs_re1;
    wnFlt* rDmrs_im_f1 = (wnFlt*)&rDmrs_im1;
    wnFlt* gDmrs_re_f1 = (wnFlt*)&gDmrs_re1;
    wnFlt* gDmrs_im_f1 = (wnFlt*)&gDmrs_im1;
    wnFlt* rDmrs_re_f2 = (wnFlt*)&rDmrs_re2;
    wnFlt* rDmrs_im_f2 = (wnFlt*)&rDmrs_im2;
    wnFlt* gDmrs_re_f2 = (wnFlt*)&gDmrs_re2;
    wnFlt* gDmrs_im_f2 = (wnFlt*)&gDmrs_im2;

    wnFlt *recDmrs_rePtr1;
    wnFlt *recDmrs_imPtr1;
    wnFlt *recDmrs_rePtr2;
    wnFlt *recDmrs_imPtr2;
    wnFlt *genDmrs_rePtr1;
    wnFlt *genDmrs_imPtr1;
    wnFlt *genDmrs_rePtr2;
    wnFlt *genDmrs_imPtr2;
    // wnFlt *lsEst_rePtr;
    // wnFlt *lsEst_imPtr;
    wnFlt *chEst_rePtr, *chEst_imPtr;
    // wnUInt8 CDMFlag0 = 0, CDMFlag1 = 0;
    wnUInt8 cdmTotal[2] = {0, 0};

    __m256 alternateSigns = _mm256_setr_ps(1, -1, 1, -1, 1, -1, 1, -1);
    wnFlt* alternateSigns_f = (wnFlt *)&alternateSigns;

    wnFlt cdm0EstRe[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][2][MAX_NUM_ALLOC*6];
    wnFlt cdm0EstIm[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][2][MAX_NUM_ALLOC*6];
    wnFlt *cdm0Est_rePtr1, *cdm0Est_imPtr1;
    wnFlt *cdm0Est_rePtr2, *cdm0Est_imPtr2;
    __m256 cdm0EstRe1, cdm0EstRe2, cdm0EstIm1, cdm0EstIm2;
    
    wnFlt cdm1EstRe[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][2][MAX_NUM_ALLOC*6];
    wnFlt cdm1EstIm[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][2][MAX_NUM_ALLOC*6];
    wnFlt *cdm1Est_rePtr1, *cdm1Est_imPtr1;
    wnFlt *cdm1Est_rePtr2, *cdm1Est_imPtr2;
    __m256 cdm1EstRe1, cdm1EstRe2, cdm1EstIm1, cdm1EstIm2;

    // cdmTotal[0] => number of ports used from CDM group 0.
    // cdmTotal[1] => number of ports used from CDM group 1.
    for(idxLyr = 0;idxLyr<pNrPuschInParams->nNrOfLayers;idxLyr++)
    {
        if(pNrPuschInParams->nPortIndex[idxLyr] == 0 || pNrPuschInParams->nPortIndex[idxLyr] == 1 ||
           pNrPuschInParams->nPortIndex[idxLyr] == 4 || pNrPuschInParams->nPortIndex[idxLyr] == 5)
        {
            cdmTotal[0] = cdmTotal[0]+1;
        }
        else
        {
            cdmTotal[1] = cdmTotal[1]+1;
        }
    }

    __m256i const1 = _mm256_set_epi32(1, 1, 6, 1, 4, 1, 2, 1);
    __m256i const2 = _mm256_set_epi32(0, 6, 5, 4, 3, 2, 1, 2);

    // wnUInt8 dmrsIdx = 0;

    // Get DM-RS Index
    if(segPart == SEG1)
    {
        dmrsSym = pNrPuschInParams->dmrs1_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG2)
    {
        dmrsSym = pNrPuschInParams->dmrs2_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG3)
    {
        dmrsSym = pNrPuschInParams->dmrs3_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG4)
    {
        dmrsSym = pNrPuschInParams->dmrs4_sym - pNrPuschInParams->nStartSymbolIndex;
    }

    nRxAnt   = pNrPuschInParams->nAntennaPorts;
    nLyr     = pNrPuschInParams->nNrOfLayers;
    nSbAlloc = pNrPuschInParams->nSbAlloc;

    if(cdmTotal[0]>=2)
    {
        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
            {
                startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
                toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);

                //Received DMRS Signal 1                          
                recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(MAX_NUM_SC/2) + startScIdx];
                //Received DMRS Signal 2
                recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][startScIdx];
                recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(MAX_NUM_SC/2) + startScIdx];
                //Generated DMRS Signal 1
                genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                //Generated DMRS Signal 2
                genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];

                cdm0Est_rePtr1 = &cdm0EstRe[idxAnt][0][startScIdx];
                cdm0Est_imPtr1 = &cdm0EstIm[idxAnt][0][startScIdx];

                cdm0Est_rePtr2 = &cdm0EstRe[idxAnt][1][startScIdx];
                cdm0Est_imPtr2 = &cdm0EstIm[idxAnt][1][startScIdx];

                for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                {
                    remSc = endScIdx - idxSc;
                    if( remSc < 8)
                    {
                        rDmrs_re1 = _mm256_setzero_ps();
                        rDmrs_im1 = _mm256_setzero_ps();
                        rDmrs_re2 = _mm256_setzero_ps();
                        rDmrs_im2 = _mm256_setzero_ps();

                        gDmrs_re1 = _mm256_setzero_ps();
                        gDmrs_im1 = _mm256_setzero_ps();
                        gDmrs_re2 = _mm256_setzero_ps();
                        gDmrs_im2 = _mm256_setzero_ps();

                        for(i = 0; i < remSc; i++)
                        {
                            rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                            rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                            rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                            rDmrs_im_f2[i] = recDmrs_imPtr2[i];

                            gDmrs_re_f1[i] = genDmrs_rePtr1[i];
                            gDmrs_im_f1[i] = genDmrs_imPtr1[i];
                            gDmrs_re_f2[i] = genDmrs_rePtr2[i];
                            gDmrs_im_f2[i] = genDmrs_imPtr2[i];
                        }

                        bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                        ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);

                        out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                        out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);

                        // out1 = ( Y1*conj(Xdmrs1) )/Beta_dmrs;
                        out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                        out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                        bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                        ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);

                        out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                        out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);

                        // out2 = ( Y2*conj(Xdmrs2) )/Beta_dmrs;
                        out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                        out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                        for(i = 0; i < remSc; i++)
                        {
                            cdm0Est_rePtr1[i] = out_re1_f[i];
                            cdm0Est_imPtr1[i] = out_im1_f[i];
                            cdm0Est_rePtr2[i] = out_re2_f[i];
                            cdm0Est_imPtr2[i] = out_im2_f[i];
                        }
                    }
                    else
                    {
                        rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                        rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                        gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                        gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);

                        rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                        rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                        gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                        gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);

                        bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                        ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);

                        out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                        out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);

                        // out1 = ( Y1*conj(Xdmrs1) )/Beta_dmrs;
                        out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                        out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                        bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                        ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);

                        out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                        out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);

                        // out2 = ( Y2*conj(Xdmrs2) )/Beta_dmrs;
                        out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                        out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                        _mm256_storeu_ps(cdm0Est_rePtr1, out_re1);
                        _mm256_storeu_ps(cdm0Est_imPtr1, out_im1);

                        _mm256_storeu_ps(cdm0Est_rePtr2, out_re2);
                        _mm256_storeu_ps(cdm0Est_imPtr2, out_im2);
                    }
                    cdm0Est_rePtr1 += 8;
                    cdm0Est_imPtr1 += 8;

                    cdm0Est_rePtr2 += 8;
                    cdm0Est_imPtr2 += 8;

                    recDmrs_rePtr1 += 8;
                    recDmrs_imPtr1 += 8;

                    recDmrs_rePtr2 += 8;
                    recDmrs_imPtr2 += 8;

                    genDmrs_rePtr1 += 8;
                    genDmrs_imPtr1 += 8;

                    genDmrs_rePtr2 += 8;
                    genDmrs_imPtr2 += 8;
                }
            }
        }
        

        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
            {
                for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                {
                    startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                    endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                    toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
                    toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);

                    cdm0Est_rePtr1 = &cdm0EstRe[idxAnt][0][startScIdx];
                    cdm0Est_imPtr1 = &cdm0EstIm[idxAnt][0][startScIdx];

                    cdm0Est_rePtr2 = &cdm0EstRe[idxAnt][1][startScIdx];
                    cdm0Est_imPtr2 = &cdm0EstIm[idxAnt][1][startScIdx];

                    //Channel Estimates after Linear Interpolation
                    chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                    chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];

                    if(pNrPuschInParams->nPortIndex[idxLyr] == 0)
                    {
                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc <= 8)
                            {
                                for(i = 0; i < remSc-1; i++)
                                {
                                    chEst_rePtr[2*i] = ((cdm0Est_rePtr1[i] + cdm0Est_rePtr1[i+1]) + (cdm0Est_rePtr2[i]+cdm0Est_rePtr2[i+1]))*0.25;
                                    chEst_rePtr[2*i+1] = 0;
                                    chEst_imPtr[2*i] = ((cdm0Est_imPtr1[i] + cdm0Est_imPtr1[i+1]) + (cdm0Est_imPtr2[i]+cdm0Est_imPtr2[i+1]))*0.25;
                                    chEst_imPtr[2*i+1] = 0;
                                }
                                chEst_rePtr[2*(remSc-1)] = chEst_rePtr[2*(remSc-2)];
                                chEst_rePtr[2*(remSc-1)+1] = 0;
                                chEst_imPtr[2*(remSc-1)] = chEst_imPtr[2*(remSc-2)];
                                chEst_imPtr[2*(remSc-1)+1] = 0;
                            }
                            else
                            {
                                cdm0EstRe1 = _mm256_loadu_ps(cdm0Est_rePtr1);
                                cdm0EstRe2 = _mm256_loadu_ps(cdm0Est_rePtr1+1);
                                cdm0EstIm1 = _mm256_loadu_ps(cdm0Est_imPtr1);
                                cdm0EstIm2 = _mm256_loadu_ps(cdm0Est_imPtr1+1);

                                out_re1 = _mm256_add_ps(cdm0EstRe1, cdm0EstRe2);
                                out_im1 = _mm256_add_ps(cdm0EstIm1, cdm0EstIm2);
                                out_re1 = _mm256_mul_ps(out_re1, div2Vec);
                                out_im1 = _mm256_mul_ps(out_im1, div2Vec);

                                cdm0EstRe1 = _mm256_loadu_ps(cdm0Est_rePtr2);
                                cdm0EstRe2 = _mm256_loadu_ps(cdm0Est_rePtr2+1);
                                cdm0EstIm1 = _mm256_loadu_ps(cdm0Est_imPtr2);
                                cdm0EstIm2 = _mm256_loadu_ps(cdm0Est_imPtr2+1);

                                out_re2 = _mm256_add_ps(cdm0EstRe1, cdm0EstRe2);
                                out_im2 = _mm256_add_ps(cdm0EstIm1, cdm0EstIm2);
                                out_re2 = _mm256_mul_ps(out_re2, div2Vec);
                                out_im2 = _mm256_mul_ps(out_im2, div2Vec);

                                // Mean of two adjacent estimates.
                                out_re = _mm256_add_ps(out_re1, out_re2);
                                out_im = _mm256_add_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                
                            }
                            cdm0Est_rePtr1 += 8;
                            cdm0Est_imPtr1 += 8;

                            cdm0Est_rePtr2 += 8;
                            cdm0Est_imPtr2 += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }
                    }
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 1)
                    {
                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc <= 8)
                            {
                                for(i = 0; i < remSc-1; i++)
                                {
                                    chEst_rePtr[2*i] = alternateSigns_f[i]*((cdm0Est_rePtr1[i] - cdm0Est_rePtr1[i+1]) + (cdm0Est_rePtr2[i] - cdm0Est_rePtr2[i+1]))*0.25;
                                    chEst_rePtr[2*i+1] = 0;
                                    chEst_imPtr[2*i] = alternateSigns_f[i]*((cdm0Est_imPtr1[i] - cdm0Est_imPtr1[i+1]) + (cdm0Est_imPtr2[i] - cdm0Est_imPtr2[i+1]))*0.25;
                                    chEst_imPtr[2*i+1] = 0;
                                }
                                chEst_rePtr[2*(remSc-1)] = chEst_rePtr[2*(remSc-2)];
                                chEst_rePtr[2*(remSc-1)+1] = 0;
                                chEst_imPtr[2*(remSc-1)] = chEst_imPtr[2*(remSc-2)];
                                chEst_imPtr[2*(remSc-1)+1] = 0;
                            }
                            else
                            {
                                cdm0EstRe1 = _mm256_loadu_ps(cdm0Est_rePtr1);
                                cdm0EstRe2 = _mm256_loadu_ps(cdm0Est_rePtr1+1);
                                cdm0EstIm1 = _mm256_loadu_ps(cdm0Est_imPtr1);
                                cdm0EstIm2 = _mm256_loadu_ps(cdm0Est_imPtr1+1);

                                out_re1 = _mm256_sub_ps(cdm0EstRe1, cdm0EstRe2);
                                out_re1 = _mm256_mul_ps(out_re1, alternateSigns);
                                out_im1 = _mm256_sub_ps(cdm0EstIm1, cdm0EstIm2);
                                out_im1 = _mm256_mul_ps(out_im1, alternateSigns);
                                out_re1 = _mm256_mul_ps(out_re1, div2Vec);
                                out_im1 = _mm256_mul_ps(out_im1, div2Vec);

                                cdm0EstRe1 = _mm256_loadu_ps(cdm0Est_rePtr2);
                                cdm0EstRe2 = _mm256_loadu_ps(cdm0Est_rePtr2+1);
                                cdm0EstIm1 = _mm256_loadu_ps(cdm0Est_imPtr2);
                                cdm0EstIm2 = _mm256_loadu_ps(cdm0Est_imPtr2+1);

                                out_re2 = _mm256_sub_ps(cdm0EstRe1, cdm0EstRe2);
                                out_re2 = _mm256_mul_ps(out_re2, alternateSigns);
                                out_im2 = _mm256_sub_ps(cdm0EstIm1, cdm0EstIm2);
                                out_im2 = _mm256_mul_ps(out_im2, alternateSigns);
                                out_re2 = _mm256_mul_ps(out_re2, div2Vec);
                                out_im2 = _mm256_mul_ps(out_im2, div2Vec);

                                // Mean of two adjacent estimates.
                                out_re = _mm256_add_ps(out_re1, out_re2);
                                out_im = _mm256_add_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                
                            }
                            cdm0Est_rePtr1 += 8;
                            cdm0Est_imPtr1 += 8;

                            cdm0Est_rePtr2 += 8;
                            cdm0Est_imPtr2 += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }
                    }
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 4)
                    {
                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc <= 8)
                            {
                                for(i = 0; i < remSc-1; i++)
                                {
                                    chEst_rePtr[2*i] = ((cdm0Est_rePtr1[i] + cdm0Est_rePtr1[i+1]) - (cdm0Est_rePtr2[i] + cdm0Est_rePtr2[i+1]))*0.25;
                                    chEst_rePtr[2*i+1] = 0;
                                    chEst_imPtr[2*i] = ((cdm0Est_imPtr1[i] + cdm0Est_imPtr1[i+1]) - (cdm0Est_imPtr2[i] + cdm0Est_imPtr2[i+1]))*0.25;
                                    chEst_imPtr[2*i+1] = 0;
                                }
                                chEst_rePtr[2*(remSc-1)] = chEst_rePtr[2*(remSc-2)];
                                chEst_rePtr[2*(remSc-1)+1] = 0;
                                chEst_imPtr[2*(remSc-1)] = chEst_imPtr[2*(remSc-2)];
                                chEst_imPtr[2*(remSc-1)+1] = 0;
                            }
                            else
                            {
                                cdm0EstRe1 = _mm256_loadu_ps(cdm0Est_rePtr1);
                                cdm0EstRe2 = _mm256_loadu_ps(cdm0Est_rePtr1+1);
                                cdm0EstIm1 = _mm256_loadu_ps(cdm0Est_imPtr1);
                                cdm0EstIm2 = _mm256_loadu_ps(cdm0Est_imPtr1+1);

                                out_re1 = _mm256_add_ps(cdm0EstRe1, cdm0EstRe2);
                                out_im1 = _mm256_add_ps(cdm0EstIm1, cdm0EstIm2);
                                out_re1 = _mm256_mul_ps(out_re1, div2Vec);
                                out_im1 = _mm256_mul_ps(out_im1, div2Vec);

                                cdm0EstRe1 = _mm256_loadu_ps(cdm0Est_rePtr2);
                                cdm0EstRe2 = _mm256_loadu_ps(cdm0Est_rePtr2+1);
                                cdm0EstIm1 = _mm256_loadu_ps(cdm0Est_imPtr2);
                                cdm0EstIm2 = _mm256_loadu_ps(cdm0Est_imPtr2+1);

                                out_re2 = _mm256_add_ps(cdm0EstRe1, cdm0EstRe2);
                                out_im2 = _mm256_add_ps(cdm0EstIm1, cdm0EstIm2);
                                out_re2 = _mm256_mul_ps(out_re2, div2Vec);
                                out_im2 = _mm256_mul_ps(out_im2, div2Vec);

                                // Mean of two adjacent estimates.
                                out_re = _mm256_sub_ps(out_re1, out_re2);
                                out_im = _mm256_sub_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                
                            }
                            cdm0Est_rePtr1 += 8;
                            cdm0Est_imPtr1 += 8;

                            cdm0Est_rePtr2 += 8;
                            cdm0Est_imPtr2 += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }
                    }
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 5)
                    {
                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc <= 8)
                            {
                                for(i = 0; i < remSc-1; i++)
                                {
                                    chEst_rePtr[2*i] = alternateSigns_f[i]*((cdm0Est_rePtr1[i] - cdm0Est_rePtr1[i+1]) - (cdm0Est_rePtr2[i] - cdm0Est_rePtr2[i+1]))*0.25;
                                    chEst_rePtr[2*i+1] = 0;
                                    chEst_imPtr[2*i] = alternateSigns_f[i]*((cdm0Est_imPtr1[i] - cdm0Est_imPtr1[i+1]) - (cdm0Est_imPtr2[i] - cdm0Est_imPtr2[i+1]))*0.25;
                                    chEst_imPtr[2*i+1] = 0;
                                }
                                chEst_rePtr[2*(remSc-1)] = chEst_rePtr[2*(remSc-2)];
                                chEst_rePtr[2*(remSc-1)+1] = 0;
                                chEst_imPtr[2*(remSc-1)] = chEst_imPtr[2*(remSc-2)];
                                chEst_imPtr[2*(remSc-1)+1] = 0;
                            }
                            else
                            {
                                cdm0EstRe1 = _mm256_loadu_ps(cdm0Est_rePtr1);
                                cdm0EstRe2 = _mm256_loadu_ps(cdm0Est_rePtr1+1);
                                cdm0EstIm1 = _mm256_loadu_ps(cdm0Est_imPtr1);
                                cdm0EstIm2 = _mm256_loadu_ps(cdm0Est_imPtr1+1);

                                out_re1 = _mm256_sub_ps(cdm0EstRe1, cdm0EstRe2);
                                out_re1 = _mm256_mul_ps(out_re1, alternateSigns);
                                out_im1 = _mm256_sub_ps(cdm0EstIm1, cdm0EstIm2);
                                out_im1 = _mm256_mul_ps(out_im1, alternateSigns);
                                out_re1 = _mm256_mul_ps(out_re1, div2Vec);
                                out_im1 = _mm256_mul_ps(out_im1, div2Vec);

                                cdm0EstRe1 = _mm256_loadu_ps(cdm0Est_rePtr2);
                                cdm0EstRe2 = _mm256_loadu_ps(cdm0Est_rePtr2+1);
                                cdm0EstIm1 = _mm256_loadu_ps(cdm0Est_imPtr2);
                                cdm0EstIm2 = _mm256_loadu_ps(cdm0Est_imPtr2+1);

                                out_re2 = _mm256_sub_ps(cdm0EstRe1, cdm0EstRe2);
                                out_re2 = _mm256_mul_ps(out_re2, alternateSigns);
                                out_im2 = _mm256_sub_ps(cdm0EstIm1, cdm0EstIm2);
                                out_im2 = _mm256_mul_ps(out_im2, alternateSigns);
                                out_re2 = _mm256_mul_ps(out_re2, div2Vec);
                                out_im2 = _mm256_mul_ps(out_im2, div2Vec);

                                // Mean of two adjacent estimates.
                                out_re = _mm256_sub_ps(out_re1, out_re2);
                                out_im = _mm256_sub_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                
                            }
                            cdm0Est_rePtr1 += 8;
                            cdm0Est_imPtr1 += 8;

                            cdm0Est_rePtr2 += 8;
                            cdm0Est_imPtr2 += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }
                    }


                    if(pNrPuschInParams->nPortIndex[idxLyr] == 0 || pNrPuschInParams->nPortIndex[idxLyr] == 1 ||\
                       pNrPuschInParams->nPortIndex[idxLyr] == 4 || pNrPuschInParams->nPortIndex[idxLyr] == 5)
                    {
                        //Frequency Interpolation.
                        //Calculate estimates for null tones using Left and Right estimates.
                        //For dead end null tone, its two adjacent estimates will be used.
                        //Inp : [ a 0 b 0 c ...]
                        //Out : [a (a+b)/2 b (b+c)/2 c...]
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];
                        AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                        AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                        for(idxSc = toneStartIdx; idxSc < toneEndIdx; idxSc += 8)
                        {
                            remSc = toneEndIdx - idxSc;
                            if(remSc<=8)
                            {
                                for(i = 0; i < remSc; i += 2)
                                {
                                    if(remSc - i == 2)
                                    {
                                        chEst_rePtr[i+1] = (chEst_rePtr[i]+chEst_rePtr[i-1])*0.5;
                                        chEst_imPtr[i+1] = (chEst_imPtr[i]+chEst_imPtr[i-1])*0.5;
                                    }
                                    else
                                    {
                                        chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                        chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                    }
                                }
                            }
                            else
                            {
                                AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                chEst_rePtr += 8;
                                chEst_imPtr += 8;
                            }
                        }


                    }
                    
                }
            }
        }

        if(cdmTotal[1] == 1)
        {
            for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
            {
                for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
                {
                    for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                    {
                        startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                        endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                        toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
                        toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);

                        // port 2 processing
                        if(pNrPuschInParams->nPortIndex[idxLyr] == 2)
                        {
                            //Received DMRS Signal 1.
                            recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC + startScIdx];
                            recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                            //DMRS Signal 1.
                            genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                            genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                            //Received DMRS Signal 2.
                            recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][MAX_NUM_SC + startScIdx];
                            recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                            //DMRS Signal 2.
                            genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                            genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];
                            //Channel Estimates after Linear Interpolation.
                            chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                            chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                remSc = endScIdx - idxSc;
                                if( remSc < 8)
                                {
                                    rDmrs_re1 = _mm256_setzero_ps();
                                    rDmrs_im1 = _mm256_setzero_ps();
                                    gDmrs_re1 = _mm256_setzero_ps();
                                    gDmrs_im1 = _mm256_setzero_ps();

                                    rDmrs_re2 = _mm256_setzero_ps();
                                    rDmrs_im2 = _mm256_setzero_ps();
                                    gDmrs_re2 = _mm256_setzero_ps();
                                    gDmrs_im2 = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                                        rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                                        gDmrs_re_f1[i] = genDmrs_rePtr1[i];
                                        gDmrs_im_f1[i] = genDmrs_imPtr1[i];

                                        rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                                        rDmrs_im_f2[i] = recDmrs_imPtr2[i];
                                        gDmrs_re_f2[i] = genDmrs_rePtr2[i];
                                        gDmrs_im_f2[i] = genDmrs_imPtr2[i];
                                    }

                                    bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                    ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                    out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                    out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                    // out1 = ( Y1*conj(Xdmrs1) )/beta_DMRS;
                                    out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                    out_im1 = _mm256_mul_ps(out_im1, dmrsScale);
                                    bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                    ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                    out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                    out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                    // out2 = ( Y2*conj(Xdmrs2) )/beta_DMRS;
                                    out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                    out_im2 = _mm256_mul_ps(out_im2, dmrsScale);
                                    out_re = _mm256_add_ps(out_re1, out_re2);
                                    out_im = _mm256_add_ps(out_im1, out_im2);
                                    out_re = _mm256_mul_ps(out_re, div2Vec);
                                    out_im = _mm256_mul_ps(out_im, div2Vec);

                                    for(i = 0; i < remSc; i++)
                                    {
                                        chEst_rePtr[2*i] = out_re_f[i];
                                        chEst_rePtr[2*i+1] = 0;

                                        chEst_imPtr[2*i] = out_im_f[i];
                                        chEst_imPtr[2*i+1] = 0;
                                    }
                                }
                                else
                                {
                                    rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                                    rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                                    gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                                    gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);

                                    rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                                    rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                                    gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                                    gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);

                                    bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                    ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                    out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                    out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                    // out1 = ( Y1*conj(Xdmrs1) )/beta_DMRS;
                                    out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                    out_im1 = _mm256_mul_ps(out_im1, dmrsScale);
                                    bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                    ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                    out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                    out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                    // out2 = ( Y2*conj(Xdmrs2) )/beta_DMRS;
                                    out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                    out_im2 = _mm256_mul_ps(out_im2, dmrsScale);
                                    out_re = _mm256_add_ps(out_re1, out_re2);
                                    out_im = _mm256_add_ps(out_im1, out_im2);
                                    out_re = _mm256_mul_ps(out_re, div2Vec);
                                    out_im = _mm256_mul_ps(out_im, div2Vec);

                                    temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                    temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                    temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                    temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                    temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                    temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                    temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                    temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                    _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                    _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                    _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                    _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                }
                                recDmrs_rePtr1 += 8;
                                recDmrs_imPtr1 += 8;

                                genDmrs_rePtr1 += 8;
                                genDmrs_imPtr1 += 8;

                                recDmrs_rePtr2 += 8;
                                recDmrs_imPtr2 += 8;

                                genDmrs_rePtr2 += 8;
                                genDmrs_imPtr2 += 8;

                                chEst_rePtr += 16;
                                chEst_imPtr += 16;
                            }

                            //Frequency Interpolation.
                            //Calculate estimates for null tones using Left and Right estimates.
                            //Inp : [0 a 0 b 0 c ...]
                            //Out : [(a+(a+b)/2)/2 a (a+b)/2 b (b+c)/2 c ...]
                            chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                            chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];
                            AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                            AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                            for(idxSc = toneStartIdx+1; idxSc < toneEndIdx; idxSc += 8)
                            {
                                remSc = toneEndIdx - idxSc;
                                if(remSc<8)
                                {
                                    for(i = 0; i < remSc; i += 2)
                                    {
                                        chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                        chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                    }

                                    pNrPuschOutParams->chEst_re[idxAnt][idxLyr][0] = \
                                    (pNrPuschOutParams->chEst_re[idxAnt][idxLyr][1]+\
                                    pNrPuschOutParams->chEst_re[idxAnt][idxLyr][2])*0.5;

                                    pNrPuschOutParams->chEst_im[idxAnt][idxLyr][0] = \
                                    (pNrPuschOutParams->chEst_im[idxAnt][idxLyr][1]+\
                                    pNrPuschOutParams->chEst_im[idxAnt][idxLyr][2])*0.5;
                                }
                                else
                                {
                                    AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                    AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                    temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                    temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                    temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                    temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                    temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                    temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                    temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                    temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                    temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                    temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                    temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                    temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                    AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                    AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                    _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                    _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);


                                    AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                    AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                    chEst_rePtr += 8;
                                    chEst_imPtr += 8;
                                }
                            }
                        }
                        //port 3 processing
                        else if(pNrPuschInParams->nPortIndex[idxLyr] == 3)
                        {
                            //Received DMRS Signal 1.
                            recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC + startScIdx];
                            recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                            //DMRS Signal 1.
                            genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                            genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                            //Received DMRS Signal 2.
                            recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][MAX_NUM_SC + startScIdx];
                            recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                            //DMRS Signal 2.
                            genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                            genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];
                            //Channel Estimates after Linear Interpolation.
                            chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                            chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                
                                if( (endScIdx - idxSc) < 8)
                                {
                                    remSc = endScIdx - idxSc;

                                    rDmrs_re1 = _mm256_setzero_ps();
                                    rDmrs_im1 = _mm256_setzero_ps();
                                    gDmrs_re1 = _mm256_setzero_ps();
                                    gDmrs_im1 = _mm256_setzero_ps();

                                    rDmrs_re2 = _mm256_setzero_ps();
                                    rDmrs_im2 = _mm256_setzero_ps();
                                    gDmrs_re2 = _mm256_setzero_ps();
                                    gDmrs_im2 = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                                        rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                                        gDmrs_re_f1[i] = alternateSigns_f[i]*genDmrs_rePtr1[i];
                                        gDmrs_im_f1[i] = alternateSigns_f[i]*genDmrs_imPtr1[i];

                                        rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                                        rDmrs_im_f2[i] = recDmrs_imPtr2[i];
                                        gDmrs_re_f2[i] = alternateSigns_f[i]*genDmrs_rePtr2[i];
                                        gDmrs_im_f2[i] = alternateSigns_f[i]*genDmrs_imPtr2[i];
                                    }

                                    bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                    ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                    out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                    out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                    // out1 = ( Y1*conj(Xdmrs1) )/beta_DMRS;
                                    out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                    out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                    bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                    ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                    out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                    out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                    // out2 = ( Y2*conj(Xdmrs2) )/beta_DMRS;
                                    out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                    out_im2 = _mm256_mul_ps(out_im2, dmrsScale);
                                    
                                    //Average estimates from two adjacent symbols.
                                    out_re = _mm256_add_ps(out_re1, out_re2);
                                    out_im = _mm256_add_ps(out_im1, out_im2);
                                    out_re = _mm256_mul_ps(out_re, div2Vec);
                                    out_im = _mm256_mul_ps(out_im, div2Vec);

                                    // Pad zeros between estimates for Linear Interpolation i.e., out = [0 h0 0 h1 0 h2...].
                                    for(i = 0; i < remSc; i++)
                                    {
                                        chEst_rePtr[2*i] = out_re_f[i];
                                        chEst_rePtr[2*i+1] = 0;

                                        chEst_imPtr[2*i] = out_im_f[i];
                                        chEst_imPtr[2*i+1] = 0;
                                    }
                                }
                                else
                                {
                                    rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                                    rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                                    //Reconstruct DMRS signal 1 for port 3.
                                    gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                                    gDmrs_re1 = _mm256_mul_ps(gDmrs_re1, alternateSigns);
                                    gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);
                                    gDmrs_im1 = _mm256_mul_ps(gDmrs_im1, alternateSigns);

                                    rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                                    rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                                    //Reconstruct DMRS signal 2 for port 3.
                                    gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                                    gDmrs_re2 = _mm256_mul_ps(gDmrs_re2, alternateSigns);
                                    gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);
                                    gDmrs_im2 = _mm256_mul_ps(gDmrs_im2, alternateSigns);

                                    bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                    ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                    out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                    out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                    // out1 = ( Y1*conj(Xdmrs1) )/beta_DMRS;
                                    out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                    out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                    bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                    ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                    out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                    out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                    // out2 = ( Y2*conj(Xdmrs2) )/beta_DMRS;
                                    out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                    out_im2 = _mm256_mul_ps(out_im2, dmrsScale);
                                    
                                    //Average estimates from two adjacent symbols.
                                    out_re = _mm256_add_ps(out_re1, out_re2);
                                    out_im = _mm256_add_ps(out_im1, out_im2);
                                    out_re = _mm256_mul_ps(out_re, div2Vec);
                                    out_im = _mm256_mul_ps(out_im, div2Vec);

                                    // Pad zeros between estimates for Linear Interpolation i.e., out = [0 h0 0 h1 0 h2...].
                                    temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                    temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                    temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                    temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                    temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                    temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                    temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                    temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                    _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                    _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                    _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                    _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                }

                                chEst_rePtr += 16;
                                chEst_imPtr += 16;
                                
                                recDmrs_rePtr1 += 8;
                                recDmrs_imPtr1 += 8;

                                genDmrs_rePtr1 += 8;
                                genDmrs_imPtr1 += 8;

                                recDmrs_rePtr2 += 8;
                                recDmrs_imPtr2 += 8;

                                genDmrs_rePtr2 += 8;
                                genDmrs_imPtr2 += 8;
                            }
                            
                            //Frequency Interpolation
                            //Calculate estimates for null tones using Left and Right estimates.
                            //For dead end null tone, its two adjacent estimates will be used.
                            //Inp : [0 a 0 b 0 c ...]
                            //Out : [(a+(a+b)/2)/2 a (a+b)/2 b (b+c)/2 c ...]
                            chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                            chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];
                            AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                            AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                            for(idxSc = toneStartIdx+1; idxSc < toneEndIdx; idxSc += 8)
                            {
                                remSc = toneEndIdx - idxSc;
                                if(remSc<8)
                                {
                                    for(i = 0; i < remSc; i += 2)
                                    {
                                        chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                        chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                    }

                                    pNrPuschOutParams->chEst_re[idxAnt][idxLyr][0] = \
                                    (pNrPuschOutParams->chEst_re[idxAnt][idxLyr][1]+\
                                    pNrPuschOutParams->chEst_re[idxAnt][idxLyr][2])*0.5;

                                    pNrPuschOutParams->chEst_im[idxAnt][idxLyr][0] = \
                                    (pNrPuschOutParams->chEst_im[idxAnt][idxLyr][1]+\
                                    pNrPuschOutParams->chEst_im[idxAnt][idxLyr][2])*0.5;
                                }
                                else
                                {
                                    AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                    AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                    temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                    temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                    temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                    temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                    temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                    temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                    temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                    temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                    temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                    temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                    temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                    temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                    AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                    AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                    _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                    _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);


                                    AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                    AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                    chEst_rePtr += 8;
                                    chEst_imPtr += 8;
                                }
                            }
                        }
                        //port 6 processing
                        else if(pNrPuschInParams->nPortIndex[idxLyr] == 6)
                        {
                            //Received DMRS Signal 1.
                            recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC + startScIdx];
                            recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                            //DMRS Signal 1.
                            genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                            genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                            //Received DMRS Signal 2.
                            recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][MAX_NUM_SC + startScIdx];
                            recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                            //DMRS Signal 2.
                            genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                            genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];
                            //Channel Estimates after Linear Interpolation.
                            chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                            chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                remSc = (endScIdx) - idxSc;
                                if( remSc < 8)
                                {
                                    rDmrs_re1 = _mm256_setzero_ps();
                                    rDmrs_im1 = _mm256_setzero_ps();
                                    gDmrs_re1 = _mm256_setzero_ps();
                                    gDmrs_im1 = _mm256_setzero_ps();

                                    rDmrs_re2 = _mm256_setzero_ps();
                                    rDmrs_im2 = _mm256_setzero_ps();
                                    gDmrs_re2 = _mm256_setzero_ps();
                                    gDmrs_im2 = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                                        rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                                        gDmrs_re_f1[i] = genDmrs_rePtr1[i];
                                        gDmrs_im_f1[i] = genDmrs_imPtr1[i];

                                        rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                                        rDmrs_im_f2[i] = recDmrs_imPtr2[i];
                                        gDmrs_re_f2[i] = genDmrs_rePtr2[i];
                                        gDmrs_im_f2[i] = genDmrs_imPtr2[i];
                                    }

                                    bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                    ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                    out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                    out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                    //out1 = ( Y1*conj(Xdmrs1) )/beta_DMRS;
                                    out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                    out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                    bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                    ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                    out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                    out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                    //out2 = ( Y2*conj(Xdmrs2) )/beta_DMRS;
                                    out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                    out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                    //Average estimates from two adjacent symbols.
                                    out_re = _mm256_sub_ps(out_re1, out_re2);
                                    out_im = _mm256_sub_ps(out_im1, out_im2);
                                    out_re = _mm256_mul_ps(out_re, div2Vec);
                                    out_im = _mm256_mul_ps(out_im, div2Vec);

                                    //Pad zeros between estimates for Linear Interpolation i.e., out = [0 h0 0 h1 0 h2 ...].
                                    for(i = 0; i < remSc; i++)
                                    {
                                        chEst_rePtr[2*i] = out_re_f[i];
                                        chEst_rePtr[2*i+1] = 0;

                                        chEst_imPtr[2*i] = out_im_f[i];
                                        chEst_imPtr[2*i+1] = 0;
                                    }
                                }
                                else
                                {
                                    rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                                    rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                                    gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                                    gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);

                                    rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                                    rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                                    gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                                    gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);

                                    bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                    ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                    out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                    out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                    //out1 = ( Y1*conj(Xdmrs1) )/beta_DMRS;
                                    out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                    out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                    bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                    ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                    out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                    out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                    //out2 = ( Y2*conj(Xdmrs2) )/beta_DMRS;
                                    out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                    out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                    //Average estimates from two adjacent symbols.
                                    out_re = _mm256_sub_ps(out_re1, out_re2);
                                    out_im = _mm256_sub_ps(out_im1, out_im2);
                                    out_re = _mm256_mul_ps(out_re, div2Vec);
                                    out_im = _mm256_mul_ps(out_im, div2Vec);

                                    //Pad zeros between estimates for Linear Interpolation i.e., out = [0 h0 0 h1 0 h2 ...].
                                    temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                    temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                    temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                    temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                    temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                    temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                    temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                    temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                    _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                    _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                    _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                    _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                }

                                recDmrs_rePtr1 += 8;
                                recDmrs_imPtr1 += 8;

                                genDmrs_rePtr1 += 8;
                                genDmrs_imPtr1 += 8;

                                recDmrs_rePtr2 += 8;
                                recDmrs_imPtr2 += 8;

                                genDmrs_rePtr2 += 8;
                                genDmrs_imPtr2 += 8;

                                chEst_rePtr += 16;
                                chEst_imPtr += 16;
                            }

                            //Frequency Interpolation
                            //Calculate estimates for null tones using Left and Right estimates.
                            //Inp : [0 a 0 b 0 c ...]
                            //Out : [(a+(a+b)/2)/2 a (a+b)/2 b (b+c)/2 c ...]
                            chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                            chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                            AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                            AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                            for(idxSc = toneStartIdx+1; idxSc < toneEndIdx; idxSc += 8)
                            {
                                remSc = toneEndIdx - idxSc;
                                if(remSc<8)
                                {
                                    for(i = 0; i < remSc; i += 2)
                                    {
                                        chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                        chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                    }

                                    pNrPuschOutParams->chEst_re[idxAnt][idxLyr][0] = \
                                    (pNrPuschOutParams->chEst_re[idxAnt][idxLyr][1]+\
                                    pNrPuschOutParams->chEst_re[idxAnt][idxLyr][2])*0.5;

                                    pNrPuschOutParams->chEst_im[idxAnt][idxLyr][0] = \
                                    (pNrPuschOutParams->chEst_im[idxAnt][idxLyr][1]+\
                                    pNrPuschOutParams->chEst_im[idxAnt][idxLyr][2])*0.5;
                                }
                                else
                                {
                                    AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                    AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                    temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                    temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                    temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                    temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                    temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                    temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                    temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                    temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                    temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                    temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                    temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                    temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                    AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                    AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                    _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                    _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                    AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                    AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                    chEst_rePtr += 8;
                                    chEst_imPtr += 8;
                                }
                            }
                        }
                        //port 7 processing
                        else if(pNrPuschInParams->nPortIndex[idxLyr] == 7)
                        {
                            //Received DMRS Signal 1.
                            recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC + startScIdx];
                            recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                            //DMRS Signal 1.
                            genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                            genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                            //Received DMRS Signal 2.
                            recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][MAX_NUM_SC + startScIdx];
                            recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                            //DMRS Signal 2.
                            genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                            genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];
                            //Channel Estimates after Linear Interpolation.
                            chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                            chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                remSc = endScIdx - idxSc;
                                if( remSc < 8)
                                {
                                    rDmrs_re1 = _mm256_setzero_ps();
                                    rDmrs_im1 = _mm256_setzero_ps();
                                    gDmrs_re1 = _mm256_setzero_ps();
                                    gDmrs_im1 = _mm256_setzero_ps();

                                    rDmrs_re2 = _mm256_setzero_ps();
                                    rDmrs_im2 = _mm256_setzero_ps();
                                    gDmrs_re2 = _mm256_setzero_ps();
                                    gDmrs_im2 = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                                        rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                                        gDmrs_re_f1[i] = alternateSigns_f[i]*genDmrs_rePtr1[i];
                                        gDmrs_im_f1[i] = alternateSigns_f[i]*genDmrs_imPtr1[i];

                                        rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                                        rDmrs_im_f2[i] = recDmrs_imPtr2[i];
                                        gDmrs_re_f2[i] = alternateSigns_f[i]*genDmrs_rePtr2[i];
                                        gDmrs_im_f2[i] = alternateSigns_f[i]*genDmrs_imPtr2[i];
                                    }

                                    bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                    ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                    out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                    out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                    // out = ( Y*conj(Xdmrs) )/beta_DMRS;
                                    out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                    out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                    bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                    ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                    out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                    out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                    // out = ( Y*conj(Xdmrs) )/beta_DMRS;
                                    out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                    out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                    //Average estimates from two adjacent symbols.
                                    out_re = _mm256_sub_ps(out_re1, out_re2);
                                    out_im = _mm256_sub_ps(out_im1, out_im2);
                                    out_re = _mm256_mul_ps(out_re, div2Vec);
                                    out_im = _mm256_mul_ps(out_im, div2Vec);

                                    // Pad zeros between estimates for Linear Interpolation i.e., out = [0 h0 0 h1 0 h2 ...].
                                    for(i = 0; i < remSc; i++)
                                    {
                                        chEst_rePtr[2*i] = out_re_f[i];
                                        chEst_rePtr[2*i+1] = 0;

                                        chEst_imPtr[2*i] = out_im_f[i];
                                        chEst_imPtr[2*i+1] = 0;
                                    }
                                }
                                else
                                {
                                    rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                                    rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                                    gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                                    gDmrs_re1 = _mm256_mul_ps(gDmrs_re1, alternateSigns);
                                    gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);
                                    gDmrs_im1 = _mm256_mul_ps(gDmrs_im1, alternateSigns);

                                    rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                                    rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                                    gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                                    gDmrs_re2 = _mm256_mul_ps(gDmrs_re2, alternateSigns);
                                    gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);
                                    gDmrs_im2 = _mm256_mul_ps(gDmrs_im2, alternateSigns);

                                    bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                    ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                    out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                    out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                    // out = ( Y*conj(Xdmrs) )/beta_DMRS;
                                    out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                    out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                    bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                    ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                    out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                    out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                    // out = ( Y*conj(Xdmrs) )/beta_DMRS;
                                    out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                    out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                    //Average estimates from two adjacent symbols.
                                    out_re = _mm256_sub_ps(out_re1, out_re2);
                                    out_im = _mm256_sub_ps(out_im1, out_im2);
                                    out_re = _mm256_mul_ps(out_re, div2Vec);
                                    out_im = _mm256_mul_ps(out_im, div2Vec);

                                    // Pad zeros between estimates for Linear Interpolation i.e., out = [0 h0 0 h1 0 h2 ...].
                                    temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                    temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                    temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                    temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                    temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                    temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                    temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                    temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                    _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                    _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                    _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                    _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                }

                                recDmrs_rePtr1 += 8;
                                recDmrs_imPtr1 += 8;

                                genDmrs_rePtr1 += 8;
                                genDmrs_imPtr1 += 8;

                                recDmrs_rePtr2 += 8;
                                recDmrs_imPtr2 += 8;

                                genDmrs_rePtr2 += 8;
                                genDmrs_imPtr2 += 8;

                                chEst_rePtr += 16;
                                chEst_imPtr += 16;
                            }

                            //Frequency Interpolation
                            //Calculate estimates for null tones using Left and Right estimates.
                            //Inp : [0 a 0 b 0 c ...]
                            //Out : [(a+(a+b)/2)/2 a (a+b)/2 b (b+c)/2 c ...]

                            chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                            chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                            AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                            AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                            for(idxSc = toneStartIdx+1; idxSc < toneEndIdx; idxSc += 8)
                            {
                                remSc = toneEndIdx - idxSc;
                                if( remSc <8)
                                {
                                    for(i = 0; i < remSc; i += 2)
                                    {
                                        chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                        chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                    }

                                    pNrPuschOutParams->chEst_re[idxAnt][idxLyr][0] = \
                                    (pNrPuschOutParams->chEst_re[idxAnt][idxLyr][1]+\
                                    pNrPuschOutParams->chEst_re[idxAnt][idxLyr][2])*0.5;

                                    pNrPuschOutParams->chEst_im[idxAnt][idxLyr][0] = \
                                    (pNrPuschOutParams->chEst_im[idxAnt][idxLyr][1]+\
                                    pNrPuschOutParams->chEst_im[idxAnt][idxLyr][2])*0.5;
                                }
                                else
                                {
                                    AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                    AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                    temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                    temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                    temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                    temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                    temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                    temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                    temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                    temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                    temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                    temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                    temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                    temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                    AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                    AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                    _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                    _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                    AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                    AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                    chEst_rePtr += 8;
                                    chEst_imPtr += 8;
                                }
                            }
                        }
                    
                    }
                }
            }

        }
    }

    if(cdmTotal[1]>=2)
    {
        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
            {
                startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
                toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);

                //Received DMRS Signal 1.
                recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC + startScIdx];
                recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                //DMRS Signal 1.
                genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                //Received DMRS Signal 2.
                recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][MAX_NUM_SC + startScIdx];
                recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                //DMRS Signal 2.
                genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];

                cdm1Est_rePtr1 = &cdm1EstRe[idxAnt][0][startScIdx];
                cdm1Est_imPtr1 = &cdm1EstIm[idxAnt][0][startScIdx];

                cdm1Est_rePtr2 = &cdm1EstRe[idxAnt][1][startScIdx];
                cdm1Est_imPtr2 = &cdm1EstIm[idxAnt][1][startScIdx];

                for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                {
                    remSc = endScIdx - idxSc;
                    if( remSc < 8)
                    {
                        rDmrs_re1 = _mm256_setzero_ps();
                        rDmrs_im1 = _mm256_setzero_ps();
                        rDmrs_re2 = _mm256_setzero_ps();
                        rDmrs_im2 = _mm256_setzero_ps();

                        gDmrs_re1 = _mm256_setzero_ps();
                        gDmrs_im1 = _mm256_setzero_ps();
                        gDmrs_re2 = _mm256_setzero_ps();
                        gDmrs_im2 = _mm256_setzero_ps();

                        for(i = 0; i < remSc; i++)
                        {
                            rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                            rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                            rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                            rDmrs_im_f2[i] = recDmrs_imPtr2[i];

                            gDmrs_re_f1[i] = genDmrs_rePtr1[i];
                            gDmrs_im_f1[i] = genDmrs_imPtr1[i];
                            gDmrs_re_f2[i] = genDmrs_rePtr2[i];
                            gDmrs_im_f2[i] = genDmrs_imPtr2[i];
                        }

                        bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                        ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);

                        out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                        out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);

                        // out1 = ( Y1*conj(Xdmrs1) )/Beta_dmrs;
                        out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                        out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                        bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                        ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);

                        out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                        out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);

                        // out2 = ( Y2*conj(Xdmrs2) )/Beta_dmrs;
                        out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                        out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                        for(i = 0; i < remSc; i++)
                        {
                            cdm1Est_rePtr1[i] = out_re1_f[i];
                            cdm1Est_imPtr1[i] = out_im1_f[i];
                            cdm1Est_rePtr2[i] = out_re2_f[i];
                            cdm1Est_imPtr2[i] = out_im2_f[i];
                        }
                    }
                    else
                    {
                        rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                        rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                        gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                        gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);

                        rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                        rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                        gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                        gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);

                        bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                        ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);

                        out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                        out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);

                        // out1 = ( Y1*conj(Xdmrs1) )/Beta_dmrs;
                        out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                        out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                        bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                        ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);

                        out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                        out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);

                        // out2 = ( Y2*conj(Xdmrs2) )/Beta_dmrs;
                        out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                        out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                        _mm256_storeu_ps(cdm1Est_rePtr1, out_re1);
                        _mm256_storeu_ps(cdm1Est_imPtr1, out_im1);

                        _mm256_storeu_ps(cdm1Est_rePtr2, out_re2);
                        _mm256_storeu_ps(cdm1Est_imPtr2, out_im2);
                    }
                    cdm1Est_rePtr1 += 8;
                    cdm1Est_imPtr1 += 8;

                    cdm1Est_rePtr2 += 8;
                    cdm1Est_imPtr2 += 8;

                    recDmrs_rePtr1 += 8;
                    recDmrs_imPtr1 += 8;

                    recDmrs_rePtr2 += 8;
                    recDmrs_imPtr2 += 8;

                    genDmrs_rePtr1 += 8;
                    genDmrs_imPtr1 += 8;

                    genDmrs_rePtr2 += 8;
                    genDmrs_imPtr2 += 8;
                }
            }
        }

        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
            {
                for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                {
                    startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                    endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                    toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
                    toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);

                    cdm1Est_rePtr1 = &cdm1EstRe[idxAnt][0][startScIdx];
                    cdm1Est_imPtr1 = &cdm1EstIm[idxAnt][0][startScIdx];

                    cdm1Est_rePtr2 = &cdm1EstRe[idxAnt][1][startScIdx];
                    cdm1Est_imPtr2 = &cdm1EstIm[idxAnt][1][startScIdx];

                    //Channel Estimates after Linear Interpolation
                    chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                    chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                    if(pNrPuschInParams->nPortIndex[idxLyr] == 2)
                    {
                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc <= 8)
                            {
                                for(i = 0; i < remSc-1; i++)
                                {
                                    chEst_rePtr[2*i] = ((cdm1Est_rePtr1[i] + cdm1Est_rePtr1[i+1]) + (cdm1Est_rePtr2[i]+cdm1Est_rePtr2[i+1]))*0.25;
                                    chEst_rePtr[2*i+1] = 0;
                                    chEst_imPtr[2*i] = ((cdm1Est_imPtr1[i] + cdm1Est_imPtr1[i+1]) + (cdm1Est_imPtr2[i]+cdm1Est_imPtr2[i+1]))*0.25;
                                    chEst_imPtr[2*i+1] = 0;
                                }
                                chEst_rePtr[2*(remSc-1)] = chEst_rePtr[2*(remSc-2)];
                                chEst_imPtr[2*(remSc-1)] = chEst_imPtr[2*(remSc-2)];
                            }
                            else
                            {
                                cdm1EstRe1 = _mm256_loadu_ps(cdm1Est_rePtr1);
                                cdm1EstRe2 = _mm256_loadu_ps(cdm1Est_rePtr1+1);
                                cdm1EstIm1 = _mm256_loadu_ps(cdm1Est_imPtr1);
                                cdm1EstIm2 = _mm256_loadu_ps(cdm1Est_imPtr1+1);

                                out_re1 = _mm256_add_ps(cdm1EstRe1, cdm1EstRe2);
                                out_im1 = _mm256_add_ps(cdm1EstIm1, cdm1EstIm2);
                                out_re1 = _mm256_mul_ps(out_re1, div2Vec);
                                out_im1 = _mm256_mul_ps(out_im1, div2Vec);

                                cdm1EstRe1 = _mm256_loadu_ps(cdm1Est_rePtr2);
                                cdm1EstRe2 = _mm256_loadu_ps(cdm1Est_rePtr2+1);
                                cdm1EstIm1 = _mm256_loadu_ps(cdm1Est_imPtr2);
                                cdm1EstIm2 = _mm256_loadu_ps(cdm1Est_imPtr2+1);

                                out_re2 = _mm256_add_ps(cdm1EstRe1, cdm1EstRe2);
                                out_im2 = _mm256_add_ps(cdm1EstIm1, cdm1EstIm2);
                                out_re2 = _mm256_mul_ps(out_re2, div2Vec);
                                out_im2 = _mm256_mul_ps(out_im2, div2Vec);

                                // Mean of two adjacent estimates.
                                out_re = _mm256_add_ps(out_re1, out_re2);
                                out_im = _mm256_add_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                
                            }
                            cdm1Est_rePtr1 += 8;
                            cdm1Est_imPtr1 += 8;

                            cdm1Est_rePtr2 += 8;
                            cdm1Est_imPtr2 += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }
                    }
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 3)
                    {
                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc <= 8)
                            {
                                for(i = 0; i < remSc-1; i++)
                                {
                                    chEst_rePtr[2*i] = alternateSigns_f[i]*((cdm1Est_rePtr1[i] - cdm1Est_rePtr1[i+1]) + (cdm1Est_rePtr2[i] - cdm1Est_rePtr2[i+1]))*0.25;
                                    chEst_rePtr[2*i+1] = 0;
                                    chEst_imPtr[2*i] = alternateSigns_f[i]*((cdm1Est_imPtr1[i] - cdm1Est_imPtr1[i+1]) + (cdm1Est_imPtr2[i] - cdm1Est_imPtr2[i+1]))*0.25;
                                    chEst_imPtr[2*i+1] = 0;
                                }
                                chEst_rePtr[2*(remSc-1)] = chEst_rePtr[2*(remSc-2)];
                                chEst_imPtr[2*(remSc-1)] = chEst_imPtr[2*(remSc-2)];
                            }
                            else
                            {
                                cdm1EstRe1 = _mm256_loadu_ps(cdm1Est_rePtr1);
                                cdm1EstRe2 = _mm256_loadu_ps(cdm1Est_rePtr1+1);
                                cdm1EstIm1 = _mm256_loadu_ps(cdm1Est_imPtr1);
                                cdm1EstIm2 = _mm256_loadu_ps(cdm1Est_imPtr1+1);

                                out_re1 = _mm256_sub_ps(cdm1EstRe1, cdm1EstRe2);
                                out_re1 = _mm256_mul_ps(out_re1, alternateSigns);
                                out_im1 = _mm256_sub_ps(cdm1EstIm1, cdm1EstIm2);
                                out_im1 = _mm256_mul_ps(out_im1, alternateSigns);
                                out_re1 = _mm256_mul_ps(out_re1, div2Vec);
                                out_im1 = _mm256_mul_ps(out_im1, div2Vec);

                                cdm1EstRe1 = _mm256_loadu_ps(cdm1Est_rePtr2);
                                cdm1EstRe2 = _mm256_loadu_ps(cdm1Est_rePtr2+1);
                                cdm1EstIm1 = _mm256_loadu_ps(cdm1Est_imPtr2);
                                cdm1EstIm2 = _mm256_loadu_ps(cdm1Est_imPtr2+1);

                                out_re2 = _mm256_sub_ps(cdm1EstRe1, cdm1EstRe2);
                                out_re2 = _mm256_mul_ps(out_re2, alternateSigns);
                                out_im2 = _mm256_sub_ps(cdm1EstIm1, cdm1EstIm2);
                                out_im2 = _mm256_mul_ps(out_im2, alternateSigns);
                                out_re2 = _mm256_mul_ps(out_re2, div2Vec);
                                out_im2 = _mm256_mul_ps(out_im2, div2Vec);

                                // Mean of two adjacent estimates.
                                out_re = _mm256_add_ps(out_re1, out_re2);
                                out_im = _mm256_add_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                
                            }
                            cdm1Est_rePtr1 += 8;
                            cdm1Est_imPtr1 += 8;

                            cdm1Est_rePtr2 += 8;
                            cdm1Est_imPtr2 += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }
                    }
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 6)
                    {
                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc <= 8)
                            {
                                for(i = 0; i < remSc-1; i++)
                                {
                                    chEst_rePtr[2*i] = ((cdm1Est_rePtr1[i] + cdm1Est_rePtr1[i+1]) - (cdm1Est_rePtr2[i] + cdm1Est_rePtr2[i+1]))*0.25;
                                    chEst_rePtr[2*i+1] = 0;
                                    chEst_imPtr[2*i] = ((cdm1Est_imPtr1[i] + cdm1Est_imPtr1[i+1]) - (cdm1Est_imPtr2[i] + cdm1Est_imPtr2[i+1]))*0.25;
                                    chEst_imPtr[2*i+1] = 0;
                                }
                                chEst_rePtr[2*(remSc-1)] = chEst_rePtr[2*(remSc-2)];
                                chEst_imPtr[2*(remSc-1)] = chEst_imPtr[2*(remSc-2)];
                            }
                            else
                            {
                                cdm1EstRe1 = _mm256_loadu_ps(cdm1Est_rePtr1);
                                cdm1EstRe2 = _mm256_loadu_ps(cdm1Est_rePtr1+1);
                                cdm1EstIm1 = _mm256_loadu_ps(cdm1Est_imPtr1);
                                cdm1EstIm2 = _mm256_loadu_ps(cdm1Est_imPtr1+1);

                                out_re1 = _mm256_add_ps(cdm1EstRe1, cdm1EstRe2);
                                out_im1 = _mm256_add_ps(cdm1EstIm1, cdm1EstIm2);
                                out_re1 = _mm256_mul_ps(out_re1, div2Vec);
                                out_im1 = _mm256_mul_ps(out_im1, div2Vec);

                                cdm1EstRe1 = _mm256_loadu_ps(cdm1Est_rePtr2);
                                cdm1EstRe2 = _mm256_loadu_ps(cdm1Est_rePtr2+1);
                                cdm1EstIm1 = _mm256_loadu_ps(cdm1Est_imPtr2);
                                cdm1EstIm2 = _mm256_loadu_ps(cdm1Est_imPtr2+1);

                                out_re2 = _mm256_add_ps(cdm1EstRe1, cdm1EstRe2);
                                out_im2 = _mm256_add_ps(cdm1EstIm1, cdm1EstIm2);
                                out_re2 = _mm256_mul_ps(out_re2, div2Vec);
                                out_im2 = _mm256_mul_ps(out_im2, div2Vec);

                                // Mean of two adjacent estimates.
                                out_re = _mm256_sub_ps(out_re1, out_re2);
                                out_im = _mm256_sub_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                
                            }
                            cdm1Est_rePtr1 += 8;
                            cdm1Est_imPtr1 += 8;

                            cdm1Est_rePtr2 += 8;
                            cdm1Est_imPtr2 += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }
                    }
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 7)
                    {
                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc <= 8)
                            {
                                for(i = 0; i < remSc-1; i++)
                                {
                                    chEst_rePtr[2*i] = alternateSigns_f[i]*((cdm1Est_rePtr1[i] - cdm1Est_rePtr1[i+1]) - (cdm1Est_rePtr2[i] - cdm1Est_rePtr2[i+1]))*0.25;
                                    chEst_rePtr[2*i+1] = 0;
                                    chEst_imPtr[2*i] = alternateSigns_f[i]*((cdm1Est_imPtr1[i] - cdm1Est_imPtr1[i+1]) - (cdm1Est_imPtr2[i] - cdm1Est_imPtr2[i+1]))*0.25;
                                    chEst_imPtr[2*i+1] = 0;
                                }
                                chEst_rePtr[2*(remSc-1)] = chEst_rePtr[2*(remSc-2)];
                                chEst_imPtr[2*(remSc-1)] = chEst_imPtr[2*(remSc-2)];                                
                            }
                            else
                            {
                                cdm1EstRe1 = _mm256_loadu_ps(cdm1Est_rePtr1);
                                cdm1EstRe2 = _mm256_loadu_ps(cdm1Est_rePtr1+1);
                                cdm1EstIm1 = _mm256_loadu_ps(cdm1Est_imPtr1);
                                cdm1EstIm2 = _mm256_loadu_ps(cdm1Est_imPtr1+1);

                                out_re1 = _mm256_sub_ps(cdm1EstRe1, cdm1EstRe2);
                                out_re1 = _mm256_mul_ps(out_re1, alternateSigns);
                                out_im1 = _mm256_sub_ps(cdm1EstIm1, cdm1EstIm2);
                                out_im1 = _mm256_mul_ps(out_im1, alternateSigns);
                                out_re1 = _mm256_mul_ps(out_re1, div2Vec);
                                out_im1 = _mm256_mul_ps(out_im1, div2Vec);

                                cdm1EstRe1 = _mm256_loadu_ps(cdm1Est_rePtr2);
                                cdm1EstRe2 = _mm256_loadu_ps(cdm1Est_rePtr2+1);
                                cdm1EstIm1 = _mm256_loadu_ps(cdm1Est_imPtr2);
                                cdm1EstIm2 = _mm256_loadu_ps(cdm1Est_imPtr2+1);

                                out_re2 = _mm256_sub_ps(cdm1EstRe1, cdm1EstRe2);
                                out_re2 = _mm256_mul_ps(out_re2, alternateSigns);
                                out_im2 = _mm256_sub_ps(cdm1EstIm1, cdm1EstIm2);
                                out_im2 = _mm256_mul_ps(out_im2, alternateSigns);
                                out_re2 = _mm256_mul_ps(out_re2, div2Vec);
                                out_im2 = _mm256_mul_ps(out_im2, div2Vec);

                                // Mean of two adjacent estimates.
                                out_re = _mm256_sub_ps(out_re1, out_re2);
                                out_im = _mm256_sub_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                
                            }
                            cdm1Est_rePtr1 += 8;
                            cdm1Est_imPtr1 += 8;

                            cdm1Est_rePtr2 += 8;
                            cdm1Est_imPtr2 += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }
                    }

                    if(pNrPuschInParams->nPortIndex[idxLyr] == 2 || pNrPuschInParams->nPortIndex[idxLyr] == 3 ||\
                       pNrPuschInParams->nPortIndex[idxLyr] == 6 || pNrPuschInParams->nPortIndex[idxLyr] == 7)
                       {
                           //Frequency Interpolation.
                           //Calculate estimates for null tones using Left and Right estimates.
                           //Inp : [0 a 0 b 0 c ...]
                           //Out : [(a+(a+b)/2)/2 a (a+b)/2 b (b+c)/2 c ...]
                           chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                           chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];
                           AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                           AVX_chEst_imPtr = (__m256*)chEst_imPtr;
                           for(idxSc = toneStartIdx+1; idxSc < toneEndIdx; idxSc += 8)
                            {
                                if(toneEndIdx - idxSc<8)
                                {
                                    remSc = toneEndIdx - idxSc;

                                    for(i = 0; i < remSc; i += 2)
                                    {
                                        chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                        chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                    }

                                    pNrPuschOutParams->chEst_re[idxAnt][idxLyr][0] = \
                                    (pNrPuschOutParams->chEst_re[idxAnt][idxLyr][1]+\
                                    pNrPuschOutParams->chEst_re[idxAnt][idxLyr][2])*0.5;

                                    pNrPuschOutParams->chEst_im[idxAnt][idxLyr][0] = \
                                    (pNrPuschOutParams->chEst_im[idxAnt][idxLyr][1]+\
                                    pNrPuschOutParams->chEst_im[idxAnt][idxLyr][2])*0.5;
                                }
                                else
                                {
                                    AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                    AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                    temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                    temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                    temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                    temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                    temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                    temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                    temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                    temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                    temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                    temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                    temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                    temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                    AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                    AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                    _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                    _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);


                                    AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                    AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                    chEst_rePtr += 8;
                                    chEst_imPtr += 8;
                                }
                            }

                       }
                    
                }
            }
        }
    
        if(cdmTotal[0] == 1)
        {
            for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
            {
                for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
                {
                    for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                    {
                        startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                        endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                        toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
                        toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);

                        // port 0 processing
                        if(pNrPuschInParams->nPortIndex[idxLyr] == 0)
                        {          
                            //Received DMRS Signal 1                          
                            recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                            recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(MAX_NUM_SC/2) + startScIdx];
                            //Received DMRS Signal 2
                            recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][startScIdx];
                            recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(MAX_NUM_SC/2) + startScIdx];
                            //Generated DMRS Signal 1
                            genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                            genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                            //Generated DMRS Signal 2
                            genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                            genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];
                            //Channel Estimates after Linear Interpolation
                            chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                            chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];

                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                remSc = endScIdx - idxSc;
                                if( remSc < 8)
                                {
                                    rDmrs_re1 = _mm256_setzero_ps();
                                    rDmrs_im1 = _mm256_setzero_ps();
                                    rDmrs_re2 = _mm256_setzero_ps();
                                    rDmrs_im2 = _mm256_setzero_ps();

                                    gDmrs_re1 = _mm256_setzero_ps();
                                    gDmrs_im1 = _mm256_setzero_ps();
                                    gDmrs_re2 = _mm256_setzero_ps();
                                    gDmrs_im2 = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                                        rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                                        rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                                        rDmrs_im_f2[i] = recDmrs_imPtr2[i];

                                        gDmrs_re_f1[i] = genDmrs_rePtr1[i];
                                        gDmrs_im_f1[i] = genDmrs_imPtr1[i];
                                        gDmrs_re_f2[i] = genDmrs_rePtr2[i];
                                        gDmrs_im_f2[i] = genDmrs_imPtr2[i];
                                    }

                                    bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                    ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);

                                    out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                    out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);

                                    // out1 = ( Y1*conj(Xdmrs1) )/Beta_dmrs;
                                    out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                    out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                    bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                    ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);

                                    out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                    out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);

                                    // out1 = ( Y2*conj(Xdmrs2) )/Beta_dmrs;
                                    out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                    out_im2 = _mm256_mul_ps(out_im2, dmrsScale);
                                    
                                    // Mean of two adjacent estimates.
                                    out_re = _mm256_add_ps(out_re1, out_re2);
                                    out_im = _mm256_add_ps(out_im1, out_im2);
                                    out_re = _mm256_mul_ps(out_re, div2Vec);
                                    out_im = _mm256_mul_ps(out_im, div2Vec);

                                    for(i = 0; i < remSc; i++)
                                    {
                                        chEst_rePtr[2*i] = out_re_f[i];
                                        chEst_rePtr[2*i+1] = 0;

                                        chEst_imPtr[2*i] = out_im_f[i];
                                        chEst_imPtr[2*i+1] = 0;
                                    }
                                }
                                else
                                {
                                    rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                                    rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                                    gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                                    gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);

                                    rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                                    rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                                    gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                                    gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);

                                    bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                    ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);

                                    out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                    out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);

                                    // out1 = ( Y1*conj(Xdmrs1) )/Beta_dmrs;
                                    out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                    out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                    bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                    ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);

                                    out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                    out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);

                                    // out1 = ( Y2*conj(Xdmrs2) )/Beta_dmrs;
                                    out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                    out_im2 = _mm256_mul_ps(out_im2, dmrsScale);
                                    
                                    // Mean of two adjacent estimates.
                                    out_re = _mm256_add_ps(out_re1, out_re2);
                                    out_im = _mm256_add_ps(out_im1, out_im2);
                                    out_re = _mm256_mul_ps(out_re, div2Vec);
                                    out_im = _mm256_mul_ps(out_im, div2Vec);

                                    temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                    temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                    temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                    temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                    temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                    temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                    temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                    temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                    _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                    _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                    _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                    _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                }
                                recDmrs_rePtr1 += 8;
                                recDmrs_imPtr1 += 8;

                                recDmrs_rePtr2 += 8;
                                recDmrs_imPtr2 += 8;

                                genDmrs_rePtr1 += 8;
                                genDmrs_imPtr1 += 8;

                                genDmrs_rePtr2 += 8;
                                genDmrs_imPtr2 += 8;

                                chEst_rePtr += 16;
                                chEst_imPtr += 16;
                            }   //for idxSc

                            //Frequency Interpolation.
                            //Calculate estimates for null tones using Left and Right estimates.
                            //For dead end null tone, its two adjacent estimates will be used.
                            //Inp : [ a 0 b 0 c ...]
                            //Out : [a (a+b)/2 b (b+c)/2 c...]
                            chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                            chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];
                            AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                            AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                            for(idxSc = toneStartIdx; idxSc < toneEndIdx; idxSc += 8)
                            {
                                remSc = toneEndIdx - idxSc;
                                if(remSc<=8)
                                {
                                    for(i = 0; i < remSc; i += 2)
                                    {
                                        if(remSc - i == 2)
                                        {
                                            chEst_rePtr[i+1] = (chEst_rePtr[i]+chEst_rePtr[i-1])*0.5;
                                            chEst_imPtr[i+1] = (chEst_imPtr[i]+chEst_imPtr[i-1])*0.5;
                                        }
                                        else
                                        {
                                            chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                            chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                        }
                                    }
                                }
                                else
                                {
                                    AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                    AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                    temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                    temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                    temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                    temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                    temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                    temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                    temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                    temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                    temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                    temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                    temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                    temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                    AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                    AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                    _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                    _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                    AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                    AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                    chEst_rePtr += 8;
                                    chEst_imPtr += 8;
                                }
                            }
                        }                  
                        // port 1 processing
                        else if(pNrPuschInParams->nPortIndex[idxLyr] == 1)
                        {
                            //Recived DMRS Signal 1
                            recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                            recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(MAX_NUM_SC/2) + startScIdx];
                            //Recived DMRS Signal 2
                            recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][startScIdx];
                            recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(MAX_NUM_SC/2) + startScIdx];
                            //DMRS Signal 1.
                            genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                            genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                            //DMRS Signal 2.
                            genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                            genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];
                            //Channel Estimates after Linear Interpolation.
                            chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                            chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];

                            
                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                remSc = endScIdx - idxSc;
                                if(remSc<8)
                                {
                                    rDmrs_re1 = _mm256_setzero_ps();
                                    rDmrs_im1 = _mm256_setzero_ps();
                                    gDmrs_re1 = _mm256_setzero_ps();
                                    gDmrs_im1 = _mm256_setzero_ps();

                                    rDmrs_re2 = _mm256_setzero_ps();
                                    rDmrs_im2 = _mm256_setzero_ps();
                                    gDmrs_re2 = _mm256_setzero_ps();
                                    gDmrs_im2 = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                                        rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                                        gDmrs_re_f1[i] = alternateSigns_f[i]*genDmrs_rePtr1[i];
                                        gDmrs_im_f1[i] = alternateSigns_f[i]*genDmrs_imPtr1[i];

                                        rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                                        rDmrs_im_f2[i] = recDmrs_imPtr2[i];
                                        gDmrs_re_f2[i] = alternateSigns_f[i]*genDmrs_rePtr2[i];
                                        gDmrs_im_f2[i] = alternateSigns_f[i]*genDmrs_imPtr2[i];
                                    }

                                    bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                    ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);

                                    out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                    out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);

                                    // out1 = ( Y1*conj(Xdmrs1) )/beta_dmrs;
                                    out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                    out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                    bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                    ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);

                                    out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                    out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);

                                    // out2 = ( Y2*conj(Xdmrs2) )/beta_dmrs;
                                    out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                    out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                    out_re = _mm256_add_ps(out_re1, out_re2);
                                    out_im = _mm256_add_ps(out_im1, out_im2);
                                    out_re = _mm256_mul_ps(out_re, div2Vec);
                                    out_im = _mm256_mul_ps(out_im, div2Vec);

                                    // Pad zeros between estimates for Linear Interpolation i.e., out = [h0 0 h1 0 h2...].
                                    for(i = 0; i < remSc; i++)
                                    {
                                        chEst_rePtr[2*i] = out_re[i];
                                        chEst_rePtr[2*i+1] = 0;

                                        chEst_imPtr[2*i] = out_im[i];
                                        chEst_imPtr[2*i+1] = 0;
                                    }
                                }
                                else
                                {
                                    rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                                    rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                                    gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                                    gDmrs_re1 = _mm256_mul_ps(gDmrs_re1, alternateSigns);
                                    gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);
                                    gDmrs_im1 = _mm256_mul_ps(gDmrs_im1, alternateSigns);

                                    rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                                    rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                                    gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                                    gDmrs_re2 = _mm256_mul_ps(gDmrs_re2, alternateSigns);
                                    gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);
                                    gDmrs_im2 = _mm256_mul_ps(gDmrs_im2, alternateSigns);

                                    bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                    ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);

                                    out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                    out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);

                                    // out1 = ( Y1*conj(Xdmrs1) )/beta_dmrs;
                                    out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                    out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                    bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                    ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);

                                    out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                    out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);

                                    // out2 = ( Y2*conj(Xdmrs2) )/beta_dmrs;
                                    out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                    out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                    out_re = _mm256_add_ps(out_re1, out_re2);
                                    out_im = _mm256_add_ps(out_im1, out_im2);
                                    out_re = _mm256_mul_ps(out_re, div2Vec);
                                    out_im = _mm256_mul_ps(out_im, div2Vec);
                                    
                                    // Pad zeros between estimates for Linear Interpolation i.e., out = [h0 0 h1 0 h2...].
                                    temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                    temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                    temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                    temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                    temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                    temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                    temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                    temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                    _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                    _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                    _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                    _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                }
                                
                                
                                recDmrs_rePtr1 += 8;
                                recDmrs_imPtr1 += 8;

                                genDmrs_rePtr1 += 8;
                                genDmrs_imPtr1 += 8;

                                recDmrs_rePtr2 += 8;
                                recDmrs_imPtr2 += 8;

                                genDmrs_rePtr2 += 8;
                                genDmrs_imPtr2 += 8;

                                chEst_rePtr += 16;
                                chEst_imPtr += 16;
                            }   //for idxSc

                            //Frequency Interpolation
                            //Calculate estimates for null tones using Left and Right estimates.
                            //Inp : [ a 0 b 0 c ...]
                            //Out : [a (a+b)/2 b (b+c)/2 c ...]
                            chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                            chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];
                            AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                            AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                            for(idxSc = toneStartIdx; idxSc < toneEndIdx; idxSc += 8)
                            {
                                remSc = toneEndIdx - idxSc;
                                if(remSc<=8)
                                {
                                    for(i = 0; i < remSc; i += 2)
                                    {
                                        if(remSc - i == 2)
                                        {
                                            chEst_rePtr[i+1] = (chEst_rePtr[i]+chEst_rePtr[i-1])*0.5;
                                            chEst_imPtr[i+1] = (chEst_imPtr[i]+chEst_imPtr[i-1])*0.5;
                                        }
                                        else
                                        {
                                            chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                            chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                        }
                                    }
                                }
                                else
                                {
                                    AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                    AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                    temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                    temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                    temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                    temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                    temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                    temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                    temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                    temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                    temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                    temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                    temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                    temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                    AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                    AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                    _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                    _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                    AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                    AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                    chEst_rePtr += 8;
                                    chEst_imPtr += 8;
                                }
                            }
                        }
                        //port 4 processing
                        else if(pNrPuschInParams->nPortIndex[idxLyr] == 4)
                        {          
                            //Received DMRS Signal 1.                         
                            recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                            recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(MAX_NUM_SC/2) + startScIdx];
                            //Received DMRS Signal 2.
                            recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][startScIdx];
                            recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(MAX_NUM_SC/2) + startScIdx];
                            //DMRS Signal 1.
                            genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                            genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                            //DMRS Signal 2.
                            genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                            genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];
                            //Channel Estimates after Linear Interpolation.
                            chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                            chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];

                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                remSc = endScIdx - idxSc;
                                if( remSc < 8)
                                {
                                    rDmrs_re1 = _mm256_setzero_ps();
                                    rDmrs_im1 = _mm256_setzero_ps();
                                    rDmrs_re2 = _mm256_setzero_ps();
                                    rDmrs_im2 = _mm256_setzero_ps();

                                    gDmrs_re1 = _mm256_setzero_ps();
                                    gDmrs_im1 = _mm256_setzero_ps();
                                    gDmrs_re2 = _mm256_setzero_ps();
                                    gDmrs_im2 = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                                        rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                                        rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                                        rDmrs_im_f2[i] = recDmrs_imPtr2[i];

                                        gDmrs_re_f1[i] = genDmrs_rePtr1[i];
                                        gDmrs_im_f1[i] = genDmrs_imPtr1[i];
                                        gDmrs_re_f2[i] = genDmrs_rePtr2[i];
                                        gDmrs_im_f2[i] = genDmrs_imPtr2[i];
                                    }

                                    bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                    ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                    out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                    out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                    //out1 = ( Y1*conj(Xdmrs1) )/Beta_dmrs;
                                    out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                    out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                    bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                    ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                    out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                    out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                    //out2 = ( Y2*conj(Xdmrs2) )/Beta_dmrs;
                                    out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                    out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                    //Average estimates from two adjacent symbols.
                                    out_re = _mm256_sub_ps(out_re1, out_re2);
                                    out_im = _mm256_sub_ps(out_im1, out_im2);
                                    out_re = _mm256_mul_ps(out_re, div2Vec);
                                    out_im = _mm256_mul_ps(out_im, div2Vec);
                                    // Pad zeros between estimates for Linear Interpolation i.e., out = [h0 0 h1 0 h2 ...];
                                    for(i = 0; i < remSc; i++)
                                    {
                                        chEst_rePtr[2*i] = out_re_f[i];
                                        chEst_rePtr[2*i+1] = 0;

                                        chEst_imPtr[2*i] = out_im_f[i];
                                        chEst_imPtr[2*i+1] = 0;
                                    }
                                }
                                else
                                {
                                    rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                                    rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                                    gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                                    gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);

                                    rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                                    rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                                    gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                                    gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);

                                    bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                    ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                    out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                    out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                    //out1 = ( Y1*conj(Xdmrs1) )/Beta_dmrs;
                                    out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                    out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                    bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                    ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                    out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                    out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                    //out2 = ( Y2*conj(Xdmrs2) )/Beta_dmrs;
                                    out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                    out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                    //Average estimates from two adjacent symbols.
                                    out_re = _mm256_sub_ps(out_re1, out_re2);
                                    out_im = _mm256_sub_ps(out_im1, out_im2);
                                    out_re = _mm256_mul_ps(out_re, div2Vec);
                                    out_im = _mm256_mul_ps(out_im, div2Vec);

                                    // Pad zeros between estimates for Linear Interpolation i.e., out = [h0 0 h1 0 h2 ...];
                                    temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                    temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                    temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                    temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                    temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                    temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                    temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                    temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                    _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                    _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                    _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                    _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                }
                                
                                recDmrs_rePtr1 += 8;
                                recDmrs_imPtr1 += 8;

                                recDmrs_rePtr2 += 8;
                                recDmrs_imPtr2 += 8;

                                genDmrs_rePtr1 += 8;
                                genDmrs_imPtr1 += 8;

                                genDmrs_rePtr2 += 8;
                                genDmrs_imPtr2 += 8;

                                chEst_rePtr += 16;
                                chEst_imPtr += 16;
                            }

                            //Frequency Interpolation 
                            //Calculate estimates for null tones using Left and Right estimates.
                            //Inp : [ a 0 b 0 c ...]
                            //Out : [a (a+b)/2 b (b+c)/2 c...]
                            chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                            chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];
                            AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                            AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                            for(idxSc = toneStartIdx; idxSc < toneEndIdx; idxSc += 8)
                            {
                                remSc = toneEndIdx - idxSc;
                                if(remSc<=8)
                                {
                                    for(i = 0; i < remSc; i += 2)
                                    {
                                        if(remSc - i == 2)
                                        {
                                            chEst_rePtr[i+1] = (chEst_rePtr[i]+chEst_rePtr[i-1])*0.5;
                                            chEst_imPtr[i+1] = (chEst_imPtr[i]+chEst_imPtr[i-1])*0.5;
                                        }
                                        else
                                        {
                                            chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                            chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                        }
                                    }
                                }
                                else
                                {
                                    AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                    AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                    temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                    temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                    temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                    temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                    temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                    temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                    temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                    temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                    temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                    temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                    temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                    temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                    AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                    AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                    _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                    _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                    AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                    AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                    chEst_rePtr += 8;
                                    chEst_imPtr += 8;
                                }
                            }
                        }
                        //port 5 processing
                        else if(pNrPuschInParams->nPortIndex[idxLyr] == 5)
                        {
                            //Received DMRS signal 1.
                            recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                            recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(MAX_NUM_SC/2) + startScIdx];
                            //Received DMRS signal 2.
                            recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][startScIdx];
                            recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(MAX_NUM_SC/2) + startScIdx];
                            //DMRS signal 1.
                            genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                            genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                            //DMRS signal 2.
                            genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                            genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];
                            //Channel Estimates after Linear Interpolation.
                            chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                            chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];
                            
                            for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                            {
                                if( (endScIdx - idxSc)<8)
                                {
                                    remSc = endScIdx - idxSc;
                                    rDmrs_re1 = _mm256_setzero_ps();
                                    rDmrs_im1 = _mm256_setzero_ps();
                                    gDmrs_re1 = _mm256_setzero_ps();
                                    gDmrs_im1 = _mm256_setzero_ps();

                                    rDmrs_re2 = _mm256_setzero_ps();
                                    rDmrs_im2 = _mm256_setzero_ps();
                                    gDmrs_re2 = _mm256_setzero_ps();
                                    gDmrs_im2 = _mm256_setzero_ps();

                                    for(i = 0; i < remSc; i++)
                                    {
                                        rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                                        rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                                        gDmrs_re_f1[i] = alternateSigns_f[i]*genDmrs_rePtr1[i];
                                        gDmrs_im_f1[i] = alternateSigns_f[i]*genDmrs_imPtr1[i];

                                        rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                                        rDmrs_im_f2[i] = recDmrs_imPtr2[i];
                                        gDmrs_re_f2[i] = alternateSigns_f[i]*genDmrs_rePtr2[i];
                                        gDmrs_im_f2[i] = alternateSigns_f[i]*genDmrs_imPtr2[i];
                                    }

                                    bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                    ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                    out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                    out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                    //out1 = ( Y1*conj(Xdmrs1) )/beta_dmrs;
                                    out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                    out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                    bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                    ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                    out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                    out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                    //out2 = ( Y2*conj(Xdmrs2) )/beta_dmrs;
                                    out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                    out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                    //Average estimates from two adjacent symbols.
                                    out_re = _mm256_sub_ps(out_re1, out_re2);
                                    out_im = _mm256_sub_ps(out_im1, out_im2);
                                    out_re = _mm256_mul_ps(out_re, div2Vec);
                                    out_im = _mm256_mul_ps(out_im, div2Vec);

                                    // Pad zeros between estimates for Linear Interpolation i.e., out = [h0 0 h1 0 h2...].
                                    for(i = 0; i < remSc; i++)
                                    {
                                        chEst_rePtr[2*i] = out_re_f[i];
                                        chEst_rePtr[2*i+1] = 0;
                                        chEst_imPtr[2*i] = out_im_f[i];
                                        chEst_imPtr[2*i+1] = 0;
                                    }
                                }
                                else
                                {
                                    rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                                    rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                                    gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                                    gDmrs_re1 = _mm256_mul_ps(gDmrs_re1, alternateSigns);
                                    gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);
                                    gDmrs_im1 = _mm256_mul_ps(gDmrs_im1, alternateSigns);

                                    rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                                    rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                                    gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                                    gDmrs_re2 = _mm256_mul_ps(gDmrs_re2, alternateSigns);
                                    gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);
                                    gDmrs_im2 = _mm256_mul_ps(gDmrs_im2, alternateSigns);

                                    bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                    ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                    out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                    out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                    //out1 = ( Y1*conj(Xdmrs1) )/beta_dmrs;
                                    out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                    out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                    bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                    ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                    out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                    out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                    //out2 = ( Y2*conj(Xdmrs2) )/beta_dmrs;
                                    out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                    out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                    //Average estimates from two adjacent symbols.
                                    out_re = _mm256_sub_ps(out_re1, out_re2);
                                    out_im = _mm256_sub_ps(out_im1, out_im2);
                                    out_re = _mm256_mul_ps(out_re, div2Vec);
                                    out_im = _mm256_mul_ps(out_im, div2Vec);

                                    // Pad zeros between estimates for Linear Interpolation i.e., out = [h0 0 h1 0 h2...].
                                    temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                    temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                    temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                    temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                    temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                    temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                    temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                    temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                    _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                    _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                    _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                    _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                                }

                                recDmrs_rePtr1 += 8;
                                recDmrs_imPtr1 += 8;

                                genDmrs_rePtr1 += 8;
                                genDmrs_imPtr1 += 8;

                                recDmrs_rePtr2 += 8;
                                recDmrs_imPtr2 += 8;

                                genDmrs_rePtr2 += 8;
                                genDmrs_imPtr2 += 8;

                                chEst_rePtr += 16;
                                chEst_imPtr += 16;
                            }

                            //Frequency Interpolation
                            //Calculate estimates for null tones using Left and Right estimates.
                            //Inp : [ a 0 b 0 c ...]
                            //Out : [a (a+b)/2 b (b+c)/2 c ...]
                            chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                            chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];
                            AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                            AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                            for(idxSc = toneStartIdx; idxSc < toneEndIdx; idxSc += 8)
                            {
                                remSc = toneEndIdx - idxSc;
                                if(remSc<=8)
                                {
                                    for(i = 0; i < remSc; i += 2)
                                    {
                                        if(remSc - i == 2)
                                        {
                                            chEst_rePtr[i+1] = (chEst_rePtr[i]+chEst_rePtr[i-1])*0.5;
                                            chEst_imPtr[i+1] = (chEst_imPtr[i]+chEst_imPtr[i-1])*0.5;
                                        }
                                        else
                                        {
                                            chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                            chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                        }
                                    }
                                }
                                else
                                {
                                    AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                    AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                    temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                    temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                    temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                    temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                    temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                    temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                    temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                    temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                    temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                    temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                    temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                    temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                    AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                    AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                    _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                    _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                    AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                    AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                    chEst_rePtr += 8;
                                    chEst_imPtr += 8;
                                }
                            }
                        }
                    }
                }
            }
    

        }

    }

    if(nLyr == 1 || (nLyr == 2 && cdmTotal[0] == 1))
    {
        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
        {
            for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
            {
                for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
                {
                    startScIdx = pNrPuschInParams->rbStart[idxAlloc] * 6;
                    endScIdx   = startScIdx + (pNrPuschInParams->nPrb[idxAlloc] * 6);

                    toneStartIdx = pNrPuschInParams->rbStart[idxAlloc] * 12;
                    toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);

                    // port 0 processing
                    if(pNrPuschInParams->nPortIndex[idxLyr] == 0)
                    {          
                        //Received DMRS Signal 1                          
                        recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                        recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(MAX_NUM_SC/2) + startScIdx];
                        //Received DMRS Signal 2
                        recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][startScIdx];
                        recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(MAX_NUM_SC/2) + startScIdx];
                        //Generated DMRS Signal 1
                        genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                        genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                        //Generated DMRS Signal 2
                        genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                        genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];
                        //Channel Estimates after Linear Interpolation
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];

                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc < 8)
                            {
                                rDmrs_re1 = _mm256_setzero_ps();
                                rDmrs_im1 = _mm256_setzero_ps();
                                rDmrs_re2 = _mm256_setzero_ps();
                                rDmrs_im2 = _mm256_setzero_ps();

                                gDmrs_re1 = _mm256_setzero_ps();
                                gDmrs_im1 = _mm256_setzero_ps();
                                gDmrs_re2 = _mm256_setzero_ps();
                                gDmrs_im2 = _mm256_setzero_ps();

                                for(i = 0; i < remSc; i++)
                                {
                                    rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                                    rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                                    rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                                    rDmrs_im_f2[i] = recDmrs_imPtr2[i];

                                    gDmrs_re_f1[i] = genDmrs_rePtr1[i];
                                    gDmrs_im_f1[i] = genDmrs_imPtr1[i];
                                    gDmrs_re_f2[i] = genDmrs_rePtr2[i];
                                    gDmrs_im_f2[i] = genDmrs_imPtr2[i];
                                }

                                bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);

                                out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);

                                // out1 = ( Y1*conj(Xdmrs1) )/Beta_dmrs;
                                out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);

                                out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);

                                // out1 = ( Y2*conj(Xdmrs2) )/Beta_dmrs;
                                out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                out_im2 = _mm256_mul_ps(out_im2, dmrsScale);
                                
                                // Mean of two adjacent estimates.
                                out_re = _mm256_add_ps(out_re1, out_re2);
                                out_im = _mm256_add_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                for(i = 0; i < remSc; i++)
                                {
                                    chEst_rePtr[2*i] = out_re_f[i];
                                    chEst_rePtr[2*i+1] = 0;

                                    chEst_imPtr[2*i] = out_im_f[i];
                                    chEst_imPtr[2*i+1] = 0;
                                }
                            }
                            else
                            {
                                rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                                rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                                gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                                gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);

                                rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                                rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                                gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                                gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);

                                bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);

                                out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);

                                // out1 = ( Y1*conj(Xdmrs1) )/Beta_dmrs;
                                out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);

                                out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);

                                // out1 = ( Y2*conj(Xdmrs2) )/Beta_dmrs;
                                out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                out_im2 = _mm256_mul_ps(out_im2, dmrsScale);
                                
                                // Mean of two adjacent estimates.
                                out_re = _mm256_add_ps(out_re1, out_re2);
                                out_im = _mm256_add_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                            }
                            recDmrs_rePtr1 += 8;
                            recDmrs_imPtr1 += 8;

                            recDmrs_rePtr2 += 8;
                            recDmrs_imPtr2 += 8;

                            genDmrs_rePtr1 += 8;
                            genDmrs_imPtr1 += 8;

                            genDmrs_rePtr2 += 8;
                            genDmrs_imPtr2 += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }   //for idxSc

                        //Frequency Interpolation.
                        //Calculate estimates for null tones using Left and Right estimates.
                        //For dead end null tone, its two adjacent estimates will be used.
                        //Inp : [ a 0 b 0 c ...]
                        //Out : [a (a+b)/2 b (b+c)/2 c...]
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];
                        AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                        AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                        for(idxSc = toneStartIdx; idxSc < toneEndIdx; idxSc += 8)
                        {
                            remSc = toneEndIdx - idxSc;
                            if(remSc<=8)
                            {
                                for(i = 0; i < remSc; i += 2)
                                {
                                    if(remSc - i == 2)
                                    {
                                        chEst_rePtr[i+1] = (chEst_rePtr[i]+chEst_rePtr[i-1])*0.5;
                                        chEst_imPtr[i+1] = (chEst_imPtr[i]+chEst_imPtr[i-1])*0.5;
                                    }
                                    else
                                    {
                                        chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                        chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                    }
                                }
                            }
                            else
                            {
                                AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                chEst_rePtr += 8;
                                chEst_imPtr += 8;
                            }
                        }
                    }                  
                    // port 1 processing
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 1)
                    {
                        //Recived DMRS Signal 1
                        recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                        recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(MAX_NUM_SC/2) + startScIdx];
                        //Recived DMRS Signal 2
                        recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][startScIdx];
                        recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(MAX_NUM_SC/2) + startScIdx];
                        //DMRS Signal 1.
                        genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                        genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                        //DMRS Signal 2.
                        genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                        genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];
                        //Channel Estimates after Linear Interpolation.
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];

                        
                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if(remSc<8)
                            {
                                rDmrs_re1 = _mm256_setzero_ps();
                                rDmrs_im1 = _mm256_setzero_ps();
                                gDmrs_re1 = _mm256_setzero_ps();
                                gDmrs_im1 = _mm256_setzero_ps();

                                rDmrs_re2 = _mm256_setzero_ps();
                                rDmrs_im2 = _mm256_setzero_ps();
                                gDmrs_re2 = _mm256_setzero_ps();
                                gDmrs_im2 = _mm256_setzero_ps();

                                for(i = 0; i < remSc; i++)
                                {
                                    rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                                    rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                                    gDmrs_re_f1[i] = alternateSigns_f[i]*genDmrs_rePtr1[i];
                                    gDmrs_im_f1[i] = alternateSigns_f[i]*genDmrs_imPtr1[i];

                                    rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                                    rDmrs_im_f2[i] = recDmrs_imPtr2[i];
                                    gDmrs_re_f2[i] = alternateSigns_f[i]*genDmrs_rePtr2[i];
                                    gDmrs_im_f2[i] = alternateSigns_f[i]*genDmrs_imPtr2[i];
                                }

                                bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);

                                out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);

                                // out1 = ( Y1*conj(Xdmrs1) )/beta_dmrs;
                                out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);

                                out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);

                                // out2 = ( Y2*conj(Xdmrs2) )/beta_dmrs;
                                out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                out_re = _mm256_add_ps(out_re1, out_re2);
                                out_im = _mm256_add_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                // Pad zeros between estimates for Linear Interpolation i.e., out = [h0 0 h1 0 h2...].
                                for(i = 0; i < remSc; i++)
                                {
                                    chEst_rePtr[2*i] = out_re[i];
                                    chEst_rePtr[2*i+1] = 0;

                                    chEst_imPtr[2*i] = out_im[i];
                                    chEst_imPtr[2*i+1] = 0;
                                }
                            }
                            else
                            {
                                rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                                rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                                gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                                gDmrs_re1 = _mm256_mul_ps(gDmrs_re1, alternateSigns);
                                gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);
                                gDmrs_im1 = _mm256_mul_ps(gDmrs_im1, alternateSigns);

                                rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                                rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                                gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                                gDmrs_re2 = _mm256_mul_ps(gDmrs_re2, alternateSigns);
                                gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);
                                gDmrs_im2 = _mm256_mul_ps(gDmrs_im2, alternateSigns);

                                bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);

                                out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);

                                // out1 = ( Y1*conj(Xdmrs1) )/beta_dmrs;
                                out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);

                                out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);

                                // out2 = ( Y2*conj(Xdmrs2) )/beta_dmrs;
                                out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                out_re = _mm256_add_ps(out_re1, out_re2);
                                out_im = _mm256_add_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);
                                
                                // Pad zeros between estimates for Linear Interpolation i.e., out = [h0 0 h1 0 h2...].
                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                            }
                            
                            
                            recDmrs_rePtr1 += 8;
                            recDmrs_imPtr1 += 8;

                            genDmrs_rePtr1 += 8;
                            genDmrs_imPtr1 += 8;

                            recDmrs_rePtr2 += 8;
                            recDmrs_imPtr2 += 8;

                            genDmrs_rePtr2 += 8;
                            genDmrs_imPtr2 += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }   //for idxSc

                        //Frequency Interpolation
                        //Calculate estimates for null tones using Left and Right estimates.
                        //Inp : [ a 0 b 0 c ...]
                        //Out : [a (a+b)/2 b (b+c)/2 c ...]
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];
                        AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                        AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                        for(idxSc = toneStartIdx; idxSc < toneEndIdx; idxSc += 8)
                        {
                            remSc = toneEndIdx - idxSc;
                            if(remSc<=8)
                            {
                                for(i = 0; i < remSc; i += 2)
                                {
                                    if(remSc - i == 2)
                                    {
                                        chEst_rePtr[i+1] = (chEst_rePtr[i]+chEst_rePtr[i-1])*0.5;
                                        chEst_imPtr[i+1] = (chEst_imPtr[i]+chEst_imPtr[i-1])*0.5;
                                    }
                                    else
                                    {
                                        chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                        chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                    }
                                }
                            }
                            else
                            {
                                AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                chEst_rePtr += 8;
                                chEst_imPtr += 8;
                            }
                        }
                    }
                    // port 2 processing
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 2)
                    {
                        //Received DMRS Signal 1.
                        recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC + startScIdx];
                        recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                        //DMRS Signal 1.
                        genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                        genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                        //Received DMRS Signal 2.
                        recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][MAX_NUM_SC + startScIdx];
                        recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                        //DMRS Signal 2.
                        genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                        genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];
                        //Channel Estimates after Linear Interpolation.
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc < 8)
                            {
                                rDmrs_re1 = _mm256_setzero_ps();
                                rDmrs_im1 = _mm256_setzero_ps();
                                gDmrs_re1 = _mm256_setzero_ps();
                                gDmrs_im1 = _mm256_setzero_ps();

                                rDmrs_re2 = _mm256_setzero_ps();
                                rDmrs_im2 = _mm256_setzero_ps();
                                gDmrs_re2 = _mm256_setzero_ps();
                                gDmrs_im2 = _mm256_setzero_ps();

                                for(i = 0; i < remSc; i++)
                                {
                                    rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                                    rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                                    gDmrs_re_f1[i] = genDmrs_rePtr1[i];
                                    gDmrs_im_f1[i] = genDmrs_imPtr1[i];

                                    rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                                    rDmrs_im_f2[i] = recDmrs_imPtr2[i];
                                    gDmrs_re_f2[i] = genDmrs_rePtr2[i];
                                    gDmrs_im_f2[i] = genDmrs_imPtr2[i];
                                }

                                bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                // out1 = ( Y1*conj(Xdmrs1) )/beta_DMRS;
                                out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                out_im1 = _mm256_mul_ps(out_im1, dmrsScale);
                                bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                // out2 = ( Y2*conj(Xdmrs2) )/beta_DMRS;
                                out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                out_im2 = _mm256_mul_ps(out_im2, dmrsScale);
                                out_re = _mm256_add_ps(out_re1, out_re2);
                                out_im = _mm256_add_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                for(i = 0; i < remSc; i++)
                                {
                                    chEst_rePtr[2*i] = out_re_f[i];
                                    chEst_rePtr[2*i+1] = 0;

                                    chEst_imPtr[2*i] = out_im_f[i];
                                    chEst_imPtr[2*i+1] = 0;
                                }
                            }
                            else
                            {
                                rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                                rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                                gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                                gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);

                                rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                                rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                                gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                                gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);

                                bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                // out1 = ( Y1*conj(Xdmrs1) )/beta_DMRS;
                                out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                out_im1 = _mm256_mul_ps(out_im1, dmrsScale);
                                bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                // out2 = ( Y2*conj(Xdmrs2) )/beta_DMRS;
                                out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                out_im2 = _mm256_mul_ps(out_im2, dmrsScale);
                                out_re = _mm256_add_ps(out_re1, out_re2);
                                out_im = _mm256_add_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                            }
                            recDmrs_rePtr1 += 8;
                            recDmrs_imPtr1 += 8;

                            genDmrs_rePtr1 += 8;
                            genDmrs_imPtr1 += 8;

                            recDmrs_rePtr2 += 8;
                            recDmrs_imPtr2 += 8;

                            genDmrs_rePtr2 += 8;
                            genDmrs_imPtr2 += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }

                        //Frequency Interpolation.
                        //Calculate estimates for null tones using Left and Right estimates.
                        //Inp : [0 a 0 b 0 c ...]
                        //Out : [(a+(a+b)/2)/2 a (a+b)/2 b (b+c)/2 c ...]
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];
                        AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                        AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                        for(idxSc = toneStartIdx+1; idxSc < toneEndIdx; idxSc += 8)
                        {
                            remSc = toneEndIdx - idxSc;
                            if(remSc<8)
                            {
                                for(i = 0; i < remSc; i += 2)
                                {
                                    chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                    chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                }

                                pNrPuschOutParams->chEst_re[idxAnt][idxLyr][0] = \
                                (pNrPuschOutParams->chEst_re[idxAnt][idxLyr][1]+\
                                pNrPuschOutParams->chEst_re[idxAnt][idxLyr][2])*0.5;

                                pNrPuschOutParams->chEst_im[idxAnt][idxLyr][0] = \
                                (pNrPuschOutParams->chEst_im[idxAnt][idxLyr][1]+\
                                pNrPuschOutParams->chEst_im[idxAnt][idxLyr][2])*0.5;
                            }
                            else
                            {
                                AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);


                                AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                chEst_rePtr += 8;
                                chEst_imPtr += 8;
                            }
                        }
                    }
                    //port 3 processing
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 3)
                    {
                        //Received DMRS Signal 1.
                        recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC + startScIdx];
                        recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                        //DMRS Signal 1.
                        genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                        genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                        //Received DMRS Signal 2.
                        recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][MAX_NUM_SC + startScIdx];
                        recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                        //DMRS Signal 2.
                        genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                        genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];
                        //Channel Estimates after Linear Interpolation.
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            
                            if( (endScIdx - idxSc) < 8)
                            {
                                remSc = endScIdx - idxSc;

                                rDmrs_re1 = _mm256_setzero_ps();
                                rDmrs_im1 = _mm256_setzero_ps();
                                gDmrs_re1 = _mm256_setzero_ps();
                                gDmrs_im1 = _mm256_setzero_ps();

                                rDmrs_re2 = _mm256_setzero_ps();
                                rDmrs_im2 = _mm256_setzero_ps();
                                gDmrs_re2 = _mm256_setzero_ps();
                                gDmrs_im2 = _mm256_setzero_ps();

                                for(i = 0; i < remSc; i++)
                                {
                                    rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                                    rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                                    gDmrs_re_f1[i] = alternateSigns_f[i]*genDmrs_rePtr1[i];
                                    gDmrs_im_f1[i] = alternateSigns_f[i]*genDmrs_imPtr1[i];

                                    rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                                    rDmrs_im_f2[i] = recDmrs_imPtr2[i];
                                    gDmrs_re_f2[i] = alternateSigns_f[i]*genDmrs_rePtr2[i];
                                    gDmrs_im_f2[i] = alternateSigns_f[i]*genDmrs_imPtr2[i];
                                }

                                bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                // out1 = ( Y1*conj(Xdmrs1) )/beta_DMRS;
                                out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                // out2 = ( Y2*conj(Xdmrs2) )/beta_DMRS;
                                out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                out_im2 = _mm256_mul_ps(out_im2, dmrsScale);
                                
                                //Average estimates from two adjacent symbols.
                                out_re = _mm256_add_ps(out_re1, out_re2);
                                out_im = _mm256_add_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                // Pad zeros between estimates for Linear Interpolation i.e., out = [0 h0 0 h1 0 h2...].
                                for(i = 0; i < remSc; i++)
                                {
                                    chEst_rePtr[2*i] = out_re_f[i];
                                    chEst_rePtr[2*i+1] = 0;

                                    chEst_imPtr[2*i] = out_im_f[i];
                                    chEst_imPtr[2*i+1] = 0;
                                }
                            }
                            else
                            {
                                rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                                rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                                //Reconstruct DMRS signal 1 for port 3.
                                gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                                gDmrs_re1 = _mm256_mul_ps(gDmrs_re1, alternateSigns);
                                gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);
                                gDmrs_im1 = _mm256_mul_ps(gDmrs_im1, alternateSigns);

                                rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                                rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                                //Reconstruct DMRS signal 2 for port 3.
                                gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                                gDmrs_re2 = _mm256_mul_ps(gDmrs_re2, alternateSigns);
                                gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);
                                gDmrs_im2 = _mm256_mul_ps(gDmrs_im2, alternateSigns);

                                bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                // out1 = ( Y1*conj(Xdmrs1) )/beta_DMRS;
                                out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                // out2 = ( Y2*conj(Xdmrs2) )/beta_DMRS;
                                out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                out_im2 = _mm256_mul_ps(out_im2, dmrsScale);
                                
                                //Average estimates from two adjacent symbols.
                                out_re = _mm256_add_ps(out_re1, out_re2);
                                out_im = _mm256_add_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                // Pad zeros between estimates for Linear Interpolation i.e., out = [0 h0 0 h1 0 h2...].
                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                            }

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                            
                            recDmrs_rePtr1 += 8;
                            recDmrs_imPtr1 += 8;

                            genDmrs_rePtr1 += 8;
                            genDmrs_imPtr1 += 8;

                            recDmrs_rePtr2 += 8;
                            recDmrs_imPtr2 += 8;

                            genDmrs_rePtr2 += 8;
                            genDmrs_imPtr2 += 8;
                        }
                        
                        //Frequency Interpolation
                        //Calculate estimates for null tones using Left and Right estimates.
                        //For dead end null tone, its two adjacent estimates will be used.
                        //Inp : [0 a 0 b 0 c ...]
                        //Out : [(a+(a+b)/2)/2 a (a+b)/2 b (b+c)/2 c ...]
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];
                        AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                        AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                        for(idxSc = toneStartIdx+1; idxSc < toneEndIdx; idxSc += 8)
                        {
                            remSc = toneEndIdx - idxSc;
                            if(remSc<8)
                            {
                                for(i = 0; i < remSc; i += 2)
                                {
                                    chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                    chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                }

                                pNrPuschOutParams->chEst_re[idxAnt][idxLyr][0] = \
                                (pNrPuschOutParams->chEst_re[idxAnt][idxLyr][1]+\
                                pNrPuschOutParams->chEst_re[idxAnt][idxLyr][2])*0.5;

                                pNrPuschOutParams->chEst_im[idxAnt][idxLyr][0] = \
                                (pNrPuschOutParams->chEst_im[idxAnt][idxLyr][1]+\
                                pNrPuschOutParams->chEst_im[idxAnt][idxLyr][2])*0.5;
                            }
                            else
                            {
                                AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);


                                AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                chEst_rePtr += 8;
                                chEst_imPtr += 8;
                            }
                        }
                    }
                    //port 4 processing
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 4)
                    {          
                        //Received DMRS Signal 1.                         
                        recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                        recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(MAX_NUM_SC/2) + startScIdx];
                        //Received DMRS Signal 2.
                        recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][startScIdx];
                        recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(MAX_NUM_SC/2) + startScIdx];
                        //DMRS Signal 1.
                        genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                        genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                        //DMRS Signal 2.
                        genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                        genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];
                        //Channel Estimates after Linear Interpolation.
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];

                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc < 8)
                            {
                                rDmrs_re1 = _mm256_setzero_ps();
                                rDmrs_im1 = _mm256_setzero_ps();
                                rDmrs_re2 = _mm256_setzero_ps();
                                rDmrs_im2 = _mm256_setzero_ps();

                                gDmrs_re1 = _mm256_setzero_ps();
                                gDmrs_im1 = _mm256_setzero_ps();
                                gDmrs_re2 = _mm256_setzero_ps();
                                gDmrs_im2 = _mm256_setzero_ps();

                                for(i = 0; i < remSc; i++)
                                {
                                    rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                                    rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                                    rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                                    rDmrs_im_f2[i] = recDmrs_imPtr2[i];

                                    gDmrs_re_f1[i] = genDmrs_rePtr1[i];
                                    gDmrs_im_f1[i] = genDmrs_imPtr1[i];
                                    gDmrs_re_f2[i] = genDmrs_rePtr2[i];
                                    gDmrs_im_f2[i] = genDmrs_imPtr2[i];
                                }

                                bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                //out1 = ( Y1*conj(Xdmrs1) )/Beta_dmrs;
                                out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                //out2 = ( Y2*conj(Xdmrs2) )/Beta_dmrs;
                                out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                //Average estimates from two adjacent symbols.
                                out_re = _mm256_sub_ps(out_re1, out_re2);
                                out_im = _mm256_sub_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);
                                // Pad zeros between estimates for Linear Interpolation i.e., out = [h0 0 h1 0 h2 ...];
                                for(i = 0; i < remSc; i++)
                                {
                                    chEst_rePtr[2*i] = out_re_f[i];
                                    chEst_rePtr[2*i+1] = 0;

                                    chEst_imPtr[2*i] = out_im_f[i];
                                    chEst_imPtr[2*i+1] = 0;
                                }
                            }
                            else
                            {
                                rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                                rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                                gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                                gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);

                                rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                                rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                                gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                                gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);

                                bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                //out1 = ( Y1*conj(Xdmrs1) )/Beta_dmrs;
                                out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                //out2 = ( Y2*conj(Xdmrs2) )/Beta_dmrs;
                                out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                //Average estimates from two adjacent symbols.
                                out_re = _mm256_sub_ps(out_re1, out_re2);
                                out_im = _mm256_sub_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                // Pad zeros between estimates for Linear Interpolation i.e., out = [h0 0 h1 0 h2 ...];
                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                            }
                            
                            recDmrs_rePtr1 += 8;
                            recDmrs_imPtr1 += 8;

                            recDmrs_rePtr2 += 8;
                            recDmrs_imPtr2 += 8;

                            genDmrs_rePtr1 += 8;
                            genDmrs_imPtr1 += 8;

                            genDmrs_rePtr2 += 8;
                            genDmrs_imPtr2 += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }

                        //Frequency Interpolation 
                        //Calculate estimates for null tones using Left and Right estimates.
                        //Inp : [ a 0 b 0 c ...]
                        //Out : [a (a+b)/2 b (b+c)/2 c...]
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];
                        AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                        AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                        for(idxSc = toneStartIdx; idxSc < toneEndIdx; idxSc += 8)
                        {
                            remSc = toneEndIdx - idxSc;
                            if(remSc<=8)
                            {
                                for(i = 0; i < remSc; i += 2)
                                {
                                    if(remSc - i == 2)
                                    {
                                        chEst_rePtr[i+1] = (chEst_rePtr[i]+chEst_rePtr[i-1])*0.5;
                                        chEst_imPtr[i+1] = (chEst_imPtr[i]+chEst_imPtr[i-1])*0.5;
                                    }
                                    else
                                    {
                                        chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                        chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                    }
                                }
                            }
                            else
                            {
                                AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                chEst_rePtr += 8;
                                chEst_imPtr += 8;
                            }
                        }
                    }
                    //port 5 processing
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 5)
                    {
                        //Received DMRS signal 1.
                        recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][startScIdx];
                        recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(MAX_NUM_SC/2) + startScIdx];
                        //Received DMRS signal 2.
                        recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][startScIdx];
                        recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(MAX_NUM_SC/2) + startScIdx];
                        //DMRS signal 1.
                        genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                        genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                        //DMRS signal 2.
                        genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                        genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];
                        //Channel Estimates after Linear Interpolation.
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];
                        
                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            if( (endScIdx - idxSc)<8)
                            {
                                remSc = endScIdx - idxSc;
                                rDmrs_re1 = _mm256_setzero_ps();
                                rDmrs_im1 = _mm256_setzero_ps();
                                gDmrs_re1 = _mm256_setzero_ps();
                                gDmrs_im1 = _mm256_setzero_ps();

                                rDmrs_re2 = _mm256_setzero_ps();
                                rDmrs_im2 = _mm256_setzero_ps();
                                gDmrs_re2 = _mm256_setzero_ps();
                                gDmrs_im2 = _mm256_setzero_ps();

                                for(i = 0; i < remSc; i++)
                                {
                                    rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                                    rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                                    gDmrs_re_f1[i] = alternateSigns_f[i]*genDmrs_rePtr1[i];
                                    gDmrs_im_f1[i] = alternateSigns_f[i]*genDmrs_imPtr1[i];

                                    rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                                    rDmrs_im_f2[i] = recDmrs_imPtr2[i];
                                    gDmrs_re_f2[i] = alternateSigns_f[i]*genDmrs_rePtr2[i];
                                    gDmrs_im_f2[i] = alternateSigns_f[i]*genDmrs_imPtr2[i];
                                }

                                bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                //out1 = ( Y1*conj(Xdmrs1) )/beta_dmrs;
                                out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                //out2 = ( Y2*conj(Xdmrs2) )/beta_dmrs;
                                out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                //Average estimates from two adjacent symbols.
                                out_re = _mm256_sub_ps(out_re1, out_re2);
                                out_im = _mm256_sub_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                // Pad zeros between estimates for Linear Interpolation i.e., out = [h0 0 h1 0 h2...].
                                for(i = 0; i < remSc; i++)
                                {
                                    chEst_rePtr[2*i] = out_re_f[i];
                                    chEst_rePtr[2*i+1] = 0;
                                    chEst_imPtr[2*i] = out_im_f[i];
                                    chEst_imPtr[2*i+1] = 0;
                                }
                            }
                            else
                            {
                                rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                                rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                                gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                                gDmrs_re1 = _mm256_mul_ps(gDmrs_re1, alternateSigns);
                                gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);
                                gDmrs_im1 = _mm256_mul_ps(gDmrs_im1, alternateSigns);

                                rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                                rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                                gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                                gDmrs_re2 = _mm256_mul_ps(gDmrs_re2, alternateSigns);
                                gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);
                                gDmrs_im2 = _mm256_mul_ps(gDmrs_im2, alternateSigns);

                                bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                //out1 = ( Y1*conj(Xdmrs1) )/beta_dmrs;
                                out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                //out2 = ( Y2*conj(Xdmrs2) )/beta_dmrs;
                                out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                //Average estimates from two adjacent symbols.
                                out_re = _mm256_sub_ps(out_re1, out_re2);
                                out_im = _mm256_sub_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                // Pad zeros between estimates for Linear Interpolation i.e., out = [h0 0 h1 0 h2...].
                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                            }

                            recDmrs_rePtr1 += 8;
                            recDmrs_imPtr1 += 8;

                            genDmrs_rePtr1 += 8;
                            genDmrs_imPtr1 += 8;

                            recDmrs_rePtr2 += 8;
                            recDmrs_imPtr2 += 8;

                            genDmrs_rePtr2 += 8;
                            genDmrs_imPtr2 += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }

                        //Frequency Interpolation
                        //Calculate estimates for null tones using Left and Right estimates.
                        //Inp : [ a 0 b 0 c ...]
                        //Out : [a (a+b)/2 b (b+c)/2 c ...]
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx];
                        AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                        AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                        for(idxSc = toneStartIdx; idxSc < toneEndIdx; idxSc += 8)
                        {
                            remSc = toneEndIdx - idxSc;
                            if(remSc<=8)
                            {
                                for(i = 0; i < remSc; i += 2)
                                {
                                    if(remSc - i == 2)
                                    {
                                        chEst_rePtr[i+1] = (chEst_rePtr[i]+chEst_rePtr[i-1])*0.5;
                                        chEst_imPtr[i+1] = (chEst_imPtr[i]+chEst_imPtr[i-1])*0.5;
                                    }
                                    else
                                    {
                                        chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                        chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                    }
                                }
                            }
                            else
                            {
                                AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                chEst_rePtr += 8;
                                chEst_imPtr += 8;
                            }
                        }
                    }
                    //port 6 processing
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 6)
                    {
                        //Received DMRS Signal 1.
                        recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC + startScIdx];
                        recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                        //DMRS Signal 1.
                        genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                        genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                        //Received DMRS Signal 2.
                        recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][MAX_NUM_SC + startScIdx];
                        recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                        //DMRS Signal 2.
                        genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                        genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];
                        //Channel Estimates after Linear Interpolation.
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = (endScIdx) - idxSc;
                            if( remSc < 8)
                            {
                                rDmrs_re1 = _mm256_setzero_ps();
                                rDmrs_im1 = _mm256_setzero_ps();
                                gDmrs_re1 = _mm256_setzero_ps();
                                gDmrs_im1 = _mm256_setzero_ps();

                                rDmrs_re2 = _mm256_setzero_ps();
                                rDmrs_im2 = _mm256_setzero_ps();
                                gDmrs_re2 = _mm256_setzero_ps();
                                gDmrs_im2 = _mm256_setzero_ps();

                                for(i = 0; i < remSc; i++)
                                {
                                    rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                                    rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                                    gDmrs_re_f1[i] = genDmrs_rePtr1[i];
                                    gDmrs_im_f1[i] = genDmrs_imPtr1[i];

                                    rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                                    rDmrs_im_f2[i] = recDmrs_imPtr2[i];
                                    gDmrs_re_f2[i] = genDmrs_rePtr2[i];
                                    gDmrs_im_f2[i] = genDmrs_imPtr2[i];
                                }

                                bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                //out1 = ( Y1*conj(Xdmrs1) )/beta_DMRS;
                                out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                //out2 = ( Y2*conj(Xdmrs2) )/beta_DMRS;
                                out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                //Average estimates from two adjacent symbols.
                                out_re = _mm256_sub_ps(out_re1, out_re2);
                                out_im = _mm256_sub_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                //Pad zeros between estimates for Linear Interpolation i.e., out = [0 h0 0 h1 0 h2 ...].
                                for(i = 0; i < remSc; i++)
                                {
                                    chEst_rePtr[2*i] = out_re_f[i];
                                    chEst_rePtr[2*i+1] = 0;

                                    chEst_imPtr[2*i] = out_im_f[i];
                                    chEst_imPtr[2*i+1] = 0;
                                }
                            }
                            else
                            {
                                rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                                rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                                gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                                gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);

                                rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                                rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                                gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                                gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);

                                bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                //out1 = ( Y1*conj(Xdmrs1) )/beta_DMRS;
                                out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                //out2 = ( Y2*conj(Xdmrs2) )/beta_DMRS;
                                out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                //Average estimates from two adjacent symbols.
                                out_re = _mm256_sub_ps(out_re1, out_re2);
                                out_im = _mm256_sub_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                //Pad zeros between estimates for Linear Interpolation i.e., out = [0 h0 0 h1 0 h2 ...].
                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                            }

                            recDmrs_rePtr1 += 8;
                            recDmrs_imPtr1 += 8;

                            genDmrs_rePtr1 += 8;
                            genDmrs_imPtr1 += 8;

                            recDmrs_rePtr2 += 8;
                            recDmrs_imPtr2 += 8;

                            genDmrs_rePtr2 += 8;
                            genDmrs_imPtr2 += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }

                        //Frequency Interpolation
                        //Calculate estimates for null tones using Left and Right estimates.
                        //Inp : [0 a 0 b 0 c ...]
                        //Out : [(a+(a+b)/2)/2 a (a+b)/2 b (b+c)/2 c ...]
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                        AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                        AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                        for(idxSc = toneStartIdx+1; idxSc < toneEndIdx; idxSc += 8)
                        {
                            remSc = toneEndIdx - idxSc;
                            if(remSc<8)
                            {
                                for(i = 0; i < remSc; i += 2)
                                {
                                    chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                    chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                }

                                pNrPuschOutParams->chEst_re[idxAnt][idxLyr][0] = \
                                (pNrPuschOutParams->chEst_re[idxAnt][idxLyr][1]+\
                                pNrPuschOutParams->chEst_re[idxAnt][idxLyr][2])*0.5;

                                pNrPuschOutParams->chEst_im[idxAnt][idxLyr][0] = \
                                (pNrPuschOutParams->chEst_im[idxAnt][idxLyr][1]+\
                                pNrPuschOutParams->chEst_im[idxAnt][idxLyr][2])*0.5;
                            }
                            else
                            {
                                AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                chEst_rePtr += 8;
                                chEst_imPtr += 8;
                            }
                        }
                    }
                    //port 7 processing
                    else if(pNrPuschInParams->nPortIndex[idxLyr] == 7)
                    {
                        //Received DMRS Signal 1.
                        recDmrs_rePtr1 = &rxFdSamples[idxAnt][dmrsSym][MAX_NUM_SC + startScIdx];
                        recDmrs_imPtr1 = &rxFdSamples[idxAnt][dmrsSym][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                        //DMRS Signal 1.
                        genDmrs_rePtr1 = &pNrPuschOutParams->genDMRS[0][startScIdx];
                        genDmrs_imPtr1 = &pNrPuschOutParams->genDMRS[0][(MAX_NUM_SC/2) + startScIdx];
                        //Received DMRS Signal 2.
                        recDmrs_rePtr2 = &rxFdSamples[idxAnt][dmrsSym+1][MAX_NUM_SC + startScIdx];
                        recDmrs_imPtr2 = &rxFdSamples[idxAnt][dmrsSym+1][(int)(MAX_NUM_SC * 1.5) + startScIdx];
                        //DMRS Signal 2.
                        genDmrs_rePtr2 = &pNrPuschOutParams->genDMRS[1][startScIdx];
                        genDmrs_imPtr2 = &pNrPuschOutParams->genDMRS[1][(MAX_NUM_SC/2) + startScIdx];
                        //Channel Estimates after Linear Interpolation.
                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                        for(idxSc = startScIdx; idxSc < endScIdx; idxSc += 8)
                        {
                            remSc = endScIdx - idxSc;
                            if( remSc < 8)
                            {
                                rDmrs_re1 = _mm256_setzero_ps();
                                rDmrs_im1 = _mm256_setzero_ps();
                                gDmrs_re1 = _mm256_setzero_ps();
                                gDmrs_im1 = _mm256_setzero_ps();

                                rDmrs_re2 = _mm256_setzero_ps();
                                rDmrs_im2 = _mm256_setzero_ps();
                                gDmrs_re2 = _mm256_setzero_ps();
                                gDmrs_im2 = _mm256_setzero_ps();

                                for(i = 0; i < remSc; i++)
                                {
                                    rDmrs_re_f1[i] = recDmrs_rePtr1[i];
                                    rDmrs_im_f1[i] = recDmrs_imPtr1[i];
                                    gDmrs_re_f1[i] = alternateSigns_f[i]*genDmrs_rePtr1[i];
                                    gDmrs_im_f1[i] = alternateSigns_f[i]*genDmrs_imPtr1[i];

                                    rDmrs_re_f2[i] = recDmrs_rePtr2[i];
                                    rDmrs_im_f2[i] = recDmrs_imPtr2[i];
                                    gDmrs_re_f2[i] = alternateSigns_f[i]*genDmrs_rePtr2[i];
                                    gDmrs_im_f2[i] = alternateSigns_f[i]*genDmrs_imPtr2[i];
                                }

                                bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                // out = ( Y*conj(Xdmrs) )/beta_DMRS;
                                out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                // out = ( Y*conj(Xdmrs) )/beta_DMRS;
                                out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                //Average estimates from two adjacent symbols.
                                out_re = _mm256_sub_ps(out_re1, out_re2);
                                out_im = _mm256_sub_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                // Pad zeros between estimates for Linear Interpolation i.e., out = [0 h0 0 h1 0 h2 ...].
                                for(i = 0; i < remSc; i++)
                                {
                                    chEst_rePtr[2*i] = out_re_f[i];
                                    chEst_rePtr[2*i+1] = 0;

                                    chEst_imPtr[2*i] = out_im_f[i];
                                    chEst_imPtr[2*i+1] = 0;
                                }
                            }
                            else
                            {
                                rDmrs_re1 = _mm256_loadu_ps(recDmrs_rePtr1);
                                rDmrs_im1 = _mm256_loadu_ps(recDmrs_imPtr1);
                                gDmrs_re1 = _mm256_loadu_ps(genDmrs_rePtr1);
                                gDmrs_re1 = _mm256_mul_ps(gDmrs_re1, alternateSigns);
                                gDmrs_im1 = _mm256_loadu_ps(genDmrs_imPtr1);
                                gDmrs_im1 = _mm256_mul_ps(gDmrs_im1, alternateSigns);

                                rDmrs_re2 = _mm256_loadu_ps(recDmrs_rePtr2);
                                rDmrs_im2 = _mm256_loadu_ps(recDmrs_imPtr2);
                                gDmrs_re2 = _mm256_loadu_ps(genDmrs_rePtr2);
                                gDmrs_re2 = _mm256_mul_ps(gDmrs_re2, alternateSigns);
                                gDmrs_im2 = _mm256_loadu_ps(genDmrs_imPtr2);
                                gDmrs_im2 = _mm256_mul_ps(gDmrs_im2, alternateSigns);

                                bd_re = _mm256_mul_ps(rDmrs_im1, gDmrs_im1);
                                ad_im = _mm256_mul_ps(rDmrs_re1, gDmrs_im1);
                                out_re1 = _mm256_fmadd_ps(rDmrs_re1, gDmrs_re1, bd_re);
                                out_im1 = _mm256_fmsub_ps(rDmrs_im1, gDmrs_re1, ad_im);
                                // out = ( Y*conj(Xdmrs) )/beta_DMRS;
                                out_re1 = _mm256_mul_ps(out_re1, dmrsScale);
                                out_im1 = _mm256_mul_ps(out_im1, dmrsScale);

                                bd_re = _mm256_mul_ps(rDmrs_im2, gDmrs_im2);
                                ad_im = _mm256_mul_ps(rDmrs_re2, gDmrs_im2);
                                out_re2 = _mm256_fmadd_ps(rDmrs_re2, gDmrs_re2, bd_re);
                                out_im2 = _mm256_fmsub_ps(rDmrs_im2, gDmrs_re2, ad_im);
                                // out = ( Y*conj(Xdmrs) )/beta_DMRS;
                                out_re2 = _mm256_mul_ps(out_re2, dmrsScale);
                                out_im2 = _mm256_mul_ps(out_im2, dmrsScale);

                                //Average estimates from two adjacent symbols.
                                out_re = _mm256_sub_ps(out_re1, out_re2);
                                out_im = _mm256_sub_ps(out_im1, out_im2);
                                out_re = _mm256_mul_ps(out_re, div2Vec);
                                out_im = _mm256_mul_ps(out_im, div2Vec);

                                // Pad zeros between estimates for Linear Interpolation i.e., out = [0 h0 0 h1 0 h2 ...].
                                temp_re1 = _mm256_permute_ps(out_re, 0b11010100);
                                temp_re1 = _mm256_blend_ps(temp_re1, zeroVec, 0b10101010);
                                temp_re2 = _mm256_permute_ps(out_re, 0b11110110);
                                temp_re2 = _mm256_blend_ps(temp_re2, zeroVec, 0b10101010);
                                temp_im1 = _mm256_permute_ps(out_im, 0b11010100);
                                temp_im1 = _mm256_blend_ps(temp_im1, zeroVec, 0b10101010);
                                temp_im2 = _mm256_permute_ps(out_im, 0b11110110);
                                temp_im2 = _mm256_blend_ps(temp_im2, zeroVec, 0b10101010);
                                _mm256_storeu_ps(chEst_rePtr, _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_rePtr+8), _mm256_permute2f128_ps(temp_re1, temp_re2, 0b00110001) );
                                _mm256_storeu_ps( chEst_imPtr, _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00100000) );
                                _mm256_storeu_ps( (chEst_imPtr+8), _mm256_permute2f128_ps(temp_im1, temp_im2, 0b00110001) );
                            }

                            recDmrs_rePtr1 += 8;
                            recDmrs_imPtr1 += 8;

                            genDmrs_rePtr1 += 8;
                            genDmrs_imPtr1 += 8;

                            recDmrs_rePtr2 += 8;
                            recDmrs_imPtr2 += 8;

                            genDmrs_rePtr2 += 8;
                            genDmrs_imPtr2 += 8;

                            chEst_rePtr += 16;
                            chEst_imPtr += 16;
                        }

                        //Frequency Interpolation
                        //Calculate estimates for null tones using Left and Right estimates.
                        //Inp : [0 a 0 b 0 c ...]
                        //Out : [(a+(a+b)/2)/2 a (a+b)/2 b (b+c)/2 c ...]

                        chEst_rePtr = &pNrPuschOutParams->chEst_re[idxAnt][idxLyr][toneStartIdx+1];
                        chEst_imPtr = &pNrPuschOutParams->chEst_im[idxAnt][idxLyr][toneStartIdx+1];

                        AVX_chEst_rePtr = (__m256*)chEst_rePtr;
                        AVX_chEst_imPtr = (__m256*)chEst_imPtr;

                        for(idxSc = toneStartIdx+1; idxSc < toneEndIdx; idxSc += 8)
                        {
                            remSc = toneEndIdx - idxSc;
                            if( remSc <8)
                            {
                                for(i = 0; i < remSc; i += 2)
                                {
                                    chEst_rePtr[i+1] = ( chEst_rePtr[i] + chEst_rePtr[i+2] )*0.5;
                                    chEst_imPtr[i+1] = ( chEst_imPtr[i] + chEst_imPtr[i+2] )*0.5;
                                }

                                pNrPuschOutParams->chEst_re[idxAnt][idxLyr][0] = \
                                (pNrPuschOutParams->chEst_re[idxAnt][idxLyr][1]+\
                                pNrPuschOutParams->chEst_re[idxAnt][idxLyr][2])*0.5;

                                pNrPuschOutParams->chEst_im[idxAnt][idxLyr][0] = \
                                (pNrPuschOutParams->chEst_im[idxAnt][idxLyr][1]+\
                                pNrPuschOutParams->chEst_im[idxAnt][idxLyr][2])*0.5;
                            }
                            else
                            {
                                AVX_chEst_re = _mm256_loadu_ps(chEst_rePtr);
                                AVX_chEst_im = _mm256_loadu_ps(chEst_imPtr);

                                temp_re1 = _mm256_permutevar8x32_ps(AVX_chEst_re, const1);
                                temp_im1 = _mm256_permutevar8x32_ps(AVX_chEst_im, const1 );
                                temp_re2 = _mm256_permute_ps( AVX_chEst_re, 0b10010001);
                                temp_im2 = _mm256_permute_ps( AVX_chEst_im, 0b10010001);
                                temp_re1 = _mm256_blend_ps(temp_re1, *(AVX_chEst_rePtr+1), 0b00000001);
                                temp_re1 = _mm256_permutevar8x32_ps(temp_re1,  const2);
                                temp_im1 = _mm256_blend_ps(temp_im1, *(AVX_chEst_imPtr+1), 0b00000001);
                                temp_im1 = _mm256_permutevar8x32_ps(temp_im1, const2 );
                                temp_re1 = _mm256_add_ps(temp_re1, temp_re2);
                                temp_re1 = _mm256_mul_ps(temp_re1, _mm256_set1_ps(0.5));
                                temp_im1 = _mm256_add_ps(temp_im1, temp_im2);
                                temp_im1 = _mm256_mul_ps(temp_im1, _mm256_set1_ps(0.5));

                                AVX_chEst_re = _mm256_add_ps( AVX_chEst_re, temp_re1);
                                AVX_chEst_im = _mm256_add_ps( AVX_chEst_im, temp_im1);

                                _mm256_storeu_ps(chEst_rePtr , AVX_chEst_re);
                                _mm256_storeu_ps(chEst_imPtr , AVX_chEst_im);

                                AVX_chEst_rePtr = AVX_chEst_rePtr + 1;
                                AVX_chEst_imPtr = AVX_chEst_imPtr + 1;

                                chEst_rePtr += 8;
                                chEst_imPtr += 8;
                            }
                        }
                    }
                
                }
            }
        }
    }
    
    
    #ifdef DEBUG
        FILE *fp1;
        // printf("segPart: %u\n", segPart);
        if( segPart == 1)
        {
            fp1 = fopen("IOs/modulationRemovalSeg1_mmse_PerTone.txt", "w");
        }
        else if( segPart == 2)
        {
            fp1 = fopen("IOs/modulationRemovalSeg2_mmse_PerTone.txt", "w");
        }
        else if( segPart == 3)
        {
            fp1 = fopen("IOs/modulationRemovalSeg3_mmse_PerTone.txt", "w");
        }
        else if( segPart == 4)
        {
            fp1 = fopen("IOs/modulationRemovalSeg4_mmse_PerTone.txt", "w");
        }

        for(int idxAnt = 0 ;idxAnt<pNrPuschInParams->nAntennaPorts;idxAnt++)
        {
            for(int idxLyr = 0;idxLyr<pNrPuschInParams->nNrOfLayers;idxLyr++)
            {
                for(int idxScTest = 0; idxScTest < 3300; idxScTest++)
                {
                    fprintf(fp1 , "%f\n",  pNrPuschOutParams->lsEst[idxAnt][idxLyr][idxScTest]);
                }
            }
        }
        fclose(fp1);

        FILE *fp2, *pf2;
        if( segPart == 1)
        {
            fp2 = fopen("IOs/estimatesRealInterpolatedSeg1_mmse_PerTone.txt", "w");
            pf2 = fopen("IOs/estimatesImagInterpolatedSeg1_mmse_PerTone.txt", "w");
        }
        else if( segPart == 2)
        {
            fp2 = fopen("IOs/estimatesRealInterpolatedSeg2_mmse_PerTone.txt", "w");
            pf2 = fopen("IOs/estimatesImagInterpolatedSeg2_mmse_PerTone.txt", "w");
        }
        else if( segPart == 3)
        {
            fp2 = fopen("IOs/estimatesRealInterpolatedSeg3_mmse_PerTone.txt", "w");
            pf2 = fopen("IOs/estimatesImagInterpolatedSeg3_mmse_PerTone.txt", "w");
        }
        else //if( segPart == 4)
        {
            fp2 = fopen("IOs/estimatesRealInterpolatedSeg4_mmse_PerTone.txt", "w");
            pf2 = fopen("IOs/estimatesImagInterpolatedSeg4_mmse_PerTone.txt", "w");
        }

        for(int idxAntTest = 0; idxAntTest < pNrPuschInParams->nAntennaPorts; idxAntTest++)
        {
            for(int idxLyrTest = 0; idxLyrTest < pNrPuschInParams->nNrOfLayers; idxLyrTest++)
            {
                for(int idxScTest = 0; idxScTest < 3300; idxScTest++)
                {
                    fprintf(fp2, "%f\n",\
                    pNrPuschOutParams->chEst_re[idxAntTest][idxLyrTest][idxScTest]);
                    fprintf(pf2, "%f\n",\
                    pNrPuschOutParams->chEst_im[idxAntTest][idxLyrTest][idxScTest]);
                }
            }
        }
        fclose(fp2);
        fclose(pf2);
    #endif
}


/*void nr_pusch_mmse_perTone_DS_equ_avx2(wnUInt8 nSeg,
                                    wnUInt8 segPart,
                                    wnUInt8* startSym,
                                    wnUInt8* endSym,
                                    commonUlConfigT* pNrUlCommonParams,
                                    puschConfigT* pNrPuschInParams,
                                    P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams)
{
    //Calculate noise variance.
    nr_pusch_mmse_perTone_nVarCalc_avx2(segPart, pNrUlCommonParams, pNrPuschInParams, pNrPuschOutParams);

    printf("nVar_overAll: %f\n", pNrPuschOutParams->nVar_overAll);

    wnUInt8 idxAlloc, idxSym, idx;
    wnUInt8 idxLyr, idxAnt;
    wnUInt8 remTones, dmrsSym;
    wnUInt8 nRxAnt = pNrPuschInParams->nAntennaPorts;
    wnUInt8 nLyr = pNrPuschInParams->nNrOfLayers;
    wnUInt8 nSbAlloc = pNrPuschInParams->nSbAlloc;
    wnUInt8 counterOf2, equaSym;
    wnUInt16 mseCnt = 0, idxTone;
    wnUInt16 total_Prb = 0;
    wnUInt16 toneStartIdx, toneEndIdx;
    wnUInt16 ntones_alloc, ntones_total, completedTones;
    wnUInt32 layerMapIdx, temp1;
    
    VECT_T H_re[nLyr][nRxAnt];
	VECT_T H_im[nLyr][nRxAnt];
	VECT_T HHh_re[nLyr][nLyr];
	VECT_T HHh_im[nLyr][nLyr];
	VECT_T invHHh_re[nLyr][nLyr];
	VECT_T invHHh_im[nLyr][nLyr];
	VECT_T weights_re[nRxAnt][nLyr];
	VECT_T weights_im[nRxAnt][nLyr];
	VECT_T Y_vect_re[1][nRxAnt];
	VECT_T Y_vect_im[1][nRxAnt];
	VECT_T Y_tilda_re[1][nLyr];
	VECT_T Y_tilda_im[1][nLyr];
	VECT_T nVar_vec = VECT_SET1(pNrPuschOutParams->nVar_overAll);
    g_negate = VECT_SET1(-0.0);

    wnFlt one_By_total_Prb = 0;
    wnFlt firstToneSNR, sum;
    
    
    //Get DMRS Symbol Index.
    if(segPart == SEG1)
    {
        dmrsSym = pNrPuschInParams->dmrs1_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG2)
    {
        dmrsSym = pNrPuschInParams->dmrs2_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG3)
    {
        dmrsSym = pNrPuschInParams->dmrs3_sym - pNrPuschInParams->nStartSymbolIndex;
    }
    else if(segPart == SEG4)
    {
        dmrsSym = pNrPuschInParams->dmrs4_sym - pNrPuschInParams->nStartSymbolIndex;
    }


    for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
    {
        total_Prb = total_Prb + pNrPuschInParams->nPrb[idxAlloc];
    }
    one_By_total_Prb = 1.0/total_Prb;

    for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
    {
        pNrPuschOutParams->oneBymse_avg_re[idxLyr] = 0;
    }

    ntones_total = 0;
    for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
    {
        completedTones = ntones_total;
        ntones_alloc = pNrPuschInParams->nPrb[idxAlloc]*12;
        ntones_total += ntones_alloc;

        toneStartIdx = pNrPuschInParams->rbStart[idxAlloc]*12;
        toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);

        counterOf2 = 1;
        for(idxTone = toneStartIdx; idxTone < toneEndIdx; idxTone += 8)
        {
            //Load Channel Estimates for MMSE weight's calculation.
            for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
            {
                for(idxAnt = 0;idxAnt < nRxAnt;idxAnt++)
                {
                    H_re[idxLyr][idxAnt] = \
                    _mm256_loadu_ps(&pNrPuschOutParams->chEst_re[idxAnt][idxLyr][idxTone]);
                    H_im[idxLyr][idxAnt] = \
                    _mm256_loadu_ps(&pNrPuschOutParams->chEst_im[idxAnt][idxLyr][idxTone]);
                }
            }

            if((idxTone + 8) >= toneEndIdx)
            {
                remTones = toneEndIdx - idxTone;

                for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
                {
                    for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                    {
                        for(idx = remTones; idx < 8; idx++)
                        {
                            H_re[idxLyr][idxAnt][idx] = 1;
                            H_im[idxLyr][idxAnt][idx] = 1;
                        }
                    }
                }
            }
            else
            {
                remTones = 8;
            }

            // HH'
            MATRIXMULTTH(H_re, H_im, HHh_re, HHh_im, nLyr, nRxAnt, nLyr);

            for(idx = 0; idx < nLyr; idx++)
            {
                HHh_re[idx][idx] = VECT_ADD(nVar_vec, HHh_re[idx][idx]);
            }


            // Inverse of (HH' + (epsilon)I)
            if(nLyr == 1)
            {
                InverseMatrix1x1HTranspose(HHh_re, HHh_im, invHHh_re, invHHh_im);
            }
            else if(nLyr == 2)
            {
                InverseMatrix2x2H(HHh_re, HHh_im, invHHh_re, invHHh_im);
            }
            else if(nLyr == 4)
            {
                InverseMatrix4x4H(HHh_re, HHh_im, invHHh_re, invHHh_im);
            }

            #ifdef SNR24 
            //SNR per two PRBs.
            if(counterOf2 == 1)
            {
                for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
                {
                    firstToneSNR = 1/invHHh_re[idxLyr][idxLyr][0];
                    pNrPuschOutParams->oneBymse_re[idxLyr][mseCnt] = firstToneSNR;

                    sum = pNrPuschOutParams->oneBymse_avg_re[idxLyr];
                    sum = sum + firstToneSNR;
                    pNrPuschOutParams->oneBymse_avg_re[idxLyr] = sum;
                }
                mseCnt += 1;
                counterOf2 += 1;
            }
            else
            {
                if(counterOf2 == 3)
                {
                    counterOf2 = 1;
                }
                else
                {
                    counterOf2 += 1;
                }
            }
            #else
                // Calculate SNR per Tone.
                for(idx = 0; idx < remTones; idx++)
                {
                    for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
                    {
                        firstToneSNR = 1/invHHh_re[idxLyr][idxLyr][idx];
                        pNrPuschOutParams->oneBymse_re[idxLyr][mseCnt] = firstToneSNR;

                        sum = pNrPuschOutParams->oneBymse_avg_re[idxLyr];
                        sum = sum + firstToneSNR;
                        pNrPuschOutParams->oneBymse_avg_re[idxLyr] = sum;

                    }
                    mseCnt += 1;
                }
            #endif
            

            // Multiplication of H' and inverse of (HH' + (epsilon)I)
            MATRIXMULTHTTRANSPOSE2(H_re, H_im, invHHh_re, invHHh_im, weights_re, weights_im, nRxAnt, nLyr, nLyr);

            #ifdef DEBUG
                // if(idxTone == toneStartIdx)
                // {
                //     printf("H:\n");
                //     for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++) 
                //     {
                //         for(idx = 0; idx < nLyr; idx++) 
                //         {
                //             printf("%f%+f*1i ", H_re[idx][idxAnt][0], H_im[idx][idxAnt][0]);
                //         }
                //         printf("\n");
                //     }

                //     printf("invHHh:\n");
                //     for(idxLyr = 0; idxLyr < nLyr; idxLyr++) 
                //     {
                //         for(idx = 0; idx < nLyr; idx++) 
                //         {
                //             printf("%f%+f*1i ", invHHh_re[idx][idxLyr][0], invHHh_im[idx][idxLyr][0]);
                //         }
                //         printf("\n");
                //     }
                // }
            #endif
            
            equaSym = 0;
            //Use same MMSE weight's across data symbols for equalisation.
            for(idxSym = startSym[segPart - 1]; idxSym <= endSym[segPart - 1]; idxSym++)
            {
                if(idxSym != dmrsSym)
                {
                    temp1 = equaSym*pNrPuschOutParams->total_nPRB*12*nLyr + completedTones*nLyr;
                    
                    //Load received data samples.
                    if( (idxTone + 8) >= toneEndIdx )
                    {
                        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                        {
                            Y_vect_re[0][idxAnt] = _mm256_set1_ps(0);
                            Y_vect_im[0][idxAnt] = _mm256_set1_ps(0);
                            remTones = toneEndIdx - idxTone;

                            for(idx = 0; idx < remTones; idx++)
                            {
                                Y_vect_re[0][idxAnt][idx] = rxFdSamples[idxAnt][idxSym][idxTone+idx];
                                Y_vect_im[0][idxAnt][idx] = rxFdSamples[idxAnt][idxSym][MAX_NUM_SC+idxTone+idx];
                            }
                        }
                    }
                    else
                    {
                        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                        {
                            Y_vect_re[0][idxAnt] = VECT_SET(rxFdSamples[idxAnt][idxSym][idxTone+7],
                                                            rxFdSamples[idxAnt][idxSym][idxTone+6],
                                                            rxFdSamples[idxAnt][idxSym][idxTone+5],
                                                            rxFdSamples[idxAnt][idxSym][idxTone+4],
                                                            rxFdSamples[idxAnt][idxSym][idxTone+3],
                                                            rxFdSamples[idxAnt][idxSym][idxTone+2],
                                                            rxFdSamples[idxAnt][idxSym][idxTone+1],
                                                            rxFdSamples[idxAnt][idxSym][idxTone] );

                            Y_vect_im[0][idxAnt] = VECT_SET(rxFdSamples[idxAnt][idxSym][MAX_NUM_SC+idxTone+7],
                                                            rxFdSamples[idxAnt][idxSym][MAX_NUM_SC+idxTone+6],
                                                            rxFdSamples[idxAnt][idxSym][MAX_NUM_SC+idxTone+5],
                                                            rxFdSamples[idxAnt][idxSym][MAX_NUM_SC+idxTone+4],
                                                            rxFdSamples[idxAnt][idxSym][MAX_NUM_SC+idxTone+3],
                                                            rxFdSamples[idxAnt][idxSym][MAX_NUM_SC+idxTone+2],
                                                            rxFdSamples[idxAnt][idxSym][MAX_NUM_SC+idxTone+1],
                                                            rxFdSamples[idxAnt][idxSym][MAX_NUM_SC+idxTone] );

                        }
                    }
                    
                    //MMSE Equalisation.
                    MATRIXMULNN(Y_vect_re, Y_vect_im, weights_re, weights_im, Y_tilda_re, Y_tilda_im, 1, nRxAnt, nLyr);

                    #ifdef DEBUG
                        // if(equaSym == 0 && idxTone == toneStartIdx)
                        // {
                        //     printf("Y_vect:\n");
                        //     for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                        //     {
                        //         printf("%f%+f*1i \n", Y_vect_re[0][idxAnt][0], Y_vect_im[0][idxAnt][0]);
                        //     }
                            
                        //     printf("weights:\n");
                        //     for(idx = 0; idx < nLyr; idx++) 
                        //     {
                        //         for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                        //         {
                        //             printf("%f%+f*1i ", weights_re[idxAnt][idx][0], weights_im[idxAnt][idx][0]);
                        //         }
                        //         printf("\n");
                        //     }

                        //     printf("Y_tilda:\n");
                        //     for(idx = 0; idx < nLyr; idx++)
                        //     {
                        //         printf("%f%+f*1i \n", Y_tilda_re[0][idx][0], Y_tilda_im[0][idx][0]);
                        //     }
                        // }
                    #endif

                    for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
                    {
                        for(idx = 0; idx < remTones; idx++)
                        {
                            if(pNrPuschInParams->FOEnable == 1)
                            {
                                pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][idxTone+idx] = Y_tilda_re[0][idxLyr][idx];
                                pNrPuschOutParams->equaOutSamples[idxLyr][equaSym][MAX_NUM_SC+idxTone+idx] = Y_tilda_im[0][idxLyr][idx];
                            }
                            else
                            {
                                // Directly stores data in their respective Layer Demapped Positions
                                layerMapIdx = temp1 + (idxTone+idx-toneStartIdx)*nLyr + idxLyr;
                                pNrPuschOutParams->layerDemapperOutReal[layerMapIdx] = Y_tilda_re[0][idxLyr][idx];
                                pNrPuschOutParams->layerDemapperOutImag[layerMapIdx] = Y_tilda_im[0][idxLyr][idx];
                            }
                        }
                    }
                    equaSym = equaSym + 1;
                }
                else
                {
                    idxSym = idxSym+1;
                }
            }
        }
    }

    totalLength = 12*pNrPuschOutParams->total_nPRB*12*nLyr;

    // Calculate mean of all SNR values.
    #ifdef SNR24
        for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
        {
            pNrPuschOutParams->oneBymse_avg_re[idxLyr] = \
                pNrPuschOutParams->oneBymse_avg_re[idxLyr]*one_By_total_Prb*0.5;
        }
    #else
        wnFlt oneByNumOfTones = one_By_total_Prb/12;
        for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
        {
            pNrPuschOutParams->oneBymse_avg_re[idxLyr] = \
            pNrPuschOutParams->oneBymse_avg_re[idxLyr]*oneByNumOfTones;
        }
    #endif



    #ifdef DEBUG
        FILE *fpMSE;
        FILE *fpMSEAvg;
        FILE *fpLD, *pfLD;
        wnUInt8 nSymbs;

        if( segPart == SEG1)
        {
            fpMSE = fopen("IOs/mse_perTone_Seg1.txt", "w");
            fpMSEAvg = fopen("IOs/mseAvg_perTone_Seg1.txt", "w");
        }
        else if(segPart == SEG2)
        {
            fpMSE = fopen("IOs/mse_perTone_Seg2.txt", "w");
            fpMSEAvg = fopen("IOs/mseAvg_perTone_Seg2.txt", "w");
        }
        else if(segPart == SEG3)
        {
            fpMSE = fopen("IOs/mse_perTone_Seg3.txt", "w");
            fpMSEAvg = fopen("IOs/mseAvg_perTone_Seg3.txt", "w");
        }
        else if(segPart == SEG4)
        {
            fpMSE = fopen("IOs/mse_perTone_Seg4.txt", "w");
            fpMSEAvg = fopen("IOs/mseAvg_perTone_Seg4.txt", "w");
        }

        #ifdef SNR24
            for(int lyrIdx = 0;lyrIdx<pNrPuschInParams->nNrOfLayers;lyrIdx++)
            {
                for(int toneIdx = 0;toneIdx<pNrPuschInParams->nRBSize;toneIdx++)
                {
                    fprintf(fpMSE, "%f\n ", pNrPuschOutParams->oneBymse_re[lyrIdx][toneIdx]);
                }
            }
        #else
            for(int lyrIdx = 0;lyrIdx<pNrPuschInParams->nNrOfLayers;lyrIdx++)
            {
                for(int toneIdx = 0;toneIdx<pNrPuschInParams->nRBSize*12;toneIdx++)
                {
                    fprintf(fpMSE, "%f\n ", pNrPuschOutParams->oneBymse_re[lyrIdx][toneIdx]);
                }
            }
        #endif

        for(int lyrIdx = 0;lyrIdx<pNrPuschInParams->nNrOfLayers;lyrIdx++)
        {
            fprintf(fpMSEAvg, "%f\n ", pNrPuschOutParams->oneBymse_avg_re[lyrIdx]);
        }

        fclose(fpMSE);
        fclose(fpMSEAvg);

        nSymbs = endSym[segPart-1] - startSym[segPart-1]-(pNrPuschInParams->nNrOfDMRSSymbols-1);
        if( segPart == SEG1)
        {
            fpLD = fopen("IOs/layerDemapperReal_mmse_perToneSeg1.txt", "w");
            pfLD = fopen("IOs/layerDemapperImag_mmse_perToneSeg1.txt", "w");
        }
        else if(segPart == SEG2)
        {
            fpLD = fopen("IOs/layerDemapperReal_mmse_perToneSeg2.txt", "w");
            pfLD = fopen("IOs/layerDemapperImag_mmse_perToneSeg2.txt", "w");
        }
        else if(segPart == SEG3)
        {
            fpLD = fopen("IOs/layerDemapperReal_mmse_perToneSeg3.txt", "w");
            pfLD = fopen("IOs/layerDemapperImag_mmse_perToneSeg3.txt", "w");
        }
        else if(segPart == SEG4)
        {
            fpLD = fopen("IOs/layerDemapperReal_mmse_perToneSeg4.txt", "w");
            pfLD = fopen("IOs/layerDemapperImag_mmse_perToneSeg4.txt", "w");
        }

        // printf("nSamples: %u\n", pNrPuschOutParams->total_nPRB*12*nLyr*nSymbs );// *12*nLyr*nSymbs
        for(int I1 = 0;I1< pNrPuschOutParams->total_nPRB*12*nLyr*nSymbs;I1++)
        {
            fprintf(fpLD, "%f\n", pNrPuschOutParams->layerDemapperOutReal[I1]);
            fprintf(pfLD, "%f\n", pNrPuschOutParams->layerDemapperOutImag[I1]);
        }
        fclose(fpLD);
        fclose(pfLD);

        if(pNrPuschInParams->FOEnable == 1)
        {
            // printf("equaOutSamples[0][0][0]: %f\n", pNrPuschOutParams->equaOutSamples[0][0][0]);
            FILE *eqRePtr, *eqImPtr;
            if(segPart == SEG1)
            {
                eqRePtr = fopen("IOs/equalisedRealSeg1.txt", "w");
                eqImPtr = fopen("IOs/equalisedImagSeg1.txt", "w");
            }
            else if(segPart == SEG2)
            {
                eqRePtr = fopen("IOs/equalisedRealSeg2.txt", "w");
                eqImPtr = fopen("IOs/equalisedImagSeg2.txt", "w");
            }
            else if(segPart == SEG3)
            {
                eqRePtr = fopen("IOs/equalisedRealSeg3.txt", "w");
                eqImPtr = fopen("IOs/equalisedImagSeg3.txt", "w");
            }
            else if(segPart == SEG4)
            {
                eqRePtr = fopen("IOs/equalisedRealSeg4.txt", "w");
                eqImPtr = fopen("IOs/equalisedImagSeg4.txt", "w");
            }

            for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
            {
                toneStartIdx = pNrPuschInParams->rbStart[idxAlloc]*12;
                toneEndIdx = toneStartIdx + (pNrPuschInParams->nPrb[idxAlloc] * 12);
                for(idxSym = 0;idxSym<nSymbs;idxSym++)
                {
                    for(idxTone = toneStartIdx; idxTone < toneEndIdx; idxTone++)
                    {
                        for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
                        {
                            fprintf(eqRePtr, "%f\n", pNrPuschOutParams->equaOutSamples[idxLyr][idxSym][idxTone]);
                            fprintf(eqImPtr, "%f\n", pNrPuschOutParams->equaOutSamples[idxLyr][idxSym][MAX_NUM_SC+idxTone]);
                        }
                    }
                }
            }
            fclose(eqRePtr);
            fclose(eqImPtr);
        }
    #endif
}//*/


wnVoid nr_pusch_mmse_TF_IFFT(wnUInt8 nSeg,
                            wnUInt8 segPart,
                            wnUInt8* startSym,
                            wnUInt8* endSym,
                            puschConfigT* pNrPuschInParams,
                            P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams)
{
    wnUInt8 dataSymbols, idxSym;
    wnFlt scaleUp;
    wnUInt16 toneIdx, nTonesPerSymb, nRBSize;
    fftwf_complex *IFFTIn, *IFFTOut;
    fftwf_plan planIFFT;

    if(segPart == SEG1)
    {
        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
        {
            dataSymbols = endSym[0]-startSym[0];
        }
        else //if(pNrPuschInParams->nNrOfDMRSSymbols == 2)
        {
            dataSymbols = endSym[0]-startSym[0]-1;
        }
    }
    else if(segPart == SEG2)
    {
        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
        {
            dataSymbols = endSym[1]-startSym[1];
        }
        else //if(pNrPuschInParams->nNrOfDMRSSymbols == 2)
        {
            dataSymbols = endSym[1]-startSym[1]-1;
        }
    }
    else if(segPart == SEG3)
    {
        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
        {
            dataSymbols = endSym[2]-startSym[2];
        }
        else //if(pNrPuschInParams->nNrOfDMRSSymbols == 2)
        {
            dataSymbols = endSym[2]-startSym[2]-1;
        }
    }
    else //if(segPart == SEG4)
    {
        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
        {
            dataSymbols = endSym[3]-startSym[3];
        }
        else //if(pNrPuschInParams->nNrOfDMRSSymbols == 2)
        {
            dataSymbols = endSym[3]-startSym[3]-1;
        }
    }
    
    nRBSize = pNrPuschInParams->nRBSize;    
    nTonesPerSymb = nRBSize*12;
    scaleUp = sqrt(nTonesPerSymb);
    IFFTIn = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * nTonesPerSymb);
    IFFTOut = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * nTonesPerSymb);
    cmplx inp;
    
    // Plan for IFFT Operation is created.
    planIFFT = fftwf_plan_dft_1d(nTonesPerSymb, IFFTIn, IFFTOut, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    for(idxSym = 0;idxSym<dataSymbols;idxSym++)
    {
        for(toneIdx = 0;toneIdx<nTonesPerSymb;toneIdx++)
        {
            inp.real = pNrPuschOutParams->layerDemapperOutReal[idxSym*nTonesPerSymb+toneIdx];
            inp.imag = pNrPuschOutParams->layerDemapperOutImag[idxSym*nTonesPerSymb+toneIdx];
            memcpy(&IFFTIn[toneIdx], &inp, sizeof(fftwf_complex));
        }

        fftwf_execute(planIFFT);
        
        for(toneIdx = 0;toneIdx<nTonesPerSymb;toneIdx++)
        {
            memcpy(&inp, &IFFTOut[toneIdx],  sizeof(fftwf_complex));
            pNrPuschOutParams->layerDemapperOutReal[idxSym*nTonesPerSymb+toneIdx] = 1/scaleUp*inp.real;
            pNrPuschOutParams->layerDemapperOutImag[idxSym*nTonesPerSymb+toneIdx] = 1/scaleUp*inp.imag;
        }
    }

    fftwf_free(IFFTIn);
    fftwf_free(IFFTOut);
    fftwf_destroy_plan(planIFFT);
}


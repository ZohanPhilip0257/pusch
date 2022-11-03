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
  *         for NR PUSCH Channel Estimation and Equalization using MMSE-IRC Algorithm.
  * @file    : nr_pusch_mmse_irc.c
  * @ingroup : nr_pusch
  * @author  : MIRZA SAMI BAIG
  **/

#include "nr_pusch_mmse_irc.h"

/*void nr_pusch_mmse_irc_perPrb_est_avx2( int symPart,
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
    __m256 Rnn_inv_re[nRxAnt][nRxAnt], Rnn_inv_im[nRxAnt][nRxAnt];
    __m256 RnnTemp_re[nRxAnt][nRxAnt], RnnTemp_im[nRxAnt][nRxAnt];

    complex h1Avg[nRxAnt][8], h2Avg[nRxAnt][8], hAvg;
    complex Rnn[nRxAnt][nRxAnt], Inf[nRxAnt];
    float prgRnn_re[8], prgRnn_im[8];

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

    g_negate = VECT_SET1(-0.0);

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

            //Time offset estimation (TOE) on Even Tones
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

            //Time offset estimation (TOE) on Odd Tones
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

            //Over all Time offset estimation (TOE)
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

            //Channel Estimation & Covariance Matrix (Rnn) computation on Even Tone
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
                if((pNrPuschInParams->nullTone_enable == 1) || (nOddTonesLyr == 0))
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

            //Channel Estimation & Covariance Matrix (Rnn) computation on Odd Tone
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

        if(COV_GRP_FACTOR == 1)
        {
            for(idxPrb = startPrbIdx; idxPrb <= endPrbIdx; idxPrb+=8)
            {
                for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                {
                    for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                    {
                        if((pNrPuschInParams->nullTone_enable == 1) || ((pNrPuschOutParams->nEvenTonesLyr > 0) && (pNrPuschOutParams->nOddTonesLyr > 0)))
                        {
                            prgRnn_re[0] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb] ) * 0.5;
                            prgRnn_re[1] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+1] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+1] ) * 0.5;
                            prgRnn_re[2] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+2] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+2] ) * 0.5;
                            prgRnn_re[3] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+3] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+3] ) * 0.5;
                            prgRnn_re[4] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+4] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+4] ) * 0.5;
                            prgRnn_re[5] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+5] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+5] ) * 0.5;
                            prgRnn_re[6] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+6] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+6] ) * 0.5;
                            prgRnn_re[7] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+7] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+7] ) * 0.5;

                            prgRnn_im[0] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb] ) * 0.5;
                            prgRnn_im[1] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+1] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+1] ) * 0.5;
                            prgRnn_im[2] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+2] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+2] ) * 0.5;
                            prgRnn_im[3] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+3] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+3] ) * 0.5;
                            prgRnn_im[4] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+4] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+4] ) * 0.5;
                            prgRnn_im[5] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+5] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+5] ) * 0.5;
                            prgRnn_im[6] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+6] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+6] ) * 0.5;
                            prgRnn_im[7] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+7] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+7] ) * 0.5;

                            Rnn_re[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_re);
                            Rnn_im[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_im);
                        }
                        else if (pNrPuschOutParams->nEvenTonesLyr > 0)
                        {
                            prgRnn_re[0] = pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb];
                            prgRnn_re[1] = pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+1];
                            prgRnn_re[2] = pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+2];
                            prgRnn_re[3] = pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+3];
                            prgRnn_re[4] = pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+4];
                            prgRnn_re[5] = pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+5];
                            prgRnn_re[6] = pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+6];
                            prgRnn_re[7] = pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+7];

                            prgRnn_im[0] = pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb];
                            prgRnn_im[1] = pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+1];
                            prgRnn_im[2] = pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+2];
                            prgRnn_im[3] = pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+3];
                            prgRnn_im[4] = pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+4];
                            prgRnn_im[5] = pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+5];
                            prgRnn_im[6] = pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+6];
                            prgRnn_im[7] = pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+7];

                            Rnn_re[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_re);
                            Rnn_im[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_im);

                        }
                        else
                        {
                            prgRnn_re[0] = pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb];
                            prgRnn_re[1] = pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+1];
                            prgRnn_re[2] = pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+2];
                            prgRnn_re[3] = pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+3];
                            prgRnn_re[4] = pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+4];
                            prgRnn_re[5] = pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+5];
                            prgRnn_re[6] = pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+6];
                            prgRnn_re[7] = pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+7];

                            prgRnn_im[0] = pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb];
                            prgRnn_im[1] = pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+1];
                            prgRnn_im[2] = pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+2];
                            prgRnn_im[3] = pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+3];
                            prgRnn_im[4] = pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+4];
                            prgRnn_im[5] = pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+5];
                            prgRnn_im[6] = pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+6];
                            prgRnn_im[7] = pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+7];

                            Rnn_re[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_re);
                            Rnn_im[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_im);
                        }
                    }
                }

                if(nRxAnt == 1)
                {
                    InverseMatrix1x1HTranspose(Rnn_re, Rnn_im, Rnn_inv_re, Rnn_inv_im);
                }
                else if(nRxAnt == 2)
                {
                    InverseMatrix2x2H(Rnn_re, Rnn_im, Rnn_inv_re, Rnn_inv_im);
                }
                else if(nRxAnt == 4)
                {
                    InverseMatrix4x4H(Rnn_re, Rnn_im, Rnn_inv_re, Rnn_inv_im);
                }
                else if(nRxAnt == 8)
                {
                    InverseMatrix8x8H(Rnn_re, Rnn_im, Rnn_inv_re, Rnn_inv_im);
                }

                for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                {
                    for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                    {
                        _mm256_storeu_ps(pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol], Rnn_inv_re[idxRow][idxCol]);
                        _mm256_storeu_ps(pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol], Rnn_inv_im[idxRow][idxCol]);
                    }
                }

            }

        }
        else if(COV_GRP_FACTOR == 2)
        {
            for(idxPrb = startPrbIdx; idxPrb <= endPrbIdx; idxPrb+=16)
            {
                for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                {
                    for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                    {
                        if((pNrPuschInParams->nullTone_enable == 1) || ((pNrPuschOutParams->nEvenTonesLyr > 0) && (pNrPuschOutParams->nOddTonesLyr > 0)))
                        {
                            prgRnn_re[0] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb] +
                                             pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+1] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+1] ) * 0.25;
                            prgRnn_re[1] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+2] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+2] +
                                             pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+3] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+3] ) * 0.25;
                            prgRnn_re[2] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+4] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+4] +
                                             pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+5] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+5] ) * 0.25;
                            prgRnn_re[3] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+6] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+6] +
                                             pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+7] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+7] ) * 0.25;
                            prgRnn_re[4] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+8] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+8] +
                                             pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+9] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+9] ) * 0.25;
                            prgRnn_re[5] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+10] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+10] +
                                             pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+11] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+11] ) * 0.25;
                            prgRnn_re[6] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+12] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+12] +
                                             pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+13] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+13] ) * 0.25;
                            prgRnn_re[7] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+14] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+14] +
                                             pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+15] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+15] ) * 0.25;

                            prgRnn_im[0] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb] +
                                             pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+1] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+1] ) * 0.25;
                            prgRnn_im[1] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+2] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+2] +
                                             pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+3] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+3] ) * 0.25;
                            prgRnn_im[2] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+4] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+4] +
                                             pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+5] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+5] ) * 0.25;
                            prgRnn_im[3] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+6] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+6] +
                                             pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+7] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+7] ) * 0.25;
                            prgRnn_im[4] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+8] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+8] +
                                             pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+9] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+9] ) * 0.25;
                            prgRnn_im[5] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+10] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+10] +
                                             pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+11] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+11] ) * 0.25;
                            prgRnn_im[6] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+12] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+12] +
                                             pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+13] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+13] ) * 0.25;
                            prgRnn_im[7] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+14] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+14] +
                                             pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+15] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+15] ) * 0.25;

                            Rnn_re[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_re);
                            Rnn_im[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_im);
                        }
                        else if (pNrPuschOutParams->nEvenTonesLyr > 0)
                        {
                            prgRnn_re[0] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb] + pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+1] ) * 0.5;
                            prgRnn_re[1] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+2] + pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+3] ) * 0.5;
                            prgRnn_re[2] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+4] + pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+5] ) * 0.5;
                            prgRnn_re[3] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+6] + pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+7] ) * 0.5;
                            prgRnn_re[4] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+8] + pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+9] ) * 0.5;
                            prgRnn_re[5] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+10] + pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+11] ) * 0.5;
                            prgRnn_re[6] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+12] + pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+13] ) * 0.5;
                            prgRnn_re[7] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+14] + pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+15] ) * 0.5;

                            prgRnn_im[0] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb] + pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+1] ) * 0.5;
                            prgRnn_im[1] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+2] + pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+3] ) * 0.5;
                            prgRnn_im[2] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+4] + pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+5] ) * 0.5;
                            prgRnn_im[3] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+6] + pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+7] ) * 0.5;
                            prgRnn_im[4] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+8] + pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+9] ) * 0.5;
                            prgRnn_im[5] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+10] + pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+11] ) * 0.5;
                            prgRnn_im[6] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+12] + pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+13] ) * 0.5;
                            prgRnn_im[7] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+14] + pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+15] ) * 0.5;

                            Rnn_re[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_re);
                            Rnn_im[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_im);

                        }
                        else
                        {
                            prgRnn_re[0] = ( pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+1] ) * 0.5;
                            prgRnn_re[1] = ( pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+2] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+3] ) * 0.5;
                            prgRnn_re[2] = ( pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+4] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+5] ) * 0.5;
                            prgRnn_re[3] = ( pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+6] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+7] ) * 0.5;
                            prgRnn_re[4] = ( pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+8] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+9] ) * 0.5;
                            prgRnn_re[5] = ( pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+10] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+11] ) * 0.5;
                            prgRnn_re[6] = ( pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+12] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+13] ) * 0.5;
                            prgRnn_re[7] = ( pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+14] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+15] ) * 0.5;

                            prgRnn_im[0] = ( pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+1] ) * 0.5;
                            prgRnn_im[1] = ( pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+2] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+3] ) * 0.5;
                            prgRnn_im[2] = ( pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+4] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+5] ) * 0.5;
                            prgRnn_im[3] = ( pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+6] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+7] ) * 0.5;
                            prgRnn_im[4] = ( pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+8] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+9] ) * 0.5;
                            prgRnn_im[5] = ( pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+10] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+11] ) * 0.5;
                            prgRnn_im[6] = ( pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+12] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+13] ) * 0.5;
                            prgRnn_im[7] = ( pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+14] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+15] ) * 0.5;

                            Rnn_re[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_re);
                            Rnn_im[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_im);

                        }

                    }
                }

                if(nRxAnt == 1)
                {
                    InverseMatrix1x1HTranspose(Rnn_re, Rnn_im, Rnn_inv_re, Rnn_inv_im);
                }
                else if(nRxAnt == 2)
                {
                    InverseMatrix2x2H(Rnn_re, Rnn_im, Rnn_inv_re, Rnn_inv_im);
                }
                else if(nRxAnt == 4)
                {
                    InverseMatrix4x4H(Rnn_re, Rnn_im, Rnn_inv_re, Rnn_inv_im);
                }
                else if(nRxAnt == 8)
                {
                    InverseMatrix8x8H(Rnn_re, Rnn_im, Rnn_inv_re, Rnn_inv_im);
                }

                for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                {
                    for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                    {
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb] = Rnn_inv_re[idxRow][idxCol][0];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+1] = Rnn_inv_re[idxRow][idxCol][0];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+2] = Rnn_inv_re[idxRow][idxCol][1];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+3] = Rnn_inv_re[idxRow][idxCol][1];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+4] = Rnn_inv_re[idxRow][idxCol][2];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+5] = Rnn_inv_re[idxRow][idxCol][2];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+6] = Rnn_inv_re[idxRow][idxCol][3];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+7] = Rnn_inv_re[idxRow][idxCol][3];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+8] = Rnn_inv_re[idxRow][idxCol][4];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+9] = Rnn_inv_re[idxRow][idxCol][4];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+10] = Rnn_inv_re[idxRow][idxCol][5];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+11] = Rnn_inv_re[idxRow][idxCol][5];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+12] = Rnn_inv_re[idxRow][idxCol][6];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+13] = Rnn_inv_re[idxRow][idxCol][6];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+14] = Rnn_inv_re[idxRow][idxCol][7];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+15] = Rnn_inv_re[idxRow][idxCol][7];

                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb] = Rnn_inv_im[idxRow][idxCol][0];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+1] = Rnn_inv_im[idxRow][idxCol][0];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+2] = Rnn_inv_im[idxRow][idxCol][1];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+3] = Rnn_inv_im[idxRow][idxCol][1];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+4] = Rnn_inv_im[idxRow][idxCol][2];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+5] = Rnn_inv_im[idxRow][idxCol][2];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+6] = Rnn_inv_im[idxRow][idxCol][3];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+7] = Rnn_inv_im[idxRow][idxCol][3];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+8] = Rnn_inv_im[idxRow][idxCol][4];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+9] = Rnn_inv_im[idxRow][idxCol][4];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+10] = Rnn_inv_im[idxRow][idxCol][5];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+11] = Rnn_inv_im[idxRow][idxCol][5];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+12] = Rnn_inv_im[idxRow][idxCol][6];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+13] = Rnn_inv_im[idxRow][idxCol][6];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+14] = Rnn_inv_im[idxRow][idxCol][7];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+15] = Rnn_inv_im[idxRow][idxCol][7];
                    }
                }

            }

        }

    }
}

void nr_pusch_mmse_irc_perPrb_equ_avx2( int symPart,
                                        commonUlConfigT* pNrUlCommonParams,
                                        puschConfigT* pNrPuschInParams,
                                        P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams)
{

    int idxAnt, idxLyr, idxAlloc, idxPrb, idxSym, idxSc, rIndx, cIndx, i, idx;
    int startScIdx, endScIdx, remSc;
    int startPrbIdx, endPrbIdx, CtPrbScStartIdx, avxPrbIdx, ctPrgIdx, lyrPort, lyrType;
    int printPrbIdx = 0, avxPrintIdx = 0;

    int startSym, endSym, dmrsSym;
    int nRxAnt, nLyr, nSbAlloc, nPrb_alloc, nPrb_total;

    int prbScStartIdx_re[8], prbScStartIdx_im[8];

    nRxAnt   = pNrPuschInParams->nAntennaPorts;
    nLyr     = pNrPuschInParams->nNrOfLayers;
    nSbAlloc = pNrPuschInParams->nSbAlloc;

    float Rnn_grp1_re, Rnn_grp1_im;
    float Rnn_grp2_re, Rnn_grp2_im;
    float Rnn_grp3_re, Rnn_grp3_im;
    float Rnn_grp4_re, Rnn_grp4_im;

    __m256 mid_re, mid_im;

    __m256 hMat_re[nRxAnt][nLyr], hMat_im[nRxAnt][nLyr];
    __m256 hHerMat_re[nLyr][nRxAnt], hHerMat_im[nLyr][nRxAnt];
    __m256 Rnn_re[nRxAnt][nRxAnt], Rnn_im[nRxAnt][nRxAnt];
    __m256 Rnn_inv_re[nRxAnt][nRxAnt], Rnn_inv_im[nRxAnt][nRxAnt];
    __m256 hHer_RnnInv_re[nLyr][nRxAnt], hHer_RnnInv_im[nLyr][nRxAnt];
    __m256 hHer_RnnInv_h_re[nLyr][nLyr], hHer_RnnInv_h_im[nLyr][nLyr];
    __m256 hHer_RnnInv_h_inv_re[nLyr][nLyr], hHer_RnnInv_h_inv_im[nLyr][nLyr];

    __m256 hHerRnnInvhInv_hHer_re[nLyr][nRxAnt], hHerRnnInvhInv_hHer_im[nLyr][nRxAnt];
    __m256 weights_re[nLyr][nRxAnt], weights_im[nLyr][nRxAnt];
    __m256 Y_vect_re[nRxAnt][1], Y_vect_im[nRxAnt][1];
    __m256 Y_tilda_re[nLyr][1], Y_tilda_im[nLyr][1];
    __m256 Y_cap_re, Y_cap_im;
    __m256 sq_toc_re[2][12];
	__m256 sq_toc_im[2][12];

	__m256 mseVec, mseAvgVec, mseAvgInvVec;

	g_negate = VECT_SET1(-0.0);

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

    float complex toc_overAll, sqrt_toc[12];
    float mse, mse_inv;
    float mseAvg, mseAvg_inv;
    int completedPrb = 0;
    nPrb_total = 0;
    for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
    {
        nPrb_alloc = pNrPuschInParams->nPrb[idxAlloc];
        completedPrb = nPrb_total;
        nPrb_total += nPrb_alloc;

        startPrbIdx = pNrPuschInParams->rbStart[idxAlloc];
        endPrbIdx   = startPrbIdx + (pNrPuschInParams->nPrb[idxAlloc]) - 1;
        CtPrbScStartIdx = startPrbIdx * 12;

        for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb+=8, CtPrbScStartIdx += 12)
        {
            ctPrgIdx = idxPrb / 8;

            prbScStartIdx_re[0] = idxPrb * 12;
            prbScStartIdx_im[0] = MAX_NUM_SC + prbScStartIdx_re[0];

            for(idx = 1; idx < 8; idx++)
            {
                prbScStartIdx_re[idx] = prbScStartIdx_re[idx-1] + 12;
                prbScStartIdx_im[idx] = prbScStartIdx_im[idx-1] + 12;
            }

            // Load H and H_her Matrix
            for(rIndx = 0; rIndx < nRxAnt; rIndx++)
            {
                for(cIndx = 0; cIndx < nLyr; cIndx++)
                {
                    lyrPort = pNrPuschInParams->nPortIndex[cIndx];
                    hMat_re[rIndx][cIndx] = _mm256_loadu_ps(&pNrPuschOutParams->chEst_re[rIndx][lyrPort][idxPrb]);
                    hMat_im[rIndx][cIndx] = _mm256_loadu_ps(&pNrPuschOutParams->chEst_im[rIndx][lyrPort][idxPrb]);

                    hHerMat_re[cIndx][rIndx] = _mm256_loadu_ps(&pNrPuschOutParams->chEst_re[rIndx][lyrPort][idxPrb]);
                    hHerMat_im[cIndx][rIndx] = VECT_NEG(_mm256_loadu_ps(&pNrPuschOutParams->chEst_im[rIndx][lyrPort][idxPrb]));
                }
            }

            // Load Rnn_inv Matrix
            for(rIndx = 0; rIndx < nRxAnt; rIndx++)
            {
                for(cIndx = 0; cIndx < nRxAnt; cIndx++)
                {
                    Rnn_inv_re[rIndx][cIndx] = _mm256_loadu_ps(&pNrPuschOutParams->RnnInv_prb_re[rIndx][cIndx][idxPrb]);
                    Rnn_inv_im[rIndx][cIndx] = _mm256_loadu_ps(&pNrPuschOutParams->RnnInv_prb_im[rIndx][cIndx][idxPrb]);
                }
            }

            // Perform [H_her * Rnn_inv]
            MATRIXMULH2(hHerMat_re, hHerMat_im, Rnn_inv_re, Rnn_inv_im, hHer_RnnInv_re, hHer_RnnInv_im, nLyr, nRxAnt, nRxAnt);

            // Perform [(H_her * Rnn_inv) * H]
            MATRIXMULH2(hHer_RnnInv_re, hHer_RnnInv_im, hMat_re, hMat_im, hHer_RnnInv_h_re, hHer_RnnInv_h_im, nLyr, nRxAnt, nLyr);

            if(nLyr == 1)
            {
                for(idx = 0; idx < 8; idx++)
                {
                    pNrPuschOutParams->msePerPrb[0][idxPrb + idx] = hHer_RnnInv_h_re[0][0][idx];
                }
            }

            // Perform [(H_her * Rnn_inv * H) + I]
            for(avxPrbIdx = 0; avxPrbIdx < 8; avxPrbIdx++)
            {
                for(rIndx = 0; rIndx < nLyr; rIndx++)
                {
                    hHer_RnnInv_h_re[rIndx][rIndx][avxPrbIdx] += 1;
                }
            }

            // Perform [inv((H_her * Rnn_inv * H) + I)]
            if(nLyr == 1)
            {
                InverseMatrix1x1HTranspose(hHer_RnnInv_h_re, hHer_RnnInv_h_im, hHer_RnnInv_h_inv_re, hHer_RnnInv_h_inv_im);
            }
            else if(nLyr == 2)
            {
                InverseMatrix2x2H(hHer_RnnInv_h_re, hHer_RnnInv_h_im, hHer_RnnInv_h_inv_re, hHer_RnnInv_h_inv_im);
            }
            else if(nLyr == 4)
            {
                InverseMatrix4x4H(hHer_RnnInv_h_re, hHer_RnnInv_h_im, hHer_RnnInv_h_inv_re, hHer_RnnInv_h_inv_im);
            }

            if(nLyr > 1)
            {
                for(rIndx = 0; rIndx < nLyr; rIndx++)
                {
                    for(idx = 0; idx < 8; idx++)
                    {
                        mse = hHer_RnnInv_h_inv_re[rIndx][rIndx][idx];
                        mse_inv = 1 / mse;
                        pNrPuschOutParams->msePerPrb[rIndx][idxPrb + idx] = mse_inv;
                    }
                }
            }

            // Perform [(inv((H_her * Rnn_inv * H) + I)) * H_her]
            MATRIXMULH2(hHer_RnnInv_h_inv_re, hHer_RnnInv_h_inv_im, hHerMat_re, hHerMat_im, hHerRnnInvhInv_hHer_re, hHerRnnInvhInv_hHer_im, nLyr, nLyr, nRxAnt);

            // Perform W = [(inv((H_her * Rnn_inv * H) + I) * H_her) * Rnn_inv]
            MATRIXMULH2(hHerRnnInvhInv_hHer_re, hHerRnnInvhInv_hHer_im, Rnn_inv_re, Rnn_inv_im, weights_re, weights_im, nLyr, nRxAnt, nRxAnt);

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

            // Perform Equalization on each data symbol
            int equaSym = 0;
            for(idxSym = startSym; idxSym <= endSym; idxSym++)
            {
                if(idxSym != dmrsSym)
                {
                    int layerMapIdx,temp1 = (6*(startSym>0)+equaSym)*pNrPuschOutParams->total_nPRB*12*nLyr + completedPrb*12*nLyr;
                    for(idxSc = 0; idxSc < 12; idxSc++)
                    {
                        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                        {
                            Y_vect_re[idxAnt][0] = VECT_SET(rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[7] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[6] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[5] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[4] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[3] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[2] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[1] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[0] + idxSc] );

                            Y_vect_im[idxAnt][0] = VECT_SET(rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[7] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[6] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[5] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[4] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[3] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[2] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[1] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[0] + idxSc] );
                        }

                        // y_tilda = W * Y
                        MATRIXMULH2(weights_re, weights_im, Y_vect_re, Y_vect_im, Y_tilda_re, Y_tilda_im, nLyr, nRxAnt, 1);

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

                            // Y_cap = Y_tilda * TOC
                            C_VECT_MUL_NA_NB(Y_tilda_re[idxLyr][0], Y_tilda_im[idxLyr][0], sq_toc_re[lyrType][idxSc], sq_toc_im[lyrType][idxSc], Y_cap_re, Y_cap_im);

                            // Store the Equalized output into memory
                            for(idx = 0; idx < 8; idx++)
                            {
                                layerMapIdx = temp1 + (prbScStartIdx_re[idx]-startPrbIdx*12+idxSc)*nLyr + idxLyr;
                                // pNrPuschOutParams->equaOutSamples_re[idxLyr][equaSym][prbScStartIdx_re[idx] + idxSc] = Y_cap_re[idx];
                                // pNrPuschOutParams->equaOutSamples_im[idxLyr][equaSym][prbScStartIdx_re[idx] + idxSc] = Y_cap_im[idx];
                                pNrPuschOutParams->layerDemapperOutReal[layerMapIdx] = Y_cap_re[idx];
                                pNrPuschOutParams->layerDemapperOutImag[layerMapIdx] = Y_cap_im[idx];
                            }
                        }
                    }
                    equaSym = equaSym + 1;
                }
            }

        }

        for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
        {
            mseAvgVec = _mm256_setzero_ps();

            for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb+=8)
            {
                mseVec = _mm256_loadu_ps(&pNrPuschOutParams->msePerPrb[idxLyr][idxPrb]);

                mseAvgVec = _mm256_add_ps(mseAvgVec, mseVec);
            }

            mid_re = _mm256_permute2f128_ps (mseAvgVec, mseAvgVec, 1);
            mseAvgVec = _mm256_add_ps(mseAvgVec, mid_re);
            mseAvgVec = _mm256_hadd_ps(mseAvgVec, mseAvgVec);
            mseAvgVec = _mm256_hadd_ps(mseAvgVec, mseAvgVec);

            mseAvg = mseAvgVec[0];
            mseAvg_inv = 1 / mseAvg;

            mseAvgInvVec = _mm256_set1_ps(mseAvg_inv);

            for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb+=8)
            {
                mseVec = _mm256_loadu_ps(&pNrPuschOutParams->msePerPrb[idxLyr][idxPrb]);

                mseVec = _mm256_mul_ps(mseVec, mseAvgInvVec);

                _mm256_storeu_ps(&pNrPuschOutParams->msePerPrb[idxLyr][idxPrb], mseVec);
            }
        }
    }
}


void nr_pusch_mmse_irc_perTone_est_avx2( int symPart,
                                         commonUlConfigT* pNrUlCommonParams,
                                         puschConfigT* pNrPuschInParams,
                                         P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams)
{
    int idxAnt, idxLyr, idxAlloc, idxPrb, idxSc, idxRow, idxCol, itr, idx;
    int startPrbIdx, endPrbIdx, startScIdx, endScIdx;
    int ctPrbIdx, ctScIdx, perToneScIdx, avxPrbIdx, ctPrgIdx;
    int printPrbIdx = 0, avxPrintIdx = 0;

    int prgScStartIdx, prbScStartIdx[8], perToneScStartIdx[8];

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
    __m256 chEst_perTone_re, chEst_perTone_im;

    __m256 bd_re, ad_im;
    __m256 mid_re, mid_im;
    __m256 mid2_re, mid2_im;
    __m256 mid3_re, mid3_im;
    __m256 out_re, out_im;

    __m256 scale_lsEst = _mm256_set1_ps(0.5012);

    __m256 div2Vec = _mm256_set1_ps(0.5);
    __m256 div4Vec = _mm256_set1_ps(0.25);
    __m256 div5Vec = _mm256_set1_ps(0.2);
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
    __m256 Rnn_inv_re[nRxAnt][nRxAnt], Rnn_inv_im[nRxAnt][nRxAnt];
    __m256 RnnTemp_re[nRxAnt][nRxAnt], RnnTemp_im[nRxAnt][nRxAnt];

    complex h1Avg[nRxAnt][8], h2Avg[nRxAnt][8], hAvg;
    complex Rnn[nRxAnt][nRxAnt], Inf[nRxAnt];
    float prgRnn_re[8], prgRnn_im[8];

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

    g_negate = VECT_SET1(-0.0);

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
        for(idxPrb = startPrbIdx; idxPrb <= endPrbIdx; idxPrb += 8, ctPrgIdx++)
        {
            prgScStartIdx = idxPrb * 6;
            ctScIdx = prgScStartIdx;
            perToneScIdx = idxPrb * 12;
            for(idx = 0; idx < 8; idx++)
            {
                prbScStartIdx[idx] = ctScIdx;
                perToneScStartIdx[idx] = perToneScIdx;

                ctScIdx += 6;
                perToneScIdx += 12;
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

                if((nEvenTonesLyr == 2) || (nOddTonesLyr == 2))
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

            //LS Estimation & Time offset estimation (TOE) on Even Tones
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
                            bd_re = _mm256_mul_ps(lsEst_evenIm[idxAnt][idxSc+2], lsEst_evenIm[idxAnt][idxSc]);
                            ad_im = _mm256_mul_ps(lsEst_evenRe[idxAnt][idxSc+2], lsEst_evenIm[idxAnt][idxSc]);

                            out_re = _mm256_fmadd_ps(lsEst_evenRe[idxAnt][idxSc+2], lsEst_evenRe[idxAnt][idxSc], bd_re);
                            out_im = _mm256_fmsub_ps(lsEst_evenIm[idxAnt][idxSc+2], lsEst_evenRe[idxAnt][idxSc], ad_im);

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
            }

            //LS Estimation & Time offset estimation (TOE) on Odd Tones
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
                            bd_re = _mm256_mul_ps(lsEst_oddIm[idxAnt][idxSc+2], lsEst_oddIm[idxAnt][idxSc]);
                            ad_im = _mm256_mul_ps(lsEst_oddRe[idxAnt][idxSc+2], lsEst_oddIm[idxAnt][idxSc]);

                            out_re = _mm256_fmadd_ps(lsEst_oddRe[idxAnt][idxSc+2], lsEst_oddRe[idxAnt][idxSc], bd_re);
                            out_im = _mm256_fmsub_ps(lsEst_oddIm[idxAnt][idxSc+2], lsEst_oddRe[idxAnt][idxSc], ad_im);

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
            }

            //Over all Time offset estimation (TOE)
            if((nEvenTonesLyr > 0) && (nOddTonesLyr > 0))
            {
                toe_overAll = ((toe_evenTones) + (toe_oddTones)) * 0.5;
//                toe_overAll = toe_oddTones;
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

            if((nEvenTonesLyr == 2) || (nOddTonesLyr == 2))
            {
                toe_overAll = csqrt(toe_overAll);
            }

            pNrPuschOutParams->toc_overAll_re[ctPrgIdx] = creal(toe_overAll);
            pNrPuschOutParams->toc_overAll_im[ctPrgIdx] = cimag(toe_overAll);

            //Channel Estimation & Covariance Matrix (Rnn) computation on Even Tone
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

                //Per Tone Channel Estimation
                if(nEvenTonesLyr == 1)
                {
                    for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                    {
                        for(idxSc = 0; idxSc < 6; idxSc++)
                        {
                            for(idx = 0; idx < 8; idx++)
                            {
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][0][perToneScStartIdx[idx] + (2*idxSc)] = lsEst_evenRe[idxAnt][idxSc][idx];
                                pNrPuschOutParams->chEstPerTone_im[idxAnt][0][perToneScStartIdx[idx] + (2*idxSc)] = lsEst_evenIm[idxAnt][idxSc][idx];
                            }

                            if(idxSc == 5)
                            {
                                for(idx = 0; idx < 7; idx++)
                                {
                                    pNrPuschOutParams->chEstPerTone_re[idxAnt][0][perToneScStartIdx[idx] + (2*idxSc+1)] = 0.5 * (lsEst_evenRe[idxAnt][idxSc][idx] + lsEst_evenRe[idxAnt][0][idx+1]);
                                    pNrPuschOutParams->chEstPerTone_im[idxAnt][0][perToneScStartIdx[idx] + (2*idxSc+1)] = 0.5 * (lsEst_evenIm[idxAnt][idxSc][idx] + lsEst_evenIm[idxAnt][0][idx+1]);
                                }

                                pNrPuschOutParams->chEstPerTone_re[idxAnt][0][perToneScStartIdx[7] + (2*idxSc+1)] = lsEst_evenRe[idxAnt][idxSc][idx];
                                pNrPuschOutParams->chEstPerTone_im[idxAnt][0][perToneScStartIdx[7] + (2*idxSc+1)] = lsEst_evenIm[idxAnt][idxSc][idx];
                            }
                            else
                            {
                                chEst_perTone_re = _mm256_add_ps(lsEst_evenRe[idxAnt][idxSc], lsEst_evenRe[idxAnt][idxSc+1]);
                                chEst_perTone_im = _mm256_add_ps(lsEst_evenIm[idxAnt][idxSc], lsEst_evenIm[idxAnt][idxSc+1]);

                                chEst_perTone_re = _mm256_mul_ps(chEst_perTone_re, div2Vec);
                                chEst_perTone_im = _mm256_mul_ps(chEst_perTone_im, div2Vec);

                                for(idx = 0; idx < 8; idx++)
                                {
                                    pNrPuschOutParams->chEstPerTone_re[idxAnt][0][perToneScStartIdx[idx] + (2*idxSc+1)] = chEst_perTone_re[idx];
                                    pNrPuschOutParams->chEstPerTone_im[idxAnt][0][perToneScStartIdx[idx] + (2*idxSc+1)] = chEst_perTone_im[idx];
                                }
                            }
                        }
                    }
                }
                else
                {
                    for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                    {
                        for(idxSc = 0; idxSc < 6; idxSc += 2)
                        {
                            chEst_perTone_re = _mm256_add_ps(lsEst_evenRe[idxAnt][idxSc], lsEst_evenRe[idxAnt][idxSc+1]);
                            chEst_perTone_im = _mm256_add_ps(lsEst_evenIm[idxAnt][idxSc], lsEst_evenIm[idxAnt][idxSc+1]);

                            chEst_perTone_re = _mm256_mul_ps(chEst_perTone_re, div2Vec);
                            chEst_perTone_im = _mm256_mul_ps(chEst_perTone_im, div2Vec);

                            for(idx = 0; idx < 8; idx++)
                            {
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][0][perToneScStartIdx[idx] + (2*idxSc  )] = chEst_perTone_re[idx];
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][0][perToneScStartIdx[idx] + (2*idxSc+1)] = chEst_perTone_re[idx];
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][0][perToneScStartIdx[idx] + (2*idxSc+2)] = chEst_perTone_re[idx];
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][0][perToneScStartIdx[idx] + (2*idxSc+3)] = chEst_perTone_re[idx];

                                pNrPuschOutParams->chEstPerTone_im[idxAnt][0][perToneScStartIdx[idx] + (2*idxSc  )] = chEst_perTone_im[idx];
                                pNrPuschOutParams->chEstPerTone_im[idxAnt][0][perToneScStartIdx[idx] + (2*idxSc+1)] = chEst_perTone_im[idx];
                                pNrPuschOutParams->chEstPerTone_im[idxAnt][0][perToneScStartIdx[idx] + (2*idxSc+2)] = chEst_perTone_im[idx];
                                pNrPuschOutParams->chEstPerTone_im[idxAnt][0][perToneScStartIdx[idx] + (2*idxSc+3)] = chEst_perTone_im[idx];
                            }

                            chEst_perTone_re = _mm256_sub_ps(lsEst_evenRe[idxAnt][idxSc], lsEst_evenRe[idxAnt][idxSc+1]);
                            chEst_perTone_im = _mm256_sub_ps(lsEst_evenIm[idxAnt][idxSc], lsEst_evenIm[idxAnt][idxSc+1]);

                            chEst_perTone_re = _mm256_mul_ps(chEst_perTone_re, div2Vec);
                            chEst_perTone_im = _mm256_mul_ps(chEst_perTone_im, div2Vec);

                            for(idx = 0; idx < 8; idx++)
                            {
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][1][perToneScStartIdx[idx] + (2*idxSc  )] = chEst_perTone_re[idx];
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][1][perToneScStartIdx[idx] + (2*idxSc+1)] = chEst_perTone_re[idx];
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][1][perToneScStartIdx[idx] + (2*idxSc+2)] = chEst_perTone_re[idx];
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][1][perToneScStartIdx[idx] + (2*idxSc+3)] = chEst_perTone_re[idx];

                                pNrPuschOutParams->chEstPerTone_im[idxAnt][1][perToneScStartIdx[idx] + (2*idxSc  )] = chEst_perTone_im[idx];
                                pNrPuschOutParams->chEstPerTone_im[idxAnt][1][perToneScStartIdx[idx] + (2*idxSc+1)] = chEst_perTone_im[idx];
                                pNrPuschOutParams->chEstPerTone_im[idxAnt][1][perToneScStartIdx[idx] + (2*idxSc+2)] = chEst_perTone_im[idx];
                                pNrPuschOutParams->chEstPerTone_im[idxAnt][1][perToneScStartIdx[idx] + (2*idxSc+3)] = chEst_perTone_im[idx];
                            }
                        }
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

            //Channel Estimation & Covariance Matrix (Rnn) computation on Odd Tone
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

                //Per Tone Channel Estimation
                if(nOddTonesLyr == 1)
                {
                    for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                    {
                        for(idxSc = 0; idxSc < 6; idxSc++)
                        {
                            for(idx = 0; idx < 8; idx++)
                            {
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][2][perToneScStartIdx[idx] + (2*idxSc+1)] = lsEst_oddRe[idxAnt][idxSc][idx];
                                pNrPuschOutParams->chEstPerTone_im[idxAnt][2][perToneScStartIdx[idx] + (2*idxSc+1)] = lsEst_oddIm[idxAnt][idxSc][idx];
                            }

                            chEst_perTone_re = _mm256_add_ps(lsEst_oddRe[idxAnt][idxSc], lsEst_oddRe[idxAnt][idxSc+1]);
                            chEst_perTone_im = _mm256_add_ps(lsEst_oddIm[idxAnt][idxSc], lsEst_oddIm[idxAnt][idxSc+1]);

                            chEst_perTone_re = _mm256_mul_ps(chEst_perTone_re, div2Vec);
                            chEst_perTone_im = _mm256_mul_ps(chEst_perTone_im, div2Vec);

                            for(idx = 0; idx < 8; idx++)
                            {
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][2][perToneScStartIdx[idx] + (2*idxSc+2)] = chEst_perTone_re[idx];
                                pNrPuschOutParams->chEstPerTone_im[idxAnt][2][perToneScStartIdx[idx] + (2*idxSc+2)] = chEst_perTone_im[idx];
                            }
                        }

                        //Adjusting 1st tone
                        pNrPuschOutParams->chEstPerTone_re[idxAnt][2][perToneScStartIdx[0]] = lsEst_oddRe[idxAnt][0][0];
                        pNrPuschOutParams->chEstPerTone_im[idxAnt][2][perToneScStartIdx[0]] = lsEst_oddIm[idxAnt][0][0];

                        for(idx = 1; idx < 8; idx++)
                        {
                            pNrPuschOutParams->chEstPerTone_re[idxAnt][2][perToneScStartIdx[idx]] = 0.5 * (lsEst_oddRe[idxAnt][0][idx] + lsEst_oddRe[idxAnt][5][idx-1]);
                            pNrPuschOutParams->chEstPerTone_im[idxAnt][2][perToneScStartIdx[idx]] = 0.5 * (lsEst_oddIm[idxAnt][0][idx] + lsEst_oddIm[idxAnt][5][idx-1]);
                        }
                    }
                }
                else
                {
                    for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                    {
                        for(idxSc = 0; idxSc < 6; idxSc += 2)
                        {
                            chEst_perTone_re = _mm256_add_ps(lsEst_oddRe[idxAnt][idxSc], lsEst_oddRe[idxAnt][idxSc+1]);
                            chEst_perTone_im = _mm256_add_ps(lsEst_oddIm[idxAnt][idxSc], lsEst_oddIm[idxAnt][idxSc+1]);

                            chEst_perTone_re = _mm256_mul_ps(chEst_perTone_re, div2Vec);
                            chEst_perTone_im = _mm256_mul_ps(chEst_perTone_im, div2Vec);

                            for(idx = 0; idx < 8; idx++)
                            {
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][2][perToneScStartIdx[idx] + (2*idxSc  )] = chEst_perTone_re[idx];
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][2][perToneScStartIdx[idx] + (2*idxSc+1)] = chEst_perTone_re[idx];
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][2][perToneScStartIdx[idx] + (2*idxSc+2)] = chEst_perTone_re[idx];
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][2][perToneScStartIdx[idx] + (2*idxSc+3)] = chEst_perTone_re[idx];

                                pNrPuschOutParams->chEstPerTone_im[idxAnt][2][perToneScStartIdx[idx] + (2*idxSc  )] = chEst_perTone_im[idx];
                                pNrPuschOutParams->chEstPerTone_im[idxAnt][2][perToneScStartIdx[idx] + (2*idxSc+1)] = chEst_perTone_im[idx];
                                pNrPuschOutParams->chEstPerTone_im[idxAnt][2][perToneScStartIdx[idx] + (2*idxSc+2)] = chEst_perTone_im[idx];
                                pNrPuschOutParams->chEstPerTone_im[idxAnt][2][perToneScStartIdx[idx] + (2*idxSc+3)] = chEst_perTone_im[idx];
                            }

                            chEst_perTone_re = _mm256_sub_ps(lsEst_oddRe[idxAnt][idxSc], lsEst_oddRe[idxAnt][idxSc+1]);
                            chEst_perTone_im = _mm256_sub_ps(lsEst_oddIm[idxAnt][idxSc], lsEst_oddIm[idxAnt][idxSc+1]);

                            chEst_perTone_re = _mm256_mul_ps(chEst_perTone_re, div2Vec);
                            chEst_perTone_im = _mm256_mul_ps(chEst_perTone_im, div2Vec);

                            for(idx = 0; idx < 8; idx++)
                            {
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][3][perToneScStartIdx[idx] + (2*idxSc  )] = chEst_perTone_re[idx];
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][3][perToneScStartIdx[idx] + (2*idxSc+1)] = chEst_perTone_re[idx];
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][3][perToneScStartIdx[idx] + (2*idxSc+2)] = chEst_perTone_re[idx];
                                pNrPuschOutParams->chEstPerTone_re[idxAnt][3][perToneScStartIdx[idx] + (2*idxSc+3)] = chEst_perTone_re[idx];

                                pNrPuschOutParams->chEstPerTone_im[idxAnt][3][perToneScStartIdx[idx] + (2*idxSc  )] = chEst_perTone_im[idx];
                                pNrPuschOutParams->chEstPerTone_im[idxAnt][3][perToneScStartIdx[idx] + (2*idxSc+1)] = chEst_perTone_im[idx];
                                pNrPuschOutParams->chEstPerTone_im[idxAnt][3][perToneScStartIdx[idx] + (2*idxSc+2)] = chEst_perTone_im[idx];
                                pNrPuschOutParams->chEstPerTone_im[idxAnt][3][perToneScStartIdx[idx] + (2*idxSc+3)] = chEst_perTone_im[idx];
                            }
                        }
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

                        if(nOddTonesLyr == 2)
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

        if(COV_GRP_FACTOR == 1)
        {
            for(idxPrb = startPrbIdx; idxPrb <= endPrbIdx; idxPrb+=8)
            {
                for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                {
                    for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                    {
                        if((pNrPuschInParams->nullTone_enable == 1) || ((pNrPuschOutParams->nEvenTonesLyr > 0) && (pNrPuschOutParams->nOddTonesLyr > 0)))
                        {
                            prgRnn_re[0] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb] ) * 0.5;
                            prgRnn_re[1] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+1] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+1] ) * 0.5;
                            prgRnn_re[2] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+2] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+2] ) * 0.5;
                            prgRnn_re[3] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+3] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+3] ) * 0.5;
                            prgRnn_re[4] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+4] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+4] ) * 0.5;
                            prgRnn_re[5] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+5] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+5] ) * 0.5;
                            prgRnn_re[6] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+6] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+6] ) * 0.5;
                            prgRnn_re[7] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+7] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+7] ) * 0.5;

                            prgRnn_im[0] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb] ) * 0.5;
                            prgRnn_im[1] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+1] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+1] ) * 0.5;
                            prgRnn_im[2] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+2] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+2] ) * 0.5;
                            prgRnn_im[3] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+3] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+3] ) * 0.5;
                            prgRnn_im[4] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+4] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+4] ) * 0.5;
                            prgRnn_im[5] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+5] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+5] ) * 0.5;
                            prgRnn_im[6] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+6] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+6] ) * 0.5;
                            prgRnn_im[7] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+7] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+7] ) * 0.5;

                            Rnn_re[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_re);
                            Rnn_im[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_im);
                        }
                        else if (pNrPuschOutParams->nEvenTonesLyr > 0)
                        {
                            prgRnn_re[0] = pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb];
                            prgRnn_re[1] = pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+1];
                            prgRnn_re[2] = pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+2];
                            prgRnn_re[3] = pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+3];
                            prgRnn_re[4] = pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+4];
                            prgRnn_re[5] = pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+5];
                            prgRnn_re[6] = pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+6];
                            prgRnn_re[7] = pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+7];

                            prgRnn_im[0] = pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb];
                            prgRnn_im[1] = pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+1];
                            prgRnn_im[2] = pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+2];
                            prgRnn_im[3] = pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+3];
                            prgRnn_im[4] = pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+4];
                            prgRnn_im[5] = pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+5];
                            prgRnn_im[6] = pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+6];
                            prgRnn_im[7] = pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+7];

                            Rnn_re[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_re);
                            Rnn_im[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_im);

                        }
                        else
                        {
                            prgRnn_re[0] = pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb];
                            prgRnn_re[1] = pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+1];
                            prgRnn_re[2] = pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+2];
                            prgRnn_re[3] = pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+3];
                            prgRnn_re[4] = pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+4];
                            prgRnn_re[5] = pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+5];
                            prgRnn_re[6] = pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+6];
                            prgRnn_re[7] = pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+7];

                            prgRnn_im[0] = pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb];
                            prgRnn_im[1] = pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+1];
                            prgRnn_im[2] = pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+2];
                            prgRnn_im[3] = pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+3];
                            prgRnn_im[4] = pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+4];
                            prgRnn_im[5] = pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+5];
                            prgRnn_im[6] = pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+6];
                            prgRnn_im[7] = pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+7];

                            Rnn_re[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_re);
                            Rnn_im[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_im);

                        }
                    }
                }

                if(nRxAnt == 1)
                {
                    InverseMatrix1x1HTranspose(Rnn_re, Rnn_im, Rnn_inv_re, Rnn_inv_im);
                }
                else if(nRxAnt == 2)
                {
                    InverseMatrix2x2H(Rnn_re, Rnn_im, Rnn_inv_re, Rnn_inv_im);
                }
                else if(nRxAnt == 4)
                {
                    InverseMatrix4x4H(Rnn_re, Rnn_im, Rnn_inv_re, Rnn_inv_im);
                }
                else if(nRxAnt == 8)
                {
                    InverseMatrix8x8H(Rnn_re, Rnn_im, Rnn_inv_re, Rnn_inv_im);
                }

                for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                {
                    for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                    {
                        _mm256_storeu_ps(pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol], Rnn_inv_re[idxRow][idxCol]);
                        _mm256_storeu_ps(pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol], Rnn_inv_im[idxRow][idxCol]);
                    }
                }
            }
        }
        else if(COV_GRP_FACTOR == 2)
        {
            for(idxPrb = startPrbIdx; idxPrb <= endPrbIdx; idxPrb+=16)
            {
                for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                {
                    for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                    {
                        if((pNrPuschInParams->nullTone_enable == 1) || ((pNrPuschOutParams->nEvenTonesLyr > 0) && (pNrPuschOutParams->nOddTonesLyr > 0)))
                        {
                            prgRnn_re[0] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb] +
                                             pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+1] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+1] ) * 0.25;
                            prgRnn_re[1] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+2] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+2] +
                                             pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+3] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+3] ) * 0.25;
                            prgRnn_re[2] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+4] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+4] +
                                             pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+5] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+5] ) * 0.25;
                            prgRnn_re[3] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+6] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+6] +
                                             pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+7] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+7] ) * 0.25;
                            prgRnn_re[4] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+8] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+8] +
                                             pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+9] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+9] ) * 0.25;
                            prgRnn_re[5] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+10] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+10] +
                                             pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+11] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+11] ) * 0.25;
                            prgRnn_re[6] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+12] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+12] +
                                             pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+13] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+13] ) * 0.25;
                            prgRnn_re[7] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+14] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+14] +
                                             pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+15] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+15] ) * 0.25;

                            prgRnn_im[0] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb] +
                                             pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+1] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+1] ) * 0.25;
                            prgRnn_im[1] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+2] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+2] +
                                             pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+3] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+3] ) * 0.25;
                            prgRnn_im[2] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+4] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+4] +
                                             pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+5] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+5] ) * 0.25;
                            prgRnn_im[3] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+6] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+6] +
                                             pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+7] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+7] ) * 0.25;
                            prgRnn_im[4] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+8] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+8] +
                                             pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+9] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+9] ) * 0.25;
                            prgRnn_im[5] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+10] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+10] +
                                             pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+11] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+11] ) * 0.25;
                            prgRnn_im[6] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+12] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+12] +
                                             pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+13] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+13] ) * 0.25;
                            prgRnn_im[7] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+14] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+14] +
                                             pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+15] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+15] ) * 0.25;

                            Rnn_re[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_re);
                            Rnn_im[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_im);
                        }
                        else if (pNrPuschOutParams->nEvenTonesLyr > 0)
                        {
                            prgRnn_re[0] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb] + pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+1] ) * 0.5;
                            prgRnn_re[1] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+2] + pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+3] ) * 0.5;
                            prgRnn_re[2] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+4] + pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+5] ) * 0.5;
                            prgRnn_re[3] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+6] + pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+7] ) * 0.5;
                            prgRnn_re[4] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+8] + pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+9] ) * 0.5;
                            prgRnn_re[5] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+10] + pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+11] ) * 0.5;
                            prgRnn_re[6] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+12] + pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+13] ) * 0.5;
                            prgRnn_re[7] = ( pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+14] + pNrPuschOutParams->RnnEven_prb_re[idxRow][idxCol][idxPrb+15] ) * 0.5;

                            prgRnn_im[0] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb] + pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+1] ) * 0.5;
                            prgRnn_im[1] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+2] + pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+3] ) * 0.5;
                            prgRnn_im[2] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+4] + pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+5] ) * 0.5;
                            prgRnn_im[3] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+6] + pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+7] ) * 0.5;
                            prgRnn_im[4] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+8] + pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+9] ) * 0.5;
                            prgRnn_im[5] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+10] + pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+11] ) * 0.5;
                            prgRnn_im[6] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+12] + pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+13] ) * 0.5;
                            prgRnn_im[7] = ( pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+14] + pNrPuschOutParams->RnnEven_prb_im[idxRow][idxCol][idxPrb+15] ) * 0.5;

                            Rnn_re[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_re);
                            Rnn_im[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_im);

                        }
                        else
                        {
                            prgRnn_re[0] = ( pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+1] ) * 0.5;
                            prgRnn_re[1] = ( pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+2] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+3] ) * 0.5;
                            prgRnn_re[2] = ( pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+4] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+5] ) * 0.5;
                            prgRnn_re[3] = ( pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+6] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+7] ) * 0.5;
                            prgRnn_re[4] = ( pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+8] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+9] ) * 0.5;
                            prgRnn_re[5] = ( pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+10] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+11] ) * 0.5;
                            prgRnn_re[6] = ( pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+12] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+13] ) * 0.5;
                            prgRnn_re[7] = ( pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+14] + pNrPuschOutParams->RnnOdd_prb_re[idxRow][idxCol][idxPrb+15] ) * 0.5;

                            prgRnn_im[0] = ( pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+1] ) * 0.5;
                            prgRnn_im[1] = ( pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+2] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+3] ) * 0.5;
                            prgRnn_im[2] = ( pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+4] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+5] ) * 0.5;
                            prgRnn_im[3] = ( pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+6] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+7] ) * 0.5;
                            prgRnn_im[4] = ( pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+8] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+9] ) * 0.5;
                            prgRnn_im[5] = ( pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+10] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+11] ) * 0.5;
                            prgRnn_im[6] = ( pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+12] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+13] ) * 0.5;
                            prgRnn_im[7] = ( pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+14] + pNrPuschOutParams->RnnOdd_prb_im[idxRow][idxCol][idxPrb+15] ) * 0.5;

                            Rnn_re[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_re);
                            Rnn_im[idxRow][idxCol] = _mm256_loadu_ps(prgRnn_im);

                        }

                    }
                }

                if(nRxAnt == 1)
                {
                    InverseMatrix1x1HTranspose(Rnn_re, Rnn_im, Rnn_inv_re, Rnn_inv_im);
                }
                else if(nRxAnt == 2)
                {
                    InverseMatrix2x2H(Rnn_re, Rnn_im, Rnn_inv_re, Rnn_inv_im);
                }
                else if(nRxAnt == 4)
                {
                    InverseMatrix4x4H(Rnn_re, Rnn_im, Rnn_inv_re, Rnn_inv_im);
                }
                else if(nRxAnt == 8)
                {
                    InverseMatrix8x8H(Rnn_re, Rnn_im, Rnn_inv_re, Rnn_inv_im);
                }

                for(idxRow = 0; idxRow < nRxAnt; idxRow++)
                {
                    for(idxCol = 0; idxCol < nRxAnt; idxCol++)
                    {
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb] = Rnn_inv_re[idxRow][idxCol][0];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+1] = Rnn_inv_re[idxRow][idxCol][0];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+2] = Rnn_inv_re[idxRow][idxCol][1];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+3] = Rnn_inv_re[idxRow][idxCol][1];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+4] = Rnn_inv_re[idxRow][idxCol][2];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+5] = Rnn_inv_re[idxRow][idxCol][2];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+6] = Rnn_inv_re[idxRow][idxCol][3];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+7] = Rnn_inv_re[idxRow][idxCol][3];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+8] = Rnn_inv_re[idxRow][idxCol][4];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+9] = Rnn_inv_re[idxRow][idxCol][4];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+10] = Rnn_inv_re[idxRow][idxCol][5];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+11] = Rnn_inv_re[idxRow][idxCol][5];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+12] = Rnn_inv_re[idxRow][idxCol][6];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+13] = Rnn_inv_re[idxRow][idxCol][6];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+14] = Rnn_inv_re[idxRow][idxCol][7];
                        pNrPuschOutParams->RnnInv_prb_re[idxRow][idxCol][idxPrb+15] = Rnn_inv_re[idxRow][idxCol][7];

                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb] = Rnn_inv_im[idxRow][idxCol][0];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+1] = Rnn_inv_im[idxRow][idxCol][0];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+2] = Rnn_inv_im[idxRow][idxCol][1];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+3] = Rnn_inv_im[idxRow][idxCol][1];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+4] = Rnn_inv_im[idxRow][idxCol][2];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+5] = Rnn_inv_im[idxRow][idxCol][2];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+6] = Rnn_inv_im[idxRow][idxCol][3];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+7] = Rnn_inv_im[idxRow][idxCol][3];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+8] = Rnn_inv_im[idxRow][idxCol][4];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+9] = Rnn_inv_im[idxRow][idxCol][4];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+10] = Rnn_inv_im[idxRow][idxCol][5];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+11] = Rnn_inv_im[idxRow][idxCol][5];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+12] = Rnn_inv_im[idxRow][idxCol][6];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+13] = Rnn_inv_im[idxRow][idxCol][6];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+14] = Rnn_inv_im[idxRow][idxCol][7];
                        pNrPuschOutParams->RnnInv_prb_im[idxRow][idxCol][idxPrb+15] = Rnn_inv_im[idxRow][idxCol][7];
                    }
                }
            }
        }
    }
}

void nr_pusch_mmse_irc_perTone_equ_avx2( int symPart,
                                         commonUlConfigT* pNrUlCommonParams,
                                         puschConfigT* pNrPuschInParams,
                                         P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams)
{
    int idxAnt, idxLyr, idxAlloc, idxPrb, idxSym, idxSc, rIndx, cIndx, i, idx;
    int startScIdx, endScIdx, remSc;
    int startPrbIdx, endPrbIdx, CtPrbScStartIdx, avxPrbIdx, ctPrgIdx, lyrPort, lyrType;
    int printPrbIdx = 0, avxPrintIdx = 0;

    int startSym, endSym, dmrsSym;
    int nRxAnt, nLyr, nSbAlloc, nPrb_alloc, nPrb_total;

    int prbScStartIdx[8], prbScStartIdx_re[8], prbScStartIdx_im[8];

    nRxAnt   = pNrPuschInParams->nAntennaPorts;
    nLyr     = pNrPuschInParams->nNrOfLayers;
    nSbAlloc = pNrPuschInParams->nSbAlloc;

    float Rnn_grp1_re, Rnn_grp1_im;
    float Rnn_grp2_re, Rnn_grp2_im;
    float Rnn_grp3_re, Rnn_grp3_im;
    float Rnn_grp4_re, Rnn_grp4_im;

    __m256 mid_re, mid_im;

    __m256 hMat_re[nRxAnt][nLyr], hMat_im[nRxAnt][nLyr];
    __m256 hHerMat_re[nLyr][nRxAnt], hHerMat_im[nLyr][nRxAnt];
    __m256 Rnn_re[nRxAnt][nRxAnt], Rnn_im[nRxAnt][nRxAnt];
    __m256 Rnn_inv_re[nRxAnt][nRxAnt], Rnn_inv_im[nRxAnt][nRxAnt];
    __m256 hHer_RnnInv_re[nLyr][nRxAnt], hHer_RnnInv_im[nLyr][nRxAnt];
    __m256 hHer_RnnInv_h_re[nLyr][nLyr], hHer_RnnInv_h_im[nLyr][nLyr];
    __m256 hHer_RnnInv_h_inv_re[nLyr][nLyr], hHer_RnnInv_h_inv_im[nLyr][nLyr];

    __m256 hHerRnnInvhInv_hHer_re[nLyr][nRxAnt], hHerRnnInvhInv_hHer_im[nLyr][nRxAnt];
    __m256 weights_re[nLyr][nRxAnt], weights_im[nLyr][nRxAnt];
    __m256 Y_vect_re[nRxAnt][1], Y_vect_im[nRxAnt][1];
    __m256 Y_tilda_re[nLyr][1], Y_tilda_im[nLyr][1];
    __m256 Y_cap_re, Y_cap_im;
    __m256 sq_toc_re[2][12];
	__m256 sq_toc_im[2][12];

	__m256 mseVec, mseAvgVec, mseAvgInvVec;

	g_negate = VECT_SET1(-0.0);

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

    float complex toc_overAll, sqrt_toc[12];
    float mse, mse_inv;
    float mseAvg, mseAvg_inv;
    int completedPrb = 0;
    nPrb_total = 0;
    for(idxAlloc = 0; idxAlloc < nSbAlloc; idxAlloc++)
    {
        nPrb_alloc = pNrPuschInParams->nPrb[idxAlloc];
        completedPrb = nPrb_total;
        nPrb_total += nPrb_alloc;

        startPrbIdx = pNrPuschInParams->rbStart[idxAlloc];
        endPrbIdx   = startPrbIdx + (pNrPuschInParams->nPrb[idxAlloc]) - 1;
        CtPrbScStartIdx = startPrbIdx * 12;

        for(idxPrb = startPrbIdx; idxPrb < (endPrbIdx); idxPrb += 8, CtPrbScStartIdx += 12)
        {
            ctPrgIdx = idxPrb / 8;

            prbScStartIdx[0] = idxPrb * 12;

            prbScStartIdx_re[0] = prbScStartIdx[0];
            prbScStartIdx_im[0] = MAX_NUM_SC + prbScStartIdx[0];

            for(idx = 1; idx < 8; idx++)
            {
                prbScStartIdx[idx] = prbScStartIdx[idx-1] + 12;

                prbScStartIdx_re[idx] = prbScStartIdx_re[idx-1] + 12;
                prbScStartIdx_im[idx] = prbScStartIdx_im[idx-1] + 12;
            }

            // Load Rnn_inv Matrix
            for(rIndx = 0; rIndx < nRxAnt; rIndx++)
            {
                for(cIndx = 0; cIndx < nRxAnt; cIndx++)
                {
                    Rnn_inv_re[rIndx][cIndx] = _mm256_loadu_ps(&pNrPuschOutParams->RnnInv_prb_re[rIndx][cIndx][idxPrb]);
                    Rnn_inv_im[rIndx][cIndx] = _mm256_loadu_ps(&pNrPuschOutParams->RnnInv_prb_im[rIndx][cIndx][idxPrb]);
                }
            }

            //TOC Vector of current PRG for Even Layers
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

            //TOC Vector of current PRG for Odd Layers
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

            for(idxSc = 0; idxSc < 12; idxSc++)
            {
                // Load H and H_her Matrix
                for(rIndx = 0; rIndx < nRxAnt; rIndx++)
                {
                    for(cIndx = 0; cIndx < nLyr; cIndx++)
                    {
                        lyrPort = pNrPuschInParams->nPortIndex[cIndx];

                        hMat_re[rIndx][cIndx] = VECT_SET(pNrPuschOutParams->chEstPerTone_re[rIndx][lyrPort][prbScStartIdx[7] + idxSc],
                                                         pNrPuschOutParams->chEstPerTone_re[rIndx][lyrPort][prbScStartIdx[6] + idxSc],
                                                         pNrPuschOutParams->chEstPerTone_re[rIndx][lyrPort][prbScStartIdx[5] + idxSc],
                                                         pNrPuschOutParams->chEstPerTone_re[rIndx][lyrPort][prbScStartIdx[4] + idxSc],
                                                         pNrPuschOutParams->chEstPerTone_re[rIndx][lyrPort][prbScStartIdx[3] + idxSc],
                                                         pNrPuschOutParams->chEstPerTone_re[rIndx][lyrPort][prbScStartIdx[2] + idxSc],
                                                         pNrPuschOutParams->chEstPerTone_re[rIndx][lyrPort][prbScStartIdx[1] + idxSc],
                                                         pNrPuschOutParams->chEstPerTone_re[rIndx][lyrPort][prbScStartIdx[0] + idxSc] );

                        hMat_im[rIndx][cIndx] = VECT_SET(pNrPuschOutParams->chEstPerTone_im[rIndx][lyrPort][prbScStartIdx[7] + idxSc],
                                                         pNrPuschOutParams->chEstPerTone_im[rIndx][lyrPort][prbScStartIdx[6] + idxSc],
                                                         pNrPuschOutParams->chEstPerTone_im[rIndx][lyrPort][prbScStartIdx[5] + idxSc],
                                                         pNrPuschOutParams->chEstPerTone_im[rIndx][lyrPort][prbScStartIdx[4] + idxSc],
                                                         pNrPuschOutParams->chEstPerTone_im[rIndx][lyrPort][prbScStartIdx[3] + idxSc],
                                                         pNrPuschOutParams->chEstPerTone_im[rIndx][lyrPort][prbScStartIdx[2] + idxSc],
                                                         pNrPuschOutParams->chEstPerTone_im[rIndx][lyrPort][prbScStartIdx[1] + idxSc],
                                                         pNrPuschOutParams->chEstPerTone_im[rIndx][lyrPort][prbScStartIdx[0] + idxSc] );

                        hHerMat_re[cIndx][rIndx] = VECT_SET(pNrPuschOutParams->chEstPerTone_re[rIndx][lyrPort][prbScStartIdx[7] + idxSc],
                                                            pNrPuschOutParams->chEstPerTone_re[rIndx][lyrPort][prbScStartIdx[6] + idxSc],
                                                            pNrPuschOutParams->chEstPerTone_re[rIndx][lyrPort][prbScStartIdx[5] + idxSc],
                                                            pNrPuschOutParams->chEstPerTone_re[rIndx][lyrPort][prbScStartIdx[4] + idxSc],
                                                            pNrPuschOutParams->chEstPerTone_re[rIndx][lyrPort][prbScStartIdx[3] + idxSc],
                                                            pNrPuschOutParams->chEstPerTone_re[rIndx][lyrPort][prbScStartIdx[2] + idxSc],
                                                            pNrPuschOutParams->chEstPerTone_re[rIndx][lyrPort][prbScStartIdx[1] + idxSc],
                                                            pNrPuschOutParams->chEstPerTone_re[rIndx][lyrPort][prbScStartIdx[0] + idxSc] );

                        hHerMat_im[cIndx][rIndx] = VECT_SET(-pNrPuschOutParams->chEstPerTone_im[rIndx][lyrPort][prbScStartIdx[7] + idxSc],
                                                            -pNrPuschOutParams->chEstPerTone_im[rIndx][lyrPort][prbScStartIdx[6] + idxSc],
                                                            -pNrPuschOutParams->chEstPerTone_im[rIndx][lyrPort][prbScStartIdx[5] + idxSc],
                                                            -pNrPuschOutParams->chEstPerTone_im[rIndx][lyrPort][prbScStartIdx[4] + idxSc],
                                                            -pNrPuschOutParams->chEstPerTone_im[rIndx][lyrPort][prbScStartIdx[3] + idxSc],
                                                            -pNrPuschOutParams->chEstPerTone_im[rIndx][lyrPort][prbScStartIdx[2] + idxSc],
                                                            -pNrPuschOutParams->chEstPerTone_im[rIndx][lyrPort][prbScStartIdx[1] + idxSc],
                                                            -pNrPuschOutParams->chEstPerTone_im[rIndx][lyrPort][prbScStartIdx[0] + idxSc] );

                    }
                }

                // Perform [H_her * Rnn_inv]
                MATRIXMULH2(hHerMat_re, hHerMat_im, Rnn_inv_re, Rnn_inv_im, hHer_RnnInv_re, hHer_RnnInv_im, nLyr, nRxAnt, nRxAnt);

                // Perform [(H_her * Rnn_inv) * H]
                MATRIXMULH2(hHer_RnnInv_re, hHer_RnnInv_im, hMat_re, hMat_im, hHer_RnnInv_h_re, hHer_RnnInv_h_im, nLyr, nRxAnt, nLyr);

                // Store MSE for LLR Calculations
                if(nLyr == 1)
                {
                    for(idx = 0; idx < 8; idx++)
                    {
                        pNrPuschOutParams->msePerPrb[0][idxPrb + idx] = hHer_RnnInv_h_re[0][0][idx];
                    }
                }

                // Perform [(H_her * Rnn_inv * H) + I]
                for(avxPrbIdx = 0; avxPrbIdx < 8; avxPrbIdx++)
                {
                    for(rIndx = 0; rIndx < nLyr; rIndx++)
                    {
                        hHer_RnnInv_h_re[rIndx][rIndx][avxPrbIdx] += 1;
                    }
                }

                // Perform [inv((H_her * Rnn_inv * H) + I)]
                if(nLyr == 1)
                {
                    InverseMatrix1x1HTranspose(hHer_RnnInv_h_re, hHer_RnnInv_h_im, hHer_RnnInv_h_inv_re, hHer_RnnInv_h_inv_im);
                }
                else if(nLyr == 2)
                {
                    InverseMatrix2x2H(hHer_RnnInv_h_re, hHer_RnnInv_h_im, hHer_RnnInv_h_inv_re, hHer_RnnInv_h_inv_im);
                }
                else if(nLyr == 4)
                {
                    InverseMatrix4x4H(hHer_RnnInv_h_re, hHer_RnnInv_h_im, hHer_RnnInv_h_inv_re, hHer_RnnInv_h_inv_im);
                }

                // Store MSE for LLR Calculations
                if(nLyr > 1)
                {
                    for(rIndx = 0; rIndx < nLyr; rIndx++)
                    {
                        for(idx = 0; idx < 8; idx++)
                        {
                            mse = hHer_RnnInv_h_inv_re[rIndx][rIndx][idx];
                            mse_inv = 1 / mse;
                            pNrPuschOutParams->msePerPrb[rIndx][idxPrb + idx] = mse_inv;
                        }
                    }
                }

                // Perform [(inv((H_her * Rnn_inv * H) + I)) * H_her]
                MATRIXMULH2(hHer_RnnInv_h_inv_re, hHer_RnnInv_h_inv_im, hHerMat_re, hHerMat_im, hHerRnnInvhInv_hHer_re, hHerRnnInvhInv_hHer_im, nLyr, nLyr, nRxAnt);

                // Perform W = [(inv((H_her * Rnn_inv * H) + I) * H_her) * Rnn_inv]
                MATRIXMULH2(hHerRnnInvhInv_hHer_re, hHerRnnInvhInv_hHer_im, Rnn_inv_re, Rnn_inv_im, weights_re, weights_im, nLyr, nRxAnt, nRxAnt);

                // Perform Equalization on each data symbol
                int equaSym = 0;
                for(idxSym = startSym; idxSym <= endSym; idxSym++)
                {
                    if(idxSym != dmrsSym)
                    {
                        int layerMapIdx,temp1 = (6*(startSym>0)+equaSym)*pNrPuschOutParams->total_nPRB*12*nLyr + completedPrb*12*nLyr;
                        for(idxAnt = 0; idxAnt < nRxAnt; idxAnt++)
                        {
                            Y_vect_re[idxAnt][0] = VECT_SET(rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[7] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[6] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[5] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[4] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[3] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[2] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[1] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_re[0] + idxSc] );

                            Y_vect_im[idxAnt][0] = VECT_SET(rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[7] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[6] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[5] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[4] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[3] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[2] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[1] + idxSc],
                                                            rxFdSamples[idxAnt][idxSym][prbScStartIdx_im[0] + idxSc] );
                        }

                        // y_tilda = W * Y
                        MATRIXMULH2(weights_re, weights_im, Y_vect_re, Y_vect_im, Y_tilda_re, Y_tilda_im, nLyr, nRxAnt, 1);

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

                            // Y_cap = Y_tilda * TOC
                            C_VECT_MUL_NA_NB(Y_tilda_re[idxLyr][0], Y_tilda_im[idxLyr][0], sq_toc_re[lyrType][idxSc], sq_toc_im[lyrType][idxSc], Y_cap_re, Y_cap_im);

                            // Store the Equalized output into memory
                            for(idx = 0; idx < 8; idx++)
                            {
                                layerMapIdx = temp1 + (prbScStartIdx_re[idx]-startPrbIdx*12+idxSc)*nLyr + idxLyr;
                                // pNrPuschOutParams->equaOutSamples_re[idxLyr][equaSym][prbScStartIdx_re[idx] + idxSc] = Y_cap_re[idx];
                                // pNrPuschOutParams->equaOutSamples_im[idxLyr][equaSym][prbScStartIdx_re[idx] + idxSc] = Y_cap_im[idx];
                                pNrPuschOutParams->layerDemapperOutReal[layerMapIdx] = Y_cap_re[idx];
                                pNrPuschOutParams->layerDemapperOutImag[layerMapIdx] = Y_cap_im[idx];
                            }
                        }

                        equaSym = equaSym + 1;
                    }
                }

            }

        }

        for(idxLyr = 0; idxLyr < nLyr; idxLyr++)
        {
            mseAvgVec = _mm256_setzero_ps();

            for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb+=8)
            {
                mseVec = _mm256_loadu_ps(&pNrPuschOutParams->msePerPrb[idxLyr][idxPrb]);

                mseAvgVec = _mm256_add_ps(mseAvgVec, mseVec);
            }

            mid_re = _mm256_permute2f128_ps (mseAvgVec, mseAvgVec, 1);
            mseAvgVec = _mm256_add_ps(mseAvgVec, mid_re);
            mseAvgVec = _mm256_hadd_ps(mseAvgVec, mseAvgVec);
            mseAvgVec = _mm256_hadd_ps(mseAvgVec, mseAvgVec);

            mseAvg = mseAvgVec[0];
            mseAvg_inv = 1 / mseAvg;

            mseAvgInvVec = _mm256_set1_ps(mseAvg_inv);

            for(idxPrb = startPrbIdx; idxPrb <= (endPrbIdx); idxPrb+=8)
            {
                mseVec = _mm256_loadu_ps(&pNrPuschOutParams->msePerPrb[idxLyr][idxPrb]);

                mseVec = _mm256_mul_ps(mseVec, mseAvgInvVec);

                _mm256_storeu_ps(&pNrPuschOutParams->msePerPrb[idxLyr][idxPrb], mseVec);
            }
        }
    }
}
*/
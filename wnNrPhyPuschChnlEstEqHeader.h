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
  * @brief This file includes all the Macros and Structure definitions
  *         for NR PUSCH Channel Estimation and Equalization.
  * @file common-def.h
  * @ingroup nr_pusch
  * @author WiSig Networks
  **/

#ifndef _COMMON_DEF_H_
#define	_COMMON_DEF_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <immintrin.h>

#include "wnNrPhyBwFftCpConfig.h"

#define MAX_SC_PER_RB       ( 12   )
#define MAX_NUM_SC          ( MAX_SUB_CARRIERS_OFDM_SYMBOL )
#define HALF_MAX_NUM_SC     ( MAX_SUB_CARRIERS_OFDM_SYMBOL/2 )
#define ONE_HALF_MAX_NUM_SC     ( 4950 )
#define TWICE_MAX_NUM_SC    ( MAX_SUB_CARRIERS_OFDM_SYMBOL*2 )

#define MAX_MOD_ORDER       ( 8 )

#define MAX_PUSCH_SYM       ( NUMBER_OF_OFDM_SYMBOLS_SLOT )
#define MAX_DMRS_SYM        ( 2 )
#define MAX_PUSCH_UE        ( 1 )
#define COV_GRP_FACTOR      ( 2 )

// #define LOW_PUSCH_SYM       ( 0 )
// #define HIGH_PUSCH_SYM      ( 1 )

#define SEG1                ( 1 )
#define SEG2                ( 2 )
#define SEG3                ( 3 )
#define SEG4                ( 4 )

#define EVEN_PUSCH_LYR      ( 0 )
#define ODD_PUSCH_LYR       ( 1 )


#define VECT_T          __m256
#define PER_VECTOR      8
#define VECT_SET1       _mm256_set1_ps
#define VECT_SET        _mm256_set_ps
#define VECT_STORE      _mm256_store_ps
#define VECT_MUL        _mm256_mul_ps
#define VECT_ADD        _mm256_add_ps
#define VECT_SUB        _mm256_sub_ps
#define VECT_RCP        _mm256_rcp_ps
#define VECT_MUL_ADD    _mm256_fmadd_ps
#define VECT_NMUL_ADD   _mm256_fnmadd_ps
#define VECT_MUL_SUB    _mm256_fmsub_ps
#define VECT_SQRT       _mm256_sqrt_ps
#define VECT_NEG(xxx)   _mm256_xor_ps(g_negate, (xxx))
#define VECT_SETZERO    _mm256_setzero_ps

// #define ALIGN_BOUNDARY  64
#define UNROLL_PRAGMA   _Pragma("unroll")


float rxFdSamples[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_PUSCH_SYM][2 * MAX_NUM_SC];
// float layerDemapperOutReal[MAX_PUSCH_LAYERS*12*MAX_NUM_SC], layerDemapperOutImag[MAX_PUSCH_LAYERS*12*MAX_NUM_SC];
int totalLength;
__attribute__((aligned(32))) float llrOut[12*MAX_NUM_PRB*13*MAX_PUSCH_LAYERS*MAX_MOD_ORDER];//12
float llrOutScaledFlt[12*MAX_NUM_PRB*13*MAX_PUSCH_LAYERS*MAX_MOD_ORDER];
__attribute__((aligned(32)))  wnInt8 llrFxd[12*MAX_NUM_PRB*13*MAX_PUSCH_LAYERS*MAX_MOD_ORDER];//12



typedef struct nrPuschOutParams_t
{
    int total_nPRB;
    wnUInt16 nEvenTonesLyr;
    wnUInt16 nOddTonesLyr;

    float genDMRS[2][MAX_SUB_CARRIERS_OFDM_SYMBOL];
    // float genDMRS2[MAX_SUB_CARRIERS_OFDM_SYMBOL];

    float lsEst[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_PUSCH_LAYERS][MAX_SUB_CARRIERS_OFDM_SYMBOL];
    // Channel estimates after per Linear Interpolation.
    float chEst_re[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_PUSCH_LAYERS][12 * MAX_NUM_PRB];// For WO IRC Per Tone
    float chEst_im[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_PUSCH_LAYERS][12 * MAX_NUM_PRB];// For WO IRC Per Tone
    float chEstAvg_re[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_PUSCH_LAYERS][6 * MAX_NUM_PRB];
    float chEstAvg_im[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_PUSCH_LAYERS][6 * MAX_NUM_PRB];


    // float chEstPerTone_re[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_PUSCH_LAYERS][MAX_NUM_SC];
    // float chEstPerTone_im[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_PUSCH_LAYERS][MAX_NUM_SC];   

    // float RnnEven_prb_re[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_SUMIMO_PUSCH_ANTENNA_PORTS][2 * MAX_NUM_PRB];
    // float RnnEven_prb_im[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_SUMIMO_PUSCH_ANTENNA_PORTS][2 * MAX_NUM_PRB];

    // float RnnOdd_prb_re[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_SUMIMO_PUSCH_ANTENNA_PORTS][2 * MAX_NUM_PRB];
    // float RnnOdd_prb_im[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_SUMIMO_PUSCH_ANTENNA_PORTS][2 * MAX_NUM_PRB];

    // float Rnn_prb_re[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_NUM_PRB];
    // float Rnn_prb_im[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_NUM_PRB];

    // float RnnInv_prb_re[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_NUM_PRB];
    // float RnnInv_prb_im[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_NUM_PRB];

    wnFlt oneBymse_re[MAX_PUSCH_LAYERS][12 * MAX_NUM_PRB];//For wo IRC Modules
    wnFlt oneBymse_avg_re[MAX_PUSCH_LAYERS];//For wo IRC Modules

    // float msePerPrb[MAX_PUSCH_LAYERS][2 * MAX_NUM_PRB];//For IRC Modules
    
    wnFlt equaOutSamples[MAX_PUSCH_LAYERS][13][2 * MAX_SUB_CARRIERS_OFDM_SYMBOL];// For WO IRC

    // float equaOutSamples_re[MAX_PUSCH_LAYERS][6][MAX_NUM_SC];// For IRC
    // float equaOutSamples_im[MAX_PUSCH_LAYERS][6][MAX_NUM_SC];// For IRC

   __attribute__((aligned(32))) float layerDemapperOutReal[MAX_PUSCH_LAYERS*13*MAX_NUM_SC];//7
   __attribute__((aligned(32))) float layerDemapperOutImag[MAX_PUSCH_LAYERS*13*MAX_NUM_SC];//7

    // float nVar[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][MAX_PUSCH_LAYERS][MAX_NUM_ALLOC];
    float nVar_overAll; // For WO IRC Per Tone
    // float nVar_overAll[MAX_NUM_PRB / 8];
    // float toc_overAll[2];
    // float toc_overAll_re[MAX_NUM_PRB / 8];
    // float toc_overAll_im[MAX_NUM_PRB / 8];

} __attribute__((__packed__)) NR_PUSCH_OUT_PARAMS, __attribute__((__packed__)) *P_NR_PUSCH_OUT_PARAMS;

/** Output Configuration Parameters from Channel Equalizer(Low = 1-7 Symbol data,High = 8-14 Symbol data) .*/
NR_PUSCH_OUT_PARAMS nrPuschOutParamsSeg1;
NR_PUSCH_OUT_PARAMS nrPuschOutParamsSeg2;
NR_PUSCH_OUT_PARAMS nrPuschOutParamsSeg3;
NR_PUSCH_OUT_PARAMS nrPuschOutParamsSeg4;



wnFlt FOEst;
wnFlt FOEst_PRG[MAX_NUM_PRB];

VECT_T g_negate;
VECT_T g_mzero2x2[2][2];
VECT_T g_mzero4x4[4][4];
VECT_T g_mzero8x8[8][8];

// typedef struct
// {
//     float real;
//     float imag;
// }cmplx;

//MACRO for Complex Vectors Multiplications
#define C_VECT_MUL_NA_NB(a_re, a_im, b_re, b_im, c_re, c_im) \
{ \
   VECT_T mid_re; \
   VECT_T mid_im; \
\
   mid_re = VECT_MUL(a_im, b_im); \
   mid_im = VECT_MUL(a_im, b_re); \
\
   c_re = VECT_MUL_SUB(a_re, b_re, mid_re); \
   c_im = VECT_MUL_ADD(a_re, b_im, mid_im); \
}

// MACRO to multiply Nomral matrix with Normal Matrix
#define MATRIXMULNN(a_re, a_im, b_re, b_im, d_re, d_im, l_sz, m_sz, n_sz) \
{ \
    int ii; \
    int ij; \
    int ik; \
    UNROLL_PRAGMA; \
    for (ii = 0; ii < l_sz; ii++) \
    { \
        for (ij = 0; ij < n_sz; ij++) \
        { \
            d_re[ii][ij] = d_im[ii][ij] = VECT_SETZERO(); \
\
            for (ik = 0; ik < m_sz; ik++) \
            { \
                d_re[ii][ij] = VECT_MUL_ADD(a_re[ii][ik], b_re[ik][ij], d_re[ii][ij]); \
                d_re[ii][ij] = VECT_NMUL_ADD(a_im[ii][ik], b_im[ik][ij], d_re[ii][ij]); \
\
                d_im[ii][ij] = VECT_MUL_ADD(a_re[ii][ik], b_im[ik][ij], d_im[ii][ij]); \
                d_im[ii][ij] = VECT_MUL_ADD(a_im[ii][ik], b_re[ik][ij], d_im[ii][ij]); \
            } \
        } \
    } \
}

// MACRO to multiply Nomral matrix with Hermitian Matrix
#define MATRIXMULNH(a_re, a_im, b_re, b_im, d_re, d_im, l_sz, m_sz, n_sz) \
{ \
    int ii; \
    int ij; \
    int ik; \
    UNROLL_PRAGMA; \
    for (ii = 0; ii < l_sz; ii++) { \
        for (ij = 0; ij < n_sz; ij++) { \
            d_re[ii][ij] = d_im[ii][ij] = VECT_SETZERO(); \
            for (ik = 0; ik < m_sz; ik++) { \
                if (ik != ij) { \
                    d_re[ii][ij] = VECT_MUL_ADD(a_re[ii][ik], \
                                                b_re[ik][ij], d_re[ii][ij]); \
                    d_re[ii][ij] = VECT_NMUL_ADD(a_im[ii][ik], \
                                                 b_im[ik][ij], d_re[ii][ij]); \
                    d_im[ii][ij] = VECT_MUL_ADD(a_re[ii][ik], \
                                                b_im[ik][ij], d_im[ii][ij]); \
                    d_im[ii][ij] = VECT_MUL_ADD(a_im[ii][ik], \
                                                b_re[ik][ij], d_im[ii][ij]); \
                } else { \
                    d_re[ii][ij] = VECT_MUL_ADD(a_re[ii][ik], \
                                                b_re[ik][ik], d_re[ii][ij]); \
                    d_im[ii][ij] = VECT_MUL_ADD(a_im[ii][ik], \
                                                b_re[ik][ik], d_im[ii][ij]); \
                } \
            } \
        } \
    } \
}

// MACRO to multiply Hermitian Matrix with Normal matrix
#define MATRIXMULHN(a_re, a_im, b_re, b_im, d_re, d_im, l_sz, m_sz, n_sz) \
{ \
    int ii; \
    int ij; \
    int ik; \
    UNROLL_PRAGMA; \
    for (ii = 0; ii < l_sz; ii++) { \
        for (ij = 0; ij < n_sz; ij++) { \
            d_re[ii][ij] = d_im[ii][ij] = VECT_SETZERO(); \
            for (ik = 0; ik < m_sz; ik++) { \
                if (ii != ik) { \
                    d_re[ii][ij] = VECT_MUL_ADD(a_re[ii][ik], \
                                                b_re[ik][ij], d_re[ii][ij]); \
                    d_re[ii][ij] = VECT_NMUL_ADD(a_im[ii][ik], \
                                                 b_im[ik][ij], d_re[ii][ij]); \
                    d_im[ii][ij] = VECT_MUL_ADD(a_re[ii][ik], \
                                                b_im[ik][ij], d_im[ii][ij]); \
                    d_im[ii][ij] = VECT_MUL_ADD(a_im[ii][ik], \
                                                b_re[ik][ij], d_im[ii][ij]); \
                } else { \
                    d_re[ii][ij] = VECT_MUL_ADD(a_re[ii][ii], \
                                                b_re[ik][ij], d_re[ii][ij]); \
                    d_im[ii][ij] = VECT_MUL_ADD(a_re[ii][ii], \
                                                b_im[ik][ij], d_im[ii][ij]); \
                } \
            } \
        } \
    } \
}

#define MATRIXMULH2(a_re, a_im, b_re, b_im, d_re, d_im, l_sz, m_sz, n_sz) \
{ \
    int ii; \
    int ij; \
    int ik; \
    UNROLL_PRAGMA; \
    for (ii = 0; ii < l_sz; ii++) { \
        for (ij = 0; ij < n_sz; ij++) { \
            d_re[ii][ij] = d_im[ii][ij] = VECT_SETZERO(); \
            for (ik = 0; ik < m_sz; ik++) { \
                d_re[ii][ij] = VECT_MUL_ADD(a_re[ii][ik], \
                                            b_re[ik][ij], d_re[ii][ij]); \
                d_re[ii][ij] = VECT_NMUL_ADD(a_im[ii][ik], \
                                             b_im[ik][ij], d_re[ii][ij]); \
                d_im[ii][ij] = VECT_MUL_ADD(a_re[ii][ik], \
                                            b_im[ik][ij], d_im[ii][ij]); \
                d_im[ii][ij] = VECT_MUL_ADD(a_im[ii][ik], \
                                            b_re[ik][ij], d_im[ii][ij]); \
            } \
        } \
    } \
}



// Result of matrix multiplication is a Hermitian matrix
// Calculate the upper triangular of the result and populate the lower triangular elements
// Diagonal values of the result will be real, so while doing complex multiplication, only calculate the real part and assign imaginary part as 0
//VECT_T temp_re;
//VECT_T temp_im; 
#define MATRIXMULTHT(a_re, a_im, b_re, b_im, d_re, d_im, l_sz, m_sz, n_sz) \
{ \
    int ii; \
    int ij; \
    int ik; \
    UNROLL_PRAGMA; \
    for (ii = 0; ii < l_sz; ii++) { \
        for (ij = ii + 1; ij < n_sz; ij++) { \
            d_re[ii][ij] = d_im[ii][ij] = VECT_SETZERO(); \
            for (ik = 0; ik < m_sz; ik++) { \
                d_re[ii][ij] = VECT_MUL_ADD(a_re[ii][ik], \
                                            b_re[ik][ij], d_re[ii][ij]); \
                d_re[ii][ij] = VECT_NMUL_ADD(a_im[ii][ik], \
                                             b_im[ik][ij], d_re[ii][ij]); \
                d_im[ii][ij] = VECT_MUL_ADD(a_re[ii][ik], \
                                            b_im[ik][ij], d_im[ii][ij]); \
                d_im[ii][ij] = VECT_MUL_ADD(a_im[ii][ik], \
                                            b_re[ik][ij], d_im[ii][ij]); \
            } \
            d_re[ij][ii] = d_re[ii][ij]; \
            d_im[ij][ii] = VECT_NEG(d_im[ii][ij]);  \
        } \
    } \
    UNROLL_PRAGMA; \
    for (ii = 0; ii < l_sz; ii++) { \
        d_re[ii][ii] = d_im[ii][ii] = VECT_SETZERO(); \
        for (ik = 0; ik < m_sz; ik++) { \
            d_re[ii][ii] = VECT_MUL_ADD(a_re[ii][ik], \
                                        b_re[ik][ii], d_re[ii][ii]); \
            d_re[ii][ii] = VECT_NMUL_ADD(a_im[ii][ik], \
                                         b_im[ik][ii], d_re[ii][ii]); \
        } \
    } \
}

// MACRO to add two matrices
#define MATRIXADD(a_re, a_im, b_re, b_im, d_re, d_im, l_sz, m_sz) \
{ \
    int ii; \
    int ij; \
    UNROLL_PRAGMA; \
    for (ii = 0; ii < l_sz; ii++) { \
        for (ij = 0; ij < m_sz; ij++) { \
            d_re[ii][ij] = VECT_ADD(a_re[ii][ij], b_re[ii][ij]); \
            d_im[ii][ij] = VECT_ADD(a_im[ii][ij], b_im[ii][ij]); \
        } \
    } \
}

// MACRO to subtract two matrices
#define MATRIXSUB(a_re, a_im, b_re, b_im, d_re, d_im, l_sz, m_sz) \
{ \
    int ii; \
    int ij; \
    UNROLL_PRAGMA; \
    for (ii = 0; ii < l_sz; ii++) { \
        for (ij = 0; ij < m_sz; ij++) { \
            d_re[ii][ij] = VECT_SUB(a_re[ii][ij], b_re[ii][ij]); \
            d_im[ii][ij] = VECT_SUB(a_im[ii][ij], b_im[ii][ij]); \
        } \
    } \
}


// Performs A*A' = A*conj(A).' (result is hermitian matrix)
// Calculate the upper triangular of the result and populate the lower triangular elements
// Diagonal values of the result will be real, so while doing complex multiplication, only calculate the real part and assign imaginary part as 0
#define MATRIXMULTTH(a_re, a_im, d_re, d_im, l_sz, m_sz, n_sz) \
{ \
    int ii; \
    int ij; \
    int ik; \
\
    UNROLL_PRAGMA; \
    for (ii = 0; ii < l_sz; ii++) \
    { \
        for (ij = ii + 1; ij < n_sz; ij++) \
        { \
            d_re[ii][ij] = VECT_MUL(a_re[ii][0], a_re[ij][0]); \
            d_re[ii][ij] = VECT_MUL_ADD(a_im[ii][0], a_im[ij][0], d_re[ii][ij]); \
            d_im[ii][ij] = VECT_MUL(a_im[ii][0], a_re[ij][0]); \
            d_im[ii][ij] = VECT_NMUL_ADD(a_re[ii][0], a_im[ij][0], d_im[ii][ij]); \
\
            for (ik = 1; ik < m_sz; ik++) \
            { \
                d_re[ii][ij] =  VECT_MUL_ADD(a_re[ii][ik], a_re[ij][ik], d_re[ii][ij]); \
                d_re[ii][ij] =  VECT_MUL_ADD(a_im[ii][ik], a_im[ij][ik], d_re[ii][ij]); \
                d_im[ii][ij] = VECT_NMUL_ADD(a_re[ii][ik], a_im[ij][ik], d_im[ii][ij]); \
                d_im[ii][ij] =  VECT_MUL_ADD(a_im[ii][ik], a_re[ij][ik], d_im[ii][ij]); \
            } \
            d_re[ij][ii] = d_re[ii][ij]; \
            d_im[ij][ii] = VECT_NEG(d_im[ii][ij]);  \
        } \
    } \
        \
    d_im[0][0] = VECT_SETZERO(); \
    d_re[0][0] = VECT_MUL(a_re[0][0], a_re[0][0]); \
    d_re[0][0] = VECT_MUL_ADD(a_im[0][0], a_im[0][0], d_re[0][0]);\
\
    for (ik = 1; ik < m_sz; ik++) \
    { \
        d_re[0][0] = VECT_MUL_ADD(a_re[0][ik], a_re[0][ik], d_re[0][0]); \
        d_re[0][0] = VECT_MUL_ADD(a_im[0][ik], a_im[0][ik], d_re[0][0]); \
    } \
\
    UNROLL_PRAGMA; \
    for (ii = 1; ii < l_sz; ii++) \
    { \
        d_im[ii][ii] = VECT_SETZERO(); \
        d_re[ii][ii] = VECT_MUL(a_re[ii][0], a_re[ii][0]); \
        d_re[ii][ii] = VECT_MUL_ADD(a_im[ii][0], a_im[ii][0], d_re[ii][ii]);\
\
        for (ik = 1; ik < m_sz; ik++) \
        { \
            d_re[ii][ii] = VECT_MUL_ADD(a_re[ii][ik], a_re[ii][ik], d_re[ii][ii]); \
            d_re[ii][ii] = VECT_MUL_ADD(a_im[ii][ik], a_im[ii][ik], d_re[ii][ii]); \
        } \
    } \
}

#define MATRIXMULTTH2(a_re, a_im, d_re, d_im, l_sz, m_sz, n_sz) \
{ \
    int ii; \
    int ij; \
    int ik; \
\
    UNROLL_PRAGMA; \
    for (ii = 0; ii < l_sz; ii++) \
    { \
        for (ij = ii + 1; ij < n_sz; ij++) \
        { \
            d_re[ii][ij] = VECT_MUL(a_re[ii][0], a_re[ij][0]); \
            d_re[ii][ij] = VECT_MUL_ADD(a_im[ii][0], a_im[ij][0], d_re[ii][ij]); \
            d_im[ii][ij] = VECT_MUL(a_im[ii][0], a_re[ij][0]); \
            d_im[ii][ij] = VECT_NMUL_ADD(a_re[ii][0], a_im[ij][0], d_im[ii][ij]); \
\
            for (ik = 1; ik < m_sz; ik++) \
            { \
                d_re[ii][ij] =  VECT_MUL_ADD(a_re[ii][ik], a_re[ij][ik], d_re[ii][ij]); \
                d_re[ii][ij] =  VECT_MUL_ADD(a_im[ii][ik], a_im[ij][ik], d_re[ii][ij]); \
                d_im[ii][ij] = VECT_NMUL_ADD(a_re[ii][ik], a_im[ij][ik], d_im[ii][ij]); \
                d_im[ii][ij] =  VECT_MUL_ADD(a_im[ii][ik], a_re[ij][ik], d_im[ii][ij]); \
            } \
            d_re[ij][ii] = d_re[ii][ij]; \
            d_im[ij][ii] = VECT_NEG(d_im[ii][ij]);  \
        } \
    } \
        \
}

// A'*B
// A and B are assumed to be transposed
#define MATRIXMULTHTTRANSPOSE(a_re, a_im, b_re, b_im, d_re, d_im, l_sz, m_sz, n_sz) \
{ \
    int ii; \
    int ij; \
    int ik; \
    UNROLL_PRAGMA; \
    for (ii = 0; ii < l_sz; ii++) \
    { \
        for (ij = 0; ij < n_sz; ij++) \
        { \
            d_re[ii][ij] = VECT_MUL(a_re[0][ii], b_re[ij][0]); \
            d_re[ii][ij] = VECT_MUL_ADD(a_im[0][ii], b_im[ij][0], d_re[ii][ij]); \
            d_im[ii][ij] = VECT_MUL(a_re[0][ii], b_im[ij][0]); \
            d_im[ii][ij] = VECT_NMUL_ADD(a_im[0][ii], b_re[ij][0], d_im[ii][ij]); \
            for (ik = 1; ik < m_sz; ik++) \
            { \
                d_re[ii][ij] = VECT_MUL_ADD(a_re[ik][ii], b_re[ij][ik], d_re[ii][ij]); \
                d_re[ii][ij] = VECT_MUL_ADD(a_im[ik][ii], b_im[ij][ik], d_re[ii][ij]); \
                d_im[ii][ij] = VECT_MUL_ADD(a_re[ik][ii], b_im[ij][ik], d_im[ii][ij]); \
                d_im[ii][ij] = VECT_NMUL_ADD(a_im[ik][ii], b_re[ij][ik], d_im[ii][ij]); \
            } \
        } \
    } \
}

#define MATRIXMULTHTTRANSPOSE2(a_re, a_im, b_re, b_im, d_re, d_im, l_sz, m_sz, n_sz) \
{ \
    int ii; \
    int ij; \
    int ik; \
    UNROLL_PRAGMA; \
    for (ii = 0; ii < l_sz; ii++) \
    { \
        for (ij = 0; ij < n_sz; ij++) \
        { \
            d_re[ii][ij] = VECT_MUL(a_re[0][ii], b_re[0][ij]); \
            d_re[ii][ij] = VECT_MUL_ADD(a_im[0][ii], b_im[0][ij], d_re[ii][ij]); \
            d_im[ii][ij] = VECT_MUL(a_re[0][ii], b_im[0][ij]); \
            d_im[ii][ij] = VECT_NMUL_ADD(a_im[0][ii], b_re[0][ij], d_im[ii][ij]); \
            for (ik = 1; ik < m_sz; ik++) \
            { \
                d_re[ii][ij] = VECT_MUL_ADD(a_re[ik][ii], b_re[ik][ij], d_re[ii][ij]); \
                d_re[ii][ij] = VECT_MUL_ADD(a_im[ik][ii], b_im[ik][ij], d_re[ii][ij]); \
                d_im[ii][ij] = VECT_MUL_ADD(a_re[ik][ii], b_im[ik][ij], d_im[ii][ij]); \
                d_im[ii][ij] = VECT_NMUL_ADD(a_im[ik][ii], b_re[ik][ij], d_im[ii][ij]); \
            } \
        } \
    } \
}

#include "matrix_invert.inc"

#endif // _COMMON_DEF_H_

/*
** Copyright (c) 2016-2018, WiSig Networks Pvt Ltd. All rights reserved.
** www.wisig.com
**
** All information contained herein is property of WiSig Networks Pvt Ltd.
** unless otherwise explicitly mentioned.
**
** The intellectual and technical concepts in this file are proprietary
** to WiSig Networks and may be covered by granted or in process national
** and international patents and are protect by trade secrets and
** copyright law.
**
** Redistribution and use in source and binary forms of the content in
** this file, with or without modification are not permitted unless
** permission is explicitly granted by WiSig Networks.
** If WiSig Networks permits this source code to be used as a part of
** open source project, the terms and conditions of CC-By-ND (No Derivative) license
** (https://creativecommons.org/licenses/by-nd/4.0/legalcode) shall apply.
*/


// NOTE --- To avoid Stack Smashing or other Compilation Errors, make sure that size of output array is a multiple of 32.

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "wnNrPhyPuschHeader.h"

/*TS 38.211 (Release 15) Section 5.2.1 - Pseudo-random sequence generation*/
wnVoid wnNrPhyPuschPnSequence(wnInt64 seedVal,    // Seed value for PN sequence generation
                  	  	  	  wnUInt32           seqLen,     // Length of the PN sequence to be generated
							  wnUInt8         *pnSeqOut)   // Output PN sequence
{
    /* Each bit of x1 and x2 represents a value corresponding to the m-sequences defined in 38.211 Sec-5.2.1*/
    wnUInt32 x1, x2, i;

    /* seq1(0) = 1, seq1(i) = 0 for 0 < i < 31 and seq1(n + 31) = (seq1(n + 3) + seq1(n))mod2*/
    /* Seed value of seq1 is fixed, so the values are pre-calculated and stored*/
    x1 =  0x54D21B24;  /* seq1(n), 1600 <= n <= 1630 (Since Nc for 5G NR is 1600)*/

    wnUInt8 *k = pnSeqOut;

    /* seq2(n), 1600 <= n <= 1631(Since Nc for 5G NR is 1600)*/
    x2 = 0;
    for (i = 0; seedVal>>i ; i++)
    {
        if ((seedVal>>i)&1)
            x2 ^= SEQUENCE_X2_MASKS[i];
    }

    wnUInt32 tmp1, tmp2;
    wnUInt32  N_iter, Ntail;

    N_iter = seqLen / 24;
    Ntail = seqLen % 24;

    for (i = 0; i < N_iter; i++)
    {
        tmp1 = x1 ^ (x1 >> 3);
        x1 >>= 24;
        x1 = (x1&2147483775)|((tmp1&16777215)<<7);

        tmp2 = x2 ^ (x2 >> 1) ^ (x2 >> 2) ^ (x2 >> 3);
        x2 >>= 24;
        x2 = (x2&2147483775)|((tmp2&16777215)<<7);

        *(wnUInt32*) pnSeqOut = tmp1 ^ tmp2;
        pnSeqOut += 3;

    }

    if (Ntail) {
        x1 &= 0x7fffffff;
        tmp1 = x1 ^ (x1 >> 3);
        x1 >>= Ntail;
        x1 |= (tmp1 << (31 - Ntail));

        x2 &= 0x7fffffff;
        tmp2 = x2 ^ (x2 >> 1) ^ (x2 >> 2) ^ (x2 >> 3);
        x2 >>= Ntail;
        x2 |= (tmp2 << (31 - Ntail));

        *(wnUInt32*) pnSeqOut = tmp1 ^ tmp2;

    }
//printf("%d\n",pnSeqOut);
#ifdef DEBUG
FILE* fid= fopen("R.txt","w");
for(i = 0; i < (seqLen>>3); i++)
    for (wnUInt32 j = 0; j < 8; j++)
           fprintf(fid,"%d\n",(k[i]>>j)&1);
#endif
return;
}

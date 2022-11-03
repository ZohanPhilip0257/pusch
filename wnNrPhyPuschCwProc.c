#include "wnNrPhyPuschHeader.h"
#include "wnNrPhygNBUlChnlConfigInit.h"

/************************************************************
 * Input Arguments:                                         *
 * payload --- Input data stacked into bytes                *
 * numBits --- Length of data in bits                       *
 *                                                          *
 * Output Arguments:                                        *
 * crc  --- Output in stacked in bytes for CRC16,CRC24A     *
 *             Output in bits for CRC24B, CRC24C            *
 ************************************************************/


/* -----------------------   CRC16    ----------------------- */
wnVoid wnNrPhyCrc16Removal(wnUInt8 *payload, wnInt32 numBits, wnUInt8 *crcRemoved, wnUInt8 *check)
{

    wnInt32 iDx,numBytes = (numBits-16)>>3;                                        //Indices used in Loops
    wnUInt16 remainder = 0;
    wnUInt8 check1, check2;

    //Start of CRC operation
    for (iDx = 0; iDx < numBytes; iDx++)
    {
        remainder = ( remainder << 8 ) ^ CRC16_TAB[( (wnInt32) (remainder >> 8) ^ payload[iDx])& 0xff];
    }
    //End of CRC operation

    // Checks whether the CRC decoded matches or not
    //CRC Remainder
    check1 = (remainder & 0x00ff) - payload[numBytes+1];
    check2 = ((remainder>>8) & 0x00ff) - payload[numBytes];

    if ((check1 == 0) & (check2 == 0))
        *check = 0;
    else
        *check = 1;
       
    memcpy (crcRemoved, payload, numBytes);

    return;
}
/*-------------------------------------------------------------*/



/* -----------------------   CRC24A    ----------------------- */
wnVoid wnNrPhyCrc24ARemoval(wnUInt8 *payload, wnInt32 numBits, wnUInt8 *crcRemoved, wnUInt8 *check)
{
    wnInt32 iDx,numBytes = (numBits-24)>>3;                                        //Indices used in Loops
    wnUInt32 remainder = 0;
    wnUInt8 check1, check2, check3;

    //Start of CRC operation
    for (iDx = 0; iDx < numBytes; iDx++)
    {
        remainder = ( remainder << 8 ) ^ CRC24A_TAB[( (wnInt32) (remainder >> 24) ^ (payload[iDx]))& 0xff];
    }
    //End of CRC operation

    // Checks whether the CRC decoded matches or not
    //CRC Remainder
    check1 = ((remainder>>8) & 0x000000ff) - payload[numBytes+2];
    check2 = ((remainder>>16) & 0x000000ff) - payload[numBytes+1];
    check3 = ((remainder>>24) & 0x000000ff) - payload[numBytes];

    if ((check1 == 0) & (check2 == 0) & (check3 == 0))
        *check = 0;
    else
        *check = 1;    

    memcpy (crcRemoved, payload, numBytes);
    //printf ("%d\n", *check);
    return;
}
/*-------------------------------------------------------------*/




/* -----------------------   CRC24B    ----------------------- */
wnVoid wnNrPhyCrc24BRemoval(wnUInt8 *payload, wnInt32 numBits, wnUInt8 *crcRemoved, wnUInt8 *check)
{
    wnInt32 iDx,extraBits,numBytes = (numBits-24)>>3;
    wnInt32 pdsch_crc_comp=16777216,pdsch_crc=25165923;
    wnUInt32 remainder = 0;
    wnUInt8 check1, check2, check3;

    //Start of CRC operation
    for (iDx = 0; iDx < numBytes; iDx++)
    {
        remainder = ( remainder << 8 ) ^ CRC24B_TAB[( (wnInt32) (remainder >> 24) ^ (payload[iDx]))& 0xff];
    }
    remainder = remainder>>8;

    extraBits = (numBits-24)&0x7;

    if(extraBits)
    {
        remainder = (remainder ^ (payload[iDx]<<16)) & 0x00ffffff;
        
        while (extraBits > 0) 
        {
            remainder = remainder<<1;
            remainder = remainder ^ (pdsch_crc * (remainder >= pdsch_crc_comp));        
            extraBits--;
        }
    }

    //End of CRC operation

    // Checks whether the CRC decoded matches or not
    //CRC Remainder
    check1 = ((remainder) & 0x000000ff) - payload[numBytes+2];
    check2 = ((remainder>>8) & 0x000000ff) - payload[numBytes+1];
    check3 = ((remainder>>16) & 0x000000ff) - payload[numBytes];

    if ((check1 == 0) & (check2 == 0) & (check3 == 0))
        *check = 0;
    else
        *check = 1;    

    memcpy (crcRemoved, payload, numBytes);
    //printf ("%d\n", *check);
  
    return;
}
/*-------------------------------------------------------------*/



/* -----------------------   CRC24C    ----------------------- */
wnVoid wnNrPhyCrc24CRemoval(wnUInt8 *payload, wnInt32 numBits, wnUInt8 *crcRemoved, wnUInt8 *check)
{
    wnInt32 iDx,extraBits,numBytes = (numBits-24)>>3;
    wnInt32 pdsch_crc_comp=16777216,pdsch_crc=28487959;
    wnInt32 remainder = 0;
    wnUInt8 check1, check2, check3;

    //Start of CRC operation
    for (iDx = 0; iDx < numBytes; iDx++)
    {
        remainder = ( remainder << 8 ) ^ CRC24C_TAB[( (wnInt32) (remainder >> 24) ^ payload[iDx])& 0xff];
    }
    remainder = remainder>>8;

    extraBits = numBits - (numBytes<<3);
    if(extraBits)
    {
        remainder = (remainder ^ (payload[iDx]<<16)) & 0x00ffffff;
        
        while (extraBits > 0) 
        {
            remainder = remainder<<1;
            remainder = remainder ^ (pdsch_crc * (remainder >= pdsch_crc_comp));        
            extraBits--;
        }
    }
    //End of CRC operation

    // Checks whether the CRC decoded matches or not
    //CRC Remainder
    check1 = ((remainder) & 0x000000ff) - payload[numBytes+2];
    check2 = ((remainder>>8) & 0x000000ff) - payload[numBytes+1];
    check3 = ((remainder>>16) & 0x000000ff) - payload[numBytes];

    if ((check1 == 0) & (check2 == 0) & (check3 == 0))
        *check = 0;
    else
        *check = 1;    

    memcpy (crcRemoved, payload, numBytes);    
  
    return;
}
/*-------------------------------------------------------------*/

// Rate DeMatching
wnVoid wnNrPhyPuschRateDematching(wnInt8 scrambleOut[], wnInt8 ldpcIn[], wnInt32 Ncb, wnInt32 start_pos_bit_select, wnInt32 fillerStart, wnInt32 fillerEnd, wnInt32 rmLength, wnInt8 modOrder)
{
    __m256i inter;//, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, inter1, inter2, inter3, excessM[9];
    wnUInt32 iDx, iDx2, iDx3, interval, interval2, interval3, interval4, interval5, interval6, interval7, excess;
    wnInt8 *deInterleaved;

    switch (modOrder)
    {
        /***************** qpsk ************************/
        case 2:
        
        iDx2 = (rmLength>>5)<<5;
        iDx3 = iDx2>>1;

        interval = rmLength>>1;

        excess = (rmLength & 0x1f)>>1;

        __m256i qpskShift1 = _mm256_set_epi32(7,6,3,2,5,4,1,0);
        __m256i qpskShift2 = _mm256_set_epi8(15,13,11,9,7,5,3,1,14,12,10,8,6,4,2,0,15,13,11,9,7,5,3,1,14,12,10,8,6,4,2,0);

        for (iDx = 0; iDx < iDx3 ; iDx=iDx+16)
        {
            inter = _mm256_lddqu_si256 ((__m256i *)&scrambleOut[iDx<<1]);

            inter = _mm256_shuffle_epi8 (inter, qpskShift2);
            inter = _mm256_permutevar8x32_epi32 (inter, qpskShift1);
            
            deInterleaved = (wnInt8 *)&inter;

            memcpy (&bitSelectIn[iDx], &deInterleaved[0], 16);
            memcpy (&bitSelectIn[interval + iDx], &deInterleaved[16], 16);
        }

        if (excess)
        {
            inter = _mm256_lddqu_si256 ((__m256i *)&scrambleOut[iDx<<1]);

            inter = _mm256_shuffle_epi8 (inter, qpskShift2);
            inter = _mm256_permutevar8x32_epi32 (inter, qpskShift1);
            
            deInterleaved = (wnInt8 *)&inter;

            memcpy (&bitSelectIn[iDx], &deInterleaved[0], excess);
            memcpy (&bitSelectIn[interval + iDx], &deInterleaved[16], excess);
        }
        
        break;

        /************************** 16 QAM *********************/
        case 4:
        
        iDx2 = (rmLength>>5)<<5;
        iDx3 = iDx2>>2;

        interval = rmLength>>2;
        interval2 = interval<<1;
        interval3 = interval + interval2;

        excess = (rmLength & 0x1f)>>2;

        __m256i q16Shift1 = _mm256_set_epi32(7,3,6,2,5,1,4,0);
        __m256i q16Shift2 = _mm256_set_epi8(15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0,15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0);

        for (iDx = 0; iDx < iDx3 ; iDx=iDx+8)
        {
            inter = _mm256_lddqu_si256 ((__m256i *)&scrambleOut[iDx<<2]);

            inter = _mm256_shuffle_epi8 (inter, q16Shift2);
            inter = _mm256_permutevar8x32_epi32 (inter, q16Shift1);
            
            deInterleaved = (wnInt8 *)&inter;

            memcpy (&bitSelectIn[iDx], &deInterleaved[0], 8);
            memcpy (&bitSelectIn[interval + iDx], &deInterleaved[8], 8);
            memcpy (&bitSelectIn[interval2 + iDx], &deInterleaved[16], 8);
            memcpy (&bitSelectIn[interval3 + iDx], &deInterleaved[24], 8);
        }

        if (excess)
        {
            inter = _mm256_lddqu_si256 ((__m256i *)&scrambleOut[iDx<<2]);

            inter = _mm256_shuffle_epi8 (inter, q16Shift2);
            inter = _mm256_permutevar8x32_epi32 (inter, q16Shift1);
            
            deInterleaved = (wnInt8 *)&inter;

            memcpy (&bitSelectIn[iDx], &deInterleaved[0], excess);
            memcpy (&bitSelectIn[interval + iDx], &deInterleaved[8], excess);
            memcpy (&bitSelectIn[interval2 + iDx], &deInterleaved[16], excess);
            memcpy (&bitSelectIn[interval3 + iDx], &deInterleaved[24], excess);
        }

        break;


        /************************** 64 QAM *********************/
        case 6:

        iDx2 = (rmLength/24)*24;
        iDx3 = iDx2/6;

        interval = rmLength/6;
        interval2 = interval<<1;
        interval3 = interval + interval2;
        interval4 = interval<<2;
        interval5 = interval + interval4;
        
        excess = (rmLength % 24)/6;

        __m256i q64Shift1 = _mm256_set_epi32(7,6,5,3,4,2,1,0);
        __m256i q64Shift2a = _mm256_set_epi8(15,14,13,12,11,10,9,8,7,6,1,0,5,4,3,2,11,5,10,4,15,9,3,14,8,2,13,7,1,12,6,0);
        __m256i q64Shift2b = _mm256_set_epi8(15,14,13,12,11,10,9,8,7,5,3,2,6,4,1,0,15,11,10,9,14,8,7,6,13,5,4,3,12,2,1,0);

        for (iDx = 0; iDx < iDx3 ; iDx=iDx+4)
        {
            inter = _mm256_lddqu_si256 ((__m256i *)&scrambleOut[iDx*6]);

            inter = _mm256_shuffle_epi8 (inter, q64Shift2a);
            inter = _mm256_permutevar8x32_epi32 (inter, q64Shift1);
            inter = _mm256_shuffle_epi8 (inter, q64Shift2b);
            
            deInterleaved = (wnInt8 *)&inter;

            memcpy (&bitSelectIn[iDx], &deInterleaved[0], 4);
            memcpy (&bitSelectIn[interval + iDx], &deInterleaved[4], 4);
            memcpy (&bitSelectIn[interval2 + iDx], &deInterleaved[8], 4);
            memcpy (&bitSelectIn[interval3 + iDx], &deInterleaved[12], 4);
            memcpy (&bitSelectIn[interval4 + iDx], &deInterleaved[16], 4);
            memcpy (&bitSelectIn[interval5 + iDx], &deInterleaved[20], 4);
        }

        if (excess)
        {
            inter = _mm256_lddqu_si256 ((__m256i *)&scrambleOut[iDx*6]);

            inter = _mm256_shuffle_epi8 (inter, q64Shift2a);
            inter = _mm256_permutevar8x32_epi32 (inter, q64Shift1);
            inter = _mm256_shuffle_epi8 (inter, q64Shift2b);
            
            deInterleaved = (wnInt8 *)&inter;

            memcpy (&bitSelectIn[iDx], &deInterleaved[0], excess);
            memcpy (&bitSelectIn[interval + iDx], &deInterleaved[4], excess);
            memcpy (&bitSelectIn[interval2 + iDx], &deInterleaved[8], excess);
            memcpy (&bitSelectIn[interval3 + iDx], &deInterleaved[12], excess);
            memcpy (&bitSelectIn[interval4 + iDx], &deInterleaved[16], excess);
            memcpy (&bitSelectIn[interval5 + iDx], &deInterleaved[20], excess);       
        }

        break;


        /************************** 256 QAM *********************/
        case 8:

        iDx2 = (rmLength>>5)<<5;
        iDx3 = iDx2>>3;

        interval = rmLength>>3;
        interval2 = interval<<1;
        interval3 = interval + interval2;
        interval4 = interval<<2;
        interval5 = interval + interval4;
        interval6 = interval + interval5;
        interval7 = interval + interval6;
        

        excess = (rmLength & 0x1f)>>3;

        __m256i q256Shift1 = _mm256_set_epi32(7,3,6,2,5,1,4,0);
        __m256i q256Shift2a = _mm256_set_epi8(15,7,14,6,13,5,12,4,11,3,10,2,9,1,8,0,15,7,14,6,13,5,12,4,11,3,10,2,9,1,8,0);
        __m256i q256Shift2b = _mm256_set_epi8(15,14,11,10,13,12,9,8,7,6,3,2,5,4,1,0,15,14,11,10,13,12,9,8,7,6,3,2,5,4,1,0);

        for (iDx = 0; iDx < iDx3 ; iDx=iDx+4)
        {
            inter = _mm256_lddqu_si256 ((__m256i *)&scrambleOut[iDx<<3]);

            inter = _mm256_shuffle_epi8 (inter, q256Shift2a);
            inter = _mm256_permutevar8x32_epi32 (inter, q256Shift1);
            inter = _mm256_shuffle_epi8 (inter, q256Shift2b);
            
            deInterleaved = (wnInt8 *)&inter;

            memcpy (&bitSelectIn[iDx], &deInterleaved[0], 4);
            memcpy (&bitSelectIn[interval + iDx], &deInterleaved[4], 4);
            memcpy (&bitSelectIn[interval2 + iDx], &deInterleaved[8], 4);
            memcpy (&bitSelectIn[interval3 + iDx], &deInterleaved[12], 4);
            memcpy (&bitSelectIn[interval4 + iDx], &deInterleaved[16], 4);
            memcpy (&bitSelectIn[interval5 + iDx], &deInterleaved[20], 4);
            memcpy (&bitSelectIn[interval6 + iDx], &deInterleaved[24], 4);
            memcpy (&bitSelectIn[interval7 + iDx], &deInterleaved[28], 4);
        }

        if (excess)
        {
            inter = _mm256_lddqu_si256 ((__m256i *)&scrambleOut[iDx<<3]);

            inter = _mm256_shuffle_epi8 (inter, q256Shift2a);
            inter = _mm256_permutevar8x32_epi32 (inter, q256Shift1);
            inter = _mm256_shuffle_epi8 (inter, q256Shift2b);
            
            deInterleaved = (wnInt8 *)&inter;

            memcpy (&bitSelectIn[iDx], &deInterleaved[0], excess);
            memcpy (&bitSelectIn[interval + iDx], &deInterleaved[4], excess);
            memcpy (&bitSelectIn[interval2 + iDx], &deInterleaved[8], excess);
            memcpy (&bitSelectIn[interval3 + iDx], &deInterleaved[12], excess);
            memcpy (&bitSelectIn[interval4 + iDx], &deInterleaved[16], excess);
            memcpy (&bitSelectIn[interval5 + iDx], &deInterleaved[20], excess);
            memcpy (&bitSelectIn[interval6 + iDx], &deInterleaved[24], excess);
            memcpy (&bitSelectIn[interval7 + iDx], &deInterleaved[28], excess);
            
        }

        break;
    }

    #ifdef DEBUG
        // FILE *fidroin;printf("::=%d\n",rmLength);
        // fidroin = fopen("RIn.txt","w");

        // for (wnUInt32 i = 0; i < rmLength; i++)
        // {
        //     fprintf(fidroin,"%d\n",bitSelectIn[i]);
        // }
        // fclose(fidroin);
    #endif

    wnUInt16   flag1=0,flag2 = 0, flag3= 0;
    wnUInt16 idx1 = 0;
    wnUInt16 idx2 = start_pos_bit_select;
    wnUInt16 offset;
    wnUInt16 factor;
    wnInt8 *deRateMatchOut, *llrCountPtr;
    __m256i ldpcInVector, bitSelectInVector,LlrCountVectorIn;//,LlrCountVectorOut;
    wnInt8 countSet = 1;
    wnUInt8 LlrCount[MAX_LDPC_DECODER_INPUT_LENGTH];
    __m256i setCount = _mm256_set1_epi8(countSet);

    while (idx1 < rmLength)
    {
        /// ******************************Section :- A *****************************************                                                                                               
        if (idx2 == start_pos_bit_select && (idx1 < rmLength))
        {
            if(start_pos_bit_select <fillerStart)
                factor = fillerStart;
            else
                factor= Ncb;

            if(flag1==0)
            {
                if (idx1 + factor -start_pos_bit_select > rmLength)
                    offset = rmLength - idx1;
                else
                    offset = factor -start_pos_bit_select;

                memcpy(&ldpcIn[idx2], &bitSelectIn[idx1], offset);
                memset(&LlrCount[idx2], 1, offset);
                flag1++;
                idx1 += factor -start_pos_bit_select;

                if(start_pos_bit_select <fillerStart)
                    idx2 = fillerStart;
                else
                    idx2 = 0;
            }
            else
            {
                while ((idx1 < rmLength) && (idx2 < factor) && (idx2 >= start_pos_bit_select))
                {
                    if((idx2+32<factor) && (idx1+32<rmLength))
                    {
                        ldpcInVector = _mm256_lddqu_si256((__m256i *)&ldpcIn[idx2]);
                        LlrCountVectorIn = _mm256_lddqu_si256((__m256i *)&LlrCount[idx2]);
                        bitSelectInVector = _mm256_lddqu_si256((__m256i *)&bitSelectIn[idx1]);
                        ldpcInVector = _mm256_add_epi8(ldpcInVector,bitSelectInVector);
                        LlrCountVectorIn = _mm256_add_epi8(LlrCountVectorIn,setCount);
                        deRateMatchOut = (wnInt8 *)&ldpcInVector;
                        llrCountPtr = (wnInt8 *)&LlrCountVectorIn;
                        memcpy(&ldpcIn[idx2],&deRateMatchOut[0],32);
                        memcpy(&LlrCount[idx2],&LlrCountVectorIn[0],32);
                        idx2+=32;
                        idx1+=32;
                    }
                    else
                    {
                        while((idx1 < rmLength) && (idx2 < factor))
                        {
                            LlrCount[idx2]++;
                            ldpcIn[idx2] = (ldpcIn[idx2] + bitSelectIn[idx1]);
                            idx2++;
                            idx1++; // */
                        }
                    }
                } // end of while
                if(start_pos_bit_select <fillerStart)
                    idx2 = fillerEnd;
                else
                    idx2 = 0;
            }// end of else
        } // end of if 

        /// ******************************Section :- A *****************************************     
        if ( idx2 == 0 && (idx1 < rmLength))
        {
            if(start_pos_bit_select <fillerStart)
                factor = start_pos_bit_select;
            else
                factor= fillerStart;


            if (flag2 == 0)
            {
                if (idx1 + factor > rmLength)
                    offset = rmLength - idx1;
                else
                    offset = factor;

                memcpy(&ldpcIn[idx2], &bitSelectIn[idx1], offset);
                memset(&LlrCount[idx2], 1, offset);
                flag2++;
                idx1 += offset;
                idx2 = factor;
            }//end of if
            else
            {                  
                while ((idx1 < rmLength) && (idx2 < factor) && (idx2 >= 0))
                {
                   
                    if((idx2+32<factor) && (idx1+32<rmLength))
                    {
                        ldpcInVector = _mm256_lddqu_si256((__m256i *)&ldpcIn[idx2]);
                        LlrCountVectorIn = _mm256_lddqu_si256((__m256i *)&LlrCount[idx2]);
                        bitSelectInVector = _mm256_lddqu_si256((__m256i *)&bitSelectIn[idx1]);
                        ldpcInVector = _mm256_add_epi8(ldpcInVector,bitSelectInVector);
                        LlrCountVectorIn = _mm256_add_epi8(LlrCountVectorIn,setCount);
                        deRateMatchOut = (wnInt8 *)&ldpcInVector;
                        llrCountPtr = (wnInt8 *)&LlrCountVectorIn;
                        memcpy(&ldpcIn[idx2],&deRateMatchOut[0],32);
                        memcpy(&LlrCount[idx2],&LlrCountVectorIn[0],32);
                        idx2+=32;
                        idx1+=32;
                    }
                    else{
                        while((idx1 < rmLength) && (idx2 < factor)){
                          LlrCount[idx2]++;
                          ldpcIn[idx2] = (ldpcIn[idx2] + bitSelectIn[idx1]);
                          idx2++;
                          idx1++; // */
                        }
                    }
                }// end of while
                if(start_pos_bit_select <fillerStart)
                    idx2 = start_pos_bit_select;
                else
                    idx2 = fillerStart;
            }//end of else
        }
        /// ******************************Section :- B *****************************************
        if (idx2 == fillerStart && (idx1 < rmLength))
        {
            memset(&ldpcIn[idx2], FILLER_BIT, fillerEnd-fillerStart);
            memset(&LlrCount[idx2], FILLER_BIT, fillerEnd-fillerStart);
            idx2 = fillerEnd;
        } //end of if
        /// ******************************Section :- C *****************************************
        if (idx2 == fillerEnd && (idx1 < rmLength))
        {
            if(start_pos_bit_select <fillerStart)
                factor = Ncb;
            else
                factor= start_pos_bit_select;


            if (flag3 == 0)
            {
                if (idx1 + factor -fillerEnd > rmLength)
                    offset = rmLength - idx1;
                else
                    offset = factor -fillerEnd;
                
                memcpy(&ldpcIn[idx2], &bitSelectIn[idx1], offset);
                memset(&LlrCount[idx2], 1, offset);
                flag3++;
                idx1 += offset;
                if(start_pos_bit_select <fillerStart)
                    idx2 = 0;
                else
                    idx2 = start_pos_bit_select;
            }//end of if
            else
            {
                while ((idx1 < rmLength) && (idx2 < factor) && (idx2 >= fillerEnd))
                {
                    if((idx2+32<factor) && (idx1+32<rmLength))
                    {
                        ldpcInVector = _mm256_lddqu_si256((__m256i *)&ldpcIn[idx2]);
                        LlrCountVectorIn = _mm256_lddqu_si256((__m256i *)&LlrCount[idx2]);
                        bitSelectInVector = _mm256_lddqu_si256((__m256i *)&bitSelectIn[idx1]);
                        ldpcInVector = _mm256_add_epi8(ldpcInVector,bitSelectInVector);
                        LlrCountVectorIn = _mm256_add_epi8(LlrCountVectorIn,setCount);
                        deRateMatchOut = (wnInt8 *)&ldpcInVector;
                        llrCountPtr = (wnInt8 *)&LlrCountVectorIn;
                        memcpy(&ldpcIn[idx2],&deRateMatchOut[0],32);
                        memcpy(&LlrCount[idx2],&LlrCountVectorIn[0],32);
                        idx2+=32;
                        idx1+=32;
                    }
                    else
                    {
                        while((idx1 < rmLength) && (idx2 < factor))
                        {
                            LlrCount[idx2]++;
                            ldpcIn[idx2] = (ldpcIn[idx2] + bitSelectIn[idx1]);
                            idx2++;
                            idx1++;
                        }
                    }
                }//end of while
                if(start_pos_bit_select <fillerStart)
                    idx2 = 0;
                else
                    idx2 = start_pos_bit_select;
            }//end of else
        }//end of if 
        /// **********************************************************************************
    }// end of while

    // If not using ICC Compiler
    #ifndef ICC
    if(flag1!=0 || flag2!=0 || flag3!=0)
    {
        idx1=0;
        while(idx1<Ncb)
        {
            if(ldpcIn[idx1]==FILLER_BIT)
                idx1++;
            else
            {
                if (LlrCount[idx1] > 1)
                    ldpcIn[idx1] = ldpcIn[idx1]/ LlrCount[idx1];
                idx1++;
            }
        }
    }
    #else
    /// ****************************** Averaging *****************************************
    if(Ncb-(fillerEnd-fillerStart)<rmLength)
    {
        idx1=0;
        wnInt8 temp1 = 1;
        __m256i setVector = _mm256_set1_epi8(temp1);
        while(idx1<Ncb)
        {
            if(idx1+32<Ncb){
               ldpcInVector = _mm256_lddqu_si256((__m256i *)&ldpcIn[idx1]);
               LlrCountVectorIn = _mm256_lddqu_si256((__m256 *)&LlrCount[idx1]);
               LlrCountVectorIn = _mm256_max_epi8(LlrCountVectorIn,setVector);
               ldpcInVector = _mm256_div_epi8(ldpcInVector,LlrCountVectorIn);
               deRateMatchOut = (wnInt8 *)&ldpcInVector;
               memcpy(&ldpcIn[idx1],&deRateMatchOut[0],32);
               idx1+=32;
            }
            else{
                ldpcIn[idx1] = ldpcIn[idx1]/LlrCount[idx1];
                idx1++;               
            }
        }
    }
    #endif

    #ifdef DEBUG
        /*FILE *fidrin;printf("::=%d\n",rmLength);
        fidrin = fopen("ROut.txt","w");

        for (wnUInt32 i = 0; i < Ncb; i++)
        {
            fprintf(fidrin,"%d\n",ldpcIn[i]);
        }
        fclose(fidrin);//*/
    #endif

    return;
}


// Descrambling
wnVoid wnNrPhyPuschDescrambling (puschConfigT* nrUlPhyPuschPDUs,\
wnInt8 *scramblingOut, wnInt8 *scramblingIn, wnUInt32 seqLength)
{
    wnInt32 iDx, iDx2 = 0;//n_id,
    // __attribute__ ((aligned(ALIGN_BOUNDARY))) wnInt8 convertedArray[32];

    // Scrambler initializer 
    wnInt64 Cinit = (nrUlPhyPuschPDUs->nRNTI<<15) + nrUlPhyPuschPDUs->nNid;

    wnUInt32 length = (seqLength+31)>>5;

    // PN Sequence
    wnNrPhyPuschPnSequence(Cinit, seqLength, pnSeqOut);

    __m256i* scrmblngIn = (__m256i*)&scramblingIn[0];
    __m256i* scrmblngOut = (__m256i*)&scramblingOut[0];
    __m256i convertedVector;
    
    // Descrambling Operation
    for (iDx = 0; iDx < length; iDx++)
    {
        #if 0
        memcpy(&convertedArray[0], &scramblerConvertor[pnSeqOut[iDx2++]], 8);
        memcpy(&convertedArray[8], &scramblerConvertor[pnSeqOut[iDx2++]], 8);
        memcpy(&convertedArray[16], &scramblerConvertor[pnSeqOut[iDx2++]], 8);
        memcpy(&convertedArray[24], &scramblerConvertor[pnSeqOut[iDx2++]], 8);
        convertedVector = _mm256_lddqu_si256 ((__m256i *)&convertedArray[0]);
        #else
        convertedVector = _mm256_setr_epi64x (*(int64_t *)&scramblerConvertor[pnSeqOut[iDx2]][0],
                                              *(int64_t *)&scramblerConvertor[pnSeqOut[iDx2+1]][0],
                                              *(int64_t *)&scramblerConvertor[pnSeqOut[iDx2+2]][0],
                                              *(int64_t *)&scramblerConvertor[pnSeqOut[iDx2+3]][0]);

        iDx2 += 4;
        #endif
        scrmblngOut[iDx] = _mm256_sign_epi8 (scrmblngIn[iDx],convertedVector);
    }

    #ifdef DEBUG
        // FILE *fidsout;
        // fidsout = fopen("SOut.txt","w");
        // for (wnUInt32 i = 0; i < seqLength; i++)
        // {
        //     fprintf (fidsout,"%d\n",scramblingOut[i]);
        // }
        // fclose (fidsout);
        // printf("Scrambling Length = %d\n",seqLength);
    #endif
    return;
}

// Used to find Rate DeMatch Length
wnInt32 wnNrPhyPuschTbsLbrmCalc(puschConfigT* nrUlPhyPuschPDUs)
{
    wnInt32 iDx;
    // Reference 5.4.2.1 Bit selection for LDPC Rate-matching in 38212-f11                                                                
             
    // Assuming that UE has not reported it's capability
    // maximum number of layers for one TB supported by the UE for the serving cell 
    wnInt32 maxLayerSuppByUE = nrUlPhyPuschPDUs->maxNumLayForOneTB;

    // Assuming that UE has not reported it's capability
    // maximum modulation order configured for the serving cell
    wnInt32 maxModSuppByUE = nrUlPhyPuschPDUs->maxModOrder;

    // Assuming that UE has not reported it's capability
    // maximum coding rate supported [Fixed Value]
    wnFlt maxCRSuppByUE = 0.9258;       //948/1024 = 0.9258

    // Assuming that UE has not reported it's capability
    // Number of PRBs configured for a UE
    // nPRB LBRM according to Table 5.4.2.1-1 in 38212-f11,
    //Maximum number of REs allocated
    wnInt32 maxREAllocated;
    if (nrUlPhyPuschPDUs->nBWPSize < 33)
        maxREAllocated = 4992;
    else if (nrUlPhyPuschPDUs->nBWPSize <= 66)
        maxREAllocated = 10296;
    else if (nrUlPhyPuschPDUs->nBWPSize <= 107)
        maxREAllocated = 16692;
    else if (nrUlPhyPuschPDUs->nBWPSize <= 135)
        maxREAllocated = 21060;
    else if (nrUlPhyPuschPDUs->nBWPSize <= 162)
        maxREAllocated = 25272;
    else if (nrUlPhyPuschPDUs->nBWPSize <= 217)
        maxREAllocated = 33852;
    else if (nrUlPhyPuschPDUs->nBWPSize > 217)
        maxREAllocated = 42588;

    // Reference 5.1.3.2 Transport block size determination for PDSCH in 38214-f10 

    //Intermediate number of information bits (Ninfo) is obtained
    wnFlt N_info = maxCRSuppByUE * maxREAllocated * maxLayerSuppByUE * maxModSuppByUE;

    wnUInt32 n,C,log2NInfo,N_desh_info,pdsch_tbs_lbrm;

    if (N_info <= 3824)
    {
        log2NInfo = log2(N_info);

        if (log2NInfo <= 9)
            n = 3;
        else
            n = log2NInfo - 6;

        if ((1<<n)*(wnInt32)(N_info/(1<<n)) <= 24)
            N_desh_info = 24;
        else
            N_desh_info = (1<<n)*(wnInt32)(N_info/(1<<n));
    
   
    
        // finding the closest TBS that is not less than N_desh_info         
        iDx = 0;         
        while (N_desh_info >= tbs_table[iDx])
            iDx++;
    
        // Final TBS size determination 
        pdsch_tbs_lbrm = tbs_table[iDx];
    }
    else
    {
        n = (wnInt32)log2(N_info-24) - 5;
    
        N_desh_info = ((int)round ((N_info-24)/(pow(2,13)))<<n);
   
        if (maxCRSuppByUE <= 0.25)
        {
            C = (N_desh_info+3839)/3816;
            pdsch_tbs_lbrm = (C<<3)*(((N_desh_info+23+(C<<3))/C)>>3) - 24;
        }
        
        else if (N_desh_info > 8424)
        {
            C = (N_desh_info+8447)/8424;
            pdsch_tbs_lbrm = (C<<3)*(((N_desh_info+23+(C<<3))/C)>>3) - 24;
        }
        else
            pdsch_tbs_lbrm = (((N_desh_info+31)>>3)<<3) - 24;
        
    }
 
    return pdsch_tbs_lbrm; 
}

// LDPC, Rate DeMatching Setup
wnVoid wnNrPhyPuschSetup (puschConfigT* nrUlPhyPuschPDUs, wnUInt32 *rmOutLength,\
wnUInt8 *m_bg, wnUInt8 *n_bg, wnUInt32 *eachCbDataLength, wnUInt8 *par_start,\
 wnUInt8 *Kb, wnUInt32 *Zc, wnUInt32 *n_parity, wnUInt32 *Pb, wnUInt32 *n_all,\
  wnUInt8 *L_max_int, wnUInt8 *Ltot_max_int, wnUInt8 *offset, wnUInt32 *baseGraphIndex,\
   wnUInt8 *maxItrs, wnUInt32 *nCodeBlocks, wnUInt32 *N_cb, wnUInt32 *start_pos_bit_select,\
    wnUInt32 *fillerStart, wnUInt32 *fillerEnd, wnUInt32 *puschRmLenCb)
{
    wnUInt32 crcLength, maxCodeBlkSize;
    wnUInt8 bg_no, bw_Ltot;

    if (nrUlPhyPuschPDUs->nTbSize < 3824)
        crcLength = 16;
    else
        crcLength = 24;

    // Maximum Code block Size
    if ((nrUlPhyPuschPDUs->nTbSize <= 292) || (nrUlPhyPuschPDUs->nTbSize <= 3824 && nrUlPhyPuschPDUs->codingRate <= 0.67) || (nrUlPhyPuschPDUs->codingRate <= 0.25))
    {
        maxCodeBlkSize = 3840;
    }
    else
    {
        maxCodeBlkSize = 8448;
    }

    // Number of Code Blocks
    *nCodeBlocks = (nrUlPhyPuschPDUs->nTbSize+maxCodeBlkSize-25)/(maxCodeBlkSize-24);

    *eachCbDataLength = ((nrUlPhyPuschPDUs->nTbSize + crcLength)/(*nCodeBlocks)) + 24*(*nCodeBlocks>1);                     

    /************ Selection of Base Graph *********************/
    if ((*eachCbDataLength<=292)||((*eachCbDataLength<=3824)&&(nrUlPhyPuschPDUs->codingRate<=0.67))||(nrUlPhyPuschPDUs->codingRate<=0.25))
    {
        bg_no = 2;                              //Type of Base Graph
        *m_bg = 42;                              //Rows of Base Matrix
        *n_bg = 52;                              //Columns of Base Matrix
        *par_start = 10;                         //Column Position from which Parity Starts
    }
    else
    {
        bg_no = 1;                              //Type of Base Graph
        *m_bg = 46;                              //Rows of Base Matrix
        *n_bg = 68;                              //Columns of Base Matrix
        *par_start = 22;                         //Column Position from which Parity Starts
    }

    /************ Imports Base Graph, Expansion Factor, Parity Start, Parity Shift Inverse ****************/
    if (bg_no == 1)
        *Kb = 22;
    else
        if (*eachCbDataLength > 640)
            *Kb = 10;
        else
            if (*eachCbDataLength > 560)
                *Kb = 9;
            else
                if (*eachCbDataLength > 192)
                    *Kb = 8;
                else
                    *Kb = 6;

    // Zc Selection
    wnUInt32 i = 0, iLS;
    
    while ((Zc_all[i]*(*Kb)) < *eachCbDataLength)
    {
        i = i + 1;
    }
    *Zc = Zc_all[i];
    iLS = iLS_all[i];
    
    // n_parity
    *n_parity = (*m_bg) * (*Zc);

    wnFlt codeRate;
    if ((nrUlPhyPuschPDUs->codingRate < (1.0f/3)) && (bg_no == 1))
        codeRate = 1.0f/3;

    else if ((nrUlPhyPuschPDUs->codingRate < 0.2) && (bg_no == 2))
        codeRate = 0.2;

    else
        codeRate = nrUlPhyPuschPDUs->codingRate;


    // Pb
    *Pb = ceil(((*eachCbDataLength)*(1/codeRate - 1)/(*Zc))) + 2; 

    // Total Length with Shortened Parity Bits
    *n_all = (*par_start + *Pb)*(*Zc);




    // L_max and Ltot_max values
    bw_Ltot = BW_L + 2; //bitwidth of Ltot
    *L_max_int = (1<<(BW_L-1));
    *Ltot_max_int = (1<<(bw_Ltot-1));

    // Offset
    *offset = nrUlPhyPuschPDUs->offset;

    // Selection of Base Graph
    *baseGraphIndex = (bg_no-1)*51+increment[iLS]+log10(Zc_all[i]/indexs[iLS])/log10(2);

    // max Iterations
    *maxItrs = MAX_ITERATIONS;

    wnInt32 ldpcDecoderInpLen, ldpcDecoderOutLen;

    if (bg_no == 1)
    {
        *fillerEnd = 20*(*Zc);
        ldpcDecoderInpLen = 66*(*Zc);
        ldpcDecoderOutLen = 22*(*Zc);
    }    
    else
    {
        *fillerEnd = (*Kb - 2)*(*Zc);
        ldpcDecoderInpLen = 50*(*Zc);
        ldpcDecoderOutLen = 10*(*Zc);
    }
    
    *fillerStart = *eachCbDataLength - 2*(*Zc);

    // wnUInt32 FillerLen = ldpcDecoderOutLen - *eachCbDataLength; 

    // FIXME : Number of CSIRS positions per PRB in PUSCH
    // CSIRS is not included || Default value = 0
    // wnInt32 no_of_csirs_per_prb = 0;
    
    // FIXME : These variables need to be removed hen DMRS is added
    nPrbAllocated = nrUlPhyPuschPDUs->nRBSize;
    nPtrsSyms = 0;
    
    // Summing up all the orthogonal CDM Groups used in the current transmission
    // "dmrsCdmG" gives whether that particular CDM Group is is being used for
    // transmission or not.
    wnInt32 orthoCdmGroup = dmrsCdmGWithoutData[0] + dmrsCdmGWithoutData[1] + dmrsCdmGWithoutData[2];

    // Number of REs allocated to NR PUSCH = Total REs for PUSCH - REs to leave for PUSCH and for orthogonal DMRS - REs to leave for PTRS
    wnInt32 no_of_REs_allocated = nPrbAllocated * 12 * nrUlPhyPuschPDUs->nNrOfSymbols - orthoCdmGroup * dmrs_pr_prb * nPrbAllocated - nPtrsSyms;//nPrbAllocated * 12 * (nrUlPhyPuschPDUs->nNrOfSymbols - dmrs_symbols);//
        
    // Consider Number of Layers and Modulation order used for PUSCH
    wnInt32 nr_pusch_ratematch_len = no_of_REs_allocated *  nrUlPhyPuschPDUs->nNrOfLayers * nrUlPhyPuschPDUs->modulationOrder;//nrUlPhyPuschPDUs->nNrOfLayers * nrUlPhyPuschPDUs->modulationOrder * 8 * ceil (nrUlPhyPuschPDUs->nTbSize/(nrUlPhyPuschPDUs->codingRate * nrUlPhyPuschPDUs->nNrOfLayers * nrUlPhyPuschPDUs->modulationOrder *8));
    //no_of_REs_allocated *  nrUlPhyPuschPDUs->nNrOfLayers * nrUlPhyPuschPDUs->modulationOrder;// +  nrUlPhyPuschPDUs->modOrder[cwordIndex]*5;//tbs_layers[cwordIndex] * nrUlPhyPuschPDUs->modOrder[cwordIndex]*ceil(nrUlPhyPuschPDUs->nTbSize[cwordIndex]/(nrUlPhyPuschPDUs->codingRate[cwordIndex]*  tbs_layers[cwordIndex] * nrUlPhyPuschPDUs->modOrder[cwordIndex]));//
    *rmOutLength = nr_pusch_ratematch_len;

    // I_LBRM is fixed in NR PUSCH Rate matching specifications
    wnFlt I_LBRM;
    wnUInt8 LbrmFbrmSelect = nrUlPhyPuschPDUs->lmtBuffRm;
    if (LbrmFbrmSelect  == 1)
        I_LBRM = 1;
    else
        I_LBRM = 0;
  

    // R_LBRM is fixed in NR PUSCH Rate matching specifications
    wnFlt R_LBRM = 0.66667;
    
    // Should be calculated according to 
    // 5.1.3.2 Transport block size determination in TS 38.214
    wnInt32 pusch_tbs_lbrm = wnNrPhyPuschTbsLbrmCalc(nrUlPhyPuschPDUs);

    if (I_LBRM == 0)
        *N_cb = ldpcDecoderInpLen;
    
    else if (I_LBRM == 1)
    {
        *N_cb = pusch_tbs_lbrm/((*nCodeBlocks)*R_LBRM);

        if (rmOutLength < (wnUInt32)*N_cb)
            *N_cb = ldpcDecoderInpLen;    
    }

    // Calculating Rate Match Lenths for each Code Block
    wnInt32 iDx1, iDx2, iDx3, iDx4, cbIdx;
    iDx3 = nrUlPhyPuschPDUs->nNrOfLayers*nrUlPhyPuschPDUs->modulationOrder;
    iDx4 = iDx3*(*nCodeBlocks);
    iDx1 = iDx3*(int)((nr_pusch_ratematch_len + iDx4 -1)/iDx4);
    iDx2 = iDx3*(int)(nr_pusch_ratematch_len/iDx4);

    iDx4 = *nCodeBlocks - ((nr_pusch_ratematch_len/iDx3)%(*nCodeBlocks));

    for (cbIdx = 0; cbIdx < iDx4; cbIdx++)
    {
        puschRmLenCb[cbIdx] = iDx2;
    }
    for (cbIdx = iDx4; cbIdx < *nCodeBlocks; cbIdx++)
    {
        puschRmLenCb[cbIdx] = iDx1;
    }

    // if( *Kb == 22 )
    // {
    //     *Pb = ceil( (puschRmLenCb[0] - 20*(*Zc) + FillerLen)/(wnFlt)(*Zc));
    //     *Pb = *Pb>4?*Pb:4;
    //     *Pb = *Pb>46?46:*Pb;
    // }
    // else
    // {
    //     *Pb = ceil( (puschRmLenCb[0] - 8*(*Zc) + FillerLen)/(wnFlt)(*Zc));
    //     *Pb = *Pb>4?*Pb:4;
    //     *Pb = *Pb>42?42:*Pb;
    // }
    
    // *n_all = (*par_start + *Pb)*(*Zc);

    // Starting Position for Rate Matching
    // if redundancy version is 0
    if (nrUlPhyPuschPDUs->nRV == 0)
    {
        if (bg_no == 1)
            *start_pos_bit_select = 0;
        else
            *start_pos_bit_select = 0;
    }
    // if redundancy version is 1
    else if (nrUlPhyPuschPDUs->nRV == 1)
    {
        if (bg_no == 1)
            *start_pos_bit_select = *Zc * ((17*(*N_cb))/(66*(*Zc)));
        else
            *start_pos_bit_select = *Zc * ((13*(*N_cb))/(50*(*Zc)));
    }
    // if redundancy version is 2
    else if (nrUlPhyPuschPDUs->nRV == 2)
    {
        if (bg_no == 1)
            *start_pos_bit_select = *Zc * ((33*(*N_cb))/(66*(*Zc)));
        else
            *start_pos_bit_select = *Zc * ((25*(*N_cb))/(50*(*Zc)));
    }
    // if redundancy version is 3
    else
    {
        if (bg_no == 1)
            *start_pos_bit_select = *Zc * ((56*(*N_cb))/(66*(*Zc)));
        else
            *start_pos_bit_select = *Zc * ((43*(*N_cb))/(50*(*Zc)));
    }    

    return;
}

wnVoid wnNrPhyPuschLdpcDecoderSetup (wnUInt32 *row_wt, wnUInt32 *col_wt, wnUInt32 *r_nz, wnUInt32 *c_nz, wnUInt32 *sh_nz, wnUInt8 m_bg, wnUInt8 n_bg,  const wnInt32 H_bg[][n_bg])
{
    wnUInt32 i, j, c = 0;

    // Rows and Columns where there are -1s in Base Matrix
    for (i = 0; i < m_bg; i++)
    {
        for (j = 0; j < n_bg; j++)
        {
            if (H_bg[i][j] != -1)
            {
                row_wt[i] += 1;
                col_wt[j] += 1;
                sh_nz[c] = H_bg[i][j];
                c_nz[c] = j;
                r_nz[c] = i;
                c += 1;
            }
        }
    }


    return;
}

/************************ Function Used in Normal C *****************************/
wnVoid mul_shifted_identity_decoder(wnInt8 *in, wnInt8 *out, wnUInt32 Zc, wnUInt32 sh)
{
    wnUInt32 i;
    wnInt8 *b;
    b = (wnInt8 *)malloc(Zc*sizeof(wnInt8));
    if (sh==-1)
        for(i=0;i<Zc;i++)
            out[i] = 0;
    else
    {
        for (i=sh;i<Zc;i++)
            b[i-sh] = in[i];
        for (i=0;i<sh;i++)
            b[Zc-sh+i] = in[i];
        for (i=0;i<Zc;i++)
            out[i] = b[i];
    }
    free(b);
}


// Operations on each Code Word
wnVoid wnNrPhyPuschCwProc (puschConfigT* nrUlPhyPuschPDUs,\
                            wnInt8 *deScramIn, wnUInt32 rmOutLength)
{
    /****************************************************************************************************
     * m_bg         = Number of rows in LDPC Base Graph (BG-1 : 46, BG-2 : 42)                          *
     * n_bg         = Number of Columns in LDPC Base Graph (BG-1 : 68, BG-2 : 52)                       *
     * par_start    = Parity Start Column in the Base Graph (BG-1 : 23, BG-2 : 11)                      *
     * Kb           = Effective Message Blocks for the LDPC configuration                               *
     * offset       = Used as the offset input for the Offset Min-Sum LDPC Decoder Algorithm            *
     * maxItrs      = Maximum Iterations for the LDPC Decoder Module(Acts as Booundary Condition)       *
     * L_max_int    = Maximum LLR for the LDPC Processing                                               *
     * Ltot_max_int = Maximum LLR Sum Posiible for the LDPC Processing                                  *
     * check        = Output of CRC Module (1 ==> CRC Passed, 0 ==> CRC Failed)                         *
     * Zc           = Expansion Factor of LDPC Base Graph                                               *
     * msgLength    = Message Length of each Code Blocks                                                *
     * nCodeBlocks  = Number of Code Blocks for the given configuration                                 *
     * N_cb         = Code Block Size after Rate Matching                                               *
     * fillerStart  = Start Position of Filler Bits in Code Block                                       *
     * fillerEnd    = End Position of Filler Bits in Code Block                                         * 
     * start_pos_bit_select = Start Position for Soft combining, based on RV                            *
     * n_parity     = Number of Parity Bits                                                             *
     * baseGraphIndex = Used for selecting LDPC Base Graph                                              *
     * Pb           = Number of Parity Blocks                                                           *
     * n_all        = Length of LDPC Decoder Input                                                      *
     * row_wt       = Number of non -1's in each row                                                    *
     * col_wt       = Number of non -1's in each column                                                 *
     * r_nz         = Row Position of all non -1's in Base Graph                                        *
     * c_nz         = Column Position of all non -1's in Base Graph                                     *
     * sh_nz        = Shift by all non -1's in Base Graph                                               *
     * puschRmLenCb = Rate Match Length of each Code Block                                              *
     ****************************************************************************************************/

    wnUInt8 m_bg, n_bg, par_start, Kb, offset, maxItrs, L_max_int, Ltot_max_int;

    wnUInt32 Zc, msgLength, nCodeBlocks, N_cb, fillerStart, fillerEnd, start_pos_bit_select, n_parity,
             cbPosIdx = 0, msgIdx = 0, baseGraphIndex, Pb, n_all, row_wt[46] = {0}, col_wt[68] = {0}, 
             r_nz[316] = {0}, c_nz[316] = {0}, sh_nz[316] = {0}, puschRmLenCb[MAX_PUSCH_CB], length, bytePackIdx;

    __m256i *bytePackIn, bytePackInter;

    /********************** LDPC Decoder and RateMatching Setup Functions **************************/
    wnNrPhyPuschSetup (nrUlPhyPuschPDUs, &rmOutLength/*LookIntoThis*/, &m_bg, &n_bg, &msgLength, &par_start,
                         &Kb, &Zc, &n_parity, &Pb, &n_all, &L_max_int, &Ltot_max_int, &offset, &baseGraphIndex,
                         &maxItrs, &nCodeBlocks, &N_cb, &start_pos_bit_select, &fillerStart, &fillerEnd, &puschRmLenCb[0]);

    wnNrPhyPuschLdpcDecoderSetup (row_wt, col_wt, r_nz, c_nz, sh_nz, m_bg, n_bg, &baseGraph[baseGraphIndex][0]);
    
    #ifdef DEBUG
        printf("m_bg = %d,n_bg = %d, msgLength = %d, par_start = %d, Kb = %d, Zc = %d, n_parity = %d,Pb = %d, n_all = %d, L_max_int = %d, Ltot_max_int = %d, offset = %d, baseGraphIndex = %d,maxItrs = %d, nCodeBlocks = %d, N_Cb = %d, start_pos_bit_select = %d, fillerStart = %d,filler End = %d\n", m_bg, n_bg, msgLength, par_start, Kb, Zc, n_parity, Pb, n_all, L_max_int,Ltot_max_int, offset, baseGraphIndex, maxItrs,nCodeBlocks, N_cb, start_pos_bit_select, fillerStart, fillerEnd);

        FILE *fRmOut = fopen ("IOs/puschRmOut.txt","w");
        FILE *fTbOut = fopen ("IOs/puschTbOut.txt","w");
        FILE *fScramIn = fopen ("IOs/puschDeScIn.txt","w");
        for (wnUInt32 i = 0; i < rmOutLength; i++)
        {
            fprintf (fScramIn, "%d\n", deScramIn[i]);
        }
        fclose (fScramIn);
    #endif

    /********************** DeScrambling ******************************************/
    wnNrPhyPuschDescrambling (nrUlPhyPuschPDUs, scramblingOut, deScramIn, rmOutLength);

    #ifdef DEBUG    
    FILE *fLdpcIn = fopen("IOs/puschLdpcInp.txt", "w");
    FILE* fLdpcOut = fopen("IOs/puschLdpcOut.txt", "w");
    FILE *fDeScOut = fopen("IOs/puschDeScOut.txt", "w");
    for(int idx = 0 ;idx<rmOutLength;idx++)
    {
        fprintf(fDeScOut, "%hhd\n", scramblingOut[idx]);
    }
    fclose(fDeScOut);
    #endif

    printf("finished Descrambling\n");
    /********************** Code Word Processing **********************************/
    for (wnUInt32 cbIdx = 0; cbIdx <  nCodeBlocks; cbIdx++)
    {
        /********************** DeInterleaver and Rate Dematcher ***********************/
        wnNrPhyPuschRateDematching (&scramblingOut[cbPosIdx], &ldpcIn[2*Zc], N_cb, start_pos_bit_select, fillerStart,
                                     fillerEnd, puschRmLenCb[cbIdx], nrUlPhyPuschPDUs->modulationOrder);
        cbPosIdx += puschRmLenCb[cbIdx];
        
        #ifdef DEBUG
            for (wnUInt32 i = 2*Zc; i < N_cb+2*Zc; i++)
            {
                fprintf (fRmOut,"%d\n",ldpcIn[i]);
            }
            printf ("rmLen = %d, mod_order = %d\n",puschRmLenCb[cbIdx],nrUlPhyPuschPDUs->modulationOrder);
            
        #endif  

        // Setting the Filler Bits
        memset(&ldpcIn[fillerStart+2*Zc], 31, fillerEnd - fillerStart);

        #ifdef DEBUG
            for(int idx = 0;idx<68*Zc;idx++)
            {
                fprintf(fLdpcIn, "%d\n", ldpcIn[idx]);
            }
        #endif

        /********************** LDPC Decoder *****************************************/
        nrldpc_decoder(ldpcIn, Ltot, Lreg, &ldpcOut[0], row_wt, c_nz, sh_nz, Zc, msgLength, Kb,\
                         Pb, par_start, n_parity, n_all, maxItrs, offset, L_max_int, Ltot_max_int);
                       

        // Data being packed for crc
        length = (msgLength+31)>>5;
        bytePackIn = (__m256i *)&ldpcOut[0];
        for (bytePackIdx = 0; bytePackIdx < length; bytePackIdx++)
        {
            bytePackInter = _mm256_slli_epi32 (bytePackIn[bytePackIdx],7);
            crcIn[bytePackIdx] = _mm256_movemask_epi8 (bytePackInter);
        }

        #ifdef DEBUG
            for (wnUInt32 i = 0; i < msgLength; i++)
            {
                fprintf(fLdpcOut,"%d\n",ldpcOut[i]);
            }
        #endif

        if (nCodeBlocks > 1)
        {
            /********************** CRC_24B after Decoding and CB Attachment *************/
            wnNrPhyCrc24BRemoval ((wnUInt8 *)&crcIn[0], msgLength, &puschData[msgIdx>>3], &check);
            msgIdx += msgLength - 24;
            #ifdef DEBUG
            printf("check: %d\n", check);
            #endif
        }
        else
        {
            msgIdx += msgLength;
        }
    }

    #ifdef DEBUG
        fclose(fLdpcOut);
        fclose(fLdpcIn);
    #endif

    

    /********************** CRC After Reconstructing CB Segments *****************/
    if (nrUlPhyPuschPDUs->nTbSize < 3824)
    {
        wnNrPhyCrc16Removal (ldpcOut, msgIdx, puschDecodedData, &check);
    }
    else
    {
        if (nrUlPhyPuschPDUs->nTbSize > 8424)
            wnNrPhyCrc24ARemoval (puschData, msgIdx, puschDecodedData, &check);
        else
            wnNrPhyCrc24ARemoval ((wnUInt8 *)&crcIn, msgIdx, puschDecodedData, &check);
    }
           
    
    #ifdef DEBUG
        for (wnUInt32 i = 0; i < msgIdx>>3; i++)
        {
            for (wnUInt32 j = 0; j < 8; j++)
                fprintf (fTbOut,"%d\n",(puschDecodedData[i]>>j)&1);
        }
        fclose (fTbOut);
        fclose (fRmOut);
    #endif   
    return;
}
    

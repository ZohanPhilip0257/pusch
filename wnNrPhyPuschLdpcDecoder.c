#include "wnNrPhyPuschHeader.h"
#include <stdlib.h>
#include <string.h>

// Check Parity Function
wnUInt32 check_parity(wnInt8 *dec_cword, wnUInt32 Kb, wnUInt32 Pb, wnUInt32 Zc, wnUInt32 par_start, wnUInt32 *row_wt, wnUInt32 *c_nz, wnUInt32 *sh_nz, wnInt8 Lreg[20][384],wnUInt32 nparity_rem_Zc)
{
    #if (AVX_2)     
    __attribute__ ((aligned(ALIGN_BOUNDARY))) wnInt8 Ltot_rearranged[7680] = {0};
    __attribute__ ((aligned(ALIGN_BOUNDARY))) wnInt8 L_rearranged[7680] = {0};  

    wnUInt32 multiple_32 = (Zc+31)>>5;
    __m256i *vLtot, *vL, check;
    wnInt8 *checkingEqual;
    
    wnUInt32 pctonok = 0, i = 0, c = 0, j, k, kmax;
    
     
    while ((i < Pb) && (pctonok == 0))
    {
        memcpy(Ltot_rearranged, dec_cword + c_nz[c]*Zc + sh_nz[c], Zc - sh_nz[c]);
        memcpy(Ltot_rearranged + Zc - sh_nz[c], dec_cword + c_nz[c]*Zc, sh_nz[c]);
        
        vLtot = (__m256i*) &Ltot_rearranged[0];
        c++;
        
        for (j = 1; j < row_wt[i]; j++)
        {
            if((c_nz[c] < Kb)||((c_nz[c] >= par_start) && (c_nz[c] < par_start + Pb)))
            {
                memcpy(L_rearranged, dec_cword + c_nz[c]*Zc + sh_nz[c], Zc - sh_nz[c]);
                memcpy(L_rearranged + Zc - sh_nz[c], dec_cword + c_nz[c]*Zc, sh_nz[c]);
                
                vL = (__m256i*) &L_rearranged[0];
                
                for ( k = 0; k < multiple_32; k++)
                {
                    *vLtot = _mm256_xor_si256( *vLtot, *vL);
                    vLtot++;
                    vL++;
                }
                vLtot = vLtot - multiple_32;
                vL = vL - multiple_32;
            }

            c++;
        }
        
        if ((i == Pb-1) && (nparity_rem_Zc > 0))
        {
            kmax = (nparity_rem_Zc+31)>>5;
        }
        else
        {
            kmax = multiple_32;
        }
        
        k = 0;
        while ((pctonok == 0) && (k < kmax))
        {
            checkingEqual = (wnInt8*)vLtot;
            
            j = 0;
            while((pctonok == 0) && (j < 32))
            {
                if (checkingEqual[j] == 1)
                    pctonok = 1;
                
                j++;
            }
            k++;
            vLtot++;
        }
        vLtot = vLtot - k;
        
        if (pctonok == 0)
        {
            i++;
        }
    }

    #else
    wnUInt32 i,j,k,kmax,c,pctonok = 0;  
    i = 0; c = 0;

    while ((pctonok==0)&(i<Pb))
        {
            for (k=0;k<Zc;k++)
                    Lreg[0][k]=0;

            for (j=0;j<row_wt[i];j++)
            {
                if ((c_nz[c]<Kb)|((c_nz[c]>=par_start)&(c_nz[c]<(Pb+par_start))))
                    {   
                         mul_shifted_identity_decoder(&dec_cword[c_nz[c]*Zc], Lreg[1], Zc, sh_nz[c]);
                    
                 for (k=0;k<Zc;k++)
                             Lreg[0][k] ^= Lreg[1][k];
                    }

                    c += 1;
            }
            
        if (i==(Pb-1))  //(&& nparity_rem_Zc > 0)
                    kmax = nparity_rem_Zc;
            else
                    kmax = Zc;
            k=0;
            while ((pctonok==0)&(k<kmax))
            {
                    if (Lreg[0][k]==1)
                        pctonok = 1;        
            k += 1;
            
            }
            if (pctonok==0)
                    i += 1;
        }
    #endif
    return pctonok;
}

#if 1
typedef struct _llrArrays
{
    __attribute__ ((aligned(ALIGN_BOUNDARY))) wnInt8 Lreg[20][384];
    __attribute__ ((aligned(ALIGN_BOUNDARY))) wnInt8 Ltot_rearranged[7680];
} llrArrays;


// LDPC Decoder Function
wnVoid nrldpc_decoder(wnInt8 *L0, wnInt8 *Ltot, wnInt8 Lreg1[20][384], wnInt8 *dec_cword, wnUInt32 *row_wt, wnUInt32 *c_nz, wnUInt32 *sh_nz, wnUInt32 Zc, wnUInt32 msg_len, wnUInt8 Kb, wnUInt32 Pb, wnUInt8 par_start, wnUInt32 nparity, wnUInt32 n_all, wnUInt8 maxitrs, wnUInt8 offset, wnUInt8 L_max, wnUInt8 Ltot_max)
{
    wnInt32 i,itrs_taken,j,k,k1,s1,c,c1, tempIndx;
    wnInt32 c_ind[20], col_ind[20], sh[20];

    wnInt32 multiple_32 = ((Zc+31) >> 5) << 5;
    wnInt8 nVecs = ((Zc+31) >> 5);

    #ifdef DEBUG
    printf("Vectors : %d\n", nVecs);
    #endif

//     mexPrintf("New Func\n");
    llrArrays llrs = {0};
    __attribute__ ((aligned(ALIGN_BOUNDARY))) wnInt8 L[316][384] = {0};
    // __attribute__ ((aligned(ALIGN_BOUNDARY))) wnInt8 Lreg[20][384] = {0};
    // printf("Zc = %d, Pb = %d, Kb = %d, par_start = %d, nparity = %d, n_all = %d, offset = %d, L_max = %d, Ltot_max = %d\n",Zc, Pb, Kb, par_start, nparity, n_all,offset,n_all,msg_len);
    #if (AVX_2)
    // __attribute__ ((aligned(ALIGN_BOUNDARY))) wnInt8 Ltot_rearranged[7680] = {0}; // max(row weight for BG1, row weight for BG2)*max(Zc) = 20*384 = 7680
    wnInt8 *temp1, *temp2;

    __m256i *curr_ptr, curr, currVec, signVec, v_min1, v_minus_min1,\
            v_min2, sign, sign_min, sign_to_be_used, curr_sign[20],\
            mask, mask1, v_min_pos, v_second_operand, to_store, *v_Ltot_rearranged_ptr,\
            *v_L_ptr, mask2;

    wnInt8 *temp3, *temp4;
    
    const __m256i all_ff    = _mm256_set1_epi8(0xff);
    const __m256i all_m1s   = _mm256_set1_epi8 (-1);
    const __m256i all_ones  = _mm256_set1_epi8(1);
    const __m256i all_zeros = _mm256_set1_epi8(0);
    const __m256i v_L_max_minus_one = _mm256_set1_epi8((wnInt8)(L_max-1));
    const __m256i v_minus_L_max = _mm256_set1_epi8((wnInt8)(-L_max));
    const __m256i v_offset_minus_one = _mm256_set1_epi8((wnInt8)(offset - 1));
    const __m256i v_offset = _mm256_set1_epi8((wnInt8)offset);
    // __m128i *vLtotR128 = (__m128i *)v_Ltot_rearranged_ptr;
    // __m128i *vL128 = (__m128i *)v_L_ptr;
    // __m256i inter[4];
    __m256i currVec_arr[12], signVec_arr[12], v_min1_arr[12];
    #else
    wnInt32 sprod, min1, min2, min1pos, t1;
    wnInt32 sgns[20];
    __m256i mask;
    const __m256i all_ones  = _mm256_set1_epi8(1);
    const __m256i all_zeros = _mm256_set1_epi8(0);
    const __m256i v_offset = _mm256_set1_epi8((wnInt8)offset);
    #endif

    memcpy (Ltot, L0, n_all);
    // _mm_prefetch (&c_nz[0], 256);
    // _mm_prefetch (row_wt, 32);

    wnInt32 itr = 0;
    wnInt32 pcnotok = 1;
    wnInt8 flag = 0;
    while((pcnotok==1)&(itr<maxitrs))
    {
        // row operation
        c = 0;
        for (i=0;i<Pb;i++)
        {
            // struct timespec tmf; clock_gettime(CLOCK_REALTIME, &tmf); 
            // stRTime[rNum] = tmf.tv_sec * G_HZ + tmf.tv_nsec;

            // read one block row and row align
            #if (AVX_2)
            c1 = 0;
            // _mm_prefetch (&llrs, sizeof (llrArrays));
            // _mm_prefetch (Ltot + c_nz[c]*384, multiple_32);
            // _mm_prefetch (L + c*384, row_wt[i]*384);
            for (j = 0; j < row_wt[i]; j++)
            {
                if ((c_nz[c]<Kb)|((c_nz[c]>=par_start)&(c_nz[c]<(Pb+par_start))))
                {
                    s1 = c_nz[c]*Zc;
                    k1 = sh_nz[c];
                    // _mm_prefetch (Ltot + c_nz[c+1]*384, 384);
                    memcpy (llrs.Ltot_rearranged + c1*multiple_32, Ltot + s1 + k1, Zc - k1);
                    memcpy (llrs.Ltot_rearranged + c1*multiple_32 + Zc - k1, Ltot + s1, k1);

                    v_Ltot_rearranged_ptr = (__m256i*) &llrs.Ltot_rearranged[0] + c1*nVecs;  // 384/32 = 12
                    v_L_ptr = (__m256i*) &L[0][0] + c*12;                    // set pointer to beginning of correct row of L matrix
                    
                    for (k = 0; k < multiple_32; k += 32)
                    {
                        *v_Ltot_rearranged_ptr = _mm256_subs_epi8 (*v_Ltot_rearranged_ptr, *v_L_ptr); // signed 8 bit saturated subtration

                        // mask = _mm256_cmpgt_epi8 (*v_Ltot_rearranged_ptr, v_L_max_minus_one);
                        // *v_L_ptr = _mm256_blendv_epi8 (*v_Ltot_rearranged_ptr, v_L_max_minus_one, mask);

                        // mask = _mm256_cmpgt_epi8 (v_minus_L_max, *v_Ltot_rearranged_ptr);
                        // *v_L_ptr = _mm256_blendv_epi8 (*v_L_ptr, v_minus_L_max, mask);

                        *v_L_ptr = _mm256_max_epi8 (_mm256_min_epi8 (*v_Ltot_rearranged_ptr, v_L_max_minus_one), v_minus_L_max);

                        v_Ltot_rearranged_ptr++;
                        v_L_ptr++;
                    } // end of k loop

                    c_ind[c1] = c;
                    col_ind[c1] = c_nz[c];
                    sh[c1] = sh_nz[c];
                    c1 += 1;
                }
                c += 1;
            } // end of j loop
            
            #else
            c1 = 0;
            for (j=0;j<row_wt[i];j++)
            {
                if ((c_nz[c]<Kb)|((c_nz[c]>=par_start)&(c_nz[c]<(Pb+par_start))))
                {
                    s1 = c_nz[c]*Zc;
                    k1 = sh_nz[c];
    
                    for (k=0;k<Zc;k++)
                    {
                        t1 = Ltot[s1+k1] - L[c][k];
                            if (t1>(Ltot_max-1))
                                Ltot[s1+k1] = Ltot_max - 1;
                            else
                                if (t1<-Ltot_max)
                                    Ltot[s1+k1] = -Ltot_max;
                                else
                                    Ltot[s1+k1] = t1;
                       
                        if (t1>(L_max-1)) //Clip to L_max for L
                            t1 = L_max - 1;
                        if (t1<-L_max)
                            t1 = -L_max;
                        Lreg[c1][k] = t1;

                        k1 = k1 + 1;
                        if (k1 == Zc)
                            k1 = 0;
                    }
                    c_ind[c1] = c;
                    col_ind[c1] = c_nz[c];
                    sh[c1] = sh_nz[c];
                    c1 += 1;
                }
                c += 1;
            }
            #endif
            // struct timespec tmg; clock_gettime(CLOCK_REALTIME, &tmg); 
            // enRTime[rNum++] = tmg.tv_sec * G_HZ + tmg.tv_nsec;

            // struct timespec tmh; clock_gettime(CLOCK_REALTIME, &tmh); 
            // stMTime[dNum] = tmh.tv_sec * G_HZ + tmh.tv_nsec;

            // minsum operation: find min1, min2, min1pos, sign
            #if (AVX_2)
            #if 0
            for (k = 0; k < multiple_32; k += 32)
            {
                v_min1 = _mm256_load_si256((__m256i*) &L[c_ind[0]][k]);
                curr_sign[0] = _mm256_cmpgt_epi8(all_zeros, v_min1);
                sign = curr_sign[0];
                sign_min =  curr_sign[0]; 
                v_min1 = _mm256_abs_epi8(v_min1);
                v_min2 = _mm256_set1_epi8(127);
                v_min_pos = all_zeros;
    
                for (j = 1; j < c1; j++)
                {
                    curr_ptr = (__m256i*) &L[c_ind[j]][k];
                    curr_sign[j] = _mm256_cmpgt_epi8(all_zeros, *curr_ptr);
                    sign = _mm256_xor_si256(sign, curr_sign[j]);
                    curr = _mm256_abs_epi8(*curr_ptr);
                    
                    mask1     = _mm256_cmpgt_epi8(v_min1, curr);
                    v_min2    = _mm256_blendv_epi8(v_min2, v_min1, mask1);
                    v_min1    = _mm256_blendv_epi8(v_min1, curr, mask1);  
                    v_min_pos = _mm256_blendv_epi8(v_min_pos, _mm256_set1_epi8(j), mask1);
        
                    sign_min = _mm256_or_si256( _mm256_and_si256(sign_min,_mm256_xor_si256(all_ff,mask1)),_mm256_and_si256(mask1, curr_sign[j]));

                    mask1  = _mm256_xor_si256(all_ff, mask1);   
                    mask   = _mm256_cmpgt_epi8(v_min2, curr);
                    mask   = _mm256_and_si256(mask, mask1);
                    v_min2 = _mm256_blendv_epi8(v_min2, curr, mask);
                    
                } // end of j loop
                                   
                mask  = _mm256_cmpgt_epi8(v_min1, v_offset);
                v_second_operand = _mm256_blendv_epi8(all_zeros, v_offset, mask);
                v_min1 = _mm256_sub_epi8(v_min1, v_second_operand);
                v_minus_min1 = _mm256_sub_epi8(all_zeros, v_min1);
                               
                mask  = _mm256_cmpgt_epi8(v_min2, v_offset);
                v_second_operand = _mm256_blendv_epi8(all_zeros, v_offset, mask);
                v_min2 = _mm256_sub_epi8(v_min2, v_second_operand);
                
                
                for (j = 0; j < c1; j++)
                {
                    sign_to_be_used = _mm256_xor_si256(sign, curr_sign[j]);
                    to_store = _mm256_blendv_epi8 (v_min1,v_minus_min1,sign_to_be_used);
                    _mm256_store_si256((__m256i*) &L[c_ind[j]][k],to_store);
                } // end of j loop
                
                temp1 = (wnInt8*) &v_min_pos;
                temp2 = (wnInt8*) &v_min2;
                sign_min = _mm256_xor_si256(sign,sign_min);
                temp3 = (wnInt8*) &sign_min;
                for (j = 0; j < 32; j++)
                {
                    L[c_ind[temp1[j]]][k+j] = ((temp3[j] == 0) ? temp2[j] : -temp2[j]);
                } // end of j loop      
            } // end of k loop
            #else
                #if 1
                    #if 1
                    for (k = 0; k < multiple_32; k += 32)
                    {
                        // 1st Vector
                        currVec = *(__m256i *) &L[c_ind[1]][k];
                        v_min1 = _mm256_abs_epi8 (currVec);
                        signVec = _mm256_blendv_epi8 (all_ones, all_ff, currVec);
                        
                        for (j = 2; j < c1; j++)
                        {
                            currVec = *(__m256i *) &L[c_ind[j]][k];
                            v_min1  = _mm256_min_epu8 (v_min1, _mm256_abs_epi8 (currVec));
                            signVec = _mm256_sign_epi8 (signVec, _mm256_blendv_epi8 (all_ones, all_ff, currVec)); 
                        }
                        
                        #if 1
                        mask2 = _mm256_cmpgt_epi8 (v_min1, v_offset);
                        mask2 = _mm256_blendv_epi8 (all_zeros, v_offset, mask2);
                        *((__m256i *) &llrs.Lreg[0][k]) = _mm256_sign_epi8(_mm256_sub_epi8(v_min1, mask2), signVec);
                        #else
                        *((__m256i *) &llrs.Lreg[0][k]) = _mm256_sign_epi8 (v_min1, signVec);
                        #endif

                        // Stores the row parity
                        currVec = *(__m256i *) &L[c_ind[0]][k];
                        mask = _mm256_sign_epi8 (signVec, _mm256_blendv_epi8 (all_ones, all_ff, currVec));

                        // 1st Column data
                        mask1 = _mm256_abs_epi8(currVec);

                        // Remaining Vectors
                        for (wnUInt8 currColIdx = 1; currColIdx < c1; currColIdx++)
                        {
                            v_min1 = mask1;

                            for (j = 1; j < c1; j++)
                                if (j != currColIdx)
                                    v_min1  = _mm256_min_epu8(v_min1, _mm256_abs_epi8(*(__m256i *) &L[c_ind[j]][k]));
                            
                            signVec = _mm256_sign_epi8 (mask, _mm256_blendv_epi8 (all_ones, all_ff, *(__m256i *) &L[c_ind[currColIdx]][k]));
                            
                            #if 1
//                             mexPrintf("Testing New \n");
                            mask2 = _mm256_cmpgt_epi8 (v_min1, v_offset);
                            mask2 = _mm256_blendv_epi8 (all_zeros, v_offset, mask2);
                            *((__m256i *) &llrs.Lreg[currColIdx][k]) = _mm256_sign_epi8(_mm256_sub_epi8(v_min1, mask2), signVec);
                            #else
                            *((__m256i *) &llrs.Lreg[currColIdx][k]) = _mm256_sign_epi8(v_min1, signVec);
                            #endif
                        }
                    }
                    #else
                    _mm_prefetch (&llrs.Lreg[1][0], 384);
                    _mm_prefetch (&L[c_ind[0]][0], 384*c1);
                    for (k = 0; k < multiple_32; k += 32)
                    {
                        wnUInt8 startColIdx = 1;
                        
                        for (wnUInt8 currColIdx = 0; currColIdx < c1; currColIdx++)
                        {
                            if (currColIdx < 2)
                            {
                                currVec = *(__m256i *) &L[c_ind[startColIdx]][k];
                                mask = _mm256_sign_epi8 (all_ones, currVec);
                                mask1 = _mm256_abs_epi8(currVec);
                            }

                            v_min1 = mask1;
                            signVec = mask;

                            for (j = 1; j < c1; j++)
                            {
                                if ((j != currColIdx) && (j != startColIdx))
                                {
                                    currVec = *(__m256i *) &L[c_ind[j]][k];
                                    v_min1  = _mm256_min_epu8(v_min1, _mm256_abs_epi8(currVec));
                                    signVec  = _mm256_sign_epi8(signVec, currVec);
                                }
                            }
                            *((__m256i *) &llrs.Lreg[currColIdx][k]) = _mm256_sign_epi8(v_min1, signVec);

                            startColIdx = 0;
                        }
                    }
                    #endif
                #else
                _mm_prefetch (currVec_arr, multiple_32);
                _mm_prefetch (signVec_arr, multiple_32);
                _mm_prefetch (v_min1_arr, multiple_32);
                _mm_prefetch (&L[c_ind[0]], 384);
                _mm_prefetch (&L[c_ind[1]], 384);
                _mm_prefetch (llrs.Lreg, 2*384);

                wnUInt8 startColIdx = 1;
                for (wnUInt8 currColIdx = 0; currColIdx < c1; currColIdx++)
                {
                    if (currColIdx < 2)
                    {
                        for (k = 0; k < nVecs; k++)
                        {
                            currVec_arr[k] = *(__m256i *) &L[c_ind[startColIdx]][k<<5];
                            signVec_arr[k] = _mm256_sign_epi8 (all_ones, currVec_arr[k]);
                            v_min1_arr[k]  = _mm256_abs_epi8(currVec_arr[k]);
                        }
                    }

                    for (j = 1; j < c1; j++)
                    {
                        if ((j != currColIdx)&&(j != startColIdx))
                        {
                            for (k = 0; k < nVecs; k++)
                            {
                                currVec_arr[k] = *(__m256i *) &L[c_ind[j]][k<<5];
                                v_min1_arr[k]  = _mm256_min_epu8(v_min1_arr[k], _mm256_abs_epi8(currVec_arr[k]));
                                signVec_arr[k]  = _mm256_sign_epi8(signVec_arr[k], currVec_arr[k]);
                            }
                        }
                    }
                    // _mm_prefetch (&L[c_ind[currColIdx+1]], 384);
                    #if 1
                    for (k = 0; k < nVecs; k++)
                        *((__m256i *) &llrs.Lreg[currColIdx][k<<5]) = _mm256_sign_epi8(v_min1_arr[k], signVec_arr[k]);
                    #else
                    mask = _mm256_cmpgt_epi8 (v_offset, v_min1);
                    mask = _mm256_blendv_epi8 (v_offset, all_zeros, mask);
                    *((__m256i *) &llrs.Lreg[currColIdx][k]) = _mm256_sign_epi8(_mm256_sub_epi8(v_min1, mask), signVec);
                    #endif

                    startColIdx = 0;
                }
                #endif
            #endif

            #else
            #if 0
            for (k=0;k<Zc; k++)
            {
                sprod = 1; min1 = L_max; min2 = L_max; min1pos = 0;
                for (j=0;j<c1;j++)
                {
                    sgns[j] = 1-2*(llrs.Lreg[j][k]<0);
                    sprod = sprod * sgns[j];
                    t1 = sgns[j] * llrs.Lreg[j][k];
                    
                    //if (t1<=min1)
                    if (t1<min1)
                    {
                        min2 = min1; min1 = t1; min1pos = j;
                    }
                    else
                        if (t1<min2) //if (t1<=min2)
                            min2 = t1;
                }
                if (min1 >= offset)         
                    min1 = min1 - offset;
                if (min2 >= offset)
                    min2 = min2 - offset;


                for (j=0;j<c1;j++)
                    if (j==min1pos)
                        llrs.Lreg[j][k] = min2*sgns[j]*sprod;
                    else
                        llrs.Lreg[j][k] = min1*sgns[j]*sprod;
            }
            #else
            __m256i v_min11;
            for (k = 0; k < multiple_32; k += 32)
            {
                int startColIdx = 1;
                for (int currColIdx = 0; currColIdx < c1; currColIdx++)
                {
                    __m256i currVec = *(__m256i *) &llrs.Lreg[startColIdx][k];
                    __m256i signVec = _mm256_sign_epi8 (all_ones, currVec);
                    v_min11  = _mm256_abs_epi8(currVec);

                    for (j = 1; (j < c1); j++)
                    {
                        if ((j != currColIdx)&&(j != startColIdx))
                        {
                            currVec = *(__m256i *) &llrs.Lreg[j][k];
                            v_min11  = _mm256_min_epu8(v_min11, _mm256_abs_epi8(currVec));
                            signVec  = _mm256_sign_epi8(signVec, currVec);
                        }
                    }
                    #if 0
                    *((__m256i *) &Lreg1[currColIdx][k]) = _mm256_sign_epi8(v_min11, signVec);
                    #else
                    mask = _mm256_cmpgt_epi8 (v_min11, v_offset);
                    mask = _mm256_blendv_epi8 (all_zeros, v_offset, mask);
                    *((__m256i *) &Lreg1[currColIdx][k]) = _mm256_sign_epi8(_mm256_sub_epi8(v_min11, mask), signVec);
                    #endif

                    startColIdx = 0;
                }
            }
            memcpy (&llrs.Lreg[0][0], &Lreg1[0][0], 20*384);
            #endif
            #endif
            // struct timespec tmj; clock_gettime(CLOCK_REALTIME, &tmj); 
            // enMTime[dNum++] = tmj.tv_sec * G_HZ + tmj.tv_nsec;

            // struct timespec tmk; clock_gettime(CLOCK_REALTIME, &tmk); 
            // stCTime[cNum] = tmk.tv_sec * G_HZ + tmk.tv_nsec;
           
            // column operation
            #if (AVX_2)
            for (j = 0; j < c1; j++)
            {
                v_Ltot_rearranged_ptr = (__m256i*) &llrs.Ltot_rearranged[0] + j*nVecs;  // 384/32 = 12
                #if 0
                v_L_ptr = (__m256i*) &L[c_ind[j]][0];                          // set pointer to beginning of correct row of L matrix
                #else
                v_L_ptr = (__m256i *) &llrs.Lreg[j][0];
                curr_ptr = (__m256i *) &L[c_ind[j]][0];
                #endif

                for (k = 0; k < multiple_32; k += 32)
                {
                    *curr_ptr = *v_L_ptr;
                    *v_Ltot_rearranged_ptr = _mm256_adds_epi8 (*v_Ltot_rearranged_ptr, *v_L_ptr);
                    v_Ltot_rearranged_ptr++;
                    v_L_ptr++;
                    curr_ptr++;
                } // end of k loop

                s1 = col_ind[j]*Zc;
                k1 = sh[j];
                memcpy(Ltot + s1 + k1, llrs.Ltot_rearranged + j*multiple_32, Zc - k1);
                memcpy(Ltot + s1, llrs.Ltot_rearranged + j*multiple_32 + Zc - k1, k1);

                // _mm_prefetch (&L[c_ind[j]][0], multiple_32);
            } // end of j loop
            
            #else
            for (j=0;j<c1;j++)
            {
                s1 = col_ind[j]*Zc;
                k1 = (sh[j] == 0) ? 0 : (Zc-sh[j]);
                for (k=0;k<Zc;k++)
                {   
                    L[c_ind[j]][k] = llrs.Lreg[j][k];

                    t1 = Ltot[s1+k] + llrs.Lreg[j][k1];
                    if (t1 > (Ltot_max-1))
                        t1 = Ltot_max-1;
                    else
                        if (t1 < -Ltot_max)
                            t1 = -Ltot_max;

                    Ltot[s1+k] = t1;
                    k1 = k1 + 1;
                    if (k1 == Zc)
                        k1 = 0;
                }
            }

            #endif
            // struct timespec tml; clock_gettime(CLOCK_REALTIME, &tml); 
            // enCTime[cNum++] = tml.tv_sec * G_HZ + tml.tv_nsec;
        }

        //make decisions
            
        #if (CHECK_PARITY)
            for (i=0;i<msg_len;i++)
                dec_cword[i] = (Ltot[i]<0);
            for (i=0;i<nparity;i++)
                dec_cword[par_start*Zc+i] = (Ltot[par_start*Zc+i]<0);
           
            //check for parity
            pcnotok = check_parity(dec_cword, Kb, Pb, Zc, par_start, row_wt, c_nz, sh_nz, llrs.Lreg, 0);//n_parity_rem_Zc); Replaced n_parity_rem_ZC by 0;
            itrs_taken = itr+1;
            //printf("%d,%d\n",itrs_taken, pcnotok);
            if (pcnotok==1) 
              itr += 1;
        #else
              itr += 1;
        #endif
    }
       // printf("\nitrs_taken = %d\n",itr);
    if (CHECK_PARITY == 0)
    {
        for (i=0;i<msg_len;i++)
            dec_cword[i] = (Ltot[i]<0);
        for (i=0;i<nparity;i++)
            dec_cword[par_start*Zc+i] = (Ltot[par_start*Zc+i]<0);
    }


    return;
}

#else
// LDPC Decoder Function
wnVoid nrldpc_decoder(wnInt8 *L0, wnInt8 *Ltot, wnInt8 Lreg[20][384], wnInt8 *dec_cword, wnUInt32 *row_wt, wnUInt32 *c_nz, wnUInt32 *sh_nz, wnUInt32 Zc, wnUInt32 msg_len, wnUInt8 Kb, wnUInt32 Pb, wnUInt8 par_start, wnUInt32 nparity, wnUInt32 n_all, wnUInt8 maxitrs, wnUInt8 offset, wnUInt8 L_max, wnUInt8 Ltot_max)
{
    wnInt32 i,itrs_taken,j,k,k1,s1,c,c1, tempIndx;
    wnInt32 c_ind[20], col_ind[20], sh[20];

    __attribute__ ((aligned(ALIGN_BOUNDARY))) wnInt8 L[316][384] = {0};
    // printf("Zc = %d, Pb = %d, Kb = %d, par_start = %d, nparity = %d, n_all = %d, offset = %d, L_max = %d, Ltot_max = %d\n",Zc, Pb, Kb, par_start, nparity, n_all,offset,n_all,msg_len);
    #if (AVX_2)
    __attribute__ ((aligned(ALIGN_BOUNDARY))) wnInt8 Ltot_rearranged[7680] = {0}; // max(row weight for BG1, row weight for BG2)*max(Zc) = 20*384 = 7680
    wnInt8 *temp1, *temp2;

    __m256i *curr_ptr, curr, v_min1, v_minus_min1,  v_min2, sign, sign_min, sign_to_be_used, curr_sign[20], mask, mask1, v_min_pos, v_second_operand, to_store, *v_Ltot_rearranged_ptr, *v_L_ptr;

    wnInt8 *temp3, *temp4;
    const __m256i all_ff    = _mm256_set1_epi8(0xff);
	const __m256i all_ones  = _mm256_set1_epi8(1);
	const __m256i all_zeros = _mm256_set1_epi8(0);
    const __m256i v_L_max_minus_one = _mm256_set1_epi8((wnInt8)(L_max-1));
    const __m256i v_minus_L_max = _mm256_set1_epi8((wnInt8)(-L_max));
    const __m256i v_offset_minus_one = _mm256_set1_epi8((wnInt8)(offset - 1));
	const __m256i v_offset = _mm256_set1_epi8((wnInt8)offset);
    
    wnInt32 multiple_32 = ((Zc+31) >> 5) << 5;

    #else
    wnInt32 sprod, min1, min2, min1pos, t1;
    wnInt32 sgns[20];
    #endif

    memcpy(Ltot, L0, n_all);
    
    wnInt32 itr = 0;
    wnInt32 pcnotok = 1;
    while((pcnotok==1)&(itr<maxitrs))
    {
        //row operation
        c = 0;
        for (i=0;i<Pb;i++)
        {
            //read one block row and row align
            #if (AVX_2)
            c1 = 0;
            for (j=0;j<row_wt[i];j++)
            {
                if ((c_nz[c]<Kb)|((c_nz[c]>=par_start)&(c_nz[c]<(Pb+par_start))))
                {
                    s1 = c_nz[c]*Zc;
                    k1 = sh_nz[c];

                    memcpy(Ltot_rearranged + c1*384, Ltot + s1 + k1, Zc - k1);
                    memcpy(Ltot_rearranged + c1*384 + Zc - k1, Ltot + s1, k1);
			
		            v_Ltot_rearranged_ptr = (__m256i*) &Ltot_rearranged[0] + c1*12;  // 384/32 = 12
                    v_L_ptr = (__m256i*) &L[0][0] + c*12;                    // set pointer to beginning of correct row of L matrix
                    
                    for (k = 0; k < multiple_32; k += 32)
                    {
                        *v_Ltot_rearranged_ptr = _mm256_subs_epi8(*v_Ltot_rearranged_ptr, *v_L_ptr); // signed 8 bit saturated subtration

                        mask = _mm256_cmpgt_epi8(*v_Ltot_rearranged_ptr, v_L_max_minus_one);
                        *v_L_ptr = _mm256_blendv_epi8(*v_Ltot_rearranged_ptr, v_L_max_minus_one, mask);

                        mask = _mm256_cmpgt_epi8(v_minus_L_max, *v_Ltot_rearranged_ptr);
                        *v_L_ptr = _mm256_blendv_epi8(*v_L_ptr, v_minus_L_max, mask);

                        v_Ltot_rearranged_ptr++;
                        v_L_ptr++;
                        
                    } // end of k loop

                    c_ind[c1] = c;
                    col_ind[c1] = c_nz[c];
                    sh[c1] = sh_nz[c];
                    c1 += 1;
                }
                c += 1;
            } // end of j loop
                  
            #else
            c1 = 0;
            for (j=0;j<row_wt[i];j++)
            {
                if ((c_nz[c]<Kb)|((c_nz[c]>=par_start)&(c_nz[c]<(Pb+par_start))))
                {
                    s1 = c_nz[c]*Zc;
                    k1 = sh_nz[c];
    
                    for (k=0;k<Zc;k++)
                    {
                        t1 = Ltot[s1+k1] - L[c][k];
                            if (t1>(Ltot_max-1))
                                Ltot[s1+k1] = Ltot_max - 1;
                            else
                                if (t1<-Ltot_max)
                                    Ltot[s1+k1] = -Ltot_max;
                                else
                                    Ltot[s1+k1] = t1;
                       
                        if (t1>(L_max-1)) //Clip to L_max for L
                            t1 = L_max - 1;
                        if (t1<-L_max)
                            t1 = -L_max;
                        Lreg[c1][k] = t1;

                        k1 = k1 + 1;
                        if (k1 == Zc)
                            k1 = 0;
                    }
                    c_ind[c1] = c;
                    col_ind[c1] = c_nz[c];
                    sh[c1] = sh_nz[c];
                    c1 += 1;
                }
                c += 1;
            }
            #endif
            

            //minsum operation: find min1, min2, min1pos, sign
            
            #if (AVX_2)
            for (k = 0; k < multiple_32; k += 32)
            {
            	v_min1 = _mm256_load_si256((__m256i*) &L[c_ind[0]][k]);
				curr_sign[0] = _mm256_cmpgt_epi8(all_zeros, v_min1);
				sign = curr_sign[0];
				sign_min =  curr_sign[0]; 
				v_min1 = _mm256_abs_epi8(v_min1);
				v_min2 = _mm256_set1_epi8(127);
				v_min_pos = all_zeros;
	
            	for (j = 1; j < c1; j++)
            	{
                    curr_ptr = (__m256i*) &L[c_ind[j]][k];
					curr_sign[j] = _mm256_cmpgt_epi8(all_zeros, *curr_ptr);
					sign = _mm256_xor_si256(sign, curr_sign[j]);
					curr = _mm256_abs_epi8(*curr_ptr);
					
					mask1     = _mm256_cmpgt_epi8(v_min1, curr);
					v_min2    = _mm256_blendv_epi8(v_min2, v_min1, mask1);
					v_min1    = _mm256_blendv_epi8(v_min1, curr, mask1);  
					v_min_pos = _mm256_blendv_epi8(v_min_pos, _mm256_set1_epi8(j), mask1);
		
					sign_min = _mm256_or_si256( _mm256_and_si256(sign_min,_mm256_xor_si256(all_ff,mask1)),_mm256_and_si256(mask1, curr_sign[j]));

					mask1  = _mm256_xor_si256(all_ff, mask1);   
					mask   = _mm256_cmpgt_epi8(v_min2, curr);
					mask   = _mm256_and_si256(mask, mask1);
					v_min2 = _mm256_blendv_epi8(v_min2, curr, mask);
					
            	} // end of j loop
            	                   
			    mask  = _mm256_cmpgt_epi8(v_min1, v_offset);
				v_second_operand = _mm256_blendv_epi8(all_zeros, v_offset, mask);
				v_min1 = _mm256_sub_epi8(v_min1, v_second_operand);
				v_minus_min1 = _mm256_sub_epi8(all_zeros, v_min1);
			                   
				 mask  = _mm256_cmpgt_epi8(v_min2, v_offset);
				v_second_operand = _mm256_blendv_epi8(all_zeros, v_offset, mask);
				v_min2 = _mm256_sub_epi8(v_min2, v_second_operand);
				
				
				for (j = 0; j < c1; j++)
            	{
            	    sign_to_be_used = _mm256_xor_si256(sign, curr_sign[j]);
					to_store = _mm256_blendv_epi8 (v_min1,v_minus_min1,sign_to_be_used);
            		_mm256_store_si256((__m256i*) &L[c_ind[j]][k],to_store);
            	} // end of j loop
            	
            	temp1 = (wnInt8*) &v_min_pos;
            	temp2 = (wnInt8*) &v_min2;
		        sign_min = _mm256_xor_si256(sign,sign_min);
		        temp3 = (wnInt8*) &sign_min;
            	for (j = 0; j < 32; j++)
            	{
			        L[c_ind[temp1[j]]][k+j] = ((temp3[j] == 0) ? temp2[j] : -temp2[j]);
            	} // end of j loop  	
            	
            } // end of k loop

            #else
            for (k=0;k<Zc; k++)
            {
                sprod = 1; min1 = L_max; min2 = L_max; min1pos = 0;
                for (j=0;j<c1;j++)
                {
                    sgns[j] = 1-2*(Lreg[j][k]<0);
                    sprod = sprod * sgns[j];
                    t1 = sgns[j] * Lreg[j][k];
                    
                    //if (t1<=min1)
                    if (t1<min1)
                    {
                        min2 = min1; min1 = t1; min1pos = j;
                    }
                    else
                        if (t1<min2) //if (t1<=min2)
                            min2 = t1;
                }
                if (min1 >= offset)         
                    min1 = min1 - offset;
                if (min2 >= offset)
                    min2 = min2 - offset;


                for (j=0;j<c1;j++)
                    if (j==min1pos)
                        Lreg[j][k] = min2*sgns[j]*sprod;
                    else
                        Lreg[j][k] = min1*sgns[j]*sprod;
            }
            #endif

            
            //column operation
            #if (AVX_2)
            for (j = 0; j < c1; j++)
            {
                v_Ltot_rearranged_ptr = (__m256i*) &Ltot_rearranged[0] + j*12;  // 384/32 = 12
		        v_L_ptr = (__m256i*) &L[c_ind[j]][0];                          // set pointer to beginning of correct row of L matrix

                for (k = 0; k < multiple_32; k += 32)
                {
                    *v_Ltot_rearranged_ptr = _mm256_adds_epi8(*v_Ltot_rearranged_ptr, *v_L_ptr);
                    v_Ltot_rearranged_ptr++;
		            v_L_ptr++;
                } // end of k loop

                s1 = col_ind[j]*Zc;
                k1 = sh[j];
                memcpy(Ltot + s1 + k1, Ltot_rearranged + j*384, Zc - k1);
                memcpy(Ltot + s1, Ltot_rearranged + j*384 + Zc - k1, k1);

            } // end of j loop
            
            #else
            for (j=0;j<c1;j++)
            {
                s1 = col_ind[j]*Zc;
                k1 = (sh[j] == 0) ? 0 : (Zc-sh[j]);
                for (k=0;k<Zc;k++)
                {   
                    L[c_ind[j]][k] = Lreg[j][k];

                    t1 = Ltot[s1+k] + Lreg[j][k1];
                    if (t1 > (Ltot_max-1))
                        t1 = Ltot_max-1;
                    else
                        if (t1 < -Ltot_max)
                            t1 = -Ltot_max;

                    Ltot[s1+k] = t1;
                    k1 = k1 + 1;
                    if (k1 == Zc)
                        k1 = 0;
                }
            }
            #endif

        }

        //make decisions
        
    #if (CHECK_PARITY)
        for (i=0;i<msg_len;i++)
            dec_cword[i] = (Ltot[i]<0);
        for (i=0;i<nparity;i++)
            dec_cword[par_start*Zc+i] = (Ltot[par_start*Zc+i]<0);
       
        //check for parity
        pcnotok = check_parity(dec_cword, Kb, Pb, Zc, par_start, row_wt, c_nz, sh_nz, Lreg, 0);//n_parity_rem_Zc); Replaced n_parity_rem_ZC by 0;
        itrs_taken = itr+1;
//printf("%d,%d\n",itrs_taken, pcnotok);
        if (pcnotok==1) 
          itr += 1;
    #else
          itr += 1;
    #endif

    }
//    printf("\nitrs_taken = %d\n",itr);
    if (CHECK_PARITY == 0)
    {
	    for (i=0;i<msg_len;i++)
            dec_cword[i] = (Ltot[i]<0);
        for (i=0;i<nparity;i++)
            dec_cword[par_start*Zc+i] = (Ltot[par_start*Zc+i]<0);
    }

return;
}

#endif












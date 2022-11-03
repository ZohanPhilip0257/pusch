#include "wnNrPhyPuschHeader.h"
#include<math.h>

wnFlt wnSinRad(wnFlt inp);
wnFlt wnCosRad(wnFlt inp);
wnUInt16 maxPrimesOfMzc[275] = {5,11,17,23,29,31,41,47,53,59,61,
                                71,73,83,89,89,101,107,113,113,
                                113,131,137,139,149,151,157,167,
                                173,179,181,191,197,199,199,211,
                                211,227,233,239,241,251,257,263,
                                269,271,281,283,293,293,293,311,
                                317,317,317,331,337,347,353,359,
                                359,367,373,383,389,389,401,401,
                                409,419,421,431,433,443,449,449,
                                461,467,467,479,479,491,491,503,
                                509,509,521,523,523,523,541,547,
                                557,563,569,571,577,587,593,599,
                                601,607,617,619,619,631,641,647,
                                653,659,661,661,677,683,683,691,
                                701,701,709,719,719,727,733,743,
                                743,751,761,761,773,773,773,787,
                                797,797,809,811,821,827,829,839,
                                839,839,857,863,863,863,881,887,
                                887,887,887,911,911,919,929,929,
                                941,947,953,953,953,971,977,983,
                                983,991,997,997,1013,1019,1021,
                                1031,1033,1039,1049,1051,1061,
                                1063,1069,1069,1069,1091,1097,
                                1103,1109,1109,1117,1123,1129,
                                1129,1129,1151,1153,1163,1163,
                                1171,1181,1187,1193,1193,1201,
                                1201,1217,1223,1229,1231,1237,
                                1237,1249,1259,1259,1259,1277,
                                1283,1289,1291,1301,1307,1307,
                                1319,1321,1327,1327,1327,1327,
                                1327,1361,1367,1373,1373,1381,
                                1381,1381,1399,1409,1409,1409,
                                1427,1433,1439,1439,1451,1453,
                                1459,1459,1471,1481,1487,1493,
                                1499,1499,1511,1511,1523,1523,
                                1531,1531,1543,1553,1559,1559,
                                1571,1571,1583,1583,1583,1601,
                                1607,1613,1619,1621,1627,1637,
                                1637,1637};

wnInt8 pi_onePRB[30][6] = {{-3, -1,  3,  3, -1, -3},
                           {-3,  3, -1, -1,  3, -3},
                           {-3, -3, -3,  3,  1, -3},
                           { 1,  1,  1,  3, -1, -3},
                           { 1,  1,  1, -3, -1,  3},
                           {-3,  1, -1, -3, -3, -3},
                           {-3,  1,  3, -3, -3, -3},
                           {-3, -1,  1, -3,  1, -1},
                           {-3, -1, -3,  1, -3, -3},
                           {-3, -3,  1, -3,  3, -3},
                           {-3,  1,  3,  1, -3, -3},
                           {-3, -1, -3,  1,  1, -3},
                           { 1,  1,  3, -1, -3,  3},
                           { 1,  1,  3,  3, -1,  3},
                           { 1,  1,  1, -3,  3, -1},
                           { 1,  1,  1, -1,  3, -3},
                           {-3, -1, -1, -1,  3, -1},
                           {-3, -3, -1,  1, -1, -3},
                           {-3, -3, -3,  1, -3, -1},
                           {-3,  1,  1, -3, -1, -3},
                           {-3,  3, -3,  1,  1, -3},
                           {-3,  1, -3, -3, -3, -1},
                           { 1,  1, -3,  3,  1,  3},
                           { 1,  1, -3, -3,  1, -3},
                           { 1,  1,  3, -1,  3,  3},
                           { 1,  1, -3,  1,  3,  3},
                           { 1,  1, -1, -1,  3, -1},
                           { 1,  1, -1,  3, -1, -1},
                           { 1,  1, -1,  3, -3, -1},
                           { 1,  1, -3,  1, -1, -1}};

wnInt8 pi_twoPRBs[30][12] = {{-3,  1, -3, -3, -3,  3, -3, -1,  1,  1,  1, -3},
                             {-3,  3,  1, -3,  1,  3, -1, -1,  1,  3,  3,  3},
                             {-3,  3,  3,  1, -3,  3, -1,  1,  3, -3,  3, -3},
                             {-3, -3, -1,  3,  3,  3, -3,  3, -3,  1, -1, -3},
                             {-3, -1, -1,  1,  3,  1,  1, -1,  1, -1, -3,  1},
                             {-3, -3,  3,  1, -3, -3, -3, -1,  3, -1,  1,  3},
                             { 1, -1,  3, -1, -1, -1, -3, -1,  1,  1,  1, -3},
                             {-1, -3,  3, -1, -3, -3, -3, -1,  1, -1,  1, -3},
                             {-3, -1,  3,  1, -3, -1, -3,  3,  1,  3,  3,  1},
                             {-3, -1, -1, -3, -3, -1, -3,  3,  1,  3, -1, -3},
                             {-3,  3, -3,  3,  3, -3, -1, -1,  3,  3,  1, -3},
                             {-3, -1, -3, -1, -1, -3,  3,  3, -1, -1,  1, -3},
                             {-3, -1,  3, -3, -3, -1, -3,  1, -1, -3,  3,  3},
                             {-3,  1, -1, -1,  3,  3, -3, -1, -1, -3, -1, -3},
                             { 1,  3, -3,  1,  3,  3,  3,  1, -1,  1, -1,  3},
                             {-3,  1,  3, -1, -1, -3, -3, -1, -1,  3,  1, -3},
                             {-1, -1, -1, -1,  1, -3, -1,  3,  3, -1, -3,  1},
                             {-1,  1,  1, -1,  1,  3,  3, -1, -1, -3,  1, -3},
                             {-3,  1,  3,  3, -1, -1, -3,  3,  3, -3,  3, -3},
                             {-3, -3,  3, -3, -1,  3,  3,  3, -1, -3,  1, -3},
                             { 3,  1,  3,  1,  3, -3, -1,  1,  3,  1, -1, -3},
                             {-3,  3,  1,  3, -3,  1,  1,  1,  1,  3, -3,  3},
                             {-3,  3,  3,  3, -1, -3, -3, -1, -3,  1,  3, -3},
                             { 3, -1, -3,  3, -3, -1,  3,  3,  3, -3, -1, -3},
                             {-3, -1,  1, -3,  1,  3,  3,  3, -1, -3,  3,  3},
                             {-3,  3,  1, -1,  3,  3, -3,  1, -1,  1, -1,  1},
                             {-1,  1,  3, -3,  1, -1,  1, -1, -1, -3,  1, -1},
                             {-3, -3,  3,  3,  3, -3, -1,  1, -3,  3,  1, -3},
                             { 1, -1,  3,  1,  1, -1, -1, -1,  1,  3, -3,  1},
                             {-3,  3, -3,  3, -3, -3,  3, -1, -1,  1,  3, -3}};

wnInt8 pi_threePRBs[30][18] = {{-1,  3, -1, -3,  3,  1, -3, -1,  3, -3, -1, -1,  1,  1,  1, -1, -1, -1},
                               { 3, -3,  3, -1,  1,  3, -3, -1, -3, -3, -1, -3,  3,  1, -1,  3, -3,  3},
                               {-3,  3,  1, -1, -1,  3, -3, -1,  1,  1,  1,  1,  1, -1,  3, -1, -3, -1},
                               {-3, -3,  3,  3,  3,  1, -3,  1,  3,  3,  1, -3, -3,  3, -1, -3, -1,  1},
                               { 1,  1, -1, -1, -3, -1,  1, -3, -3, -3,  1, -3, -1, -1,  1, -1,  3,  1},
                               { 3, -3,  1,  1,  3, -1,  1, -1, -1, -3,  1,  1, -1,  3,  3, -3,  3, -1},
                               {-3,  3, -1,  1,  3,  1, -3, -1,  1,  1, -3,  1,  3,  3, -1, -3, -3, -3},
                               { 1,  1, -3,  3,  3,  1,  3, -3,  3, -1,  1,  1, -1,  1, -3, -3, -1,  3},
                               {-3,  1, -3, -3,  1, -3, -3,  3,  1, -3, -1, -3, -3, -3, -1,  1,  1,  3},
                               { 3, -1,  3,  1, -3, -3, -1,  1, -3, -3,  3,  3,  3,  1,  3, -3,  3, -3},
                               {-3, -3, -3,  1, -3,  3,  1,  1,  3, -3, -3,  1,  3, -1,  3, -3, -3,  3},
                               {-3, -3,  3,  3,  3, -1, -1, -3, -1, -1, -1,  3,  1, -3, -3, -1,  3, -1},
                               {-3, -1, -3, -3,  1,  1, -1, -3, -1, -3, -1, -1,  3,  3, -1,  3,  1,  3},
                               { 1,  1, -3, -3, -3, -3,  1,  3, -3,  3,  3,  1, -3, -1,  3, -1, -3,  1},
                               {-3,  3, -1, -3, -1, -3,  1,  1, -3, -3, -1, -1,  3, -3,  1,  3,  1,  1},
                               { 3,  1, -3,  1, -3,  3,  3, -1, -3, -3, -1, -3, -3,  3, -3, -1,  1,  3},
                               {-3, -1, -3, -1, -3,  1,  3, -3, -1,  3,  3,  3,  1, -1, -3,  3, -1, -3},
                               {-3, -1,  3,  3, -1,  3, -1, -3, -1,  1, -1, -3, -1, -1, -1,  3,  3,  1},
                               {-3,  1, -3, -1, -1,  3,  1, -3, -3, -3, -1, -3, -3,  1,  1,  1, -1, -1},
                               { 3,  3,  3, -3, -1, -3, -1,  3, -1,  1, -1, -3,  1, -3, -3, -1,  3,  3},
                               {-3,  1,  1, -3,  1,  1,  3, -3, -1, -3, -1,  3, -3,  3, -1, -1, -1, -3},
                               { 1, -3, -1, -3,  3,  3, -1, -3,  1, -3, -3, -1, -3, -1,  1,  3,  3,  3},
                               {-3, -3,  1, -1, -1,  1,  1, -3, -1,  3,  3,  3,  3, -1,  3,  1,  3,  1},
                               { 3, -1, -3,  1, -3, -3, -3,  3,  3, -1,  1, -3, -1,  3,  1,  1,  3,  3},
                               { 3, -1, -1,  1, -3, -1, -3, -1, -3, -3, -1, -3,  1,  1,  1, -3, -3,  3},
                               {-3, -3,  1, -3,  3,  3,  3, -1,  3,  1,  1, -3, -3, -3,  3, -3, -1, -1},
                               {-3, -1, -1, -3,  1, -3,  3, -1, -1, -3,  3,  3, -3, -1,  3, -1, -1, -1},
                               {-3, -3,  3,  3, -3,  1,  3, -1, -3,  1, -1, -3,  3, -3, -1, -1, -1,  3},
                               {-1, -3,  1, -3, -3, -3,  1,  1,  3,  3, -3,  3,  3, -3, -1,  3, -3,  1},
                               {-3,  3,  1, -1, -1, -1, -1,  1, -1,  3,  3, -3, -1,  1,  3, -1,  3, -1}};
wnInt8 pi_fourPRBs[30][24] = {{-1, -3,  3, -1,  3,  1,  3, -1,  1, -3, -1, -3, -1,  1,  3, -3, -1, -3,  3,  3,  3, -3, -3, -3},
                              {-1, -3,  3,  1,  1, -3,  1, -3, -3,  1, -3, -1, -1,  3, -3,  3,  3,  3, -3,  1,  3,  3, -3, -3},
                              {-1, -3, -3,  1, -1, -1, -3,  1,  3, -1, -3, -1, -1, -3,  1,  1,  3,  1, -3, -1, -1,  3, -3, -3},
                              { 1, -3,  3, -1, -3, -1,  3,  3,  1, -1,  1,  1,  3, -3, -1, -3, -3, -3, -1,  3, -3, -1, -3, -3},
                              {-1,  3, -3, -3, -1,  3, -1, -1,  1,  3,  1,  3, -1, -1, -3,  1,  3,  1, -1, -3,  1, -1, -3, -3},
                              {-3, -1,  1, -3, -3,  1,  1, -3,  3, -1, -1, -3,  1,  3,  1, -1, -3, -1, -3,  1, -3, -3, -3, -3},
                              {-3,  3,  1,  3, -1,  1, -3,  1, -3,  1, -1, -3, -1, -3, -3, -3, -3, -1, -1, -1,  1,  1, -3, -3},
                              {-3,  1,  3, -1,  1, -1,  3, -3,  3, -1, -3, -1, -3,  3, -1, -1, -1, -3, -1, -1, -3,  3,  3, -3},
                              {-3,  1, -3,  3, -1, -1, -1, -3,  3,  1, -1, -3, -1,  1,  3, -1,  1, -1,  1, -3, -3, -3, -3, -3},
                              { 1,  1, -1, -3, -1,  1,  1, -3,  1, -1,  1, -3,  3, -3, -3,  3, -1, -3,  1,  3, -3,  1, -3, -3},
                              {-3, -3, -3, -1,  3, -3,  3,  1,  3,  1, -3, -1, -1, -3,  1,  1,  3,  1, -1, -3,  3,  1,  3, -3},
                              {-3,  3, -1,  3,  1, -1, -1, -1,  3,  3,  1,  1,  1,  3,  3,  1, -3, -3, -1,  1, -3,  1,  3, -3},
                              { 3, -3,  3, -1, -3,  1,  3,  1, -1, -1, -3, -1,  3, -3,  3, -1, -1,  3,  3, -3, -3,  3, -3, -3},
                              {-3,  3, -1,  3, -1,  3,  3,  1,  1, -3,  1,  3, -3,  3, -3, -3, -1,  1,  3, -3, -1, -1, -3, -3},
                              {-3,  1, -3, -1, -1,  3,  1,  3, -3,  1, -1,  3,  3, -1, -3,  3, -3, -1, -1, -3, -3, -3,  3, -3},
                              {-3, -1, -1, -3,  1, -3, -3, -1, -1,  3, -1,  1, -1,  3,  1, -3, -1,  3,  1,  1, -1, -1, -3, -3},
                              {-3, -3,  1, -1,  3,  3, -3, -1,  1, -1, -1,  1,  1, -1, -1,  3, -3,  1, -3,  1, -1, -1, -1, -3},
                              { 3, -1,  3, -1,  1, -3,  1,  1, -3, -3,  3, -3, -1, -1, -1, -1, -1, -3, -3, -1,  1,  1, -3, -3},
                              {-3,  1, -3,  1, -3, -3,  1, -3,  1, -3, -3, -3, -3, -3,  1, -3, -3,  1,  1, -3,  1,  1, -3, -3},
                              {-3, -3,  3,  3,  1, -1, -1, -1,  1, -3, -1,  1, -1,  3, -3, -1, -3, -1, -1,  1, -3,  3, -1, -3},
                              {-3, -3, -1, -1, -1, -3,  1, -1, -3, -1,  3, -3,  1, -3,  3, -3,  3,  3,  1, -1, -1,  1, -3, -3},
                              { 3, -1,  1, -1,  3, -3,  1,  1,  3, -1, -3,  3,  1, -3,  3, -1, -1, -1, -1,  1, -3, -3, -3, -3},
                              {-3,  1, -3,  3, -3,  1, -3,  3,  1, -1, -3, -1, -3, -3, -3, -3,  1,  3, -1,  1,  3,  3,  3, -3},
                              {-3, -1,  1, -3, -1, -1,  1,  1,  1,  3,  3, -1,  1, -1,  1, -1, -1, -3, -3, -3,  3,  1, -1, -3},
                              {-3,  3, -1, -3, -1, -1, -1,  3, -1, -1,  3, -3, -1,  3, -3,  3, -3, -1,  3,  1,  1, -1, -3, -3},
                              {-3,  1, -1, -3, -3, -1,  1, -3, -1, -3,  1,  1, -1,  1,  1,  3,  3,  3, -1,  1, -1,  1, -1, -3},
                              {-1,  3, -1, -1,  3,  3, -1, -1, -1,  3, -1, -3,  1,  3,  1,  1, -3, -3, -3, -1, -3, -1, -3, -3},
                              { 3, -3, -3, -1,  3,  3, -3, -1,  3,  1,  1,  1,  3, -1,  3, -3, -1,  3, -1,  3,  1, -1, -3, -3},
                              {-3,  1, -3,  1, -3,  1,  1,  3,  1, -3, -3, -1,  1,  3, -1, -3,  3,  1, -1, -3, -3, -3, -3, -3},
                              { 3, -3, -1,  1,  3, -1, -1, -3, -1,  3, -1, -3, -1, -3,  3, -1,  3,  1,  1, -3,  3, -3, -3, -3}};

cmplx nr_pusch_dmrs_data[MAX_PUSCH_INFO_PER_SYMBOL];

// Modulation
wnVoid wnNrPhyPuschModulation(cmplx *mod_out, wnUInt8 *mod_inp, wnUInt32 inp_len, wnUInt32 modulation_type)
{
    wnUInt32 i, j, l;

    if (modulation_type == 2) /* qpsk Modulation */
    {
        j = (inp_len+7)>>3;

        for (i = 0; i < j; i++)
        {
            memcpy(&mod_out[i<<2].real,&qpsk[mod_inp[i]][0],32);
        }
    }
    else if (modulation_type == 4) /* 16QAM Modulation */
    {
        j = (inp_len+7)>>3;

        for (i = 0; i < j; i++)
        {
            memcpy(&mod_out[i<<1].real,&q16[mod_inp[i]][0],16);
        }
    }
    else if (modulation_type == 6) /* 64QAM Modulation */
    {
        j = (inp_len+7)>>3;
        l = 0;

        for (i = 0; i < j; i = i + 3)
        {
            memcpy(&mod_out[l].real,&q64[mod_inp[i]&63][0],8);
            memcpy(&mod_out[l+1].real,&q64[((mod_inp[i]>>6)^(mod_inp[i+1]<<2))&63][0],8);
            memcpy(&mod_out[l+2].real,&q64[((mod_inp[i+1]>>4)^(mod_inp[i+2]<<4))&63][0],8);
            memcpy(&mod_out[l+3].real,&q64[(mod_inp[i+2]>>2)&63][0],8);
            
            l += 4;
        }
    }
    else if (modulation_type == 8) /* 256QAM Modulation */
    {
        j = (inp_len+7)>>3;

        for (i = 0; i < j; i++)
        {
            memcpy(&mod_out[i].real,&q256[mod_inp[i]][0],8);
        }
    }
    else
    {
        printf("ERROR :: Modulation not supported\n");
    }
}

/********************************* Start of DMRS Function ********************************/

wnVoid wnNrPhyPuschDmrs(puschConfigT* nrUlPhyPuschPDUs)
{
    // Variables used in the function
    wnUInt8 l_bar[4], l_bar_length;
    wnUInt32 iDx, iDx2, dmrs_sym_1prb_1sym, pnSeqOffset, dmrs_sequence_length, l_note, n_id_scid;//, C, dmrs_f, portNumber, portIndex,
    wnUInt64 seedVal;
    wnFlt* dmrsData;
    
    dmrsCdmg[0] = 0; dmrsCdmg[1] = 0; dmrsCdmg[2] = 0;

     if (nrUlPhyPuschPDUs->nSCID == 0)
        n_id_scid = nrUlPhyPuschPDUs->nNIDnSCID0;
     else
        n_id_scid = nrUlPhyPuschPDUs->nNIDnSCID1;
    
 
    // PDSCH DMRS signals transmission power 
    wnFlt puschDmrsPowerBoostFactor;

    if (nrUlPhyPuschPDUs->nDmrsWithoutData == 0)
        puschDmrsPowerBoostFactor = 1;
    else if (nrUlPhyPuschPDUs->nDmrsWithoutData == 1)
        puschDmrsPowerBoostFactor = pow(10,0.15);
    else if (nrUlPhyPuschPDUs->nDmrsWithoutData == 2)
        puschDmrsPowerBoostFactor = pow(10,0.2385);

    //puschDmrsPowerBoostFactor = 1;
    // Number of DMRS symbols in 1 OFDM symbol in 1 PRB (1 PRB in frequency)
    if (nrUlPhyPuschPDUs->nDMRSConfigType == 0)
    {
        dmrs_sym_1prb_1sym = 3 * nrUlPhyPuschPDUs->nBWPSize;
        pnSeqOffset = 6 * (nrUlPhyPuschPDUs->nRBStart + nrUlPhyPuschPDUs->nBWPStart);
        // C = 2;
    }    
    else if (nrUlPhyPuschPDUs->nDMRSConfigType == 1)
    {
        dmrs_sym_1prb_1sym = nrUlPhyPuschPDUs->nBWPSize<<1;
        pnSeqOffset = (nrUlPhyPuschPDUs->nRBStart + nrUlPhyPuschPDUs->nBWPStart)<<2;
        // C = 3;
    }


    // Total number of DMRS symbols produced
    dmrs_sequence_length = (dmrs_sym_1prb_1sym<<1) + pnSeqOffset;


    // Reference point for l and the position 0 l of the first DM-RS symbol  
    // If NR PDSCH Mapping type A
    if (nrUlPhyPuschPDUs->nMappingType == 0)
    {
        if (nrUlPhyPuschPDUs->ulDmrsTypeAPos == 3)
            l_note = 3;
        else 
            l_note = 2;
    }
    // If NR PDSCH Mapping type B
    else if (nrUlPhyPuschPDUs->nMappingType == 1)
        l_note = 0;


    // Calculating the position of DMRS (time scale) (look-up table for l_bar) 
    // If it is a single symbol DMRS
    if (nrUlPhyPuschPDUs->nNrOfDMRSSymbols == 1)
    {
        // Temporary variable to move in table from left to right 
        wnUInt32 temp_var = (nrUlPhyPuschPDUs->nMappingType<<2) + nrUlPhyPuschPDUs->nDMRSAddPos;                                                     

        l_bar_length = puschSingleSymTbl[nrUlPhyPuschPDUs->nNrOfSymbols-1][temp_var][0];
        memcpy(&l_bar[0],&puschSingleSymTbl[nrUlPhyPuschPDUs->nNrOfSymbols-1][temp_var][1],l_bar_length);
        
        if (l_bar[0] == x)
            l_bar[0] = l_note;
   
    }
    // If it is a double symbol DMRS    
    else 
    {
        // Temporary variable to move in table from left to right 
        wnUInt32 temp_var = nrUlPhyPuschPDUs->nMappingType*4 + nrUlPhyPuschPDUs->nDMRSAddPos;                                                     
  
        l_bar_length = puschDoubleSymTbl[nrUlPhyPuschPDUs->nNrOfSymbols-1][temp_var][0];
        memcpy(&l_bar[0],&puschDoubleSymTbl[nrUlPhyPuschPDUs->nNrOfSymbols-1][temp_var][1],l_bar_length);
       if (l_bar[0] == x)
            l_bar[0] = l_note;
    }


    // Following loop calculate all the OFDM symbol positions where DMRS will be transmitted
    // Array 'dmrs_s' contains all the positions
    iDx2 = 0;
    for (iDx = 0; iDx < l_bar_length; iDx++)
    {        
        dmrs_s[iDx2] = l_bar[iDx] + nrUlPhyPuschPDUs->nMappingType*nrUlPhyPuschPDUs->nStartSymbolIndex;
        iDx2++;
        // IF double symbol DMRS
        if (nrUlPhyPuschPDUs->nNrOfDMRSSymbols == 2)
        {        
            dmrs_s[iDx2] = l_bar[iDx] + nrUlPhyPuschPDUs->nMappingType*nrUlPhyPuschPDUs->nStartSymbolIndex + 1;
            iDx2++;
        }
    }


    // wnUInt32 imagOffset = TOTAL_RES_IN_SYMBOL>>1;
    wnUInt8 DoubleDMRSIdx = 0;
    if (nrUlPhyPuschPDUs->nDMRSConfigType == 0)
    {
        //If transform precoding is disabled.
        if(nrUlPhyPuschPDUs->transformPrecode == 0)
        {
            for (iDx = 0; iDx < l_bar_length; iDx++)
            {
                for(DoubleDMRSIdx = 0;DoubleDMRSIdx<nrUlPhyPuschPDUs->nNrOfDMRSSymbols;DoubleDMRSIdx++)
                {
                    // DMRS sequence generation
                    // DMRS sequence is being generated for each symbols separately 
                    seedVal = (((14*nrUlPhyPuschPDUs->nSlotNumber + dmrs_s[iDx*nrUlPhyPuschPDUs->nNrOfDMRSSymbols+DoubleDMRSIdx] + 1)*((n_id_scid<<1)+1)<<17) + (n_id_scid<<1) + nrUlPhyPuschPDUs->nSCID)&0x7FFFFFFF;

                    // PN sequence generation function
                    wnNrPhyPuschPnSequence(seedVal, dmrs_sequence_length<<1, pnSeqOut);

                    // DMRS for PUSCH
                    wnNrPhyPuschModulation(nr_pusch_dmrs_data,pnSeqOut,dmrs_sequence_length<<1,2);
            
                    if (iDx == 0)
                        dmrsData = (wnFlt*)&nrPuschOutParamsSeg1.genDMRS[DoubleDMRSIdx][0];
                    else if(iDx == 1)
                        dmrsData = (wnFlt*)&nrPuschOutParamsSeg2.genDMRS[DoubleDMRSIdx][0];
                    else if (iDx == 2)
                        dmrsData = (wnFlt*)&nrPuschOutParamsSeg3.genDMRS[DoubleDMRSIdx][0];
                    else
                        dmrsData = (wnFlt*)&nrPuschOutParamsSeg4.genDMRS[DoubleDMRSIdx][0];
                        
                    // This is for DMRS symbols in each OFDM symbol.
                    for (iDx2 = 0; iDx2 < (dmrs_sym_1prb_1sym<<1); iDx2 += 2)
                    {  
                        dmrsData[pnSeqOffset+iDx2] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2].real;
                        dmrsData[HALF_MAX_NUM_SC+pnSeqOffset+iDx2] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2].imag;

                        dmrsData[pnSeqOffset+iDx2+1] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2 + 1].real;
                        dmrsData[HALF_MAX_NUM_SC+pnSeqOffset+iDx2+1] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2 + 1].imag;
                    }
                }
            }
        }
        else
        {
            if(nrUlPhyPuschPDUs->Rel15 == 1)
            {
                wnUInt8 gh, sh, f_gh, v, u;
                wnUInt16 N_zc, q, M_zc;
                wnUInt32 seqLen, tempMem, temp, itr;
                wnFlt q_bar;
                cmplx xq[1637];//maximum possible value of N_zc = 1637.
                wnFlt oneBy31 = 1/(wnFlt)31;

                gh = nrUlPhyPuschPDUs->ghFlag;//base group hopping enable(1)/disable(0).
                sh = nrUlPhyPuschPDUs->shFlag;//base sequence hopping enable(1)/disable(0).

                dmrs_sequence_length = nrUlPhyPuschPDUs->nRBSize*6;

                M_zc = dmrs_sequence_length;

                if(gh == 0 && sh == 0)
                {
                    f_gh = 0;
                    v = 0;

                    u = (f_gh+nrUlPhyPuschPDUs->NIdRS)%30;
                    // printf("Fix DMRS Generation for DFT-S-OFDM\n");
                    
                    // Generate Zad-off chu sequence.    
                    if(nrUlPhyPuschPDUs->nRBSize>=6)
                    {
                        N_zc = maxPrimesOfMzc[M_zc/6-1];
                        
                        q_bar = N_zc*(u+1)/(wnFlt)31;
                        q = (wnUInt16)(q_bar+0.5)+v*(1-2*( ((wnUInt16)(2*q_bar))&1 ));
                        // FILE *fpcos = fopen("testCos&sin/cos.txt", "w");
                        // FILE *fpsin = fopen("testCos&sin/sin.txt", "w");

                        for(itr = 0;itr<N_zc;itr++)
                        {
                            xq[itr].real = cos((double)-1*PI*q*itr*(itr+1)/N_zc);
                            xq[itr].imag = sin((double)-1*PI*q*itr*(itr+1)/N_zc);
                            
                            // xq[itr].real = wnCosRad((double)-1*PI*q*itr*(itr+1)/N_zc);
                            // xq[itr].imag = wnSinRad((double)-1*PI*q*itr*(itr+1)/N_zc);

                            // printf("xq[itr].imag: %f\n", xq[itr].imag );
                            // fprintf(fpcos, "%f\n", xq[itr].real );
                            // fprintf(fpsin, "%f\n", xq[itr].imag );
                        }
                        // fclose(fpcos);
                        // fclose(fpsin);
                        
                        for(itr = 0;itr<M_zc;itr++)
                        {
                            nr_pusch_dmrs_data[itr].real = xq[itr%N_zc].real;
                            nr_pusch_dmrs_data[itr].imag = xq[itr%N_zc].imag;
                        }
                    }
                    else if(nrUlPhyPuschPDUs->nRBSize==1)
                    {
                        //Computer Generated Sequence.
                        for(itr = 0; itr<6;itr++)
                        {
                            nr_pusch_dmrs_data[itr].real = wnCosRad(pi_onePRB[u][itr] * PI*0.25);
                            nr_pusch_dmrs_data[itr].imag = wnSinRad(pi_onePRB[u][itr] * PI*0.25);
                        }
                    }
                    else if(nrUlPhyPuschPDUs->nRBSize==2)
                    {
                        //Computer Generated Sequence.
                        for(itr = 0; itr<12;itr++)
                        {
                            nr_pusch_dmrs_data[itr].real = wnCosRad(pi_twoPRBs[u][itr] * PI*0.25);
                            nr_pusch_dmrs_data[itr].imag = wnSinRad(pi_twoPRBs[u][itr] * PI*0.25);
                        }
                    }
                    else if(nrUlPhyPuschPDUs->nRBSize==3)
                    {
                        //Computer Generated Sequence.
                        for(itr = 0; itr<18;itr++)
                        {
                            nr_pusch_dmrs_data[itr].real = wnCosRad(pi_threePRBs[u][itr] * PI*0.25);
                            nr_pusch_dmrs_data[itr].imag = wnSinRad(pi_threePRBs[u][itr] * PI*0.25);
                        }
                    }
                    else if(nrUlPhyPuschPDUs->nRBSize==4)
                    {
                        //Computer Generated Sequence.
                        for(itr = 0; itr<24;itr++)
                        {
                            nr_pusch_dmrs_data[itr].real = wnCosRad(pi_fourPRBs[u][itr] * PI*0.25);
                            nr_pusch_dmrs_data[itr].imag = wnSinRad(pi_fourPRBs[u][itr] * PI*0.25);
                        }
                    }
                    else //if(nrUlPhyPuschPDUs->nRBSize==5)
                    {
                        for(itr = 0; itr<30;itr++)
                        {
                            nr_pusch_dmrs_data[itr].real = wnCosRad(-1*(u+1)*(itr+1)*(itr+2)* PI*oneBy31);
                            nr_pusch_dmrs_data[itr].imag = wnSinRad(-1*(u+1)*(itr+1)*(itr+2)* PI*oneBy31);
                        }
                    }

                    //Base Sequence will be same across symbols and across slots i.e., base sequence is 
                    //independent of DMRS symbol index and slot number within a frame in this case.
                    for (iDx = 0; iDx < l_bar_length; iDx++)
                    {
                        for(DoubleDMRSIdx = 0;DoubleDMRSIdx<nrUlPhyPuschPDUs->nNrOfDMRSSymbols;DoubleDMRSIdx++)
                        {
                            if (iDx == 0)
                                dmrsData = (wnFlt*)&nrPuschOutParamsSeg1.genDMRS[DoubleDMRSIdx][0];
                            else if(iDx == 1)
                                dmrsData = (wnFlt*)&nrPuschOutParamsSeg2.genDMRS[DoubleDMRSIdx][0];
                            else if (iDx == 2)
                                dmrsData = (wnFlt*)&nrPuschOutParamsSeg3.genDMRS[DoubleDMRSIdx][0];
                            else
                                dmrsData = (wnFlt*)&nrPuschOutParamsSeg4.genDMRS[DoubleDMRSIdx][0];
                                
                            // This is for DMRS symbols in each OFDM symbol.
                            for (iDx2 = 0; iDx2 < (dmrs_sym_1prb_1sym<<1); iDx2 += 2)
                            {  
                                dmrsData[pnSeqOffset+iDx2] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[iDx2].real;
                                dmrsData[HALF_MAX_NUM_SC+pnSeqOffset+iDx2] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[iDx2].imag;

                                dmrsData[pnSeqOffset+iDx2+1] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[iDx2 + 1].real;
                                dmrsData[HALF_MAX_NUM_SC+pnSeqOffset+iDx2+1] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[iDx2 + 1].imag;
                            }
                        }
                    }
                }
                else if(gh == 1 && sh == 0)
                {
                    seedVal = (wnUInt64)nrUlPhyPuschPDUs->NIdRS/30;
                    for (iDx = 0; iDx < l_bar_length; iDx++)
                    {
                        for(DoubleDMRSIdx = 0;DoubleDMRSIdx<nrUlPhyPuschPDUs->nNrOfDMRSSymbols;DoubleDMRSIdx++)
                        {
                            tempMem = 8*(14*nrUlPhyPuschPDUs->nSlotNumber+ dmrs_s[iDx*nrUlPhyPuschPDUs->nNrOfDMRSSymbols+DoubleDMRSIdx] );
                            seqLen = tempMem+8;
                            // PN sequence generation function
                            wnNrPhyPuschPnSequence(seedVal, seqLen, pnSeqOut);
                            temp = 0;
                            for(itr = 0 ;itr<8;itr++)
                            {
                                temp = temp+(pnSeqOut[tempMem+itr]*(1<<itr));
                            }

                            f_gh = temp%30;
                            v = 0;

                            u = (f_gh+nrUlPhyPuschPDUs->NIdRS)%30;
                            // Generate Zadd-off chu sequence.    
                            if(nrUlPhyPuschPDUs->nRBSize>=6)
                            {
                                N_zc = maxPrimesOfMzc[M_zc/6];
                                q_bar = N_zc*(u+1)/(wnFlt)31;
                                q = (wnUInt16)(q_bar+0.5)+v*(1-2*( ((wnUInt16)(2*q_bar))&1 ));
                                
                                for(itr = 0;itr<N_zc;itr++)
                                {
                                    xq[itr].real = wnCosRad(-1*PI*q*itr*(itr+1)/N_zc);
                                    xq[itr].imag = wnSinRad(-1*PI*q*itr*(itr+1)/N_zc);
                                }
                                
                                for(itr = 0;itr<M_zc;itr++)
                                {
                                    nr_pusch_dmrs_data[itr].real = xq[itr%N_zc].real;
                                    nr_pusch_dmrs_data[itr].imag = xq[itr%N_zc].imag;
                                }
                            }
                            //PNsequence offset????.
                            else if(nrUlPhyPuschPDUs->nRBSize==1)
                            {
                                //Computer Generated Sequence.
                                for(itr = 0; itr<6;itr++)
                                {
                                    nr_pusch_dmrs_data[itr].real = wnCosRad(pi_onePRB[u][itr] * PI*0.25);
                                    nr_pusch_dmrs_data[itr].imag = wnSinRad(pi_onePRB[u][itr] * PI*0.25);
                                }
                            }
                            else if(nrUlPhyPuschPDUs->nRBSize==2)
                            {
                                //Computer Generated Sequence.
                                for(itr = 0; itr<12;itr++)
                                {
                                    nr_pusch_dmrs_data[itr].real = wnCosRad(pi_twoPRBs[u][itr] * PI*0.25);
                                    nr_pusch_dmrs_data[itr].imag = wnSinRad(pi_twoPRBs[u][itr] * PI*0.25);
                                }
                            }
                            else if(nrUlPhyPuschPDUs->nRBSize==3)
                            {
                                //Computer Generated Sequence.
                                for(itr = 0; itr<18;itr++)
                                {
                                    nr_pusch_dmrs_data[itr].real = wnCosRad(pi_threePRBs[u][itr] * PI*0.25);
                                    nr_pusch_dmrs_data[itr].imag = wnSinRad(pi_threePRBs[u][itr] * PI*0.25);
                                }
                            }
                            else if(nrUlPhyPuschPDUs->nRBSize==4)
                            {
                                //Computer Generated Sequence.
                                for(itr = 0; itr<24;itr++)
                                {
                                    nr_pusch_dmrs_data[itr].real = wnCosRad(pi_fourPRBs[u][itr] * PI*0.25);
                                    nr_pusch_dmrs_data[itr].imag = wnSinRad(pi_fourPRBs[u][itr] * PI*0.25);
                                }
                            }
                            else //if(nrUlPhyPuschPDUs->nRBSize==5)
                            {
                                for(itr = 0; itr<30;itr++)
                                {
                                    nr_pusch_dmrs_data[itr].real = wnCosRad(-1*(u+1)*(itr+1)*(itr+2)* PI*oneBy31);
                                    nr_pusch_dmrs_data[itr].imag = wnSinRad(-1*(u+1)*(itr+1)*(itr+2)* PI*oneBy31);
                                }
                            }
                    
                            if (iDx == 0)
                                dmrsData = (wnFlt*)&nrPuschOutParamsSeg1.genDMRS[DoubleDMRSIdx][0];
                            else if(iDx == 1)
                                dmrsData = (wnFlt*)&nrPuschOutParamsSeg2.genDMRS[DoubleDMRSIdx][0];
                            else if (iDx == 2)
                                dmrsData = (wnFlt*)&nrPuschOutParamsSeg3.genDMRS[DoubleDMRSIdx][0];
                            else
                                dmrsData = (wnFlt*)&nrPuschOutParamsSeg4.genDMRS[DoubleDMRSIdx][0];
                                
                            // This is for DMRS symbols in each OFDM symbol
                            for (iDx2 = 0; iDx2 < (dmrs_sym_1prb_1sym<<1); iDx2 += 2)
                            {  
                                dmrsData[pnSeqOffset+iDx2] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2].real;
                                dmrsData[HALF_MAX_NUM_SC+pnSeqOffset+iDx2] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2].imag;

                                dmrsData[pnSeqOffset+iDx2+1] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2 + 1].real;
                                dmrsData[HALF_MAX_NUM_SC+pnSeqOffset+iDx2+1] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2 + 1].imag;
                            }
                        
                        }
                    }
                }
                else if(gh == 0 && sh == 1)
                {
                    f_gh = 0;
                    v = 0;
                    seedVal = (wnUInt64)nrUlPhyPuschPDUs->NIdRS;
                    
                    if(nrUlPhyPuschPDUs->nRBSize>=12)
                    {
                        for (iDx = 0; iDx < l_bar_length; iDx++)
                        {
                            for(DoubleDMRSIdx = 0;DoubleDMRSIdx<nrUlPhyPuschPDUs->nNrOfDMRSSymbols;DoubleDMRSIdx++)
                            {
                                seqLen = 14*nrUlPhyPuschPDUs->nSlotNumber + dmrs_s[iDx*nrUlPhyPuschPDUs->nNrOfDMRSSymbols+DoubleDMRSIdx];
                                // PN sequence generation function
                                wnNrPhyPuschPnSequence(seedVal, seqLen, pnSeqOut);
                                v = pnSeqOut[14*nrUlPhyPuschPDUs->nSlotNumber + dmrs_s[iDx*nrUlPhyPuschPDUs->nNrOfDMRSSymbols+DoubleDMRSIdx]];
                            
                                u = (f_gh+nrUlPhyPuschPDUs->NIdRS)%30;
                                // Generate Zadd-off chu sequence.    
                                if(nrUlPhyPuschPDUs->nRBSize>=6)
                                {
                                    N_zc = maxPrimesOfMzc[M_zc/6];
                                    q_bar = N_zc*(u+1)/(wnFlt)31;
                                    q = (wnUInt16)(q_bar+0.5)+v*(1-2*( ((wnUInt16)(2*q_bar))&1 ));
                                    
                                    for(itr = 0;itr<N_zc;itr++)
                                    {
                                        xq[itr].real = wnCosRad(-1*PI*q*itr*(itr+1)/N_zc);
                                        xq[itr].imag = wnSinRad(-1*PI*q*itr*(itr+1)/N_zc);
                                    }
                                    
                                    for(itr = 0;itr<M_zc;itr++)
                                    {
                                        nr_pusch_dmrs_data[itr].real = xq[itr%N_zc].real;
                                        nr_pusch_dmrs_data[itr].imag = xq[itr%N_zc].imag;
                                    }
                                }
                                //PNsequence offset????.
                                //Computer Generated Sequence.
                                else if(nrUlPhyPuschPDUs->nRBSize==1)
                                {
                                    for(itr = 0; itr<6;itr++)
                                    {
                                        nr_pusch_dmrs_data[itr].real = wnCosRad(pi_onePRB[u][itr] * PI*0.25);
                                        nr_pusch_dmrs_data[itr].imag = wnSinRad(pi_onePRB[u][itr] * PI*0.25);
                                    }
                                }
                                //Computer Generated Sequence.
                                else if(nrUlPhyPuschPDUs->nRBSize==2)
                                {
                                    for(itr = 0; itr<12;itr++)
                                    {
                                        nr_pusch_dmrs_data[itr].real = wnCosRad(pi_twoPRBs[u][itr] * PI*0.25);
                                        nr_pusch_dmrs_data[itr].imag = wnSinRad(pi_twoPRBs[u][itr] * PI*0.25);
                                    }
                                }
                                //Computer Generated Sequence.
                                else if(nrUlPhyPuschPDUs->nRBSize==3)
                                {
                                    for(itr = 0; itr<18;itr++)
                                    {
                                        nr_pusch_dmrs_data[itr].real = wnCosRad(pi_threePRBs[u][itr] * PI*0.25);
                                        nr_pusch_dmrs_data[itr].imag = wnSinRad(pi_threePRBs[u][itr] * PI*0.25);
                                    }
                                }
                                //Computer Generated Sequence.
                                else if(nrUlPhyPuschPDUs->nRBSize==4)
                                {
                                    for(itr = 0; itr<24;itr++)
                                    {
                                        nr_pusch_dmrs_data[itr].real = wnCosRad(pi_fourPRBs[u][itr] * PI*0.25);
                                        nr_pusch_dmrs_data[itr].imag = wnSinRad(pi_fourPRBs[u][itr] * PI*0.25);
                                    }
                                }
                                else //if(nrUlPhyPuschPDUs->nRBSize==5)
                                {
                                    for(itr = 0; itr<30;itr++)
                                    {
                                        nr_pusch_dmrs_data[itr].real = wnCosRad(-1*(u+1)*(itr+1)*(itr+2)* PI*oneBy31);
                                        nr_pusch_dmrs_data[itr].imag = wnSinRad(-1*(u+1)*(itr+1)*(itr+2)* PI*oneBy31);
                                    }
                                }
                        
                                if (iDx == 0)
                                    dmrsData = (wnFlt*)&nrPuschOutParamsSeg1.genDMRS[DoubleDMRSIdx][0];
                                else if(iDx == 1)
                                    dmrsData = (wnFlt*)&nrPuschOutParamsSeg2.genDMRS[DoubleDMRSIdx][0];
                                else if (iDx == 2)
                                    dmrsData = (wnFlt*)&nrPuschOutParamsSeg3.genDMRS[DoubleDMRSIdx][0];
                                else
                                    dmrsData = (wnFlt*)&nrPuschOutParamsSeg4.genDMRS[DoubleDMRSIdx][0];
                                    
                                // This is for DMRS symbols in each OFDM symbol.
                                for (iDx2 = 0; iDx2 < (dmrs_sym_1prb_1sym<<1); iDx2 += 2)
                                {  
                                    dmrsData[pnSeqOffset+iDx2] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2].real;
                                    dmrsData[HALF_MAX_NUM_SC+pnSeqOffset+iDx2] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2].imag;

                                    dmrsData[pnSeqOffset+iDx2+1] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2 + 1].real;
                                    dmrsData[HALF_MAX_NUM_SC+pnSeqOffset+iDx2+1] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2 + 1].imag;
                                }
                
                            }
                        }
                    }
                    else
                    {
                        v = 0;

                        u = (f_gh+nrUlPhyPuschPDUs->NIdRS)%30;
                        // Generate Zadd-off chu sequence.    
                        if(nrUlPhyPuschPDUs->nRBSize>=6)
                        {
                            N_zc = maxPrimesOfMzc[M_zc/6];
                            q_bar = N_zc*(u+1)/(wnFlt)31;
                            q = (wnUInt16)(q_bar+0.5)+v*(1-2*( ((wnUInt16)(2*q_bar))&1 ));
                            
                            for(itr = 0;itr<N_zc;itr++)
                            {
                                xq[itr].real = wnCosRad(-1*PI*q*itr*(itr+1)/N_zc);
                                xq[itr].imag = wnSinRad(-1*PI*q*itr*(itr+1)/N_zc);
                            }
                            
                            for(itr = 0;itr<M_zc;itr++)
                            {
                                nr_pusch_dmrs_data[itr].real = xq[itr%N_zc].real;
                                nr_pusch_dmrs_data[itr].imag = xq[itr%N_zc].imag;
                            }
                        }
                        //PNsequence offset????.
                        //Computer Generated Sequence.
                        else if(nrUlPhyPuschPDUs->nRBSize==1)
                        {
                            for(itr = 0; itr<6;itr++)
                            {
                                nr_pusch_dmrs_data[itr].real = wnCosRad(pi_onePRB[u][itr] * PI*0.25);
                                nr_pusch_dmrs_data[itr].imag = wnSinRad(pi_onePRB[u][itr] * PI*0.25);
                            }
                        }
                        //Computer Generated Sequence.
                        else if(nrUlPhyPuschPDUs->nRBSize==2)
                        {
                            for(itr = 0; itr<12;itr++)
                            {
                                nr_pusch_dmrs_data[itr].real = wnCosRad(pi_twoPRBs[u][itr] * PI*0.25);
                                nr_pusch_dmrs_data[itr].imag = wnSinRad(pi_twoPRBs[u][itr] * PI*0.25);
                            }
                        }
                        //Computer Generated Sequence.
                        else if(nrUlPhyPuschPDUs->nRBSize==3)
                        {
                            for(itr = 0; itr<18;itr++)
                            {
                                nr_pusch_dmrs_data[itr].real = wnCosRad(pi_threePRBs[u][itr] * PI*0.25);
                                nr_pusch_dmrs_data[itr].imag = wnSinRad(pi_threePRBs[u][itr] * PI*0.25);
                            }
                        }
                        //Computer Generated Sequence.
                        else if(nrUlPhyPuschPDUs->nRBSize==4)
                        {
                            for(itr = 0; itr<24;itr++)
                            {
                                nr_pusch_dmrs_data[itr].real = wnCosRad(pi_fourPRBs[u][itr] * PI*0.25);
                                nr_pusch_dmrs_data[itr].imag = wnSinRad(pi_fourPRBs[u][itr] * PI*0.25);
                            }
                        }
                        else //if(nrUlPhyPuschPDUs->nRBSize==5)
                        {
                            for(itr = 0; itr<30;itr++)
                            {
                                nr_pusch_dmrs_data[itr].real = wnCosRad(-1*(u+1)*(itr+1)*(itr+2)* PI*oneBy31);
                                nr_pusch_dmrs_data[itr].imag = wnSinRad(-1*(u+1)*(itr+1)*(itr+2)* PI*oneBy31);
                            }
                        }

                        for (iDx = 0; iDx < l_bar_length; iDx++)
                        {
                            for(DoubleDMRSIdx = 0;DoubleDMRSIdx<nrUlPhyPuschPDUs->nNrOfDMRSSymbols;DoubleDMRSIdx++)
                            {
                                                      
                                if (iDx == 0)
                                    dmrsData = (wnFlt*)&nrPuschOutParamsSeg1.genDMRS[DoubleDMRSIdx][0];
                                else if(iDx == 1)
                                    dmrsData = (wnFlt*)&nrPuschOutParamsSeg2.genDMRS[DoubleDMRSIdx][0];
                                else if (iDx == 2)
                                    dmrsData = (wnFlt*)&nrPuschOutParamsSeg3.genDMRS[DoubleDMRSIdx][0];
                                else
                                    dmrsData = (wnFlt*)&nrPuschOutParamsSeg4.genDMRS[DoubleDMRSIdx][0];
                                    
                                // This is for DMRS symbols in each OFDM symbol.
                                for (iDx2 = 0; iDx2 < (dmrs_sym_1prb_1sym<<1); iDx2 += 2)
                                {  
                                    dmrsData[pnSeqOffset+iDx2] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2].real;
                                    dmrsData[HALF_MAX_NUM_SC+pnSeqOffset+iDx2] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2].imag;

                                    dmrsData[pnSeqOffset+iDx2+1] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2 + 1].real;
                                    dmrsData[HALF_MAX_NUM_SC+pnSeqOffset+iDx2+1] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2 + 1].imag;
                                }
                            }
                        }
                    }                          
                }
            }
        }
    }
    else
    {
        for (iDx = 0; iDx < l_bar_length*nrUlPhyPuschPDUs->nNrOfDMRSSymbols; iDx++)
        {
            for(DoubleDMRSIdx = 0;DoubleDMRSIdx<nrUlPhyPuschPDUs->nNrOfDMRSSymbols;DoubleDMRSIdx++)
            {
                // DMRS sequence generation
                // DMRS sequence is being generated for each symbols separately      
                seedVal = (((14*nrUlPhyPuschPDUs->nSlotNumber + dmrs_s[iDx*nrUlPhyPuschPDUs->nNrOfDMRSSymbols+DoubleDMRSIdx] + 1)*((n_id_scid<<1)+1)<<17) + (n_id_scid<<1) + nrUlPhyPuschPDUs->nSCID)&0x7FFFFFFF;

                // PN sequence generation function
                wnNrPhyPuschPnSequence(seedVal, dmrs_sequence_length<<1, pnSeqOut);

                // DMRS for PDSCH
                wnNrPhyPuschModulation(nr_pusch_dmrs_data,pnSeqOut,dmrs_sequence_length<<1,2);                               
               
                if (iDx == 0)
                    dmrsData = (wnFlt*)&nrPuschOutParamsSeg1.genDMRS[DoubleDMRSIdx][0];
                else if(iDx == 1)
                    dmrsData = (wnFlt*)&nrPuschOutParamsSeg2.genDMRS[DoubleDMRSIdx][0];
                else if(iDx == 2)
                    dmrsData = (wnFlt*)&nrPuschOutParamsSeg3.genDMRS[DoubleDMRSIdx][0];
                else if(iDx == 3)
                    dmrsData = (wnFlt*)&nrPuschOutParamsSeg4.genDMRS[DoubleDMRSIdx][0];

                
                // This is for DMRS symbols in each OFDM symbol
                for (iDx2 = 0; iDx2 < (dmrs_sym_1prb_1sym<<1); iDx2 += 2)
                { 
                    dmrsData[pnSeqOffset+iDx2] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2].real;
                    dmrsData[HALF_MAX_NUM_SC+pnSeqOffset+iDx2] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2].imag;

                    dmrsData[pnSeqOffset+iDx2+1] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2 + 1].real;
                    dmrsData[HALF_MAX_NUM_SC+pnSeqOffset+iDx2+1] = puschDmrsPowerBoostFactor*nr_pusch_dmrs_data[pnSeqOffset + iDx2 + 1].imag;
                }
            }
        }
    }



    if (nrUlPhyPuschPDUs->nDmrsWithoutData == 0)
    {
        dmrsCdmGWithoutData[0] = 1;dmrsCdmGWithoutData[1] = 0;dmrsCdmGWithoutData[2] = 0;
        if (nrUlPhyPuschPDUs->nDMRSConfigType == 0)
            dmrsCdmGWithData = 1;
        else
            dmrsCdmGWithData = 2;
    }
    else if (nrUlPhyPuschPDUs->nDmrsWithoutData == 1)
    {
        dmrsCdmGWithoutData[0] = 1;dmrsCdmGWithoutData[1] = 1;dmrsCdmGWithoutData[2] = 0;
        if (nrUlPhyPuschPDUs->nDMRSConfigType == 0)
            dmrsCdmGWithData = 0;
        else
            dmrsCdmGWithData = 1;
    }
    else if (nrUlPhyPuschPDUs->nDmrsWithoutData == 2)
    {
        dmrsCdmGWithoutData[0] = 1;dmrsCdmGWithoutData[1] = 1;dmrsCdmGWithoutData[2] = 1;
        if (nrUlPhyPuschPDUs->nDMRSConfigType == 0)
            dmrsCdmGWithData = 0;
        else
            dmrsCdmGWithData = 0;
    }
    
    // printf("l_bar_length: %u\n", l_bar_length);
    dmrs_pr_prb = (dmrs_sym_1prb_1sym<<1)*l_bar_length*nrUlPhyPuschPDUs->nNrOfDMRSSymbols/ nrUlPhyPuschPDUs->nBWPSize;
    dmrs_symbols = l_bar_length*nrUlPhyPuschPDUs->nNrOfDMRSSymbols;

    if(nrUlPhyPuschPDUs->nNrOfDMRSSymbols == 1)
    {
        for(iDx = 0;iDx<l_bar_length;iDx++)
        {
            if(iDx == 0)
                nrUlPhyPuschPDUs->dmrs1_sym = dmrs_s[0];
            else if(iDx == 1)
                nrUlPhyPuschPDUs->dmrs2_sym = dmrs_s[1];
            else if (iDx == 2)
                nrUlPhyPuschPDUs->dmrs3_sym = dmrs_s[2];
            else if (iDx == 3)
                nrUlPhyPuschPDUs->dmrs4_sym = dmrs_s[3];
        }
    }
    else
    {
        for(iDx = 0;iDx<l_bar_length;iDx++)
        {
            if(iDx == 0)
                nrUlPhyPuschPDUs->dmrs1_sym = dmrs_s[0];
            else if(iDx == 1)
                nrUlPhyPuschPDUs->dmrs2_sym = dmrs_s[2];
        }
    }
    
    
    
    //printf ("%d,%d,%d,%d\n",(dmrs_sym_1prb_1sym<<1),dmrs_pr_prb,nrUlPhyPuschPDUs->nNrOfDMRSSymbols,nrUlPhyPuschPDUs->nBWPSize);
return;
}



wnVoid wnSegBoundaries(wnUInt8 *nSeg, wnUInt8* startSym, wnUInt8* endSym, puschConfigT* pNrPuschInParams)
{
    wnUInt8 numOfSymbs = pNrPuschInParams->nNrOfSymbols;
    // wnUInt8 symbStartIdx = pNrPuschInParams->nStartSymbolIndex;
    wnUInt8 MappingType = pNrPuschInParams->nMappingType;
    wnUInt8 NrOfDMRSSymbols = pNrPuschInParams->nNrOfDMRSSymbols;
    wnUInt8 addPos = pNrPuschInParams->nDMRSAddPos;

    if(NrOfDMRSSymbols == 1)// Single Symbol DM-RS.
    {
        if(MappingType == 0)//Mapping Type A.
        {
            if(numOfSymbs<=7 || ((numOfSymbs<=14)&& addPos == 0) )
            {
                *nSeg = 1;
            }
            else if(numOfSymbs<=9 || ((numOfSymbs<=14)&& addPos == 1))
            {
                *nSeg = 2;
            }
            else if(numOfSymbs<=11 || ((numOfSymbs<=14)&& addPos == 2))
            {
                *nSeg = 3;
            }
            else
            {
                *nSeg = 4;
            }
        }
        else//Mapping Type B.
        {
            if(numOfSymbs<=4 || ((numOfSymbs<=14)&& addPos == 0) )
            {
                *nSeg = 1;
            }
            else if(numOfSymbs<=7 || ((numOfSymbs<=14)&& addPos == 1))
            {
                *nSeg = 2;
            }
            else if(numOfSymbs<=9 || ((numOfSymbs<=14)&& addPos == 2))
            {
                *nSeg = 3;
            }
            else
            {
                *nSeg = 4;
            }
        }
    }
    else//Double Symbol DM-RS.
    {
        if(MappingType == 0)//Mapping Type A.
        {
            if(numOfSymbs<= 9 || ((numOfSymbs<=14)&& addPos == 0) )
            {
                *nSeg = 1;
            }
            else //if( ((numOfSymbs<=14) && addPos == 1))
            {
                *nSeg = 2;
            }
        }
        else//Mapping Type B.
        {
            if(numOfSymbs<=7 || ((numOfSymbs<=14)&& addPos == 0) )
            {
                *nSeg = 1;
            }
            else //if( ((numOfSymbs<=14) && addPos == 1) )
            {
                *nSeg = 2;
            }
        }
    }

    
    if(*nSeg == 1)
    {
        startSym[0] = 0;
        endSym[0] = numOfSymbs-1;
        /*if(NrOfDMRSSymbols == 1 || NrOfDMRSSymbols == 2)//Single Front Loaded.
        {
            if(MappingType == 0)//Mapping Type A.
            {
                startSym[0] = 0;
                endSym[0] = numOfSymbs-1;
            }
            else//Mapping Type B
            {
                startSym[0] = 0;
                endSym[0] = numOfSymbs-1;
            }
        }*/
    }
    else if(*nSeg == 2)
    {
        if(NrOfDMRSSymbols == 1)//Single Front Loaded.
        {
            if(MappingType == 0)//Mapping Type A
            {
                if(numOfSymbs <= 9)
                {
                    startSym[0] = 0;endSym[0] = 4;
                    startSym[1] = 5;endSym[1] = numOfSymbs-1;
                }
                else if(numOfSymbs <= 12)
                {
                    startSym[0] = 0;endSym[0] = 5;
                    startSym[1] = 6;endSym[1] = numOfSymbs-1;
                }
                else if(numOfSymbs <= 14)
                {
                    startSym[0] = 0;endSym[0] = 6;
                    startSym[1] = 7;endSym[1] = numOfSymbs-1;
                }
            }
            else//Mapping Type B
            {
                if(numOfSymbs <= 7)
                {
                    startSym[0] = 0;endSym[0] = 2;
                    startSym[1] = 3;endSym[1] = numOfSymbs-1;
                }
                else if(numOfSymbs <= 9)
                {
                    startSym[0] = 0;endSym[0] = 3;
                    startSym[1] = 4;endSym[1] = numOfSymbs-1;
                }
                else if(numOfSymbs <= 11)
                {
                    startSym[0] = 0;endSym[0] = 4;
                    startSym[1] = 5;endSym[1] = numOfSymbs-1;
                }
                else if(numOfSymbs <= 14)
                {
                    startSym[0] = 0;endSym[0] = 5;
                    startSym[1] = 6;endSym[1] = numOfSymbs-1;
                }
            }
        }
        else //if(NrOfDMRSSymbols == 2)//Double Front Loaded.
        {
            if(MappingType == 0)//Mapping Type A
            {
                if(numOfSymbs <= 12)
                {
                    startSym[0] = 0;endSym[0] = 5;
                    startSym[1] = 6;endSym[1] = numOfSymbs-1;
                }
                else //if(numOfSymbs <= 14)
                {
                    startSym[0] = 0;endSym[0] = 6;
                    startSym[1] = 7;endSym[1] = numOfSymbs-1;
                }
            }
            else//Mapping Type B
            {
                if(numOfSymbs <= 9)
                {
                    startSym[0] = 0;endSym[0] = 3;
                    startSym[1] = 4;endSym[1] = numOfSymbs-1;
                }
                else if(numOfSymbs <= 11)
                {
                    startSym[0] = 0;endSym[0] = 4;
                    startSym[1] = 5;endSym[1] = numOfSymbs-1;
                }
                else //if(numOfSymbs <= 14)
                {
                    startSym[0] = 0;endSym[0] = 5;
                    startSym[1] = 6;endSym[1] = numOfSymbs-1;
                }
            }
        }
    }
    else if(*nSeg == 3)
    {
        if(NrOfDMRSSymbols == 1)//Single Front Loaded.
        {
            if(MappingType == 0)//Mapping Type A
            {
                if(numOfSymbs <= 12)
                {
                    if(pNrPuschInParams->dmrs1_sym == 2)
                    {
                        startSym[0] = 0;endSym[0] = 3;
                        startSym[1] = 4;endSym[1] = 7;
                        startSym[2] = 8;endSym[2] = numOfSymbs-1;
                    }
                    else if(pNrPuschInParams->dmrs1_sym == 3)
                    {
                        startSym[0] = 0;endSym[0] = 4;
                        startSym[1] = 5;endSym[1] = 7;
                        startSym[2] = 8;endSym[2] = numOfSymbs-1;
                    }
                }
                else if(numOfSymbs <= 14)
                {
                    startSym[0] = 0;endSym[0] = 4;
                    startSym[1] = 5;endSym[1] = 8;
                    startSym[2] = 9;endSym[2] = numOfSymbs-1;
                }
            }
            else // Mapping Type B
            {
                if(numOfSymbs <= 9)
                {
                    startSym[0] = 0;endSym[0] = 1;
                    startSym[1] = 2;endSym[1] = 4;
                    startSym[2] = 5;endSym[2] = numOfSymbs-1;
                }
                else if(numOfSymbs <= 11)
                {
                    startSym[0] = 0;endSym[0] = 2;
                    startSym[1] = 3;endSym[1] = 6;
                    startSym[2] = 7;endSym[2] = numOfSymbs-1; 
                }
                else if(numOfSymbs <= 14)
                {
                    startSym[0] = 0;endSym[0] = 3;
                    startSym[1] = 4;endSym[1] = 7;
                    startSym[2] = 8;endSym[2] = numOfSymbs-1;
                }
            }
            
        }
    }
    else if(*nSeg == 4)
    {
        if(NrOfDMRSSymbols == 1)//Single Front Loaded.
        {
            if(MappingType == 0)//Mapping Type A
            {
                startSym[0] = 0;endSym[0] = 3;
                startSym[1] = 4;endSym[1] = 6;
                startSym[2] = 7;endSym[2] = 9;
                startSym[3] = 10;endSym[3] = numOfSymbs-1;
            }
            else// Mapping Type B.
            {
                startSym[0] = 0;endSym[0] = 1;
                startSym[1] = 2;endSym[1] = 4;
                startSym[2] = 5;endSym[2] = 7;
                startSym[3] = 8;endSym[3] = numOfSymbs-1;
            }
        }
    }
}
    

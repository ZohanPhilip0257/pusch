#ifndef WNNRPHYPUSCHHEADER_H
#define WNNRPHYPUSCHHEADER_H

#include "wnNrPhyBwFftCpConfig.h"
#include "wnNrPhyPuschChnlEstEqHeader.h"
#include "immintrin.h"
#include "string.h"
#include<stdio.h>
#include<math.h>
#include<stdint.h>

#define TOTAL_RES_IN_SYMBOL                 (MAX_SUB_CARRIERS_OFDM_SYMBOL)
#define MAX_PUSCH_DMRS_PER_SYMBOL           (MAX_SUB_CARRIERS_OFDM_SYMBOL/2)
#define MAX_PUSCH_INFO_PER_SYMBOL           (MAX_SUB_CARRIERS_OFDM_SYMBOL)
#define MAX_PUSCH_DATA                      (1270000)
#define MAX_LDPC_DECODER_INPUT_LENGTH       (26112)
#define MAX_CRC_INPUT_LENGTH_IN_4BYTES      (265)
#define MAX_PUSCH_CB                        (150)
#define AVX_2                               (1)        // For LDPC
#define CHECK_PARITY                        (0)        // For LDPC
#define BW_L                                (6)        // For LDPC
#define OFFSET                              (0)        // For LDPC  
#define MAX_ITERATIONS                      (32)        // For LDPC
#define FILLER_BIT                          (120)  
#define ALIGN_BOUNDARY 	                    (32)
#define x                                   (200)

#define PER_PRB_EST_EQU         ( 1 )
#define PER_HALF_PRB_EST_EQU    ( 2 )
#define PER_TONE_EST_EQU        ( 3 )

//cmplx pusch_dmrs_generated[NUMBER_OF_OFDM_SYMBOLS_SLOT][TOTAL_RES_IN_SYMBOL][MAX_SUMIMO_PUSCH_ANTENNA_PORTS];
wnFlt rxFdSamples[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][NUMBER_OF_OFDM_SYMBOLS_SLOT][2 * TOTAL_RES_IN_SYMBOL];
wnFlt pusch_dmrs_generated[MAX_SUMIMO_PUSCH_ANTENNA_PORTS][NUMBER_OF_OFDM_SYMBOLS_SLOT][TOTAL_RES_IN_SYMBOL*2];
__attribute__ ((aligned(ALIGN_BOUNDARY))) wnInt8 llrScaledOut[MAX_PUSCH_DATA];
__attribute__ ((aligned(ALIGN_BOUNDARY))) wnFlt softDemodOut[MAX_PUSCH_DATA];
__attribute__ ((aligned(ALIGN_BOUNDARY))) wnInt8 scramblingOut[MAX_PUSCH_DATA];
__attribute__ ((aligned(ALIGN_BOUNDARY))) wnInt8 ldpcIn[MAX_LDPC_DECODER_INPUT_LENGTH]; // max(nbg1, nbg2)*max(Zc) = max(68,52)*384 = 68*384 = 26112   
__attribute__ ((aligned(ALIGN_BOUNDARY))) wnInt8 Ltot[MAX_LDPC_DECODER_INPUT_LENGTH]; // max(nbg1, nbg2)*max(Zc) = max(68,52)*384 = 68*384 = 26112
__attribute__ ((aligned(ALIGN_BOUNDARY))) wnInt8 Lreg[20][384];
wnUInt8 dmrsCdmg[3], dmrs_pr_prb, dmrsCdmGWithData, dmrsCdmGWithoutData[3], dmrs_symbols, dmrs_s[16], nPtrsSyms;
wnUInt32 nPrbAllocated;
__attribute__ ((aligned(ALIGN_BOUNDARY))) wnUInt8 ldpcOut[MAX_PUSCH_DATA];
__attribute__ ((aligned(32))) wnInt8 bitSelectIn[MAX_PUSCH_DATA];
__attribute__ ((aligned(ALIGN_BOUNDARY))) wnUInt32 crcIn[MAX_CRC_INPUT_LENGTH_IN_4BYTES];
extern const wnInt32 baseGraph[102][3128];
wnUInt8 pnSeqOut[MAX_PUSCH_DATA];
wnUInt8 puschData[MAX_PUSCH_DATA];
wnUInt8 cbSeg[MAX_PUSCH_DATA];
wnUInt8 puschDecodedData[MAX_PUSCH_DATA];
wnUInt8 nSeg;
wnUInt8 startSym[4], endSym[4];

// CRC_16 Table
extern const wnUInt16 CRC16_TAB[256];

// CRC_24A Table
extern const wnUInt32 CRC24A_TAB[256];

// CRC_24B Table
extern const wnUInt32 CRC24B_TAB[256];

// CRC_24C Table
extern const wnUInt32 CRC24C_TAB[256];

// Used for Descrambling
extern const wnInt8 scramblerConvertor[256][8];

// Used to find Zc
extern const wnUInt32 Zc_all[51];

// Used to find Base Graph
extern const wnUInt32 iLS_all[51];
extern const wnUInt32 increment[8];
extern const wnUInt32 indexs[8];

// qpsk Modulation
extern const wnFlt qpsk[256][8];

// 16 QAM Modulation
extern const wnFlt q16[256][4];

// 64 QAM Modulation
extern const wnFlt q64[64][2];

// 256 QAM Modulation
extern const wnFlt q256[256][2];

// To access this :: pusch_dmrs_config01_cdm_group(mod(port_number,1000))
extern const wnUInt8 pusch_dmrs_config01_cdm_group[8];

// To access this :: pusch_dmrs_config01_delta(mod(port_number,1000))
extern const wnUInt8 pusch_dmrs_config01_delta[8];

// To access this :: pusch_dmrs_config01_wf_k(mod(port_number,1000),k_desh)
extern const wnInt8 pusch_dmrs_config01_wf_k[2][8];

// To access this :: pusch_dmrs_config01_wt_l(mod(port_number,1000), l_desh)
extern const wnInt8 pusch_dmrs_config01_wt_l[2][8];

// To access this :: pusch_dmrs_config01_cdm_group(mod(port_number,1000))
extern const wnUInt8 pusch_dmrs_config02_cdm_group[12];

// To access this :: pusch_dmrs_config01_delta(mod(port_number,1000))
extern const wnUInt8 pusch_dmrs_config02_delta[12];

// To access this :: pusch_dmrs_config01_wf_k(mod(port_number,1000),k_desh)
extern const wnInt8 pusch_dmrs_config02_wf_k[2][12];

// To access this :: pusch_dmrs_config01_wt_l(mod(port_number,1000), l_desh)
extern const wnInt8 pusch_dmrs_config02_wt_l[2][12];

// Used to determine DMRS Symbol Locations in case of Single Symbol per slot
extern const wnUInt8 puschSingleSymTbl[14][8][5];

// Used to determine DMRS Symbol Locations in case of Double Symbol per slot
extern const wnUInt8 puschDoubleSymTbl[14][8][3];

// Table 5.1.3.2-2: TBS for Ninfo <= 3824
extern const wnUInt32 tbs_table[93];

// Used in PRS Generation
extern const wnUInt32 SEQUENCE_X2_MASKS[31];


wnVoid wnNrPhyPuschPnSequence(wnInt64 seedVal, wnUInt32 seqLen, wnUInt8* pnSeqOut);
wnVoid wnNrPhyPuschDmrs(puschConfigT* nrUlPhyPuschPDUs);
wnVoid wnNrPhyPuschdemapping(puschConfigT* nrUlPhyPuschPDUs);
wnVoid wnSegBoundaries(wnUInt8 *nSeg, wnUInt8* startSym, wnUInt8* endSym, puschConfigT* pNrPuschInParams);
wnVoid wnNrPhyPuschInit(puschConfigT* commonVar);
wnVoid wnNrPhyPuschDescrambling (puschConfigT* nrUlPhyPuschPDUs, wnInt8 *scramblingOut, wnInt8 *scramblingIn, wnUInt32 seqLength);
wnVoid wnNrPhyPuschCwProc (puschConfigT* nrUlPhyPuschPDUs, wnInt8 *layerMappedOut, wnUInt32 rmOutLength);

wnVoid wnNrPhyPuschSetup (puschConfigT* nrUlPhyPuschPDUs, wnUInt32 *rmOutLength, wnUInt8 *m_bg, wnUInt8 *n_bg, wnUInt32 *eachCbDataLength, wnUInt8 *par_start, wnUInt8 *Kb, wnUInt32 *Zc, wnUInt32 *n_parity, wnUInt32 *Pb, wnUInt32 *n_all, wnUInt8 *L_max_int, wnUInt8 *Ltot_max_int, wnUInt8 *offset, wnUInt32 *baseGraphIndex, wnUInt8 *maxItrs, wnUInt32 *nCodeBlocks, wnUInt32 *N_cb, wnUInt32 *start_pos_bit_select, wnUInt32 *fillerStart, wnUInt32 *fillerEnd, wnUInt32 *puschRmLenCb);
wnVoid wnNrPhyPuschLdpcDecoderSetup (wnUInt32 *row_wt, wnUInt32 *col_wt, wnUInt32 *r_nz, wnUInt32 *c_nz, wnUInt32 *sh_nz, wnUInt8 m_bg, wnUInt8 n_bg, const wnInt32 H_bg[][n_bg]);
wnVoid nrldpc_decoder(wnInt8 *L0, wnInt8 *Ltot, wnInt8 Lreg[20][384], wnInt8 *dec_cword, wnUInt32 *row_wt, wnUInt32 *c_nz, wnUInt32 *sh_nz, wnUInt32 Zc, wnUInt32 msg_len, wnUInt8 Kb, wnUInt32 Pb, wnUInt8 par_start, wnUInt32 nparity, wnUInt32 n_all, wnUInt8 maxitrs, wnUInt8 offset, wnUInt8 L_max, wnUInt8 Ltot_max);
wnVoid wnNrPhyCrc16Removal(wnUInt8 *payload, wnInt32 numBits, wnUInt8 *crcRemoved, wnUInt8 *check);
wnVoid wnNrPhyCrc24ARemoval(wnUInt8 *payload, wnInt32 numBits, wnUInt8 *crcRemoved, wnUInt8 *check);
wnVoid wnNrPhyCrc24BRemoval(wnUInt8 *payload, wnInt32 numBits, wnUInt8 *crcRemoved, wnUInt8 *check);
#endif
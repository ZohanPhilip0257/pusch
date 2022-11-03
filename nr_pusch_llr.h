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
  * @brief This file includes all the function prototypes
  *         for NR PUSCH Channel Estimation and Equalization.
  * @file nr_pusch.h
  * @ingroup nr_pusch
  * @author WiSig Networks
  **/

#ifndef NR_PUSCH_LLR_H
#define NR_PUSCH_LLR_H

// #include "common_def.h"
#include "nr_pusch_mmse_irc.h"
#include "nr_pusch_mmse.h"

void nr_pusch_llr_compute_avx2( wnUInt8 segPart,
                                wnUInt8* startSym,
                                wnUInt8* endSym,
                                commonUlConfigT* pNrUlCommonParams,
                                puschConfigT* pNrPuschInParams,
                                P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams );

void nr_pusch_llr_scaling_avx2( wnUInt8 segPart,
                                wnUInt8* startSym,
                                wnUInt8* endSym,
                                commonUlConfigT* pNrUlCommonParams,
                                puschConfigT* pNrPuschInParams,
                                P_NR_PUSCH_OUT_PARAMS pNrPuschOutParams,
                                wnUInt8 bwl,
                                wnFlt Lmax );

wnVoid wnNrPhyPuschCwProc (puschConfigT* nrUlPhyPuschPDUs,\
                            wnInt8 *layerMappedOut, wnUInt32 rmOutLength);

#endif // NR_PUSCH_LLR_H

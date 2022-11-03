#include "wnNrPhyPuschHeader.h"
#include "nr_pusch_mmse_irc.h"
#include "nr_pusch_mmse.h"
#include "nr_pusch_llr.h"
#define SNR24
#include <stdio.h>
#include <stdlib.h>
#include "wnNrExeTimeX86Intel.h"


wnVoid wnNrPhyPuschMain(puschConfigT* pNrPuschInParams,
                        commonUlConfigT* pNrUlCommonParams)
{ 
	// Generates DMRS Data for the Channel Estimation Function
	wnNrPhyPuschDmrs(pNrPuschInParams);
    
    
    #ifdef DEBUG
        FILE *fp, *pf, *fp1, *pf1;
        wnUInt32 dmrsSyms;
        if (pNrPuschInParams->nDMRSConfigType == 0)
        {
            dmrsSyms = pNrPuschInParams->nRBSize*6;
        }
        else
        {
            dmrsSyms = pNrPuschInParams->nRBSize*4;
        }
        
        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
        {
            fp = fopen("IOs/DMRSReal1.txt", "w");
            pf = fopen("IOs/DMRSImag1.txt", "w");
            
            for(int i = 0;i<dmrsSyms;i++)
            {
                fprintf(fp, "%f\n", nrPuschOutParamsSeg1.genDMRS[0][i]);
                fprintf(pf, "%f\n", nrPuschOutParamsSeg1.genDMRS[0][275*6+i]);
            }
            fclose(fp);
            fclose(pf);
        }
        else
        {
            fp = fopen("IOs/DMRSReal1.txt", "w");
            pf = fopen("IOs/DMRSImag1.txt", "w");
            fp1 = fopen("IOs/DMRSReal2.txt", "w");
            pf1 = fopen("IOs/DMRSImag2.txt", "w");
            for(int i = 0;i<dmrsSyms;i++)
            {
                fprintf(fp, "%f\n", nrPuschOutParamsSeg1.genDMRS[0][i]);
                fprintf(pf, "%f\n", nrPuschOutParamsSeg1.genDMRS[0][275*6+i]);
                fprintf(fp1, "%f\n", nrPuschOutParamsSeg1.genDMRS[1][i]);
                fprintf(pf1, "%f\n", nrPuschOutParamsSeg1.genDMRS[1][275*6+i]);
            }
            fclose(fp);
            fclose(pf);
            fclose(fp1);
            fclose(pf1);
        }
    #endif
	
    //Remove location offset in Time and Frequency domain.
    wnNrPhyPuschdemapping(pNrPuschInParams);

    // Computes number of segments and symbol boundaries for each segment.
    wnSegBoundaries(&nSeg, startSym, endSym, pNrPuschInParams);

    // printf("nSeg: %u\n", nSeg);
    // printf("startSym: %u\n", startSym[0]);
    // printf("endSym: %u\n", endSym[0]);
    
    wnUInt8 bwl = 6;
    wnFlt Lmax = pNrPuschInParams->Lmax;

    if(pNrPuschInParams->irc_enable)
    {
        switch(pNrPuschInParams->equType)
        {
            case PER_PRB_EST_EQU:
                /*nr_pusch_mmse_irc_perPrb_est_avx2(LOW_PUSCH_SYM,\
                //                                   pNrUlCommonParams,\
                //                                   pNrPuschInParams,\
                //                                   &nrPuschOutParamsLow);
				// nr_pusch_mmse_irc_perPrb_est_avx2(HIGH_PUSCH_SYM,\
                //                                   pNrUlCommonParams,\
                //                                   pNrPuschInParams,\
                //                                   &nrPuschOutParamsHigh);
                // nr_pusch_mmse_irc_perPrb_equ_avx2(LOW_PUSCH_SYM,\
                //                                   pNrUlCommonParams,\
                //                                   pNrPuschInParams,\
                //                                   &nrPuschOutParamsLow);
				// nr_pusch_mmse_irc_perPrb_equ_avx2(HIGH_PUSCH_SYM,\
                //                                   pNrUlCommonParams,\
                //                                   pNrPuschInParams,\
                //                                   &nrPuschOutParamsHigh);*/
            break;

            case PER_HALF_PRB_EST_EQU:
                // Yet to recieve code from Sami.
            break;

            case PER_TONE_EST_EQU:
                /*nr_pusch_mmse_irc_perTone_est_avx2(LOW_PUSCH_SYM,\
                //                                    pNrUlCommonParams,\
                //                                    pNrPuschInParams,\
                //                                    &nrPuschOutParamsLow);
				// nr_pusch_mmse_irc_perTone_est_avx2(HIGH_PUSCH_SYM,
                //                                    pNrUlCommonParams,\
                //                                    pNrPuschInParams,\
                //                                    &nrPuschOutParamsHigh);
                // nr_pusch_mmse_irc_perTone_equ_avx2(LOW_PUSCH_SYM,\
                //                                    pNrUlCommonParams,\
                //                                    pNrPuschInParams,\
                //                                    &nrPuschOutParamsLow);
				// nr_pusch_mmse_irc_perTone_equ_avx2(HIGH_PUSCH_SYM,\
                //                                    pNrUlCommonParams,\
                //                                    pNrPuschInParams,\
                //                                    &nrPuschOutParamsHigh);*/
            break;
        }
    }
    else
    {
        switch(pNrPuschInParams->equType)
        {
            case PER_PRB_EST_EQU:
                /* nr_pusch_mmse_perPrb_est_avx2(LOW_PUSCH_SYM,\
                //                               pNrUlCommonParams,\
                //                               pNrPuschInParams,\
                //                               &nrPuschOutParamsLow);
				// nr_pusch_mmse_perPrb_est_avx2(HIGH_PUSCH_SYM,\
                //                               pNrUlCommonParams,\
                //                               pNrPuschInParams,\
                //                               &nrPuschOutParamsHigh);
                // nr_pusch_mmse_perPrb_equ_avx2(LOW_PUSCH_SYM,\
                //                             pNrUlCommonParams,\
                //                             pNrPuschInParams,\
                //                               &nrPuschOutParamsLow);
				// nr_pusch_mmse_perPrb_equ_avx2(HIGH_PUSCH_SYM,\
                //                               pNrUlCommonParams,\
                //                               pNrPuschInParams,\
                //                               &nrPuschOutParamsHigh);*/
            break;

            case PER_HALF_PRB_EST_EQU:
                /* nr_pusch_mmse_perHPrb_est_avx2(LOW_PUSCH_SYM,\
                //                               pNrUlCommonParams,\
                //                               pNrPuschInParams,\
                //                               &nrPuschOutParamsLow);
				// nr_pusch_mmse_perHPrb_est_avx2(HIGH_PUSCH_SYM,\
                //                               pNrUlCommonParams,\
                //                               pNrPuschInParams,\
                //                               &nrPuschOutParamsHigh);
                // nr_pusch_mmse_perHPrb_equ_avx2(LOW_PUSCH_SYM,\
                //                             pNrUlCommonParams,\
                //                             pNrPuschInParams,\
                //                               &nrPuschOutParamsLow);
				// nr_pusch_mmse_perHPrb_equ_avx2(HIGH_PUSCH_SYM,\
                //                               pNrUlCommonParams,\
                //                               pNrPuschInParams,\
                //                               &nrPuschOutParamsHigh);*/
            break;

            case PER_TONE_EST_EQU:
                if(nSeg == 1)
                {
                    if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
                    {
                        // printf("test 1\n");
                        nr_pusch_mmse_perTone_est_avx2(SEG1,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg1);
                        // printf("estimation finished\n");
                        nr_pusch_mmse_perTone_equ_avx2(nSeg,\
                                                SEG1,\
                                                startSym,\
                                                endSym,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg1);//*/
                    }
                    else
                    {
                        nr_pusch_mmse_perTone_DS_est_avx2(SEG1,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg1);

                        nr_pusch_mmse_perTone_equ_avx2(nSeg,\
                                                    SEG1,\
                                                    startSym,\
                                                    endSym,\
                                                    pNrUlCommonParams,\
                                                    pNrPuschInParams,\
                                                    &nrPuschOutParamsSeg1);
                    }

                    //If trasform precoding is enabled, perform IFFT on each data symbol.
                    if(pNrPuschInParams->transformPrecode == 1)
                    {
                        nr_pusch_mmse_TF_IFFT(nSeg,
                                              SEG1,
                                              startSym,
                                              endSym,
                                              pNrPuschInParams,
                                              &nrPuschOutParamsSeg1);
                    }

                    #ifdef DEBUG
                        wnUInt8 nSymbs;
                        wnUInt8 nLyr = pNrPuschInParams->nNrOfLayers;
                        wnUInt32 Idx;

                        fp = fopen("IOs/layerDemapperReal_mmse_perToneSeg1.txt", "w");
                        pf = fopen("IOs/layerDemapperImag_mmse_perToneSeg1.txt", "w");

                        
                        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
                        {
                            nSymbs = endSym[0] - startSym[0];
                        }
                        else
                        {
                            nSymbs = endSym[0] - startSym[0] - 1;
                        }
                        
                        for(Idx = 0;Idx< nrPuschOutParamsSeg1.total_nPRB*12*nLyr*nSymbs;Idx++)
                        {
                            fprintf(fp, "%f\n", nrPuschOutParamsSeg1.layerDemapperOutReal[Idx]);
                            fprintf(pf, "%f\n", nrPuschOutParamsSeg1.layerDemapperOutImag[Idx]);
                        }

                        fclose(fp);
                        fclose(pf);
                    #endif

                    // Soft demodulation.
                    nr_pusch_llr_compute_avx2(SEG1,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg1);  

                    nr_pusch_llr_scaling_avx2(SEG1,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg1,\
                                              bwl,\
                                              Lmax);
                }
                else if(nSeg == 2)
                {
                    if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
                    {
                        nr_pusch_mmse_perTone_est_avx2(SEG1,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg1);

                        nr_pusch_mmse_perTone_est_avx2(SEG2,\
                                                    pNrUlCommonParams,\
                                                    pNrPuschInParams,\
                                                    &nrPuschOutParamsSeg2);

                        nr_pusch_mmse_perTone_equ_avx2(nSeg,\
                                                    SEG1,\
                                                    startSym,\
                                                    endSym,\
                                                    pNrUlCommonParams,\
                                                    pNrPuschInParams,\
                                                    &nrPuschOutParamsSeg1);

                        nr_pusch_mmse_perTone_equ_avx2(nSeg,\
                                                    SEG2,\
                                                    startSym,\
                                                    endSym,\
                                                    pNrUlCommonParams,\
                                                    pNrPuschInParams,\
                                                    &nrPuschOutParamsSeg2);
                    }
                    else
                    {
                        nr_pusch_mmse_perTone_DS_est_avx2(SEG1,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg1);

                        nr_pusch_mmse_perTone_DS_est_avx2(SEG2,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg2);

                        nr_pusch_mmse_perTone_equ_avx2(nSeg,\
                                                    SEG1,\
                                                    startSym,\
                                                    endSym,\
                                                    pNrUlCommonParams,\
                                                    pNrPuschInParams,\
                                                    &nrPuschOutParamsSeg1);

                        nr_pusch_mmse_perTone_equ_avx2(nSeg,\
                                                    SEG2,\
                                                    startSym,\
                                                    endSym,\
                                                    pNrUlCommonParams,\
                                                    pNrPuschInParams,\
                                                    &nrPuschOutParamsSeg2);
                    }

                    //If trasform precoding is enabled, perform IFFT on each data symbol.
                    if(pNrPuschInParams->transformPrecode == 1)
                    {
                        nr_pusch_mmse_TF_IFFT(nSeg,
                                              SEG1,
                                              startSym,
                                              endSym,
                                              pNrPuschInParams,
                                              &nrPuschOutParamsSeg1);

                        nr_pusch_mmse_TF_IFFT(nSeg,
                                              SEG2,
                                              startSym,
                                              endSym,
                                              pNrPuschInParams,
                                              &nrPuschOutParamsSeg2);
                    }

                    #ifdef DEBUG
                        // FILE *fp, *pf;
                        wnUInt8 nSymbs;
                        wnUInt8 nLyr = pNrPuschInParams->nNrOfLayers;
                        wnUInt32 Idx;

                        fp = fopen("IOs/layerDemapperReal_mmse_perToneSeg1.txt", "w");
                        pf = fopen("IOs/layerDemapperImag_mmse_perToneSeg1.txt", "w");
                        
                        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
                        {
                            nSymbs = endSym[0] - startSym[0];
                        }
                        else
                        {
                            nSymbs = endSym[0] - startSym[0] - 1;
                        }
                        
                        for(Idx = 0;Idx< nrPuschOutParamsSeg1.total_nPRB*12*nLyr*nSymbs;Idx++)
                        {
                            fprintf(fp, "%f\n", nrPuschOutParamsSeg1.layerDemapperOutReal[Idx]);
                            fprintf(pf, "%f\n", nrPuschOutParamsSeg1.layerDemapperOutImag[Idx]);
                        }

                        fclose(fp);
                        fclose(pf);

                        fp = fopen("IOs/layerDemapperReal_mmse_perToneSeg2.txt", "w");
                        pf = fopen("IOs/layerDemapperImag_mmse_perToneSeg2.txt", "w");
                        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
                        {
                            nSymbs = endSym[1] - startSym[1];
                        }
                        else
                        {
                            nSymbs = endSym[1] - startSym[1] - 1;
                        }
                        
                        for(Idx = 0;Idx< nrPuschOutParamsSeg2.total_nPRB*12*nLyr*nSymbs;Idx++)
                        {
                            fprintf(fp, "%f\n", nrPuschOutParamsSeg2.layerDemapperOutReal[Idx]);
                            fprintf(pf, "%f\n", nrPuschOutParamsSeg2.layerDemapperOutImag[Idx]);
                        }
                        fclose(fp);
                        fclose(pf);
                    #endif
                    
                    if(pNrPuschInParams->FOEnable)
                    {
                        #ifdef PROFILE_GNB
                        for (int iterIdx = 0; iterIdx < PROF_NUM_OF_ITER; iterIdx++)
                        {
                            /* To get the current time */
                            clock_gettime(CLOCK_REALTIME, &estStTime);
                        #endif

                        nr_pusch_mmse_perTone_FO_est(pNrPuschInParams,\
                                                    &nrPuschOutParamsSeg1,\
                                                    &nrPuschOutParamsSeg2);

                        #ifdef PROFILE_GNB
                        /* To get the current time */
                            clock_gettime(CLOCK_REALTIME, &estEndTime);
                            exeTimeEst[iterIdx] = ((estEndTime.tv_sec - estStTime.tv_sec) * 1E9 + (estEndTime.tv_nsec - estStTime.tv_nsec));
    
                        }
                        #endif

                        #ifdef PROFILE_GNB
                        for (int iterIdx = 0; iterIdx < PROF_NUM_OF_ITER; iterIdx++)
                        {
                            /* To get the current time */
                            clock_gettime(CLOCK_REALTIME, &corrStTime);
                        #endif
                        
                        nr_pusch_mmse_perTone_FO_Corr(SEG1,\
                                                    startSym,\
                                                    endSym,\
                                                    pNrUlCommonParams,\
                                                    pNrPuschInParams,\
                                                    &nrPuschOutParamsSeg1);

                        nr_pusch_mmse_perTone_FO_Corr(SEG2,\
                                                    startSym,\
                                                    endSym,\
                                                    pNrUlCommonParams,\
                                                    pNrPuschInParams,\
                                                    &nrPuschOutParamsSeg2);
                        #ifdef PROFILE_GNB
                        /* To get the current time */
                            clock_gettime(CLOCK_REALTIME, &corrEndTime);
                            exeTimeCorr[iterIdx] = ((corrEndTime.tv_sec - corrStTime.tv_sec) * 1E9 + (corrEndTime.tv_nsec - corrStTime.tv_nsec));
                        }
                        wnNrExeTimeX86Intel(PROF_NUM_OF_ITER);
                        #endif
                    }    

                    nr_pusch_llr_compute_avx2(SEG1,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg1);

                    nr_pusch_llr_compute_avx2(SEG2,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg2);   

                    nr_pusch_llr_scaling_avx2(SEG1,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg1,\
                                              bwl,\
                                              Lmax);

                    nr_pusch_llr_scaling_avx2(SEG2,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg2,\
                                              bwl,\
                                              Lmax);   

                }
                else if(nSeg == 3)       
                {
                    nr_pusch_mmse_perTone_est_avx2(SEG1,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg1);

                    nr_pusch_mmse_perTone_est_avx2(SEG2,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg2);

                    nr_pusch_mmse_perTone_est_avx2(SEG3,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg3);

                    nr_pusch_mmse_perTone_equ_avx2(nSeg,\
                                                SEG1,\
                                                startSym,\
                                                endSym,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg1);

                    nr_pusch_mmse_perTone_equ_avx2(nSeg,\
                                                SEG2,\
                                                startSym,\
                                                endSym,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg2);

                    nr_pusch_mmse_perTone_equ_avx2(nSeg,\
                                                SEG3,\
                                                startSym,\
                                                endSym,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg3);

                    #ifdef DEBUG
                        wnUInt8 nSymbs;
                        wnUInt8 nLyr = pNrPuschInParams->nNrOfLayers;
                        wnUInt32 Idx;

                        fp = fopen("IOs/layerDemapperReal_mmse_perToneSeg1.txt", "w");
                        pf = fopen("IOs/layerDemapperImag_mmse_perToneSeg1.txt", "w");
                        
                        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
                        {
                            nSymbs = endSym[0] - startSym[0];
                        }
                        else
                        {
                            nSymbs = endSym[0] - startSym[0] - 1;
                        }
                        
                        for(Idx = 0;Idx< nrPuschOutParamsSeg1.total_nPRB*12*nLyr*nSymbs;Idx++)
                        {
                            fprintf(fp, "%f\n", nrPuschOutParamsSeg1.layerDemapperOutReal[Idx]);
                            fprintf(pf, "%f\n", nrPuschOutParamsSeg1.layerDemapperOutImag[Idx]);
                        }

                        fclose(fp);
                        fclose(pf);

                        fp = fopen("IOs/layerDemapperReal_mmse_perToneSeg2.txt", "w");
                        pf = fopen("IOs/layerDemapperImag_mmse_perToneSeg2.txt", "w");
                        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
                        {
                            nSymbs = endSym[1] - startSym[1];
                        }
                        else
                        {
                            nSymbs = endSym[1] - startSym[1] - 1;
                        }
                        
                        for(Idx = 0;Idx< nrPuschOutParamsSeg2.total_nPRB*12*nLyr*nSymbs;Idx++)
                        {
                            fprintf(fp, "%f\n", nrPuschOutParamsSeg2.layerDemapperOutReal[Idx]);
                            fprintf(pf, "%f\n", nrPuschOutParamsSeg2.layerDemapperOutImag[Idx]);
                        }
                        fclose(fp);
                        fclose(pf);

                        fp = fopen("IOs/layerDemapperReal_mmse_perToneSeg3.txt", "w");
                        pf = fopen("IOs/layerDemapperImag_mmse_perToneSeg3.txt", "w");
                        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
                        {
                            nSymbs = endSym[2] - startSym[2];
                        }
                        else
                        {
                            nSymbs = endSym[2] - startSym[2] - 1;
                        }
                        
                        for(Idx = 0;Idx< nrPuschOutParamsSeg3.total_nPRB*12*nLyr*nSymbs;Idx++)
                        {
                            fprintf(fp, "%f\n", nrPuschOutParamsSeg3.layerDemapperOutReal[Idx]);
                            fprintf(pf, "%f\n", nrPuschOutParamsSeg3.layerDemapperOutImag[Idx]);
                        }
                        fclose(fp);
                        fclose(pf);
                    #endif

                    //If trasform precoding is enabled, perform IFFT on each data symbol.
                    //Note: With valid PDU, #layers with transform precoding enabled will be always equal to 1.
                    if(pNrPuschInParams->transformPrecode == 1 && pNrPuschInParams->nNrOfLayers == 1)
                    {
                        nr_pusch_mmse_TF_IFFT(nSeg,
                                              SEG1,
                                              startSym,
                                              endSym,
                                              pNrPuschInParams,
                                              &nrPuschOutParamsSeg1);

                        nr_pusch_mmse_TF_IFFT(nSeg,
                                              SEG2,
                                              startSym,
                                              endSym,
                                              pNrPuschInParams,
                                              &nrPuschOutParamsSeg2);

                        nr_pusch_mmse_TF_IFFT(nSeg,
                                              SEG3,
                                              startSym,
                                              endSym,
                                              pNrPuschInParams,
                                              &nrPuschOutParamsSeg3);
                    }

                    nr_pusch_llr_compute_avx2(SEG1,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg1);

                    nr_pusch_llr_compute_avx2(SEG2,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg2);

                    nr_pusch_llr_compute_avx2(SEG3,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg3);   

                    nr_pusch_llr_scaling_avx2(SEG1,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg1,\
                                              bwl,\
                                              Lmax);

                    nr_pusch_llr_scaling_avx2(SEG2,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg2,\
                                              bwl,\
                                              Lmax); 

                    nr_pusch_llr_scaling_avx2(SEG3,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg3,\
                                              bwl,\
                                              Lmax); 
                }
                else if(nSeg == 4)       
                {
                    nr_pusch_mmse_perTone_est_avx2(SEG1,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg1);

                    nr_pusch_mmse_perTone_est_avx2(SEG2,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg2);

                    nr_pusch_mmse_perTone_est_avx2(SEG3,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg3);

                    nr_pusch_mmse_perTone_est_avx2(SEG4,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg4);

                     nr_pusch_mmse_perTone_equ_avx2(nSeg,\
                                                SEG1,\
                                                startSym,\
                                                endSym,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg1);

                    nr_pusch_mmse_perTone_equ_avx2(nSeg,\
                                                SEG2,\
                                                startSym,\
                                                endSym,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg2);

                    nr_pusch_mmse_perTone_equ_avx2(nSeg,\
                                                SEG3,\
                                                startSym,\
                                                endSym,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg3);

                    nr_pusch_mmse_perTone_equ_avx2(nSeg,\
                                                SEG4,\
                                                startSym,\
                                                endSym,\
                                                pNrUlCommonParams,\
                                                pNrPuschInParams,\
                                                &nrPuschOutParamsSeg4);

                    #ifdef DEBUG
                        wnUInt8 nSymbs;
                        wnUInt8 nLyr = pNrPuschInParams->nNrOfLayers;
                        wnUInt32 Idx;

                        fp = fopen("IOs/layerDemapperReal_mmse_perToneSeg1.txt", "w");
                        pf = fopen("IOs/layerDemapperImag_mmse_perToneSeg1.txt", "w");
                        
                        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
                        {
                            nSymbs = endSym[0] - startSym[0];
                        }
                        else
                        {
                            nSymbs = endSym[0] - startSym[0] - 1;
                        }
                        
                        for(Idx = 0;Idx< nrPuschOutParamsSeg1.total_nPRB*12*nLyr*nSymbs;Idx++)
                        {
                            fprintf(fp, "%f\n", nrPuschOutParamsSeg1.layerDemapperOutReal[Idx]);
                            fprintf(pf, "%f\n", nrPuschOutParamsSeg1.layerDemapperOutImag[Idx]);
                        }

                        fclose(fp);
                        fclose(pf);

                        fp = fopen("IOs/layerDemapperReal_mmse_perToneSeg2.txt", "w");
                        pf = fopen("IOs/layerDemapperImag_mmse_perToneSeg2.txt", "w");
                        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
                        {
                            nSymbs = endSym[1] - startSym[1];
                        }
                        else
                        {
                            nSymbs = endSym[1] - startSym[1] - 1;
                        }
                        
                        for(Idx = 0;Idx< nrPuschOutParamsSeg2.total_nPRB*12*nLyr*nSymbs;Idx++)
                        {
                            fprintf(fp, "%f\n", nrPuschOutParamsSeg2.layerDemapperOutReal[Idx]);
                            fprintf(pf, "%f\n", nrPuschOutParamsSeg2.layerDemapperOutImag[Idx]);
                        }
                        fclose(fp);
                        fclose(pf);

                        fp = fopen("IOs/layerDemapperReal_mmse_perToneSeg3.txt", "w");
                        pf = fopen("IOs/layerDemapperImag_mmse_perToneSeg3.txt", "w");
                        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
                        {
                            nSymbs = endSym[2] - startSym[2];
                        }
                        else
                        {
                            nSymbs = endSym[2] - startSym[2] - 1;
                        }
                        
                        for(Idx = 0;Idx< nrPuschOutParamsSeg3.total_nPRB*12*nLyr*nSymbs;Idx++)
                        {
                            fprintf(fp, "%f\n", nrPuschOutParamsSeg3.layerDemapperOutReal[Idx]);
                            fprintf(pf, "%f\n", nrPuschOutParamsSeg3.layerDemapperOutImag[Idx]);
                        }
                        fclose(fp);
                        fclose(pf);

                        fp = fopen("IOs/layerDemapperReal_mmse_perToneSeg4.txt", "w");
                        pf = fopen("IOs/layerDemapperImag_mmse_perToneSeg4.txt", "w");
                        if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
                        {
                            nSymbs = endSym[3] - startSym[3];
                        }
                        else
                        {
                            nSymbs = endSym[3] - startSym[3] - 1;
                        }
                        
                        for(Idx = 0;Idx< nrPuschOutParamsSeg4.total_nPRB*12*nLyr*nSymbs;Idx++)
                        {
                            fprintf(fp, "%f\n", nrPuschOutParamsSeg4.layerDemapperOutReal[Idx]);
                            fprintf(pf, "%f\n", nrPuschOutParamsSeg4.layerDemapperOutImag[Idx]);
                        }
                        fclose(fp);
                        fclose(pf);
                    #endif

                    //If trasform precoding is enabled, perform IFFT on each data symbol.
                    //Note: With valid PDU, #layers with transform precoding enabled will be always equal to 1.
                    if(pNrPuschInParams->transformPrecode == 1 && pNrPuschInParams->nNrOfLayers == 1)
                    {
                        nr_pusch_mmse_TF_IFFT(nSeg,
                                              SEG1,
                                              startSym,
                                              endSym,
                                              pNrPuschInParams,
                                              &nrPuschOutParamsSeg1);

                        nr_pusch_mmse_TF_IFFT(nSeg,
                                              SEG2,
                                              startSym,
                                              endSym,
                                              pNrPuschInParams,
                                              &nrPuschOutParamsSeg2);

                        nr_pusch_mmse_TF_IFFT(nSeg,
                                              SEG3,
                                              startSym,
                                              endSym,
                                              pNrPuschInParams,
                                              &nrPuschOutParamsSeg3);

                        nr_pusch_mmse_TF_IFFT(nSeg,
                                              SEG4,
                                              startSym,
                                              endSym,
                                              pNrPuschInParams,
                                              &nrPuschOutParamsSeg4);
                    }

                    nr_pusch_llr_compute_avx2(SEG1,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg1);

                    nr_pusch_llr_compute_avx2(SEG2,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg2);

                    nr_pusch_llr_compute_avx2(SEG3,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg3);   

                    nr_pusch_llr_compute_avx2(SEG4,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg4);  

                    nr_pusch_llr_scaling_avx2(SEG1,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg1,\
                                              bwl,\
                                              Lmax);

                    nr_pusch_llr_scaling_avx2(SEG2,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg2,\
                                              bwl,\
                                              Lmax);

                    nr_pusch_llr_scaling_avx2(SEG3,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg3,\
                                              bwl,\
                                              Lmax);

                    nr_pusch_llr_scaling_avx2(SEG4,\
                                              startSym,\
                                              endSym,\
                                              pNrUlCommonParams,\
                                              pNrPuschInParams,\
                                              &nrPuschOutParamsSeg4,\
                                              bwl,\
                                              Lmax); 
                }
            break;
        }
    }

    wnUInt8 dataSyms = 0;
    if(pNrPuschInParams->nNrOfDMRSSymbols == 1)
    {
        for(wnUInt8 segIdx = 0;segIdx<nSeg;segIdx++)
        {
            dataSyms = dataSyms + (endSym[segIdx] - startSym[segIdx]); 
        }
    }
    else
    {
        for(wnUInt8 segIdx = 0;segIdx<nSeg;segIdx++)
        {
            dataSyms = dataSyms + (endSym[segIdx] - startSym[segIdx]-1); 
        }
    }
    
    
    // This variable will be helpful for mexing.
    totalLength = dataSyms*pNrPuschInParams->nRBSize*12*pNrPuschInParams->nNrOfLayers*pNrPuschInParams->modulationOrder;
    
    #ifdef DEBUG
        FILE *fpsoft = fopen("IOs/softDemodOut.txt", "w");
		for(int i = 0;i<(dataSyms*pNrPuschInParams->nRBSize*12*pNrPuschInParams->nNrOfLayers*pNrPuschInParams->modulationOrder);i++)
		{
			fprintf(fpsoft, "%f\n", llrOut[i]);
		}
		fclose(fpsoft);
    #endif         
    
    #ifdef DEBUG
        FILE *fpsoftFxd = fopen("IOs/LLRFxd.txt", "w");
		for(int i = 0;i<(dataSyms*pNrPuschInParams->nRBSize*12*pNrPuschInParams->nNrOfLayers*pNrPuschInParams->modulationOrder);i++)
		{
			fprintf(fpsoftFxd, "%hhd\n", llrFxd[i]);
		}
		fclose(fpsoftFxd);
    #endif

    wnUInt32 LLRLen = pNrPuschInParams->nRBSize*12*dataSyms*\
                      pNrPuschInParams->modulationOrder*pNrPuschInParams->nNrOfLayers;
    
    
    // Processes the output of LLR Scaling Module, and decodes
    // the PUSCH Data of the UE
    // DeScrambling to CRC Removal
	wnNrPhyPuschCwProc(pNrPuschInParams, llrFxd, LLRLen);
    return;
}
    

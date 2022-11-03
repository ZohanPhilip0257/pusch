#include <stdio.h>
#include "wnNrPhygNBUlChnlConfigInit.h"

wnVoid* wnNrPhygNBUlChnlConfigInit()
{
	// Initialising all the GNB Dl Common Configurations
	wnNrPhycommonUlConfigInit ();

    // Initialising all the gNB Common Tables
    // LUTs related to LDPC, CRC
    wnNrPhygNBCommonTables ();

	// Initialising IFFT CP UL Configurations of GNB
	wnNrPhyUlBwFftCpConfig(&gnbUlChnlSt.commonUlConfig);

	// Initalising PUSCH Data
	for (wnInt32 puschIdx = 0; puschIdx < 1; puschIdx++)
   		wnNrPhyPuschInit(&gnbUlChnlSt.puschConfig[puschIdx]);

	return NULL;
}


wnVoid wnNrPhycommonUlConfigInit()
{
	// Numerology
	gnbUlChnlSt.commonUlConfig.numerology = 1;

	// Allocated Band Width
	gnbUlChnlSt.commonUlConfig.allocatedBw = 50;

	// Operating Frequency Band
	gnbUlChnlSt.commonUlConfig.nrOperatingFreqBand = 0;

    // Slot Number
    gnbUlChnlSt.commonUlConfig.nSlotNumber = 0;

    // Number of Physical Antenna
    gnbUlChnlSt.commonUlConfig.nPhysicalAnt = 1;

	return;
}


wnVoid wnNrPhyPuschInit(puschConfigT* commonVar)
{
    /*************************** INITIALISATIONS *************************/
    /*********** Parameters Initialisation ****************/
    commonVar->allocatedBw = 50;
    commonVar->nBWPSize = 133;
    commonVar->nBWPStart = 0;
    commonVar->nRBStart = 15;
    commonVar->nRBSize = 6;

    commonVar->nMappingType = 0;
    commonVar->ulDmrsTypeAPos = 3;
    commonVar->nNrOfDMRSSymbols = 1;
    commonVar->nDMRSAddPos = 2;
    commonVar->nDMRSConfigType = 0;

    commonVar->nNrOfLayers = 1;
    commonVar->nAntennaPorts = gnbUlChnlSt.commonUlConfig.nPhysicalAnt;
    commonVar->nPortIndex[0] = 0;
    commonVar->nPortIndex[1] = 1; 
    commonVar->nPortIndex[2] = 2; 
    commonVar->nPortIndex[3] = 3; 
	commonVar->nPortIndex[4] = 0; 
	commonVar->nPortIndex[5] = 0; 
	commonVar->nPortIndex[6] = 0; 
	commonVar->nPortIndex[7] = 0;

    commonVar->nNrOfSymbols = 14;
    commonVar->nStartSymbolIndex = 0;//15
    
    commonVar->nDmrsWithoutData = 1;
    commonVar->numerology = 1;
    commonVar->nSlotNumber = 0;
    commonVar->nrOperatingFreqBand = 0;
    commonVar->nSCID = 0;
    commonVar->nNIDnSCID0 = 50;
    commonVar->nNIDnSCID1 = 50;

    commonVar->modulationOrder = 2;
    commonVar->nTbSize = 184;
    commonVar->codingRate = (float)120/1024;
    commonVar->nRNTI = 20000;
    commonVar->nNid = 50;
    commonVar->nRV = 0;

    /*********** Channel Estimator Initialisation **********/
    commonVar->nSbAlloc = 1;
    for(wnUInt8 sbIdx = 0; sbIdx < commonVar->nSbAlloc; sbIdx++)
    {
        commonVar->rbStart[sbIdx] = 15;
        commonVar->nPrb[sbIdx]    = 6;
    }

    commonVar->irc_enable = 0;
    commonVar->nullTone_enable = 1;
    commonVar->equType = 3;

    commonVar->FOEnable = 0;
    commonVar->FOFullAlloc = 1;
    commonVar->FOPRG = 2;

    //1 => enable.
    //0 => disable.
    commonVar->transformPrecode = 0;
    commonVar->ghFlag = 0;
    commonVar->shFlag = 0;
    commonVar->lmtBuffRm = 0;
    commonVar->maxNumLayForOneTB = 2;
    commonVar->maxModOrder = 6;

    commonVar->Rel15 = 1;
    commonVar->NIdRS = 0;
    commonVar->Lmax = 2;
    commonVar->offset = 1;
    commonVar->ZF = 0;
    /************************* END OF INITIALISATIONS *******************/

return;
}

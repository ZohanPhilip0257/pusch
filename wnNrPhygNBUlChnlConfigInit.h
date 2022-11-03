/**
 * @file wnNrPhygNBChnlConfig.h
 * @author Sakshama,Chandra Sekhar
 * @brief Phy gNB channel configuration file
 */

#ifndef WN_NR_PHY_GNB_UL_CHNL_CONFIG_INIT_H
#define WN_NR_PHY_GNB_UL_CHNL_CONFIG_INIT_H

#include "wnDataTypes.h"

#define MAX_SUMIMO_PUSCH_ANTENNA_PORTS              (8)
#define MAX_PUSCH_LAYERS                            (4)
#define MAX_NUM_ALLOC                              (275)
#define MAX_NUM_PRB                                (275)

#define MAX_UL_UE_PER_SLOT_SUPP                     (12)
#define MAX_SRS_UE                                  (12)

#define MAX_SUB_CARRIERS_OFDM_SYMBOL                (3300)
#define MAX_ZERO_PADDED_SUB_CARRIERS_OFDM_SYMBOL    (4096)
#define NUMBER_OF_OFDM_SYMBOLS_SLOT                 (14)
#define MAX_IFFT_CP_DATA_IN_SLOT                    (61440)
#define ALIGN_BOUNDARY                              (32)
#define MAX_SRS_SEQ_LEN                             (1632)

#define PI                  ( (float) 3.14159265358979 )
wnUInt8 check;
float snr;
double BER[10000][50];

// Common Uplink Configurations Structure
typedef struct commonUlConfigT
{
    /* Frame duplex type
       Value:
       0 = FDD
       1 = TDD*/
    wnUInt8 frameDuplexType;

    /** Cyclic prefix type
      * 0: normal
      * 1: extended */
    wnUInt8 cpType;
    
    /** Operating frequency band.
     * 0 : FR1,
     * 1 : FR2.*/
    wnUInt8 nrOperatingFreqBand;
    
    /* Subcarrier spacing
     * Values: 0 to 4
     */
    wnUInt8 numerology;
    
    /* Number of Physical antennas
    */
    wnUInt8 nPhysicalAnt;

    /* Slot Number
    */
    wnUInt8 nSlotNumber;

    /* Number of Cyclic Prefix Symbols
    */
    wnUInt16 cpSyms;

    /** IFFT size to be used
      */
    wnUInt16 fftSize;

    /** Number of PRBs in BW
    */
    wnUInt16 bwInPrbs;
    
    /** BW in MHz
      * Values: 5, 10, 15, 20, 25, 30, 40,50, 60, 70,
      * 80,90,100,200,400
      */
    wnUInt16 allocatedBw;
    
    /*Absolute frequency of point A in KHz
    Value: 450000 -> 52600000*/
    wnUInt16 carrierFreq;

    /*Inverse of fftSize
    */
    wnFlt inverseOfFftSize;

        /* Number of scheduled Pusch in current slot
     * Value 0 means no Pusch scheduled
     * Max value: MAX_PUSCH_PER_SLOT
    */
    wnUInt8 scheduledPusch;
    
    /* Number of scheduled PUCCH in current slot
     * Value 0 means no PUCCH scheduled
     * Max value: MAX_UE_PER_SLOT_SUPP
    */
    wnUInt8 scheduledPucchF0;
    wnUInt8 scheduledPucchF2;

    /* Number of scheduled SRS in current slot
     * Value 0 means no SRS scheduled
     * Max value: MAX_UE_PER_SLOT_SUPP
    */
    wnUInt8 scheduledSrs;
    
    /* Number of scheduled PRACH in current slot
     * Value 0 means no PRACH scheduled
     * Max value: MAX_UE_PER_SLOT_SUPP
    */
    wnUInt8 scheduledPrach;

    wnInt32 cellId;
    
}__attribute__((__packed__)) commonUlConfigT;


// PRACH Parameters
typedef struct prachConfigParams
{
  wnChar preambleFormat;

  /* 1 word */
  wnUInt8 Repetitions;
  wnUInt8 TotalPreambles;
  wnUInt8 Ncs;
  wnUInt8 FilterFactor;
  wnUInt8 rxAntennas;

  /* 1 Word */
  wnInt16 LRA;
  wnInt16 nFFT;
  wnInt16 pusch_FFTsize;
  
  /* 4 word */
  wnInt32 CpLength;
  wnInt32 GpLength;
  wnInt32 totalLength;
  wnInt32 rootSeqIndex;
  wnInt32 IFFTSize;
  wnInt32 SubcarrierStartIndex;
  wnInt32 prachScs;
  wnInt32 pusch_scs;
  wnInt32 kappaTcFs;
  wnUInt64 samplingFreq;
  wnUInt64 nSamplesSubframe;

}prachConfigParamsT;

// PUCCH Format 0 Parameters
typedef struct
{
    // Operating frequency band
    // 0 : FR1
    // 1 : FR2
    wnUInt8 nrOperatingFreqBand;

    // Allocated Bandwidth size in MHz
    wnUInt16 allocatedBandwidth;

    // Start of the resource grid for numerology mu
    // Higher-layer parameter 'offsetToCarrier' in 'SCS-SpecificCarrier' IE
    // Range: 0 to 2199
    wnUInt16  NStartMuGrid;
    wnUInt16  NStartMu0Grid;

    // Size of the resource grid for numerology mu
    // Higher-layer parameter 'carrierBandwidth' in 'SCS-SpecificCarrier' IE
    // Range: 1 to 275
    wnUInt16  NSizeMuGrid;
    wnUInt16  NSizeMu0Grid;

    // Physical-layer cell identity
    // Range: 0 to 1007
    // NCellID = 3*N(1)ID + N(2)ID
    // where N(1)ID has range 0 to 335 and N(2)ID has range 0 to 2
    wnUInt16  NCellID;

    // Higher-layer parameter 'hoppingId' in 'PUCCH-ConfigCommon' IE
    // Range: 0 to 1023
    wnUInt16  hoppingId;

    // Subcarrier spacing configuration (Numerology)
    // Higher-layer parameter 'subcarrierSpacing' in 'BWP' IE
    ////////////////////////////////////////////////////////////////
    ////////////  mu  Subcarrier Spacing = 2^mu*15 kHz  ////////////
    ////////////   0                  15                ////////////
    ////////////   1                  30                ////////////
    ////////////   2                  60                ////////////
    ////////////   3                 120                ////////////
    ////////////   4                 240                ////////////
    ////////////////////////////////////////////////////////////////
    wnUInt8   mu;

    // Largest mu value among the subcarrier spacing configurations
    // Higher-layer parameter 'scs-SpecificCarrierList' in 'FrequencyInfoUL' IE
    wnUInt8   mu0;

    // Number of sub-carriers in a resource block
    // There are 12 sub-carriers in a resource block
    wnUInt8   NRBsc;

    // Number of OFDM Symbol per slot
    // There are 14 symbols in a slot
    wnUInt8   NSlotSymb;

    // Slot number within a frame for numerology mu
    //////////////////////////////////////////////
    ////////////  mu  Range of nmusf  ////////////
    ////////////   0      0 to 9      ////////////
    ////////////   1      0 to 19     ////////////
    ////////////   2      0 to 39     ////////////
    ////////////   3      0 to 79     ////////////
    ////////////   4      0 to 159    ////////////
    //////////////////////////////////////////////
    wnUInt8   nmusf;

    // Number of UCI bits to be transmitted
    // Can be 1 or 2 for PUCCH Format 0
    wnUInt8   noOfBits;
    
    // SR Falg to check if it is being transmitted
    wnUInt8 SRFlag;

    // Index of the first PUCCH OFDM symbol in the slot (lDash)
    // Higher layer parameter 'startingSymbolIndex' in 'PUCCH-format0'
    wnUInt8   startSymbIndx;

    // Number of PUCCH symbols to be transmitted
    // Can be 1 or 2 for PUCCH Format 0
    // Obtained from higher layer parameter 'nrOfSymbols' in 'PUCCH-format0'
    wnUInt8   numbPucchSymb;

    // Number of PRBs
    // For PUCCH Format 0, it is always 1
    wnUInt8   numbPucchPrb;

    // Initial cyclic shift, mo
    // Higher-layer parameter 'initialCyclicShift' in 'PUCCH-format0'
    // Range: 0 to 11
    wnUInt8   mo;

    // Enables/disables group hop and/or sequence hop
    // Higher-layer parameter 'pucch-GroupHopping' in 'PUCCH-ConfigCommon' IE
    // grpHopFlag = {0, 1, 2}, where
    // 0 -> disable (group hopping disabled and sequence hopping enabled)
    // 1 -> enable  (group hopping enabled and sequence hopping disabled)
    // 2 -> neither (group hopping disabled and sequence hopping disabled)
    wnUInt8  grpHopFlag;

    // Enables/disables intra-slot frequency hopping
    // Higher-layer parameter 'intraSlotFrequencyHopping' in 'PUCCH-Resource'
    // 0 -> intra-slot frequency hopping disabled
    // 1 -> intra-slot frequency hopping enabled
    // For Format 0, it is always disabled i.e 0
    wnUInt8  intraSlotFreqHopEnable;

    // Antennas
    wnUInt8 numAntennas;

}__attribute__((__packed__)) pucchConfigF0T;

// PUCCH Format 2 Configurations Structure
typedef struct
{
    // Operating frequency band
    // 0 : FR1
    // 1 : FR2
    wnUInt8 nrOperatingFreqBand;

    // Allocated Bandwidth size in MHz
    wnUInt16 allocatedBandwidth;

    // Start of the resource grid for numerology mu
    // Higher-layer parameter 'offsetToCarrier' in 'SCS-SpecificCarrier' IE
    // Range: 0 to 2199
    wnUInt16  NStartMuGrid;
    wnUInt16  NStartMu0Grid;

    // Size of the resource grid for numerology mu
    // Higher-layer parameter 'carrierBandwidth' in 'SCS-SpecificCarrier' IE
    // Range: 1 to 275
    wnUInt16  NSizeMuGrid;
    wnUInt16  NSizeMu0Grid;

    // Physical-layer cell identity
    // Range: 0 to 1007
    // NCellID = 3*N(1)ID + N(2)ID
    // where N(1)ID has range 0 to 335 and N(2)ID has range 0 to 2
    wnUInt16  NCellID;

    // Higher-layer parameter 'dataScramblingIdentityPUSCH' in 'PUSCH-Config' IE
    // Range: 0 to 1023
    wnUInt16  dataScramblingIdentityPUSCH;

    // Radio Network Temporary Identifier (RNTI)
    // Given by the C-RNTI (Cell RNTI)
    // Range: 1 to 65519
    wnUInt16  nRNTI;

    // Given by higher-layer parameter 'scramblingID0' in 'DMRS-UplinkConfig' IE if provided
    // and by NCellID otherwise
    // If UE is configured with both 'dmrs-UplinkForPUSCH-MappingTypeA' and 'dmrs-UplinkForPUSCHMappingTypeB',
    // 'scramblingID0' is obtained from 'dmrs-UplinkForPUSCH-MappingTypeB'
    // Range: 0 to 65535
    wnUInt16 NID0;

    // Subcarrier spacing configuration (Numerology)
    // Higher-layer parameter 'subcarrierSpacing' in 'BWP' IE
    ////////////////////////////////////////////////////////////////
    ////////////  mu  Subcarrier Spacing = 2^mu*15 kHz  ////////////
    ////////////   0                  15                ////////////
    ////////////   1                  30                ////////////
    ////////////   2                  60                ////////////
    ////////////   3                 120                ////////////
    ////////////   4                 240                ////////////
    ////////////////////////////////////////////////////////////////
    wnUInt8   mu;

    // Largest mu value among the subcarrier spacing configurations
    // Higher-layer parameter 'scs-SpecificCarrierList' in 'FrequencyInfoUL' IE
    wnUInt8   mu0;

    // Number of sub-carriers in a resource block
    // There are 12 sub-carriers in a resource block
    wnUInt8   NRBsc;

    // Number of OFDM Symbol per slot
    // There are 14 symbols in a slot
    wnUInt8   NSlotSymb;

    // Slot number within a frame for numerology mu
    //////////////////////////////////////////////
    ////////////  mu  Range of nmusf  ////////////
    ////////////   0      0 to 9      ////////////
    ////////////   1      0 to 19     ////////////
    ////////////   2      0 to 39     ////////////
    ////////////   3      0 to 79     ////////////
    ////////////   4      0 to 159    ////////////
    //////////////////////////////////////////////
    wnUInt8   nmusf;

    // Number of UCI bits to be transmitted
    // Has to be greater than 2 for PUCCH Format 2
    // 3 to 11 bits use RM encoding
    // More than 11 bits use Polar encoding
    wnUInt16  noOfBits;

    // Index of the first PUCCH OFDM symbol in the slot (lDash)
    // Higher layer parameter 'startingSymbolIndex' in 'PUCCH-format2'
    wnUInt8   startSymbIndx;

    // Number of PUCCH symbols to be transmitted
    // Can be 1 or 2 for PUCCH Format 2
    // Obtained from higher layer parameter 'nrOfSymbols' in 'PUCCH-format2'
    wnUInt8   numbPucchSymb;

    // Number of PRBs
    // Range: 1 to 16
    // Obtained from higher layer parameter 'nrOfPRBs' in 'PUCCH-format2'
    wnUInt8   numbPucchPrb;

    wnUInt8 RbOffset;

    // Modulation order
    // Qm = 1 -> pi/2 BPSK
    // Qm = 2 -> QPSK
    // For Format 2, it is always 2 i.e. QPSK
    wnUInt8   Qm;

    // Enables/disables intra-slot frequency hopping
    // Higher-layer parameter 'intraSlotFrequencyHopping' in 'PUCCH-Resource'
    // 0 -> intra-slot frequency hopping disabled
    // 1 -> intra-slot frequency hopping enabled
    // For Format 2, it is always disabled i.e 0
    wnUInt8   intraSlotFreqHopEnable;

    // Antennas
    wnUInt8 numAntennas;

}__attribute__((__packed__)) pucchConfigF2T;

// PUSCH Configurations Structure
typedef struct puschConfigT
{
    /**Numerology.
     * Value: 0->4.
     */
    wnUInt8 numerology;

    /** Operating Frequency Band .*/
    wnUInt8 nrOperatingFreqBand;

    /** Allocated Band Width .*/
    wnUInt8 allocatedBw;

    /**USCH duration in symbols.
     * Value: 0->14.
     */
    wnUInt8 nNrOfSymbols;  

    /**Number of DM-RS CDM groups without data .
     *It determines the ratio of PDSCH EPRE to DM-RS EPRE.
     *Value: 1->3.
     */
    wnUInt8 nDmrsWithoutData;

    /**Number of Antenna Ports.
     * Value: 1->8.
     */
    wnUInt8 nAntennaPorts;

    /**Antenna port index to get DMRS ports, number_of_CDMGroups_withoutData and nNrOfDMRSSymbols
     * Value: 0->31.
     */
    wnUInt8 antPortValue;

    /**Antenna port index.
     *0: port 1000,
     *1: port 1001,
     *....
     *11: port 1011.
     *Value : 0->11.
     *Make sure the antenna ports are in the increasing order.*/
    wnUInt16 nPortIndex[MAX_SUMIMO_PUSCH_ANTENNA_PORTS];

    /* Indicates maximum number DMRS font-load symbols
        1 => single front loaded.
        2 => single front loaded or double front loaded.
    */
    wnUInt8 maxLen;

    /**UL DMRS symbol number.
     *Value: 1: single symbol   2: double symbol.
     */
    wnUInt8 nNrOfDMRSSymbols;

    /**UL DMRS config type.
     *0: type1,
     *1: type2.*/
    wnUInt8 nDMRSConfigType;

    /**UL DMRS Additional Positions.
     *Value: 0->3.*/
    wnUInt8 nDMRSAddPos;

   /**PDSCH mapping Type .
     * 0: mapping type A,
     * 1: mapping type B*/
    wnUInt8 nMappingType;

    /**higher-layer parameter DL-DMRS-typeA-pos.*/
    wnUInt8 ulDmrsTypeAPos; 

	 /** Starting OFDM symbol for the PUSCH.
	 *Value: 0->13 .*/
    wnUInt8 nStartSymbolIndex; 

    /**Slot Number.
     * Value: 0->319
     */
    wnUInt8 nSlotNumber;

    /**DMRS sequence initialization.
     *Should match what is sent in DCI 1_1,otherwise set to 0.
     *Value: 0->1.*/
    wnUInt8 nSCID;

    /**Number of Layers.
     * Value: 1->8.
     */      
    wnUInt8 nNrOfLayers;

    wnUInt16 lyrPortIndex[MAX_PUSCH_LAYERS];

    /** Modulation Order for Code Word.
     * Value: 1, 2, 4, 6, 8 .*/
    wnUInt8 modulationOrder;

    /** Redundancy version for Code Word.
     * Value:0, 1, 2, 3 .*/
    wnUInt8 nRV;

    /**bandwidth part start RB index from reference CRB.
     * Value: 0->274.
     */
    wnUInt16 nBWPStart;

    /**For resource allocation type 1.
     *The starting resource block within the BWP for this PUSCH.
	 *Value: 0->274
     */
    wnUInt16 nRBStart;

    /**For resource allocation type 1.
     *The Size resource block within the BWP for this PUSCH.
     *Value: 1->275
     */
    wnUInt16 nRBSize;

     /**Size of the BandwidthPart.
     * Value: 1->275.
     */
    wnUInt16 nBWPSize;

    /**The RNTI used for identifying the UE.
     *Value: 1 -> 65535.
     */
    wnUInt16 nRNTI;

    /**Data-scrambling-Identity.
     * Value : 0->1023.
     */
    wnUInt16 nNid;

    /**Physical Cell ID.
     * Value: 0->1007.
     */
    wnUInt16 nPhyCellId;

    /**DL-DMRS-Scrambling-ID 0.
     * If provided by the higher-layer and the PUSCH is scheduled by PDCCH with 
       CRC scrambled by C-RNTI or CS-RNTI, otherwise, L2 should set this to physical cell id.
     * Value: 0->65535.
     */
    wnUInt16 nNIDnSCID0;

    /**DL-DMRS-Scrambling-ID 1.
     * If provided by the higher-layer and the PUSCH is scheduled by PDCCH with 
       CRC scrambled by C-RNTI or CS-RNTI, otherwise, L2 should set this to physical cell id.
     * Value: 0->65535.
     */
    wnUInt16 nNIDnSCID1;

    /**PT-RS-to-PUSCH EPRE ratio.
     *DMRS EPRE Value :1->20000 : 0.001 dB step, -6dB to 14 dB
     */
    wnUInt16 nEpreRatioOfDmrsToSSB;

    /**TB Size for codeword.*/
    wnUInt32 nTbSize;

    /**Coderate for codeword.*/
    wnFlt codingRate;

    /** DMRS Symbol Positions .*/
    wnUInt8 dmrs1_sym;
    wnUInt8 dmrs2_sym;
    wnUInt8 dmrs3_sym;
    wnUInt8 dmrs4_sym;

    //0: Half-PRB, 1: 1 PRB, 2: 2 PRB
    wnUInt8 prMethod;

    //Data Positions or Sub Block Allocations and their information
    wnUInt16 nSbAlloc;

    wnUInt16 rbStart[MAX_NUM_ALLOC];

    wnUInt16 nPrb[MAX_NUM_ALLOC];

    wnUInt16 irc_enable;
    wnUInt16 nullTone_enable;
    wnUInt16 equType;
    wnUInt8 FOEnable;
    wnUInt8 FOFullAlloc;
    wnFlt FOPRG;
    wnUInt16 modOrder;

    /*Flag to indicate whether transform precoding is enabled or not.
    *0-> disabled
    *1-> enabled. */
    wnUInt8 transformPrecode;

    /*Flag to indicate whether group hoping is enabled or not. 
     * 0-> disabled
     * 1-> enabled*/
    wnUInt8 ghFlag;

    /*Flag to indicate whether sequence hoping is enabled or not. 
     * 0-> disabled
     * 1-> enabled*/
    wnUInt8 shFlag;

    /*0: Full buffer ratematching
     *1: Limited buffer ratematching*/
    wnUInt8 lmtBuffRm;

    wnUInt8 maxNumLayForOneTB;

    wnUInt8 maxModOrder;

    wnUInt8 Rel15;
    wnUInt32 NIdRS;
    wnFlt Lmax;
    wnUInt8 offset;
    wnUInt8 ZF;

}__attribute__((__packed__)) puschConfigT;

/*********** SRS Structures **************/ 
/* Hopping Type */
typedef enum
{
    neither      = 0,
    groupHopping = 1,
    seqHopping   = 2
}hoppingtype;

/* SRS Resource Type */
typedef enum
{
    aperiodic      = 0,
    semiPersistent = 1,
    Periodic       = 2,
}SRS_resource_Type;



typedef struct srsConfigT
{
  /* Number of UEs */
  wnUInt16          nUEs;    

  /* SRS-Sequence-Id; Value :: 0->1023  */
  wnUInt16          n_ID_SRS; 

  /* Transmission comb Value ; Value --> K_TC = {2,4} */
  wnUInt16          K_TC;   

  /* Transmission comb offset ; Value --> k_bar_TC = {0,1,...K_TC-1}                         */
  wnUInt16          k_bar_TC; 

  /* starting position; Value :: 0-->5  */
  wnUInt16          l_offset;  

  /* Number of SRS Symbols ; Value :: {1,2,4}  */
  wnUInt16          nNrofSRS_Symbols;  

  /* Repetition Factor; Value :: {1,2,4}  */
  wnUInt16          repetitionFactor;   

  /* Frequency Domain Position ; Value :: 0->67 */
  wnUInt16          nRRC;     

  /* Frequency Domain Shift; Value :: 0->268   */
  wnUInt16          n_shift; 

  /* SRS Bandwidth Configuration Index; C_SRS Value :: 0->63  */
  wnUInt16          Csrs;  

  /* SRS Bandwidth; Value :: 0->3  */
  wnUInt16          bsrs; 

  /* SRS Hopping Bandwidth; Value :: 0->3  */
  wnUInt16          bHop;  

  /* Cyclic Shifts for each UE  */
  wnUInt16          nCS_UE[MAX_SRS_UE]; 

  /* Amplitude Scaling factor */
  wnFlt             beta_srs;  

  /* Sub carrier Spacing :Numerology :: 0->4 */
  wnUInt16          numerology;   

  /* Value :0->2^NrUlSrsPDUs.nSubcSpacing*10 */
  wnUInt16          slot_number; 

  /* No. of PRB Allocated */
  wnUInt16          nPRBAllocated;   

  /* IFFT Size */
  wnUInt16          ifftSize; 

  /* 0->neither ; 1->Group Hopping Enable ; 2->Sequence Hopping Enable */
  hoppingtype       groupORSeqHopping;   

  /* 0->Aperiodic; 1->Semi-persistent ; 2-> Periodic  */
  SRS_resource_Type resourceType;   

  /* Number of SRS_Antenna Ports--> N_SRS_Ports  */
  wnUInt16          nNrOfAntennaPorts; 

  wnUInt16          Tsrs;

  wnUInt16          Toffset;

  /* Frame Number */
  wnUInt16          nf;                             

}__attribute__((__packed__)) srsConfigT;


// This is the final Ul structure
typedef struct gnbUlChnlStT
{
    commonUlConfigT commonUlConfig;

    pucchConfigF0T pucchConfigF0[MAX_UL_UE_PER_SLOT_SUPP];

    pucchConfigF2T pucchConfigF2[MAX_UL_UE_PER_SLOT_SUPP];

    puschConfigT   puschConfig[MAX_UL_UE_PER_SLOT_SUPP];

    srsConfigT     srsConfig;

    prachConfigParamsT prachConfig;

}__attribute__((__packed__)) gnbUlChnlStT;


gnbUlChnlStT gnbUlChnlSt;


wnVoid* wnNrPhygNBReDeMapMemPoolInit();
wnVoid* wnNrPhygNBPuschMemPoolInit();
wnVoid* wnNrPhygNBUlChnlConfigInit();
wnVoid* wnNrPhyFftCpRemoval(commonUlConfigT* commonUlConfig);
wnVoid wnNrPhyPuschMain(puschConfigT* pNrPuschInParams, commonUlConfigT* pNrUlCommonParams);
wnVoid* wnNrPhyUlBwFftCpConfig (commonUlConfigT* bwFftCp);
wnVoid wnNrPhycommonUlConfigInit();
wnVoid wnNrPhygNBCommonTables();
wnVoid wnNrPhyPuschInit(puschConfigT* commonVar);



#endif
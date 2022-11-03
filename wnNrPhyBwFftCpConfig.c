#include <stdio.h>
#include "wnNrPhyBwFftCpConfig.h"

void *wnNrPhyUlBwFftCpConfig (commonUlConfigT* bwFftCp)
{

    //printf("wnNrPhyBwFftCpConfig: Inside the FFT-CP config Init function\n");    

    if (bwFftCp->nrOperatingFreqBand == 0) // When the Frequency band is FR1
      {        
        if (bwFftCp->numerology == 0){            // 15 kHz Subcarrier spacing
            if(bwFftCp->allocatedBw == 5){
                bwFftCp->bwInPrbs = 25;
                bwFftCp->fftSize = 512;
                bwFftCp->cpSyms = 36;
                bwFftCp->inverseOfFftSize = 0.001953125;
            }else if(bwFftCp->allocatedBw == 10){
                bwFftCp->bwInPrbs = 52;
                bwFftCp->fftSize = 1024;
                bwFftCp->cpSyms = 72;
                bwFftCp->inverseOfFftSize = 0.0009765625;
            }else if(bwFftCp->allocatedBw == 15){
                bwFftCp->bwInPrbs = 79;
                bwFftCp->fftSize = 1536;
                bwFftCp->cpSyms = 108;
                bwFftCp->inverseOfFftSize = 0.00065104166;
            }else if(bwFftCp->allocatedBw == 20){
                bwFftCp->bwInPrbs = 106;
                bwFftCp->fftSize = 2048;
                bwFftCp->cpSyms = 144;
                bwFftCp->inverseOfFftSize = 0.00048828125;
            }else if(bwFftCp->allocatedBw == 25){
                bwFftCp->bwInPrbs = 133;
                bwFftCp->fftSize = 2048;
                bwFftCp->cpSyms = 144;
                bwFftCp->inverseOfFftSize = 0.00048828125;
            }else if(bwFftCp->allocatedBw == 30){
                bwFftCp->bwInPrbs = 160;
                bwFftCp->fftSize = 3072;
                bwFftCp->cpSyms = 216;
                bwFftCp->inverseOfFftSize = 0.000325520833;
            }else if(bwFftCp->allocatedBw == 40){
                bwFftCp->bwInPrbs = 216;
                bwFftCp->fftSize = 4096;
                bwFftCp->cpSyms = 288;
                bwFftCp->inverseOfFftSize = 0.000244140625;
            }else if(bwFftCp->allocatedBw == 50){
                bwFftCp->bwInPrbs = 270;
                bwFftCp->fftSize = 4096;
                bwFftCp->cpSyms = 288;
                bwFftCp->inverseOfFftSize = 0.000244140625;
            }else
                printf("Allocated bandwidth is not supported for this configuration");

        }else if(bwFftCp->numerology == 1){          // 30 kHz Subcarrier spacing
            if (bwFftCp->allocatedBw == 5){
                bwFftCp->bwInPrbs = 11;
                bwFftCp->fftSize = 256;
                bwFftCp->cpSyms = 18;
                bwFftCp->inverseOfFftSize = 0.00390625;
            }
            else if(bwFftCp->allocatedBw == 10){
                bwFftCp->bwInPrbs = 24;
                bwFftCp->fftSize = 512;
                bwFftCp->cpSyms = 36;
                bwFftCp->inverseOfFftSize = 0.001953125;
            }
            else if(bwFftCp->allocatedBw == 15){
                bwFftCp->bwInPrbs = 38;
                bwFftCp->fftSize = 768;
                bwFftCp->cpSyms = 54;
                bwFftCp->inverseOfFftSize = 0.0013020833;
            }
            else if(bwFftCp->allocatedBw == 20){
                bwFftCp->bwInPrbs = 51;
                bwFftCp->fftSize = 1024;
                bwFftCp->cpSyms = 72;
                bwFftCp->inverseOfFftSize = 0.0009765625;
            }else if(bwFftCp->allocatedBw == 25){
                bwFftCp->bwInPrbs = 65;
                bwFftCp->fftSize = 1024;
                bwFftCp->cpSyms = 72;
                bwFftCp->inverseOfFftSize = 0.0009765625;
            }else if(bwFftCp->allocatedBw == 30){
                bwFftCp->bwInPrbs = 78;
                bwFftCp->fftSize = 1536;
                bwFftCp->cpSyms = 108;
                bwFftCp->inverseOfFftSize = 0.00065104166;
            }else if (bwFftCp->allocatedBw == 40){
                bwFftCp->bwInPrbs = 106;
                bwFftCp->fftSize = 2048;
                bwFftCp->cpSyms = 144;
                bwFftCp->inverseOfFftSize = 0.00048828125;
            }else if(bwFftCp->allocatedBw == 50){
                bwFftCp->bwInPrbs = 133;
                bwFftCp->fftSize = 2048;
                bwFftCp->cpSyms = 144;
                bwFftCp->inverseOfFftSize = 0.00048828125;
            }else if (bwFftCp->allocatedBw == 60){
                bwFftCp->bwInPrbs = 162;
                bwFftCp->fftSize = 3072;
                bwFftCp->cpSyms = 216;
                bwFftCp->inverseOfFftSize = 0.000325520833;
            }else if (bwFftCp->allocatedBw == 80){
                bwFftCp->bwInPrbs = 217;
                bwFftCp->fftSize = 4096;
                bwFftCp->cpSyms = 288;
                bwFftCp->inverseOfFftSize = 0.000244140625;
            }else if(bwFftCp->allocatedBw == 90){
                bwFftCp->bwInPrbs = 245;
                bwFftCp->fftSize = 4096;
                bwFftCp->cpSyms = 288;
                bwFftCp->inverseOfFftSize = 0.000244140625;
            }else if (bwFftCp->allocatedBw == 100){
                bwFftCp->bwInPrbs = 273;
                bwFftCp->fftSize = 4096;
                bwFftCp->cpSyms = 288;
                bwFftCp->inverseOfFftSize = 0.000244140625;
            }
            else
                printf("Allocated bandwidth is not supported for this configuration");
        }else if (bwFftCp->numerology == 2 ) {        // 60 kHz Subcarrier spacing
            if (bwFftCp->allocatedBw == 10){
                bwFftCp->bwInPrbs = 11;
                bwFftCp->fftSize = 256;
                bwFftCp->cpSyms = 18;
                bwFftCp->inverseOfFftSize = 0.00390625;
            }else if(bwFftCp->allocatedBw == 15){
                bwFftCp->bwInPrbs = 18;
                bwFftCp->fftSize = 384;
                bwFftCp->cpSyms = 27;
                bwFftCp->inverseOfFftSize = 0.0026041666;
            }else if (bwFftCp->allocatedBw == 20){
                bwFftCp->bwInPrbs = 24;
                bwFftCp->fftSize = 512;
                bwFftCp->cpSyms = 36;
                bwFftCp->inverseOfFftSize = 0.001953125;
            }else if (bwFftCp->allocatedBw == 25){
                bwFftCp->bwInPrbs = 31;
                bwFftCp->fftSize = 512;
                bwFftCp->cpSyms = 36;
                bwFftCp->inverseOfFftSize = 0.001953125;
            }else if (bwFftCp->allocatedBw == 30){
                bwFftCp->bwInPrbs = 38;
                bwFftCp->fftSize = 768;
                bwFftCp->cpSyms = 54;
                bwFftCp->inverseOfFftSize = 0.0013020833;
            }else if (bwFftCp->allocatedBw == 40){
                bwFftCp->bwInPrbs = 51;
                bwFftCp->fftSize = 1024;
                bwFftCp->cpSyms = 72;
                bwFftCp->inverseOfFftSize = 0.0009765625;
            }else if (bwFftCp->allocatedBw == 50){
                bwFftCp->bwInPrbs = 65;
                bwFftCp->fftSize = 1024;
                bwFftCp->cpSyms = 72;
                bwFftCp->inverseOfFftSize = 0.0009765625;
            }else if (bwFftCp->allocatedBw == 60){
                bwFftCp->bwInPrbs = 79;
                bwFftCp->fftSize = 1536;
                bwFftCp->cpSyms = 108;
                bwFftCp->inverseOfFftSize = 0.00065104166;
            }else if (bwFftCp->allocatedBw == 80){
                bwFftCp->bwInPrbs = 107;
                bwFftCp->fftSize = 2048;
                bwFftCp->cpSyms = 144;
                bwFftCp->inverseOfFftSize = 0.00048828125;
            }else if (bwFftCp->allocatedBw == 90){
                bwFftCp->bwInPrbs = 121;
                bwFftCp->fftSize = 2048;
                bwFftCp->cpSyms = 144;
                bwFftCp->inverseOfFftSize = 0.00048828125;
            }else if(bwFftCp->allocatedBw == 100){
                bwFftCp->bwInPrbs = 135;
                bwFftCp->fftSize = 2048;
                bwFftCp->cpSyms = 144;
                bwFftCp->inverseOfFftSize = 0.00048828125;
            }else
                printf("Allocated bandwidth is not supported for this configuration");
        }

    }else if (bwFftCp->nrOperatingFreqBand == 1){     // When the Frequency band is FR2
        if(bwFftCp->numerology == 2 )  {          // 60 kHz Subcarrier spacing
            if (bwFftCp->allocatedBw == 50){
                bwFftCp->bwInPrbs = 66;
                bwFftCp->fftSize = 1024;
                bwFftCp->cpSyms = 72;
                bwFftCp->inverseOfFftSize = 0.0009765625;
            }else if (bwFftCp->allocatedBw == 100){
                bwFftCp->bwInPrbs = 132;
                bwFftCp->fftSize = 2048;
                bwFftCp->cpSyms = 144;
                bwFftCp->inverseOfFftSize = 0.00048828125;
            }else if (bwFftCp->allocatedBw == 200){
                bwFftCp->bwInPrbs = 264;
                bwFftCp->fftSize = 4096;
                bwFftCp->cpSyms = 288;
                bwFftCp->inverseOfFftSize = 0.000244140625;
            }else
                printf("Allocated bandwidth is not supported for this configuration");
        }else if( bwFftCp->numerology == 3){        // 120 kHz Subcarrier spacing
            if(bwFftCp->allocatedBw == 50){
                bwFftCp->bwInPrbs = 32;
                bwFftCp->fftSize = 512;
                bwFftCp->cpSyms = 36;
                bwFftCp->inverseOfFftSize = 0.001953125;
            }else if (bwFftCp->allocatedBw == 100){
                bwFftCp->bwInPrbs = 66;
                bwFftCp->fftSize = 1024;
                bwFftCp->cpSyms = 72;
                bwFftCp->inverseOfFftSize = 0.0009765625;
            }else if (bwFftCp->allocatedBw == 200){
                bwFftCp->bwInPrbs = 132;
                bwFftCp->fftSize = 2048;
                bwFftCp->cpSyms = 144;
                bwFftCp->inverseOfFftSize = 0.00048828125;
            }else if(bwFftCp->allocatedBw == 400){
                bwFftCp->bwInPrbs = 264;
                bwFftCp->fftSize = 4096;
                bwFftCp->cpSyms = 288;
                bwFftCp->inverseOfFftSize = 0.000244140625;
            }
            else
                printf("Allocated bandwidth is not supported for this configuration");
        }
    }
    zerosPaddedForIfft = (bwFftCp->fftSize - bwFftCp->bwInPrbs*12)>>1;

    return NULL;
}

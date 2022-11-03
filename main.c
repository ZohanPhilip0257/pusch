#include "wnNrPhygNBUlChnlConfigInit.h"
#include<stdio.h>

int main()
{
	// Memory Init
	wnNrPhygNBReDeMapMemPoolInit();

	// pool of memory elements few assigned with 0's and constants (DM-RS related tables).
	wnNrPhygNBPuschMemPoolInit();

	// Configure the UL PHY channels
	wnNrPhygNBUlChnlConfigInit();

	// Calls Rx IFFT and CP Removal Function
	wnNrPhyFftCpRemoval(&gnbUlChnlSt.commonUlConfig);

	// PUSCH Chain
	wnNrPhyPuschMain(&gnbUlChnlSt.puschConfig[0], &gnbUlChnlSt.commonUlConfig);
	
	printf("cBLER: %d\n", check);
	
	return 0;
}
#include "Coherence5LEWrapper.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*           Program to test wrapper
 *
 * to compile :
 *   gcc -o testwrapper testwrapper.c -L. -lCoherence5LEWrapper
 *   export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH
 */

int main(int argc, char *argv[]) {

    TCoh3 FileInfo;

    Eeg3_DebugFileSwitch(0);

// Initialize
    if (Eeg3_Initialisation()!=0)
        fprintf(stderr,"Eeg3_Initialisation() failed\n");	
   

// Unlock
    TUnlock3LE unlockStruct;
    unlockStruct.int1=1; // insert your Unlock Code here
	unlockStruct.int2=2;
	unlockStruct.int3=3;
	unlockStruct.int4=4;

    if (Eeg3_Unlock(unlockStruct)!=0){
    	fprintf(stdout,"Eeg3_Unlock failed\n");
        exit(0);
    }


// Version
    TVersion version;
    Eeg3_Version(&version);
    fprintf(stdout,"Eeg3_Version: %d.%d.%d.%d\n", version.major, version.minor, version.compile, version.number);

// Open File
    if (Eeg3_OpenFile("test_0001.eeg",&FileInfo)!=0){
        fprintf(stderr,"Eeg3_OpenFile failed\n");
        exit(1);
    }
    else{
       fprintf(stdout,"Eeg3_OpenFile: File opened successfully\n\tDuration: %d seconds\n\tFrequency: %d Hz\n\tElectrodes: %d\n",FileInfo.duration,FileInfo.frequency,FileInfo.electrodes);
    }

// Test Eeg3_GetEeg
    int res; 
    int start = 0;
    int blocksize = 20000;
    short* HBuf1 = malloc(blocksize*sizeof(short)*FileInfo.electrodes);

    res = 0; // return code of Eeg3_GetEeg
    while (res>=0){ // EOF
        res = Eeg3_GetEeg(start,blocksize,HBuf1);
        if (res>=0)        
            start +=res; // increase start for next block
    }
    fprintf(stdout,"Eeg3_GetEeg: %d samples read (end with code: %i)\n",start,res);
    res = Eeg3_GetEeg(0,blocksize,HBuf1);
    fprintf(stdout,"\tfirst 5 values: %d %d %d %d %d\n",HBuf1[0],HBuf1[1],HBuf1[2],HBuf1[3],HBuf1[4]);
    free(HBuf1);
    int lastsample=start;
    
// GetMarkersNumber
    int markerNumber = Eeg3_GetMarkersNumber(1,lastsample);  
    fprintf(stdout,"Eeg3_GetMarkersNumber: returned %d\n",markerNumber);

 

// Get Markers
    struct TMarker* markerBuf = malloc(markerNumber*sizeof(struct TMarker));
    if ( (res=Eeg3_GetMarkers(1,lastsample,markerBuf))<0)
        fprintf(stderr,"Eeg3_GetMarkers failed with code: %d\n",res);
    else
        fprintf(stdout,"Eeg3_GetMarkers: %d markers found\n\tcomment of first marker at position %d : %s\n",res,markerBuf[0].pos,markerBuf[0].text);

// Put Marker
   struct TMarker newMarker;
    newMarker.pos = 1000;
    newMarker.duration = 1;
    newMarker.evttype= 2;
    strcpy(newMarker.text,"Test Marker");

    if ( (res=Eeg3_PutMarker(&newMarker))<0)
        fprintf(stderr,"Eeg3_PutMarker failed with code: %d\n",res);
    else
        fprintf(stdout,"Eeg3_PutMarker: success with code: %d\n",res);
         
    

// Get Impedances 
    struct TImpedances imp;
    res = Eeg3_GetImpedances(lastsample,&imp); // start at last sample
    if (res<0)
        fprintf(stderr,"Eeg3_GetImpedances failed with Code %d\n",res);
    else
        fprintf(stdout,"Eeg3_GetImpedances successfull with code %d\n\tposition %d\n\ttext: %s\n",res, imp.pos,imp.text);


// Next File
    char* nextfile = malloc(260);
    if ( (res=Eeg3_NextFile(1,nextfile,&FileInfo))!=0)
        fprintf(stderr,"Eeg3_NextFile failed with code %d %s\n",res,nextfile);
    else{
       fprintf(stdout,"Eeg3_NextFile: File %s opened successfully\n\tDuration: %d seconds\n\tFrequency: %d Hz\n\tElectrodes: %d\n",nextfile, FileInfo.duration,FileInfo.frequency,FileInfo.electrodes);
    }

    if ( (res=Eeg3_NextFile(1,nextfile,&FileInfo))!=0)
        fprintf(stderr,"Eeg3_NextFile failed with code %d %s\n",res,nextfile);
    else{
       fprintf(stdout,"Eeg3_NextFile: File %s opened successfully\n\tDuration: %d seconds\n\tFrequency: %d Hz\n\tElectrodes: %d\n",nextfile, FileInfo.duration,FileInfo.frequency,FileInfo.electrodes);
    }


// Close File
    if (Eeg3_CloseFile()!=0)
        fprintf(stderr,"Eeg3_CloseFile failed\n");
    else   
        fprintf(stdout,"Eeg3_CloseFile(): closed\n");


// Terminate
    Eeg3_Termination();    


}

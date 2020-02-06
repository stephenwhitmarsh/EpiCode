/* Coherence5LEWrapper.c
 * Wrapper for Coherence5LE
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unistd.h>
#include "Coherence5LEstructs.h"

FILE *tochild, *fromchild;
int debug;
int initialized;
TCoh3 FileInfo;

/*
 * for testing purposes ... rename into main
 
int test_main(int argc, char *argv[]) {


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
        fprintf(stderr,"Eeg3_NextFile failed with code %d\n",res);
    else{
       fprintf(stdout,"Eeg3_NextFile: File opened successfully\n\tDuration: %d seconds\n\tFrequency: %d Hz\n\tElectrodes: %d\n",FileInfo.duration,FileInfo.frequency,FileInfo.electrodes);
    }

// Close File
    if (Eeg3_CloseFile()!=0)
        fprintf(stderr,"Eeg3_CloseFile failed\n");
    else   
        fprintf(stdout,"Eeg3_CloseFile(): closed\n");


// Terminate
    Eeg3_Termination();    
    

}  */

/*
 *  Shared Object Functions
 */

/*

void _init()
{
	//printf("Inside _init()\n"); 
} 
void _fini()
{
	//printf("Inside _fini()\n"); 
}
*/

/*
 *  Wrapper Functions
 */


int Eeg3_DebugFileSwitch(int val){
    debug=val;

    return 1;
}

int	Eeg3_OpenFile(char* filename, TCoh3* tch){
    int res=-166; // if something goes wrong
    char command[5] = "opfi";
    fwrite(command,1,5,tochild);
    fprintf(tochild,"%s\n",filename);

    fread(&res,1,sizeof(int),fromchild);
    if (res==0){
        fread(tch,1,sizeof(TCoh3),fromchild);
    }
    FileInfo = *tch;
   

    return res;
}

int	Eeg3_Version(TVersion* version){
    int res=-166; // if something goes wrong
    char command[5] = "gver";
    fwrite(command,1,5,tochild);
    
    fread(&res,1,sizeof(int),fromchild);
    if (res==0){
        fread(version,1,sizeof(TVersion),fromchild);
    }

    return res;
}

int	Eeg3_GetImpedances(int startpos, TImpedances* imps){
    int res=-166; // if something goes wrong
    char command[5] = "gimp";
    fwrite(command,1,5,tochild);

    fwrite(&startpos,1,sizeof(int),tochild);

    fread(&res,1,sizeof(int),fromchild);

    if (res>=0){
        fread(imps,1,sizeof(TImpedances),fromchild);
    }

    return res;
}

int	Eeg3_PutMarker(TMarker* marker){
    int res=-166; // if something goes wrong
    char command[5] = "pmar";
    fwrite(command,1,5,tochild);
    fwrite(marker,1,sizeof(TMarker),tochild);
    fread(&res,1,sizeof(int),fromchild);
    return res;
}

int Eeg3_NextFile(int direction, char* filename, TCoh3* FileInfo){
    int res=-166; // if something goes wrong
    char command[5] = "nefi";
    fwrite(command,1,5,tochild);
    fwrite(&direction,1,sizeof(int),tochild);

    // read filename

    char* tfilename = malloc(260);
    fgets(tfilename,260,fromchild);
    //strncpy(filename,tfilename,strlen(tfilename)-1);  
    // read 
    fread(&res,1,sizeof(int),fromchild);
    if (res==0){
        fread(FileInfo,1,sizeof(struct TCoh3),fromchild);
    }
    free(tfilename);
    return res;

}

int Eeg3_GetEeg(int start, int duration, short* buffer){
    int res=-166; // if something goes wrong
    char command[5] = "geeg";
    fwrite(command,1,5,tochild);

    int msg[2];
    msg[0]=start;
    msg[1]=duration;
    fwrite(msg,1,2*sizeof(int),tochild);      
    fread(&res,1,sizeof(int),fromchild);
    int count =-1;    
    if (res>0)
        count=fread(buffer,1,res*sizeof(short)*FileInfo.electrodes,fromchild);   
    return res;
}

int Eeg3_GetMarkers(int start, int end, TMarker* buffer){
    int res=-166; // if something goes wrong
    char command[5] = "gmar";
    fwrite(command,1,5,tochild);

    int msg[2];
    msg[0]=start;
    msg[1]=end;
    fwrite(msg,1,2*sizeof(int),tochild);      
    fread(&res,1,sizeof(int),fromchild);
    if (res>0)
        fread(buffer,1,res*sizeof(TMarker),fromchild);   
    return res;
}


int Eeg3_CloseFile(void){
    int res=-166; // if something goes wrong
    char command[5] = "clfi";
    fwrite(command,1,5,tochild);

    fread(&res,1,sizeof(int),fromchild);
    return res;
}


int Eeg3_GetMarkersNumber(int begin, int end){
    int res=-166; // if something goes wrong
    char command[5] = "gman";
    fwrite(command,1,5,tochild);
    

    int msg[2];
    msg[0]=begin;
    msg[1]=end;
    fwrite(msg,1,2*sizeof(int),tochild);    

    fread(&res,1,sizeof(int),fromchild);
    return res;
}

int Eeg3_Initialisation(){

    if (initialized)
        return -1; // allready initialized
    pid_t pid;

    int p2c[2],c2p[2];
   
    pipe(p2c);
    pipe(c2p);
     
    pid = fork();

    if (pid == 0) { /* child */

        close(p2c[1]);

        // connect pipe output to stdin
        dup2( p2c[0], 0 );

        close(c2p[0]);

        dup2( c2p[1], 1 );
          

        // loading image of cat application
        int res = execl("./Coherence5LE", "Coherence5LE", "wrapper", NULL);
        /*system("./Coherence5LE wrapper &");
        int res =1;*/
        fprintf(stderr,"execl failed with code %d\n",res);

    }  else if (pid < 0){           // failed to fork
        fprintf(stderr,"Fork failed\n");
        exit(1);
    }else {

        tochild = fdopen(p2c[1], "w");
        fromchild = fdopen(c2p[0], "r");
        
        // disable buffering, write data to buffer after each operation
        setvbuf(tochild,(char*)NULL,_IONBF,0);
        setvbuf(fromchild,(char*)NULL,_IONBF,0);

        // set DebugFileSwitch first ... if necessary
        int res=-166; // if something goes wrong
        char command[5] = "debu";
        fwrite(command,1,5,tochild);

        fwrite(&debug,1,sizeof(int),tochild);

        fread(&res,1,sizeof(int),fromchild);

        res=-166; // if something goes wrong
        char nextcommand[5] = "init";
        fwrite(nextcommand,1,5,tochild);
        fread(&res,1,sizeof(int),fromchild);
        if (res==0)
           	initialized = 1;
        return res;
    }
}

int Eeg3_Unlock(TUnlock3LE unlockstruct){
    int res=-166; // if something goes wrong
    char command[5] = "unlk";
    fwrite(command,1,5,tochild);
    fwrite(&unlockstruct,1,sizeof(unlockstruct),tochild);

    fread(&res,1,sizeof(int),fromchild);
    return res;
}

int Eeg3_Termination(){
    int res=-166; // if something goes wrong
    char command[5] = "term";
    fwrite(command,1,5,tochild);
    
    fread(&res,1,sizeof(int),fromchild);
  	initialized = 0;
    char nextcommand[5] = "exit";
    fwrite(command,1,5,tochild);
    return res;
}




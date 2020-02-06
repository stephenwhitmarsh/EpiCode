#include <stdio.h>
#include <windows.h>
#include "Coherence5LE.h"


HANDLE h; 
int usedebug;
TCoh3 FileInfo;

int main(int argc, char *argv[]) {


    if (argc!=2 || strcmp(argv[1],"wrapper")){
        fprintf(stdout,"program can only be used by wrapper %d %s\n",argc, argv[1]);
        exit(0);
    }
    setvbuf(stdout,(char*)NULL,_IONBF,0);

    char command[5];
    
    int exit =0;    

    while (!exit && !feof(stdin)) {
        
        int bytesread = fread(&command, 1, 5,stdin);
           
        if (strcmp(command, "init") == 0){
            int res = dll_Eeg3_Initialisation();
            fwrite(&res,1,sizeof(int),stdout);
        }else if (strcmp(command, "exit")==0){
            exit = 1;
        }else if (strcmp(command, "unlk") == 0){
            internal_Eeg3_Unlock();
        }else if (strcmp(command, "opfi") == 0){
            internal_Eeg3_OpenFile();
        }else if (strcmp(command, "term")==0){
            int res = dll_Eeg3_Termination();
            fwrite(&res,1,sizeof(int),stdout);
        }else if (strcmp(command, "geeg")==0){
            internal_Eeg3_GetEeg();
        }else if (strcmp(command, "gman")==0){
            internal_Eeg3_GetMarkersNumber();
        }else if (strcmp(command, "clfi")==0){
            internal_Eeg3_CloseFile();
        }else if (strcmp(command, "gver")==0){
            internal_Eeg3_Version();
        }else if (strcmp(command, "gmar")==0){
            internal_Eeg3_GetMarkers();
        }else if (strcmp(command, "pmar")==0){
            internal_Eeg3_PutMarker();
        }else if (strcmp(command, "gimp")==0){
            internal_Eeg3_GetImpedances();
        }else if (strcmp(command, "nefi")==0){
            internal_Eeg3_NextFile();
        }else if (strcmp(command, "debu")==0){
            internal_Eeg3_DebugFileSwitch();
        }
    }
}

int internal_Eeg3_Unlock(){
    TUnlock3LE unlockstruct;
    fread(&unlockstruct,1,sizeof( TUnlock3LE),stdin);    
    int res = dll_Eeg3_Unlock(unlockstruct);
    fwrite(&res,1,sizeof(int),stdout);
    return res;
}

int internal_Eeg3_CloseFile(){
    int res = dll_Eeg3_CloseFile();
    fwrite(&res,1,sizeof(int),stdout);
    return res;
}


int internal_Eeg3_DebugFileSwitch(){
    int val= -1;

    fread(&val,1,sizeof(int),stdin);
    int res = dll_Eeg3_DebugFileSwitch(val);
    fwrite(&res,1,sizeof(int),stdout);
    return res;
}

int internal_Eeg3_GetEeg(){
    int ints[2];
    fread(&ints,1,2*sizeof(int),stdin);    
    int buffersize = ints[1]*sizeof(short)*FileInfo.electrodes;

    short* HBuf = malloc(buffersize);
    int res = -1;
    res = dll_Eeg3_GetEeg(ints[0],ints[1],HBuf);
    fwrite(&res,1,sizeof(int),stdout); // write return value


    if (res>0)
        fwrite(HBuf,1,res*sizeof(short)*FileInfo.electrodes,stdout);


    free(HBuf);
    return res;
}

int internal_Eeg3_GetMarkers(){
    int ints[2];
    fread(&ints,1,2*sizeof(int),stdin);    

    int markerNumber = dll_Eeg3_GetMarkersNumber(ints[0],ints[1]);
    
    TMarker* markerBuf = malloc(markerNumber*sizeof(TMarker));
    int res = -1;
    res = dll_Eeg3_GetMarkers2(ints[0],ints[1],markerBuf);
    fwrite(&res,1,sizeof(int),stdout); // write return value

    if (res>0)
        fwrite(markerBuf,1,markerNumber*sizeof(TMarker),stdout);
    free(markerBuf);
    return res;
}

int internal_Eeg3_GetImpedances(){    
    TImpedances imps;
    int res = -1;
    int startpos = -1;
    fread(&startpos,1,sizeof(int),stdin);    // get startpos

    res = dll_Eeg3_GetImpedances(startpos,&imps);
    fwrite(&res,1,sizeof(int),stdout); // write return value
    if (res>=0)
        fwrite(&imps,1,sizeof(TImpedances),stdout);

    return res;
}

int internal_Eeg3_PutMarker(){
    TMarker marker;
    fread(&marker,1,sizeof(TMarker),stdin);    

    int res = -1;
    res = dll_Eeg3_PutMarker(&marker);
    fwrite(&res,1,sizeof(int),stdout); // write return value

    return res;
}

int internal_Eeg3_GetMarkersNumber(){
    int ints[2];
    fread(&ints,1,2*sizeof(int),stdin);    

    int res = -1;
    res = dll_Eeg3_GetMarkersNumber(ints[0],ints[1]);
    fwrite(&res,1,sizeof(int),stdout); // write return value

    return res;
}


int internal_Eeg3_OpenFile(){
    // read in filename
    char filename[260]="";
    fgets(filename,260,stdin);
    char cleanfilename[260]="";  // we have to remove \n at end
    strncpy(cleanfilename,filename,strlen(filename)-1);    
    char fullname[260];
    GetFullPathName(cleanfilename,260,fullname,NULL);
       

    int res = dll_Eeg3_OpenFile(fullname,&FileInfo);
    
    fwrite(&res,1,sizeof(int),stdout);
    if (res==0){ // OpenFile succeeded, transfer TCoh3
        fwrite(&FileInfo,1,sizeof(struct TCoh3),stdout);
    }

    return res;
}

int internal_Eeg3_NextFile(){
    // read direction
    int direction = 0;
    fread(&direction,1,sizeof(int),stdin);


    char filename[260]=""; 

    //NextFile  
    int res = dll_Eeg3_NextFile(direction,filename,&FileInfo);
    fprintf(stdout,"%s\n",filename);
    fwrite(&res,1,sizeof(int),stdout);
    if (res==0){ // OpenFile succeeded, transfer TCoh3
        fwrite(&FileInfo,1,sizeof(struct TCoh3),stdout);
    }

    return res;
}

int internal_Eeg3_Version(){
    TVersion version;
    int res = dll_Eeg3_Version(&version);
    fwrite(&res,1,sizeof(int),stdout);
    if (res==0){ // OpenFile succeeded, transfer TCoh3
        fwrite(&version,1,sizeof(version),stdout);
    }

    return res;
}


int dll_Eeg3_DebugFileSwitch(int val){
    usedebug = val;
    return 1;
}

int dll_Eeg3_Initialisation(){ 

    int (/*WINAPI?*/ *local_Eeg3_Initialisation)(void);     
    int (/*WINAPI?*/ *local_Eeg3_DebugFileSwitch)(int);
    h = LoadLibrary("Coherence5LE.dll"); 
	if (h==NULL)
        fprintf(stderr,"Loading DLL failed");

    dll_Eeg3_Version = (void *) GetProcAddress(h, (LPCSTR)3);
    dll_Eeg3_OpenFile = (void *) GetProcAddress(h, (LPCSTR)4);
    dll_Eeg3_CloseFile = (void *) GetProcAddress(h, (LPCSTR)5);
    dll_Eeg3_GetEeg = (void *) GetProcAddress(h, (LPCSTR)6);
    dll_Eeg3_PutMarker = (void *) GetProcAddress(h, (LPCSTR)7);
    dll_Eeg3_GetMarkers = (void *) GetProcAddress(h, (LPCSTR)8);
    dll_Eeg3_GetImpedances = (void *) GetProcAddress(h, (LPCSTR)9);
	dll_Eeg3_Unlock = (void *) GetProcAddress(h, (LPCSTR)10);
    dll_Eeg3_NextFile = (void *) GetProcAddress(h, (LPCSTR)12);
    dll_Eeg3_GetEeg2 =(void *) GetProcAddress(h, (LPCSTR)13);
    dll_Eeg3_GetMarkers2 = (void *) GetProcAddress(h, (LPCSTR)14);
    dll_Eeg3_GetMarkersNumber = (void *) GetProcAddress(h, (LPCSTR)15);

    local_Eeg3_Initialisation = (void *) GetProcAddress(h, (LPCSTR)1);
    local_Eeg3_DebugFileSwitch = (void *) GetProcAddress(h, (LPCSTR)11);     

    // Debug
    int res;
    if ((res=local_Eeg3_DebugFileSwitch(usedebug))<0)
        fprintf(stderr,"Eeg3_DebugFileSwitch failed to create debug file with code: %d\n",res);   
    res = local_Eeg3_Initialisation();
    return res;
}

int dll_Eeg3_Termination(){ 
    int (/*WINAPI?*/ *local_Eeg3_Termination)(void);
    local_Eeg3_Termination = (void *)  GetProcAddress(h, (LPCSTR)2);
    int res = local_Eeg3_Termination();
    //FreeLibrary(h);
    return res;
}



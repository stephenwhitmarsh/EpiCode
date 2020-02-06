#include <stdlib.h>
#pragma pack(push,1)
#include <windows.h>

/* Only define EXTERN_C if it hasn't been defined already. This allows
 * individual modules to have more control over managing their exports.
 */
#ifndef EXTERN_C

#ifdef __cplusplus
  #define EXTERN_C extern "C"
#else
  #define EXTERN_C extern
#endif

#endif


#include "Coherence5LEstructs.h"

EXTERN_C int Eeg3_NextFile(int,char*,TCoh3*); // 12
EXTERN_C int Eeg3_Initialisation(void);  //1
EXTERN_C int Eeg3_Termination(void);    // 2
EXTERN_C int Eeg3_Unlock(TUnlock3LE);    // 10
EXTERN_C int Eeg3_OpenFile(char*,TCoh3*); // 4  
EXTERN_C int Eeg3_CloseFile(void); // 5
EXTERN_C int Eeg3_GetEeg(int,int,HGLOBAL*); // 6
EXTERN_C int Eeg3_GetMarkers(int,int,HGLOBAL*);	    // 8
EXTERN_C int Eeg3_GetImpedances(int,TImpedances*);  // 9
EXTERN_C int Eeg3_Version(TVersion*);               // 3
EXTERN_C int Eeg3_PutMarker(TMarker*);              // 7
EXTERN_C int Eeg3_GetEeg2(int,int,short*); // 6
EXTERN_C int Eeg3_GetMarkersNumber(int,int);
EXTERN_C int Eeg3_GetMarkers2(int,int,int*); // hack, because Matlab doesn't support struct arrays
EXTERN_C int Eeg3_DebugFileSwitch(int);


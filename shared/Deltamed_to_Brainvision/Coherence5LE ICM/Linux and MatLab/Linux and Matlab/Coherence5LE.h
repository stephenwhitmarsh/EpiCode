#include <windows.h>
#include "Coherence5LEstructs.h"
// FUNCTIONS


int (/*WINAPI?*/ *dll_Eeg3_Version)(TVersion*);
int	(/*WINAPI?*/ *dll_Eeg3_OpenFile) (char*, struct TCoh3*);
int (/*WINAPI?*/ *dll_Eeg3_CloseFile)(void);
int (/*WINAPI?*/ *dll_Eeg3_Unlock)(TUnlock3LE);
int (/*WINAPI?*/ *dll_Eeg3_GetEeg)(int, int, short*);
int (/*WINAPI?*/ *dll_Eeg3_PutMarker)(TMarker*);
int (/*WINAPI?*/ *dll_Eeg3_GetMarkers)(int, int, HGLOBAL*);
int (/*WINAPI?*/ *dll_Eeg3_GetImpedances)(int,  TImpedances*);
int (/*WINAPI?*/ *dll_Eeg3_NextFile)(short, char*, TCoh3*);
int (/*WINAPI?*/ *dll_Eeg3_GetEeg2)(int, int, short*);
int (/*WINAPI?*/ *dll_Eeg3_GetMarkers2)(int, int, TMarker*);
int (/*WINAPI?*/ *dll_Eeg3_GetMarkersNumber)(int, int);





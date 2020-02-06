/* Coherence5LEWrapper.h */
#include "Coherence5LEstructs.h"

int Eeg3_Initialisation();
int Eeg3_Termination();
int Eeg3_Unlock(TUnlock3LE);
int Eeg3_GetEeg(int, int, short*);
int	Eeg3_OpenFile(char*, TCoh3*);
int Eeg3_CloseFile(void);
int Eeg3_Version(TVersion*);
int Eeg3_GetMarkers(int, int, int*); // hack, because Matlab doesn't support struct arrays
int	Eeg3_PutMarker(TMarker*);
int	Eeg3_GetImpedances(int, TImpedances*);
int Eeg3_NextFile(int, char*, TCoh3*);
int Eeg3_DebugFileSwitch(int);
int Eeg3_GetMarkersNumber(int,int);

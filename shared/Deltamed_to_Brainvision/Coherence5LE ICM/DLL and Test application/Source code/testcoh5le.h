//#include <stdlib.h>
#define MAXELEC 1024
#pragma pack(push,1)

typedef struct						// objet événement de Coherence5LE
{
	int
		typeevt,					// type de l'événement
		pos,						  // position début
		duree;						  // durée de l'événement
	char
		buf[80];
} TEventLe;

typedef struct						// objet événement de Coherence5LE
{
	int
		typeevt,					// type de l'événement
		pos,						  // position début
		duree;						  // durée de l'événement
	char
		buf[252];
} TEventLeLong;

typedef struct
{
	int
		duree,						  // duree de l'enregistrement (en secondes)
		frequence,					// fréquence d'echantillonage de chaque groupe
		basetemps,					// base de tps de la fréq d'échan de chaque groupe (sec)
		totalvoies;					// nombre total de voies
	char							    // (entre parenthèses le nb d'octets)
		date[20],					  // date et heure de début d'enregistrement
                        // (hh:mm:ss jj/mm/aaa)
		nom[MAXELEC][8],		// nom de la voie (8)
		type[MAXELEC];			// type de voie (1)
	int
		theta[MAXELEC],			// coordonnée angulaire theta (4, -200 à +200 grades)
		phi[MAXELEC],				// coordonnée angulaire phi (4, -100 à +100 grades)
		r[MAXELEC],					// rayon (4)
		minanal[MAXELEC],		// valeur min analogique (4 int stocke=1000*float ori)
		maxanal[MAXELEC],		// valeur max analogique (4 int stocke=1000*float ori)
		minconv[MAXELEC],		// valeur min du convertisseur (4)
		maxconv[MAXELEC];		// valeur max du convertisseur (4)
	char
		unite[MAXELEC][4];	// unite (µV,mmHg...) (4 char)
} TCoh3Le;

typedef struct
{
	unsigned char
		texte[20];				  // texte de l'événement impédance
	int
		pos,					      // position de l'événement impédance
		nbvoies,				    // nombre de voies utilisées
		imped[MAXELEC];			// valeurs des impédances (en Kohm, max=250)
} TImpedances3;

typedef struct				  // Version de la DLL
{
	short
		major,				      // version majeure
		minor,				      // sous-version (sur 100)
		compile,				    // compilation (jour)
		number;				      // numero dans la journee
} Tversion;

// YA 26/09/2007 : structure permettant de débloquer la DLL si les bonnes
// valeurs sont entrées
typedef struct
{
  int int1;
  int int2;
  int int3;
  int int4;
} TUnlock3LE;

//WR 28/01/09
typedef struct
{
  char
    name[50],        // patient name
    firstname[30],   // firstname
    date[11],        // date of birth
    sex,             // sex (M/F)
    file[20],        // file number
    center[39],      // origin center of the recording
    commentary[256]; // commentary
}TPatientInfo3LE;

typedef struct
{
	int
		pos,
                realtime;
}TRealTimeMarker;

typedef struct
{
  int Shift;
  short Drift;
}TVideoShiftAndDrift;

typedef int (*TEeg3_Initialisation)(void);              // 1
typedef int (*TEeg3_Termination)(void);                 // 2
typedef int (*TEeg3_Version)(Tversion*);                // 3
typedef int (*TEeg3_OpenFile)(char*,TCoh3Le*);          // 4
typedef int (*TEeg3_CloseFile)(void);                   // 5
typedef int (*TEeg3_GetEeg)(int,int,HGLOBAL);           // 6
typedef int (*TEeg3_PutMarker)(TEventLe*);              // 7
typedef int (*TEeg3_GetMarkers)(int,int,HGLOBAL*);      // 8
typedef int (*TEeg3_GetImpedances)(int,TImpedances3*);  // 9
typedef int (*TEeg3_Unlock)(TUnlock3LE);                // 10
typedef int (*TEeg3_DebugFileSwitch)(bool);             // 11
typedef int (*TEeg3_NextFile)(int,char*,TCoh3Le*);      // 12
typedef int (*TEeg3_GetEeg2)(int,int,short*);           // 13
typedef int (*TEeg3_GetMarkers2)(int,int,TEventLe*);    // 14
typedef int (*TEeg3_GetMarkersNumber)(int,int);         // 15
typedef int (*TEeg3_GetPatientInfo)(TPatientInfo3LE *infopat); //16
typedef int (*TEeg3_GetNumberOfBlocFiles)(char*);       // 17
typedef int (*TEeg3_GetAvailableBlocFiles)(char*, String*); //18
typedef int (*TEeg3_GetNextMarker)(int,TEventLe*);      //20
typedef int (*TEeg3_TempFolderSwitch)(bool);            //21
typedef int (*TEeg3_GetMarkersLong)(int,int,HGLOBAL*);  //22
typedef int (*TEeg3_GetMarkersLong2)(int,int,TEventLeLong*); //23
typedef int (*TEeg3_GetNextMarkerLong)(int,TEventLeLong*); // 24
typedef int (*TEeg3_GetRealTimeMarkers)(HGLOBAL*); // 25  //WR 04/11/2011
typedef int (*TEeg3_GetNumberOfRealTimeMarkers)(void); // 26  //WR 07/11/2011
typedef int (*TEeg3_GetRealTimeMarkers2)(TRealTimeMarker *realTimeMrk); // 27  //WR 07/11/2011
typedef int (*TEeg3_GetVideoShiftAndDrift)(TVideoShiftAndDrift *shiftAndDrift); //28 //WR 22/12/11


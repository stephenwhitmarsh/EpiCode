// STRUCTS

typedef struct TCoh3
{
	int duration;	// duration of the record in seconds
	int frequency;	// sampling rate
	int timebase;	// sampling time base in seconds
	int electrodes;	// number of electrodes
	char date[20];	// date and time of the record (‘hh:mm:ss dd/mm/yyyy’)
	char name[8192];	// electrode names
	char type[1024];	// electrode types
    int theta[1024];	// theta angular coordinate
	int phi[1024]; // phi angular coordinate
	int r[1024];	// radius
	int minanal[1024];	// min analog value
	int maxanal[1024];	// max analog value
	int minconv[1024];	// min value of the converter
	int maxconv[1024];	// max value of the converter
	char unit[4096];	// display unit (‘µV’, ‘mmHg’…)
} TCoh3;


typedef struct TUnlock3LE
{
  int int1,		// first integer		
     int2,		// second integer	
     int3, 		// third integer		
     int4; 		// fourth integer	
} TUnlock3LE;

typedef struct TVersion
{
    short
		major,	// major
		minor,	// minor version
		compile,	// compilation
		number;	// number
} TVersion;


typedef struct TMarker
{
	int evttype;	// type of the marker
    int pos;	// position
	int duration;	// duration
	char text[80];	// comment
} TMarker;


typedef struct TImpedances
{
	char text[20];	// text associated to the impedance test
	int pos;	// position of the impedance test
	int nbchannel;	// number of channels used in the impedance test
	int imped[128];	// impedances values (in Kohm, range : 0 to 250)
} TImpedances;

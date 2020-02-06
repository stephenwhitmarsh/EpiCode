//---------------------------------------------------------------------------
#ifndef mainletestH
#define mainletestH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include"testcoh5le.h"
#include <Dialogs.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <ExtCtrls.hpp>
#include <time.h>
#include <ComCtrls.hpp>
#include <Buttons.hpp>

extern HANDLE hlib;

extern TEeg3_Version Eeg3_Version;
extern TEeg3_OpenFile Eeg3_OpenFile;
extern TEeg3_CloseFile Eeg3_CloseFile;
extern TEeg3_GetEeg Eeg3_GetEeg;
extern TEeg3_GetMarkers Eeg3_GetMarkers;
extern TEeg3_GetImpedances Eeg3_GetImpedances;
extern TEeg3_PutMarker Eeg3_PutMarker;
extern TEeg3_Unlock Eeg3_Unlock;
extern TEeg3_DebugFileSwitch Eeg3_DebugFileSwitch;
extern TEeg3_NextFile Eeg3_NextFile;
extern TEeg3_GetEeg2 Eeg3_GetEeg2;
extern TEeg3_GetMarkers2 Eeg3_GetMarkers2;
extern TEeg3_GetMarkersNumber Eeg3_GetMarkersNumber;
extern TEeg3_GetPatientInfo Eeg3_GetPatientInfo;
extern TEeg3_GetNumberOfBlocFiles Eeg3_GetNumberOfBlocFiles;
extern TEeg3_GetAvailableBlocFiles Eeg3_GetAvailableBlocFiles;
extern TEeg3_Initialisation Eeg3_Initialisation; //WR 26/02/09 : pour test init et termination dll
extern TEeg3_Termination Eeg3_Termination; //WR 26/02/09 : pour test init et termination dll
extern TEeg3_GetNextMarker Eeg3_GetNextMarker;
extern TEeg3_GetMarkersLong Eeg3_GetMarkersLong; //WR 04/11/11
extern TEeg3_GetMarkersLong2 Eeg3_GetMarkersLong2; //WR 04/11/11
extern TEeg3_GetNextMarkerLong Eeg3_GetNextMarkerLong; //WR 04/11/11
extern TEeg3_GetRealTimeMarkers Eeg3_GetRealTimeMarkers; //WR 04/11/11
extern TEeg3_GetNumberOfRealTimeMarkers Eeg3_GetNumberOfRealTimeMarkers; // 26  //WR 07/11/2011
extern TEeg3_GetRealTimeMarkers2 Eeg3_GetRealTimeMarkers2; // 27  //WR 07/11/2011
extern TEeg3_GetVideoShiftAndDrift Eeg3_GetVideoShiftAndDrift; //28 //WR 22/12/11

#define INC 64
#define clPerso (TColor)0xf4f8f7
//---------------------------------------------------------------------------
class TForm3le : public TForm
{
__published:	// Composants gérés par l'EDI
    TOpenDialog *OpD;
   TButton *ButtonMk;
   TButton *ButtonEv;
   TImage *Image;
   TButton *ButtonGo;
   TTimer *Timer;
   TEdit *Edit;
   TLabel *LabelT;
   TLabel *Label1;
   TLabel *Label2;
   TEdit *EditChan;
   TUpDown *UpDChan;
  TStatusBar *StatusBar;
  TSpeedButton *SpeedButtonOnline;
  TLabel *NbSecs;
  TTimer *Timer2;
  TLabel *NbPts;
        TButton *btPrev;
        TButton *btNext;
        TCheckBox *cbGetEeg2;
  TButton *btGetEvent;
  TLabel *Label3;
  TEdit *edDeb;
  TLabel *Label4;
  TEdit *edFin;
  TButton *btGetEvent2;
  TLabel *Label5;
  TEdit *edDeb2;
  TLabel *Label6;
  TEdit *edFin2;
  TButton *btMax;
  TButton *btMax2;
  TScrollBox *ScrollBox1;
  TMemo *Memo;
        TButton *TermBt;
        TCheckBox *UseMarkersLongCb;
        TCheckBox *UseMarkers2Cb;
        TButton *GetRealTimeBt;
        TButton *GetSAndDBt;
    void __fastcall ButtonMkClick(TObject *Sender);
    void __fastcall ButtonEvClick(TObject *Sender);
   void __fastcall ButtonGoClick(TObject *Sender);
   void __fastcall FormDestroy(TObject *Sender);
   void __fastcall TimerTimer(TObject *Sender);
   void __fastcall EditChange(TObject *Sender);
   void __fastcall FormCloseQuery(TObject *Sender, bool &CanClose);
   void __fastcall FormActivate(TObject *Sender);
   void __fastcall ImageMouseDown(TObject *Sender, TMouseButton Button,
          TShiftState Shift, int X, int Y);
   void __fastcall EditClick(TObject *Sender);
   void __fastcall UpDChanClick(TObject *Sender, TUDBtnType Button);
   void __fastcall ImageMouseMove(TObject *Sender, TShiftState Shift,
          int X, int Y);
  void __fastcall FormCreate(TObject *Sender);
  void __fastcall SpeedButtonOnlineClick(TObject *Sender);
  void __fastcall Timer2Timer(TObject *Sender);
        void __fastcall btPrevClick(TObject *Sender);
        void __fastcall btNextClick(TObject *Sender);
  void __fastcall btGetEventClick(TObject *Sender);
  void __fastcall btGetEvent2Click(TObject *Sender);
  void __fastcall btMaxClick(TObject *Sender);
  void __fastcall btMax2Click(TObject *Sender);
        void __fastcall FirstBtClick(TObject *Sender);
        void __fastcall lastBtClick(TObject *Sender);
        void __fastcall TermBtClick(TObject *Sender);
        void __fastcall Button1Click(TObject *Sender);
        void __fastcall GetRealTimeBtClick(TObject *Sender);
        void __fastcall GetSAndDBtClick(TObject *Sender);
private:	// Déclarations de l'utilisateur
 TCoh3Le Dsc;
// TEventLe Evt;
 Tversion Ver;
 TPatientInfo3LE infoPat;
 Graphics::TBitmap *Bitmap;
 int y0[8][512],
     s0,
     Start,
     Chanmax,
     Chan0,
     Nbmark,
     Mark,
     Duree;
 char Sbuf[256];
 bool _run;
 HGLOBAL Heeg;
 HGLOBAL Hevt;
 short *eeg;
 TEventLe *evt;
 TEventLeLong *evtlg;
 int DataRead(int);
 int DataRead2(int); //WR 18/11/2008
 int ScanMark(void);
 int ScanMark2(void); //WR 20/11/2008
public:		// Déclarations de l'utilisateur
    __fastcall TForm3le(TComponent* Owner);
    void __fastcall updateAffichage(); //WR 17/11/2008
    int __fastcall nextFile(int direction); //WR 18/11/2008
 int NbPointsRead; // YA
 char* fileName;
};

//---------------------------------------------------------------------------
extern PACKAGE TForm3le *Form3le;
//---------------------------------------------------------------------------
#endif
 
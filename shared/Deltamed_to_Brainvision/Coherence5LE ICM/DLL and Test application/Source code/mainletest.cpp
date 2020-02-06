//---------------------------------------------------------------------------
#include <vcl.h>
#include <stdio.h>
#include <values.h>
#pragma hdrstop
USERES("testcoh3le.res");
#include "mainletest.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"

//YA 25/09/2013 : unlock codes for ICM
#define UNLOCK1	654125
#define UNLOCK2	894265
#define UNLOCK3	512743
#define UNLOCK4	554190

HANDLE hlib;

TEeg3_Initialisation Eeg3_Initialisation;
TEeg3_Termination Eeg3_Termination;
TEeg3_Version Eeg3_Version;
TEeg3_OpenFile Eeg3_OpenFile;
TEeg3_CloseFile Eeg3_CloseFile;
TEeg3_GetEeg Eeg3_GetEeg;
TEeg3_PutMarker Eeg3_PutMarker;
TEeg3_GetMarkers Eeg3_GetMarkers;
TEeg3_GetImpedances Eeg3_GetImpedances;
TEeg3_Unlock Eeg3_Unlock; // YA 26/09/2007
TEeg3_DebugFileSwitch Eeg3_DebugFileSwitch; //WR 26/11/2008
TEeg3_NextFile Eeg3_NextFile; //WR 17/11/2008
TEeg3_GetEeg2 Eeg3_GetEeg2; //WR 18/11/2008
TEeg3_GetMarkers2 Eeg3_GetMarkers2; //WR 20/11/2008
TEeg3_GetMarkersNumber Eeg3_GetMarkersNumber; //WR 20/11/2008
TEeg3_GetPatientInfo Eeg3_GetPatientInfo; //WR 29/01/09
TEeg3_GetNumberOfBlocFiles Eeg3_GetNumberOfBlocFiles; //WR 29/01/09
TEeg3_GetAvailableBlocFiles Eeg3_GetAvailableBlocFiles; //Wr 29/01/09
TEeg3_GetNextMarker Eeg3_GetNextMarker; //WR 31/05/10
TEeg3_TempFolderSwitch Eeg3_TempFolderSwitch; //WR 07/06/10
TEeg3_GetMarkersLong Eeg3_GetMarkersLong; //WR 04/11/11
TEeg3_GetMarkersLong2 Eeg3_GetMarkersLong2; //WR 04/11/11
TEeg3_GetNextMarkerLong Eeg3_GetNextMarkerLong; //WR 04/11/2011
TEeg3_GetRealTimeMarkers Eeg3_GetRealTimeMarkers; //WR 04/11/2011
TEeg3_GetNumberOfRealTimeMarkers Eeg3_GetNumberOfRealTimeMarkers; // 26  //WR 07/11/2011
TEeg3_GetRealTimeMarkers2 Eeg3_GetRealTimeMarkers2; // 27  //WR 07/11/2011
TEeg3_GetVideoShiftAndDrift Eeg3_GetVideoShiftAndDrift; //28 //WR 22/12/11

HANDLE DLLinit(void)
{
 HANDLE hlib;

 hlib=LoadLibrary("Coherence5le.dll");
 if(!hlib)
 {
  return NULL;
 }
 Eeg3_Initialisation=(TEeg3_Initialisation)GetProcAddress(hlib,(LPCSTR) 1);
 Eeg3_Termination=   (TEeg3_Termination)   GetProcAddress(hlib,(LPCSTR) 2);
 Eeg3_Version=       (TEeg3_Version)       GetProcAddress(hlib,(LPCSTR) 3);
 Eeg3_OpenFile=      (TEeg3_OpenFile)      GetProcAddress(hlib,(LPCSTR) 4);
 Eeg3_CloseFile=     (TEeg3_CloseFile)     GetProcAddress(hlib,(LPCSTR) 5);
 Eeg3_GetEeg=        (TEeg3_GetEeg)        GetProcAddress(hlib,(LPCSTR) 6);
 Eeg3_PutMarker=     (TEeg3_PutMarker)     GetProcAddress(hlib,(LPCSTR) 7);
 Eeg3_GetMarkers=    (TEeg3_GetMarkers)    GetProcAddress(hlib,(LPCSTR) 8);
 Eeg3_GetImpedances= (TEeg3_GetImpedances) GetProcAddress(hlib,(LPCSTR) 9);
 Eeg3_Unlock=        (TEeg3_Unlock)        GetProcAddress(hlib,(LPCSTR) 10); // YA 26/09/2007
 Eeg3_DebugFileSwitch= (TEeg3_DebugFileSwitch) GetProcAddress(hlib,(LPCSTR) 11); // WR 26/11/2008
 Eeg3_NextFile=      (TEeg3_NextFile)      GetProcAddress(hlib,(LPCSTR) 12); // WR 17/11/2008
 Eeg3_GetEeg2=       (TEeg3_GetEeg2)       GetProcAddress(hlib,(LPCSTR) 13); // WR 18/11/2008
 Eeg3_GetMarkers2=   (TEeg3_GetMarkers2)   GetProcAddress(hlib,(LPCSTR) 14); // WR 20/11/2008
 Eeg3_GetMarkersNumber= (TEeg3_GetMarkersNumber) GetProcAddress(hlib,(LPCSTR) 15); // WR 20/11/2008
 Eeg3_GetPatientInfo=(TEeg3_GetPatientInfo)GetProcAddress(hlib,(LPCSTR) 16); //WR 29/01/09
 Eeg3_GetNumberOfBlocFiles= (TEeg3_GetNumberOfBlocFiles) GetProcAddress(hlib,(LPCSTR) 17); //WR 29/01/09
 Eeg3_GetAvailableBlocFiles= (TEeg3_GetAvailableBlocFiles) GetProcAddress(hlib,(LPCSTR) 18); //WR 29/01/09
 Eeg3_GetNextMarker= (TEeg3_GetNextMarker) GetProcAddress(hlib,(LPCSTR) 20); //WR 31/05/10
 Eeg3_TempFolderSwitch=(TEeg3_TempFolderSwitch) GetProcAddress(hlib,(LPCSTR) 21); //WR 07/06/10
 Eeg3_GetMarkersLong=(TEeg3_GetMarkersLong)GetProcAddress(hlib,(LPCSTR) 22); //WR 04/11/11
 Eeg3_GetMarkersLong2=(TEeg3_GetMarkersLong2)GetProcAddress(hlib,(LPCSTR) 23); // WR 04/11/11
 Eeg3_GetNextMarkerLong= (TEeg3_GetNextMarkerLong) GetProcAddress(hlib,(LPCSTR) 24); //WR 04/11/11
 Eeg3_GetRealTimeMarkers= (TEeg3_GetRealTimeMarkers) GetProcAddress(hlib,(LPCSTR) 25); //WR 04/11/11
 Eeg3_GetNumberOfRealTimeMarkers= (TEeg3_GetNumberOfRealTimeMarkers) GetProcAddress(hlib,(LPCSTR) 26); //WR 07/11/2011
 Eeg3_GetRealTimeMarkers2= (TEeg3_GetRealTimeMarkers2) GetProcAddress(hlib,(LPCSTR) 27); //WR 07/11/2011
 Eeg3_GetVideoShiftAndDrift=(TEeg3_GetVideoShiftAndDrift) GetProcAddress(hlib,(LPCSTR) 28); //WR 22/12/11
 return hlib;
}

TForm3le *Form3le;
extern int result;
//---------------------------------------------------------------------------

int TForm3le::ScanMark(void)
{
  GlobalUnlock(Hevt);

  /*if(UseMarkersLongCb->Checked)
  {
    Nbmark=Eeg3_GetMarkersLong(0,Dsc.duree*Dsc.frequence,&Hevt);
    Nbmark=GlobalSize(Hevt)/sizeof(TEventLeLong);
    evtlg=(TEventLeLong*)GlobalLock(Hevt);
    for(Mark=0;Mark<Nbmark;Mark++)
      if(evtlg[Mark].pos*INC/Dsc.frequence>Start)
        break;
  }
  else
  {*/
    Nbmark=Eeg3_GetMarkers(0,Dsc.duree*Dsc.frequence,&Hevt);
    Nbmark=GlobalSize(Hevt)/sizeof(TEventLe);
    evt=(TEventLe*)GlobalLock(Hevt);
    for(Mark=0;Mark<Nbmark;Mark++)
      if(evt[Mark].pos*INC/Dsc.frequence>Start)
        break;
  //}

 return Nbmark;
}

//---------------------------------------------------------------------------

int TForm3le::ScanMark2(void)
{

 Nbmark=Eeg3_GetMarkersNumber(0,Dsc.duree*Dsc.frequence);

 if(Nbmark>0){
   /*if(UseMarkersLongCb->Checked)
   {
     evtlg = new TEventLeLong[Nbmark];

     Eeg3_GetMarkersLong2(0,Dsc.duree*Dsc.frequence,evtlg);

     for(Mark=0;Mark<Nbmark;Mark++)
      if(evtlg[Mark].pos*INC/Dsc.frequence>Start)
       break;
     }
   else
   {*/
     evt = new TEventLe[Nbmark];

     Eeg3_GetMarkers2(0,Dsc.duree*Dsc.frequence,evt);

     for(Mark=0;Mark<Nbmark;Mark++)
      if(evt[Mark].pos*INC/Dsc.frequence>Start)
       break;
   //}
 }

 return Nbmark;
}

//---------------------------------------------------------------------------

int TForm3le::DataRead(int start)
{
 int chan,
     smp,
     ret,
     x0=start&511;

 ret=Eeg3_GetEeg(Dsc.frequence/INC*start,Dsc.frequence*2,Heeg); //two seconds readed

 eeg=(short*)GlobalLock(Heeg);
 for(chan=0;chan<8;chan++)
  for(smp=0;smp<128;smp++)
   y0[chan][x0+smp]=ret<0?0:
                    eeg[Dsc.frequence/INC*smp*Dsc.totalvoies+Chan0+chan]/20;
 GlobalUnlock(Heeg);

 //if (ret<=0) ShowMessage("oups");

 return ret;
}

//---------------------------------------------------------------------------

int TForm3le::DataRead2(int start)
{
 int chan,
     smp,
     ret,
     x0=start&511;

 ret=Eeg3_GetEeg2(Dsc.frequence/INC*start,Dsc.frequence*2,eeg); //two seconds readed

 for(chan=0;chan<8;chan++){
   for(smp=0;smp<128;smp++){
     y0[chan][x0+smp]=ret<0?0:eeg[Dsc.frequence/INC*smp*Dsc.totalvoies+Chan0+chan]/20;
   }
 }

 //if (ret<=0) ShowMessage("oups");

 return ret;
}


//---------------------------------------------------------------------------
__fastcall TForm3le::TForm3le(TComponent* Owner)
    : TForm(Owner)
{

  hlib=DLLinit();
  if(!hlib)                        // DLL introuvable
  {
    Application->MessageBox("Can't load coherence5le.dll","Error",
                MB_SYSTEMMODAL|MB_OK|MB_ICONHAND);
    Application->Terminate();
    return;
  }
  else {
    //Eeg3_DebugFileSwitch(true); //WR 26/11/2008
    Eeg3_TempFolderSwitch(false); //WR 07/06/10
    result=Eeg3_Initialisation();
  }
  if(result<0)
  {
    char buf[80];

    sprintf(buf,"Error %d",result);
    Application->MessageBox("Can't use coherence5le.dll",buf,
                MB_SYSTEMMODAL|MB_OK|MB_ICONHAND);
    Application->Terminate();
    return;
  }


 int j;
 TUnlock3LE UL3LE;

 Eeg3_Version(&Ver);
 Memo->Lines->Add(AnsiString("Coherence5LE.dll : ")+Ver.major+"."+Ver.minor+"."+
                                       Ver.compile+"."+Ver.number+"\r\n");
 if(!OpD->Execute())
 {
  Application->Terminate();
  return;
 }

 fileName = OpD->FileName.c_str();

 UL3LE.int1=UNLOCK1;
 UL3LE.int2=UNLOCK2;
 UL3LE.int3=UNLOCK3;
 UL3LE.int4=UNLOCK4;

 if (Eeg3_Unlock(UL3LE)<0)
 {
  Application->MessageBox("TestCoh5LE.exe has failed to unlock Coherence5LE.dll. Execution will stop.","Information",
               MB_SYSTEMMODAL|MB_OK|MB_ICONHAND);
  Application->Terminate();
  return;
 }

 j=Eeg3_OpenFile(fileName,&Dsc);
 if (j==-1)
 {
  Application->MessageBox("Coherence5LE.dll is locked, you must unlock it.",fileName,
               MB_SYSTEMMODAL|MB_OK|MB_ICONHAND);
  Application->Terminate();
  return;
 }
 if (j<0)
 {
  Application->MessageBox(fileName,"Can't open",
               MB_SYSTEMMODAL|MB_OK|MB_ICONHAND);
  Application->Terminate();
  return;
 }

 ScanMark();
 Bitmap=new Graphics::TBitmap;
 Bitmap->Width=512;
 Bitmap->Height=512;
 Image->Canvas->Draw(0,0,Bitmap);

 updateAffichage();
 edFin->Text = IntToStr(Dsc.duree);
 edFin2->Text = IntToStr(Dsc.duree*Dsc.frequence);
 Heeg=GlobalAlloc(GMEM_MOVEABLE|GMEM_ZEROINIT,Dsc.totalvoies*Dsc.frequence*sizeof(short)*2); //two seconds readed
 eeg=new short[Dsc.totalvoies*Dsc.frequence*2]; //two seconds readed
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::FormCloseQuery(TObject *Sender, bool &CanClose)
{
 if(Timer->Enabled)
  CanClose=false;
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::FormDestroy(TObject *Sender)
{
 delete Bitmap;
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::ButtonMkClick(TObject *Sender)
{
 TEventLe evt={0,0,0,"Marker"};

 evt.typeevt=1;
 evt.pos=Start*Dsc.frequence/INC;
 Mark=Eeg3_PutMarker(&evt);
 Start-=256;
 Start=Start/128*128;
 ScanMark2();
 if(Start<0)
  Start=0;
 Edit->Text=AnsiString(Start/INC)+" s";
 sprintf(Sbuf,"%02d:%02d:%02d",Start/INC/3600,Start/INC/60%60,Start/INC%60);
 LabelT->Caption=Sbuf;
 _run=false;
 ActiveControl=Edit;
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::ButtonEvClick(TObject *Sender)
{
 TEventLe evt={0,0,0,"Event"};

 evt.typeevt=0;
 evt.pos=Start*Dsc.frequence/INC;
 Mark=Eeg3_PutMarker(&evt);
 Start-=256;
 Start=Start/128*128;
 ScanMark2();
 if(Start<0)
  Start=0;
 Edit->Text=AnsiString(Start/INC)+" s";
 sprintf(Sbuf,"%02d:%02d:%02d",Start/INC/3600,Start/INC/60%60,Start/INC%60);
 LabelT->Caption=Sbuf;
 _run=false;
 ActiveControl=Edit;
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::ButtonGoClick(TObject *Sender)
{
 Eeg3_CloseFile();

 Eeg3_OpenFile(fileName,&Dsc);
 if(Duree!=Dsc.duree)
 {
  Memo->Lines->Add(AnsiString(Duree=Dsc.duree)+" seconds    "+Dsc.frequence+" Hz");
  Start-=512;
  if(Start<0)
   Start=0;
  Start=Start/128*128;
  ScanMark2();
  _run=false;
 }
 if(!Timer->Enabled)
 {
  ClientHeight=850;
  Timer->Enabled=true;

  s0=clock()-Start*1000/INC;

  if(!_run)
  {
   Bitmap->Canvas->Brush->Color=clPerso;
   Bitmap->Canvas->FillRect(Rect(0,0,512,512));
   Image->Canvas->Draw(0,0,Bitmap);
   s0-=8000;
  }
  ButtonGo->Caption="Stop";
  Edit->Enabled=false;
  ButtonMk->Enabled=false;
  ButtonEv->Enabled=false;
 }
 else
 {
  Timer->Enabled=false;
  ButtonGo->Caption="Start";
  Edit->Enabled=true;
  ButtonMk->Enabled=true;
  ButtonEv->Enabled=true;
  ActiveControl=Edit;
 }
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::TimerTimer(TObject *Sender)
{
 static AnsiString S;
 static int dw=-1,      // précédent
            w;
 bool _read=false;
 int j,
     chan,
     x0,      // point courant
     x1;

 for(j=0;Start<(clock()-s0)*INC/1000&&Timer->Enabled;j++)
 {
  if((Start&127)==0)
   if(!_run)
   {
     if(cbGetEeg2->Checked){
       if((j=DataRead2(Start))<0);
     } else {
       if((j=DataRead(Start))<0);
     }
   }
   else
    _read=true;
  x0=Start%512;
  //StatusBar->SimpleText="a";//IntToStr(x0);
  x1=(Start+511)%512;
  if(Start%INC==0)
  {
   int t=(Start/INC)%60;

   sprintf(Sbuf,"%02d:%02d:%02d",Start/INC/3600,Start/INC/60%60,Start/INC%60);
   if(t)
    S=t;
   else
    S=Sbuf;
   dw=(w=Bitmap->Canvas->TextWidth(S)+4)+1;
   Bitmap->Canvas->Pen->Color=(TColor)0x00e0ff;
  }
  else
   Bitmap->Canvas->Pen->Color=clPerso;
  if(Mark<Nbmark&&evt[Mark].pos*INC/Dsc.frequence<=Start)
  {
   AnsiString E=AnsiString(evt[Mark].buf+evt[Mark].typeevt*4);

   Bitmap->Canvas->Pen->Color=clRed;
   Bitmap->Canvas->Font->Color=clRed;
   Bitmap->Canvas->TextOut(x0-Bitmap->Canvas->TextWidth(E)-4,500,E);
   Bitmap->Canvas->TextOut(x0+512-Bitmap->Canvas->TextWidth(E)-4,500,E);
   Mark++;
  }
  if(!dw)
  {
   LabelT->Caption=Sbuf;
   Edit->Text=AnsiString((Start)/INC)+" s";
   Bitmap->Canvas->Font->Color=clBlack;
   Bitmap->Canvas->TextOut(x0-w,Memo->Top/2,S);
   Bitmap->Canvas->TextOut(x0+512-w,Memo->Top/2,S);
  }
  dw--;
  Bitmap->Canvas->MoveTo(x0,0);
  Bitmap->Canvas->LineTo(x0,512);
  Bitmap->Canvas->Pen->Color=clBlack;
  for(chan=Chan0;chan<Chanmax;chan++)
  {
   int o=256/(Chanmax-Chan0)+512/(Chanmax-Chan0)*(chan-Chan0);

   Bitmap->Canvas->MoveTo(x1,o+y0[chan&7][x1]);
   if(x0==0)
   {
    Bitmap->Canvas->LineTo(512,o+y0[chan&7][x0]);
    Bitmap->Canvas->MoveTo(-1,o+y0[chan&7][x1]);
   }
   Bitmap->Canvas->LineTo(x0,o+y0[chan&7][x0]);
  }
  Start++;
  if(Start>Duree*INC)
  {
   _read=0;
   break;
  }
 }
 if(!_run)
  _read=true;
 _run=true;

 Image->Canvas->Draw(-x0-1,0,Bitmap);
 Image->Canvas->Draw(-x0-1+512,0,Bitmap);

 if (Start>Duree*INC)
   ButtonGo->Click();
 if(_read)
 {
   if(cbGetEeg2->Checked){
     if((j=DataRead2(((Start+128)>>7)<<7))<0)
     _read=0;
   } else {
     if((j=DataRead(((Start+128)>>7)<<7))<0)
     _read=0;
   }
 }
  StatusBar->SimpleText=LabelT->Caption;//IntToStr(x0);
}


//---------------------------------------------------------------------------
void __fastcall TForm3le::EditClick(TObject *Sender)
{
 Edit->SelectAll();
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::EditChange(TObject *Sender)
{
 if(!Timer->Enabled)
 {
  Start=atof(Edit->Text.c_str())*INC;
  if(Start>(Duree-1)*INC)
  {
   Start=(Duree-1)*INC;
   Edit->Text=AnsiString(Duree-1)+" s";
  }
  sprintf(Sbuf,"%02d:%02d:%02d",Start/INC/3600,Start/INC/60%60,Start/INC%60);
  LabelT->Caption=Sbuf;
  for(Mark=0;Mark<Nbmark;Mark++)
   if(evt[Mark].pos*INC/Dsc.frequence>Start)
    break;
  _run=false;
 }
}
//---------------------------------------------------------------------------


void __fastcall TForm3le::FormActivate(TObject *Sender)
{

 Application->BringToFront();
 ActiveControl=Edit;
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::ImageMouseMove(TObject *Sender,
      TShiftState Shift, int X, int Y)
{
 if(Timer->Enabled)
  Image->ShowHint=false;
 else
  Image->ShowHint=true;
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::ImageMouseDown(TObject *Sender,
      TMouseButton Button, TShiftState Shift, int X, int Y)
{
 int j;

 if(Timer->Enabled)
  return;
 for(j=0;j<Nbmark;j++)
 {
  if(abs(evt[j].pos*INC/Dsc.frequence-Start+512-X)<2)
  {
   evt[j].typeevt-=2;
   Mark=Eeg3_PutMarker(&evt[j]);
   Start-=512;
   Start=Start/128*128;
   ScanMark2();
   _run=0;
  }
 }
 ActiveControl=Edit;
}
//---------------------------------------------------------------------------


void __fastcall TForm3le::UpDChanClick(TObject *Sender, TUDBtnType Button)
{

 Chan0=UpDChan->Position*8;
 Chanmax=Chan0+8;
 if(Chanmax>Dsc.totalvoies)
  Chanmax=Dsc.totalvoies;
 Start-=512;
 if(Start<0)
  Start=0;
 Start=Start/128*128;
 ScanMark2();
 _run=0;
 EditChan->Text=AnsiString(Chan0+1)+"-"+Chanmax;
}


//---------------------------------------------------------------------------
void __fastcall TForm3le::FormCreate(TObject *Sender)
{
  DoubleBuffered=true; // YA 25/08/2007 : pour virer le scintillement
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::SpeedButtonOnlineClick(TObject *Sender)
{
  NbPointsRead=0;
  Timer2->Enabled=SpeedButtonOnline->Down;
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::Timer2Timer(TObject *Sender)
{
  int ret=-1, NbPointsToRead, NbS;

  NbPointsToRead=32;

  Eeg3_CloseFile();
  Eeg3_OpenFile(fileName,&Dsc);

  if(cbGetEeg2->Checked)
    ret=Eeg3_GetEeg2(NbPointsRead,NbPointsToRead,eeg);
  else
    ret=Eeg3_GetEeg(NbPointsRead,NbPointsToRead,Heeg);

  if (ret>0)
  {
    NbPointsRead+=ret;
    NbS=NbPointsRead / Dsc.frequence;
    NbPts->Caption=IntToStr(NbPointsRead);
    NbSecs->Caption=IntToStr(NbS);
  } else {
    if(ret==-106){ //fin de fichier
      nextFile(1);
    }
  }
}

//---------------------------------------------------------------------------

void __fastcall TForm3le::updateAffichage(){ //WR 17/11/2008
  int i;
  String* tabBloc;

  Memo->Clear();

  //Finding and displaying patient informations
  Memo->Lines->Add("Patient informations :");
  int ret = Eeg3_GetPatientInfo(&infoPat);
  if(ret<0){
    Memo->Lines->Add("Error reading patient informations : "+ret);
  } else {
    Memo->Lines->Add("Name : " + AnsiString(infoPat.firstname) + " " + AnsiString(infoPat.name));
    Memo->Lines->Add("Date of birth : " + AnsiString(infoPat.date));
    Memo->Lines->Add("Sex : " + AnsiString(infoPat.sex));
    Memo->Lines->Add("Origin and file number : " + AnsiString(infoPat.center) + ", " + AnsiString(infoPat.file));
    Memo->Lines->Add("Commentary : " + AnsiString(infoPat.commentary));
  }

  //Finding and displaying number of bloc
  // + the name of all available blocs
  /*int NumberOfBloc = Eeg3_GetNumberOfBlocFiles(fileName);

  if (NumberOfBloc>0){

    Memo->Lines->Add(IntToStr(NumberOfBloc)+" blocs in the file : ");

    tabBloc = new String[NumberOfBloc];

    for(i=0;i<NumberOfBloc;i++){
      tabBloc[i]="";
    }

    ret = Eeg3_GetAvailableBlocFiles(fileName, tabBloc);

    for(i=0;i<NumberOfBloc;i++){
      Memo->Lines->Add(tabBloc[i]);
    }

  }
  Memo->Lines->Add(" "); */

  //Displaying file informations
  Memo->Lines->Add(fileName);

  //Memo->Lines->Add(AnsiString(Dsc.date));

  int Hour, Min, Sec, Day, Month, Year;
  char* szTestDate;
  char* szTestTime;
  szTestDate=new char[20];
  szTestTime=new char[20];
  if (6 == sscanf(Dsc.date, "%d:%d:%d %d/%d/%d", &Hour, &Min, &Sec, &Day, &Month, &Year))
  {
    sprintf(szTestDate, "%02d/%02d/%04d", Day, Month, Year);
    sprintf(szTestTime, "%02d:%02d:%02d", Hour, Min, Sec);
  }

  Memo->Lines->Add(szTestDate);
  Memo->Lines->Add(szTestTime);

  Memo->Lines->Add(AnsiString(Duree=Dsc.duree)+" seconds    "+Dsc.frequence+" Hz");
  Memo->Lines->Add(AnsiString(Dsc.totalvoies)+" chan");

  UpDChan->Max=(short)((Dsc.totalvoies-1)/8);
  Chanmax=Dsc.totalvoies>8?8:Dsc.totalvoies;
  EditChan->Text=AnsiString(Chan0+1)+"-"+Chanmax;

  for(int j=0;j<Dsc.totalvoies;j++)
  {
    Dsc.nom[j][7]=0;
    while(Dsc.nom[j][strlen(Dsc.nom[j])-1]==' ')
      Dsc.nom[j][strlen(Dsc.nom[j])-1]=0;
    Memo->Lines->Text=Memo->Lines->Text+(j?",":"")+Dsc.nom[j];
  }

  //displaying events

  //the function should go through all the blocs file that can be opened to read their events

  Memo->Lines->Add(" ");
  Memo->Lines->Add("Events :");

  int NumberOfFiles = Eeg3_GetNumberOfBlocFiles(fileName); // = eeg3NumberOfFiles();

  char* temp = new char[260];
  //memcpy(oldFileName, fileName, strlen(fileName)+1);
  int startPosForCurrentBloc = 0;

  for (int j=0; j<NumberOfFiles; j++)
  {

    int Count = Eeg3_GetMarkersNumber(0, MAXINT); //INT_MAX
    //if (Count <= 0)
      //return 0

    if (Count>0){

      TEventLe* pEvents = new TEventLe[Count];

      Eeg3_GetMarkers2(0, MAXINT, pEvents); //INT_MAX

      for (int i = 0; i < Count; i++)
      {
        TEventLe* p = &pEvents[i];
        Memo->Lines->Add(IntToStr( (p->pos+startPosForCurrentBloc)/Dsc.frequence)+" "+IntToStr(p->duree/Dsc.frequence)+" "+p->buf);
      }

      delete pEvents;
      
    }

    if (j<NumberOfFiles-1)
    {
      startPosForCurrentBloc += Dsc.duree*Dsc.frequence;
      Eeg3_NextFile(1,temp,&Dsc);
    }

  }

  // the first bloc is re-opened
  //memcpy(fileName, oldFileName, strlen(oldFileName)+1);
  Eeg3_OpenFile(fileName, &Dsc);

  // return Count

}

//---------------------------------------------------------------------------

int __fastcall TForm3le::nextFile(int direction){ //WR 18/11/2008

  char* temp = new char[260];
  int ret = Eeg3_NextFile(direction, temp, &Dsc);
  if(ret==0){
    fileName = temp;
    NbPointsRead = 0;
    updateAffichage();
  }
  return ret;
}

//---------------------------------------------------------------------------

void __fastcall TForm3le::btPrevClick(TObject *Sender)
{
  nextFile(-1);
}

//---------------------------------------------------------------------------

void __fastcall TForm3le::btNextClick(TObject *Sender)
{
  nextFile(1);
}

//---------------------------------------------------------------------------

//WR 26/11/2008
void __fastcall TForm3le::btGetEventClick(TObject *Sender)
{
  //to be sure that it's the last version of the file
  Eeg3_CloseFile();
  Eeg3_OpenFile(fileName,&Dsc);

  //update display
  updateAffichage();

  //recuperation of the values and exeptions
  int debut, fin;

  try {
    debut = StrToInt(edDeb->Text);
  } catch(Exception &e) {
    debut = 0;
    edDeb->Text = "0";
  }
  if(debut<0){
    debut = 0;
    edDeb->Text = "0";
  }

  try {
    fin = StrToInt(edFin->Text);
  } catch(Exception &e) {
    fin = Dsc.duree;
    edFin->Text = IntToStr(Dsc.duree);
  }
  if(fin>Dsc.duree){
    fin = Dsc.duree;
    edFin->Text = IntToStr(Dsc.duree);
  }

  //recuperation of the events
  int nbEvt = Eeg3_GetMarkersNumber(debut*Dsc.frequence,fin*Dsc.frequence);

  Memo->Lines->Add(" ");

  Memo->Lines->Add(IntToStr(nbEvt) + " marker(s) found between point " + IntToStr(debut*Dsc.frequence)
    + " (" + IntToStr(debut) + "s) and " + IntToStr(fin*Dsc.frequence) + " (" + IntToStr(fin) + "s)." );

  if(nbEvt>0){

    if(UseMarkersLongCb->Checked)
    {
      TEventLeLong* myEvt = new TEventLeLong[nbEvt];

      Eeg3_GetMarkersLong2(debut*Dsc.frequence,fin*Dsc.frequence,myEvt);

      for(int i=0;i<nbEvt;i++){
        Memo->Lines->Add("Marker n°" + IntToStr(i+1) + " : " + myEvt[i].buf + ", pos : "
          + myEvt[i].pos + " (" + IntToStr(int(myEvt[i].pos/Dsc.frequence)) + "s)." );
      }
    }
    else
    {
      TEventLe* myEvt = new TEventLe[nbEvt];

      Eeg3_GetMarkers2(debut*Dsc.frequence,fin*Dsc.frequence,myEvt);

      for(int i=0;i<nbEvt;i++){
        Memo->Lines->Add("Marker n°" + IntToStr(i+1) + " : " + myEvt[i].buf + ", pos : "
          + myEvt[i].pos + " (" + IntToStr(int(myEvt[i].pos/Dsc.frequence)) + "s)." );
      }
    }
  }

 return;

}
//---------------------------------------------------------------------------

//WR 26/11/2008
void __fastcall TForm3le::btGetEvent2Click(TObject *Sender)
{
  //to be sure that its the last version of the file
  Eeg3_CloseFile();
  Eeg3_OpenFile(fileName,&Dsc);

  //update display
  updateAffichage();

  //Récupération of the value in the field and exceptions
  int debut, fin;

  try {
    debut = StrToInt(edDeb2->Text);
  } catch(Exception &e) {
    debut = 0;
    edDeb2->Text = "0";
  }
  if(debut<0){
    debut = 0;
    edDeb2->Text = "0";
  }

  try {
    fin = StrToInt(edFin2->Text);
  } catch(Exception &e) {
    fin = Dsc.duree*Dsc.frequence;
    edFin2->Text = IntToStr(Dsc.duree*Dsc.frequence);
  }
  if(fin>Dsc.duree*Dsc.frequence){
    fin = Dsc.duree*Dsc.frequence;
    edFin2->Text = IntToStr(Dsc.duree*Dsc.frequence);
  }

  //event recuperation
  int nbEvt = Eeg3_GetMarkersNumber(debut,fin);

  Memo->Lines->Add(" ");

  Memo->Lines->Add(IntToStr(nbEvt) + " marker(s) found between point " + IntToStr(debut) + " ("
    + IntToStr(int(debut/Dsc.frequence)) + "s) and " + IntToStr(fin)
    + " (" + IntToStr(int(fin/Dsc.frequence)) + "s) :" );

  if(nbEvt>0){

    if(UseMarkersLongCb->Checked)
    {
      TEventLeLong* myEvt = new TEventLeLong[nbEvt];

      Eeg3_GetMarkersLong2(debut,fin,myEvt);

      for(int i=0;i<nbEvt;i++){
        Memo->Lines->Add("- Marker n°" + IntToStr(i+1) + " : " + myEvt[i].buf + ", pos : "
          + myEvt[i].pos + " (" + IntToStr(int(myEvt[i].pos/Dsc.frequence)) + "s)." );
      }
    }
    else
    {
      TEventLe* myEvt = new TEventLe[nbEvt];

      Eeg3_GetMarkers2(debut,fin,myEvt);

      for(int i=0;i<nbEvt;i++){
        Memo->Lines->Add("- Marker n°" + IntToStr(i+1) + " : " + myEvt[i].buf + ", pos : "
          + myEvt[i].pos + " (" + IntToStr(int(myEvt[i].pos/Dsc.frequence)) + "s)." );
      }
    }
  }

 return;
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::btMaxClick(TObject *Sender)
{
  //pour être sûr d'avoir la dernière version du fichier
  Eeg3_CloseFile();
  Eeg3_OpenFile(fileName,&Dsc);

  edFin->Text = IntToStr(Dsc.duree);
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::btMax2Click(TObject *Sender)
{
  //pour être sûr d'avoir la dernière version du fichier
  Eeg3_CloseFile();
  Eeg3_OpenFile(fileName,&Dsc);

  edFin2->Text = IntToStr(Dsc.duree*Dsc.frequence);
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::FirstBtClick(TObject *Sender)
{
  int NumberOfBloc = Eeg3_GetNumberOfBlocFiles(fileName);
  String* tabBloc;
  int i,ret;

  if (NumberOfBloc>0){

    tabBloc = new String[NumberOfBloc];

    for(i=0;i<NumberOfBloc;i++){
      tabBloc[i]="";
    }

    ret = Eeg3_GetAvailableBlocFiles(fileName, tabBloc);

    Eeg3_CloseFile();
    fileName = tabBloc[0].c_str();
    Eeg3_OpenFile(fileName,&Dsc);
    NbPointsRead = 0;
    updateAffichage();
  }
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::lastBtClick(TObject *Sender)
{
  int NumberOfBloc = Eeg3_GetNumberOfBlocFiles(fileName);
  String* tabBloc;
  int i,ret;

  if (NumberOfBloc>0){

    tabBloc = new String[NumberOfBloc];

    for(i=0;i<NumberOfBloc;i++){
      tabBloc[i]="";
    }

    ret = Eeg3_GetAvailableBlocFiles(fileName, tabBloc);

    Eeg3_CloseFile();
    fileName = tabBloc[NumberOfBloc-1].c_str();
    Eeg3_OpenFile(fileName,&Dsc);
    NbPointsRead = 0;
    updateAffichage();
  }
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::TermBtClick(TObject *Sender)
{
  if(hlib)
  {
    //close
    Eeg3_CloseFile();
    Eeg3_Termination();
    FreeLibrary(hlib);
    hlib=0;

    TermBt->Caption = "Activate library";

    //disable buttons
    btPrev->Enabled=false;
    btNext->Enabled=false;
    SpeedButtonOnline->Enabled=false;
    ButtonMk->Enabled=false;
    ButtonEv->Enabled=false;
    ButtonGo->Enabled=false;
    btGetEvent->Enabled=false;
    btMax->Enabled=false;
    btGetEvent2->Enabled=false;
    btMax2->Enabled=false;
  }
  else
  {
    //open
    hlib=DLLinit();
    if(!hlib)                     // DLL introuvable
    {
      Application->MessageBox("Can't load coherence5le.dll","Error",
                  MB_SYSTEMMODAL|MB_OK|MB_ICONHAND);
      Application->Terminate();
      return;
    }
    else {
      Eeg3_DebugFileSwitch(true); // Display debug file
      Eeg3_TempFolderSwitch(true); //WR 07/06/10
      result=Eeg3_Initialisation();
    }
    if(result<0)                  // Problem with the DLL
    {
      char buf[80];

      sprintf(buf,"Error %d",result);
      Application->MessageBox("Can't use coherence5le.dll",buf,
                  MB_SYSTEMMODAL|MB_OK|MB_ICONHAND);
      Application->Terminate();
      return;
    }

    TermBt->Caption = "Terminate library";

    //enable buttons
    btPrev->Enabled=true;
    btNext->Enabled=true;
    SpeedButtonOnline->Enabled=true;
    ButtonMk->Enabled=true;
    ButtonEv->Enabled=true;
    ButtonGo->Enabled=true;
    btGetEvent->Enabled=true;
    btMax->Enabled=true;
    btGetEvent2->Enabled=true;
    btMax2->Enabled=true;

    //load a file
    int j;
    TUnlock3LE UL3LE;

    Eeg3_Version(&Ver);
    Memo->Lines->Add(AnsiString("Coherence5LE.dll : ")+Ver.major+"."+Ver.minor+"."+
                                     Ver.compile+"."+Ver.number+"\r\n");
    if(!OpD->Execute())
    {
    Application->Terminate();
    return;
    }

    fileName = OpD->FileName.c_str();

    UL3LE.int1=UNLOCK1;
    UL3LE.int2=UNLOCK2;
    UL3LE.int3=UNLOCK3;
    UL3LE.int4=UNLOCK4;

    if (Eeg3_Unlock(UL3LE)<0)
    {
    Application->MessageBox("Information","TestCoh5LE.exe has failed to unlock Coherence5LE.dll. Execution will stop.",
             MB_SYSTEMMODAL|MB_OK|MB_ICONHAND);
    Application->Terminate();
    return;
    }

    j=Eeg3_OpenFile(fileName,&Dsc);

    if (j==-1)
    {
    Application->MessageBox("Coherence5LE.dll is locked, you must unlock it.",fileName,
             MB_SYSTEMMODAL|MB_OK|MB_ICONHAND);
    Application->Terminate();
    return;
    }
    if (j<0)
    {
    Application->MessageBox(fileName,"Can't open",
             MB_SYSTEMMODAL|MB_OK|MB_ICONHAND);
    Application->Terminate();
    return;
    }

    ScanMark();
    Bitmap=new Graphics::TBitmap;
    Bitmap->Width=512;
    Bitmap->Height=512;
    Image->Canvas->Draw(0,0,Bitmap);

    updateAffichage();
    edFin->Text = IntToStr(Dsc.duree);
    edFin2->Text = IntToStr(Dsc.duree*Dsc.frequence);
    Heeg=GlobalAlloc(GMEM_MOVEABLE|GMEM_ZEROINIT,Dsc.totalvoies*Dsc.frequence*sizeof(short)*2); //two seconds readed
    eeg=new short[Dsc.totalvoies*Dsc.frequence*2]; //two seconds readed
  }

}
//---------------------------------------------------------------------------

//Sample code to use GetNextMarker
/*void __fastcall TForm3le::Button1Click(TObject *Sender)
{
  //to be sure that its the last version of the file
  Eeg3_CloseFile();
  Eeg3_OpenFile(fileName,&Dsc);

  //update display
  updateAffichage();

  TEventLe evt;
  int start=0;
  int retour=0;

  while (retour>=0) {
    retour = Eeg3_GetNextMarker(start,&evt);
    start=evt.pos+1;
  }
}*/
//---------------------------------------------------------------------------


void __fastcall TForm3le::Button1Click(TObject *Sender)
{
  int ret=-1;

  Eeg3_CloseFile();
  Eeg3_OpenFile(fileName,&Dsc);

  ret=Eeg3_GetEeg(68188-100,32,Heeg);
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::GetRealTimeBtClick(TObject *Sender)
{

  if(UseMarkers2Cb->Checked)
  {
    int NbRTE=Eeg3_GetNumberOfRealTimeMarkers();
    if(NbRTE>0){
      TRealTimeMarker *myRTM=new TRealTimeMarker[NbRTE];
      Eeg3_GetRealTimeMarkers2(myRTM);
      for(int i=0;i<NbRTE;i++)
      {
        Memo->Lines->Add("Real time "+IntToStr(i)+" : at pos "+IntToStr(myRTM[i].pos)+" the time is "+IntToStr(myRTM[i].realtime));
      }
    }
  }
  else
  {
    HGLOBAL HRTM;
    TRealTimeMarker *myRTM;
    int NbRTE;

    GlobalUnlock(HRTM);

    NbRTE=Eeg3_GetRealTimeMarkers(&HRTM);
    NbRTE=GlobalSize(HRTM)/sizeof(TRealTimeMarker);
    myRTM=(TRealTimeMarker*)GlobalLock(HRTM);
    for(int i=0;i<NbRTE;i++)
    {
      Memo->Lines->Add("Real time "+IntToStr(i)+" : at pos "+IntToStr(myRTM[i].pos)+" the time is "+IntToStr(myRTM[i].realtime));
    }
  }
}
//---------------------------------------------------------------------------

void __fastcall TForm3le::GetSAndDBtClick(TObject *Sender)
{
  TVideoShiftAndDrift myVSAD;
  myVSAD.Shift=0;
  myVSAD.Drift=0;

  int retour=Eeg3_GetVideoShiftAndDrift(&myVSAD);

  if(retour>=0)
  {
    Memo->Lines->Add("-Video Shift="+IntToStr(myVSAD.Shift));
    Memo->Lines->Add("-Video Drift="+IntToStr(myVSAD.Drift));
  }
}
//---------------------------------------------------------------------------


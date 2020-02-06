//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop
USEFORM("mainletest.cpp", Form3le);
//---------------------------------------------------------------------------

#include "mainletest.h"

int result;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{
    try
    {
        Application->Initialize();
        Application->CreateForm(__classid(TForm3le), &Form3le);
        Application->Run();

        if(hlib)
        {
         Eeg3_Termination();
         FreeLibrary(hlib);
        }
    }
    catch (Exception &exception)
    {
        Application->ShowException(&exception);
    }
    return 0;
}
//---------------------------------------------------------------------------

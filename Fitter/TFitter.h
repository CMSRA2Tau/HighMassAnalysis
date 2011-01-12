/*********************************************************
********************************************************** 
***  April 23,2010					**
***  Indara Suarez					**
***	define Class TFitter  				**
***  version 1.0			     		**
**********************************************************
*********************************************************/
#include <string>

using namespace std;

#include "TString.h"
#include "TH1F.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "TPrProp.h"
#include "TChProp.h"
#ifndef TFitter_h
#define TFitter_h


class TFitter {
private:
TPrProp *fProcess[15];
TChProp *fChannel;

public :

TFitter(){
  _theNExp = 300;
  _theLumi = 50.;
  _theNMCInt = 10;
  _theMass = 500.;
  _theXSection = 1.914;
  
  GetTheFiles();
  BookCLHistos();
};
virtual ~TFitter(){};

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   virtual void     Loop();
   void SetNExp(int);
   void SetLumi(float);
   void SetLumiErr(float);
   void SetSystSettings(int);
   void NMCInt(int);
   void SetSignalMassXsection(float, float);
   //void AddTree(std::string);

private:
  int _theNExp;
  float _theLumi;
  float _theLumiErr;
  int _theNMCInt;
  int _SystSettings;
  float _theMass;
  float _theXSection;
  std::map<int, TH1F*> _theNameHistoMap;
  std::vector<std::string> _theRootFiles;
  void BookCLHistos();
  void GetTheFiles();

};

#endif


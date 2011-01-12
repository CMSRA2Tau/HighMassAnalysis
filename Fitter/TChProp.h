/*********************************************************
********************************************************** 
***  October 10,2010					**
***  Indara Suarez					**
***  Texas A&M University				**
***	   TChProp constructor				**
***	version 2.0        				**
**********************************************************
*********************************************************/
#include <string>

using namespace std;
#include "TString.h"
#include "TH1F.h"
#include "TPrProp.h"
#include <TFile.h>

#ifndef __TChProp_h__
#define __TChProp_h__


class TChProp{
 private:  
  
  string fchannelt;
  
  Int_t fNPr;
  Float_t Lumi2;
  Float_t Lumi2err;
  
  TH1F* fSignalSh;
  TH1F* fTotBGSh;
  
  TPrProp fProcessN[15];
  TPrProp fSignal[1];
  
  
  TH1F* fLikelihood;
  TH1F* f95CL;
  
 public:
  
  TChProp(){};
  virtual ~TChProp(){};
  
  void SetSignalShape(TH1F* th0){fSignalSh = th0;};
  void SetTotBGShape(TH1F* th1){fTotBGSh = th1;};
  void SetNProcess(Int_t tch0){fNPr = tch0;};
  void SetChannelTitle(string fpr1){fchannelt=fpr1;};
  void SetLimit95CLH(TH1F* th2){f95CL = th2;};
  void SetChLumi(Float_t tch3){Lumi2 = tch3;};
  void SetChLumiErr(Float_t tch5){Lumi2err = tch5;};

  void SetPrProp(int ichnl, TFile* fntuple);
  void SetProcessProp(TPrProp* fprocess4, int ipr);
  void SetSignalProp(TPrProp* fprocess3);
  
  Float_t GetChLumi(){return Lumi2;};
  Float_t GetChLumiErr(){return Lumi2err;};
  TH1F* GetSignalShape(){return fSignalSh;};
  TH1F* GetTotBGShape(){return fTotBGSh;};
  TH1F* GetLimit95CLH(){return f95CL;};
  
  TPrProp GetSignalProp(){return fSignal[0];};
  TPrProp GetProcessProp(int p){return fProcessN[p];};
  
  Int_t GetNProcess(){return fNPr;};
  string GetChannelTitle(){return fchannelt;};
  
  TH1F* Likelihood(int setting0, TChProp fChannelPro, int mcintg,TH1F* pseudo, TH1F* LogLVsSigma1);
  Double_t Limit95CL(TH1F* LogLVsSigma);
  
  
};

#endif

//
// TChProp
//
#include "TH1F.h"
#include "TPrProp.h"

#ifndef __TChProp_h__
#define __TChProp_h__


class TChProp{
private:  
  
  string fchannelt;
  
  Int_t fNPr;

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
  
  void SetProcessProp(TPrProp* fprocess, int ipr);
  void SetSignalProp(TPrProp* fprocess3);
  

  TH1F* GetSignalShape(){return fSignalSh;};
  TH1F* GetTotBGShape(){return fTotBGSh;};
  TH1F* GetLimit95CLH(){return f95CL;};
  
  TPrProp GetSignalProp(){return fSignal;};
  TPrProp GetProcessProp(int p){return fProcessN[p];};
  
  Int_t GetNProcess(){return fNPr;};
  string GetChannelTitle(){return fchannelt;};
  
  TH1F* Likelihood(int ichn, TChProp fChannelPro, int mcintg, TH1F* LogLVsSigma1);
  Double_t Limit95CL(TH1F* LogLVsSigma);
  
  
};

#endif

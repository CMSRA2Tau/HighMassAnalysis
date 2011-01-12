/*********************************************************
**********************************************************
***  October 10,2010					**
***  Indara Suarez					**
***  Texas A&M University				**
***	   TPrProp constructor  			**
***	version 2.0    	    				**
**********************************************************
*********************************************************/
#include <string>


using namespace std;
#include <vector>

#include "TString.h"
#include "TH1F.h"
#ifndef __TPrProp_h__
#define __TPrProp_h__


class TPrProp{
 private:  
  
  Float_t	fTotEff;
  Float_t	fPreScl;
  
  Double_t	fLumi;
  Double_t	fLumierr;
  Double_t	fSigma;
  Float_t	fBRatio;
  Float_t	fBRatioerr;
  
  TH1F* fMDTShapeProcess;
  TH1F* fSystematicShape[6];
  TH1F* fPDFShapeSys;
  TH1F* fJetShapeSys;
  TH1F* fTauEtaShapeSys; 

  Float_t ftotscl;
  
  Float_t fFilterEff;
  Float_t fAccEff;
  Float_t feIDEff;
  Float_t fmuIDEff;
  Float_t ftauIDEff;
  Float_t fTopoEff;
  
  Float_t fBFracStatErr ;
  Float_t fAccEffStatErr;
  Float_t feIDEffStatErr;
  Float_t fmuIDEffStatErr;
  Float_t ftauIDStatErr ;
  Float_t fTopoEffStatErr;
  Float_t fBFracSysErr ;
  Float_t fAccEffSysErr;
  Float_t feIDEffSysErr;
  Float_t fmuIDEffSysErr;
  Float_t ftauIDSysErr ;
  Float_t fTopoEffSysErr;
  Float_t fCmErrors[10];
  vector<float> *fSystematics; 
  
  
  string fprocessname;

public:
  TPrProp(){};
  virtual ~TPrProp(){};
  
  void SetProcessTitle(string fprname){fprocessname =(string)fprname;};
  void SetTotEff(Float_t tch0){fTotEff = tch0;};
  void SetPreScl(Float_t tch1){fPreScl = tch1;};
  
  void SetFilterEfficiency(Float_t eff0){fFilterEff = eff0;};
  void SetacceptanceEfficiency(Float_t eff1){fAccEff = eff1;};
  void SetelectronIDEfficiency(Float_t eff2){feIDEff = eff2;};
  void SetmuonIDEfficiency(Float_t eff3){fmuIDEff = eff3;};
  void SettauIDEfficiency(Float_t eff4){ftauIDEff = eff4;};
  void SettopologyEfficiency(Float_t eff5){fTopoEff = eff5;};
  
  void SetSystematics(vector<float> *systerrs){fSystematics=systerrs;};

  void SetbranchingFractionStatError(Float_t eff01){fBFracStatErr = eff01;};  
  void SetacceptanceEfficiencyStatError(Float_t eff02){fAccEffStatErr = eff02;}; 
  void SetelectronIDEfficiencyStatError(Float_t eff03){feIDEffStatErr = eff03;}; 
  void SetmuonIDEfficiencyStatError(Float_t eff04){fmuIDEffStatErr = eff04;};
  void SettauIDEfficiencyStatError(Float_t eff05){ftauIDStatErr = eff05;};
  void SettopologyEfficiencyStatError(Float_t eff06){fTopoEffStatErr = eff06;};
  
  void SetbranchingFractionSystError(Float_t eff07){fBFracSysErr = eff07;};  
  void SetacceptanceEfficiencySystError(Float_t eff08){fAccEffSysErr = eff08;};
  void SetelectronIDEfficiencySystError(Float_t eff09){feIDEffSysErr = eff09;};
  void SetmuonIDEfficiencySystError(Float_t eff010){fmuIDEffSysErr = eff010;};
  void SettauIDEfficiencySystError(Float_t eff011){ftauIDSysErr = eff011;};
  void SettopologyEfficiencySystError(Float_t eff012){fTopoEffSysErr = eff012;};
  
  void SetLuminosity(Double_t lumi){fLumi = lumi;};
  void SetLuminosityErr(Double_t lumir){fLumierr = lumir;};
  void SetSigma(Double_t sigma){fSigma = sigma;};
  void SetBR(Float_t bratio){fBRatio=bratio;};
  
  void SetProcessShape(TH1F* th4){fMDTShapeProcess = th4;};
  void SetSystematicShape(int m, TH1F* th5){fSystematicShape[m] = th5;};
  void SetPDFSysShape(TH1F* th01){fPDFShapeSys = th01;};
  void SetJETSysShape(TH1F* th02){fJetShapeSys = th02;};
  void SetTauEtaSysShape(TH1F* th03){fTauEtaShapeSys = th03;};
  
  string GetProcessTitle(){return fprocessname;};
  

  Float_t GetTotEff(){return fTotEff;};
  Float_t GetPreScl(){return fPreScl;};
  
  Float_t GetFilterEfficiency(){return fFilterEff;};
  Float_t GetacceptanceEfficiency(){return fAccEff;};
  Float_t GetelectronIDEfficiency(){return feIDEff;};
  Float_t GetmuonIDEfficiency(){return fmuIDEff;};
  Float_t GettauIDEfficiency(){return ftauIDEff;};
  Float_t GettopologyEfficiency(){return fTopoEff;};
  vector<float> * GetSystematics(){return fSystematics;};  
  
  Float_t GetbranchingFractionStatError(){return fBFracStatErr;};  
  Float_t GetacceptanceEfficiencyStatError(){return fAccEffStatErr;}; 
  Float_t GetelectronIDEfficiencyStatError(){return feIDEffStatErr;}; 
  Float_t GetmuonIDEfficiencyStatError(){return fmuIDEffStatErr;};
  Float_t GettauIDEfficiencyStatError(){return ftauIDStatErr;};
  Float_t GettopologyEfficiencyStatError(){return fTopoEffStatErr;};
  
  Float_t GetbranchingFractionSystError(){return fBFracSysErr;};  
  Float_t GetacceptanceEfficiencySystError(){return fAccEffSysErr;};
  Float_t GetelectronIDEfficiencySystError(){return feIDEffSysErr;};
  Float_t GetmuonIDEfficiencySystError(){return fmuIDEffSysErr;};
  Float_t GettauIDEfficiencySystError(){return ftauIDSysErr;};
  Float_t GettopologyEfficiencySystError(){return fTopoEffSysErr;};
  

  Double_t GetLuminosity(){return fLumi;};
  Double_t GetLuminosityErr(){return fLumierr;};
  Double_t GetSigma(){return fSigma;};
  Float_t GetBR(){return fBRatio;};
  
  
  TH1F* GetProcessShape(){return fMDTShapeProcess;};
  TH1F* GetSystematicShape(int n){return fSystematicShape[n];};
  TH1F* GetPDFSysShape(){return fPDFShapeSys;};
  TH1F* GetJETSysShape(){return fJetShapeSys;};
  TH1F* GettauetaSysShape(){return fTauEtaShapeSys;};
  

   

 
};

#endif

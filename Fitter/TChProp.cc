/*********************************************************
********************************************************** 
***  October 10,2010					**
***  Indara Suarez					**
***  Texas A&M University				**
***	runs with TFitter.cc			    	**
*** 	fills TPrProp objects 			   	**
*** 	Calculates single channel likelihood 		**
***	Calculates joint likelihood			**
***  version 2.0				 	**
**********************************************************
*********************************************************/

#include <math.h>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <string.h>

#include <fstream>
using namespace std;
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TKey.h" 
#include "TH1.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TRandom.h"
#include "TVector.h"
#include "TLorentzVector.h" 
#include "TList.h" 
#include "TChain.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRint.h"

#include "TChProp.h"
#include "TPrProp.h"
#include "TSystematicStateGenerator.h"


#include <vector>


void TChProp::SetPrProp(int ichnl, TFile* fntuple){

  gROOT->ProcessLine("#include <vector>");
  
  Float_t  fbranchingFraction;
  Float_t  fFilterEfficiency;
  Float_t  fxsection;
  Float_t  facceptanceEfficiency;
  Float_t  felectronIDEfficiency;
  Float_t  fmuonIDEfficiency;
  Float_t  ftauIDEfficiency;
  Float_t  ftopologyEfficiency;
  vector<float>* fsystematics0 = new vector<float>;
  Float_t  fbranchingFractionStatError;
  Float_t  facceptanceEfficiencyStatError;
  Float_t  felectronIDEfficiencyStatError;
  Float_t  fmuonIDEfficiencyStatError;
  Float_t  ftauIDEfficiencyStatError;
  Float_t  ftopologyEfficiencyStatError;
  Float_t  fbranchingFractionSystError;
  Float_t  facceptanceEfficiencySystError;
  Float_t  felectronIDEfficiencySystError;
  Float_t  fmuonIDEfficiencySystError;
  Float_t  ftauIDEfficiencySystError;
  Float_t  ftopologyEfficiencySystError;

  vector<float> systv;
  vector<float> systv2[10];
  TKey *key;
  TKey *keyh;
  int nb=0;
  int ch=0;
  int ipr=0;  
  int RebinFactor =3;
  Float_t binerr=0.;
  Float_t peak=0.;
  Int_t nbinz2=0;
  
  TPrProp fProcess[15];
  TPrProp fProcess2[1];
  
  Float_t toteff=0;
  
  if ( fntuple->IsOpen() ) printf("File opened successfully\n"); 
  
  TDirectory *dirsav = gDirectory;
  TIter next(dirsav->GetListOfKeys());
  
  while ((key = (TKey*)next())) {
    string histname1 = key->GetName();
    cout<<"title1  "<<histname1<<"   ";
    
    TObject *obj = key->ReadObj();
    string histname2 = obj->ClassName();
    cout<<"type  "<<histname2<<endl;
    
    if ( obj->IsA()->InheritsFrom( "TTree" ) ) {
      string teffs = obj->GetName();
      TTree *tobj = (TTree*)obj;
      fsystematics0 = new vector<float>;
      Int_t Nentries = (Int_t)tobj->GetEntries();     
      tobj->SetBranchAddress("branchingFraction",	      &fbranchingFraction	       );   //&b_fbranchingFraction);
      tobj->SetBranchAddress("FilterEfficiency",	      &fFilterEfficiency	       );   //&b_fFilterEfficiency);
      tobj->SetBranchAddress("xsection",		      &fxsection		       );   //&b_fxsection);
      tobj->SetBranchAddress("acceptanceEfficiency",	      &facceptanceEfficiency	       );   //&b_facceptanceEfficiency);
      tobj->SetBranchAddress("electronIDEfficiency",	      &felectronIDEfficiency	       );   //&b_felectronIDEfficiency);
      tobj->SetBranchAddress("muonIDEfficiency",	      &fmuonIDEfficiency	       );   //&b_fmuonIDEfficiency);
      tobj->SetBranchAddress("tauIDEfficiency", 	      &ftauIDEfficiency 	       );   //&b_ftauIDEfficiency);
      tobj->SetBranchAddress("topologyEfficiency",	      &ftopologyEfficiency	       );   //&b_ftopologyEfficiency);
      tobj->SetBranchAddress("systematics",		      &fsystematics0                   );
      tobj->SetBranchAddress("branchingFractionStatError",    &fbranchingFractionStatError     );   //&b_fbranchingFractionStatError);
      tobj->SetBranchAddress("acceptanceEfficiencyStatError", &facceptanceEfficiencyStatError  );   //&b_facceptanceEfficiencyStatError);
      tobj->SetBranchAddress("electronIDEfficiencyStatError", &felectronIDEfficiencyStatError  );   //&b_felectronIDEfficiencyStatError);
      tobj->SetBranchAddress("muonIDEfficiencyStatError",     &fmuonIDEfficiencyStatError      );   //&b_fmuonIDEfficiencyStatError);
      tobj->SetBranchAddress("tauIDEfficiencyStatError",      &ftauIDEfficiencyStatError       );   //&b_ftauIDEfficiencyStatError);
      tobj->SetBranchAddress("topologyEfficiencyStatError",   &ftopologyEfficiencyStatError    );   //&b_ftopologyEfficiencyStatError);
      tobj->SetBranchAddress("branchingFractionSystError",    &fbranchingFractionSystError     );   //&b_fbranchingFractionSystError);
      tobj->SetBranchAddress("acceptanceEfficiencySystError", &facceptanceEfficiencySystError  );   //&b_facceptanceEfficiencySystError);
      tobj->SetBranchAddress("electronIDEfficiencySystError", &felectronIDEfficiencySystError  );   //&b_felectronIDEfficiencySystError);
      tobj->SetBranchAddress("muonIDEfficiencySystError",     &fmuonIDEfficiencySystError      );   //&b_fmuonIDEfficiencySystError);
      tobj->SetBranchAddress("tauIDEfficiencySystError",      &ftauIDEfficiencySystError       );   //&b_ftauIDEfficiencySystError);
      tobj->SetBranchAddress("topologyEfficiencySystError",   &ftopologyEfficiencySystError    );   //&b_ftopologyEfficiencySystError);
      
      for ( Int_t iev=0; iev<Nentries; iev++ ) {
	nb=tobj->GetEntry(iev);
	if(histname1=="ZprimeTree"){
	  toteff=fFilterEfficiency;
	  toteff = toteff*(facceptanceEfficiency * felectronIDEfficiency * fmuonIDEfficiency * ftauIDEfficiency * ftopologyEfficiency*fbranchingFraction);
	  fProcess2[0].SetTotEff(toteff);
	  fProcess2[0].SetSigma(fxsection);
	  fProcess2[0].SetBR(fbranchingFraction);
	  fProcess2[0].SetFilterEfficiency(fFilterEfficiency);
	  fProcess2[0].SetacceptanceEfficiency(facceptanceEfficiency);
	  fProcess2[0].SetelectronIDEfficiency(felectronIDEfficiency);
	  fProcess2[0].SetmuonIDEfficiency(fmuonIDEfficiency);
	  fProcess2[0].SettauIDEfficiency(ftauIDEfficiency);
	  fProcess2[0].SettopologyEfficiency(ftopologyEfficiency);
	  
	  fProcess2[0].SetSystematics(fsystematics0);
	  fProcess2[0].SetbranchingFractionStatError(fbranchingFractionStatError);  
	  fProcess2[0].SetacceptanceEfficiencyStatError(facceptanceEfficiencyStatError);  
	  fProcess2[0].SetelectronIDEfficiencyStatError(felectronIDEfficiencyStatError);  
	  fProcess2[0].SetmuonIDEfficiencyStatError(fmuonIDEfficiencyStatError);	 
	  fProcess2[0].SettauIDEfficiencyStatError(ftauIDEfficiencyStatError);	
	  fProcess2[0].SettopologyEfficiencyStatError(ftopologyEfficiencyStatError);    
	  fProcess2[0].SetbranchingFractionSystError(fbranchingFractionSystError);	  
	  fProcess2[0].SetacceptanceEfficiencySystError(facceptanceEfficiencySystError);  
	  fProcess2[0].SetelectronIDEfficiencySystError(felectronIDEfficiencySystError);  
	  fProcess2[0].SetmuonIDEfficiencySystError(fmuonIDEfficiencySystError);	 
	  fProcess2[0].SettauIDEfficiencySystError(ftauIDEfficiencySystError);	
	  fProcess2[0].SettopologyEfficiencySystError(ftopologyEfficiencySystError);
	}
	
	if(histname1!="ZprimeTree"){
	  toteff=fFilterEfficiency;
	  toteff = toteff*(facceptanceEfficiency * felectronIDEfficiency * fmuonIDEfficiency * ftauIDEfficiency * ftopologyEfficiency*fbranchingFraction);
	  fProcess[ch].SetTotEff(toteff);
	  fProcess[ch].SetSigma(fxsection);
	  fProcess[ch].SetBR(fbranchingFraction);
	  fProcess[ch].SetFilterEfficiency(fFilterEfficiency);
	  fProcess[ch].SetacceptanceEfficiency(facceptanceEfficiency);
	  fProcess[ch].SetelectronIDEfficiency(felectronIDEfficiency);
	  fProcess[ch].SetmuonIDEfficiency(fmuonIDEfficiency);
	  fProcess[ch].SettauIDEfficiency(ftauIDEfficiency);
	  fProcess[ch].SettopologyEfficiency(ftopologyEfficiency);
	  
	  fProcess[ch].SetSystematics(fsystematics0);
	  
	  fProcess[ch].SetbranchingFractionStatError(fbranchingFractionStatError);  
	  fProcess[ch].SetacceptanceEfficiencyStatError(facceptanceEfficiencyStatError);  
	  fProcess[ch].SetelectronIDEfficiencyStatError(felectronIDEfficiencyStatError);  
	  fProcess[ch].SetmuonIDEfficiencyStatError(fmuonIDEfficiencyStatError);	   
	  fProcess[ch].SettauIDEfficiencyStatError(ftauIDEfficiencyStatError);	  
	  fProcess[ch].SettopologyEfficiencyStatError(ftopologyEfficiencyStatError);    
	  fProcess[ch].SetbranchingFractionSystError(fbranchingFractionSystError);     
	  fProcess[ch].SetacceptanceEfficiencySystError(facceptanceEfficiencySystError);  
	  fProcess[ch].SetelectronIDEfficiencySystError(felectronIDEfficiencySystError);  
	  fProcess[ch].SetmuonIDEfficiencySystError(fmuonIDEfficiencySystError);	   
	  fProcess[ch].SettauIDEfficiencySystError(ftauIDEfficiencySystError);	  
	  fProcess[ch].SettopologyEfficiencySystError(ftopologyEfficiencySystError);
	  ch=ch+1;
	}
	SetNProcess(ch);
      }
    }  
    
    int histn=0;
    int histn2=0;  
    int ipr2=0; 
    
    if (obj->IsA()->InheritsFrom("TDirectoryFile")) {
      TDirectory *dir = (TDirectory*)obj;
      dir->cd(dir->GetPath());
      TIter nexth(dir->GetListOfKeys());
      histn=0;
      
      while ((keyh = (TKey*)nexth())) {
	TObject *objh = keyh->ReadObj();
	
	if ( objh->IsA()->InheritsFrom( "TH1F" ) ) {
	  string histnameh = objh->GetName();
	  TH1F* hobj = (TH1F*)objh->Clone("hobj");
	  fProcess[ipr].SetProcessTitle(histname1);
	  
	  if(histname1=="ZprimeTemplateDirectory"){ 
	    
	    if(histnameh=="Template"){
	      fProcess2[0].SetProcessShape(hobj);
	      if (fProcess2[0].GetProcessShape()->Integral() > 0.) fProcess2[0].GetProcessShape()->Scale(1./fProcess2[0].GetProcessShape()->Integral());
	      fProcess2[0].GetProcessShape()->Rebin(RebinFactor);
	    }
	    if(histnameh!="Template"){
	      histn2=histn-1;           
	      ///artificial shifting
	      /*
		
	      Float_t distIntegral=0;
	      
	      float binvalue =0;
	      
	      TH1F* Shifted = (TH1F*)hobj->Clone("Shifted");
	      
	      Shifted->Reset();
	      Int_t nbinz = Shifted->GetNbinsX();
	      for (int in= 1; in <= nbinz; in++) {
	      
	      binvalue =  0;
	      binvalue =hobj->GetBinContent(in);
	      //cout<<" bin "<<binvalue<<endl;
	      // if(in==0){
	      for(int i=1;i<10;i++){
	      Shifted->SetBinContent(i, 0.0);
	      }
	      if(binvalue>0.) Shifted->SetBinContent(in+3+histn2, binvalue);
	      
	      //}
	      
	      }
	      //TH1F* default_process10 = (TH1F*)fChannel[0].GetProcessProp(7).GetProcessShape()->Clone("default_process10");
	      //fProcess[10].GetProcessShape()->Reset();
	      //fProcess[ipr2].SetSystematicShape(histn2,Shifted);
	      //fProcess2[0].SetSystematicShape(histn2,Shifted);
	      
	      ///////////*/
	      fProcess2[0].SetSystematicShape(histn2,hobj);
	      if (fProcess2[0].GetSystematicShape(histn2)->Integral() > 0.) fProcess2[0].GetProcessShape()->Scale(1./fProcess2[0].GetSystematicShape(histn2)->Integral());
	      fProcess2[0].GetSystematicShape(histn2)->Rebin(RebinFactor);
	    }
	    histn=histn+1;
	  }
	  
	  if(histname1!="ZprimeTemplateDirectory"){
	    
	    if(histnameh=="Template"){
	      fProcess[ipr].SetProcessShape(hobj);
	      if (fProcess[ipr].GetProcessShape()->Integral() > 0.) fProcess[ipr].GetProcessShape()->Scale(1./fProcess[ipr].GetProcessShape()->Integral());
	      fProcess[ipr].GetProcessShape()->Rebin(RebinFactor);
	      ipr =ipr+1;     
	    }
	    if(histnameh!="Template"){
	      histn2=histn-1;
	      
	      ipr2 =ipr-1;
	      fProcess[ipr2].SetSystematicShape(histn2,hobj);
	      if (fProcess[ipr2].GetSystematicShape(histn2)->Integral() > 0.) fProcess[ipr2].GetProcessShape()->Scale(1./fProcess[ipr2].GetSystematicShape(histn2)->Integral());
	      fProcess[ipr2].GetSystematicShape(histn2)->Rebin(RebinFactor);
	      
	    }
	    
	    histn=histn+1;
	    
	  }
	  
	}
	
      }
    }
  }  
  
  SetProcessProp(fProcess, 4);
  SetSignalProp(fProcess2);
}

void TChProp::SetProcessProp(TPrProp* fprocess, int ipr0){
  
  for(int i=0;i<ipr0;i++){
    
    
    fProcessN[i].SetTotEff(fprocess[i].GetTotEff());
    fProcessN[i].SetPreScl(fprocess[i].GetPreScl());
    
    fProcessN[i].SetFilterEfficiency(fprocess[i].GetFilterEfficiency());
    fProcessN[i].SetacceptanceEfficiency(fprocess[i].GetacceptanceEfficiency());
    fProcessN[i].SetelectronIDEfficiency(fprocess[i].GetelectronIDEfficiency());
    fProcessN[i].SetmuonIDEfficiency(fprocess[i].GetmuonIDEfficiency());
    fProcessN[i].SettauIDEfficiency(fprocess[i].GettauIDEfficiency());
    fProcessN[i].SettopologyEfficiency(fprocess[i].GettopologyEfficiency());
    
    fProcessN[i].SetSystematics(fprocess[i].GetSystematics());
    
    fProcessN[i].SetbranchingFractionStatError(fprocess[i].GetbranchingFractionStatError());  
    fProcessN[i].SetacceptanceEfficiencyStatError(fprocess[i].GetacceptanceEfficiencyStatError());  
    fProcessN[i].SetelectronIDEfficiencyStatError(fprocess[i].GetelectronIDEfficiencyStatError());  
    fProcessN[i].SetmuonIDEfficiencyStatError(fprocess[i].GetmuonIDEfficiencyStatError());	
    fProcessN[i].SettauIDEfficiencyStatError(fprocess[i].GettauIDEfficiencyStatError());       
    fProcessN[i].SettopologyEfficiencyStatError(fprocess[i].GettopologyEfficiencyStatError());    
    fProcessN[i].SetbranchingFractionSystError(fprocess[i].GetbranchingFractionSystError());	 
    fProcessN[i].SetacceptanceEfficiencySystError(fprocess[i].GetacceptanceEfficiencySystError());  
    fProcessN[i].SetelectronIDEfficiencySystError(fprocess[i].GetelectronIDEfficiencySystError());  
    fProcessN[i].SetmuonIDEfficiencySystError(fprocess[i].GetmuonIDEfficiencySystError());	
    fProcessN[i].SettauIDEfficiencySystError(fprocess[i].GettauIDEfficiencySystError());       
    fProcessN[i].SettopologyEfficiencySystError(fprocess[i].GettopologyEfficiencySystError());    
    fProcessN[i].SetLuminosity(fprocess[i].GetLuminosity());
    fProcessN[i].SetSigma(fprocess[i].GetSigma());
    fProcessN[i].SetBR(fprocess[i].GetBR());
    
    TH1F* default_process = (TH1F*)fprocess[i].GetProcessShape()->Clone("default_process");
    fProcessN[i].SetProcessShape(default_process);
    int sizevect = fProcessN[i].GetSystematics()->size();
    for(int hists=0;hists<sizevect;hists++){
      TH1F* syst_processh = (TH1F*)fprocess[i].GetSystematicShape(hists)->Clone("syst_processh");
      fProcessN[i].SetSystematicShape(hists,syst_processh);
    }    
    
  }
}

void TChProp::SetSignalProp(TPrProp* fprocess3){
  //gROOT->ProcessLine("#include <vector>");
  
  
  fSignal[0].SetTotEff(fprocess3[0].GetTotEff());
  fSignal[0].SetPreScl(fprocess3[0].GetPreScl());
  
  fSignal[0].SetFilterEfficiency(fprocess3[0].GetFilterEfficiency());
  fSignal[0].SetacceptanceEfficiency(fprocess3[0].GetacceptanceEfficiency());
  fSignal[0].SetelectronIDEfficiency(fprocess3[0].GetelectronIDEfficiency());
  fSignal[0].SetmuonIDEfficiency(fprocess3[0].GetmuonIDEfficiency());
  fSignal[0].SettauIDEfficiency(fprocess3[0].GettauIDEfficiency());
  fSignal[0].SettopologyEfficiency(fprocess3[0].GettopologyEfficiency());
  
  fSignal[0].SetSystematics(fprocess3[0].GetSystematics());
  
  fSignal[0].SetbranchingFractionStatError(fprocess3[0].GetbranchingFractionStatError());  
  fSignal[0].SetacceptanceEfficiencyStatError(fprocess3[0].GetacceptanceEfficiencyStatError());  
  fSignal[0].SetelectronIDEfficiencyStatError(fprocess3[0].GetelectronIDEfficiencyStatError());  
  fSignal[0].SetmuonIDEfficiencyStatError(fprocess3[0].GetmuonIDEfficiencyStatError());      
  fSignal[0].SettauIDEfficiencyStatError(fprocess3[0].GettauIDEfficiencyStatError());	    
  fSignal[0].SettopologyEfficiencyStatError(fprocess3[0].GettopologyEfficiencyStatError());    
  fSignal[0].SetbranchingFractionSystError(fprocess3[0].GetbranchingFractionSystError());     
  fSignal[0].SetacceptanceEfficiencySystError(fprocess3[0].GetacceptanceEfficiencySystError());  
  fSignal[0].SetelectronIDEfficiencySystError(fprocess3[0].GetelectronIDEfficiencySystError());  
  fSignal[0].SetmuonIDEfficiencySystError(fprocess3[0].GetmuonIDEfficiencySystError());      
  fSignal[0].SettauIDEfficiencySystError(fprocess3[0].GettauIDEfficiencySystError());	    
  fSignal[0].SettopologyEfficiencySystError(fprocess3[0].GettopologyEfficiencySystError());    
  
  fSignal[0].SetLuminosity(fprocess3[0].GetLuminosity());
  fSignal[0].SetSigma(fprocess3[0].GetSigma());
  fSignal[0].SetBR(fprocess3[0].GetBR());
  
  TH1F* default_process = (TH1F*)fprocess3[0].GetProcessShape()->Clone("default_process");
  fSignal[0].SetProcessShape(default_process);
  int sizevectsig = fSignal[0].GetSystematics()->size();
  for(int hists=0;hists<sizevectsig;hists++){
    TH1F* syst_process = (TH1F*)fprocess3[0].GetSystematicShape(hists)->Clone("syst_process");
    fSignal[0].SetSystematicShape(hists,syst_process);
  }
}


TH1F* TChProp::Likelihood(int setting0, TChProp fChannelPro, int mcintg, TH1F* pseudo, TH1F* LogLVsSigma1){
  
  Float_t fitValue = 0; Float_t distValue_num=0; Float_t lkl=0; 
  
  
  double lumi = 50.0;
  
  TPrProp fProcessd;
  TSystematicStateGenerator fSystemacticState;
  TH1F* default_temp_bgtot    = (TH1F*)fChannelPro.GetProcessProp(0).GetProcessShape()->Clone("default_temp_totbg");
  TH1F* default_temp_signal   = (TH1F*)fChannelPro.GetSignalProp().GetProcessShape()  ->Clone("default_temp_signal");
  
  Int_t nBins = default_temp_signal->GetNbinsX();
  
  TH1F* pseudodata_0 = (TH1F*)pseudo->Clone("pseudodata_0");
  Int_t nbins = LogLVsSigma1->GetNbinsX();

  if(mcintg==0){
    default_temp_signal->Scale(fChannelPro.GetSignalProp().GetTotEff()*fChannelPro.GetChLumi()); 
    default_temp_bgtot->Scale(fChannelPro.GetProcessProp(0).GetTotEff()*fChannelPro.GetChLumi()*fChannelPro.GetProcessProp(0).GetSigma());
    for(int k=1; k<fChannelPro.GetNProcess(); k++){
      default_temp_bgtot->Add(fChannelPro.GetProcessProp(k).GetProcessShape(),fChannelPro.GetProcessProp(k).GetTotEff()*fChannelPro.GetChLumi()*fChannelPro.GetProcessProp(k).GetSigma()); 		
    }		    
    LogLVsSigma1->Reset();
    for (int in= 1; in <= nbins; in++) {
      Double_t currentSigma = LogLVsSigma1->GetBinCenter(in);// Min + (float(i)-0.5)*step;
      lkl=0;
      for (int i = 1; i <= nBins; i++){
        fitValue =(currentSigma*default_temp_signal->GetBinContent(i)+default_temp_bgtot->GetBinContent(i));
        if(fitValue <1e-50) continue;	      
        lkl = lkl + pseudodata_0->GetBinContent(i)*log(fitValue) - fitValue - TMath::LnGamma(pseudodata_0->GetBinContent(i) + 1);
      }
      LogLVsSigma1->SetBinContent(in, exp(lkl));
    }
    return(LogLVsSigma1);    
  }
  
  ////////////////////////////////////Integration Method///////////////////////////////////////////////
  Float_t Lumis = 0.0;
  if(mcintg>0){
    TH1F* LikelVsSigmaMCtot = (TH1F*)LogLVsSigma1->Clone();
    for( int mcin=0; mcin<mcintg; mcin++){
      
      default_temp_signal->Reset();
      default_temp_bgtot->Reset();
      
      Lumis = (Float_t)fSystemacticState.SmearLumi(fChannelPro.GetChLumi(),fChannelPro.GetChLumiErr());
      if(setting0==3){
	//cout<<"smearing lumi only"<<endl;
	for( int k=0; k<fChannelPro.GetNProcess(); k++){
	  // TH1F* Morph0_BG=(TH1F*)fSystemacticState.GetNormalizedDistribution(ichn, fChannelPro.GetProcessProp(k));
	  default_temp_bgtot->Add(fChannelPro.GetProcessProp(k).GetProcessShape(),fChannelPro.GetProcessProp(k).GetTotEff()*Lumis*fChannelPro.GetProcessProp(k).GetSigma());		
	}
	
	default_temp_signal->Reset();
	//TH1F* Morph0=(TH1F*)fSystemacticState.GetNormalizedDistribution(ichn, fChannelPro.GetSignalProp());
	//Normalizing Histogram to nominal product of effiencies
	default_temp_signal->Add(fChannelPro.GetSignalProp().GetProcessShape(),fChannelPro.GetSignalProp().GetTotEff()*Lumis); 
	
      }
      if(setting0<3){
	// cout<<"fun stuff"<<setting0<<endl;
	for( int k=0; k<fChannelPro.GetNProcess(); k++){
	  TH1F* Morph0_BG=(TH1F*)fSystemacticState.GetNormalizedDistribution(setting0, fChannelPro.GetProcessProp(k));
	  default_temp_bgtot->Add(Morph0_BG,fChannelPro.GetProcessProp(k).GetTotEff()*Lumis*fChannelPro.GetProcessProp(k).GetSigma());		
	}
	
	default_temp_signal->Reset();
	TH1F* Morph0=(TH1F*)fSystemacticState.GetNormalizedDistribution(setting0, fChannelPro.GetSignalProp());
	//Normalizing Histogram to nominal product of effiencies
	default_temp_signal->Add(Morph0,fChannelPro.GetSignalProp().GetTotEff()*Lumis);
	
      }
      
      LogLVsSigma1->Reset();
      for (int in= 1; in <= nbins; in++) {
	Double_t currentSigma =LogLVsSigma1 ->GetBinCenter(in);// Min + (float(i)-0.5)*step;
	lkl=0;
	for (int i = 1; i <= nBins; i++){
	  fitValue =(currentSigma*default_temp_signal->GetBinContent(i)+default_temp_bgtot->GetBinContent(i));
	  distValue_num = pseudodata_0->GetBinContent(i);
	  
	  if(fitValue <1e-50) continue;	  
	  lkl = lkl + distValue_num*log(fitValue) - fitValue - TMath::LnGamma(distValue_num + 1);
	}
	LogLVsSigma1->SetBinContent(in, exp(lkl));
      }
      LikelVsSigmaMCtot->Add(LogLVsSigma1);
    }
    double scl=0.0;
    scl = (double)(1.0/mcintg);
    LikelVsSigmaMCtot->Scale(scl);
    return(LikelVsSigmaMCtot);    
  } 
  
}

Double_t TChProp::Limit95CL(TH1F* LogLVsSigma){
  Double_t totalIntegral95 = 0.95*LogLVsSigma->Integral();
  Double_t currentIntegral = 0;
  int i = 1;
  
  while(currentIntegral < totalIntegral95){
    currentIntegral = LogLVsSigma->Integral(1, i);
    i++;
  }
  
  Double_t sigma95 = LogLVsSigma->GetBinLowEdge(i+1);
  return(sigma95);
  
}

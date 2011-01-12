/*********************************************************
********************************************************** 
***  October 10,2010					**
***  Indara Suarez					**
***  Texas A&M University				**
***	runs with TFitter.cc			    	**
*** 	calculates a systematic state vector   		**
*** 	performs morphing of mass templates		**
***	version 1.0   	   				**
**********************************************************
*********************************************************/
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <TStyle.h>
#include <TCanvas.h>
#include <cassert> 
#include <TFile.h>
#include <TTree.h>
#include <TKey.h> 
#include <TH1.h>

#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <THStack.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TRandom.h>
#include <TVector.h>
#include <TLorentzVector.h>
#include <TList.h>
#include <TChain.h>
#include <TMath.h>

using namespace std;
#include "TSystematicStateGenerator.h"
#include "TRandom3.h"

vector<float> TSystematicStateGenerator::GenerateSystematicState(int n){
	float alpha;
	vector<float> galphas2;
	gRandom->SetSeed(0);
	for(int i=0;i<n;i++){
	  alpha = (float)gRandom->Gaus(0,1.0);
	  // cout<<i<<"  "<<alpha<<endl;
  	  galphas2.push_back(alpha);

  	}
	return(galphas2);
}

Float_t TSystematicStateGenerator::SmearLumi(float lumi_mean, float lumi_err){
  float alphalumi=0;
  Float_t smearedlumi=0.;
  vector<float> sysvectorlumi = (vector<float>)GenerateSystematicState(1);
  alphalumi  = sysvectorlumi.at(0);
  smearedlumi = lumi_mean+(lumi_err*alphalumi);
  return(smearedlumi);
}	

TH1F* TSystematicStateGenerator::GetNormalizedDistribution(int setting, TPrProp ProcessProps){  
  float alph=0.;
  float err =0.;
  int o = 0;
  Float_t Norm = 1.;
  
  TH1F* fNormalizedZprime = (TH1F*)ProcessProps.GetProcessShape()->Clone("fNormalizedZprime");
  TH1F* Corrections = (TH1F*)ProcessProps.GetProcessShape()->Clone("Corrections");
  TH1F* Corrections1 = (TH1F*)ProcessProps.GetProcessShape()->Clone("Corrections1");
  Corrections->Reset();
  Corrections1->Reset();
  
  vector<float>* cummulativerrors =(vector<float>*)ProcessProps.GetSystematics();  
  int size = cummulativerrors->size();
  vector<float> sysvector = (vector<float>)GenerateSystematicState(size);
  vector<float>::iterator it;
  
  //Histograms
  for (it=sysvector.begin(); it<sysvector.end(); it++){
    if(setting != 0){
      //Morphing histogram
      TH1F* fNormalizedZprime1 = (TH1F*)ProcessProps.GetProcessShape()->Clone("fNormalizedZprime1");
      alph = 0.;
      alph = *it;
      Corrections1->Reset();
      Corrections1->Add(fNormalizedZprime1,-1);
      Corrections1->Add(ProcessProps.GetSystematicShape(o));
      Corrections->Add(Corrections1, alph); 
    } 
    //Normalization factor
    err = cummulativerrors->at(o);
    Norm = Norm+(alph*err);
    o++; 
  }
  
  if(setting!=0) fNormalizedZprime->Add(Corrections);
  //morphed Histogram normalized to 1
  if (fNormalizedZprime->Integral() > 0.){
    fNormalizedZprime->Scale(1./fNormalizedZprime->Integral());
  }
  
  //Normalizing morphed histogram to modified product of efficiency errors
  if(setting!=1) fNormalizedZprime->Scale(Norm);
  
  return(fNormalizedZprime);
  
}

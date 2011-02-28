/*********************************************************
********************************************************** 
***  April 23,2010					**
***  I. Suarez						**
***	ROOT > .L TFitter.cc  				**
*** 	ROOT > TFitter t      				**
*** 	ROOT > t.Loop()       				**
**********************************************************
*********************************************************/


//#define TFitter_cxx
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <TStyle.h>
#include <TCanvas.h>

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

#include "TFitter.h"
#include "TPrProp.h"
#include "TChProp.h"
#include "TSystematicStateGenerator.h"

using namespace std;

void TFitter::Loop(){


//Declaration of Variables
  
  TChProp fChannel[_theRootFiles.size()];
  int ichannel = 0;
  Double_t sigma95=0;  
  Double_t sigma950=0;
  Double_t joint_sigma95, joint_sigma95_mc = 0;
  
  TH1F* Limit_95CL_0 =new TH1F("Limit_95CL_0","Limit_95CL_0", 60,0,20);
  TH1F* Limit_95CL_1 =new TH1F("Limit_95CL_1","Limit_95CL_1", 60,0,20);
  TH1F* Joint_Limit_95CL =new TH1F("Joint_Limit_95CL","Joint_Limit_95CL", 60,0,20);
  TH1F* Joint_Limit_95CL_MC =new TH1F("Joint_Limit_95CL_MC","Joint_Limit_95CL_MC", 60,0,20);
  
 //***********************Settings*****************************//

  Int_t npseudo = _theNExp;
  Int_t systematicset = _SystSettings;
  Int_t rebinfactor = _theRebinFactor;
    
  Float_t lumi_error = _theLumiErr;
  Float_t lumi = _theLumi; //in ipb
  Double_t mass = _theMass; //z-prime mass, 
  Double_t sigma_input = _theXSection; //ipb
  string zprimetemplate = _theZprimetemp;  //name of zprime template directory you wish to analyze
  string zprimetree = _theZprimetree;  //name of zprime tree you wish to analyze
  
  int mcint = _theNMCInt;
  
//***********************Creating Objects*****************************//   
  for(unsigned int vIt = 0; vIt < _theRootFiles.size(); vIt++){
     TFile* fchntuple = new TFile(_theRootFiles.at(vIt).c_str(),"read");
     fChannel[vIt].SetPrProp(ichannel, fchntuple, zprimetree, zprimetemplate, rebinfactor);
     fChannel[vIt].SetChLumi(lumi);
     fChannel[vIt].SetChLumiErr(lumi_error);
  }
  //Test if your objects were filled correctly
 /* for(int l=0;l<1;l++){
    cout<<"Channel: "<<l<<"  BR for signal : "<<fChannel[l].GetSignalProp().GetmuonIDEfficiency()<<endl;
    cout<<"Channel: "<<l<<"  XSEC for signal : "<<fChannel[l].GetSignalProp().GetSigma()<<endl;
    cout<<"Channel: "<<l<<"  N Process : "<<(int)fChannel[l].GetNProcess()<<endl;
    cout<<"Channel: "<<l<<"  Shape for signal : "<<fChannel[l].GetSignalProp().GetProcessShape()->Integral()<<endl;
    }*/
  
  

  //********************************Calculations*************************//
  int nBins2 = 1400;									 						
  double Min = 0;									 						
  double Max = 35;
  
  vector<float> sysvectorlumi;
  
  vector< vector<float> > systematicvector;
  
  TH1F* pseudodata_MassDistribution_0;
  TH1F* LogLVsSigma1 = new TH1F("LogLVsSigma1","LogLVsSigma1", nBins2, Min, Max);	 						
  TH1F* JointlklvSigma= new TH1F("JointlklvSigma","JointlklvSigma", nBins2, Min, Max);   
  TH1F* JointlklvSigma_MC= new TH1F("JointlklvSigma_MC","JointlklvSigma_MC", nBins2, Min, Max);
  
  TSystematicStateGenerator fSystemacticState_0;
  int vectorsize = ( (vector<float>*)fChannel[0].GetSignalProp().GetSystematics() )->size();
  
  if(npseudo==0){
    JointlklvSigma->Reset();
    
    /////////////////////generating luminosity smearing vector////////////////////////////////////////////
    if(mcint>0){
      sysvectorlumi = (vector<float>)fSystemacticState_0.GenerateSystematicState(mcint);
      for(int i=0;i<mcint;i++){
        systematicvector.push_back((vector<float>)fSystemacticState_0.GenerateSystematicState(vectorsize));
      }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    for(int i=0;i<mcint;i++){
        for(int j=0;j<vectorsize;j++){
      //  cout<<"mcint "<<i<<" systematic "<<j<<" has nuisance paramete = "<<( (vector<float>)systematicvector.at(i) ).at(j)<<endl;
      }
    }
	
	
    
    for(unsigned int ichannel = 0; ichannel < _theRootFiles.size(); ichannel++){
      TH1F* data_MassDistribution = (TH1F*)fChannel[ichannel].GetDataHistogram()->Clone("data_MassDistribution");
      cout<<"DATA INTEGRAL =  "<<data_MassDistribution->Integral()<<endl;
      //calculating likehood distributions for case with no systematics  											  
      LogLVsSigma1->Reset();	       
      TH1F* LogLVsSigma0 =(TH1F*)fChannel[ichannel].Likelihood(systematicset, sysvectorlumi,  systematicvector, fChannel[ichannel],0,data_MassDistribution, LogLVsSigma1)->Clone("LogLVsSigma0");
      sigma950 = fChannel[ichannel].Limit95CL(LogLVsSigma0);  //calculating 95%CL Limit
      Limit_95CL_0->Fill(sigma950); //This is the 95%CL Limit distribution without systematics
      cout<<"For Channel "<<ichannel<<"  Limit with no systematics for data "<<sigma950<<endl;
      //Filling output histogram
      if(mcint ==0){	
        _theNameHistoMap[ichannel]->Fill(sigma950);  
      }
      //calculating likelihood distribution with systematics
      TH1F* LogLVsSigmaMC0;
      if(mcint>0){
	LogLVsSigma1->Reset();
	LogLVsSigmaMC0=(TH1F*)fChannel[ichannel].Likelihood(systematicset, sysvectorlumi, systematicvector, fChannel[ichannel], mcint,data_MassDistribution, LogLVsSigma1)->Clone("LogLVsSigmaMC0");
	sigma95 = fChannel[ichannel].Limit95CL(LogLVsSigmaMC0);
	cout<<" 95%CL Limit for channel "<<ichannel<<" for "<<mcint<<" mc itterations fitting data is = "<<sigma95<<endl;
	Limit_95CL_1->Fill(sigma95); //This is the 95%CL Limit distribution with systematics
	//_theNameHistoMap[ichannel]->Fill(sigma95);
      }
      
      //joint likelihood
      if(ichannel ==0){
	JointlklvSigma->Reset();
	if(mcint>0) JointlklvSigma_MC->Reset();
	JointlklvSigma->Add(LogLVsSigma0);
	if(mcint>0) JointlklvSigma_MC->Add(LogLVsSigmaMC0);
      }
      
      if(ichannel>0){
	JointlklvSigma->Multiply(LogLVsSigma0);
	if(mcint>0) JointlklvSigma_MC->Multiply(LogLVsSigmaMC0);  //joint likelihood
      }
      
      pseudodata_MassDistribution_0 = (TH1F*)data_MassDistribution->Clone("pseudodata_MassDistribution_0");
      
    }
    //joint likelihood
    joint_sigma95 = fChannel[0].Limit95CL(JointlklvSigma);
    Joint_Limit_95CL->Fill(joint_sigma95);
    
    cout<<" Joint Limit without systematics "<<joint_sigma95<<endl;
    
    if(mcint>0){
    joint_sigma95_mc = fChannel[0].Limit95CL(JointlklvSigma_MC);
    Joint_Limit_95CL_MC->Fill(joint_sigma95_mc);   
    cout<<" Joint Limit with systematics "<<joint_sigma95_mc<<endl;
    }
    
  
  
  }
  //if want to run pseudoexperiments
  if(npseudo>0){
  
  for(int lim=0;lim<npseudo;lim++){	
  //cout<<"Pseudoexperiment  "<<lim<<"   ";
    //generating luminosity smearing vector
    if(mcint>0)   sysvectorlumi = (vector<float>)fSystemacticState_0.GenerateSystematicState(mcint);
    
    //													    

    JointlklvSigma->Reset();		
    
    for(unsigned int ichannel = 0; ichannel < _theRootFiles.size(); ichannel++){
    
  //  cout<<"Channel  "<<ichannel<<endl;
      
      //Generating Pseudodata
      TH1F* default_temp_bgtots = (TH1F*)fChannel[ichannel].GetProcessProp(0).GetProcessShape()->Clone("default_temp_bgtots");
      TH1F* pseudodata_MassDistribution = (TH1F*)default_temp_bgtots->Clone("pseudodata_MassDistribution");
      default_temp_bgtots->Reset();
      for( int k=0; k<fChannel[ichannel].GetNProcess(); k++){
    	default_temp_bgtots->Add(fChannel[ichannel].GetProcessProp(k).GetProcessShape(),fChannel[ichannel].GetProcessProp(k).GetTotEff()*fChannel[ichannel].GetChLumi()*fChannel[ichannel].GetProcessProp(k).GetSigma());
      }
      gRandom->SetSeed(0);
      int distIntegral = int(gRandom->Poisson(default_temp_bgtots->Integral()));
      pseudodata_MassDistribution->Reset();
      pseudodata_MassDistribution->FillRandom(default_temp_bgtots, distIntegral);
      
      //calculating likehood distributions for case with no systematics  											  
      LogLVsSigma1->Reset();	       
      TH1F* LogLVsSigma0 =(TH1F*)fChannel[ichannel].Likelihood(systematicset, sysvectorlumi, systematicvector, fChannel[ichannel],0,pseudodata_MassDistribution, LogLVsSigma1)->Clone("LogLVsSigma0");
      sigma950 = fChannel[ichannel].Limit95CL(LogLVsSigma0);  //calculating 95%CL Limit
      Limit_95CL_0->Fill(sigma950); //This is the 95%CL Limit distribution without systematics
      //Filling output histogram
      if(mcint ==0){	
        _theNameHistoMap[ichannel]->Fill(sigma950);  
      }
      //calculating likelihood distribution with systematics
      TH1F* LogLVsSigmaMC0;
      if(mcint>0){
	LogLVsSigma1->Reset();
	
	LogLVsSigmaMC0=(TH1F*)fChannel[ichannel].Likelihood(systematicset, sysvectorlumi, systematicvector, fChannel[ichannel], mcint,pseudodata_MassDistribution, LogLVsSigma1)->Clone("LogLVsSigmaMC0");

	sigma95 = fChannel[ichannel].Limit95CL(LogLVsSigmaMC0);
	Limit_95CL_1->Fill(sigma95); //This is the 95%CL Limit distribution with systematics
	
	//_theNameHistoMap[ichannel]->Fill(sigma95);
      }
      
      //joint likelihood
      if(ichannel ==0){
	JointlklvSigma->Reset();
	if(mcint>0) JointlklvSigma_MC->Reset();
	JointlklvSigma->Add(LogLVsSigma0);
	if(mcint>0) JointlklvSigma_MC->Add(LogLVsSigmaMC0);
      }
      
      if(ichannel>0){
	JointlklvSigma->Multiply(LogLVsSigma0);
	if(mcint>0) JointlklvSigma_MC->Multiply(LogLVsSigmaMC0);  //joint likelihood
      }
      
      if(lim==npseudo-1) pseudodata_MassDistribution_0 = (TH1F*)pseudodata_MassDistribution->Clone("pseudodata_MassDistribution_0");
      
    }
    //joint likelihood
    joint_sigma95 = fChannel[0].Limit95CL(JointlklvSigma);
    Joint_Limit_95CL->Fill(joint_sigma95);
    
    if(mcint>0){
    joint_sigma95_mc = fChannel[0].Limit95CL(JointlklvSigma_MC);
    Joint_Limit_95CL_MC->Fill(joint_sigma95_mc);   
    }
    }
    
  }
  
  ////////////OUTPUT HISTOGRAMS//////////////////// 
  
  //The 95%CL Limit with and without systematics output as a root file
  TCanvas *theLimitCanvas = new TCanvas("theLimitCanvas","theLimitCanvas");
  theLimitCanvas->cd();
  Limit_95CL_0->SetTitle("muTau 95%CL Limit; 95%CL Limit on #sigma_{(pp->Z')} (pb); ");
  Limit_95CL_1->SetTitle("Limit w/ systemtics");
  Limit_95CL_0->SetFillColor(30);
  Limit_95CL_1->SetLineWidth(2);
  Limit_95CL_0->Draw("hist");
  Limit_95CL_1->Draw("histsames");
  theLimitCanvas->BuildLegend();
  theLimitCanvas->Update();
  theLimitCanvas->Print("SensitivityLim95CL.root");
  
//This is the band plot
  TGraphAsymmErrors* xsec_band = new TGraphAsymmErrors();
  TGraphAsymmErrors* xsec = new TGraphAsymmErrors();
  TGraphAsymmErrors* txsec= new TGraphAsymmErrors();
  TCanvas *theCLCanvas = new TCanvas("theCLCanvas","theCLCanvas");				 						
  theCLCanvas->Clear();  									 						
  theCLCanvas->Divide(1,1);									 						
  for(unsigned int vIt = 0; vIt < _theRootFiles.size(); vIt++){
    
    theCLCanvas->cd(vIt + 1);
    std::cout << "No Systematics"<<endl;
    cout<<"Limit = \t" << (Limit_95CL_0->GetMean())
	<< "  Normalization  "<<fChannel[vIt].GetSignalProp().GetTotEff()*fChannel[vIt].GetChLumi() <<"  \t +- \t" << Limit_95CL_0->GetRMS() << std::endl;
    
    std::cout << "** Systematics **"<<endl;
    cout<<"Morphing and/or Smearing due to systematics defined in ntuple, and Smearing of Luminosity **"<<endl;
    cout<<"Limit = \t" << (Limit_95CL_1->GetMean())
	<<"  \t +- \t" << Limit_95CL_1->GetRMS() << std::endl;      

    
    fChannel[vIt].SetLimit95CLH(_theNameHistoMap[vIt]);
    fChannel[vIt].GetLimit95CLH()->SetTitle(" 95%CL Limit");
    fChannel[vIt].GetLimit95CLH()->GetXaxis()->SetTitle("Sigma");
    fChannel[vIt].GetLimit95CLH()->DrawCopy();

    xsec_band->SetPoint(vIt, vIt, _theNameHistoMap[vIt]->GetMean());
    xsec_band->SetPointError(vIt, 0., 0.,_theNameHistoMap[vIt]->GetRMS(),_theNameHistoMap[vIt]->GetRMS());
    xsec->SetPoint(vIt, vIt, _theNameHistoMap[vIt]->GetMean());
    txsec->SetPoint(vIt, vIt, 1.914);
  }
  cout<<"******************** JOINT LIMIT = \t" << (Joint_Limit_95CL_MC->GetMean())<<"  *******************************"<<endl;
  
  theCLCanvas->Print("95PerCentCL.eps");

  TCanvas *theCLBandCanvas = new TCanvas("theCLBandCanvas","theCLBandCanvas");                                                         
  theCLBandCanvas->Clear();    
    theCLBandCanvas->cd();
    xsec->SetLineColor(2);
    txsec->SetLineColor(4);
    xsec_band->SetFillColor(kCyan);
    xsec_band->GetYaxis()->SetLabelSize(0.04);
    xsec_band->GetXaxis()->SetLabelSize(0.04);
    xsec_band->GetXaxis()->SetTitleSize(0.05);
    xsec_band->GetYaxis()->SetTitleSize(0.05);
    xsec_band->GetYaxis()->SetTitleOffset(1.25);
    xsec_band->GetXaxis()->SetTitleOffset(1.25);
    xsec_band->GetXaxis()->SetTitle("L_{int}(pb^{-1})");
    xsec_band->GetYaxis()->SetTitle("#sigma_{95%}");
    xsec_band->GetYaxis()->SetRangeUser(0,10);
    xsec_band->Draw("Ae3");
    xsec->Draw("same");
    txsec->Draw("same");

    TLatex* tex1 = new TLatex(30.,4.5,"CMS Preliminary 2010");
    TLatex* tex2 = new TLatex(30.,4.2,"Z' #rightarrow #tau#tau");
    tex1->SetTextSize(0.03);
    tex2->SetTextSize(0.03);
    tex1->SetLineWidth(2); tex2->SetLineWidth(2);
    tex1->Draw(); tex2->Draw();

    TLegend *legend = new TLegend(0.65,0.77,0.92,0.88,NULL,"brNDC");
    legend->AddEntry(xsec,"Experiment", "l");
    legend->AddEntry(txsec,"Theory", "l");
    legend->SetTextFont(42);
    legend->SetLineColor(1);
    legend->SetLineStyle(1);
    legend->SetLineWidth(1);
    legend->SetFillColor(0);
    legend->SetFillStyle(1001);
    legend->SetBorderSize(0);
    legend->SetFillColor(kWhite);
    legend->Draw();

  theCLBandCanvas->Print("95PerCentCLBand.eps");

  TCanvas* theJointCLCanvas = new TCanvas("theJointCLCanvas","theJointCLCanvas");
  theJointCLCanvas->cd();
  Joint_Limit_95CL->Draw();
  Joint_Limit_95CL_MC->Draw("same");
  theJointCLCanvas->Print("Joint95PerCentCL.eps");
  
  //For Stacked Mass Histogram, make sure you modify the names and number of processes for your channel; this one is an example for mutau
  TCanvas* StackedMass = new TCanvas("StackedMass","StackedMass");
  StackedMass->cd();
  double default_background_scale = 1.; //0.000001;
  double ymin = 0.001;
  TH1F* default_temp_signal_0=  (TH1F*)fChannel[0].GetSignalProp().GetProcessShape()->Clone("default_temp_signal_0");
  THStack *hs = new THStack("hs","MuTau Channel");
  for( int k=0; k<fChannel[0].GetNProcess(); k++){
    
    if(k==0){
      TH1F *Ztautau = (TH1F*)fChannel[0].GetProcessProp(k).GetProcessShape()->Clone("Ztautau");
      
      Ztautau->Scale(fChannel[0].GetProcessProp(k).GetTotEff()*fChannel[0].GetChLumi()*fChannel[0].GetProcessProp(k).GetSigma()*default_background_scale);
      cout<<"ZTAUTAU integral  = "<<Ztautau->Integral()<<endl;
     // cout<<"ZTAUTAU  "<<fChannel[0].GetProcessProp(k).GetTotEff()*fChannel[0].GetChLumi()*fChannel[0].GetProcessProp(k).GetSigma()*default_background_scale<<endl;
      Ztautau->SetFillColor(kBlue);
      Ztautau->SetTitle("Ztautau");
      //c2->cd();
      //Ztautau->DrawCopy();
      hs->Add(Ztautau);
    }
    if(k==2){
      TH1F *TTBar =  (TH1F*)fChannel[0].GetProcessProp(k).GetProcessShape()->Clone("TTBar");
      TTBar->Scale(fChannel[0].GetProcessProp(k).GetTotEff()*fChannel[0].GetChLumi()*fChannel[0].GetProcessProp(k).GetSigma()*default_background_scale);
      cout<<"TTBar integral  = "<<TTBar->Integral()<<endl;
      TTBar->SetFillColor(kRed-9);
      TTBar->SetTitle("TTBar");
      hs->Add(TTBar);
    }
    if(k==3){
      TH1F *WJets =  (TH1F*)fChannel[0].GetProcessProp(k).GetProcessShape()->Clone("WJets");
      WJets->Scale(fChannel[0].GetProcessProp(k).GetTotEff()*fChannel[0].GetChLumi()*fChannel[0].GetProcessProp(k).GetSigma()*default_background_scale);
      cout<<"WJets integral  = "<<WJets->Integral()<<endl;
      WJets->SetFillColor(kCyan);
      WJets->SetTitle("WJets");
      hs->Add(WJets);
    }
    
    if(k==1){
      TH1F *Zmumu =  (TH1F*)fChannel[0].GetProcessProp(k).GetProcessShape()->Clone("Zmumu");
      Zmumu->Scale(fChannel[0].GetProcessProp(k).GetTotEff()*fChannel[0].GetChLumi()*fChannel[0].GetProcessProp(k).GetSigma()*default_background_scale);
    //   cout<<"ZMuMu  "<<fChannel[0].GetProcessProp(k).GetTotEff()<<endl;
      cout<<"Zmumu integral  = "<<Zmumu->Integral()<<endl;
      Zmumu->SetFillColor(kOrange);
      Zmumu->SetTitle("Zmumu");
      hs->Add(Zmumu);   
      }   
      if(k==4){
      TH1F *IncMu15 =  (TH1F*)fChannel[0].GetProcessProp(k).GetProcessShape()->Clone("IncMu15");
      IncMu15->Scale(fChannel[0].GetProcessProp(k).GetTotEff()*fChannel[0].GetChLumi()*fChannel[0].GetProcessProp(k).GetSigma()*default_background_scale);
      cout<<"IncMu15 integral  = "<<IncMu15->Integral()<<endl;
      IncMu15->SetFillColor(kYellow);
      IncMu15->SetTitle("IncMu15");
      hs->Add(IncMu15);
    
      
      default_temp_signal_0->Scale(fChannel[0].GetSignalProp().GetTotEff()*fChannel[0].GetChLumi()*fChannel[0].GetSignalProp().GetSigma());
      cout<<"MCSIGNAL  "<<fChannel[0].GetSignalProp().GetTotEff()<<endl;
      cout<<"MCSignal integral  = "<<default_temp_signal_0->Integral()<<endl;;
      default_temp_signal_0->SetTitle("Signal, Zprime");
      
      // default_temp_signal_0->SetFillColor(kRed-7);
      //  default_temp_signal_0->SetFillStyle(3325);
      //hs->Add(default_temp_signal_0);     
    }
  } 
  
  for( int k=0; k<fChannel[1].GetNProcess(); k++){
   (fChannel[1].GetProcessProp(k).GetProcessShape())->Scale(fChannel[1].GetProcessProp(k).GetTotEff()*fChannel[1].GetChLumi()*fChannel[1].GetProcessProp(k).GetSigma()*default_background_scale);
   cout<<"  Process "<<k<<" has "<<(fChannel[1].GetProcessProp(k).GetProcessShape())->Integral()<<" events "<<endl;
    }    
  
  
  StackedMass->cd();
  StackedMass->SetLogy();
  
  hs->SetTitle("BG Distributions & Pseudodata for muTau Channel; muTau Mass (GeV);Events per 50GeV Bins at 35ipb");
  hs->Draw("hist");
  hs->SetMinimum(ymin);
  hs->SetMaximum(10.);
  pseudodata_MassDistribution_0->SetLineColor(kRed+3);
  pseudodata_MassDistribution_0->SetLineWidth(2);
  pseudodata_MassDistribution_0->SetMarkerColor(kRed+3);
  pseudodata_MassDistribution_0->SetMarkerStyle(0);
  cout<<"Data Integral / pseudodata example =  "<<pseudodata_MassDistribution_0->Integral()<<endl;
  //pseudodata_MassDistribution_0->SetMarkerSize(2);
  pseudodata_MassDistribution_0->SetTitle("Data");
  pseudodata_MassDistribution_0->Draw("same");
  default_temp_signal_0->Draw("samehist");
  default_temp_signal_0->SetLineWidth(5);
 /* TSystematicStateGenerator fSystemacticState_0;
  TH1F* Morphing=(TH1F*)fSystemacticState_0.GetNormalizedDistribution(0, fChannel[0].GetSignalProp());
  Morphing->Scale(fChannel[0].GetSignalProp().GetTotEff()*fChannel[0].GetChLumi()*fChannel[0].GetSignalProp().GetSigma());
  Morphing->SetLineWidth(2);
  Morphing->SetLineColor(8);
  Morphing->SetTitle("Morphed Distribution");
  Morphing->Draw("same");*/
  StackedMass->BuildLegend();
  StackedMass->Update();      
  
  //hs->GetXaxis()->SetTitle("MuTau Mass (GeV)");
  // hs->GetYaxis()->SetLimits(1., 5.);
  
  
  StackedMass->Print("StackedMass.root"); 
  
  //These are the systematic histograms used for morphing of signal shape
  TCanvas* ShapeSyst = new TCanvas("ShapeSyst","ShapeSyst");
  ShapeSyst->cd();
  default_temp_signal_0->SetLineWidth(3);
  default_temp_signal_0->Draw();
  ShapeSyst->BuildLegend();
  ShapeSyst->Update();
  for(int kp=0;kp<6;kp++){
    fChannel[0].GetSignalProp().GetSystematicShape(kp)->SetLineColor(kBlue+2*kp);
    fChannel[0].GetSignalProp().GetSystematicShape(kp)->Scale(fChannel[0].GetSignalProp().GetTotEff()*fChannel[0].GetChLumi()*fChannel[0].GetSignalProp().GetSigma());
    fChannel[0].GetSignalProp().GetSystematicShape(kp)->Draw("same");
  }
  ShapeSyst->Print("ShapeSystematics.eps"); 
  
}

void TFitter::SetNExp(int theNExp){
  _theNExp = theNExp;
}

void TFitter::SetRebinFactor(int theRebinFactor){
  _theRebinFactor = theRebinFactor;
}

void TFitter::SetLumi(float theLumi){
  _theLumi = theLumi;
}

void TFitter::SetLumiErr(float theLumiErr){
  _theLumiErr = theLumiErr;
}

void TFitter::SetSystSettings(int SystSettings){
  _SystSettings = SystSettings;
}

void TFitter::SetZprimeTempName(string zprimetemp){
  _theZprimetemp = zprimetemp;
}

void TFitter::SetZprimeTreeName(string zprimetree){
  _theZprimetree = zprimetree;
}

void TFitter::NMCInt(int theMCInt){
  _theNMCInt = theMCInt;
}

void TFitter::SetSignalMassXsection(float theMass, float theXSection){
  _theMass = theMass;
  _theXSection = theXSection;
}

void TFitter::GetTheFiles(){
  char line[120],FileName[120];	 
  FILE *fTreeFiles = fopen("channels.dat","read");
  while ( fgets( line, 120, fTreeFiles ) ) {
    sscanf(line,"%s",FileName);
    _theRootFiles.push_back(FileName);
  }
  fclose(fTreeFiles);
}

void TFitter::BookCLHistos(){
  for(unsigned int vIt = 0; vIt < _theRootFiles.size(); vIt++){
    char* theHistoName = Form("histo95CL%d",vIt);
    _theNameHistoMap[vIt] = new TH1F(theHistoName,theHistoName,150,0,15);
  }
}

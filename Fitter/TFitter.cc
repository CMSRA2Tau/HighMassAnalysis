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
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TKey.h" 
#include "TH1.h"

#include "THStack.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom.h"
#include "TVector.h"
#include "TLorentzVector.h" 
#include "TList.h" 
#include "TChain.h"
#include "TMath.h"

#include "TFitter.h"
#include "TPrProp.h"
#include "TPrProp.cc"
#include "TChProp.h"
#include "TChProp.cc"

using namespace std;

void TFitter::Loop(){


	Float_t         fbranchingFraction;
   	Float_t         fFilterEfficiency;
   	Float_t         fxsection;
   	Float_t         facceptanceEfficiency;
   	Float_t         felectronIDEfficiency;
   	Float_t         fmuonIDEfficiency;
   	Float_t         ftauIDEfficiency;
   	Float_t         ftopologyEfficiency;
   	Float_t         fbranchingFractionStatError;
   	Float_t         facceptanceEfficiencyStatError;
   	Float_t         felectronIDEfficiencyStatError;
   	Float_t         fmuonIDEfficiencyStatError;
   	Float_t         ftauIDEfficiencyStatError;
   	Float_t         ftopologyEfficiencyStatError;
   	Float_t         fbranchingFractionSystError;
   	Float_t         facceptanceEfficiencySystError;
   	Float_t         felectronIDEfficiencySystError;
   	Float_t         fmuonIDEfficiencySystError;
   	Float_t         ftauIDEfficiencySystError;
   	Float_t         ftopologyEfficiencySystError;

//Declaration of Variables
  TPrProp fProcess[15];
  TPrProp fProcess2[1];
  TChProp fChannel[4];
  TFile *fFileRead;
  TKey *key;
  TKey *keyh;
  TList *h2;
  TList *hl;
  int nb=0;
  int ch=0;
  //int ch=0;
  int ipr=0;
  int num=0;
  int syst_vector=0;
  Float_t PreScale[4][10];
  Float_t TotEff[4][10];
  Double_t distValue_num;
  Double_t fitValue;
  Float_t prescl=0;
  Float_t toteff=0;
  Double_t lkl=0;
  Float_t totscl1=0;
  Float_t totscl0=0;
  int ichannel = 0;
  
  Float_t JointLikelihood = 0;

  Double_t sigma95=0;
  Double_t sigma951=0;
  Double_t sigma952=0;
  Double_t sigma953=0;
  
  Double_t joint_sigma95 = 0;
  TH1F* Limit_95CL_0 =new TH1F("Limit_95CL_0","Limit_95CL_0", 150,0,15);
  TH1F* Limit_95CL_1 =new TH1F("Limit_95CL_1","Limit_95CL_1", 150,0,15);
  TH1F* Limit_95CL_2 =new TH1F("Limit_95CL_2","Limit_95CL_2", 150,0,30);
  TH1F* Limit_95CL_01 =new TH1F("Limit_95CL_01","Limit_95CL_0", 6000,0,15);
  TH1F* Limit_95CL_11 =new TH1F("Limit_95CL_11","Limit_95CL_1", 6000,0,15);
  TH1F* Limit_95CL_21 =new TH1F("Limit_95CL_21","Limit_95CL_2", 6000,0,15);
  TH1F* Limit_95CL_02 =new TH1F("Limit_95CL_02","Limit_95CL_0", 6000,0,15);
  TH1F* Limit_95CL_12 =new TH1F("Limit_95CL_12","Limit_95CL_1", 6000,0,15);
  TH1F* Limit_95CL_22 =new TH1F("Limit_95CL_22","Limit_95CL_2", 6000,0,15);
  TH1F* Limit_95CL_03 =new TH1F("Limit_95CL_03","Limit_95CL_0", 6000,0,15);
  TH1F* Limit_95CL_13 =new TH1F("Limit_95CL_13","Limit_95CL_1", 6000,0,15);
  TH1F* Limit_95CL_23 =new TH1F("Limit_95CL_23","Limit_95CL_2", 6000,0,15);
  TH1F* Joint_Limit_95CL =new TH1F("Joint_Limit_95CL","Joint_Limit_95CL", 150,0,15);
  
  
    /*TCanvas *cn3 = new TCanvas("cn3","", 100, 100, 1000, 1000);
  cn3->Clear();
  cn3->Divide(1,2); 
  
  TCanvas *cn2 = new TCanvas("cn2","", 100, 100, 1000, 1000);
  cn2->Clear();
  cn2->Divide(1,2);
  

  TCanvas *cn1 = new TCanvas("cn1","", 100, 100, 1000, 1000);
  cn1->Clear();
  cn1->Divide(1,2);
  */
  
  
  
 //***********************Settings*****************************//

  Int_t npseudo=300;
    
  double lumi_error = 1;
  Double_t lumi = 50.0; //in pb
  Double_t mass=500.; //z-prime mass
  Double_t sigma_input = 1.94;

  int input_sigma_signal = 0;
  double RebinFactor =1;
  Double_t background_scale = 1;
  Double_t signal_scale =1;
  //0.000001;
  //99999+20999+2+942; //for cross checks
  // 0.0000000000000001; //for cross checks
   //more debugging
   
   //mc integrations
   
  int mcint = 0;
//***********************Settings*****************************//   

char line[120],FileName[120];	 
FILE *fTreeFiles = fopen("channels.dat","read");

	while ( fgets( line, 120, fTreeFiles ) ) {
	
	nb=0;
	ch=0;
	ipr=0;
	num=0;
 	prescl=0;
	toteff=0;
 	lkl=0;
	totscl1=0;
	 totscl0=0;

    	sscanf(line,"%s",FileName);

   	TFile* fntuple = new TFile(FileName,"read");
   	if ( fntuple->IsOpen() ) printf("File opened successfully\n"); 
	//cout<<FileName<<endl;
   
   	TDirectory *dirsav = gDirectory;

   	TIter next(dirsav->GetListOfKeys());
   	while ((key = (TKey*)next())) {
      		string histname1 = key->GetName();
		//cout<<"title1  "<<histname1<<"   ";
      
      		TObject *obj = key->ReadObj();
        	string histname2 = obj->ClassName();
		//cout<<"type  "<<histname2<<endl;
		
		if ( obj->IsA()->InheritsFrom( "TTree" ) ) {
			string teffs = obj->GetName();
			TTree *tobj = (TTree*)obj;
				
			Int_t Nentries = (Int_t)tobj->GetEntries();	
			tobj->SetBranchAddress("branchingFraction", 		&fbranchingFraction		       );   //&b_fbranchingFraction);
   			tobj->SetBranchAddress("FilterEfficiency", 		&fFilterEfficiency		       );   //&b_fFilterEfficiency);
   			tobj->SetBranchAddress("xsection", 			&fxsection			       );   //&b_fxsection);
   			tobj->SetBranchAddress("acceptanceEfficiency", 		&facceptanceEfficiency		       );   //&b_facceptanceEfficiency);
   			tobj->SetBranchAddress("electronIDEfficiency", 		&felectronIDEfficiency		       );   //&b_felectronIDEfficiency);
   			tobj->SetBranchAddress("muonIDEfficiency", 		&fmuonIDEfficiency		       );   //&b_fmuonIDEfficiency);
   			tobj->SetBranchAddress("tauIDEfficiency", 		&ftauIDEfficiency		       );   //&b_ftauIDEfficiency);
   			tobj->SetBranchAddress("topologyEfficiency", 		&ftopologyEfficiency		       );   //&b_ftopologyEfficiency);
   			tobj->SetBranchAddress("branchingFractionStatError", 	&fbranchingFractionStatError	       );   //&b_fbranchingFractionStatError);
   			tobj->SetBranchAddress("acceptanceEfficiencyStatError", &facceptanceEfficiencyStatError         );   //&b_facceptanceEfficiencyStatError);
   			tobj->SetBranchAddress("electronIDEfficiencyStatError", &felectronIDEfficiencyStatError         );   //&b_felectronIDEfficiencyStatError);
   			tobj->SetBranchAddress("muonIDEfficiencyStatError", 	&fmuonIDEfficiencyStatError	       );   //&b_fmuonIDEfficiencyStatError);
   			tobj->SetBranchAddress("tauIDEfficiencyStatError", 	&ftauIDEfficiencyStatError	       );   //&b_ftauIDEfficiencyStatError);
   			tobj->SetBranchAddress("topologyEfficiencyStatError", 	&ftopologyEfficiencyStatError	       );   //&b_ftopologyEfficiencyStatError);
   			tobj->SetBranchAddress("branchingFractionSystError", 	&fbranchingFractionSystError	       );   //&b_fbranchingFractionSystError);
   			tobj->SetBranchAddress("acceptanceEfficiencySystError", &facceptanceEfficiencySystError         );   //&b_facceptanceEfficiencySystError);
   			tobj->SetBranchAddress("electronIDEfficiencySystError", &felectronIDEfficiencySystError         );   //&b_felectronIDEfficiencySystError);
   			tobj->SetBranchAddress("muonIDEfficiencySystError", 	&fmuonIDEfficiencySystError	       );   //&b_fmuonIDEfficiencySystError);
   			tobj->SetBranchAddress("tauIDEfficiencySystError", 	&ftauIDEfficiencySystError	       );   //&b_ftauIDEfficiencySystError);
   			tobj->SetBranchAddress("topologyEfficiencySystError", 	&ftopologyEfficiencySystError	       );   //&b_ftopologyEfficiencySystError);
	
	  		for ( Int_t iev=0; iev<Nentries; iev++ ) {
				nb=tobj->GetEntry(iev);
	  		/*logfile<<"Channel  "<<fntuple->GetName()<<endl;
			logfile<<"**********************************************"<<endl;
			logfile<<"Process  "<<teffs<<endl;
			logfile<<"**********************************************"<<endl;
			logfile<<"branchingFraction\t FilterEfficiency\t xsection  "<<endl;
			logfile<<fbranchingFraction<<"\t\t\t "<< fFilterEfficiency<<"\t\t\t "<<fxsection  <<endl;
			logfile<<"acceptanceEfficiency\t electronIDEfficiency\t muonIDEfficiency\t tauIDEfficiency\t topologyEfficiency "<<endl; 	
			logfile<<facceptanceEfficiency<<"\t\t\t " <<felectronIDEfficiency<<"\t\t "<< fmuonIDEfficiency<<"\t\t\t "<<ftauIDEfficiency<<"\t\t "<< ftopologyEfficiency <<endl;
			logfile<<"branchingFractionStatError\t acceptanceEfficiencyStatError\t electronIDEfficiencyStatError\t muonIDEfficiencyStatError\t tauIDEfficiencyStatError\t topologyEfficiencyStatError "<<endl; 
			logfile<<fbranchingFractionSystError<<"\t\t\t\t "<<facceptanceEfficiencySystError<<"\t\t\t\t "<<felectronIDEfficiencySystError<<"\t\t\t\t "<< fmuonIDEfficiencySystError<<"\t\t\t\t "<<ftauIDEfficiencySystError<<"\t\t\t\t "<< ftopologyEfficiencySystError<<endl; 
     			*/      
			if(histname1=="ZprimeTree"){
				toteff=fFilterEfficiency;
				toteff = toteff*(facceptanceEfficiency * felectronIDEfficiency * fmuonIDEfficiency * ftauIDEfficiency * ftopologyEfficiency);
			//Channel[ichannel]<<"  Total Efficiencies  = "<<toteff<<endl;
				fProcess2[0].SetTotEff(toteff);
				fProcess2[0].SetSigma(fxsection);
				if(ichannel==0)fProcess2[0].SetSigma(fxsection*0.000000000001);
				fProcess2[0].SetBR(fbranchingFraction);
				fProcess2[0].SetFilterEfficiency(fFilterEfficiency);
				fProcess2[0].SetacceptanceEfficiency(facceptanceEfficiency);
				fProcess2[0].SetelectronIDEfficiency(felectronIDEfficiency);
				fProcess2[0].SetmuonIDEfficiency(fmuonIDEfficiency);
				fProcess2[0].SettauIDEfficiency(ftauIDEfficiency);
				fProcess2[0].SettopologyEfficiency(ftopologyEfficiency);
				
				
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
						

				prescl=lumi*toteff*(fbranchingFraction*fxsection);
				fProcess2[0].SetPreScl(prescl);
			}
			if(histname1!="ZprimeTree"){
			
				toteff=fFilterEfficiency;
				toteff = toteff*(facceptanceEfficiency * felectronIDEfficiency * fmuonIDEfficiency * ftauIDEfficiency * ftopologyEfficiency);
			//Channel[ichannel]<<"  Total Efficiencies  = "<<toteff<<endl;
				fProcess[ch].SetTotEff(toteff);
				fProcess[ch].SetSigma(fxsection);
				if(ichannel==0) fProcess[ch].SetSigma(fxsection*0.000000000001);
				fProcess[ch].SetBR(fbranchingFraction);
				fProcess[ch].SetFilterEfficiency(fFilterEfficiency);
				fProcess[ch].SetacceptanceEfficiency(facceptanceEfficiency);
				fProcess[ch].SetelectronIDEfficiency(felectronIDEfficiency);
				fProcess[ch].SetmuonIDEfficiency(fmuonIDEfficiency);
				fProcess[ch].SettauIDEfficiency(ftauIDEfficiency);
				fProcess[ch].SettopologyEfficiency(ftopologyEfficiency);
				
				
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
						

				prescl=lumi*toteff*(fbranchingFraction*fxsection);
				fProcess[ch].SetPreScl(prescl);
			//logfile<<"  Total Efficiencies  = "<<fProcess[ch].GetTotEff()<<"  Prescale (luminosity*tot_eff*BR*sigma) = "<<fProcess[ch].GetPreScl()<<endl;
			//closing scanning through tree
  			ch = ch+1;
			}
			}
		//cout<<"N Process "<<ch<<endl;
	 	//end looking at trees
		}
		
		//cout<<"N Process "<<ch<<endl;
	
		fChannel[ichannel].SetNProcess(ch);
		
		if (obj->IsA()->InheritsFrom("TDirectoryFile")) {
	  		TDirectory *dir = (TDirectory*)obj;
	  		dir->cd(dir->GetPath());
	  		TIter nexth(dir->GetListOfKeys());
	  
			while ((keyh = (TKey*)nexth())) {
     	      
        			TObject *objh = keyh->ReadObj();
				//cout<<"inside TDirFl  "<<endl;
       
				if ( objh->IsA()->InheritsFrom( "TH1F" ) ) {
	
      					//cout<<"its a histogram   ";
      
        				string histnameh = objh->GetName();
	 				//cout<<histnameh<<endl;
	 				TH1F* hobj = (TH1F*)objh->Clone("hobj");
				
					fProcess[ipr].SetProcessTitle(histname1);
				
					//histograms Alexei's get_DefaultTemplates()
	 				if(histname1=="ZprimeTemplateDirectory"){
						
						fProcess2[0].SetProcessShape(hobj);
				 		if (fProcess2[0].GetProcessShape()->Integral() > 0.){
 	  						fProcess2[0].GetProcessShape()->Scale(1./fProcess2[0].GetProcessShape()->Integral());
							//cout<<"Process :"<<fProcess2[0].GetProcessTitle()<<"  integral  =  "<< fProcess2[0].GetProcessShape()->Integral()<<endl;
						}
				
						fProcess2[0].GetProcessShape()->Rebin(RebinFactor);
						
					}
				if(histname1!="ZprimeTemplateDirectory"){
				
					TH1F* hobj2 = (TH1F*)objh->Clone("hobj2");
					//fProcess[ipr].SetProcessTitle(histname1);
					fProcess[ipr].SetProcessShape(hobj2);
				 	if (fProcess[ipr].GetProcessShape()->Integral() > 0.){
 	  					fProcess[ipr].GetProcessShape()->Scale(1./fProcess[ipr].GetProcessShape()->Integral());
						//cout<<"Process :"<<fProcess[ipr].GetProcessTitle()<<"  integral  =  "<< fProcess[ipr].GetProcessShape()->Integral()<<endl;
					}
				
					fProcess[ipr].GetProcessShape()->Rebin(RebinFactor);
				ipr=ipr+1;
				}
			//	logfile<<"  Process "<<histname1<<"  Number of events  = "<<fProcess[ipr].GetProcessShape()->Integral()<<endl;
				} 
				
			}

		}

	
	}
	//cout<<"IPR "<<ipr<<endl;
	//for(int p=0; p<ipr;p++){
	//cout<<fProcess[p].GetProcessShape()->Integral()<<endl;
	//}
	fChannel[ichannel].SetProcessProp(fProcess, ipr);
	fChannel[ichannel].SetSignalProp(fProcess2);
	ichannel=ichannel+1;
	//fntuple->Close();
	
	}
	
	//ichannel=ichannel+1;
	//fntuple->Close();
	
	fclose(fTreeFiles);
	
	  //TEST
	//fProcess2 = fChannel[1].GetProcessProp(0);
	//cout<<"Title "<<fProcess2.GetProcessTitle()<<endl;
	/*TH1F* default_temp_bgtot = (TH1F*)fProcess2.GetProcessShape()->Clone("default_temp_totbg");
	cout<<"Testing"<<endl;
	cout<<"Integral "<<default_temp_bgtot->Integral()<<endl;
	TCanvas *c3 = new TCanvas("c3","", 100, 100, 1000, 1000);
	default_temp_bgtot->DrawCopy();
	*/
	TCanvas *c1 = new TCanvas("c1","", 100, 100, 1000, 1000);
	c1->Divide(1,3);
	TCanvas *c2 = new TCanvas("c2","", 100, 100, 1000, 1000);
	
	TCanvas *c3 = new TCanvas("c3","", 100, 100, 1000, 1000);
	TCanvas *c4 = new TCanvas("c4","", 100, 100, 1000, 1000);
	
	TCanvas *cn3 = new TCanvas("cn3","", 100, 100, 1000, 1000);
  	cn3->Clear();
  	cn3->Divide(1,2); 
  
  	TCanvas *cn2 = new TCanvas("cn2","", 100, 100, 1000, 1000);
  	cn2->Clear();
  	cn2->Divide(1,2);
  

  	TCanvas *cn1 = new TCanvas("cn1","", 100, 100, 1000, 1000);
  	cn1->Clear();
  	cn1->Divide(1,2);
	
		TCanvas *c5 = new TCanvas("c5","c5", 100, 100, 1000, 1000);
  	c5->Clear();
  	c5->Divide(2,3);
	TCanvas *c6 = new TCanvas("c6","c6", 100, 100, 1000, 1000);
  	c6->Clear();
  	c6->Divide(1,3);
	
	THStack *hs = new THStack("hs","MuTau Channel");
	THStack *ht = new THStack("ht","eTau Channel");
	THStack *hu = new THStack("hu","eMu Channel");
	
	//Calculations
	//TPrProp fProcess2[15];
	
	int nBins2 = 250;
	double Min = 0;
	double Max = 25;
	
	TH1F* LogLVsSigma1 = new TH1F("LogLVsSigma1","LogLVsSigma1", nBins2, Min, Max);
	TH1F* JointlklvSigma= new TH1F("JointlklvSigma","JointlklvSigma", nBins2, Min, Max);
		for(int lim=0;lim<npseudo;lim++){
			
		if(lim==100) cout<<"******************   PseudoExperiment "<<lim<<"   ******************"<<endl;
		if(lim==150) cout<<"******************   PseudoExperiment "<<lim<<"   ******************"<<endl;
		if(lim==250) cout<<"******************   PseudoExperiment "<<lim<<"   ******************"<<endl;
		for(ichannel=0;ichannel<3;ichannel++){
		LogLVsSigma1->Reset();
			
		//TH1F* LogLVsSigmaMC=(TH1F*)fChannel[ichannel].Likelihood(ichannel, fChannel[ichannel], mcint, LogLVsSigma1)->Clone("LogLVsSigmaMC");
		TH1F* LogLVsSigmaMC0=(TH1F*)fChannel[ichannel].Likelihood(ichannel, fChannel[ichannel], 0, LogLVsSigma1)->Clone("LogLVsSigmaMC0");
		/*TH1F* LogLVsSigmaMC1=(TH1F*)fChannel[ichannel].Likelihood(ichannel, fChannel[ichannel], 100, LogLVsSigma1)->Clone("LogLVsSigmaMC1");
		TH1F* LogLVsSigmaMC2=(TH1F*)fChannel[ichannel].Likelihood(ichannel, fChannel[ichannel], 500, LogLVsSigma1)->Clone("LogLVsSigmaMC2");
		TH1F* LogLVsSigmaMC3=(TH1F*)fChannel[ichannel].Likelihood(ichannel, fChannel[ichannel], 1000, LogLVsSigma1)->Clone("LogLVsSigmaMC3");
		//cout<<"Integral "<<LogLVsSigmaMC->Integral()<<endl;
		
		if(ichannel==0)c1->cd(1);
		if(ichannel==1)c1->cd(2);
		if(ichannel==2)c1->cd(3);
		 LogLVsSigmaMC->DrawCopy();
		 LogLVsSigmaMC0->DrawCopy("same");
		 LogLVsSigmaMC1->DrawCopy("same");
		 LogLVsSigmaMC2->DrawCopy("same");
		 LogLVsSigmaMC3->DrawCopy("same");
	
		if(lim<5){
			 if(ichannel==0){
			  	cn1->cd(1);
			 	if(lim==0){
			 	  //LogLVsSigma->GetYaxis()->SetRangeUser(0.,1);
				  LogLVsSigmaMC->GetXaxis()->SetTitle("Sigma");
				  LogLVsSigmaMC->GetYaxis()->SetTitle("Likelihood");
				  LogLVsSigmaMC->SetTitle("Likelihood");
				  //LogLVsSigmaMC->DrawCopy();
				 // LogLVsSigmaMC0->SetFillColor(kCyan);
				 LogLVsSigmaMC0->SetLineColor(kBlue);
				 LogLVsSigmaMC0->DrawCopy();
				  LogLVsSigmaMC->DrawCopy("same");
				  cn1->SetLogy();
				  cn1->BuildLegend();
				  cn1->Update();
				
				}
			 if(lim>0)LogLVsSigmaMC->DrawCopy("same");
			 }
			 if(ichannel==1){
			  	cn2->cd(1);
			 	if(lim==0){
			 	  //LogLVsSigma->GetYaxis()->SetRangeUser(0.,1);
				  LogLVsSigmaMC->GetXaxis()->SetTitle("Sigma");
				  LogLVsSigmaMC->GetYaxis()->SetTitle("Likelihood");
				  LogLVsSigmaMC->SetTitle("Likelihood");
				  
				  LogLVsSigmaMC0->SetLineColor(kBlue);
				 // LogLVsSigmaMC0->SetFillColor(kCyan);
				  LogLVsSigmaMC0->DrawCopy();
				  LogLVsSigmaMC->DrawCopy("same");
				  cn2->SetLogy();
				  cn2->BuildLegend();
				  cn2->Update();
				  
				}
			 if(lim>0)LogLVsSigmaMC->DrawCopy("same");
			 }
 			 if(ichannel==2){
			  	cn3->cd(1);
			 	if(lim==0){
			 	  //LogLVsSigma->GetYaxis()->SetRangeUser(0.,1);
				  LogLVsSigmaMC->GetXaxis()->SetTitle("Sigma");
				  LogLVsSigmaMC->GetYaxis()->SetTitle("Likelihood");
				  LogLVsSigmaMC->SetTitle("Likelihood");
				 // LogLVsSigmaMC->DrawCopy();
				  LogLVsSigmaMC0->SetLineColor(kBlue);
				 // LogLVsSigmaMC0->SetFillColor(kCyan);
				  LogLVsSigmaMC0->DrawCopy();
				  LogLVsSigmaMC->DrawCopy("same");
				  cn3->SetLogy();
				  cn3->BuildLegend();
				  cn3->Update();
				 
				  
				}
			 if(lim>0)LogLVsSigmaMC->DrawCopy("same");
			 } 
			 }*/	
		
			if(ichannel==0) JointlklvSigma->Add(LogLVsSigmaMC0);
			if(ichannel>0) JointlklvSigma->Multiply(LogLVsSigmaMC0);
		 	
			sigma95 = fChannel[ichannel].Limit95CL(LogLVsSigmaMC0);
			/*sigma951 = fChannel[ichannel].Limit95CL(LogLVsSigmaMC1);
			sigma952 = fChannel[ichannel].Limit95CL(LogLVsSigmaMC2);
			sigma953 = fChannel[ichannel].Limit95CL(LogLVsSigmaMC3);
			*/
			if(ichannel==0){
			Limit_95CL_0->Fill(sigma95);
			/*Limit_95CL_01->Fill(sigma951);
			Limit_95CL_02->Fill(sigma952);
			Limit_95CL_03->Fill(sigma953);*/
			}
			if(ichannel==1){
			Limit_95CL_1->Fill(sigma95);
			/*Limit_95CL_12->Fill(sigma952);
			Limit_95CL_13->Fill(sigma953);
			Limit_95CL_1->Fill(sigma95);*/
			}
			if(ichannel==2){
			Limit_95CL_2->Fill(sigma95);
			/*Limit_95CL_21->Fill(sigma951);
			Limit_95CL_22->Fill(sigma952);
			Limit_95CL_23->Fill(sigma953);*/
			}
			
			
		}
			joint_sigma95 = fChannel[0].Limit95CL(JointlklvSigma);
			Joint_Limit_95CL->Fill(joint_sigma95);
	}
			fChannel[0].SetLimit95CLH(Limit_95CL_0);
			fChannel[1].SetLimit95CLH(Limit_95CL_1);
			fChannel[2].SetLimit95CLH(Limit_95CL_2);
			
			c6->cd(1);
			Limit_95CL_0->Draw();
			/*Limit_95CL_01->Draw("sames");
			Limit_95CL_02->Draw("sames");
			Limit_95CL_03->Draw("sames");*/
			c6->cd(2);
			Limit_95CL_1->Draw();
			/*Limit_95CL_11->Draw("sames");
			Limit_95CL_12->Draw("sames");
			Limit_95CL_13->Draw("sames");*/
			c6->cd(3);
			Limit_95CL_2->Draw();
			/*Limit_95CL_21->Draw("sames");
			Limit_95CL_22->Draw("sames");
			Limit_95CL_23->Draw("sames");*/
			
			c5->cd(2);
			fChannel[0].GetLimit95CLH()->SetTitle("muTau 95%CL Limit, suppressed BG");
			fChannel[0].GetLimit95CLH()->GetXaxis()->SetTitle("Sigma");
			fChannel[0].GetLimit95CLH()->DrawCopy();
			c5->cd(4);
			fChannel[1].GetLimit95CLH()->SetTitle("eTau 95%CL Limit, suppressed BG");
			fChannel[1].GetLimit95CLH()->GetXaxis()->SetTitle("Sigma");
			fChannel[1].GetLimit95CLH()->DrawCopy();
	     		c5->cd(6);
			fChannel[2].GetLimit95CLH()->SetTitle("eMu 95%CL Limit, suppressed BG");
			fChannel[2].GetLimit95CLH()->GetXaxis()->SetTitle("Sigma");
			fChannel[2].GetLimit95CLH()->DrawCopy();
			c1->cd();
			Joint_Limit_95CL->Draw();

				
	
}

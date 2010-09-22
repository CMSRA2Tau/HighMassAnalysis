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

#include "TH1F.h"
#include "TH2F.h"
#include "TRandom.h"
#include "TVector.h"
#include "TLorentzVector.h" 
#include "TList.h" 
#include "TChain.h"
#include "TMath.h"

#include "TChProp.h"
#include "TPrProp.h"

TH1F* TChProp::Likelihood(int ichn, TChProp fChannelPro, int mcintg, TH1F* LogLVsSigma1){

	
	Float_t distIntegral=0;
	Float_t fitValue = 0;
	Float_t distValue_num=0;
	Float_t lkl=0;
	Float_t lkl2 = 0;
	char hist_name[100];
	
	
	float default_background_scale=1.0;
	Float_t totscl0=0;
	Float_t totscl1=0;
	
	double lumi = 50.0;
	
	TPrProp fProcessd;
	
	TH1F* default_temp_bgtot = (TH1F*)fChannelPro.GetProcessProp(1).GetProcessShape()->Clone("default_temp_totbg");
	TH1F* default_temp_bgtot_0 = (TH1F*)fChannelPro.GetProcessProp(1).GetProcessShape()->Clone("default_temp_totbg_0");
	TH1F* def_bg = (TH1F*)fChannelPro.GetProcessProp(1).GetProcessShape()->Clone("def_bg");
	TH1F* default_temp_signal =(TH1F*)fChannelPro.GetSignalProp().GetProcessShape()->Clone("default_temp_signal");
	TH1F* default_temp_signal_0 =(TH1F*)fChannelPro.GetSignalProp().GetProcessShape()->Clone("default_temp_signal_0");
	
	TH1F* LogLVsSigma2 = (TH1F*)LogLVsSigma1->Clone("LogLVsSigma2");
	Int_t nBins = default_temp_signal->GetNbinsX();
	
	TH1F* pseudodata_MassDistribution = (TH1F*)default_temp_bgtot_0->Clone("pseudodata_MassDistribution");
	TH1F* pseudodata_0 = (TH1F*)default_temp_bgtot->Clone("pseudodata_0");
	TH1F* pseudodata_1 = (TH1F*)default_temp_bgtot->Clone("pseudodata_1");
	TH1F* pseudodata_2 = (TH1F*)default_temp_bgtot->Clone("pseudodata_2");
	pseudodata_MassDistribution->Reset();
	pseudodata_0->Reset();
	pseudodata_1->Reset();
	pseudodata_2->Reset();
	default_temp_bgtot->Reset();
	default_temp_signal->Reset();
	default_temp_bgtot_0->Reset();
	default_temp_signal_0->Reset();
	Double_t step = LogLVsSigma1->GetBinWidth(1); //(Max - Min)/nBins2;
	Int_t nbins = LogLVsSigma1->GetNbinsX();
		
	for( int k=0; k<fChannelPro.GetNProcess(); k++){
		totscl0 =lumi*fChannelPro.GetProcessProp(k).GetBR()*fChannelPro.GetProcessProp(k).GetFilterEfficiency()*fChannelPro.GetProcessProp(k).GetacceptanceEfficiency()*fChannelPro.GetProcessProp(k).GetelectronIDEfficiency()*fChannelPro.GetProcessProp(k).GetmuonIDEfficiency()*fChannelPro.GetProcessProp(k).GettopologyEfficiency()*fChannelPro.GetProcessProp(k).GettauIDEfficiency();
		//cout<<" xsection  "<<fChannelPro.GetProcessProp(k).GetSigma()<<endl;
		//if(k>0) 
		default_temp_bgtot_0->Add(fChannelPro.GetProcessProp(k).GetProcessShape(),totscl0*fChannelPro.GetProcessProp(k).GetSigma());			
	}
		//default_temp_bgtot_0->Rebin(rebinfactor);
		default_temp_bgtot_0->Scale(default_background_scale);
		//cout<<"background Integral = "<<default_temp_bgtot_0->Integral()<<endl;
		distIntegral = int(gRandom->Poisson(default_temp_bgtot_0->Integral()));
		pseudodata_MassDistribution->Reset();
	 	pseudodata_MassDistribution->FillRandom(default_temp_bgtot_0, distIntegral);
		//pseudodata_MassDistribution->Scale(0.0);
	
		
		if(mcintg==0){
			default_temp_signal->Reset();
	        	totscl1=lumi*fChannelPro.GetSignalProp().GetBR()*fChannelPro.GetSignalProp().GetFilterEfficiency()*fChannelPro.GetSignalProp().GetacceptanceEfficiency()*fChannelPro.GetSignalProp().GetelectronIDEfficiency()*fChannelPro.GetSignalProp().GetmuonIDEfficiency()*fChannelPro.GetSignalProp().GettopologyEfficiency()*fChannelPro.GetSignalProp().GettauIDEfficiency();
			default_temp_signal->Add(fChannelPro.GetSignalProp().GetProcessShape(),totscl1);		
			
			
			LogLVsSigma1->Reset();
			for (int in= 1; in <= nbins; in++) {
	        		double value[5];
				Double_t currentSigma = LogLVsSigma1->GetBinCenter(in);// Min + (float(i)-0.5)*step;
				lkl=0;
				for (int i = 1; i <= nBins; i++){
				fitValue =(currentSigma*default_temp_signal->GetBinContent(i)+default_temp_bgtot_0->GetBinContent(i));
				distValue_num = pseudodata_MassDistribution->GetBinContent(i);

                		if(fitValue <1e-50) continue;		
		  		lkl = lkl + distValue_num*log(fitValue) - fitValue - TMath::LnGamma(distValue_num + 1);
				}
				//cout<<" lkl "<<lkl<<" exp(lkl) "<<exp(lkl)<<endl;
			
			LogLVsSigma1->SetBinContent(in, exp(lkl));
			
				
			}
			
		/*	if(ichn==0){
			  c5->cd(1);
			   LogLVsSigma1->SetTitle("muTau, no smearing");
			  LogLVsSigma1->SetLineColor(kCyan);
			  LogLVsSigma1->DrawCopy("same");
			
			}
			if(ichn==1){
			  c5->cd(3);
			  LogLVsSigma1->SetTitle("eTau, no smearing");
			  LogLVsSigma1->SetLineColor(kCyan);
			  LogLVsSigma1->DrawCopy("same");
			 
			}
			if(ichn==2){
			  c5->cd(5);
			  LogLVsSigma1->SetTitle("eMu, no smearing");
			  LogLVsSigma1->SetLineColor(kCyan);
			  LogLVsSigma1->DrawCopy("same");
			  
			}*/
			TH1F* fLikelihood = (TH1F*)LogLVsSigma1->Clone("fLikelihood");
			return(fLikelihood);	
			
		}
		
		if(mcintg>0){
		for( int mcin=0; mcin<mcintg; mcin++){
			default_temp_signal->Reset();
			default_temp_bgtot->Reset();
			/*if(mcin==100) cout<<"MC 100"<<endl;
			if(mcin==200) cout<<"MC 200"<<endl;
			if(mcin==300) cout<<"MC 300"<<endl;
			if(mcin==400) cout<<"MC 400"<<endl;
			if(mcin==500) cout<<"MC 500"<<endl;
			if(mcin==700) cout<<"MC 700"<<endl;
			if(mcin==1000) cout<<"MC 1000"<<endl;*/
			for( int k=0; k<fChannelPro.GetNProcess(); k++){
	        		totscl1=fProcessd.Scaling(fChannelPro.GetProcessProp(k));
				totscl0=fProcessd.Scaling(fChannelPro.GetSignalProp());
				//if(k==0) 
				default_temp_signal->Add(fChannelPro.GetSignalProp().GetProcessShape(),totscl0);
				//if(k>0) 
				default_temp_bgtot->Add(fChannelPro.GetProcessProp(k).GetProcessShape(),totscl1*fChannelPro.GetProcessProp(k).GetSigma());		
			}
			// default_temp_signal->Rebin(rebinfactor);
			// default_temp_bgtot->Rebin(rebinfactor);
			
			default_temp_bgtot->Scale(default_background_scale);
			LogLVsSigma1->Reset();
			for (int in= 1; in <= nbins; in++) {
	        		double value[5];
				Double_t currentSigma = LogLVsSigma1->GetBinCenter(in);// Min + (float(i)-0.5)*step;
				lkl=0;
				for (int i = 1; i <= nBins; i++){
				fitValue =(currentSigma*default_temp_signal->GetBinContent(i)+default_temp_bgtot->GetBinContent(i));
				distValue_num = pseudodata_MassDistribution->GetBinContent(i);

                		if(fitValue <1e-50) continue;		
		  		lkl = lkl + distValue_num*log(fitValue) - fitValue - TMath::LnGamma(distValue_num + 1);
				}
				//cout<<" lkl "<<lkl<<" exp(lkl) "<<exp(lkl)<<endl;
			
			LogLVsSigma1->SetBinContent(in, exp(lkl));	
			}
			if(mcin<20){
			if(ichn==0){
			  	c5->cd(1);
			 	if(mcin==0){
			 	  //LogLVsSigma->GetYaxis()->SetRangeUser(0.,1);
				  LogLVsSigma1->GetXaxis()->SetTitle("Sigma");
				  LogLVsSigma1->GetYaxis()->SetTitle("Likelihood");
				  LogLVsSigma1->SetTitle("muTau Likelihood");
				LogLVsSigma1->SetLineColor(kBlack);
				 LogLVsSigma1->DrawCopy();
				}
			 if(mcin>0)LogLVsSigma1->DrawCopy("same");
			 }
			 if(ichn==1){
			  	c5->cd(3);
			 	if(mcin==0){
			 	  //LogLVsSigma->GetYaxis()->SetRangeUser(0.,1);
				  LogLVsSigma1->GetXaxis()->SetTitle("Sigma");
				  LogLVsSigma1->GetYaxis()->SetTitle("Likelihood");
				  LogLVsSigma1->SetTitle("eTau Likelihood");
				LogLVsSigma1->SetLineColor(kBlack);
				 LogLVsSigma1->DrawCopy();
				}
			 if(mcin>0)LogLVsSigma1->DrawCopy("same");
			 }
			 if(ichn==2){
			  	c5->cd(5);
			 	if(mcin==0){
			 	  //LogLVsSigma->GetYaxis()->SetRangeUser(0.,1);
				  LogLVsSigma1->GetXaxis()->SetTitle("Sigma");
				  LogLVsSigma1->GetYaxis()->SetTitle("Likelihood");
				  LogLVsSigma1->SetTitle("eMu Likelihood");
				LogLVsSigma1->SetLineColor(kBlack);
				 LogLVsSigma1->DrawCopy();
				}
			 if(mcin>0)LogLVsSigma1->DrawCopy("same");
			 }
			 
			 
			 }
			//cout<<"one integration "<<LogLVsSigma1->Integral()<<endl;
			//if(mcin==0) 
			
			//if(mcin>0) 
			
			LogLVsSigma2->Add(LogLVsSigma1);
			
			
			
			//cout<<"adding integrations "<<LogLVsSigma2->Integral()<<endl;
		}
		//cout<<"total integrations integration "<<LogLVsSigma2->Integral()<<endl;
		
		float scl=0;
		if(mcintg==10) scl = 0.1;
		if(mcintg==100) scl = 0.01;
		if(mcintg==500) scl = 0.002;
		if(mcintg==1000) scl = 0.001;
		LogLVsSigma2->Scale(scl);
		//cout<<LogLVsSigma2->Integral()<<endl;
		
		TH1F* fLikelihood = (TH1F*)LogLVsSigma2->Clone("fLikelihood");
		return(fLikelihood);	
		}	
}


void TChProp::SetProcessProp(TPrProp* fprocess, int ipr){
	//string testtitle;
	//fProcessN = fprocess;
	for(int i=0;i<ipr;i++){
	//testtitle = 
	//cout<<"fun "<<fprocess[i].GetProcessTitle()<<endl;;
	//cout<<testtitle<<endl;
	//fSignal.SetProcessTitle("fun");
	//cout<<"TEST"<<endl;
	//fSignal.SetProcessTitle("normal");
	fProcessN[i].SetTotEff(fprocess[i].GetTotEff());
	fProcessN[i].SetPreScl(fprocess[i].GetPreScl());
	
	fProcessN[i].SetFilterEfficiency(fprocess[i].GetFilterEfficiency());
	fProcessN[i].SetacceptanceEfficiency(fprocess[i].GetacceptanceEfficiency());
	fProcessN[i].SetelectronIDEfficiency(fprocess[i].GetelectronIDEfficiency());
	fProcessN[i].SetmuonIDEfficiency(fprocess[i].GetmuonIDEfficiency());
	fProcessN[i].SettauIDEfficiency(fprocess[i].GettauIDEfficiency());
	fProcessN[i].SettopologyEfficiency(fprocess[i].GettopologyEfficiency());
	
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
	
	}
	//cout<<"Integral "<<fProcessN[0].GetProcessShape()->Integral()<<endl;
	//break;
	//}
	}

void TChProp::SetSignalProp(TPrProp* fprocess3){
	//string testtitle;
	//fProcessN = fprocess;
	

	fSignal[0].SetTotEff(fprocess3[0].GetTotEff());
	fSignal[0].SetPreScl(fprocess3[0].GetPreScl());
	
	fSignal[0].SetFilterEfficiency(fprocess3[0].GetFilterEfficiency());
	fSignal[0].SetacceptanceEfficiency(fprocess3[0].GetacceptanceEfficiency());
	fSignal[0].SetelectronIDEfficiency(fprocess3[0].GetelectronIDEfficiency());
	fSignal[0].SetmuonIDEfficiency(fprocess3[0].GetmuonIDEfficiency());
	fSignal[0].SettauIDEfficiency(fprocess3[0].GettauIDEfficiency());
	fSignal[0].SettopologyEfficiency(fprocess3[0].GettopologyEfficiency());
	
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
	
	
	//cout<<"Integral "<<fProcessN[0].GetProcessShape()->Integral()<<endl;
	//break;
	//}
	}


Double_t TChProp::Limit95CL(TH1F* LogLVsSigma){
	Double_t totalIntegral95 = 0.95*LogLVsSigma->Integral();
	//cout<<"Integral95 = "<<totalIntegral95 <<endl;
	Double_t currentIntegral = 0;
	int i = 1;
		
		while(currentIntegral < totalIntegral95){
				currentIntegral = LogLVsSigma->Integral(1, i);
				i++;
	  	}
	  	
		Double_t sigma95 = LogLVsSigma->GetBinLowEdge(i+1);
		return(sigma95);
			
}

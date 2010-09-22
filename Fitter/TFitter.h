/*********************************************************
********************************************************** 
*** 	 April 23,2010					**
***  I. Suarez						**
***	define Class TFitter  				**
*** 			     				**
**********************************************************
*********************************************************/

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TPrProp.h>
#include <TChProp.h>
#ifndef TFitter_h
#define TFitter_h


class TFitter {
private:
TPrProp *fProcess[15];
TChProp *fChannel;

public :

TFitter(){};
virtual ~TFitter(){};

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

//index glossary
// {i}{j}
//i = channel id
//   0=mu+tau  j=0,1,2,9
//   1=e+tau   j=3,4,5,6,7,8,9,1
//   2=e+mu    j=0,1,2,3,10
//   3=tautau  j=1,3,4,11,10,12
//j = process id
//   0=Ztt, 1=WJets, 2=IncMu, 3=Zprime, 4=Ztautau, 5=Zee
//   6=PhotonJets, 7=QCD20to30, 8=QCD30to80, 9=QCD80to170
//   10=TTBar, 11=Zmumu, 12=QCD

	/* Float_t         branchingFraction;
   	Float_t         FilterEfficiency;
   	Float_t         xsection;
   	Float_t         acceptanceEfficiency;
   	Float_t         electronIDEfficiency;
   	Float_t         muonIDEfficiency;
   	Float_t         tauIDEfficiency;
   	Float_t         topologyEfficiency;
   	Float_t         branchingFractionStatError;
   	Float_t         acceptanceEfficiencyStatError;
   	Float_t         electronIDEfficiencyStatError;
   	Float_t         muonIDEfficiencyStatError;
   	Float_t         tauIDEfficiencyStatError;
   	Float_t         topologyEfficiencyStatError;
   	Float_t         branchingFractionSystError;
   	Float_t         acceptanceEfficiencySystError;
   	Float_t         electronIDEfficiencySystError;
   	Float_t         muonIDEfficiencySystError;
   	Float_t         tauIDEfficiencySystError;
   	Float_t         topologyEfficiencySystError;

   TBranch        *b_branchingFraction;   //!
   TBranch        *b_FilterEfficiency;   //!
   TBranch        *b_xsection;   //!
   TBranch        *b_acceptanceEfficiency;   //!
   TBranch        *b_electronIDEfficiency;   //!
   TBranch        *b_muonIDEfficiency;   //!
   TBranch        *b_tauIDEfficiency;   //!
   TBranch        *b_topologyEfficiency;   //!
   TBranch        *b_branchingFractionStatError;   //!
   TBranch        *b_acceptanceEfficiencyStatError;   //!
   TBranch        *b_electronIDEfficiencyStatError;   //!
   TBranch        *b_muonIDEfficiencyStatError;   //!
   TBranch        *b_tauIDEfficiencyStatError;   //!
   TBranch        *b_topologyEfficiencyStatError;   //!
   TBranch        *b_branchingFractionSystError;   //!
   TBranch        *b_acceptanceEfficiencySystError;   //!
   TBranch        *b_electronIDEfficiencySystError;   //!
   TBranch        *b_muonIDEfficiencySystError;   //!
   TBranch        *b_tauIDEfficiencySystError;   //!
   TBranch        *b_topologyEfficiencySystError;   //!
  */
   virtual void     Loop();
   
  

};

#endif

	/*tobj->SetBranchAddress("branchingFraction", 		&branchingFraction, 	&b_branchingFraction);
   	tobj->SetBranchAddress("FilterEfficiency", 		&FilterEfficiency, 	&b_FilterEfficiency);
	tobj->SetBranchAddress("xsection", 			&xsection, 		&b_xsection);
   	tobj->SetBranchAddress("acceptanceEfficiency", 		&acceptanceEfficiency, 	&b_acceptanceEfficiency);
   	tobj->SetBranchAddress("electronIDEfficiency", 		&electronIDEfficiency, 	&b_electronIDEfficiency);
   	tobj->SetBranchAddress("muonIDEfficiency", 		&muonIDEfficiency, 	&b_muonIDEfficiency);
   	tobj->SetBranchAddress("tauIDEfficiency", 		&tauIDEfficiency, 	&b_tauIDEfficiency);
   	tobj->SetBranchAddress("topologyEfficiency", 		&topologyEfficiency, 	&b_topologyEfficiency);
   	tobj->SetBranchAddress("branchingFractionStatError", 	&branchingFractionStatError, &b_branchingFractionStatError);
   	tobj->SetBranchAddress("acceptanceEfficiencyStatError", &acceptanceEfficiencyStatError, &b_acceptanceEfficiencyStatError);
   	tobj->SetBranchAddress("electronIDEfficiencyStatError", &electronIDEfficiencyStatError, &b_electronIDEfficiencyStatError);
   	tobj->SetBranchAddress("muonIDEfficiencyStatError", 	&muonIDEfficiencyStatError, &b_muonIDEfficiencyStatError);
   	tobj->SetBranchAddress("tauIDEfficiencyStatError", 	&tauIDEfficiencyStatError, &b_tauIDEfficiencyStatError);
   	tobj->SetBranchAddress("topologyEfficiencyStatError", 	&topologyEfficiencyStatError, &b_topologyEfficiencyStatError);
   	tobj->SetBranchAddress("branchingFractionSystError", 	&branchingFractionSystError, &b_branchingFractionSystError);
   	tobj->SetBranchAddress("acceptanceEfficiencySystError", &acceptanceEfficiencySystError, &b_acceptanceEfficiencySystError);
   	tobj->SetBranchAddress("electronIDEfficiencySystError", &electronIDEfficiencySystError, &b_electronIDEfficiencySystError);
   	tobj->SetBranchAddress("muonIDEfficiencySystError", 	&muonIDEfficiencySystError, &b_muonIDEfficiencySystError);
   	tobj->SetBranchAddress("tauIDEfficiencySystError", 	&tauIDEfficiencySystError, &b_tauIDEfficiencySystError);
   	tobj->SetBranchAddress("topologyEfficiencySystError", 	&topologyEfficiencySystError, &b_topologyEfficiencySystError);	
	*/

// Authors: Andres Florez, Alfredo Gurrola, Eduardo Luiggi, Chi Nhan Nguyen

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/LorentzVectorFwd.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/normalizedPhi.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "CLHEP/Random/RandGauss.h"
#include "CommonTools/CandUtils/interface/Booster.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include <Math/VectorUtil.h>

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Common/interface/Handle.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"


#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>

using namespace std;
using namespace edm;

typedef reco::Candidate::LorentzVector LorentzVector;
typedef std::vector< reco::Candidate::LorentzVector > LVCollection;

class HiMassTauAnalysis : public EDAnalyzer {
public:
  explicit HiMassTauAnalysis(const ParameterSet&);
  ~HiMassTauAnalysis();


private:
  virtual void beginJob() ;
  virtual void analyze(const Event&, const EventSetup&);
  virtual void endJob() ;

  void getCollections(const Event&, const EventSetup&);
  void fillHistograms();
  void fillNtuple();
  void getEventFlags(const Event&);
  bool passEventSelectionSequence();
  //bool passGenTauCuts(const LorentzVector&);
  bool passRecoTriggerCuts(const Event&);
  bool passRecoVertexCuts(const reco::Vertex&);
  bool passRecoTauCuts(const pat::Tau&,int);
  bool passRecoMuonCuts(const pat::Muon&,int);
  bool passRecoElectronCuts(const pat::Electron&,int);
  bool passRecoJetCuts(const pat::Jet&,int);
  bool passRecoFirstLeadingJetCuts(const pat::Jet&,int);
  bool passRecoSecondLeadingJetCuts(const pat::Jet&,int);
  bool passTopologyCuts(const pat::Tau&, int, const pat::Muon&, int);
  bool passTopologyCuts(const pat::Tau&, int, const pat::Electron&, int);
  bool passTopologyCuts(const pat::Electron&, int,const pat::Muon&, int);
  bool passTopologyCuts(const pat::Muon&, int,const pat::Muon&, int);
  bool passTopologyCuts(const pat::Tau&, int, const pat::Tau&, int);
  bool passTopologyCuts(const pat::Electron&, int,const pat::Electron&, int);
  bool passSusyTopologyCuts(int, int);
  pair<bool, reco::Candidate::LorentzVector> matchToGen(const pat::Electron&);
  pair<bool, reco::Candidate::LorentzVector> matchToGen(const pat::Muon&);
  pair<bool, reco::Candidate::LorentzVector> matchToGen(const pat::Tau&);
  double CalculatePZeta(const pat::Tau&, int, const pat::Muon&, int);
  double CalculatePZeta(const pat::Tau&, int, const pat::Electron&, int);
  double CalculatePZeta(const pat::Electron&, int,const pat::Muon&, int);
  double CalculatePZeta(const pat::Tau&, int, const pat::Tau&, int);
  double CalculatePZeta(const pat::Muon&, int,const pat::Muon&, int);
  double CalculatePZeta(const pat::Electron&, int,const pat::Electron&, int);
  double CalculatePZetaVis(const pat::Tau&, int, const pat::Muon&, int);
  double CalculatePZetaVis(const pat::Tau&, int, const pat::Electron&, int);
  double CalculatePZetaVis(const pat::Electron&, int, const pat::Muon&, int);
  double CalculatePZetaVis(const pat::Tau&, int, const pat::Tau&, int);
  double CalculatePZetaVis(const pat::Muon&, int,const pat::Muon&, int);
  double CalculatePZetaVis(const pat::Electron&, int, const pat::Electron&, int);
  std::pair<bool, reco::Candidate::LorentzVector> CalculateThe4Momentum(const pat::Tau&, int, const pat::Electron&, int);
  std::pair<bool, reco::Candidate::LorentzVector> CalculateThe4Momentum(const pat::Tau&, int, const pat::Muon&, int);
  std::pair<bool, reco::Candidate::LorentzVector> CalculateThe4Momentum(const pat::Electron&, int,const pat::Muon&, int);
  std::pair<bool, reco::Candidate::LorentzVector> CalculateThe4Momentum(const pat::Tau&, int, const pat::Tau&, int);
  std::pair<bool, reco::Candidate::LorentzVector> CalculateThe4Momentum(const pat::Muon&, int, const pat::Muon&, int);
  std::pair<bool, reco::Candidate::LorentzVector> CalculateThe4Momentum(const pat::Electron&, int, const pat::Electron&, int);
  double CalculateLeptonMetMt(const pat::Muon&, int);
  double CalculateLeptonMetMt(const pat::Electron&, int);
  double CalculateLeptonMetMt(const pat::Tau&, int);
  double alphaRatio(double);
  std::pair<int, double> CalculateTauTrackIsolation(const pat::Tau&);
  std::pair<int, double> CalculateTauTrackIsolation(const pat::Tau&, float, float);
  std::pair<int, double> CalculateTauEcalIsolation(const pat::Tau&);
  std::pair<int, double> CalculateTauEcalIsolation(const pat::Tau&, float, float);
  int CalculateNumberSignalTauGammas(const pat::Tau&);
  reco::Candidate::LorentzVector CalculateTauSignalTracksMass(const pat::Tau&);
  reco::Candidate::LorentzVector CalculateTauSignalTracksAndGammasMass(const pat::Tau&);
  std::pair<reco::Candidate::LorentzVector,int> CalculateTauSignalTracksAndPiZerosMass(const pat::Tau&);
  std::pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> SmearLightLepton(const pat::Muon&);
  std::pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> SmearLightLepton(const pat::Electron&);
  std::pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> SmearTau(const pat::Tau&);
  std::pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> SmearJet(const pat::Jet&);
  std::pair<bool, std::pair<float, float> > isZee(reco::Candidate::LorentzVector);
  void InitializeInfoForPDFSystematicUncertaintites();
  void bookHistograms();
  void setMapSelectionAlgoIDs();
  void initMapSelectionCounters();
  void printEfficiency();
  void setupBranches();
  bool isInTheCracks(float);
  std::pair<unsigned int, unsigned int> getMatchedPdgId(float, float, float, int);
  void clearVectors();
  void initializeVectors();

  // ----------member data ---------------------------

  //-----Generator level Inputs 
  InputTag _GenParticleSource;

  //-----Inputs to determine which channel to analyze
  bool _AnalyzeTauForLeg1;
  bool _AnalyzeMuonForLeg1;
  bool _AnalyzeElectronForLeg1;
  bool _AnalyzeTauForLeg2;
  bool _AnalyzeMuonForLeg2;
  bool _AnalyzeElectronForLeg2;

  //-----Reco Tau Inputs 
  InputTag _RecoTauSource;
  double _RecoTauPtMinCut;
  double _RecoTauPtMaxCut;
  double _RecoTauEtaCut;
  double _RecoTauLeadTrackThreshold;
  double _RecoTauSigGamThreshold;
  double _RecoTauIsoDeltaRCone;
  double _RecoTauTrackIsoTrkThreshold;
  double _RecoTauGammaIsoGamThreshold;
  double _RecoTauTrackIsoSumPtMinCutValue;
  double _RecoTauTrackIsoSumPtMaxCutValue;
  double _RecoTauEcalIsoSumPtMinCutValue;
  double _RecoTauEcalIsoSumPtMaxCutValue;
  double _RecoTauEcalIsoRphiForEllipse;
  double _RecoTauEcalIsoRetaForEllipse;
  double _RecoTauSignal3ProngAndGammasMassMinCutValue;
  double _RecoTauSignal3ProngAndGammasMassMaxCutValue;
  double _RecoTauSignal1ProngAndGammasMassForPionMinCutValue;
  double _RecoTauSignal1ProngAndGammasMassForPionMaxCutValue;
  double _RecoTauSignal1ProngAndGammasMassForKaonVetoMinCutValue;
  double _RecoTauSignal1ProngAndGammasMassForKaonVetoMaxCutValue;
  int _RecoTauTrackNisoMax;
  int _RecoTauEcalNisoMax;
  int _RecoTauLeadTrackMinHits;
  bool _DoRecoTauDiscrByIsolation;
  bool _UseRecoTauDiscrByIsolationFlag;
  bool _UseRecoTauIsoSumPtInsteadOfNiso;
  bool _UseRecoTauEllipseForEcalIso;
  bool _DoRecoTauDiscrByLeadTrack;
  bool _DoRecoTauDiscrByLeadTrackNhits;
  bool _UseRecoTauDiscrByLeadTrackFlag;
  bool _DoRecoTauDiscrAgainstElectron;
  bool _DoRecoTauDiscrAgainstMuon;
  bool _DoRecoTauDiscrByCrackCut;
  bool _DoRecoTauDiscrBySignalTracksAndGammasMass;
  bool _SetTANC;
  string _RecoTauDiscrByIsolation;
  string _RecoTauDiscrByLeadTrack;
  string _RecoTauDiscrAgainstElectron;
  string _RecoTauDiscrAgainstMuon;
  string _RecoTauDiscrByProngType;
  bool _DoRecoTauDiscrByH3x3OverP;
  double _RecoTauH3x3OverP;

  //-----Reco Muon Inputs
  InputTag _RecoMuonSource;
  double _RecoMuonPtMinCut;
  double _RecoMuonPtMaxCut;
  double _RecoMuonEtaCut;
  double _RecoMuonTrackIsoSumPtMinCutValue;
  double _RecoMuonTrackIsoSumPtMaxCutValue;
  double _RecoMuonEcalIsoSumPtMinCutValue;
  double _RecoMuonEcalIsoSumPtMaxCutValue;
  double _RecoMuonIsoDeltaRCone;
  double _RecoMuonTrackIsoTrkThreshold;
  double _RecoMuonEcalIsoRecHitThreshold;
  double _RecoMuonIpCut;
  double _RecoMuonCaloCompCoefficient;
  double _RecoMuonSegmCompCoefficient;
  double _RecoMuonAntiPionCut;
  bool _DoRecoMuonDiscrByGlobal;
  bool _DoRecoMuonDiscrByIsolation;
  bool _DoRecoMuonDiscrByIp;
  bool _DoRecoMuonDiscrByPionVeto;

  //-----Reco Electron Inputs
  InputTag _RecoElectronSource;
  double _RecoElectronPtMinCut;
  double _RecoElectronPtMaxCut;
  double _RecoElectronEtaCut;
  double _RecoElectronTrackIsoSumPtMaxCutValue;
  double _RecoElectronTrackIsoSumPtMinCutValue;
  double _RecoElectronTrackIsoDeltaRCone;
  double _RecoElectronTrackIsoTrkThreshold;
  double _RecoElectronEcalIsoSumPtMaxCutValue;
  double _RecoElectronEcalIsoSumPtMinCutValue;
  double _RecoElectronEcalIsoDeltaRCone;
  double _RecoElectronEcalIsoRecHitThreshold;
  double _RecoElectronIpCut;  
  double _RecoElectronEoverPMax;
  double _RecoElectronEoverPMin;
  double _RecoElectronHoverEmCut;
  bool _UseHeepInfo;
  bool _DoRecoElectronDiscrByTrackIsolation;
  bool _DoRecoElectronDiscrByEcalIsolation;
  bool _DoRecoElectronDiscrByIp;
  bool _DoRecoElectronDiscrByEoverP;
  bool _DoRecoElectronDiscrByHoverEm;
  bool _UseElectronUserIsolation;
  bool _DoRecoElectronDiscrByEcalDrivenSeed;
  bool _DoRecoElectronDiscrByTrackerDrivenSeed;
  bool _DoRecoElectronDiscrBySigmaIEtaIEta;
  double _RecoElectronSigmaIEtaIEta;
  bool _DoRecoElectronDiscrByDEtaIn;
  double _RecoElectronEEDEtaIn;
  double _RecoElectronEBDEtaIn;
  bool _DoRecoElectronDiscrByDPhiIn;
  double _RecoElectronEEDPhiIn;
  double _RecoElectronEBDPhiIn;
  bool _DoRecoElectronDiscrBySCE2by5Over5by5;
  double _RecoElectronEBscE1by5Over5by5;
  double _RecoElectronEBscE2by5Over5by5;
  bool _DoRecoElectronDiscrByMissingHits;
  int _RecoElectronMissingHits;

  //-----Reco ParticleFlow Inputs
  InputTag _RecoParticleFlowSource;
  bool _UsePFlowBasedIsolationInsteadOfStandard;

  //-----Reco Jet Inputs
  InputTag _RecoJetSource;
  double _RecoJetPtCut;
  double _RecoJetEtaMinCut;
  double _RecoJetEtaMaxCut;
  double _JetMuonMatchingDeltaR;
  double _JetElectronMatchingDeltaR;
  double _JetTauMatchingDeltaR;
  double _SecondLeadingJetMuonMatchingDeltaR;
  double _SecondLeadingJetElectronMatchingDeltaR;
  double _SecondLeadingJetTauMatchingDeltaR;
  double _FirstLeadingJetMuonMatchingDeltaR;
  double _FirstLeadingJetElectronMatchingDeltaR;
  double _FirstLeadingJetTauMatchingDeltaR;
  double _JetBTaggingTCHEcut;
  double _RecoFirstLeadingJetPt;
  double _RecoFirstLeadingJetEtaMinCut;
  double _RecoFirstLeadingJetEtaMaxCut;
  double _RecoSecondLeadingJetPt;
  double _RecoSecondLeadingJetEtaMinCut;
  double _RecoSecondLeadingJetEtaMaxCut;
  bool _UseCorrectedJet;
  bool _RemoveJetOverlapWithMuons;
  bool _RemoveJetOverlapWithElectrons;
  bool _RemoveJetOverlapWithTaus;
  bool _RemoveSecondLeadingJetOverlapWithMuons;
  bool _RemoveSecondLeadingJetOverlapWithElectrons;
  bool _RemoveSecondLeadingJetOverlapWithTaus;
  bool _RemoveFirstLeadingJetOverlapWithMuons;
  bool _RemoveFirstLeadingJetOverlapWithElectrons;
  bool _RemoveFirstLeadingJetOverlapWithTaus;
  bool _ApplyJetBTagging;
  bool _DoDiscrByFirstLeadingJet;
  bool _DoDiscrBySecondLeadingJet;

  double sumpxForMht;
  double sumpyForMht;
  double phiForMht;
  double sumptForHt;
  double leadingjetpt;
  double secondleadingjetpt;
  int theLeadingJetIndex;
  int theSecondLeadingJetIndex;

  //-----Vertex Inputs
  InputTag _RecoVertexSource;
  double _RecoVertexMaxZposition;
  double _RecoVertexTrackWeight;
  int _RecoVertexMinTracks;
  bool _DoPUrewieghting;
  string _DataHistos;
  string  _MCHistos;


  //-----Trigger Inputs
  InputTag _RecoTriggerSource;
  std::vector<std::string> _TriggerRequirements;

  //-----Topology Inputs
//  InputTag _RecoDiTauSource;
  InputTag _RecoMetSource;
  double _RecoMetCut;
  double _DiTauDeltaRCut;
  double _DiTauCosDphiMinCut;
  double _DiTauCosDphiMaxCut;
  double _MassMinCut;
  double _MassMaxCut;
  double _PZetaCutCoefficient;
  double _PZetaVisCutCoefficient;
  double _CDFzeta2DCutValue;
  double _DiTauDeltaPtDivSumPtMinCutValue;
  double _DiTauDeltaPtDivSumPtMaxCutValue;
  double _DiTauDeltaPtMinCutValue;
  double _DiTauDeltaPtMaxCutValue;
  double _Leg1MetDphiMinCut;
  double _Leg1MetDphiMaxCut;
  double _Leg2MetDphiMinCut;
  double _Leg2MetDphiMaxCut;
  bool _DoDiscrByMet;
  bool _CalculateMetUsingOnlyLeg1AndLeg2;
  bool _DoDiTauDiscrByDeltaR;
  bool _UseTauSeedTrackForDiTauDiscrByOSLS;
  bool _DoDiTauDiscrByCosDphi;
  bool _DoDiscrByMassReco;
  bool _UseVectorSumOfVisProductsAndMetMassReco;
  bool _UseCollinerApproxMassReco;
  bool _DoDiTauDiscrByCDFzeta2D;
  bool _DoDiTauDiscrByDeltaPtDivSumPt;
  bool _DoDiTauDiscrByDeltaPt;
  bool _DoDiscrByLeg1MetDphi;
  bool _DoDiscrByLeg2MetDphi;
  bool _DoTauDiscrByIsZeeCut;
  string _DiTauDiscrByOSLSType;

  //-----SUSY Specific Topology Inputs
  bool _DoSUSYDiscrByMHT;
  double _MhtCut;
  bool _DoSUSYDiscrByR1;
  double _R1MinCut;
  double _R1MaxCut;
  bool _DoSUSYDiscrByR2;
  double _R2MinCut;
  double _R2MaxCut;
  bool _DoSUSYDiscrByAlpha;
  double _AlphaMinCut;
  double _AlphaMaxCut;
//  bool _DoSUSYDiscrByDphi1;
//  double _Dphi1MinCut;
//  double _Dphi1MaxCut;
//  bool _DoSUSYDiscrByDphi2;
//  double _Dphi2MinCut;
//  double _Dphi2MaxCut;
  bool _DoSUSYDiscrByDphiMhtJet;
  double _DphiMhtJet1;
  double _DphiMhtJet2;

  //-----do matching to gen?
  bool _MatchTauToGen;
  bool _UseTauMotherId;
  bool _UseTauGrandMotherId;
  bool _MatchLeptonToGen;
  bool _UseLeptonMotherId;
  bool _UseLeptonGrandMotherId;
  int _TauMotherId;
  int _TauGrandMotherId;
  int _LeptonMotherId;
  int _LeptonGrandMotherId;
  double _TauToGenMatchingDeltaR;
  std::vector<reco::GenParticleRef> associatedGenParticles;
  reco::Candidate::LorentzVector MChadtau;
  const reco::Candidate * daughterCand;
  const reco::Candidate * motherCand;
  const reco::Candidate * grandMotherCand;

  //-----ntuple Inputs
  TTree *_HMTTree;
  std::string _NtupleTreeName;
  bool _DoProduceNtuple;

  //-----Fill Histograms?
  bool _FillRecoVertexHists;
  bool _FillGenTauHists;
  bool _FillRecoTauHists;
  bool _FillRecoMuonHists;
  bool _FillRecoElectronHists;
  bool _FillRecoJetHists;
  bool _FillTopologyHists;

  //-----histogram that keeps track of the number of analyzed events & the number of
  //-----events passing the user defined cuts
  std::map<unsigned int, TH1*> _hEvents;

  //-----vertex histograms
  std::map<unsigned int, TH1*> _hVertexZposition;
  std::map<unsigned int, TH1*> _hVertexNTracks;
  std::map<unsigned int, TH1*> _hNVertices;

  //-----generator level histograms
  std::map<unsigned int, TH1*> _hNGenTau;
  std::map<unsigned int, TH1*> _hGenTauEnergy;
  std::map<unsigned int, TH1*> _hGenTauPt;
  std::map<unsigned int, TH1*> _hGenTauEta;
  std::map<unsigned int, TH1*> _hGenTauPhi;
  std::map<unsigned int, TH1*> _hGenTauMotherEnergy;
  std::map<unsigned int, TH1*> _hGenTauMotherPt;
  std::map<unsigned int, TH1*> _hGenTauMotherEta;
  std::map<unsigned int, TH1*> _hGenTauMotherPhi;
  std::map<unsigned int, TH1*> _hGenTauGrandMotherEnergy;
  std::map<unsigned int, TH1*> _hGenTauGrandMotherPt;
  std::map<unsigned int, TH1*> _hGenTauGrandMotherEta;
  std::map<unsigned int, TH1*> _hGenTauGrandMotherPhi;

  //-----reconstruction level tau histograms
  std::map<unsigned int, TH1*> _hNTau;
  std::map<unsigned int, TH1*> _hTauJetEnergy;
  std::map<unsigned int, TH1*> _hTauJetPt;
  std::map<unsigned int, TH1*> _hTauJetEta;
  std::map<unsigned int, TH1*> _hBestTauJetPt;
  std::map<unsigned int, TH1*> _hBestTauJetEta;
  std::map<unsigned int, TH1*> _hTauJetPhi;
  std::map<unsigned int, TH1*> _hTauJetNumSignalTracks;
  std::map<unsigned int, TH1*> _hTauJetNumSignalGammas;
  std::map<unsigned int, TH1*> _hTauJetSeedTrackPt;
  std::map<unsigned int, TH1*> _hBestTauJetSeedTrackPt;
  std::map<unsigned int, TH1*> _hTauJetSeedTrackIpSignificance;
  std::map<unsigned int, TH1*> _hTauJetSeedTrackNhits;
  std::map<unsigned int, TH1*> _hBestTauJetSeedTrackNhits;
  std::map<unsigned int, TH1*> _hTauJetSeedTrackChi2;
  std::map<unsigned int, TH1*> _hTauJetCharge;
  std::map<unsigned int, TH1*> _hTauJetSignalTracksMass;
  std::map<unsigned int, TH1*> _hTauJetSignalTracksAndGammasMass;
  std::map<unsigned int, TH1*> _hTauJetSignalTracksMass1prong;
  std::map<unsigned int, TH1*> _hTauJetSignalTracksAndGammasMass1prong;
  std::map<unsigned int, TH1*> _hTauJetSignalTracksAndPiZerosMass1prong;
  std::map<unsigned int, TH1*> _hTauJetNumSignalPiZeros1prong;
  std::map<unsigned int, TH1*> _hTauJetSignalTracksMass3prong;
  std::map<unsigned int, TH1*> _hTauJetSignalTracksAndGammasMass3prong;
  std::map<unsigned int, TH1*> _hTauJetMass1Prong0Gamma;
  std::map<unsigned int, TH1*> _hTauJetMass1Prong1Gamma;
  std::map<unsigned int, TH1*> _hTauJetMass1Prong2orMoreGamma;
  std::map<unsigned int, TH1*> _hTauJetMass3Prong0Gamma;
  std::map<unsigned int, TH1*> _hTauJetMass3Prong1Gamma;
  std::map<unsigned int, TH1*> _hTauJetMass3Prong2orMoreGamma;
  std::map<unsigned int, TH1*> _hTauJetSignalTracksChargeFraction;
  std::map<unsigned int, TH1*> _hTauJetNumIsoTracks;
  std::map<unsigned int, TH1*> _hTauJetNumIsoGammas;
  std::map<unsigned int, TH1*> _hTauJetNumIsoCands;
  std::map<unsigned int, TH1*> _hTauJetSumPtIsoTracks;
  std::map<unsigned int, TH1*> _hTauJetSumPtIsoGammas;
  std::map<unsigned int, TH1*> _hBestTauJetSumPtIsoTracks;
  std::map<unsigned int, TH1*> _hBestTauJetSumPtIsoGammas;
  std::map<unsigned int, TH1*> _hTauJetSumPtIso;
  std::map<unsigned int, TH1*> _hTauJetNumberDensity;
  std::map<unsigned int, TH1*> _hTauJetGenTauDeltaPhi;
  std::map<unsigned int, TH1*> _hTauJetGenTauDeltaEta;
  std::map<unsigned int, TH1*> _hTauJetGenTauDeltaPt;
  std::map<unsigned int, TH1*> _hTauJetH3x3OverP;
  std::map<unsigned int, TH1*> _hBestTauJetH3x3OverP;

  //-----reconstruction level muon histograms
  std::map<unsigned int, TH1*> _hNMuon;
  std::map<unsigned int, TH1*> _hMuonEnergy;
  std::map<unsigned int, TH1*> _hMuonPt;
  std::map<unsigned int, TH1*> _hMuonEta;
  std::map<unsigned int, TH1*> _hBestMuonPt;
  std::map<unsigned int, TH1*> _hBestMuonEta;
  std::map<unsigned int, TH1*> _hMuonPhi;
  std::map<unsigned int, TH1*> _hMuonTrackIso;
  std::map<unsigned int, TH1*> _hMuonEcalIso;
  std::map<unsigned int, TH1*> _hBestMuonTrackIso;
  std::map<unsigned int, TH1*> _hBestMuonEcalIso;
  std::map<unsigned int, TH1*> _hMuonIp;
  std::map<unsigned int, TH1*> _hBestMuonIp;
  std::map<unsigned int, TH1*> _hMuonIpSignificance;
  std::map<unsigned int, TH1*> _hMuonIso;
  std::map<unsigned int, TH1*> _hMuonGenMuonDeltaPhi;
  std::map<unsigned int, TH1*> _hMuonGenMuonDeltaEta;
  std::map<unsigned int, TH1*> _hMuonGenMuonDeltaPt;
  std::map<unsigned int, TH1*> _hMuonCaloCompatibility;
  std::map<unsigned int, TH1*> _hMuonSegmentCompatibility;
  std::map<unsigned int, TH1*> _hBestMuonCaloCompatibility;
  std::map<unsigned int, TH1*> _hBestMuonSegmentCompatibility;
  std::map<unsigned int, TH1*> _hMuonAntiPion;
  std::map<unsigned int, TH1*> _hBestMuonAntiPion;
  std::map<unsigned int, TH2*> _hMuonCaloCompatibilityVsSegmentCompatibility;

  //-----reconstruction level electron histograms  
  std::map<unsigned int, TH1*> _hNElectron;
  std::map<unsigned int, TH1*> _hElectronEnergy;
  std::map<unsigned int, TH1*> _hElectronPt;
  std::map<unsigned int, TH1*> _hElectronEta;
  std::map<unsigned int, TH1*> _hBestElectronPt;
  std::map<unsigned int, TH1*> _hBestElectronEta;
  std::map<unsigned int, TH1*> _hElectronPhi;
  std::map<unsigned int, TH1*> _hElectronTrackIso;
  std::map<unsigned int, TH1*> _hElectronEcalIso;
  std::map<unsigned int, TH1*> _hBestElectronTrackIso;
  std::map<unsigned int, TH1*> _hBestElectronEcalIso;
  std::map<unsigned int, TH1*> _hElectronIp;
  std::map<unsigned int, TH1*> _hElectronEoverP;
  std::map<unsigned int, TH1*> _hElectronClassification;
  std::map<unsigned int, TH1*> _hElectronGenElectronDeltaPhi;
  std::map<unsigned int, TH1*> _hElectronGenElectronDeltaEta;
  std::map<unsigned int, TH1*> _hElectronGenElectronDeltaPt;
  std::map<unsigned int, TH1*> _hElectronEcalDriven;
  std::map<unsigned int, TH1*> _hElectronTrackerDriven;
  std::map<unsigned int, TH1*> _hElectronHoverEm;
  std::map<unsigned int, TH1*> _hElectronEESigmaIEtaIEta;
  std::map<unsigned int, TH1*> _hElectronEEDEta;
  std::map<unsigned int, TH1*> _hElectronEEDPhi;
  std::map<unsigned int, TH1*> _hElectronEBSigmaIEtaIEta;
  std::map<unsigned int, TH1*> _hElectronEBDEta;
  std::map<unsigned int, TH1*> _hElectronEBDPhi;
  std::map<unsigned int, TH1*> _hElectronEB2by5Over5by5;
  std::map<unsigned int, TH1*> _hElectronEB1by5Over5by5;
  std::map<unsigned int, TH1*> _hElectronMissingHits;
  std::map<unsigned int, TH1*> _hBestElectronHoverEm;
  std::map<unsigned int, TH1*> _hBestElectronEESigmaIEtaIEta;
  std::map<unsigned int, TH1*> _hBestElectronEEDEta;
  std::map<unsigned int, TH1*> _hBestElectronEEDPhi;
  std::map<unsigned int, TH1*> _hBestElectronEBSigmaIEtaIEta;
  std::map<unsigned int, TH1*> _hBestElectronEBDEta;
  std::map<unsigned int, TH1*> _hBestElectronEBDPhi;
  std::map<unsigned int, TH1*> _hBestElectronEB2by5Over5by5;
  std::map<unsigned int, TH1*> _hBestElectronEB1by5Over5by5;
  std::map<unsigned int, TH1*> _hBestElectronMissingHits;

  //-----reconstruction level jet histograms  
  std::map<unsigned int, TH1*> _hNJet;
  std::map<unsigned int, TH1*> _hJetEnergy;
  std::map<unsigned int, TH1*> _hJetPt;
  std::map<unsigned int, TH1*> _hJetEta;
  std::map<unsigned int, TH1*> _hJetPhi;
  std::map<unsigned int, TH1*> _hBJetDiscrByTrackCounting;
  std::map<unsigned int, TH1*> _hBJetDiscrBySimpleSecondaryV;
  std::map<unsigned int, TH1*> _hBJetDiscrByCombinedSecondaryV;
  std::map<unsigned int, TH1*> _hFirstLeadingJetPt;
  std::map<unsigned int, TH1*> _hSecondLeadingJetPt;
  std::map<unsigned int, TH1*> _hMHT;
  std::map<unsigned int, TH1*> _hHT;

  //-----reconstruction level topology histograms
  std::map<unsigned int, TH2*> _hMuonPtVsTauPt;
  std::map<unsigned int, TH1*> _hMuonTauDeltaR;
  std::map<unsigned int, TH1*> _hBestMuonTauDeltaR;
  std::map<unsigned int, TH1*> _hMuonTauDeltaPtDivSumPt;
  std::map<unsigned int, TH1*> _hMuonTauDeltaPt;
  std::map<unsigned int, TH1*> _hMuonTauOSLS;
  std::map<unsigned int, TH1*> _hBestMuonTauOSLS;
  std::map<unsigned int, TH1*> _hMuonTauCosDphi;
  std::map<unsigned int, TH1*> _hBestMuonTauCosDphi;
  std::map<unsigned int, TH1*> _hMuonMetDeltaPhi;
  std::map<unsigned int, TH1*> _hTauMetDeltaPhi;
  std::map<unsigned int, TH2*> _hMuonMetDeltaPhiVsMuonTauCosDphi;
  std::map<unsigned int, TH2*> _hElectronPtVsTauPt;
  std::map<unsigned int, TH1*> _hElectronTauDeltaR;
  std::map<unsigned int, TH1*> _hBestElectronTauDeltaR;
  std::map<unsigned int, TH1*> _hElectronTauDeltaPtDivSumPt;
  std::map<unsigned int, TH1*> _hElectronTauDeltaPt;
  std::map<unsigned int, TH1*> _hElectronTauOSLS;
  std::map<unsigned int, TH1*> _hBestElectronTauOSLS;
  std::map<unsigned int, TH1*> _hElectronTauCosDphi;
  std::map<unsigned int, TH1*> _hBestElectronTauCosDphi;
  std::map<unsigned int, TH1*> _hElectronMetDeltaPhi;
  std::map<unsigned int, TH2*> _hElectronMetDeltaPhiVsElectronTauCosDphi;
  std::map<unsigned int, TH2*> _hElectronPtVsMuonPt;
  std::map<unsigned int, TH1*> _hElectronMuonDeltaR;
  std::map<unsigned int, TH1*> _hElectronMuonDeltaPtDivSumPt;
  std::map<unsigned int, TH1*> _hElectronMuonDeltaPt;
  std::map<unsigned int, TH1*> _hElectronMuonOSLS;
  std::map<unsigned int, TH1*> _hElectronMuonCosDphi;
  std::map<unsigned int, TH2*> _hElectronMetDeltaPhiVsElectronMuonCosDphi;
  std::map<unsigned int, TH2*> _hTau1PtVsTau2Pt;
  std::map<unsigned int, TH1*> _hTau1Tau2DeltaR;
  std::map<unsigned int, TH1*> _hTau1Tau2DeltaPtDivSumPt;
  std::map<unsigned int, TH1*> _hTau1Tau2DeltaPt;
  std::map<unsigned int, TH1*> _hTau1MetMt;
  std::map<unsigned int, TH1*> _hTau2MetMt;
  std::map<unsigned int, TH1*> _hTau1Tau2OSLS;
  std::map<unsigned int, TH1*> _hTau1Tau2CosDphi;
  std::map<unsigned int, TH1*> _hTau1MetDeltaPhi;
  std::map<unsigned int, TH1*> _hTau2MetDeltaPhi;
  std::map<unsigned int, TH2*> _hTau1MetDeltaPhiVsTau1Tau2CosDphi;
  std::map<unsigned int, TH2*> _hMuon1PtVsMuon2Pt;
  std::map<unsigned int, TH1*> _hMuon1Muon2DeltaR;
  std::map<unsigned int, TH1*> _hMuon1Muon2DeltaPtDivSumPt;
  std::map<unsigned int, TH1*> _hMuon1Muon2DeltaPt;
  std::map<unsigned int, TH1*> _hMuon1MetMt;
  std::map<unsigned int, TH1*> _hMuon2MetMt;
  std::map<unsigned int, TH1*> _hMuon1Muon2OSLS;
  std::map<unsigned int, TH1*> _hMuon1Muon2CosDphi;
  std::map<unsigned int, TH1*> _hMuon1MetDeltaPhi;
  std::map<unsigned int, TH1*> _hMuon2MetDeltaPhi;
  std::map<unsigned int, TH2*> _hMuon1MetDeltaPhiVsMuon1Muon2CosDphi;
  std::map<unsigned int, TH2*> _hElectron1PtVsElectron2Pt;
  std::map<unsigned int, TH1*> _hElectron1Electron2DeltaR;
  std::map<unsigned int, TH1*> _hElectron1Electron2DeltaPtDivSumPt;
  std::map<unsigned int, TH1*> _hElectron1Electron2DeltaPt;
  std::map<unsigned int, TH1*> _hElectron1MetMt;
  std::map<unsigned int, TH1*> _hElectron2MetMt;
  std::map<unsigned int, TH1*> _hElectron1Electron2OSLS;
  std::map<unsigned int, TH1*> _hElectron1Electron2CosDphi;
  std::map<unsigned int, TH1*> _hElectron1MetDeltaPhi;
  std::map<unsigned int, TH1*> _hElectron2MetDeltaPhi;
  std::map<unsigned int, TH2*> _hElectron1MetDeltaPhiVsElectron1Electron2CosDphi;
  std::map<unsigned int, TH1*> _hTauMetMt;
  std::map<unsigned int, TH1*> _hElectronMetMt;
  std::map<unsigned int, TH1*> _hBestElectronMetMt;
  std::map<unsigned int, TH1*> _hMuonMetMt;
  std::map<unsigned int, TH1*> _hBestMuonMetMt;
  std::map<unsigned int, TH1*> _hNotReconstructableMass;
  std::map<unsigned int, TH1*> _hReconstructableMass;
  std::map<unsigned int, TH1*> _hNotReconstructableMassOS;
  std::map<unsigned int, TH1*> _hReconstructableMassOS;
  std::map<unsigned int, TH1*> _hNotReconstructableMassLS;
  std::map<unsigned int, TH1*> _hReconstructableMassLS;
  std::map<unsigned int, TH1*> _hPZeta;
  std::map<unsigned int, TH1*> _hPZetaVis;
  std::map<unsigned int, TH2*> _hZeta2D;
  std::map<unsigned int, TH1*> _hZeta1D;
  std::map<unsigned int, TH1*> _hBestZeta1D;
  std::map<unsigned int, TH1*> _hMet;
  std::map<unsigned int, TH1*> _hMetResolution;
  std::map<unsigned int, TH1*> _hElectronIsZee;

  std::map<unsigned int, TH1*> _hR1;
  std::map<unsigned int, TH1*> _hR2;
  std::map<unsigned int, TH1*> _hDphi1;
  std::map<unsigned int, TH1*> _hDphiMhtJet1;
  std::map<unsigned int, TH1*> _hDphi2;
  std::map<unsigned int, TH1*> _hDphiMhtJet2;
  std::map<unsigned int, TH1*> _hAlpha;

/* 
  std::map<unsigned int, TH1*> _hTauJetSumPtIso_JetOS;
  std::map<unsigned int, TH1*> _hTauJetSumPtIso_JetLS;
  std::map<unsigned int, TH1*> _hTauJetSumPtIso_SeedOS;
  std::map<unsigned int, TH1*> _hTauJetSumPtIso_SeedLS;
*/

  //-----Handle to tree collections
  Handle< reco::GenParticleCollection > _genParticles;
  Handle< pat::TauCollection >    _patTaus;
  Handle< pat::MuonCollection >    _patMuons;
  Handle< pat::ElectronCollection > _patElectrons;
  Handle< pat::JetCollection > _patJets;
  Handle< pat::METCollection > _patMETs;
//  Handle< reco::CompositeCandidateCollection > _patDiTaus;
  Handle< reco::VertexCollection > _primaryEventVertexCollection;
  edm::Handle< edm::TriggerResults > _triggerResults;
  edm::Handle< trigger::TriggerEvent > handleTriggerEvent;
  Handle< reco::PFCandidateCollection > _pflow;

  Handle< reco::PFTauCollection > _hpsTau;
  Handle< reco::GenMETCollection > genTrue;

/*
  Handle< edm::View<reco::Muon> >   _recoMuonsForMetCorrections;
  Handle< ValueMap<reco::MuonMETCorrectionData > > vm_muCorrData_h;
  Handle< edm::View<pat::Muon> >    _patMuonsForMetCorrections;
  const edm::View<reco::Muon>& muonsrecoForMetCorrections;
  const edm::ValueMap<reco::MuonMETCorrectionData>& vm_muCorrData;
  const edm::View<pat::Muon>& muonsForMetCorrections;
*/

  //-----Event Sequence inputTags
  int _RecoTriggersNmin;
  int _RecoVertexNmin;
  int _RecoVertexNmax;
  int _RecoLeg1Nmin;
  int _RecoLeg1Nmax;
  int _RecoLeg2Nmin;
  int _RecoLeg2Nmax;
  int _RecoJetNmin;
  int _RecoJetNmax;
  int _RecoFirstLeadingJetNmin;
  int _RecoSecondLeadingJetNmin;
  int _SusyCombinationsNmin;
  int _CombinationsNmin;
  int _CombinationsNmax;
  vector<string> _EventSelectionSequence;

  //-----Event flags
  vector<bool> _EventFlag;

  //-----Event counters
  int _totalEvents;
  int _totalEventsPassingCuts;

  //-----Mapping of cut names and cut id
  std::map<string,int> _mapSelectionAlgoID;
  std::map<string,int> _mapSelectionCounter;
  std::map<string,int> _mapSelectionCounterCumul;

  //-----Inputs for systematic uncertainties
  bool _CalculatePdfSystematicUncertanties;
  bool _CalculateFSRSystematics;
  bool _CalculateISRGluonSystematics;
  bool _CalculateISRGammaSystematics;
  std::vector<edm::InputTag> pdfWeightTags_;
  InputTag FSRWeightTag_;
  InputTag ISRWeightTag_;
  std::vector<int> pdfStart_Denominator_;
  std::vector<double> weightedEvents_Denominator_;
  std::vector<int> pdfStart_Numerator_;
  std::vector<double> weightedEvents_Numerator_;
  double isrgluon_weight;
  double isrgamma_weight;
  double fsr_weight;
  double pu_weight;
  double trig_weight;
  double pdfcentral_weight;
  bool _SmearTheMuon;
  double _MuonPtScaleOffset;
  double _MuonPtSigmaOffset;
  double _MuonEtaScaleOffset;
  double _MuonEtaSigmaOffset;
  double _MuonPhiScaleOffset;
  double _MuonPhiSigmaOffset;
  bool _SmearTheElectron;
  double _ElectronPtScaleOffset;
  double _ElectronPtSigmaOffset;
  double _ElectronEtaScaleOffset;
  double _ElectronEtaSigmaOffset;
  double _ElectronPhiScaleOffset;
  double _ElectronPhiSigmaOffset;
  bool _SmearTheTau;
  double _TauPtScaleOffset;
  double _TauPtSigmaOffset;
  double _TauEtaScaleOffset;
  double _TauEtaSigmaOffset;
  double _TauPhiScaleOffset;
  double _TauPhiSigmaOffset;
  bool _SmearTheJet;
  double _JetEnergyScaleOffset;
  bool _SmearThePt;
  bool _SmearTheEta;
  bool _SmearThePhi;
  bool _CalculatePUSystematics;
  double _PUConstantWeightFactor;
  double _PUConstantWeightFactorRMS;
  double _PUSlopeWeightParameter;
  double _PUSlopeWeightParameterRMS;

  bool _ApplyMuonTriggerScaleFactors;
  bool _ApplyElectronTriggerScaleFactors;
  bool _ApplyTauTriggerScaleFactors;


  std::vector<reco::Candidate::LorentzVector> smearedMuonMomentumVector;
  std::vector<math::PtEtaPhiMLorentzVector> smearedMuonPtEtaPhiMVector;
  std::vector<reco::Candidate::LorentzVector> smearedElectronMomentumVector;
  std::vector<math::PtEtaPhiMLorentzVector> smearedElectronPtEtaPhiMVector;
  std::vector<reco::Candidate::LorentzVector> smearedTauMomentumVector;
  std::vector<math::PtEtaPhiMLorentzVector> smearedTauPtEtaPhiMVector;
  std::vector<reco::Candidate::LorentzVector> smearedJetMomentumVector;
  reco::Candidate::LorentzVector maxPtMuonVector;
  reco::Candidate::LorentzVector maxEtElectronVector;
  std::vector<math::PtEtaPhiMLorentzVector> smearedJetPtEtaPhiMVector;
  std::vector<double> bosonPtBinEdges_;
  std::vector<double> ptWeights_;
  std::vector<double> MuonTrigPtBinEdges_;
  std::vector<double> MuonTrigptWeights_;
  std::vector<double> MuonTrigEtaBinEdges_;
  std::vector<double> MuonTrigetaWeights_;
  std::vector<double> ElectronTrigPtBinEdges_;
  std::vector<double> ElectronTrigptWeights_;
  std::vector<double> ElectronTrigEtaBinEdges_;
  std::vector<double> ElectronTrigetaWeights_;
  std::vector<double> TauTrigPtBinEdges_;
  std::vector<double> TauTrigptWeights_;
  std::vector<double> TauTrigEtaBinEdges_;
  std::vector<double> TauTrigetaWeights_;

  double deltaForMEx;
  double deltaForMEy;
  reco::Candidate::LorentzVector theMETVector;
  
  //-----For Ntuple
  vector<double> pdfWeightVector;
  vector<float> *_jetPt;
  vector<float> *_jetEt;
  vector<float> *_jetEta;
  vector<float> *_jetPhi;
  vector<float> *_jetEmFraction;
  vector<float> *_bJetDiscrByTrackCounting;
  vector<float> *_bJetDiscrBySimpleSecondaryV;
  vector<float> *_bJetDiscrByCombinedSecondaryV;
  vector<float> *_tauTrkIsoPat;
  vector<float> *_muTrkIsoPat;
  vector<float> *_mEt;
  vector<float> *_zeta;
  vector<float> *_muTauMetMass;
  vector<float> *_muMetMass;
  vector<float> *_OSLS;

};

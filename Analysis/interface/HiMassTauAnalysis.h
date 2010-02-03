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
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "CLHEP/Random/RandGauss.h"

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <string>
#include <vector>
#include <map>

using namespace std;
using namespace edm;

typedef math::XYZTLorentzVector LorentzVector;
typedef std::vector< math::XYZTLorentzVector > LVCollection;

class HiMassTauAnalysis : public EDAnalyzer {
public:
  explicit HiMassTauAnalysis(const ParameterSet&);
  ~HiMassTauAnalysis();


private:
  virtual void beginJob(const EventSetup&) ;
  virtual void analyze(const Event&, const EventSetup&);
  virtual void endJob() ;

  void getCollections(const Event&, const EventSetup&);
  void fillHistograms();
  void fillNtuple();
  void getEventFlags();
  bool passEventSelectionSequence();
  //bool passGenTauCuts(const LorentzVector&);
  bool passRecoTriggerCuts();
  bool passRecoVertexCuts(const reco::Vertex&);
  bool passRecoTauCuts(const pat::Tau&);
  bool passRecoMuonCuts(const pat::Muon&,bool,reco::Candidate::LorentzVector);
  bool passRecoElectronCuts(const pat::Electron&,bool,reco::Candidate::LorentzVector);
  bool passRecoJetCuts(const pat::Jet&);
  bool passTopologyCuts(const pat::Tau&, const pat::Muon&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  bool passTopologyCuts(const pat::Tau&, const pat::Electron&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  bool passTopologyCuts(const pat::Electron&, bool,reco::Candidate::LorentzVector,const pat::Muon&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  bool passTopologyCuts(const pat::Muon&, bool,reco::Candidate::LorentzVector,const pat::Muon&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  bool passTopologyCuts(const pat::Tau&, const pat::Tau&, const pat::MET&);
  bool passTopologyCuts(const pat::Electron&, bool,reco::Candidate::LorentzVector,const pat::Electron&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  pair<bool, math::XYZTLorentzVector> matchToGen(const pat::Electron&);
  pair<bool, math::XYZTLorentzVector> matchToGen(const pat::Muon&);
  pair<bool, math::XYZTLorentzVector> matchToGen(const pat::Tau&);
  double CalculatePZeta(const pat::Tau&, const pat::Muon&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  double CalculatePZeta(const pat::Tau&, const pat::Electron&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  double CalculatePZeta(const pat::Electron&, bool,reco::Candidate::LorentzVector,const pat::Muon&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  double CalculatePZeta(const pat::Tau&, const pat::Tau&, const pat::MET&);
  double CalculatePZeta(const pat::Muon&, bool,reco::Candidate::LorentzVector,const pat::Muon&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  double CalculatePZeta(const pat::Electron&, bool,reco::Candidate::LorentzVector,const pat::Electron&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  double CalculatePZetaVis(const pat::Tau&, const pat::Muon&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  double CalculatePZetaVis(const pat::Tau&, const pat::Electron&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  double CalculatePZetaVis(const pat::Electron&, bool,reco::Candidate::LorentzVector,const pat::Muon&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  double CalculatePZetaVis(const pat::Tau&, const pat::Tau&, const pat::MET&);
  double CalculatePZetaVis(const pat::Muon&, bool,reco::Candidate::LorentzVector,const pat::Muon&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  double CalculatePZetaVis(const pat::Electron&, bool,reco::Candidate::LorentzVector,const pat::Electron&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  std::pair<bool, reco::Candidate::LorentzVector> CalculateThe4Momentum(const pat::Tau&, const pat::Electron&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  std::pair<bool, reco::Candidate::LorentzVector> CalculateThe4Momentum(const pat::Tau&, const pat::Muon&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  std::pair<bool, reco::Candidate::LorentzVector> CalculateThe4Momentum(const pat::Electron&, bool,reco::Candidate::LorentzVector,const pat::Muon&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  std::pair<bool, reco::Candidate::LorentzVector> CalculateThe4Momentum(const pat::Tau&, const pat::Tau&, const pat::MET&);
  std::pair<bool, reco::Candidate::LorentzVector> CalculateThe4Momentum(const pat::Muon&, bool,reco::Candidate::LorentzVector,const pat::Muon&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  std::pair<bool, reco::Candidate::LorentzVector> CalculateThe4Momentum(const pat::Electron&, bool,reco::Candidate::LorentzVector,const pat::Electron&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  double CalculateLeptonMetMt(const pat::Muon&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  double CalculateLeptonMetMt(const pat::Electron&, bool,reco::Candidate::LorentzVector,const pat::MET&);
  double CalculateLeptonMetMt(const pat::Tau&, const pat::MET&);
  std::pair<int, double> CalculateTauTrackIsolation(const pat::Tau&);
  std::pair<int, double> CalculateTauTrackIsolation(const pat::Tau&, float, float);
  std::pair<int, double> CalculateTauEcalIsolation(const pat::Tau&);
  std::pair<int, double> CalculateTauEcalIsolation(const pat::Tau&, float, float);
  int CalculateNumberSignalTauGammas(const pat::Tau&);
  reco::Candidate::LorentzVector CalculateTauSignalTracksMass(const pat::Tau&);
  reco::Candidate::LorentzVector SmearLightLepton(const pat::Muon&);
  reco::Candidate::LorentzVector SmearLightLepton(const pat::Electron&);
  std::pair<bool, std::pair<float, float> > isZee(const pat::Electron&);
  void InitializeInfoForPDFSystematicUncertaintites();
  void bookHistograms();
  void setMapSelectionAlgoIDs();
  void initMapSelectionCounters();
  void printEfficiency();
  void setupBranches();
  std::pair<const pat::Tau*, const pat::Electron*> getPATComponents(const reco::CompositeCandidate&);
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
  double _RecoTauIsoSumPtMinCutValue;
  double _RecoTauIsoSumPtMaxCutValue;
  double _RecoTauEcalIsoRphiForEllipse;
  double _RecoTauEcalIsoRetaForEllipse;
  int _RecoTauNisoMax;
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
  string _RecoTauDiscrByIsolation;
  string _RecoTauDiscrByLeadTrack;
  string _RecoTauDiscrAgainstElectron;
  string _RecoTauDiscrAgainstMuon;
  string _RecoTauDiscrByProngType;

  //-----Reco Muon Inputs
  InputTag _RecoMuonSource;
  double _RecoMuonPtMinCut;
  double _RecoMuonPtMaxCut;
  double _RecoMuonEtaCut;
  double _RecoMuonIsoSumPtMinCutValue;
  double _RecoMuonIsoSumPtMaxCutValue;
  double _RecoMuonIsoDeltaRCone;
  double _RecoMuonTrackIsoTrkThreshold;
  double _RecoMuonEcalIsoRecHitThreshold;
  double _RecoMuonIpCut;
  bool _DoRecoMuonDiscrByGlobal;
  bool _DoRecoMuonDiscrByIsolation;
  bool _DoRecoMuonDiscrByIp;

  //-----Reco Electron Inputs
  InputTag _RecoElectronSource;
  double _RecoElectronPtMinCut;
  double _RecoElectronPtMaxCut;
  double _RecoElectronEtaCut;
  double _RecoElectronTrackIsoSumPtCutValue;
  double _RecoElectronTrackIsoDeltaRCone;
  double _RecoElectronTrackIsoTrkThreshold;
  double _RecoElectronEcalIsoSumPtCutValue;
  double _RecoElectronEcalIsoDeltaRCone;
  double _RecoElectronEcalIsoRecHitThreshold;
  double _RecoElectronIpCut;  
  double _RecoElectronEoverPMax;
  double _RecoElectronEoverPMin;
  double _RecoElectronHoverEmCut;
  double _RecoElectronSigmaEtaEtaCut;
  double _RecoElectronSigmaIEtaIEtaCut;
  double _RecoElectronSCE5by5Cut;
  bool _DoRecoElectronDiscrByTrackIsolation;
  bool _DoRecoElectronDiscrByEcalIsolation;
  bool _DoRecoElectronDiscrByIp;
  bool _DoRecoElectronDiscrByEoverP;
  bool _DoRecoElectronDiscrByHoverEm;
  bool _DoRecoElectronDiscrBySigmaEtaEta;
  bool _DoRecoElectronDiscrBySigmaIEtaIEta;
  bool _DoRecoElectronDiscrBySCE5by5;

  //-----Reco Jet Inputs
  InputTag _RecoJetSource;
  double _RecoJetPtCut;
  double _RecoJetEtaMinCut;
  double _RecoJetEtaMaxCut;
  double _JetMuonMatchingDeltaR;
  double _JetElectronMatchingDeltaR;
  double _JetTauMatchingDeltaR;
  bool _UseCorrectedJet;
  bool _RemoveJetOverlapWithMuons;
  bool _RemoveJetOverlapWithElectrons;
  bool _RemoveJetOverlapWithTaus;

  //-----Vertex Inputs
  InputTag _RecoVertexSource;
  double _RecoVertexMaxZposition;
  double _RecoVertexTrackWeight;
  int _RecoVertexMinTracks;

  //-----Trigger Inputs
  InputTag _RecoTriggerSource;
  std::vector<std::string> _TriggerRequirements;

  //-----Topology Inputs
  InputTag _RecoDiTauSource;
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
  double _Leg1MetDphiMinCut;
  double _Leg1MetDphiMaxCut;
  double _Leg2MetDphiMinCut;
  double _Leg2MetDphiMaxCut;
  bool _DoDiscrByMet;
  bool _DoDiTauDiscrByDeltaR;
  bool _UseTauSeedTrackForDiTauDiscrByOSLS;
  bool _DoDiTauDiscrByCosDphi;
  bool _DoDiscrByMassReco;
  bool _UseVectorSumOfVisProductsAndMetMassReco;
  bool _UseCollinerApproxMassReco;
  bool _DoDiTauDiscrByCDFzeta2D;
  bool _DoDiTauDiscrByDeltaPtDivSumPt;
  bool _DoDiscrByLeg1MetDphi;
  bool _DoDiscrByLeg2MetDphi;
  string _DiTauDiscrByOSLSType;

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
  math::XYZTLorentzVector MChadtau;
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
  TH1* _hEvents;

  //-----vertex histograms
  TH1* _hVertexZposition;
  TH1* _hVertexNTracks;
  TH1* _hNVertices;

  //-----generator level histograms
  TH1* _hNGenTau;
  TH1* _hGenTauEnergy;
  TH1* _hGenTauPt;
  TH1* _hGenTauEta;
  TH1* _hGenTauPhi;

  //-----reconstruction level tau histograms
  TH1* _hNTau;
  TH1* _hTauJetEnergy;
  TH1* _hTauJetPt;
  TH1* _hTauJetEta;
  TH1* _hTauJetPhi;
  TH1* _hTauJetNumSignalTracks;
  TH1* _hTauJetNumSignalGammas;
  TH1* _hTauJetSeedTrackPt;
  TH1* _hTauJetSeedTrackIpSignificance;
  TH1* _hTauJetSeedTrackNhits;
  TH1* _hTauJetSeedTrackChi2;
  TH1* _hTauJetCharge;
  TH1* _hTauJetSignalTracksMass;
  TH1* _hTauJetSignalTracksChargeFraction;
  TH1* _hTauJetNumIsoTracks;
  TH1* _hTauJetNumIsoGammas;
  TH1* _hTauJetNumIsoCands;
  TH1* _hTauJetSumPtIsoTracks;
  TH1* _hTauJetSumPtIsoGammas;
  TH1* _hTauJetSumPtIso;
  TH1* _hTauJetGenTauDeltaPhi;
  TH1* _hTauJetGenTauDeltaEta;
  TH1* _hTauJetGenTauDeltaPt;

  //-----reconstruction level muon histograms
  TH1* _hNMuon;
  TH1* _hMuonEnergy;
  TH1* _hMuonPt;
  TH1* _hMuonEta;
  TH1* _hMuonPhi;
  TH1* _hMuonTrackIso;
  TH1* _hMuonEcalIso;
  TH1* _hMuonIp;
  TH1* _hMuonIpSignificance;
  TH1* _hMuonIso;
  TH1* _hMuonGenMuonDeltaPhi;
  TH1* _hMuonGenMuonDeltaEta;
  TH1* _hMuonGenMuonDeltaPt;

  //-----reconstruction level electron histograms  
  TH1* _hNElectron;
  TH1* _hElectronEnergy;
  TH1* _hElectronPt;
  TH1* _hElectronEta;
  TH1* _hElectronPhi;
  TH1* _hElectronTrackIso;
  TH1* _hElectronEcalIso;
  TH1* _hElectronIp;
  TH1* _hElectronEoverP;
  TH1* _hElectronHoverEm;
  TH1* _hElectronClassification;
  TH1* _hElectronGenElectronDeltaPhi;
  TH1* _hElectronGenElectronDeltaEta;
  TH1* _hElectronGenElectronDeltaPt;

  //-----reconstruction level jet histograms  
  TH1* _hNJet;
  TH1* _hJetEnergy;
  TH1* _hJetPt;
  TH1* _hJetEta;
  TH1* _hJetPhi;
  TH1* _hBJetDiscrByTrackCounting;
  TH1* _hBJetDiscrBySimpleSecondaryV;
  TH1* _hBJetDiscrByCombinedSecondaryV;

  //-----reconstruction level topology histograms
  TH2* _hMuonPtVsTauPt;
  TH1* _hMuonTauDeltaR;
  TH1* _hMuonTauDeltaPtDivSumPt;
  TH1* _hMuonTauOSLS;
  TH1* _hMuonTauCosDphi;
  TH1* _hMuonMetDeltaPhi;
  TH1* _hTauMetDeltaPhi;
  TH2* _hMuonMetDeltaPhiVsMuonTauCosDphi;
  TH2* _hElectronPtVsTauPt;
  TH1* _hElectronTauDeltaR;
  TH1* _hElectronTauDeltaPtDivSumPt;
  TH1* _hElectronTauOSLS;
  TH1* _hElectronTauCosDphi;
  TH1* _hElectronMetDeltaPhi;
  TH2* _hElectronMetDeltaPhiVsElectronTauCosDphi;
  TH2* _hElectronPtVsMuonPt;
  TH1* _hElectronMuonDeltaR;
  TH1* _hElectronMuonDeltaPtDivSumPt;
  TH1* _hElectronMuonOSLS;
  TH1* _hElectronMuonCosDphi;
  TH2* _hElectronMetDeltaPhiVsElectronMuonCosDphi;
  TH2* _hTau1PtVsTau2Pt;
  TH1* _hTau1Tau2DeltaR;
  TH1* _hTau1Tau2DeltaPtDivSumPt;
  TH1* _hTau1MetMt;
  TH1* _hTau2MetMt;
  TH1* _hTau1Tau2OSLS;
  TH1* _hTau1Tau2CosDphi;
  TH1* _hTau1MetDeltaPhi;
  TH1* _hTau2MetDeltaPhi;
  TH2* _hTau1MetDeltaPhiVsTau1Tau2CosDphi;
  TH2* _hMuon1PtVsMuon2Pt;
  TH1* _hMuon1Muon2DeltaR;
  TH1* _hMuon1Muon2DeltaPtDivSumPt;
  TH1* _hMuon1MetMt;
  TH1* _hMuon2MetMt;
  TH1* _hMuon1Muon2OSLS;
  TH1* _hMuon1Muon2CosDphi;
  TH1* _hMuon1MetDeltaPhi;
  TH1* _hMuon2MetDeltaPhi;
  TH2* _hMuon1MetDeltaPhiVsMuon1Muon2CosDphi;
  TH2* _hElectron1PtVsElectron2Pt;
  TH1* _hElectron1Electron2DeltaR;
  TH1* _hElectron1Electron2DeltaPtDivSumPt;
  TH1* _hElectron1MetMt;
  TH1* _hElectron2MetMt;
  TH1* _hElectron1Electron2OSLS;
  TH1* _hElectron1Electron2CosDphi;
  TH1* _hElectron1MetDeltaPhi;
  TH1* _hElectron2MetDeltaPhi;
  TH2* _hElectron1MetDeltaPhiVsElectron1Electron2CosDphi;
  TH1* _hTauMetMt;
  TH1* _hElectronMetMt;
  TH1* _hMuonMetMt;
  TH1* _hNotReconstructableMass;
  TH1* _hReconstructableMass;
  TH1* _hPZeta;
  TH1* _hPZetaVis;
  TH2* _hZeta2D;
  TH1* _hZeta1D;
  TH1* _hMet;

/* 
  TH1* _hTauJetSumPtIso_JetOS;
  TH1* _hTauJetSumPtIso_JetLS;
  TH1* _hTauJetSumPtIso_SeedOS;
  TH1* _hTauJetSumPtIso_SeedLS;
*/

  //-----Handle to tree collections
  Handle< reco::GenParticleCollection > _genParticles;
  Handle< pat::TauCollection >    _patTaus;
  Handle< pat::MuonCollection >    _patMuons;
  Handle< pat::ElectronCollection > _patElectrons;
  Handle< pat::JetCollection > _patJets;
  Handle< pat::METCollection > _patMETs;
  Handle< reco::CompositeCandidateCollection > _patDiTaus;
  Handle< reco::VertexCollection > _primaryEventVertexCollection;
  edm::Handle< edm::TriggerResults > _triggerResults;
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
  std::vector<edm::InputTag> pdfWeightTags_;
  std::vector<int> pdfStart_Denominator_;
  std::vector<double> weightedEvents_Denominator_;
  std::vector<int> pdfStart_Numerator_;
  std::vector<double> weightedEvents_Numerator_;
  bool _SmearTheMuon;
  double _RelativeMuonPtOffset;
  double _RelativeMuonPtSigma;
  double _AbsoluteMuonEtaOffset;
  double _AbsoluteMuonEtaSigma;
  double _AbsoluteMuonPhiOffset;
  double _AbsoluteMuonPhiSigma;
  bool _SmearTheElectron;
  double _RelativeElectronPtOffset;
  double _RelativeElectronPtSigma;
  double _AbsoluteElectronEtaOffset;
  double _AbsoluteElectronEtaSigma;
  double _AbsoluteElectronPhiOffset;
  double _AbsoluteElectronPhiSigma;
  std::vector<reco::Candidate::LorentzVector> smearedMuonMomentumVector;
  std::vector<reco::Candidate::LorentzVector> smearedElectronMomentumVector;
  
  //-----For Ntuple
  vector<unsigned int> *_tauMotherId;
  vector<int> *_tauMatched;
  vector<int> *_tauPdgId;
  vector<int> *_tauMotherPdgId;
  vector<float> *_tauE;
  vector<float> *_tauEt; 
  vector<float> *_tauPt; 
  vector<float> *_tauCharge;
  vector<float> *_tauEta;
  vector<float> *_tauPhi;
  vector<unsigned int> *_tauNProngs;
  vector<float> *_tauIsoTrackPtSum;
  vector<float> *_tauIsoTrkPtSumDR1_0MinPt1_0;
  vector<float> *_tauIsoTrkPtSumDR1_0MinPt0_5;
  vector<float> *_tauIsoTrkPtSumDR0_75MinPt1_0;
  vector<float> *_tauIsoTrkPtSumDR0_75MinPt0_5;
  vector<float> *_tauIsoGammaEtSum;
  vector<float> *_tauIsoGammaEtSumDR1_0MinPt1_5;
  vector<float> *_tauIsoGammaEtSumDR1_0MinPt1_0;
  vector<float> *_tauIsoGammaEtSumDR0_75MinPt1_5;
  vector<float> *_tauIsoGammaEtSumDR0_75MinPt1_0;
  vector<float> *_tauLTPt;
  vector<float> *_tauLTChi2;
  vector<int> *_tauLTRecHitsSize;
  vector<float> *_tauEmFraction;
  vector<float> *_tauHcalTotOverPLead;
  vector<float> *_tauHcalMaxOverPLead;
  vector<float> *_tauHcal3x3OverPLead;
  vector<float> *_tauElectronPreId;
  vector<float> *_tauModifiedEOverP;
  vector<float> *_tauBremsRecoveryEOverPLead;
  vector<float> *_tauDiscAgainstElec;
  vector<int> *_tauLTCharge;
  vector<float> *_tauLTSignedIp;
  vector<int> *_tauIsInTheCraks;
     
  vector<unsigned int> *_eventIsZee;
  vector<unsigned int> *_eMotherId;
  vector<float> *_zeeMass;
  vector<float> *_zeePtAsymm;
  vector<int> *_eMatched;
  vector<int> *_ePdgId;
  vector<int> *_eMotherPdgId;
  vector<float> *_eE;
  vector<float> *_eEt;
  vector<float> *_ePt;
  vector<float> *_eCharge;
  vector<float> *_eEta;
  vector<float> *_ePhi;
  vector<float> *_eSigmaEtaEta;
  vector<float> *_eSigmaIEtaIEta;
  vector<float> *_eEOverP;
  vector<float> *_eHOverEm;
  vector<float> *_eDeltaPhiIn;
  vector<float> *_eDeltaEtaIn;
  
  vector<float> *_eEcalIsoPat;
  vector<float> *_eHcalIsoPat;
  vector<float> *_eTrkIsoPat;
  vector<float> *_eIsoPat;
  
  vector<float> *_eSCE1x5;
  vector<float> *_eSCE2x5;
  vector<float> *_eSCE5x5;
  vector<float> *_eIp;
  vector<int> *_eClass;
  
  vector<float> *_mEt;
  
  vector<float> *_eTauMass;
  vector<float> *_eTauCosDPhi;
  vector<float> *_eTauDelatR;
  vector<float> *_eTauMetMass;
  vector<float> *_eTauCollMass;
  vector<float> *_eTauPZeta;
  vector<float> *_eTauPZetaVis;
  
  vector<float> *_diTauMass;
  vector<float> *_diTauPt;
  vector<float> *_diTauEt;

};

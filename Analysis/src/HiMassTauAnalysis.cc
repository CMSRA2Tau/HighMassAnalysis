// Authors: Andres Florez, Alfredo Gurrola, Eduardo Luiggi, Chi Nhan Nguyen

#include "HighMassAnalysis/Analysis/interface/HiMassTauAnalysis.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include <TMath.h>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace edm;
using namespace reco;

// constructors and destructor
HiMassTauAnalysis::HiMassTauAnalysis(const ParameterSet& iConfig) {

  //-----Generator level Inputs 
  _GenParticleSource = iConfig.getUntrackedParameter<InputTag>("GenParticleSource");

  //-----Inputs to determine which channel to analyze
  _AnalyzeTauForLeg1 = iConfig.getParameter<bool>("AnalyzeTauForLeg1");
  _AnalyzeMuonForLeg1 = iConfig.getParameter<bool>("AnalyzeMuonForLeg1");
  _AnalyzeElectronForLeg1 = iConfig.getParameter<bool>("AnalyzeElectronForLeg1");
  _AnalyzeTauForLeg2 = iConfig.getParameter<bool>("AnalyzeTauForLeg2");
  _AnalyzeMuonForLeg2 = iConfig.getParameter<bool>("AnalyzeMuonForLeg2");
  _AnalyzeElectronForLeg2 = iConfig.getParameter<bool>("AnalyzeElectronForLeg2");

  //-----Reco Tau Inputs 
  _RecoTauSource = iConfig.getParameter<InputTag>("RecoTauSource");
  _RecoTauEtaCut = iConfig.getParameter<double>("RecoTauEtaCut");
  _RecoTauPtMinCut = iConfig.getParameter<double>("RecoTauPtMinCut");
  _RecoTauPtMaxCut = iConfig.getParameter<double>("RecoTauPtMaxCut");
  _DoRecoTauDiscrByLeadTrack = iConfig.getParameter<bool>("DoRecoTauDiscrByLeadTrack");
  _UseRecoTauDiscrByLeadTrackFlag = iConfig.getParameter<bool>("UseRecoTauDiscrByLeadTrackFlag");
  _RecoTauDiscrByLeadTrack = iConfig.getUntrackedParameter<string>("RecoTauDiscrByLeadTrack");
  _DoRecoTauDiscrByLeadTrackNhits = iConfig.getParameter<bool>("DoRecoTauDiscrByLeadTrackNhits");
  _RecoTauLeadTrackMinHits = iConfig.getParameter<int>("RecoTauLeadTrackMinHits");
  _DoRecoTauDiscrByIsolation = iConfig.getParameter<bool>("DoRecoTauDiscrByIsolation");
  _UseRecoTauDiscrByIsolationFlag = iConfig.getParameter<bool>("UseRecoTauDiscrByIsolationFlag");
  _UseRecoTauIsoSumPtInsteadOfNiso = iConfig.getParameter<bool>("UseRecoTauIsoSumPtInsteadOfNiso");
  _UseRecoTauEllipseForEcalIso = iConfig.getParameter<bool>("UseRecoTauEllipseForEcalIso");
  _RecoTauEcalIsoRphiForEllipse = iConfig.getParameter<double>("RecoTauEcalIsoRphiForEllipse");
  _RecoTauEcalIsoRetaForEllipse = iConfig.getParameter<double>("RecoTauEcalIsoRetaForEllipse");
  _RecoTauNisoMax = iConfig.getParameter<int>("RecoTauNisoMax");
  _RecoTauIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoTauIsoSumPtMinCutValue");
  _RecoTauIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoTauIsoSumPtMaxCutValue");
  _RecoTauDiscrByIsolation = iConfig.getUntrackedParameter<string>("RecoTauDiscrByIsolation");
  _RecoTauDiscrByProngType = iConfig.getParameter<string>("RecoTauDiscrByProngType");
  _RecoTauLeadTrackThreshold = iConfig.getParameter<double>("RecoTauLeadTrackThreshold");
  _RecoTauSigGamThreshold = iConfig.getParameter<double>("RecoTauSigGamThreshold");
  _RecoTauIsoDeltaRCone = iConfig.getParameter<double>("RecoTauIsoDeltaRCone");
  _RecoTauTrackIsoTrkThreshold = iConfig.getParameter<double>("RecoTauTrackIsoTrkThreshold");
  _RecoTauGammaIsoGamThreshold = iConfig.getParameter<double>("RecoTauGammaIsoGamThreshold");
  _DoRecoTauDiscrAgainstElectron = iConfig.getParameter<bool>("DoRecoTauDiscrAgainstElectron");
  _RecoTauDiscrAgainstElectron = iConfig.getUntrackedParameter<string>("RecoTauDiscrAgainstElectron");
  _DoRecoTauDiscrByCrackCut = iConfig.getParameter<bool>("DoRecoTauDiscrByCrackCut");
  _DoRecoTauDiscrAgainstMuon = iConfig.getParameter<bool>("DoRecoTauDiscrAgainstMuon");
  _RecoTauDiscrAgainstMuon = iConfig.getUntrackedParameter<string>("RecoTauDiscrAgainstMuon");

  //-----Reco Muon Inputs
  _RecoMuonSource = iConfig.getParameter<InputTag>("RecoMuonSource");
  _RecoMuonEtaCut = iConfig.getParameter<double>("RecoMuonEtaCut");
  _RecoMuonPtMinCut = iConfig.getParameter<double>("RecoMuonPtMinCut");
  _RecoMuonPtMaxCut = iConfig.getParameter<double>("RecoMuonPtMaxCut");
  _DoRecoMuonDiscrByGlobal = iConfig.getParameter<bool>("DoRecoMuonDiscrByGlobal");
  _DoRecoMuonDiscrByIsolation = iConfig.getParameter<bool>("DoRecoMuonDiscrByIsolation");
  _RecoMuonIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoMuonIsoSumPtMinCutValue");
  _RecoMuonIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoMuonIsoSumPtMaxCutValue");
  _RecoMuonIsoDeltaRCone = iConfig.getParameter<double>("RecoMuonIsoDeltaRCone");
  _RecoMuonTrackIsoTrkThreshold = iConfig.getParameter<double>("RecoMuonTrackIsoTrkThreshold");
  _RecoMuonEcalIsoRecHitThreshold = iConfig.getParameter<double>("RecoMuonEcalIsoRecHitThreshold");
  _DoRecoMuonDiscrByIp = iConfig.getParameter<bool>("DoRecoMuonDiscrByIp");
  _RecoMuonIpCut = iConfig.getParameter<double>("RecoMuonIpCut");

  //-----Reco Electron Inputs
  _RecoElectronSource = iConfig.getParameter<InputTag>("RecoElectronSource");
  _RecoElectronEtaCut = iConfig.getParameter<double>("RecoElectronEtaCut");
  _RecoElectronPtMinCut = iConfig.getParameter<double>("RecoElectronPtMinCut");
  _RecoElectronPtMaxCut = iConfig.getParameter<double>("RecoElectronPtMaxCut");
  _DoRecoElectronDiscrByTrackIsolation = iConfig.getParameter<bool>("DoRecoElectronDiscrByTrackIsolation");
  _RecoElectronTrackIsoSumPtCutValue = iConfig.getParameter<double>("RecoElectronTrackIsoSumPtCutValue");
  _RecoElectronTrackIsoDeltaRCone = iConfig.getParameter<double>("RecoElectronTrackIsoDeltaRCone");
  _RecoElectronTrackIsoTrkThreshold = iConfig.getParameter<double>("RecoElectronTrackIsoTrkThreshold");
  _DoRecoElectronDiscrByEcalIsolation = iConfig.getParameter<bool>("DoRecoElectronDiscrByEcalIsolation");
  _RecoElectronEcalIsoSumPtCutValue = iConfig.getParameter<double>("RecoElectronEcalIsoSumPtCutValue");
  _RecoElectronEcalIsoDeltaRCone = iConfig.getParameter<double>("RecoElectronEcalIsoDeltaRCone");
  _RecoElectronEcalIsoRecHitThreshold = iConfig.getParameter<double>("RecoElectronEcalIsoRecHitThreshold");
  _DoRecoElectronDiscrByIp = iConfig.getParameter<bool>("DoRecoElectronDiscrByIp");
  _RecoElectronIpCut = iConfig.getParameter<double>("RecoElectronIpCut");
  _DoRecoElectronDiscrByEoverP = iConfig.getParameter<bool>("DoRecoElectronDiscrByEoverP");
  _RecoElectronEoverPMax = iConfig.getParameter<double>("RecoElectronEoverPMax");
  _RecoElectronEoverPMin = iConfig.getParameter<double>("RecoElectronEoverPMin");
  _DoRecoElectronDiscrByHoverEm = iConfig.getParameter<bool>("DoRecoElectronDiscrByHoverEm");
  _RecoElectronHoverEmCut = iConfig.getParameter<double>("RecoElectronHoverEmCut");
  _DoRecoElectronDiscrBySigmaEtaEta = iConfig.getParameter<bool>("DoRecoElectronDiscrBySigmaEtaEta");
  _RecoElectronSigmaEtaEtaCut = iConfig.getParameter<double>("RecoElectronSigmaEtaEtaCut");
  _DoRecoElectronDiscrBySigmaIEtaIEta = iConfig.getParameter<bool>("DoRecoElectronDiscrBySigmaIEtaIEta");
  _RecoElectronSigmaIEtaIEtaCut = iConfig.getParameter<double>("RecoElectronSigmaIEtaIEtaCut");
  _DoRecoElectronDiscrBySCE5by5 = iConfig.getParameter<bool>("DoRecoElectronDiscrBySCE5by5");
  _RecoElectronSCE5by5Cut = iConfig.getParameter<double>("RecoElectronSCE5by5Cut");

  //-----Reco Jet Inputs
  _RecoJetSource = iConfig.getParameter<InputTag>("RecoJetSource");
  _RecoJetEtaMinCut = iConfig.getParameter<double>("RecoJetEtaMinCut");
  _RecoJetEtaMaxCut = iConfig.getParameter<double>("RecoJetEtaMaxCut");
  _RecoJetPtCut = iConfig.getParameter<double>("RecoJetPtCut");
  _UseCorrectedJet = iConfig.getParameter<bool>("UseCorrectedJet");
  _RemoveJetOverlapWithMuons = iConfig.getParameter<bool>("RemoveJetOverlapWithMuons");
  _JetMuonMatchingDeltaR = iConfig.getParameter<double>("JetMuonMatchingDeltaR");
  _RemoveJetOverlapWithElectrons = iConfig.getParameter<bool>("RemoveJetOverlapWithElectrons");
  _JetElectronMatchingDeltaR = iConfig.getParameter<double>("JetElectronMatchingDeltaR");
  _RemoveJetOverlapWithTaus = iConfig.getParameter<bool>("RemoveJetOverlapWithTaus");
  _JetTauMatchingDeltaR = iConfig.getParameter<double>("JetTauMatchingDeltaR");

  //-----Vertex Inputs
  _RecoVertexSource = iConfig.getParameter<InputTag>("RecoVertexSource");
  _RecoVertexMaxZposition = iConfig.getParameter<double>("RecoVertexMaxZposition");
  _RecoVertexMinTracks = iConfig.getParameter<int>("RecoVertexMinTracks");
  _RecoVertexTrackWeight = iConfig.getParameter<double>("RecoVertexTrackWeight");

  //-----Trigger Inputs
  _RecoTriggerSource = iConfig.getParameter<InputTag>("RecoTriggerSource");
  _TriggerRequirements = iConfig.getParameter<std::vector<std::string> >("TriggerRequirements");

  //-----Topology Inputs
  _RecoDiTauSource = iConfig.getParameter<InputTag>("RecoDiTauSource");
  _RecoMetSource = iConfig.getParameter<InputTag>("RecoMetSource");
  _DoDiscrByMet = iConfig.getParameter<bool>("DoDiscrByMet");
  _RecoMetCut = iConfig.getParameter<double>("RecoMetCut");
  _DoDiTauDiscrByDeltaR = iConfig.getParameter<bool>("DoDiTauDiscrByDeltaR");
  _DiTauDeltaRCut = iConfig.getParameter<double>("DiTauDeltaRCut");
  _UseTauSeedTrackForDiTauDiscrByOSLS = iConfig.getParameter<bool>("UseTauSeedTrackForDiTauDiscrByOSLS");
  _DiTauDiscrByOSLSType = iConfig.getParameter<string>("DiTauDiscrByOSLSType");
  _DoDiTauDiscrByCosDphi = iConfig.getParameter<bool>("DoDiTauDiscrByCosDphi");
  _DiTauCosDphiMinCut = iConfig.getParameter<double>("DiTauCosDphiMinCut");
  _DiTauCosDphiMaxCut = iConfig.getParameter<double>("DiTauCosDphiMaxCut");
  _DoDiscrByMassReco = iConfig.getParameter<bool>("DoDiscrByMassReco");
  _UseVectorSumOfVisProductsAndMetMassReco = iConfig.getParameter<bool>("UseVectorSumOfVisProductsAndMetMassReco");
  _UseCollinerApproxMassReco = iConfig.getParameter<bool>("UseCollinerApproxMassReco");
  _MassMinCut = iConfig.getParameter<double>("MassMinCut");
  _MassMaxCut = iConfig.getParameter<double>("MassMaxCut");
  _DoDiTauDiscrByCDFzeta2D = iConfig.getParameter<bool>("DoDiTauDiscrByCDFzeta2D");
  _PZetaCutCoefficient = iConfig.getParameter<double>("PZetaCutCoefficient");
  _PZetaVisCutCoefficient = iConfig.getParameter<double>("PZetaVisCutCoefficient");
  _CDFzeta2DCutValue = iConfig.getParameter<double>("CDFzeta2DCutValue");
  _DoDiTauDiscrByDeltaPtDivSumPt = iConfig.getParameter<bool>("DoDiTauDiscrByDeltaPtDivSumPt");
  _DiTauDeltaPtDivSumPtMinCutValue = iConfig.getParameter<double>("DiTauDeltaPtDivSumPtMinCutValue");
  _DiTauDeltaPtDivSumPtMaxCutValue = iConfig.getParameter<double>("DiTauDeltaPtDivSumPtMaxCutValue");
  _DoDiscrByLeg1MetDphi = iConfig.getParameter<bool>("DoDiscrByLeg1MetDphi");
  _Leg1MetDphiMinCut = iConfig.getParameter<double>("Leg1MetDphiMinCut");
  _Leg1MetDphiMaxCut = iConfig.getParameter<double>("Leg1MetDphiMaxCut");
  _DoDiscrByLeg2MetDphi = iConfig.getParameter<bool>("DoDiscrByLeg2MetDphi");
  _Leg2MetDphiMinCut = iConfig.getParameter<double>("Leg2MetDphiMinCut");
  _Leg2MetDphiMaxCut = iConfig.getParameter<double>("Leg2MetDphiMaxCut");

  //-----do matching to gen?
  _MatchLeptonToGen = iConfig.getParameter<bool>("MatchLeptonToGen");
  _UseLeptonMotherId = iConfig.getParameter<bool>("UseLeptonMotherId");
  _UseLeptonGrandMotherId = iConfig.getParameter<bool>("UseLeptonGrandMotherId");
  _LeptonMotherId = iConfig.getParameter<int>("LeptonMotherId");
  _LeptonGrandMotherId = iConfig.getParameter<int>("LeptonGrandMotherId");
  _MatchTauToGen = iConfig.getParameter<bool>("MatchTauToGen");
  _UseTauMotherId = iConfig.getParameter<bool>("UseTauMotherId");
  _UseTauGrandMotherId = iConfig.getParameter<bool>("UseTauGrandMotherId");
  _TauMotherId = iConfig.getParameter<int>("TauMotherId");
  _TauGrandMotherId = iConfig.getParameter<int>("TauGrandMotherId");
  _TauToGenMatchingDeltaR = iConfig.getParameter<double>("TauToGenMatchingDeltaR");

  //-----ntuple Inputs
  _DoProduceNtuple = iConfig.getParameter<bool>("DoProduceNtuple");
  _NtupleTreeName = (iConfig.getUntrackedParameter<std::string>("NtupleTreeName"));

  //-----Fill Histograms?
  _FillRecoVertexHists = iConfig.getParameter<bool>("FillRecoVertexHists");
  _FillGenTauHists = iConfig.getParameter<bool>("FillGenTauHists");
  _FillRecoTauHists = iConfig.getParameter<bool>("FillRecoTauHists");
  _FillRecoMuonHists = iConfig.getParameter<bool>("FillRecoMuonHists");
  _FillRecoElectronHists = iConfig.getParameter<bool>("FillRecoElectronHists");
  _FillRecoJetHists = iConfig.getParameter<bool>("FillRecoJetHists");
  _FillTopologyHists = iConfig.getParameter<bool>("FillTopologyHists");

  //-----Event Sequence inputs
  _RecoTriggersNmin = iConfig.getParameter<int>("RecoTriggersNmin");
  _RecoVertexNmin = iConfig.getParameter<int>("RecoVertexNmin");
  _RecoVertexNmax = iConfig.getParameter<int>("RecoVertexNmax");
  _RecoLeg1Nmin = iConfig.getParameter<int>("RecoLeg1Nmin");
  _RecoLeg1Nmax = iConfig.getParameter<int>("RecoLeg1Nmax");
  _RecoLeg2Nmin = iConfig.getParameter<int>("RecoLeg2Nmin");
  _RecoLeg2Nmax = iConfig.getParameter<int>("RecoLeg2Nmax");
  _RecoJetNmin = iConfig.getParameter<int>("RecoJetNmin");
  _RecoJetNmax = iConfig.getParameter<int>("RecoJetNmax");
  _CombinationsNmin = iConfig.getParameter<int>("CombinationsNmin");
  _CombinationsNmax = iConfig.getParameter<int>("CombinationsNmax");
  _EventSelectionSequence = iConfig.getParameter< vector<string> >("EventSelectionSequence");

  //-----Inputs for systematic uncertainties
  _CalculatePdfSystematicUncertanties = iConfig.getParameter<bool>("CalculatePdfSystematicUncertanties");
  pdfWeightTags_ = iConfig.getUntrackedParameter<vector<InputTag> >("PdfWeightTags");
  _SmearTheMuon = iConfig.getParameter<bool>("SmearTheMuon");
  _RelativeMuonPtOffset = iConfig.getParameter<double>("RelativeMuonPtOffset");
  _RelativeMuonPtSigma = iConfig.getParameter<double>("RelativeMuonPtSigma");
  _AbsoluteMuonEtaOffset = iConfig.getParameter<double>("AbsoluteMuonEtaOffset");
  _AbsoluteMuonEtaSigma = iConfig.getParameter<double>("AbsoluteMuonEtaSigma");
  _AbsoluteMuonPhiOffset = iConfig.getParameter<double>("AbsoluteMuonPhiOffset");
  _AbsoluteMuonPhiSigma = iConfig.getParameter<double>("AbsoluteMuonPhiSigma");
  _SmearTheElectron = iConfig.getParameter<bool>("SmearTheElectron");
  _RelativeElectronPtOffset = iConfig.getParameter<double>("RelativeElectronPtOffset");
  _RelativeElectronPtSigma = iConfig.getParameter<double>("RelativeElectronPtSigma");
  _AbsoluteElectronEtaOffset = iConfig.getParameter<double>("AbsoluteElectronEtaOffset");
  _AbsoluteElectronEtaSigma = iConfig.getParameter<double>("AbsoluteElectronEtaSigma");
  _AbsoluteElectronPhiOffset = iConfig.getParameter<double>("AbsoluteElectronPhiOffset");
  _AbsoluteElectronPhiSigma = iConfig.getParameter<double>("AbsoluteElectronPhiSigma");
}

// ------------ method called once each job just before starting event loop  ------------
void  HiMassTauAnalysis::beginJob(const EventSetup&) {
  _totalEvents = 0;  
  _totalEventsPassingCuts = 0;  
  if(_CalculatePdfSystematicUncertanties) {InitializeInfoForPDFSystematicUncertaintites();}
  setMapSelectionAlgoIDs();
  initMapSelectionCounters();
  bookHistograms();  
  if (_DoProduceNtuple) {
    initializeVectors();
    setupBranches();
  }
}

void  HiMassTauAnalysis::setupBranches() {

  _HMTTree = new TTree(_NtupleTreeName.c_str(), "HiMassDiTau Tree");
  //_HMTTree->Branch("Taus",&_totalEvents,"Taus/I");
  _HMTTree->Branch("tauMatched",&_tauMatched);
  _HMTTree->Branch("tauMotherId",&_tauMotherId);
  _HMTTree->Branch("tauPdgId",&_tauPdgId);
  _HMTTree->Branch("tauMotherPdgId",&_tauMotherPdgId);
  _HMTTree->Branch("tauE",&_tauE);
  _HMTTree->Branch("tauEt",&_tauEt);
  _HMTTree->Branch("tauPt",&_tauPt);
  _HMTTree->Branch("tauCharge",&_tauCharge);
  _HMTTree->Branch("tauEta",&_tauEta);
  _HMTTree->Branch("tauPhi",&_tauPhi);
  _HMTTree->Branch("tauNProngs",&_tauNProngs);
  _HMTTree->Branch("tauLTPt",&_tauLTPt);
  //_HMTTree->Branch("tauLTChi2",&_tauLTChi2);
  //_HMTTree->Branch("tauLTRecHitsSize",&_tauLTRecHitsSize);
  _HMTTree->Branch("tauIsoTrackSumPt",&_tauIsoTrackPtSum);
  _HMTTree->Branch("tauIsoTrackSumPtDR1_0MinPt1_0",&_tauIsoTrkPtSumDR1_0MinPt1_0);
  _HMTTree->Branch("tauIsoTrackSumPtDR1_0MinPt0_5",&_tauIsoTrkPtSumDR1_0MinPt0_5);
  _HMTTree->Branch("tauIsoTrackSumPtDR0_75MinPt1_0",&_tauIsoTrkPtSumDR0_75MinPt1_0);
  _HMTTree->Branch("tauIsoTrackSumPtDR0_75MinPt0_5",&_tauIsoTrkPtSumDR0_75MinPt0_5);
  
  _HMTTree->Branch("tauIsoGammaEtSum",&_tauIsoGammaEtSum);
  _HMTTree->Branch("tauIsoGammaEtSumDR1_0MinPt1_5",&_tauIsoGammaEtSumDR1_0MinPt1_5);
  _HMTTree->Branch("tauIsoGammaEtSumDR1_0MinPt1_0",&_tauIsoGammaEtSumDR1_0MinPt1_0);
  _HMTTree->Branch("tauIsoGammaEtSumDR0_75MinPt1_5",&_tauIsoGammaEtSumDR0_75MinPt1_5);
  _HMTTree->Branch("tauIsoGammaEtSumDR0_75MinPt1_0",&_tauIsoGammaEtSumDR0_75MinPt1_0);
  
  _HMTTree->Branch("tauEmFraction",&_tauEmFraction);
  _HMTTree->Branch("tauHcalTotOverPLead",&_tauHcalTotOverPLead);
  _HMTTree->Branch("tauHcalMaxOverPLead",&_tauHcalMaxOverPLead);
  _HMTTree->Branch("tauHcal3x3OverPLead",&_tauHcal3x3OverPLead);
  _HMTTree->Branch("tauElectronPreId",&_tauElectronPreId);
  _HMTTree->Branch("tauModifiedEOverP",&_tauModifiedEOverP);
  _HMTTree->Branch("tauBremsRecoveryEOverPLead",&_tauBremsRecoveryEOverPLead);
  _HMTTree->Branch("tauDiscAgainstElec",&_tauDiscAgainstElec);
  _HMTTree->Branch("tauLTCharge",&_tauLTCharge);
  _HMTTree->Branch("tauLTSignedIp",&_tauLTSignedIp);
  _HMTTree->Branch("tauIsInTheCraks",&_tauIsInTheCraks);
 
  _HMTTree->Branch("eventIsZee",&_eventIsZee);
  _HMTTree->Branch("zeeMass",&_zeeMass);
  _HMTTree->Branch("zeePtAsymm",&_zeePtAsymm);
  _HMTTree->Branch("eMatched",&_eMatched);
  _HMTTree->Branch("eMotherId",&_eMotherId);
  _HMTTree->Branch("ePdgId",&_ePdgId);
  _HMTTree->Branch("eMotherPdgId",&_eMotherPdgId);
  _HMTTree->Branch("eE",&_eE);
  _HMTTree->Branch("eEt",&_eEt);
  _HMTTree->Branch("ePt",&_ePt);
  _HMTTree->Branch("eCharge",&_eCharge);
  _HMTTree->Branch("eEta",&_eEta);
  _HMTTree->Branch("ePhi",&_ePhi);
  _HMTTree->Branch("eSigmaEtaEta",&_eSigmaEtaEta);
  _HMTTree->Branch("eSigmaIEtaIEta",&_eSigmaIEtaIEta);
  _HMTTree->Branch("eEOverP",&_eEOverP);
  _HMTTree->Branch("eHOverEm",&_eHOverEm);
  _HMTTree->Branch("eDeltaPhiIn", &_eDeltaPhiIn);
  _HMTTree->Branch("eDeltaEtaIn", &_eDeltaEtaIn);
 
  _HMTTree->Branch("eEcalIsoPat", &_eEcalIsoPat);
  _HMTTree->Branch("eHcalIsoPat", &_eHcalIsoPat);
  _HMTTree->Branch("eTrkIsoPat", &_eTrkIsoPat);
  _HMTTree->Branch("eIsoPat", &_eIsoPat);

  _HMTTree->Branch("eSCE1x5", &_eSCE1x5);
  _HMTTree->Branch("eSCE2x5", &_eSCE2x5);
  _HMTTree->Branch("eSCE5x5", &_eSCE5x5);
  _HMTTree->Branch("eIp", &_eIp);
  _HMTTree->Branch("eClass", &_eClass);
 
  _HMTTree->Branch("mEt", &_mEt);
    
  _HMTTree->Branch("eTauMass",&_eTauMass);
  _HMTTree->Branch("eTauCosDPhi",&_eTauCosDPhi);
  _HMTTree->Branch("eTauDeltaR",&_eTauDelatR);
  _HMTTree->Branch("eTauMetMass",&_eTauMetMass);
  _HMTTree->Branch("eTauCollMass",&_eTauCollMass);
 
  _HMTTree->Branch("eTauPZeta", &_eTauPZeta);
  _HMTTree->Branch("eTauPZetaVis", &_eTauPZetaVis);

  _HMTTree->Branch("diTauMass",&_diTauMass);
  _HMTTree->Branch("diTauPt",&_diTauPt);
  _HMTTree->Branch("diTauEt",&_diTauEt);  

}

void HiMassTauAnalysis::initializeVectors(){

  _tauMatched = 0;
  _tauMotherId = 0;
  _tauPdgId = 0;
  _tauMotherPdgId = 0;
  _tauE = 0;
  _tauEt = 0; 
  _tauPt = 0; 
  _tauCharge = 0;
  _tauEta = 0;
  _tauPhi = 0;
  _tauNProngs = 0;
  _tauLTPt = 0;
  //_tauLTChi2 = 0;
  //_tauLTRecHitsSize = 0;
  _tauIsoTrackPtSum = 0;
  _tauIsoTrkPtSumDR1_0MinPt1_0 = 0;
  _tauIsoTrkPtSumDR1_0MinPt0_5 = 0;
  _tauIsoTrkPtSumDR0_75MinPt1_0 = 0;
  _tauIsoTrkPtSumDR0_75MinPt0_5 = 0;
  _tauIsoGammaEtSum = 0;
  _tauIsoGammaEtSumDR1_0MinPt1_5 = 0;
  _tauIsoGammaEtSumDR1_0MinPt1_0 = 0;
  _tauIsoGammaEtSumDR0_75MinPt1_5 = 0;
  _tauIsoGammaEtSumDR0_75MinPt1_0 = 0;
  _tauEmFraction = 0;
  _tauHcalTotOverPLead = 0;
  _tauHcalMaxOverPLead = 0;
  _tauHcal3x3OverPLead = 0;
  _tauElectronPreId = 0;
  _tauModifiedEOverP = 0;
  _tauBremsRecoveryEOverPLead = 0;
  _tauDiscAgainstElec = 0;
  _tauLTCharge = 0;
  _tauLTSignedIp = 0;
  _tauIsInTheCraks = 0;
  
  _eventIsZee = 0;
  _zeeMass = 0;
  _zeePtAsymm = 0;
  _eMatched = 0;
  _eMotherId = 0;
  _ePdgId = 0;
  _eMotherPdgId = 0;
  _eE = 0;
  _eEt = 0;
  _ePt = 0;
  _eCharge = 0;
  _eEta = 0;
  _ePhi = 0;
  _eSigmaEtaEta = 0;
  _eSigmaIEtaIEta = 0;
  _eEOverP = 0;
  _eHOverEm = 0;
  _eDeltaPhiIn = 0;
  _eDeltaEtaIn = 0;
  
  _eEcalIsoPat = 0;
  _eHcalIsoPat = 0;
  _eTrkIsoPat = 0;
  _eIsoPat = 0;
  
  _eSCE1x5 = 0;
  _eSCE2x5 = 0;
  _eSCE5x5 = 0;
  _eIp = 0;
  _eClass = 0;
  
  _mEt =0;
  
  _eTauMass = 0;
  _eTauCosDPhi = 0;
  _eTauDelatR = 0;
  _eTauMetMass = 0;
  _eTauCollMass = 0;
  _eTauPZeta = 0;
  _eTauPZetaVis = 0;
    
  _diTauMass = 0;
  _diTauPt = 0;
  _diTauEt = 0;

}

void HiMassTauAnalysis::clearVectors(){

  _tauMotherId->clear();
  _tauMatched->clear();
  _tauPdgId->clear();
  _tauMotherPdgId->clear();
  _tauE->clear();
  _tauEt->clear(); 
  _tauPt->clear(); 
  _tauCharge->clear();
  _tauEta->clear();
  _tauPhi->clear();
  _tauNProngs->clear();
  _tauIsoTrackPtSum->clear();
  _tauIsoTrkPtSumDR1_0MinPt1_0->clear();
  _tauIsoTrkPtSumDR1_0MinPt0_5->clear();
  _tauIsoTrkPtSumDR0_75MinPt1_0->clear();
  _tauIsoTrkPtSumDR0_75MinPt0_5->clear();
  
  _tauIsoGammaEtSum->clear();
  _tauIsoGammaEtSumDR1_0MinPt1_5->clear();
  _tauIsoGammaEtSumDR1_0MinPt1_0->clear();
  _tauIsoGammaEtSumDR0_75MinPt1_5->clear();
  _tauIsoGammaEtSumDR0_75MinPt1_0->clear();
  _tauLTPt->clear();
  //_tauLTChi2->clear();
  //_tauLTRecHitsSize->clear();
  _tauEmFraction->clear();
  _tauHcalTotOverPLead->clear();
  _tauHcalMaxOverPLead->clear();
  _tauHcal3x3OverPLead->clear();
  _tauElectronPreId->clear();
  _tauModifiedEOverP->clear();
  _tauBremsRecoveryEOverPLead->clear();
  _tauDiscAgainstElec->clear();
  _tauLTCharge->clear();
  _tauLTSignedIp->clear();
  _tauIsInTheCraks->clear();
  
  _eventIsZee->clear();
  _zeeMass->clear();
  _zeePtAsymm->clear();
  _eMotherId->clear();
  _eMatched->clear();
  _ePdgId->clear();
  _eMotherPdgId->clear();
  _eE->clear();
  _eEt->clear();
  _ePt->clear();
  _eCharge->clear();
  _eEta->clear();
  _ePhi->clear();
  _eSigmaEtaEta->clear();
  _eSigmaIEtaIEta->clear();
  _eEOverP->clear();
  _eHOverEm->clear();
  _eDeltaPhiIn->clear();
  _eDeltaEtaIn->clear();
  
  _eEcalIsoPat->clear();
  _eHcalIsoPat->clear();
  _eTrkIsoPat->clear();
  _eIsoPat->clear();
 
  _eSCE1x5->clear();
  _eSCE2x5->clear();
  _eSCE5x5->clear();
  _eIp->clear();
  _eClass->clear();
  
  _mEt->clear();
    
  _eTauMass->clear();
  _eTauCosDPhi->clear();
  _eTauDelatR->clear();
  _eTauMetMass->clear();
  _eTauCollMass->clear();
  _eTauPZeta->clear();
  _eTauPZetaVis->clear();
 
  _diTauMass->clear();
  _diTauPt->clear();
  _diTauEt->clear();

}

// ------------ method called to for each event  ------------
void HiMassTauAnalysis::analyze(const Event& iEvent, const EventSetup& iSetup) {

  //------Number of events analyzed (denominator)
  _totalEvents++;
  _hEvents->Fill(0);

  //-----Get weights for the calculation of pdf systematic uncertainties for the denominator
  if(_CalculatePdfSystematicUncertanties) {
    for (unsigned int i=0; i<pdfWeightTags_.size(); ++i) {
//      std::cout << "pdf tag = " << pdfWeightTags_[i] << std::endl;
      edm::Handle<std::vector<double> > weightHandle;
      iEvent.getByLabel(pdfWeightTags_[i], weightHandle);
      std::vector<double> weights = (*weightHandle);
      unsigned int nmembers = weights.size();
      if (pdfStart_Denominator_[i]<0) {
        pdfStart_Denominator_[i] = weightedEvents_Denominator_.size();
        for (unsigned int j=0; j<nmembers; ++j) {weightedEvents_Denominator_.push_back(0.);}
      }
      for (unsigned int j=0; j<nmembers; ++j) {
        weightedEvents_Denominator_[pdfStart_Denominator_[i]+j] += weights[j];
//        std::cout << "weight member " << j << " = " << weights[j] << std::endl;
      }
    }
  }

  //------Grab the handle to the relevant collections
//  std::cout << "getting collections ..." << std::endl;
  getCollections(iEvent,iSetup);

  //-----Smearing light lepton (e/mu) momentum and position for systematic uncertanties
//  std::cout << "smearing muon ..." << std::endl;
  smearedMuonMomentumVector.clear();
/*
  if(((_AnalyzeMuonForLeg1) || (_AnalyzeMuonForLeg2)) && (_SmearTheMuon)) {
    for(pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end();++patMuon) {
      smearedMuonMomentumVector.push_back(SmearLightLepton(*patMuon));
    }
  }
*/
  if(((_AnalyzeMuonForLeg1) || (_AnalyzeMuonForLeg2))) {
    if(_SmearTheMuon) {
      for(pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end();++patMuon) {
        smearedMuonMomentumVector.push_back(SmearLightLepton(*patMuon));
      }
    } else {
      for(pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end();++patMuon) {
        smearedMuonMomentumVector.push_back(patMuon->p4());
      }
    }
  }
//  std::cout << "smearing electron ..." << std::endl;
  smearedElectronMomentumVector.clear();
/*
  if(((_AnalyzeElectronForLeg1) || (_AnalyzeElectronForLeg2)) && (_SmearTheElectron)) {
    for(pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end();++patElectron) {
      smearedElectronMomentumVector.push_back(SmearLightLepton(*patElectron));
    }
  }
*/
  if(((_AnalyzeElectronForLeg1) || (_AnalyzeElectronForLeg2))) {
    if(_SmearTheElectron) {
      for(pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end();++patElectron) {
        smearedElectronMomentumVector.push_back(SmearLightLepton(*patElectron));
      }
    } else {
      for(pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end();++patElectron) {
        smearedElectronMomentumVector.push_back(patElectron->p4());
      }
    }
  }

  //------Get the event flags (did the event pass the cuts?)
  getEventFlags();
  if (_DoProduceNtuple){
    fillNtuple();
    clearVectors();
  }
  if (!passEventSelectionSequence()) return;  

  //------Number of events passing cuts (numerator)
  _totalEventsPassingCuts++;  
  _hEvents->Fill(1);

  //-----Get weights for the calculation of pdf systematic uncertainties for the numerator
  if(_CalculatePdfSystematicUncertanties) {
    for (unsigned int i=0; i<pdfWeightTags_.size(); ++i) {
//      std::cout << "pdf tag = " << pdfWeightTags_[i] << std::endl;
      edm::Handle<std::vector<double> > weightHandle;
      iEvent.getByLabel(pdfWeightTags_[i], weightHandle);
      std::vector<double> weights = (*weightHandle);
      unsigned int nmembers = weights.size();
      if (pdfStart_Numerator_[i]<0) {
        pdfStart_Numerator_[i] = weightedEvents_Numerator_.size();
        for (unsigned int j=0; j<nmembers; ++j) {weightedEvents_Numerator_.push_back(0.);}
      }
      for (unsigned int j=0; j<nmembers; ++j) {
        weightedEvents_Numerator_[pdfStart_Numerator_[i]+j] += weights[j];
//        std::cout << "weight member " << j << " = " << weights[j] << std::endl;
      }
    }
  }

  //------If the event passed the cut criteria, then fill histograms and ntuple
  fillHistograms();
}

void  HiMassTauAnalysis::initMapSelectionCounters() {
  for (unsigned int i=0;i<_EventSelectionSequence.size();i++) {
    _mapSelectionCounter[_EventSelectionSequence[i]] = 0;
    _mapSelectionCounterCumul[_EventSelectionSequence[i]] = 0;
  }
}

void HiMassTauAnalysis::printEfficiency() {
  cout.setf(ios::floatfield,ios::fixed);
  cout<<setprecision(3);
  cout << "\n";
  cout << "Selection Efficiency " << "\n";
  cout << "Total events: " << _totalEvents << "\n";
  cout << "         Name                     Indiv.      Cumulative\n";
  cout << "---------------------------------------------------------------------------\n";
  for (unsigned int i=0;i<_EventSelectionSequence.size();i++) {
    cout<<setw(24)<<_EventSelectionSequence[i]<<" "
	<<setw(6)<<_mapSelectionCounter[_EventSelectionSequence[i]]<<" ("
	<<setw(8)<<(float)_mapSelectionCounter[_EventSelectionSequence[i]]/(float)_totalEvents<<") "
	<<setw(6)<<_mapSelectionCounterCumul[_EventSelectionSequence[i]]<<"( "
	<<setw(8)<<(float)_mapSelectionCounterCumul[_EventSelectionSequence[i]]/(float)_totalEvents<<") "
	<<endl;
  }
  cout << "---------------------------------------------------------------------------\n";  
}

void HiMassTauAnalysis::setMapSelectionAlgoIDs() {
  for (unsigned int i=0;i<_EventSelectionSequence.size();i++) {
    _mapSelectionAlgoID[_EventSelectionSequence[i]] = i;
  }
}

void HiMassTauAnalysis::getEventFlags() {

  //-----init event flags
  _EventFlag.clear();
  for (unsigned int i=0;i<_EventSelectionSequence.size();i++) { _EventFlag.push_back(false); }

  int nLeg1 = 0;
  if(_AnalyzeTauForLeg1) {nLeg1++;}
  if(_AnalyzeMuonForLeg1) {nLeg1++;}
  if(_AnalyzeElectronForLeg1) {nLeg1++;}
  int nLeg2 = 0;
  if(_AnalyzeTauForLeg2) {nLeg2++;}
  if(_AnalyzeMuonForLeg2) {nLeg2++;}
  if(_AnalyzeElectronForLeg2) {nLeg2++;}
  if(nLeg1 > 1) {std::cerr << "### HiMassTauAnalysis - CONFIGURATION ERROR:  Cannot analyze more than 1 object for leg1!!!! " << std::endl;exit(1);}
  if(nLeg2 > 1) {std::cerr << "### HiMassTauAnalysis - CONFIGURATION ERROR:  Cannot analyze more than 1 object for leg2!!!! " << std::endl;exit(1);}
  if(nLeg1 == 0) {std::cerr << "### HiMassTauAnalysis - CONFIGURATION ERROR:  ZERO objects for leg1!!!! " << std::endl;exit(1);}
  if(nLeg2 == 0) {std::cerr << "### HiMassTauAnalysis - CONFIGURATION ERROR:  ZERO objects for leg2!!!! " << std::endl;exit(1);}

  // ------Does the event pass trigger requirements?
//  std::cout << "trigger selections ..." << std::endl;
  int nTriggersSatisfied = 0;
  if(passRecoTriggerCuts()) {nTriggersSatisfied++;}
  if (nTriggersSatisfied>=_RecoTriggersNmin) _EventFlag[_mapSelectionAlgoID["RecoTriggersNmin"]] = true;

  // ------Number of Good Vertices
//  std::cout << "vertex selections ..." << std::endl;
  int nGoodVertices = 0;
  for ( reco::VertexCollection::const_iterator primaryVertex = _primaryEventVertexCollection->begin();
        primaryVertex != _primaryEventVertexCollection->end(); ++primaryVertex ) {
    if (!passRecoVertexCuts(*primaryVertex)) continue;
    nGoodVertices++;
  }
  if (nGoodVertices>=_RecoVertexNmin) _EventFlag[_mapSelectionAlgoID["RecoVertexNmin"]] = true;
  if (nGoodVertices<=_RecoVertexNmax) _EventFlag[_mapSelectionAlgoID["RecoVertexNmax"]] = true;

  //------Number of Good Candidates for leg1
//  std::cout << "leg1 selections ..." << std::endl;
  int nGoodCandidatesLeg1 = 0;
  if(_AnalyzeMuonForLeg1) {
    int theNumberOfMuons = 0;
    for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); 
  	  patMuon != _patMuons->end(); ++patMuon ) {
      theNumberOfMuons++;
      if (!passRecoMuonCuts((*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1))) continue;
      nGoodCandidatesLeg1++;
    }
  }
  if(_AnalyzeElectronForLeg1) {
    int theNumberOfElectrons = 0;
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();
          patElectron != _patElectrons->end(); ++patElectron ) {
      theNumberOfElectrons++;
      if (!passRecoElectronCuts((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1))) continue;
      nGoodCandidatesLeg1++;
    }
  }
  if(_AnalyzeTauForLeg1) {
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); 
	  patTau != _patTaus->end(); ++patTau ) {
      if (!passRecoTauCuts(*patTau)) continue;
      nGoodCandidatesLeg1++;
    }
  }
  if (nGoodCandidatesLeg1>=_RecoLeg1Nmin) _EventFlag[_mapSelectionAlgoID["RecoLeg1Nmin"]] = true;
  if (nGoodCandidatesLeg1<=_RecoLeg1Nmax) _EventFlag[_mapSelectionAlgoID["RecoLeg1Nmax"]] = true;

  //------Number of Good Candidates for leg2
//  std::cout << "leg2 selections ..." << std::endl;
  int nGoodCandidatesLeg2 = 0;
  if(_AnalyzeMuonForLeg2) {
    int theNumberOfMuons = 0;
    for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); 
  	  patMuon != _patMuons->end(); ++patMuon ) {
      theNumberOfMuons++;
      if (!passRecoMuonCuts((*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1))) continue;
      nGoodCandidatesLeg2++;
    }
  }
  if(_AnalyzeElectronForLeg2) {
    int theNumberOfElectrons = 0;
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();
          patElectron != _patElectrons->end(); ++patElectron ) {
      theNumberOfElectrons++;
      if (!passRecoElectronCuts((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1))) continue;
      nGoodCandidatesLeg2++;
    }
  }
  if(_AnalyzeTauForLeg2) {
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); 
	  patTau != _patTaus->end(); ++patTau ) {
      if (!passRecoTauCuts(*patTau)) continue;
      nGoodCandidatesLeg2++;
    }
  }
  if (nGoodCandidatesLeg2>=_RecoLeg2Nmin) _EventFlag[_mapSelectionAlgoID["RecoLeg2Nmin"]] = true;
  if (nGoodCandidatesLeg2<=_RecoLeg2Nmax) _EventFlag[_mapSelectionAlgoID["RecoLeg2Nmax"]] = true;

  // ------Number of Good Jets   
//  std::cout << "jet selections ..." << std::endl;
  int nGoodJets = 0;
  for ( pat::JetCollection::const_iterator patJet = _patJets->begin(); 
	patJet != _patJets->end(); ++patJet ) {
    if (!passRecoJetCuts(*patJet)) continue;
    nGoodJets++;
  }
  if (nGoodJets>=_RecoJetNmin) _EventFlag[_mapSelectionAlgoID["RecoJetNmin"]] = true;
  if (nGoodJets<=_RecoJetNmax) _EventFlag[_mapSelectionAlgoID["RecoJetNmax"]] = true;  

  // ------Number of Good Combinations (leg1+leg2+met combinations)
//  std::cout << "ditau selections ..." << std::endl;
  int nGoodCombinations = 0;
  if( ((_AnalyzeMuonForLeg1) && (_AnalyzeTauForLeg2)) || ((_AnalyzeMuonForLeg2) && (_AnalyzeTauForLeg1)) ) {
    int theNumberOfMuons = 0;
    for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
      theNumberOfMuons++;
      for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
        if ((passRecoTauCuts(*patTau)) && 
            (passRecoMuonCuts((*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1))) && 
            (passTopologyCuts((*patTau),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))))) {
          nGoodCombinations++;
        }
      }
    }
  }
  if( ((_AnalyzeElectronForLeg1) && (_AnalyzeTauForLeg2)) || ((_AnalyzeElectronForLeg2) && (_AnalyzeTauForLeg1)) ) {
    int theNumberOfElectrons = 0;
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
      theNumberOfElectrons++;
      for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
        if ((passRecoTauCuts(*patTau)) && 
            (passRecoElectronCuts((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1))) &&
            (passTopologyCuts((*patTau),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))))) {
          nGoodCombinations++;
        }
      }
    }
  }
  if( ((_AnalyzeMuonForLeg1) && (_AnalyzeElectronForLeg2)) || ((_AnalyzeMuonForLeg2) && (_AnalyzeElectronForLeg1)) ) {
    int theNumberOfMuons = 0;
    for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
      theNumberOfMuons++;
      int theNumberOfElectrons = 0;
      for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
        theNumberOfElectrons++;
        if ((passRecoElectronCuts((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1))) && 
            (passRecoMuonCuts((*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1))) && 
            (passTopologyCuts((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))))) {
          nGoodCombinations++;
        }
      }
    }
  }
  if( ((_AnalyzeTauForLeg1) && (_AnalyzeTauForLeg2)) ) {
    for ( pat::TauCollection::const_iterator patTau1 = _patTaus->begin();patTau1 != _patTaus->end(); ++patTau1 ) {
      for ( pat::TauCollection::const_iterator patTau2 = _patTaus->begin();patTau2 != _patTaus->end(); ++patTau2 ) {
        if ((passRecoTauCuts(*patTau1)) && (passRecoTauCuts(*patTau2)) && (passTopologyCuts((*patTau1),(*patTau2),(*(_patMETs->begin()))))) {
          nGoodCombinations++;
        }
      }
    }
  }
  if( ((_AnalyzeMuonForLeg1) && (_AnalyzeMuonForLeg2)) ) {
    int theNumberOfMuons1 = 0;
    for ( pat::MuonCollection::const_iterator patMuon1 = _patMuons->begin();patMuon1 != _patMuons->end(); ++patMuon1 ) {
      theNumberOfMuons1++;
      int theNumberOfMuons2 = 0;
      for ( pat::MuonCollection::const_iterator patMuon2 = _patMuons->begin();patMuon2 != _patMuons->end(); ++patMuon2 ) {
        theNumberOfMuons2++;
        if ((passRecoMuonCuts((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1))) && 
            (passRecoMuonCuts((*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1))) && 
            (passTopologyCuts((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))))) {
          nGoodCombinations++;
        }
      }
    }
  }
  if( ((_AnalyzeElectronForLeg1) && (_AnalyzeElectronForLeg2)) ) {
    int theNumberOfElectrons1 = 0;
    for ( pat::ElectronCollection::const_iterator patElectron1 = _patElectrons->begin();patElectron1 != _patElectrons->end(); ++patElectron1 ) {
      theNumberOfElectrons1++;
      int theNumberOfElectrons2 = 0;
      for ( pat::ElectronCollection::const_iterator patElectron2 = _patElectrons->begin();patElectron2 != _patElectrons->end(); ++patElectron2 ) {
        theNumberOfElectrons2++;
        if ((passRecoElectronCuts((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1))) && 
            (passRecoElectronCuts((*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1))) && 
            (passTopologyCuts((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))))) {
          nGoodCombinations++;
        }
      }
    }
  }
  if (nGoodCombinations>=_CombinationsNmin) _EventFlag[_mapSelectionAlgoID["CombinationsNmin"]] = true;
  if (nGoodCombinations<=_CombinationsNmax) _EventFlag[_mapSelectionAlgoID["CombinationsNmax"]] = true;

}

// --------Count number of events passing the selection criteria
bool HiMassTauAnalysis::passEventSelectionSequence() {
  bool cumulDecision = true;
  for (unsigned int i=0;i<_EventSelectionSequence.size();i++) {
    cumulDecision = cumulDecision && _EventFlag[i]; 
    if (_EventFlag[i]) {
      (_mapSelectionCounter[_EventSelectionSequence[i]])++;
    }
    if (cumulDecision) {
      (_mapSelectionCounterCumul[_EventSelectionSequence[i]])++;
    }
  }  
  return cumulDecision;
}

// -------------Apply Trigger Requirements
bool HiMassTauAnalysis::passRecoTriggerCuts() {
  edm::TriggerNames TheTriggerNames;
  TheTriggerNames.init(*(_triggerResults));
  for(std::vector<std::string>::const_iterator TheTriggerPath = _TriggerRequirements.begin();
      TheTriggerPath != _TriggerRequirements.end(); ++TheTriggerPath ) {
    unsigned int index = TheTriggerNames.triggerIndex(*TheTriggerPath);
    if(index < TheTriggerNames.size()) {if(_triggerResults->accept(index)) {return true;}}
    else {
      std::cout << "### HiMassTauAnalysis - CONFIGURATION ERROR:  Specified trigger " << (*TheTriggerPath) << " is not found/defined!!!!" << std::endl;
      std::cout << "Please use one of the following triggers :" << std::endl;
      for(edm::TriggerNames::Strings::const_iterator triggerName = TheTriggerNames.triggerNames().begin();
          triggerName != TheTriggerNames.triggerNames().end(); ++triggerName ) {
        unsigned int index = TheTriggerNames.triggerIndex(*triggerName);
        if(index < TheTriggerNames.size()) {
          std::string triggerDecision = (_triggerResults->accept(index)) ? "passed" : "failed";
          std::cout << " Trigger Name = " << (*triggerName) << " " << triggerDecision << std::endl;
        }
      }
//      exit(1);
    }
  }
  return false;
}

// -------------Apply Vertex Cuts
bool HiMassTauAnalysis::passRecoVertexCuts(const reco::Vertex& theVertex) {
  // ----remove fakes
  if(theVertex.isFake()) {return false;}
  // ----require vertex to be close to pp interaction point
  if(fabs(theVertex.z()) >= _RecoVertexMaxZposition) {return false;}
  // ----require vertex to have a minimum number of "good" tracks used for the fit
  if(_RecoVertexMinTracks >= 0) {
    int pvntrk = 0;
    reco::Vertex::trackRef_iterator pvtrk;
    for(pvtrk=theVertex.tracks_begin();pvtrk!=theVertex.tracks_end();++pvtrk) {
     if(theVertex.trackWeight(*pvtrk) > _RecoVertexTrackWeight) {pvntrk++;}
    }
    if(pvntrk < _RecoVertexMinTracks) {return false;}
  }
  return true;
}

// -------------Apply Tau Cuts
bool HiMassTauAnalysis::passRecoTauCuts(const pat::Tau& patTau) {
  // ----Matching to gen
  if(_MatchTauToGen) {if(!(matchToGen(patTau).first)) {return false;}}
  // ----Acceptance cuts
  if (abs(patTau.eta())>_RecoTauEtaCut) {return false;}
  if (patTau.pt()<_RecoTauPtMinCut) {return false;}
  if (patTau.pt()>_RecoTauPtMaxCut) {return false;}
  // ----Lead track requirement
  if (_DoRecoTauDiscrByLeadTrack) {
    if (_UseRecoTauDiscrByLeadTrackFlag) { if ( (patTau.tauID(_RecoTauDiscrByLeadTrack.data()) < 0.5) ) {return false;} }
    else {
      if (patTau.isCaloTau()) { if( (!(patTau.leadTrack().isNonnull())) || (patTau.leadTrack()->pt() < _RecoTauLeadTrackThreshold) ) {return false;} }
      else { if( (!(patTau.leadPFChargedHadrCand().isNonnull())) || (patTau.leadPFChargedHadrCand()->pt() < _RecoTauLeadTrackThreshold) ) {return false;} }
    }
  }
  // ----Lead track minimimum hits requirement
  if (_DoRecoTauDiscrByLeadTrackNhits) {
    if (patTau.isCaloTau()) {
      if( (!(patTau.leadTrack().isNonnull())) || ((int)(patTau.leadTrack()->recHitsSize()) < _RecoTauLeadTrackMinHits) ) {return false;}
    } else {
      if( (!(patTau.leadPFChargedHadrCand().isNonnull())) || (!(patTau.leadPFChargedHadrCand()->trackRef().isNonnull())) ) {return false;}
      if( (int)(patTau.leadPFChargedHadrCand()->trackRef()->recHitsSize()) < _RecoTauLeadTrackMinHits ) {return false;}
    }
  }
  // ----Isolation requirement
  if (_DoRecoTauDiscrByIsolation) {
    if (_UseRecoTauDiscrByIsolationFlag) {if ( (patTau.tauID(_RecoTauDiscrByIsolation.data()) < 0.5) ) {return false;}}
    else {
      if(_UseRecoTauIsoSumPtInsteadOfNiso) {
        if( (CalculateTauTrackIsolation(patTau).second + CalculateTauEcalIsolation(patTau).second) >= _RecoTauIsoSumPtMaxCutValue ) {return false;}
        if( (CalculateTauTrackIsolation(patTau).second + CalculateTauEcalIsolation(patTau).second) < _RecoTauIsoSumPtMinCutValue ) {return false;}
      } else { if( (CalculateTauTrackIsolation(patTau).first + CalculateTauEcalIsolation(patTau).first) > _RecoTauNisoMax ) {return false;} }
    }
  }
  // ----Require 1 or 3 prongs
  if (_RecoTauDiscrByProngType == "1or3") {
    if (patTau.isCaloTau()) {
      if((patTau.signalTracks().size() == 1) ||(patTau.signalTracks().size() == 3)) {}
      else {return false;}
    } else {
      if((patTau.signalPFChargedHadrCands().size() == 1) ||(patTau.signalPFChargedHadrCands().size() == 3)) {}
      else {return false;}
    }
  } else if (_RecoTauDiscrByProngType == "1") {
    if (patTau.isCaloTau()) {
      if(patTau.signalTracks().size() == 1) {}
      else {return false;}
    } else {
      if(patTau.signalPFChargedHadrCands().size() == 1) {}
      else {return false;}
    }
  } else if (_RecoTauDiscrByProngType == "3") {
    if (patTau.isCaloTau()) {
      if(patTau.signalTracks().size() == 3) {}
      else {return false;}
    } else {
      if(patTau.signalPFChargedHadrCands().size() == 3) {}
      else {return false;}
    }
  } else {}
  // ----Electron and Muon vetos
  if (_DoRecoTauDiscrAgainstElectron) { if ( (patTau.tauID(_RecoTauDiscrAgainstElectron.data()) < 0.5) ) {return false;} }
  if (_DoRecoTauDiscrByCrackCut) {
    if(isInTheCracks(patTau.eta())) {return false;}
  }
  if (_DoRecoTauDiscrAgainstMuon) { if ( (patTau.tauID(_RecoTauDiscrAgainstMuon.data()) < 0.5) ) {return false;} }
  return true;
}

// ---------------Apply Muon Cuts
bool HiMassTauAnalysis::passRecoMuonCuts(const pat::Muon& patMuon,bool smear,reco::Candidate::LorentzVector smearedLV) {
  // ----Matching to gen
  if(_MatchLeptonToGen) {if(!(matchToGen(patMuon).first)) {return false;}}
  // ----Require maching tracks in silicon tracker and muon chamber
  if (_DoRecoMuonDiscrByGlobal) {if (!patMuon.isGlobalMuon()) {return false;}}
  // ----Acceptance cuts
  if(smear) {
    if (fabs(smearedLV.eta())>_RecoMuonEtaCut) {return false;}
    if (smearedLV.pt()<_RecoMuonPtMinCut) {return false;}
    if (smearedLV.pt()>_RecoMuonPtMaxCut) {return false;}
  } else {
    if (fabs(patMuon.eta())>_RecoMuonEtaCut) {return false;}
    if (patMuon.pt()<_RecoMuonPtMinCut) {return false;}
    if (patMuon.pt()>_RecoMuonPtMaxCut) {return false;}
  }
/*
  if (_DoRecoMuonDiscrByTrackIsolation) {
    if ( (patMuon.trackIsoDeposit()->depositAndCountWithin(_RecoMuonTrackIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonTrackIsoTrkThreshold).first >= _RecoMuonTrackIsoSumPtCutValue) ) {return false;}
  }
  if (_DoRecoMuonDiscrByEcalIsolation) {
    if ( (patMuon.ecalIsoDeposit()->depositAndCountWithin(_RecoMuonEcalIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonEcalIsoRecHitThreshold).first >= _RecoMuonEcalIsoSumPtCutValue) ) {return false;}
  }
*/
  // ----Isolation requirement
  if (_DoRecoMuonDiscrByIsolation) {
    if( (patMuon.trackIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonTrackIsoTrkThreshold).first +
         patMuon.ecalIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonEcalIsoRecHitThreshold).first) 
         >= _RecoMuonIsoSumPtMaxCutValue) {return false;}
    if( (patMuon.trackIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonTrackIsoTrkThreshold).first +
         patMuon.ecalIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonEcalIsoRecHitThreshold).first) 
         < _RecoMuonIsoSumPtMinCutValue) {return false;}
/*
    if( (patMuon.trackIso() + patMuon.ecalIso()) 
	 >= _RecoMuonIsoSumPtMaxCutValue) {return false;}
    if( (patMuon.trackIso() + patMuon.ecalIso()) 
	 < _RecoMuonIsoSumPtMinCutValue) {return false;}
*/
  }
  // ----Impact parameter requirement
  if (_DoRecoMuonDiscrByIp) {
    const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
    if ( patMuon.track().isNonnull() ) {if ( (fabs(patMuon.track()->dxy(thePrimaryEventVertex.position())) >= _RecoMuonIpCut) ) {return false;}}
    else {return false;}
  }
  return true;
}

//--------------Apply Electron Cuts
bool HiMassTauAnalysis::passRecoElectronCuts(const pat::Electron& patElectron,bool smear,reco::Candidate::LorentzVector smearedLV) {
  // ----Matching to gen
  if(_MatchLeptonToGen) {if(!(matchToGen(patElectron).first)) {return false;}}
  // ----Acceptance cuts
  if(smear) {
    if (fabs(smearedLV.eta())>_RecoElectronEtaCut) {return false;}
    if (smearedLV.pt()<_RecoElectronPtMinCut) {return false;}
    if (smearedLV.pt()>_RecoElectronPtMaxCut) {return false;}
  } else {
    if (fabs(patElectron.eta())>_RecoElectronEtaCut) {return false;}
    if (patElectron.pt()<_RecoElectronPtMinCut) {return false;}
    if (patElectron.pt()>_RecoElectronPtMaxCut) {return false;}
  }
  // ----Isolation requirement
  if (_DoRecoElectronDiscrByTrackIsolation) {
//    if ( (patElectron.trackIsoDeposit()->depositAndCountWithin(_RecoElectronTrackIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoElectronTrackIsoTrkThreshold).first >= _RecoElectronTrackIsoSumPtCutValue) ) {return false;}
    if(patElectron.trackIso() > _RecoElectronTrackIsoSumPtCutValue) {return false;}
  }
  if (_DoRecoElectronDiscrByEcalIsolation) {
//    if ( (patElectron.ecalIsoDeposit()->depositAndCountWithin(_RecoElectronEcalIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoElectronEcalIsoRecHitThreshold).first >= _RecoElectronEcalIsoSumPtCutValue) ) {return false;}
    if(patElectron.ecalIso() > _RecoElectronEcalIsoSumPtCutValue) {return false;}
  }
  // ----Impact parameter requirement
  if (_DoRecoElectronDiscrByIp) {
    const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
//    if ( patElectron.gsfTrack().isNonnull() ) {if ( (fabs(patElectron.gsfTrack()->dxy(thePrimaryEventVertex.position())) >= _RecoElectronIpCut) ) {return false;}}
    if ( patElectron.track().isNonnull() ) {if ( (fabs(patElectron.track()->dxy(thePrimaryEventVertex.position())) >= _RecoElectronIpCut) ) {return false;}}
    else {return false;}
  }
  // ----E over p requirement
  if (_DoRecoElectronDiscrByEoverP) { if((patElectron.eSuperClusterOverP()>_RecoElectronEoverPMax) || (patElectron.eSuperClusterOverP()<_RecoElectronEoverPMin)) {return false;} }
  // ----Electromagnetic energy fraction requirement
  if (_DoRecoElectronDiscrByHoverEm) { if(patElectron.hadronicOverEm() > _RecoElectronHoverEmCut) {return false;} }
  if (_DoRecoElectronDiscrBySigmaEtaEta) {}
  if (_DoRecoElectronDiscrBySigmaIEtaIEta) {}
  if (_DoRecoElectronDiscrBySCE5by5) {}
  return true;
}

//--------------Apply Jet Cuts
bool HiMassTauAnalysis::passRecoJetCuts(const pat::Jet& patJet) {
  // use corrected jet?
  if(_UseCorrectedJet) {
    // ----Acceptance cuts
    if (fabs(patJet.eta())>_RecoJetEtaMaxCut) {return false;}
    if (fabs(patJet.eta())<_RecoJetEtaMinCut) {return false;}
    if (patJet.pt()<_RecoJetPtCut) {return false;}
    // ----anti-overlap requirements
    if (_RemoveJetOverlapWithMuons) {
      int theNumberOfMuons = 0;
      for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
        theNumberOfMuons++;
        if( (passRecoMuonCuts((*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1))) ) {
          if(_SmearTheMuon) {
            if(reco::deltaR(patJet.p4(), smearedMuonMomentumVector.at(theNumberOfMuons-1)) < _JetMuonMatchingDeltaR) {return false;}
          } else {
            if(reco::deltaR(patJet.p4(), patMuon->p4()) < _JetMuonMatchingDeltaR) {return false;}
          }
        }
      }
    }
    if (_RemoveJetOverlapWithElectrons) {
      int theNumberOfElectrons = 0;
      for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron ) {
        theNumberOfElectrons++;
        if( (passRecoElectronCuts((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1))) ) {
          if(_SmearTheElectron) {
            if(reco::deltaR(patJet.p4(), smearedElectronMomentumVector.at(theNumberOfElectrons-1)) < _JetElectronMatchingDeltaR) {return false;}
          } else {
            if(reco::deltaR(patJet.p4(), patElectron->p4()) < _JetElectronMatchingDeltaR) {return false;}
          }
        }
      }
    }
    if (_RemoveJetOverlapWithTaus) {
      for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
        if( (passRecoTauCuts(*patTau)) && (reco::deltaR(patJet.p4(), patTau->p4()) < _JetTauMatchingDeltaR) ) {return false;}
      }
    }
  } else {
    // ----Acceptance cuts
    if (fabs(patJet.correctedJet("raw","").eta())>_RecoJetEtaMaxCut) {return false;}
    if (fabs(patJet.correctedJet("raw","").eta())<_RecoJetEtaMinCut) {return false;}
    if (patJet.correctedJet("raw","").pt()<_RecoJetPtCut) {return false;}
    // ----anti-overlap requirements
    if (_RemoveJetOverlapWithMuons) {
      int theNumberOfMuons = 0;
      for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
        theNumberOfMuons++;
        if( (passRecoMuonCuts((*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1))) ) {
          if(_SmearTheMuon) {
            if(reco::deltaR(patJet.correctedJet("raw","").p4(), smearedMuonMomentumVector.at(theNumberOfMuons-1)) < _JetMuonMatchingDeltaR) {return false;}
          } else {
            if(reco::deltaR(patJet.correctedJet("raw","").p4(), patMuon->p4()) < _JetMuonMatchingDeltaR) {return false;}
          }
        }
      }
    }
    if (_RemoveJetOverlapWithElectrons) {
      int theNumberOfElectrons = 0;
      for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron ) {
        theNumberOfElectrons++;
        if( (passRecoElectronCuts((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1))) ) {
          if(_SmearTheElectron) {
            if(reco::deltaR(patJet.correctedJet("raw","").p4(), smearedElectronMomentumVector.at(theNumberOfElectrons-1)) < _JetElectronMatchingDeltaR) {return false;}
          } else {
            if(reco::deltaR(patJet.correctedJet("raw","").p4(), patElectron->p4()) < _JetElectronMatchingDeltaR) {return false;}
          }
        }
      }
    }
    if (_RemoveJetOverlapWithTaus) {
      for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
        if( (passRecoTauCuts(*patTau)) && (reco::deltaR(patJet.correctedJet("raw","").p4(), patTau->p4()) < _JetTauMatchingDeltaR) ) {return false;}
      }
    }
  }
  return true;
}

// ---------------Apply Topology Cuts
bool HiMassTauAnalysis::passTopologyCuts(const pat::Tau& patTau, const pat::Muon& patMuon, bool smear, reco::Candidate::LorentzVector smearedLV, const pat::MET& patMET) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if(smear) {
    if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(patTau.p4(), smearedLV) < _DiTauDeltaRCut) {return false;} }
  } else {
    if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(patTau.p4(), patMuon.p4()) < _DiTauDeltaRCut) {return false;} }
  }
  // ----Opposite sign - Like sign requirement
  if (_DiTauDiscrByOSLSType == "OS") {
    if(_UseTauSeedTrackForDiTauDiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patMuon.charge() * patTau.leadTrack()->charge()) >= 0) {return false;} }
      else { if((patMuon.charge() * patTau.leadPFChargedHadrCand()->charge()) >= 0) {return false;} }
    } else { if((patMuon.charge() * patTau.charge()) >= 0) {return false;} }
  } else if (_DiTauDiscrByOSLSType == "LS") {
    if(_UseTauSeedTrackForDiTauDiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patMuon.charge() * patTau.leadTrack()->charge()) <= 0) {return false;} }
      else { if((patMuon.charge() * patTau.leadPFChargedHadrCand()->charge()) <= 0) {return false;} }
    } else { if((patMuon.charge() * patTau.charge()) <= 0) {return false;} }
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if(smear) {
    if (_DoDiTauDiscrByCosDphi) {
      if(cos(TMath::Abs(normalizedPhi(smearedLV.phi() - patTau.phi()))) > _DiTauCosDphiMaxCut) {return false;}
      if(cos(TMath::Abs(normalizedPhi(smearedLV.phi() - patTau.phi()))) < _DiTauCosDphiMinCut) {return false;}
    }
  } else {
    if (_DoDiTauDiscrByCosDphi) {
      if(cos(TMath::Abs(normalizedPhi(patMuon.phi() - patTau.phi()))) > _DiTauCosDphiMaxCut) {return false;}
      if(cos(TMath::Abs(normalizedPhi(patMuon.phi() - patTau.phi()))) < _DiTauCosDphiMinCut) {return false;}
    }
  }
  // ----Mass window requirement
  if (_DoDiscrByMassReco) { if( (CalculateThe4Momentum(patTau,patMuon,smear,smearedLV,patMET).second.M() < _MassMinCut) || (CalculateThe4Momentum(patTau,patMuon,smear,smearedLV,patMET).second.M() > _MassMaxCut) ) {return false;} }
  // ----Zeta requirement
  if (_DoDiTauDiscrByCDFzeta2D) {
    if( ((_PZetaCutCoefficient * CalculatePZeta(patTau,patMuon,smear,smearedLV,patMET)) + 
         (_PZetaVisCutCoefficient * CalculatePZetaVis(patTau,patMuon,smear,smearedLV,patMET))) < _CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) { if (patMET.pt()<_RecoMetCut) {return false;} }
  if(smear) {
    if (_DoDiTauDiscrByDeltaPtDivSumPt) {
      if ( ((patTau.pt() - smearedLV.pt()) / (patTau.pt() + smearedLV.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
      if ( ((patTau.pt() - smearedLV.pt()) / (patTau.pt() + smearedLV.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
    }
    if (_DoDiscrByLeg1MetDphi) {
      if(_AnalyzeMuonForLeg1) {
        if(TMath::Abs(normalizedPhi(smearedLV.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(smearedLV.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
      }
      if(_AnalyzeTauForLeg1) {
        if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
      }
    }
    if (_DoDiscrByLeg2MetDphi) {
      if(_AnalyzeMuonForLeg2) {
        if(TMath::Abs(normalizedPhi(smearedLV.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(smearedLV.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
      }
      if(_AnalyzeTauForLeg2) {
        if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
      }
    }
  } else {
    if (_DoDiTauDiscrByDeltaPtDivSumPt) {
      if ( ((patTau.pt() - patMuon.pt()) / (patTau.pt() + patMuon.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
      if ( ((patTau.pt() - patMuon.pt()) / (patTau.pt() + patMuon.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
    }
    if (_DoDiscrByLeg1MetDphi) {
      if(_AnalyzeMuonForLeg1) {
        if(TMath::Abs(normalizedPhi(patMuon.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patMuon.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
      }
      if(_AnalyzeTauForLeg1) {
        if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
      }
    }
    if (_DoDiscrByLeg2MetDphi) {
      if(_AnalyzeMuonForLeg2) {
        if(TMath::Abs(normalizedPhi(patMuon.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patMuon.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
      }
      if(_AnalyzeTauForLeg2) {
        if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
      }
    }
  }
  return true;
}
bool HiMassTauAnalysis::passTopologyCuts(const pat::Tau& patTau, const pat::Electron& patElectron, bool smear, reco::Candidate::LorentzVector smearedLV, const pat::MET& patMET) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if(smear) {
    if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(patTau.p4(), smearedLV) < _DiTauDeltaRCut) {return false;} }
  } else {
    if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(patTau.p4(), patElectron.p4()) < _DiTauDeltaRCut) {return false;} }
  }
  // ----Opposite sign - Like sign requirement
  if (_DiTauDiscrByOSLSType == "OS") {
    if(_UseTauSeedTrackForDiTauDiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patElectron.charge() * patTau.leadTrack()->charge()) >= 0) {return false;} }
      else { if((patElectron.charge() * patTau.leadPFChargedHadrCand()->charge()) >= 0) {return false;} }
    } else { if((patElectron.charge() * patTau.charge()) >= 0) {return false;} }
  } else if (_DiTauDiscrByOSLSType == "LS") {
    if(_UseTauSeedTrackForDiTauDiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patElectron.charge() * patTau.leadTrack()->charge()) <= 0) {return false;} }
      else { if((patElectron.charge() * patTau.leadPFChargedHadrCand()->charge()) <= 0) {return false;} }
    } else { if((patElectron.charge() * patTau.charge()) <= 0) {return false;} }
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if(smear) {
    if (_DoDiTauDiscrByCosDphi) {
      if(cos(TMath::Abs(normalizedPhi(smearedLV.phi() - patTau.phi()))) > _DiTauCosDphiMaxCut) {return false;}
      if(cos(TMath::Abs(normalizedPhi(smearedLV.phi() - patTau.phi()))) < _DiTauCosDphiMinCut) {return false;}
    }
  } else {
    if (_DoDiTauDiscrByCosDphi) {
      if(cos(TMath::Abs(normalizedPhi(patElectron.phi() - patTau.phi()))) > _DiTauCosDphiMaxCut) {return false;}
      if(cos(TMath::Abs(normalizedPhi(patElectron.phi() - patTau.phi()))) < _DiTauCosDphiMinCut) {return false;}
    }
  }
  // ----Mass window requirement
  if (_DoDiscrByMassReco) { if( (CalculateThe4Momentum(patTau,patElectron,smear,smearedLV,patMET).second.M() < _MassMinCut) || (CalculateThe4Momentum(patTau,patElectron,smear,smearedLV,patMET).second.M() > _MassMaxCut) ) {return false;} }
  // ----Zeta requirement
  if (_DoDiTauDiscrByCDFzeta2D) {
    if( ((_PZetaCutCoefficient * CalculatePZeta(patTau,patElectron,smear,smearedLV,patMET)) +
         (_PZetaVisCutCoefficient * CalculatePZetaVis(patTau,patElectron,smear,smearedLV,patMET))) < _CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) { if (patMET.pt()<_RecoMetCut) {return false;} }
  if(smear) {
    if (_DoDiTauDiscrByDeltaPtDivSumPt) {
      if ( ((patTau.pt() - smearedLV.pt()) / (patTau.pt() + smearedLV.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
      if ( ((patTau.pt() - smearedLV.pt()) / (patTau.pt() + smearedLV.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
    }
    if (_DoDiscrByLeg1MetDphi) {
      if(_AnalyzeElectronForLeg1) {
        if(TMath::Abs(normalizedPhi(smearedLV.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(smearedLV.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
      }
      if(_AnalyzeTauForLeg1) {
        if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
      }
    }
    if (_DoDiscrByLeg2MetDphi) {
      if(_AnalyzeElectronForLeg2) {
        if(TMath::Abs(normalizedPhi(smearedLV.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(smearedLV.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
      }
      if(_AnalyzeTauForLeg2) {
        if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
      }
    }
  } else {
    if (_DoDiTauDiscrByDeltaPtDivSumPt) {
      if ( ((patTau.pt() - patElectron.pt()) / (patTau.pt() + patElectron.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
      if ( ((patTau.pt() - patElectron.pt()) / (patTau.pt() + patElectron.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
    }
    if (_DoDiscrByLeg1MetDphi) {
      if(_AnalyzeElectronForLeg1) {
        if(TMath::Abs(normalizedPhi(patElectron.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patElectron.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
      }
      if(_AnalyzeTauForLeg1) {
        if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
      }
    }
    if (_DoDiscrByLeg2MetDphi) {
      if(_AnalyzeElectronForLeg2) {
        if(TMath::Abs(normalizedPhi(patElectron.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patElectron.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
      }
      if(_AnalyzeTauForLeg2) {
        if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
      }
    }
  }
  return true;
}
bool HiMassTauAnalysis::passTopologyCuts(const pat::Electron& patElectron, bool smearE, reco::Candidate::LorentzVector smearedLVE, const pat::Muon& patMuon, bool smearM, reco::Candidate::LorentzVector smearedLVM, const pat::MET& patMET) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if(smearE) {
    if(smearM) {
      if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(smearedLVM, smearedLVE) < _DiTauDeltaRCut) {return false;} }
    } else {
      if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(patMuon.p4(), smearedLVE) < _DiTauDeltaRCut) {return false;} }
    }
  } else {
    if(smearM) {
      if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(smearedLVM, patElectron.p4()) < _DiTauDeltaRCut) {return false;} }
    } else {
      if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(patMuon.p4(), patElectron.p4()) < _DiTauDeltaRCut) {return false;} }
    }
  }
  // ----Opposite sign - Like sign requirement
  if (_DiTauDiscrByOSLSType == "OS") {
    if((patElectron.charge() * patMuon.charge()) >= 0) {return false;}
  } else if (_DiTauDiscrByOSLSType == "LS") {
    if((patElectron.charge() * patMuon.charge()) <= 0) {return false;}
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if(smearE) {
    if(smearM) {
      if (_DoDiTauDiscrByCosDphi) {
        if(cos(TMath::Abs(normalizedPhi(smearedLVE.phi() - smearedLVM.phi()))) > _DiTauCosDphiMaxCut) {return false;}
        if(cos(TMath::Abs(normalizedPhi(smearedLVE.phi() - smearedLVM.phi()))) < _DiTauCosDphiMinCut) {return false;}
      }
    } else {
      if (_DoDiTauDiscrByCosDphi) {
        if(cos(TMath::Abs(normalizedPhi(smearedLVE.phi() - patMuon.phi()))) > _DiTauCosDphiMaxCut) {return false;}
        if(cos(TMath::Abs(normalizedPhi(smearedLVE.phi() - patMuon.phi()))) < _DiTauCosDphiMinCut) {return false;}
      }
    }
  } else {
    if(smearM) {
      if (_DoDiTauDiscrByCosDphi) {
        if(cos(TMath::Abs(normalizedPhi(patElectron.phi() - smearedLVM.phi()))) > _DiTauCosDphiMaxCut) {return false;}
        if(cos(TMath::Abs(normalizedPhi(patElectron.phi() - smearedLVM.phi()))) < _DiTauCosDphiMinCut) {return false;}
      }
    } else {
      if (_DoDiTauDiscrByCosDphi) {
        if(cos(TMath::Abs(normalizedPhi(patElectron.phi() - patMuon.phi()))) > _DiTauCosDphiMaxCut) {return false;}
        if(cos(TMath::Abs(normalizedPhi(patElectron.phi() - patMuon.phi()))) < _DiTauCosDphiMinCut) {return false;}
      }
    }
  }
  // ----Mass window requirement
  if (_DoDiscrByMassReco) { if( (CalculateThe4Momentum(patElectron,smearE,smearedLVE,patMuon,smearM,smearedLVM,patMET).second.M() < _MassMinCut) || (CalculateThe4Momentum(patElectron,smearE,smearedLVE,patMuon,smearM,smearedLVM,patMET).second.M() > _MassMaxCut) ) {return false;} }
  // ----Zeta requirement
  if (_DoDiTauDiscrByCDFzeta2D) {
    if( ((_PZetaCutCoefficient * CalculatePZeta(patElectron,smearE,smearedLVE,patMuon,smearM,smearedLVM,patMET)) +
         (_PZetaVisCutCoefficient * CalculatePZetaVis(patElectron,smearE,smearedLVE,patMuon,smearM,smearedLVM,patMET))) < _CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) { if (patMET.pt()<_RecoMetCut) {return false;} }
  if(smearE) {
    if(smearM) {
      if (_DoDiTauDiscrByDeltaPtDivSumPt) {
        if ( ((smearedLVM.pt() - smearedLVE.pt()) / (smearedLVM.pt() + smearedLVE.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
        if ( ((smearedLVM.pt() - smearedLVE.pt()) / (smearedLVM.pt() + smearedLVE.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
      }
      if (_DoDiscrByLeg1MetDphi) {
        if(_AnalyzeElectronForLeg1) {
          if(TMath::Abs(normalizedPhi(smearedLVE.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVE.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
        if(_AnalyzeMuonForLeg1) {
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
      }
      if (_DoDiscrByLeg2MetDphi) {
        if(_AnalyzeElectronForLeg2) {
          if(TMath::Abs(normalizedPhi(smearedLVE.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVE.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
        }
        if(_AnalyzeMuonForLeg2) {
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
        }
      }
    } else {
      if (_DoDiTauDiscrByDeltaPtDivSumPt) {
        if ( ((patMuon.pt() - smearedLVE.pt()) / (patMuon.pt() + smearedLVE.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
        if ( ((patMuon.pt() - smearedLVE.pt()) / (patMuon.pt() + smearedLVE.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
      }
      if (_DoDiscrByLeg1MetDphi) {
        if(_AnalyzeElectronForLeg1) {
          if(TMath::Abs(normalizedPhi(smearedLVE.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVE.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
        if(_AnalyzeMuonForLeg1) {
          if(TMath::Abs(normalizedPhi(patMuon.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(patMuon.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
      }
      if (_DoDiscrByLeg2MetDphi) {
        if(_AnalyzeElectronForLeg2) {
          if(TMath::Abs(normalizedPhi(smearedLVE.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVE.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
        }
        if(_AnalyzeMuonForLeg2) {
          if(TMath::Abs(normalizedPhi(patMuon.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(patMuon.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
        }
      }
    }
  } else {
    if(smearM) {
      if (_DoDiTauDiscrByDeltaPtDivSumPt) {
        if ( ((smearedLVM.pt() - patElectron.pt()) / (smearedLVM.pt() + patElectron.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
        if ( ((smearedLVM.pt() - patElectron.pt()) / (smearedLVM.pt() + patElectron.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
      }
      if (_DoDiscrByLeg1MetDphi) {
        if(_AnalyzeElectronForLeg1) {
          if(TMath::Abs(normalizedPhi(patElectron.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(patElectron.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
        if(_AnalyzeMuonForLeg1) {
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
      }
      if (_DoDiscrByLeg2MetDphi) {
        if(_AnalyzeElectronForLeg2) {
          if(TMath::Abs(normalizedPhi(patElectron.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(patElectron.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
        }
        if(_AnalyzeMuonForLeg2) {
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
        }
      }
    } else {
      if (_DoDiTauDiscrByDeltaPtDivSumPt) {
        if ( ((patMuon.pt() - patElectron.pt()) / (patMuon.pt() + patElectron.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
        if ( ((patMuon.pt() - patElectron.pt()) / (patMuon.pt() + patElectron.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
      }
      if (_DoDiscrByLeg1MetDphi) {
        if(_AnalyzeElectronForLeg1) {
          if(TMath::Abs(normalizedPhi(patElectron.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(patElectron.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
        if(_AnalyzeMuonForLeg1) {
          if(TMath::Abs(normalizedPhi(patMuon.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(patMuon.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
      }
      if (_DoDiscrByLeg2MetDphi) {
        if(_AnalyzeElectronForLeg2) {
          if(TMath::Abs(normalizedPhi(patElectron.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(patElectron.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
        }
        if(_AnalyzeMuonForLeg2) {
          if(TMath::Abs(normalizedPhi(patMuon.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(patMuon.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
        }
      }
    }
  }
  return true;
}
bool HiMassTauAnalysis::passTopologyCuts(const pat::Muon& patMuon1, bool smear1, reco::Candidate::LorentzVector smearedLV1, const pat::Muon& patMuon2, bool smear2, reco::Candidate::LorentzVector smearedLV2, const pat::MET& patMET) {
  // ----Separation cut between muon legs (remove double counting)
  if(smear1) {
    if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(smearedLV1, smearedLV2) < _DiTauDeltaRCut) {return false;} }
  } else {
    if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(patMuon1.p4(), patMuon2.p4()) < _DiTauDeltaRCut) {return false;} }
  }
  // ----Opposite sign - Like sign requirement
  if (_DiTauDiscrByOSLSType == "OS") {
    if((patMuon1.charge() * patMuon2.charge()) >= 0) {return false;}
  } else if (_DiTauDiscrByOSLSType == "LS") {
    if((patMuon1.charge() * patMuon2.charge()) <= 0) {return false;}
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if(smear1) {
    if (_DoDiTauDiscrByCosDphi) {
      if(cos(TMath::Abs(normalizedPhi(smearedLV1.phi() - smearedLV2.phi()))) > _DiTauCosDphiMaxCut) {return false;}
      if(cos(TMath::Abs(normalizedPhi(smearedLV1.phi() - smearedLV2.phi()))) < _DiTauCosDphiMinCut) {return false;}
    }
  } else {
    if (_DoDiTauDiscrByCosDphi) {
      if(cos(TMath::Abs(normalizedPhi(patMuon1.phi() - patMuon2.phi()))) > _DiTauCosDphiMaxCut) {return false;}
      if(cos(TMath::Abs(normalizedPhi(patMuon1.phi() - patMuon2.phi()))) < _DiTauCosDphiMinCut) {return false;}
    }
  }
  // ----Mass window requirement
  if (_DoDiscrByMassReco) { if( (CalculateThe4Momentum(patMuon1,smear1,smearedLV1,patMuon2,smear2,smearedLV2,patMET).second.M() < _MassMinCut) || (CalculateThe4Momentum(patMuon1,smear1,smearedLV1,patMuon2,smear2,smearedLV2,patMET).second.M() > _MassMaxCut) ) {return false;} }
  // ----Zeta requirement
  if (_DoDiTauDiscrByCDFzeta2D) {
    if( ((_PZetaCutCoefficient * CalculatePZeta(patMuon1,smear1,smearedLV1,patMuon2,smear2,smearedLV2,patMET)) +
         (_PZetaVisCutCoefficient * CalculatePZetaVis(patMuon1,smear1,smearedLV1,patMuon2,smear2,smearedLV2,patMET))) < _CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) { if (patMET.pt()<_RecoMetCut) {return false;} }
  if(smear1) {
    if (_DoDiTauDiscrByDeltaPtDivSumPt) {
      if ( ((smearedLV1.pt() - smearedLV2.pt()) / (smearedLV1.pt() + smearedLV2.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
      if ( ((smearedLV1.pt() - smearedLV2.pt()) / (smearedLV1.pt() + smearedLV2.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
    }
    if (_DoDiscrByLeg1MetDphi) {
      if(_AnalyzeMuonForLeg1) {
        if(TMath::Abs(normalizedPhi(smearedLV1.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(smearedLV1.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
      }
    }
    if (_DoDiscrByLeg2MetDphi) {
      if(_AnalyzeMuonForLeg2) {
        if(TMath::Abs(normalizedPhi(smearedLV2.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(smearedLV2.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
      }
    }
  } else {
    if (_DoDiTauDiscrByDeltaPtDivSumPt) {
      if ( ((patMuon1.pt() - patMuon2.pt()) / (patMuon1.pt() + patMuon2.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
      if ( ((patMuon1.pt() - patMuon2.pt()) / (patMuon1.pt() + patMuon2.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
    }
    if (_DoDiscrByLeg1MetDphi) {
      if(_AnalyzeMuonForLeg1) {
        if(TMath::Abs(normalizedPhi(patMuon1.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patMuon1.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
      }
    }
    if (_DoDiscrByLeg2MetDphi) {
      if(_AnalyzeMuonForLeg2) {
        if(TMath::Abs(normalizedPhi(patMuon2.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patMuon2.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
      }
    }
  }
  return true;
}
bool HiMassTauAnalysis::passTopologyCuts(const pat::Electron& patElectron1, bool smear1, reco::Candidate::LorentzVector smearedLV1, const pat::Electron& patElectron2, bool smear2, reco::Candidate::LorentzVector smearedLV2, const pat::MET& patMET) {
  // ----Separation cut between electron legs (remove double counting)
  if(smear1) {
    if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(smearedLV1, smearedLV2) < _DiTauDeltaRCut) {return false;} }
  } else {
    if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(patElectron1.p4(), patElectron2.p4()) < _DiTauDeltaRCut) {return false;} }
  }
  // ----Opposite sign - Like sign requirement
  if (_DiTauDiscrByOSLSType == "OS") {
    if((patElectron1.charge() * patElectron2.charge()) >= 0) {return false;}
  } else if (_DiTauDiscrByOSLSType == "LS") {
    if((patElectron1.charge() * patElectron2.charge()) <= 0) {return false;}
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if(smear1) {
    if (_DoDiTauDiscrByCosDphi) {
      if(cos(TMath::Abs(normalizedPhi(smearedLV1.phi() - smearedLV2.phi()))) > _DiTauCosDphiMaxCut) {return false;}
      if(cos(TMath::Abs(normalizedPhi(smearedLV1.phi() - smearedLV2.phi()))) < _DiTauCosDphiMinCut) {return false;}
    }
  } else {
    if (_DoDiTauDiscrByCosDphi) {
      if(cos(TMath::Abs(normalizedPhi(patElectron1.phi() - patElectron2.phi()))) > _DiTauCosDphiMaxCut) {return false;}
      if(cos(TMath::Abs(normalizedPhi(patElectron1.phi() - patElectron2.phi()))) < _DiTauCosDphiMinCut) {return false;}
    }
  }
  // ----Mass window requirement
  if (_DoDiscrByMassReco) { if( (CalculateThe4Momentum(patElectron1,smear1,smearedLV1,patElectron2,smear2,smearedLV2,patMET).second.M() < _MassMinCut) || (CalculateThe4Momentum(patElectron1,smear1,smearedLV1,patElectron2,smear2,smearedLV2,patMET).second.M() > _MassMaxCut) ) {return false;} }
  // ----Zeta requirement
  if (_DoDiTauDiscrByCDFzeta2D) {
    if( ((_PZetaCutCoefficient * CalculatePZeta(patElectron1,smear1,smearedLV1,patElectron2,smear2,smearedLV2,patMET)) +
         (_PZetaVisCutCoefficient * CalculatePZetaVis(patElectron1,smear1,smearedLV1,patElectron2,smear2,smearedLV2,patMET))) < _CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) { if (patMET.pt()<_RecoMetCut) {return false;} }
  if(smear1) {
    if (_DoDiTauDiscrByDeltaPtDivSumPt) {
      if ( ((smearedLV1.pt() - smearedLV2.pt()) / (smearedLV1.pt() + smearedLV2.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
      if ( ((smearedLV1.pt() - smearedLV2.pt()) / (smearedLV1.pt() + smearedLV2.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
    }
    if (_DoDiscrByLeg1MetDphi) {
      if(_AnalyzeElectronForLeg1) {
        if(TMath::Abs(normalizedPhi(smearedLV1.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(smearedLV1.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
      }
    }
    if (_DoDiscrByLeg2MetDphi) {
      if(_AnalyzeElectronForLeg2) {
        if(TMath::Abs(normalizedPhi(smearedLV2.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(smearedLV2.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
      }
    }
  } else {
    if (_DoDiTauDiscrByDeltaPtDivSumPt) {
      if ( ((patElectron1.pt() - patElectron2.pt()) / (patElectron1.pt() + patElectron2.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
      if ( ((patElectron1.pt() - patElectron2.pt()) / (patElectron1.pt() + patElectron2.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
    }
    if (_DoDiscrByLeg1MetDphi) {
      if(_AnalyzeElectronForLeg1) {
        if(TMath::Abs(normalizedPhi(patElectron1.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patElectron1.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
      }
    }
    if (_DoDiscrByLeg2MetDphi) {
      if(_AnalyzeElectronForLeg2) {
        if(TMath::Abs(normalizedPhi(patElectron2.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patElectron2.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
      }
    }
  }
  return true;
}
bool HiMassTauAnalysis::passTopologyCuts(const pat::Tau& patTau1, const pat::Tau& patTau2, const pat::MET& patMET) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(patTau1.p4(), patTau2.p4()) < _DiTauDeltaRCut) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_DiTauDiscrByOSLSType == "OS") {
    if(_UseTauSeedTrackForDiTauDiscrByOSLS) {
      if (patTau1.isCaloTau()) { if((patTau1.leadTrack()->charge() * patTau2.leadTrack()->charge()) >= 0) {return false;} }
      else { if((patTau1.leadPFChargedHadrCand()->charge() * patTau2.leadPFChargedHadrCand()->charge()) >= 0) {return false;} }
    } else { if((patTau1.charge() * patTau2.charge()) >= 0) {return false;} }
  } else if (_DiTauDiscrByOSLSType == "LS") {
    if(_UseTauSeedTrackForDiTauDiscrByOSLS) {
      if (patTau1.isCaloTau()) { if((patTau1.leadTrack()->charge() * patTau2.leadTrack()->charge()) <= 0) {return false;} }
      else { if((patTau1.leadPFChargedHadrCand()->charge() * patTau2.leadPFChargedHadrCand()->charge()) <= 0) {return false;} }
    } else { if((patTau1.charge() * patTau2.charge()) <= 0) {return false;} }
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoDiTauDiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(patTau1.phi() - patTau2.phi()))) > _DiTauCosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(patTau1.phi() - patTau2.phi()))) < _DiTauCosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByMassReco) { if( (CalculateThe4Momentum(patTau1,patTau2,patMET).second.M() < _MassMinCut) || (CalculateThe4Momentum(patTau1,patTau2,patMET).second.M() > _MassMaxCut) ) {return false;} }
  // ----Zeta requirement
  if (_DoDiTauDiscrByCDFzeta2D) {
    if( ((_PZetaCutCoefficient * CalculatePZeta(patTau1, patTau2, patMET)) +
         (_PZetaVisCutCoefficient * CalculatePZetaVis(patTau1, patTau2, patMET))) < _CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) { if (patMET.pt()<_RecoMetCut) {return false;} }
  if (_DoDiTauDiscrByDeltaPtDivSumPt) {
    if ( ((patTau1.pt() - patTau2.pt()) / (patTau1.pt() + patTau2.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((patTau1.pt() - patTau2.pt()) / (patTau1.pt() + patTau2.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByLeg1MetDphi) {
    if(_AnalyzeTauForLeg1) {
      if(TMath::Abs(normalizedPhi(patTau1.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(patTau1.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
    }
  }
  if (_DoDiscrByLeg2MetDphi) {
    if(_AnalyzeTauForLeg2) {
      if(TMath::Abs(normalizedPhi(patTau2.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(patTau2.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
    }
  }
  return true;
}

// ---------------Fill Ntuple
void HiMassTauAnalysis::fillNtuple() {

  for(CompositeCandidateCollection::const_iterator theCand = _patDiTaus->begin(); theCand != _patDiTaus->end(); ++theCand){
    pair<const pat::Tau*, const pat::Electron*> theCandDaughters = getPATComponents(*theCand);
    const pat::Tau* theTau = theCandDaughters.first;
    const pat::Electron* theElectron = theCandDaughters.second;    
    if(fabs(theElectron->eta()) > _RecoElectronEtaCut || theElectron->pt() < _RecoElectronPtMinCut) continue;
    if(fabs(theTau->eta()) > _RecoTauEtaCut || theTau->tauID(_RecoTauDiscrByLeadTrack.data()) < 0.5) continue;
    _diTauMass->push_back(theCand->mass());
  }

  int theNumberOfElectrons = 0;    
  for(pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron) {
    theNumberOfElectrons++;
    if(fabs(patElectron->eta()) > _RecoElectronEtaCut || patElectron->pt() < _RecoElectronPtMinCut) continue;
    for(pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau) {
      if(reco::deltaR(patTau->p4(), patElectron->p4()) < _DiTauDeltaRCut) continue; 				     
      if(patTau->tauID(_RecoTauDiscrByLeadTrack.data()) < 0.5) continue;
            
      const GenParticleRef theGenTauRef = patTau->genParticleRef();
      _tauMatched->push_back(int(matchToGen(*patTau).first));							       
      if(theGenTauRef.isNonnull()){
        _tauMotherId->push_back(abs(theGenTauRef->mother()->mother()->pdgId()));
      }			       
      else{
        _tauMotherId->push_back(0);
      }												       
      _tauE->push_back(patTau->energy()); 											       
      _tauEt->push_back(patTau->et());												       
      _tauPt->push_back(patTau->pt());												       
      _tauCharge->push_back(patTau->charge());											       
      _tauEta->push_back(patTau->eta());  											       
      _tauPhi->push_back(patTau->phi());
      
      TrackRef theLeadTrack = patTau->leadTrack();
      if(theLeadTrack.isNonnull()){
         _tauLTPt->push_back(theLeadTrack->pt());
         _tauLTCharge->push_back(theLeadTrack->charge());
         //_tauLTChi2->push_back(theLeadTrack->chi2());
         //_tauLTRecHitsSize->push_back(theLeadTrack->recHitsSize());
      }
      else {
        _tauLTPt->push_back(-1.);
        _tauLTCharge->push_back(-10);
        //_tauLTChi2->push_back(100.);
        //_tauLTRecHitsSize->push_back(0);
      } 									
      
      if (patTau->isCaloTau()){
        _tauNProngs->push_back(patTau->signalTracks().size());									       
        _tauEmFraction->push_back(-1.);												       
      	_tauHcalTotOverPLead->push_back(-1.);		  									    
      	_tauHcalMaxOverPLead->push_back(-1.);		  									    
      	_tauHcal3x3OverPLead->push_back(-1.);		  									    
      	_tauElectronPreId->push_back(0);		  									    
      	_tauModifiedEOverP->push_back(-1.);		  									    
      	_tauBremsRecoveryEOverPLead->push_back(-1.);
      	_tauDiscAgainstElec->push_back(-1.);
	_tauLTSignedIp->push_back(patTau->leadTracksignedSipt());
      } 														       
      else {
        const PFCandidateRef& theLeadPFCand = patTau->leadPFChargedHadrCand();
	if(theLeadPFCand.isNonnull()){
	  _tauPdgId->push_back(getMatchedPdgId(theLeadPFCand->pt(), theLeadPFCand->eta(), theLeadPFCand->phi(), theLeadPFCand->charge()).first);
	  _tauMotherPdgId->push_back(getMatchedPdgId(theLeadPFCand->pt(), theLeadPFCand->eta(), theLeadPFCand->phi(), theLeadPFCand->charge()).second);
	}
	else{
	  _tauPdgId->push_back(getMatchedPdgId(patTau->pt(), patTau->eta(), patTau->phi(), patTau->charge()).first);
	  _tauMotherPdgId->push_back(getMatchedPdgId(patTau->pt(), patTau->eta(), patTau->phi(), patTau->charge()).second);
	}
        _tauNProngs->push_back(patTau->signalPFChargedHadrCands().size());							       
        _tauEmFraction->push_back(patTau->emFraction());  									       
        _tauHcalTotOverPLead->push_back(patTau->hcalTotOverPLead());								       
        _tauHcalMaxOverPLead->push_back(patTau->hcalMaxOverPLead());								       
        _tauHcal3x3OverPLead->push_back(patTau->hcal3x3OverPLead());								       
        _tauElectronPreId->push_back(patTau->electronPreIDDecision());								       
        _tauModifiedEOverP->push_back(patTau->ecalStripSumEOverPLead());  							       
        _tauBremsRecoveryEOverPLead->push_back(patTau->bremsRecoveryEOverPLead());						       
        _tauDiscAgainstElec->push_back(patTau->tauID(_RecoTauDiscrAgainstElectron.data()));
	_tauLTSignedIp->push_back(patTau->leadPFChargedHadrCandsignedSipt());
	_tauIsInTheCraks->push_back(isInTheCracks(patTau->eta()));
      } 														       
      _tauIsoTrackPtSum->push_back(CalculateTauTrackIsolation(*patTau, 0.5, 1.0).second);
      _tauIsoTrkPtSumDR1_0MinPt1_0->push_back(CalculateTauTrackIsolation(*patTau, 1.0, 1.0).second);
      _tauIsoTrkPtSumDR1_0MinPt0_5->push_back(CalculateTauTrackIsolation(*patTau, 1.0, 0.5).second);
      _tauIsoTrkPtSumDR0_75MinPt1_0->push_back(CalculateTauTrackIsolation(*patTau, 0.75, 1.0).second);
      _tauIsoTrkPtSumDR0_75MinPt0_5->push_back(CalculateTauTrackIsolation(*patTau, 0.75, 0.5).second);
      _tauIsoGammaEtSum->push_back(CalculateTauEcalIsolation(*patTau, 0.5, 1.5).second);
      _tauIsoGammaEtSumDR1_0MinPt1_5->push_back(CalculateTauEcalIsolation(*patTau, 1.0, 1.5).second);
      _tauIsoGammaEtSumDR1_0MinPt1_0->push_back(CalculateTauEcalIsolation(*patTau, 1.0, 1.0).second);
      _tauIsoGammaEtSumDR0_75MinPt1_5->push_back(CalculateTauEcalIsolation(*patTau, 0.75, 1.5).second);
      _tauIsoGammaEtSumDR0_75MinPt1_0->push_back(CalculateTauEcalIsolation(*patTau, 0.75, 1.0).second);
      _eventIsZee->push_back(int(isZee(*patElectron).first));
      _zeeMass->push_back(isZee(*patElectron).second.first);
      _zeePtAsymm->push_back(isZee(*patElectron).second.second);
      _eMatched->push_back(int(matchToGen(*patElectron).first));
      const GenParticleRef theGenElectronRef = patElectron->genParticleRef();						       
      if(theGenElectronRef.isNonnull()){ 
        _eMotherId->push_back(abs(theGenElectronRef->mother()->pdgId()));
	} 			       
      else{ 
        _eMotherId->push_back(0);
	_ePdgId->push_back(0);
      }											       
      _ePdgId->push_back(getMatchedPdgId(patElectron->pt(), patElectron->eta(), patElectron->phi(), patElectron->charge()).first);
      _eMotherPdgId->push_back(getMatchedPdgId(patElectron->pt(), patElectron->eta(), patElectron->phi(), patElectron->charge()).second);
      _eE->push_back(patElectron->energy());											       
      _ePt->push_back(patElectron->pt()); 											       
      _eEt->push_back(patElectron->et()); 											       
      _eCharge->push_back(patElectron->charge()); 										       
      _eEta->push_back(patElectron->eta());											       
      _ePhi->push_back(patElectron->phi());											       
      _eSigmaEtaEta->push_back(patElectron->scSigmaEtaEta());									       
      _eSigmaIEtaIEta->push_back(patElectron->scSigmaIEtaIEta()); 								       
      _eEOverP->push_back(patElectron->eSuperClusterOverP());									       
      _eHOverEm->push_back(patElectron->hadronicOverEm());									       
      _eDeltaPhiIn->push_back(patElectron->deltaEtaSuperClusterTrackAtVtx());							       
      _eDeltaEtaIn->push_back(patElectron->deltaPhiSuperClusterTrackAtVtx());							       
     
      _eEcalIsoPat->push_back(patElectron->ecalIso());
      _eHcalIsoPat->push_back(patElectron->hcalIso());										       
      _eIsoPat->push_back(patElectron->caloIso());											       
      _eTrkIsoPat->push_back(patElectron->trackIso());										       
     
      _eSCE1x5->push_back(patElectron->scE1x5()); 										       
      _eSCE2x5->push_back(patElectron->scE2x5Max());										       
      _eSCE5x5->push_back(patElectron->scE5x5()); 
      const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
      if(patElectron->track().isNonnull()) _eIp->push_back(patElectron->track()->dxy(thePrimaryEventVertex.position()));
      else _eIp->push_back(-100.);
      _eClass->push_back(patElectron->classification());
      
      _mEt->push_back((*(_patMETs->begin())).pt());
      								       
      // Di-tau variables												       
      math::XYZTLorentzVector The_LorentzVect = patTau->p4() + patElectron->p4();					       
      _eTauMass->push_back(The_LorentzVect.M());  										       
      _eTauCosDPhi->push_back(cos(TMath::Abs(normalizedPhi(patElectron->phi() - patTau->phi()))));				       
      _eTauDelatR->push_back(reco::deltaR(patElectron->p4(), patTau->p4()));
      _UseVectorSumOfVisProductsAndMetMassReco = true;
      _eTauMetMass->push_back(CalculateThe4Momentum((*patTau),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))).second.M());
      _UseVectorSumOfVisProductsAndMetMassReco = false;
      _UseCollinerApproxMassReco = true;
      _eTauCollMass->push_back(CalculateThe4Momentum((*patTau),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))).second.M());

      _eTauPZeta->push_back(CalculatePZeta((*patTau),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))));
      _eTauPZetaVis->push_back(CalculatePZetaVis((*patTau),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))));
    }										  
  }
  //}
  _HMTTree->Fill();

}

// ---------------Fill Histograms
void HiMassTauAnalysis::fillHistograms() {

  // ------Vertices
  if (_FillRecoVertexHists) {
    int nVertices = 0;
    for(reco::VertexCollection::const_iterator primaryVertex = _primaryEventVertexCollection->begin();
        primaryVertex != _primaryEventVertexCollection->end(); ++primaryVertex ) {
      if (!passRecoVertexCuts(*primaryVertex)) continue;
      _hVertexZposition->Fill(primaryVertex->z());
      int pvntrk = 0;
      reco::Vertex::trackRef_iterator pvtrk;
      for(pvtrk=primaryVertex->tracks_begin();pvtrk!=primaryVertex->tracks_end();++pvtrk) {
       if(primaryVertex->trackWeight(*pvtrk) > _RecoVertexTrackWeight) {pvntrk++;}
      }
      _hVertexNTracks->Fill(pvntrk);
      nVertices++;
    }
    _hNVertices->Fill(nVertices);
  }

  // ------Generated Taus
  if (_FillGenTauHists) {
    int nGenTaus = 0;
    for(GenParticleCollection::const_iterator genParticle = _genParticles->begin();genParticle != _genParticles->end();++genParticle) {
      if((abs(genParticle->pdgId()) == 15) && (genParticle->status() != 3)) {
        int neutrinos = 0;
        MChadtau = genParticle->p4();
        for(int ii=0; ii<(int)(genParticle->numberOfDaughters()); ii++) {
          daughterCand = genParticle->daughter(ii);
          if( (abs(daughterCand->pdgId()) == 12) || (abs(daughterCand->pdgId()) == 14) || (abs(daughterCand->pdgId()) == 16) ) {
            neutrinos++;
            MChadtau = MChadtau - daughterCand->p4();
          }
        }
        if(neutrinos == 1) {
          _hGenTauEnergy->Fill(MChadtau.energy());
          _hGenTauPt->Fill(MChadtau.pt());
          _hGenTauEta->Fill(MChadtau.eta());
          _hGenTauPhi->Fill(MChadtau.phi());
          nGenTaus++;
        }
      }
    }
    _hNGenTau->Fill(nGenTaus);
  }

  // ------Reco Tau Histograms
  if (_FillRecoTauHists) {
    int nTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); 
	  patTau != _patTaus->end(); ++patTau ) {
      if (!passRecoTauCuts(*patTau)) continue;
      _hTauJetEnergy->Fill(patTau->energy());
      _hTauJetPt->Fill(patTau->pt());
      _hTauJetEta->Fill(patTau->eta());
      _hTauJetPhi->Fill(patTau->phi());
      if (patTau->isCaloTau()) {
        _hTauJetNumSignalTracks->Fill(patTau->signalTracks().size());
        _hTauJetNumSignalGammas->Fill(0);
        if( (patTau->leadTrack().isNonnull()) ) {
          _hTauJetSeedTrackPt->Fill(patTau->leadTrack()->pt());
          _hTauJetSeedTrackIpSignificance->Fill(patTau->leadTracksignedSipt());
//          _hTauJetSeedTrackNhits->Fill(patTau->leadTrack()->recHitsSize());
//          _hTauJetSeedTrackChi2->Fill(patTau->leadTrack()->chi2());
        }
      } else {
        _hTauJetNumSignalTracks->Fill(patTau->signalPFChargedHadrCands().size());
        _hTauJetNumSignalGammas->Fill(CalculateNumberSignalTauGammas(*patTau));
        if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
          _hTauJetSeedTrackPt->Fill(patTau->leadPFChargedHadrCand()->pt());
          _hTauJetSeedTrackIpSignificance->Fill(patTau->leadPFChargedHadrCandsignedSipt());
          if( (patTau->leadPFChargedHadrCand()->trackRef().isNonnull()) ) {
//            _hTauJetSeedTrackNhits->Fill(patTau->leadPFChargedHadrCand()->trackRef()->recHitsSize());
//            _hTauJetSeedTrackChi2->Fill(patTau->leadPFChargedHadrCand()->trackRef()->chi2());
          }
        }
      }
      _hTauJetCharge->Fill(fabs(patTau->charge()));
      _hTauJetSignalTracksMass->Fill(CalculateTauSignalTracksMass(*patTau).M());
      _hTauJetSignalTracksChargeFraction->Fill(CalculateTauSignalTracksMass(*patTau).pt() / patTau->pt());
      if(matchToGen(*patTau).first) {
        _hTauJetGenTauDeltaPhi->Fill(normalizedPhi(patTau->phi() - matchToGen(*patTau).second.phi()));
        _hTauJetGenTauDeltaEta->Fill(patTau->eta() - matchToGen(*patTau).second.eta());
        _hTauJetGenTauDeltaPt->Fill((patTau->pt() - matchToGen(*patTau).second.pt()) / matchToGen(*patTau).second.pt());
      }
      if (_UseRecoTauDiscrByIsolationFlag) {
        if (patTau->isCaloTau()) {
          _hTauJetNumIsoTracks->Fill(patTau->isolationTracks().size());
          _hTauJetNumIsoGammas->Fill(0);
          _hTauJetNumIsoCands->Fill(patTau->isolationTracks().size() + 0);
          _hTauJetSumPtIsoTracks->Fill(patTau->isolationTracksPtSum());
          _hTauJetSumPtIsoGammas->Fill(patTau->isolationECALhitsEtSum());
          _hTauJetSumPtIso->Fill(patTau->isolationTracksPtSum() + patTau->isolationECALhitsEtSum());
        } else {
          _hTauJetNumIsoTracks->Fill(patTau->isolationPFChargedHadrCands().size());
          _hTauJetNumIsoGammas->Fill(patTau->isolationPFGammaCands().size());
          _hTauJetNumIsoCands->Fill(patTau->isolationPFChargedHadrCands().size() + patTau->isolationPFGammaCands().size());
          _hTauJetSumPtIsoTracks->Fill(patTau->isolationPFChargedHadrCandsPtSum());
          _hTauJetSumPtIsoGammas->Fill(patTau->isolationPFGammaCandsEtSum());
          _hTauJetSumPtIso->Fill(patTau->isolationPFChargedHadrCandsPtSum() + patTau->isolationPFGammaCandsEtSum());
        }
      } else {
        _hTauJetNumIsoTracks->Fill(CalculateTauTrackIsolation(*patTau).first);
        _hTauJetNumIsoGammas->Fill(CalculateTauEcalIsolation(*patTau).first);
        _hTauJetNumIsoCands->Fill(CalculateTauTrackIsolation(*patTau).first + CalculateTauEcalIsolation(*patTau).first);
        _hTauJetSumPtIsoTracks->Fill(CalculateTauTrackIsolation(*patTau).second);
        _hTauJetSumPtIsoGammas->Fill(CalculateTauEcalIsolation(*patTau).second);
        _hTauJetSumPtIso->Fill(CalculateTauTrackIsolation(*patTau).second + CalculateTauEcalIsolation(*patTau).second);
      }
      nTaus++;
    }
    _hNTau->Fill(nTaus);
  }

  // ------Reco Muon Histograms
  if (_FillRecoMuonHists) {
    int nMuons = 0;
    int theNumberOfMuons = 0;
    for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); 
	  patMuon != _patMuons->end(); ++patMuon ) {
      theNumberOfMuons++;
      if (!passRecoMuonCuts((*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1))) continue;
      if(_SmearTheMuon) {
        _hMuonEnergy->Fill(smearedMuonMomentumVector.at(theNumberOfMuons-1).energy());
        _hMuonPt->Fill(smearedMuonMomentumVector.at(theNumberOfMuons-1).pt());
        _hMuonEta->Fill(smearedMuonMomentumVector.at(theNumberOfMuons-1).eta());
        _hMuonPhi->Fill(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi());
        if(matchToGen(*patMuon).first) {
          _hMuonGenMuonDeltaPhi->Fill(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi() - matchToGen(*patMuon).second.phi()));
          _hMuonGenMuonDeltaEta->Fill(smearedMuonMomentumVector.at(theNumberOfMuons-1).eta() - matchToGen(*patMuon).second.eta());
          _hMuonGenMuonDeltaPt->Fill((smearedMuonMomentumVector.at(theNumberOfMuons-1).pt() - matchToGen(*patMuon).second.pt()) / matchToGen(*patMuon).second.pt());
        }
      } else {
        _hMuonEnergy->Fill(patMuon->energy());
        _hMuonPt->Fill(patMuon->pt());
        _hMuonEta->Fill(patMuon->eta());
        _hMuonPhi->Fill(patMuon->phi());
        if(matchToGen(*patMuon).first) {
          _hMuonGenMuonDeltaPhi->Fill(normalizedPhi(patMuon->phi() - matchToGen(*patMuon).second.phi()));
          _hMuonGenMuonDeltaEta->Fill(patMuon->eta() - matchToGen(*patMuon).second.eta());
          _hMuonGenMuonDeltaPt->Fill((patMuon->pt() - matchToGen(*patMuon).second.pt()) / matchToGen(*patMuon).second.pt());
        }
      }
      _hMuonTrackIso->Fill(patMuon->trackIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonTrackIsoTrkThreshold).first);
      _hMuonEcalIso->Fill(patMuon->ecalIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonEcalIsoRecHitThreshold).first);
      _hMuonIso->Fill(patMuon->trackIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonTrackIsoTrkThreshold).first +
                      patMuon->ecalIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonEcalIsoRecHitThreshold).first);
/*
      _hMuonTrackIso->Fill(patMuon->trackIso());
      _hMuonEcalIso->Fill(patMuon->ecalIso());
      _hMuonIso->Fill(patMuon->trackIso() + patMuon->ecalIso());
*/
      const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
      if ( patMuon->track().isNonnull() ) {
        _hMuonIp->Fill( patMuon->track()->dxy(thePrimaryEventVertex.position()) );
        if(fabs(patMuon->track()->dxyError()) != 0) {
          _hMuonIpSignificance->Fill( fabs(patMuon->track()->dxy(thePrimaryEventVertex.position()) / patMuon->track()->dxyError()) );
        } else {_hMuonIpSignificance->Fill(-1);}
      }
      nMuons++;
    }
    _hNMuon->Fill(nMuons);
  }

  //-----Fill Reco Electron Histograms
  if (_FillRecoElectronHists) {
    int nElectrons = 0;
    int theNumberOfElectrons = 0;
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();
          patElectron != _patElectrons->end(); ++patElectron ) {
      theNumberOfElectrons++;
      if (!passRecoElectronCuts((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1))) continue;
      if(_SmearTheElectron) {
        _hElectronEnergy->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).energy());
        _hElectronPt->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt());
        _hElectronEta->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).eta());
        _hElectronPhi->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi());
        if(matchToGen(*patElectron).first) {
          _hElectronGenElectronDeltaPhi->Fill(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - matchToGen(*patElectron).second.phi()));
          _hElectronGenElectronDeltaEta->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).eta() - matchToGen(*patElectron).second.eta());
          _hElectronGenElectronDeltaPt->Fill((smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt() - matchToGen(*patElectron).second.pt()) / matchToGen(*patElectron).second.pt());
        }
      } else {
        _hElectronEnergy->Fill(patElectron->energy());
        _hElectronPt->Fill(patElectron->pt());
        _hElectronEta->Fill(patElectron->eta());
        _hElectronPhi->Fill(patElectron->phi());
        if(matchToGen(*patElectron).first) {
          _hElectronGenElectronDeltaPhi->Fill(normalizedPhi(patElectron->phi() - matchToGen(*patElectron).second.phi()));
          _hElectronGenElectronDeltaEta->Fill(patElectron->eta() - matchToGen(*patElectron).second.eta());
          _hElectronGenElectronDeltaPt->Fill((patElectron->pt() - matchToGen(*patElectron).second.pt()) / matchToGen(*patElectron).second.pt());
        }
      }
      _hElectronTrackIso->Fill(patElectron->trackIso());
      _hElectronEcalIso->Fill(patElectron->ecalIso());
/*
      _hElectronTrackIso->Fill(patElectron->trackIsoDeposit()->depositAndCountWithin(_RecoElectronTrackIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoElectronTrackIsoTrkThreshold).first);
      _hElectronEcalIso->Fill(patElectron->ecalIsoDeposit()->depositAndCountWithin(_RecoElectronEcalIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoElectronEcalIsoRecHitThreshold).first);
*/
      const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
      if ( patElectron->track().isNonnull() ) { _hElectronIp->Fill( patElectron->track()->dxy(thePrimaryEventVertex.position()) ); }
//      if ( patElectron->gsfTrack().isNonnull() ) { _hElectronIp->Fill( patElectron->gsfTrack()->dxy(thePrimaryEventVertex.position()) ); }
      _hElectronEoverP->Fill(patElectron->eSuperClusterOverP());
      _hElectronHoverEm->Fill(patElectron->hadronicOverEm());
      _hElectronClassification->Fill(patElectron->classification());
      nElectrons++;
    }
    _hNElectron->Fill(nElectrons);
  }

  // ------Reco Jet Histograms
  if (_FillRecoJetHists) {
    int nJets = 0;
    for ( pat::JetCollection::const_iterator patJet = _patJets->begin(); 
	  patJet != _patJets->end(); ++patJet ) {
      if (!passRecoJetCuts(*patJet)) continue;
      if(_UseCorrectedJet) {
        _hJetEnergy->Fill(patJet->energy());
        _hJetPt->Fill(patJet->pt());
        _hJetEta->Fill(patJet->eta());
        _hJetPhi->Fill(patJet->phi());
        nJets++;
      } else {
        _hJetEnergy->Fill(patJet->correctedJet("raw","").energy());
        _hJetPt->Fill(patJet->correctedJet("raw","").pt());
        _hJetEta->Fill(patJet->correctedJet("raw","").eta());
        _hJetPhi->Fill(patJet->correctedJet("raw","").phi());
        nJets++;
      }
      _hBJetDiscrByTrackCounting->Fill(patJet->bDiscriminator("trackCountingHighEffBJetTags"));
      _hBJetDiscrBySimpleSecondaryV->Fill(patJet->bDiscriminator("simpleSecondaryVertexBJetTags"));
      _hBJetDiscrByCombinedSecondaryV->Fill(patJet->bDiscriminator("combinedSecondaryVertexBJetTags"));
    }
    _hNJet->Fill(nJets);
  }

  // ------Topology Histograms
  if (_FillTopologyHists) {
    _hMet->Fill((*(_patMETs->begin())).pt());
    if( ((_AnalyzeMuonForLeg1) && (_AnalyzeTauForLeg2)) || ((_AnalyzeMuonForLeg2) && (_AnalyzeTauForLeg1)) ) {
      int theNumberOfMuons = 0;
      for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
        theNumberOfMuons++;
        for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
          if ((passRecoTauCuts(*patTau)) && 
              (passRecoMuonCuts((*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1))) &&
              (passTopologyCuts((*patTau),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))))) {
            if(_SmearTheMuon) {
              _hMuonPtVsTauPt->Fill(smearedMuonMomentumVector.at(theNumberOfMuons-1).pt(),patTau->pt());
              _hMuonTauDeltaR->Fill(reco::deltaR(patTau->p4(), smearedMuonMomentumVector.at(theNumberOfMuons-1)));
              _hMuonTauDeltaPtDivSumPt->Fill((patTau->pt() - smearedMuonMomentumVector.at(theNumberOfMuons-1).pt()) / (patTau->pt() + smearedMuonMomentumVector.at(theNumberOfMuons-1).pt()));
              _hMuonTauCosDphi->Fill(cos(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi() - patTau->phi()))));
              _hMuonMetDeltaPhi->Fill(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi() - (*(_patMETs->begin())).phi())));
              _hMuonMetDeltaPhiVsMuonTauCosDphi->Fill(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi() - patTau->phi()))));
            } else {
              _hMuonPtVsTauPt->Fill(patMuon->pt(),patTau->pt());
              _hMuonTauDeltaR->Fill(reco::deltaR(patTau->p4(), patMuon->p4()));
              _hMuonTauDeltaPtDivSumPt->Fill((patTau->pt() - patMuon->pt()) / (patTau->pt() + patMuon->pt()));
              _hMuonTauCosDphi->Fill(cos(TMath::Abs(normalizedPhi(patMuon->phi() - patTau->phi()))));
              _hMuonMetDeltaPhi->Fill(TMath::Abs(normalizedPhi(patMuon->phi() - (*(_patMETs->begin())).phi())));
              _hMuonMetDeltaPhiVsMuonTauCosDphi->Fill(TMath::Abs(normalizedPhi(patMuon->phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(patMuon->phi() - patTau->phi()))));
            }
            if(CalculateThe4Momentum((*patTau),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))).first) {_hReconstructableMass->Fill(CalculateThe4Momentum((*patTau),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))).second.M());}
            else {_hNotReconstructableMass->Fill(CalculateThe4Momentum((*patTau),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))).second.M());}
            _hMuonMetMt->Fill(CalculateLeptonMetMt((*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))));
            _hTauMetMt->Fill(CalculateLeptonMetMt((*patTau),(*(_patMETs->begin()))));
            if(_UseTauSeedTrackForDiTauDiscrByOSLS) {
              if (patTau->isCaloTau()) {
                if( (patTau->leadTrack().isNonnull()) ) {
                  _hMuonTauOSLS->Fill(patMuon->charge() * patTau->leadTrack()->charge());
                }
              } else {
                if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
                  _hMuonTauOSLS->Fill(patMuon->charge() * patTau->leadPFChargedHadrCand()->charge());
                }
              }
            } else {_hMuonTauOSLS->Fill(patMuon->charge() * patTau->charge());}
            _hTauMetDeltaPhi->Fill(TMath::Abs(normalizedPhi(patTau->phi() - (*(_patMETs->begin())).phi())));
            _hPZeta->Fill(CalculatePZeta((*patTau),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))));
            _hPZetaVis->Fill(CalculatePZetaVis((*patTau),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))));
            _hZeta2D->Fill(CalculatePZetaVis((*patTau),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))),CalculatePZeta((*patTau),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))));
            _hZeta1D->Fill((_PZetaCutCoefficient * CalculatePZeta((*patTau),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin())))) + 
                           (_PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin())))));
          }
        }
      }
    }
    if( ((_AnalyzeElectronForLeg1) && (_AnalyzeTauForLeg2)) || ((_AnalyzeElectronForLeg2) && (_AnalyzeTauForLeg1)) ) {
      int theNumberOfElectrons = 0;
      for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
        theNumberOfElectrons++;
        for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
          if ((passRecoTauCuts(*patTau)) && 
              (passRecoElectronCuts((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1))) &&
              (passTopologyCuts((*patTau),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))))) {
            if(_SmearTheElectron) {
              _hElectronPtVsTauPt->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt(),patTau->pt());
              _hElectronTauDeltaR->Fill(reco::deltaR(patTau->p4(), smearedElectronMomentumVector.at(theNumberOfElectrons-1)));
              _hElectronTauDeltaPtDivSumPt->Fill((patTau->pt() - smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt()) / (patTau->pt() + smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt()));
              _hElectronTauCosDphi->Fill(cos(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - patTau->phi()))));
              _hElectronMetDeltaPhi->Fill(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - (*(_patMETs->begin())).phi())));
              _hElectronMetDeltaPhiVsElectronTauCosDphi->Fill(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - patTau->phi()))));
            } else {
              _hElectronPtVsTauPt->Fill(patElectron->pt(),patTau->pt());
              _hElectronTauDeltaR->Fill(reco::deltaR(patTau->p4(), patElectron->p4()));
              _hElectronTauDeltaPtDivSumPt->Fill((patTau->pt() - patElectron->pt()) / (patTau->pt() + patElectron->pt()));
              _hElectronTauCosDphi->Fill(cos(TMath::Abs(normalizedPhi(patElectron->phi() - patTau->phi()))));
              _hElectronMetDeltaPhi->Fill(TMath::Abs(normalizedPhi(patElectron->phi() - (*(_patMETs->begin())).phi())));
              _hElectronMetDeltaPhiVsElectronTauCosDphi->Fill(TMath::Abs(normalizedPhi(patElectron->phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(patElectron->phi() - patTau->phi()))));
            }
            if(CalculateThe4Momentum((*patTau),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))).first) {_hReconstructableMass->Fill(CalculateThe4Momentum((*patTau),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))).second.M());}
            else {_hNotReconstructableMass->Fill(CalculateThe4Momentum((*patTau),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))).second.M());}
            _hElectronMetMt->Fill(CalculateLeptonMetMt((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))));
            _hTauMetMt->Fill(CalculateLeptonMetMt((*patTau),(*(_patMETs->begin()))));
            if(_UseTauSeedTrackForDiTauDiscrByOSLS) {
              if (patTau->isCaloTau()) {
                if( (patTau->leadTrack().isNonnull()) ) {
                  _hElectronTauOSLS->Fill(patElectron->charge() * patTau->leadTrack()->charge());
                }
              } else {
                if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
                  _hElectronTauOSLS->Fill(patElectron->charge() * patTau->leadPFChargedHadrCand()->charge());
                }
              }
            } else {_hElectronTauOSLS->Fill(patElectron->charge() * patTau->charge());}
            _hTauMetDeltaPhi->Fill(TMath::Abs(normalizedPhi(patTau->phi() - (*(_patMETs->begin())).phi())));
            _hPZeta->Fill(CalculatePZeta((*patTau),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))));
            _hPZetaVis->Fill(CalculatePZetaVis((*patTau),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))));
            _hZeta2D->Fill(CalculatePZetaVis((*patTau),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))),CalculatePZeta((*patTau),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))));
            _hZeta1D->Fill((_PZetaCutCoefficient * CalculatePZeta((*patTau),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin())))) + 
                           (_PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin())))));
          }
        }
      }
    }
    if( ((_AnalyzeMuonForLeg1) && (_AnalyzeElectronForLeg2)) || ((_AnalyzeMuonForLeg2) && (_AnalyzeElectronForLeg1)) ) {
      int theNumberOfMuons = 0;
      for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
        theNumberOfMuons++;
        int theNumberOfElectrons = 0;
        for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
          theNumberOfElectrons++;
          if ((passRecoElectronCuts((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1))) && 
              (passRecoMuonCuts((*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1))) && 
              (passTopologyCuts((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))))) {
            if(_SmearTheElectron) {
              if(_SmearTheMuon) {
                _hElectronPtVsMuonPt->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt(),smearedMuonMomentumVector.at(theNumberOfMuons-1).pt());
                _hElectronMuonDeltaR->Fill(reco::deltaR(smearedMuonMomentumVector.at(theNumberOfMuons-1), smearedElectronMomentumVector.at(theNumberOfElectrons-1)));
                _hElectronMuonDeltaPtDivSumPt->Fill((smearedMuonMomentumVector.at(theNumberOfMuons-1).pt() - smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt()) / (smearedMuonMomentumVector.at(theNumberOfMuons-1).pt() + smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt()));
                _hElectronMuonCosDphi->Fill(cos(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - smearedMuonMomentumVector.at(theNumberOfMuons-1).phi()))));
                _hElectronMetDeltaPhi->Fill(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - (*(_patMETs->begin())).phi())));
                _hMuonMetDeltaPhi->Fill(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi() - (*(_patMETs->begin())).phi())));
                _hElectronMetDeltaPhiVsElectronMuonCosDphi->Fill(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - smearedMuonMomentumVector.at(theNumberOfMuons-1).phi()))));
              } else {
                _hElectronPtVsMuonPt->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt(),patMuon->pt());
                _hElectronMuonDeltaR->Fill(reco::deltaR(patMuon->p4(), smearedElectronMomentumVector.at(theNumberOfElectrons-1)));
                _hElectronMuonDeltaPtDivSumPt->Fill((patMuon->pt() - smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt()) / (patMuon->pt() + smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt()));
                _hElectronMuonCosDphi->Fill(cos(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - patMuon->phi()))));
                _hElectronMetDeltaPhi->Fill(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - (*(_patMETs->begin())).phi())));
                _hMuonMetDeltaPhi->Fill(TMath::Abs(normalizedPhi(patMuon->phi() - (*(_patMETs->begin())).phi())));
                _hElectronMetDeltaPhiVsElectronMuonCosDphi->Fill(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - patMuon->phi()))));
              }
            } else {
              if(_SmearTheMuon) {
                _hElectronPtVsMuonPt->Fill(patElectron->pt(),smearedMuonMomentumVector.at(theNumberOfMuons-1).pt());
                _hElectronMuonDeltaR->Fill(reco::deltaR(smearedMuonMomentumVector.at(theNumberOfMuons-1), patElectron->p4()));
                _hElectronMuonDeltaPtDivSumPt->Fill((smearedMuonMomentumVector.at(theNumberOfMuons-1).pt() - patElectron->pt()) / (smearedMuonMomentumVector.at(theNumberOfMuons-1).pt() + patElectron->pt()));
                _hElectronMuonCosDphi->Fill(cos(TMath::Abs(normalizedPhi(patElectron->phi() - smearedMuonMomentumVector.at(theNumberOfMuons-1).phi()))));
                _hElectronMetDeltaPhi->Fill(TMath::Abs(normalizedPhi(patElectron->phi() - (*(_patMETs->begin())).phi())));
                _hMuonMetDeltaPhi->Fill(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi() - (*(_patMETs->begin())).phi())));
                _hElectronMetDeltaPhiVsElectronMuonCosDphi->Fill(TMath::Abs(normalizedPhi(patElectron->phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(patElectron->phi() - smearedMuonMomentumVector.at(theNumberOfMuons-1).phi()))));
              } else {
                _hElectronPtVsMuonPt->Fill(patElectron->pt(),patMuon->pt());
                _hElectronMuonDeltaR->Fill(reco::deltaR(patMuon->p4(), patElectron->p4()));
                _hElectronMuonDeltaPtDivSumPt->Fill((patMuon->pt() - patElectron->pt()) / (patMuon->pt() + patElectron->pt()));
                _hElectronMuonCosDphi->Fill(cos(TMath::Abs(normalizedPhi(patElectron->phi() - patMuon->phi()))));
                _hElectronMetDeltaPhi->Fill(TMath::Abs(normalizedPhi(patElectron->phi() - (*(_patMETs->begin())).phi())));
                _hMuonMetDeltaPhi->Fill(TMath::Abs(normalizedPhi(patMuon->phi() - (*(_patMETs->begin())).phi())));
                _hElectronMetDeltaPhiVsElectronMuonCosDphi->Fill(TMath::Abs(normalizedPhi(patElectron->phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(patElectron->phi() - patMuon->phi()))));
              }
            }
            if(CalculateThe4Momentum((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))).first) {_hReconstructableMass->Fill(CalculateThe4Momentum((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))).second.M());}
            else {_hNotReconstructableMass->Fill(CalculateThe4Momentum((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))).second.M());}
            _hElectronMetMt->Fill(CalculateLeptonMetMt((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))));
            _hMuonMetMt->Fill(CalculateLeptonMetMt((*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))));
            _hElectronMuonOSLS->Fill(patElectron->charge() * patMuon->charge());
            _hPZeta->Fill(CalculatePZeta((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))));
            _hPZetaVis->Fill(CalculatePZetaVis((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))));
            _hZeta2D->Fill(CalculatePZetaVis((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))),CalculatePZeta((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))));
            _hZeta1D->Fill((_PZetaCutCoefficient * CalculatePZeta((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin())))) + 
                           (_PZetaVisCutCoefficient * CalculatePZetaVis((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin())))));
          }
        }
      }
    }
    if( ((_AnalyzeTauForLeg1) && (_AnalyzeTauForLeg2)) ) {
      for ( pat::TauCollection::const_iterator patTau1 = _patTaus->begin();patTau1 != _patTaus->end(); ++patTau1 ) {
        for ( pat::TauCollection::const_iterator patTau2 = _patTaus->begin();patTau2 != _patTaus->end(); ++patTau2 ) {
          if ((passRecoTauCuts(*patTau1)) && (passRecoTauCuts(*patTau2)) && (passTopologyCuts((*patTau1),(*patTau2),(*(_patMETs->begin()))))) {
            _hTau1PtVsTau2Pt->Fill(patTau1->pt(),patTau2->pt());
            if(CalculateThe4Momentum((*patTau1),(*patTau2),(*(_patMETs->begin()))).first) {_hReconstructableMass->Fill(CalculateThe4Momentum((*patTau1),(*patTau2),(*(_patMETs->begin()))).second.M());}
            else {_hNotReconstructableMass->Fill(CalculateThe4Momentum((*patTau1),(*patTau2),(*(_patMETs->begin()))).second.M());}
            _hTau1Tau2DeltaR->Fill(reco::deltaR(patTau1->p4(), patTau2->p4()));
            _hTau1Tau2DeltaPtDivSumPt->Fill((patTau1->pt() - patTau2->pt()) / (patTau1->pt() + patTau2->pt()));
            _hTau1MetMt->Fill(CalculateLeptonMetMt((*patTau1),(*(_patMETs->begin()))));
            _hTau2MetMt->Fill(CalculateLeptonMetMt((*patTau2),(*(_patMETs->begin()))));
            if(_UseTauSeedTrackForDiTauDiscrByOSLS) {
              if (patTau1->isCaloTau()) {
                if( (patTau1->leadTrack().isNonnull()) && (patTau2->leadTrack().isNonnull()) ) {
                  _hTau1Tau2OSLS->Fill(patTau1->leadTrack()->charge() * patTau2->leadTrack()->charge());
                }
              } else {
                if( (patTau1->leadPFChargedHadrCand().isNonnull()) && (patTau2->leadPFChargedHadrCand().isNonnull()) ) {
                  _hTau1Tau2OSLS->Fill(patTau1->leadPFChargedHadrCand()->charge() * patTau2->leadPFChargedHadrCand()->charge());
                }
              }
            } else {_hTau1Tau2OSLS->Fill(patTau1->charge() * patTau2->charge());}
            _hTau1Tau2CosDphi->Fill(cos(TMath::Abs(normalizedPhi(patTau1->phi() - patTau2->phi()))));
            _hTau1MetDeltaPhi->Fill(TMath::Abs(normalizedPhi(patTau1->phi() - (*(_patMETs->begin())).phi())));
            _hTau2MetDeltaPhi->Fill(TMath::Abs(normalizedPhi(patTau2->phi() - (*(_patMETs->begin())).phi())));
            _hTau1MetDeltaPhiVsTau1Tau2CosDphi->Fill(TMath::Abs(normalizedPhi(patTau1->phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(patTau1->phi() - patTau2->phi()))));
            _hPZeta->Fill(CalculatePZeta((*patTau1),(*patTau2),(*(_patMETs->begin()))));
            _hPZetaVis->Fill(CalculatePZetaVis((*patTau1),(*patTau2),(*(_patMETs->begin()))));
            _hZeta2D->Fill(CalculatePZetaVis((*patTau1),(*patTau2),(*(_patMETs->begin()))),CalculatePZeta((*patTau1),(*patTau2),(*(_patMETs->begin()))));
            _hZeta1D->Fill((_PZetaCutCoefficient * CalculatePZeta((*patTau1),(*patTau2),(*(_patMETs->begin())))) + 
                           (_PZetaVisCutCoefficient * CalculatePZetaVis((*patTau1),(*patTau2),(*(_patMETs->begin())))));
          }
        }
      }
    }
    if( ((_AnalyzeMuonForLeg1) && (_AnalyzeMuonForLeg2)) ) {
      int theNumberOfMuons1 = 0;
      for ( pat::MuonCollection::const_iterator patMuon1 = _patMuons->begin();patMuon1 != _patMuons->end(); ++patMuon1 ) {
        theNumberOfMuons1++;
        int theNumberOfMuons2 = 0;
        for ( pat::MuonCollection::const_iterator patMuon2 = _patMuons->begin();patMuon2 != _patMuons->end(); ++patMuon2 ) {
          theNumberOfMuons2++;
          if ((passRecoMuonCuts((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1))) && 
              (passRecoMuonCuts((*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1))) && 
              (passTopologyCuts((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))))) {
            if(_SmearTheMuon) {
              _hMuon1PtVsMuon2Pt->Fill(smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).pt(),smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).pt());
              _hMuon1Muon2DeltaR->Fill(reco::deltaR(smearedMuonMomentumVector.at(theNumberOfMuons1 - 1), smearedMuonMomentumVector.at(theNumberOfMuons2 - 1)));
              _hMuon1Muon2DeltaPtDivSumPt->Fill((smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).pt() - smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).pt()) / (smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).pt() + smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).pt()));
              _hMuon1Muon2CosDphi->Fill(cos(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).phi() - smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).phi()))));
              _hMuon1MetDeltaPhi->Fill(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).phi() - (*(_patMETs->begin())).phi())));
              _hMuon2MetDeltaPhi->Fill(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).phi() - (*(_patMETs->begin())).phi())));
              _hMuon1MetDeltaPhiVsMuon1Muon2CosDphi->Fill(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).phi() - smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).phi()))));
            } else {
              _hMuon1PtVsMuon2Pt->Fill(patMuon1->pt(),patMuon2->pt());
              _hMuon1Muon2DeltaR->Fill(reco::deltaR(patMuon1->p4(), patMuon2->p4()));
              _hMuon1Muon2DeltaPtDivSumPt->Fill((patMuon1->pt() - patMuon2->pt()) / (patMuon1->pt() + patMuon2->pt()));
              _hMuon1Muon2CosDphi->Fill(cos(TMath::Abs(normalizedPhi(patMuon1->phi() - patMuon2->phi()))));
              _hMuon1MetDeltaPhi->Fill(TMath::Abs(normalizedPhi(patMuon1->phi() - (*(_patMETs->begin())).phi())));
              _hMuon2MetDeltaPhi->Fill(TMath::Abs(normalizedPhi(patMuon2->phi() - (*(_patMETs->begin())).phi())));
              _hMuon1MetDeltaPhiVsMuon1Muon2CosDphi->Fill(TMath::Abs(normalizedPhi(patMuon1->phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(patMuon1->phi() - patMuon2->phi()))));
            }
            if(CalculateThe4Momentum((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))).first) {_hReconstructableMass->Fill(CalculateThe4Momentum((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))).second.M());}
            else {_hNotReconstructableMass->Fill(CalculateThe4Momentum((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))).second.M());}
            _hMuon1MetMt->Fill(CalculateLeptonMetMt((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*(_patMETs->begin()))));
            _hMuon2MetMt->Fill(CalculateLeptonMetMt((*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))));
            _hMuon1Muon2OSLS->Fill(patMuon1->charge() * patMuon2->charge());
            _hPZeta->Fill(CalculatePZeta((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))));
            _hPZetaVis->Fill(CalculatePZetaVis((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))));
            _hZeta2D->Fill(CalculatePZetaVis((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))),CalculatePZeta((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))));
            _hZeta1D->Fill((_PZetaCutCoefficient * CalculatePZeta((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin())))) + 
                           (_PZetaVisCutCoefficient * CalculatePZetaVis((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin())))));
          }
        }
      }
    }
    if( ((_AnalyzeElectronForLeg1) && (_AnalyzeElectronForLeg2)) ) {
      int theNumberOfElectrons1 = 0;
      for ( pat::ElectronCollection::const_iterator patElectron1 = _patElectrons->begin();patElectron1 != _patElectrons->end(); ++patElectron1 ) {
        theNumberOfElectrons1++;
        int theNumberOfElectrons2 = 0;
        for ( pat::ElectronCollection::const_iterator patElectron2 = _patElectrons->begin();patElectron2 != _patElectrons->end(); ++patElectron2 ) {
          theNumberOfElectrons2++;
          if ((passRecoElectronCuts((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1))) && 
              (passRecoElectronCuts((*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1))) && 
              (passTopologyCuts((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))))) {
            if(_SmearTheElectron) {
              _hElectron1PtVsElectron2Pt->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).pt(),smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).pt());
              _hElectron1Electron2DeltaR->Fill(reco::deltaR(smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1), smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1)));
              _hElectron1Electron2DeltaPtDivSumPt->Fill((smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).pt() - smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).pt()) / (smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).pt() + smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).pt()));
              _hElectron1Electron2CosDphi->Fill(cos(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).phi() - smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).phi()))));
              _hElectron1MetDeltaPhi->Fill(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).phi() - (*(_patMETs->begin())).phi())));
              _hElectron2MetDeltaPhi->Fill(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).phi() - (*(_patMETs->begin())).phi())));
              _hElectron1MetDeltaPhiVsElectron1Electron2CosDphi->Fill(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).phi() - smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).phi()))));
            } else {
              _hElectron1PtVsElectron2Pt->Fill(patElectron1->pt(),patElectron2->pt());
              _hElectron1Electron2DeltaR->Fill(reco::deltaR(patElectron1->p4(), patElectron2->p4()));
              _hElectron1Electron2DeltaPtDivSumPt->Fill((patElectron1->pt() - patElectron2->pt()) / (patElectron1->pt() + patElectron2->pt()));
              _hElectron1Electron2CosDphi->Fill(cos(TMath::Abs(normalizedPhi(patElectron1->phi() - patElectron2->phi()))));
              _hElectron1MetDeltaPhi->Fill(TMath::Abs(normalizedPhi(patElectron1->phi() - (*(_patMETs->begin())).phi())));
              _hElectron2MetDeltaPhi->Fill(TMath::Abs(normalizedPhi(patElectron2->phi() - (*(_patMETs->begin())).phi())));
              _hElectron1MetDeltaPhiVsElectron1Electron2CosDphi->Fill(TMath::Abs(normalizedPhi(patElectron1->phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(patElectron1->phi() - patElectron2->phi()))));
            }
            if(CalculateThe4Momentum((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))).first) {_hReconstructableMass->Fill(CalculateThe4Momentum((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))).second.M());}
            else {_hNotReconstructableMass->Fill(CalculateThe4Momentum((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))).second.M());}
            _hElectron1MetMt->Fill(CalculateLeptonMetMt((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*(_patMETs->begin()))));
            _hElectron2MetMt->Fill(CalculateLeptonMetMt((*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))));
            _hElectron1Electron2OSLS->Fill(patElectron1->charge() * patElectron2->charge());
            _hPZeta->Fill(CalculatePZeta((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))));
            _hPZetaVis->Fill(CalculatePZetaVis((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))));
            _hZeta2D->Fill(CalculatePZetaVis((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))),CalculatePZeta((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))));
            _hZeta1D->Fill((_PZetaCutCoefficient * CalculatePZeta((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin())))) + 
                           (_PZetaVisCutCoefficient * CalculatePZetaVis((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin())))));
          }
        }
      }
    }
  }
}

pair<const pat::Tau*, const pat::Electron*> HiMassTauAnalysis::getPATComponents(const reco::CompositeCandidate& theCandiate){
  const Candidate* theTauCand = theCandiate.daughter("tau");
  const Candidate* theElectronCand = theCandiate.daughter("electron");
  const pat::Tau* theTau = dynamic_cast<const pat::Tau*>(theTauCand);
  const pat::Electron* theElectron = dynamic_cast<const pat::Electron*>(theElectronCand);
  if(theTau == 0 && theTauCand->hasMasterClonePtr())  theTau = dynamic_cast<const pat::Tau*>(&*(theTauCand->masterClonePtr()));
  if(theElectron == 0 && theElectronCand->hasMasterClonePtr()) theElectron = dynamic_cast<const pat::Electron*>(&*(theElectronCand->masterClonePtr()));
  pair<const pat::Tau*, const pat::Electron*> thePatComponents = make_pair<const pat::Tau*, const pat::Electron*>(theTau, theElectron);
  return thePatComponents;
}

// ---------------
void HiMassTauAnalysis::getCollections(const Event& iEvent, const EventSetup& iSetup) {
  iEvent.getByLabel(_RecoTauSource, _patTaus);
  iEvent.getByLabel(_RecoMuonSource, _patMuons);
  iEvent.getByLabel(_GenParticleSource, _genParticles);
  iEvent.getByLabel(_RecoElectronSource, _patElectrons);
  iEvent.getByLabel(_RecoJetSource, _patJets);
  iEvent.getByLabel(_RecoMetSource, _patMETs);
  iEvent.getByLabel(_RecoDiTauSource, _patDiTaus);
  iEvent.getByLabel(_RecoVertexSource, _primaryEventVertexCollection);
  iEvent.getByLabel(_RecoTriggerSource, _triggerResults);
/*
  iEvent.getByLabel("muons", _recoMuonsForMetCorrections);
  iEvent.getByLabel("muonMETValueMapProducer", "muCorrData", vm_muCorrData_h);
  iEvent.getByLabel("selectedLayer1Muons", _patMuonsForMetCorrections);
*/
}

pair<bool, pair<float, float> > HiMassTauAnalysis::isZee(const pat::Electron& theObject) {
  pair<bool, pair<float, float> > theOutPutPair;
  bool eventIsZee = false;
  bool massWindow = false;
  bool ptAsymmWindow = false;
  const float zeeMass = 90.1876;
  const float zeeWidht = 2.4952;
  float zeePtAsymmetry = -10.;
  pair<float, float> theMassPtAsymmPair;
  // if mass is within 3 sigmas of z or pt asymmetry is small set to true.				     
  for(pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron){
    if(reco::deltaR(theObject.p4(), patElectron->p4()) < _DiTauDeltaRCut) continue;
    if(theObject.p4() == patElectron->p4())continue;
    math::XYZTLorentzVector The_LorentzVect = theObject.p4() + patElectron->p4();
    zeePtAsymmetry = (theObject.pt() - patElectron->pt())/(theObject.pt() + patElectron->pt());
    theMassPtAsymmPair = make_pair<float, float>(The_LorentzVect.M(), zeePtAsymmetry);
    
    if(The_LorentzVect.M() < (zeeMass + 3.*zeeWidht) && The_LorentzVect.M() > (zeeMass - 3.*zeeWidht))massWindow = true;
    if(fabs(zeePtAsymmetry) < 0.20) ptAsymmWindow = true;
    if(massWindow || ptAsymmWindow){
      eventIsZee = true;
      break;
    }
  }
  theOutPutPair = make_pair<bool, pair<float, float> >(eventIsZee, theMassPtAsymmPair);
  return theOutPutPair;
}

//-----Matching to generator level objects
// get the pdg id of the gen particle closest in DeltaR to the reco/pat object
pair<unsigned int, unsigned int> HiMassTauAnalysis::getMatchedPdgId(float pt, float eta, float phi, int charge){
//pair<unsigned int, unsigned int> HiMassTauAnalysis::getMatchedPdgId(float pt, int charge){
  pair<unsigned int, unsigned int> theTrackAndMotherPdgId;
  float minDeltaPt = 1000.;
  float minDeltaR = 0.2;
  unsigned int thePdgId = 0;
  unsigned int theMotherPdgId = 0;
  for(GenParticleCollection::const_iterator genParticle = _genParticles->begin();genParticle != _genParticles->end();++genParticle){
    if(charge != genParticle->charge() || genParticle->status() != 1)continue;  // match only to final states...
    if(reco::deltaR(eta, phi, genParticle->eta(), genParticle->phi()) > minDeltaR) continue ;
    float theDeltaPt = fabs(pt - genParticle->pt());
    if(theDeltaPt < minDeltaPt){
      minDeltaPt = theDeltaPt;
      thePdgId = abs(genParticle->pdgId());
      theMotherPdgId = abs(genParticle->mother()->pdgId());
    }
  }
  theTrackAndMotherPdgId = make_pair<unsigned int, unsigned int>(thePdgId, theMotherPdgId);
  return theTrackAndMotherPdgId;
}

//-----Matching to generator level objects
//-----Electrons
pair<bool, math::XYZTLorentzVector> HiMassTauAnalysis::matchToGen(const pat::Electron& theObject) {
  bool isGenMatched = false;
  math::XYZTLorentzVector theGenObject(0,0,0,0);
  associatedGenParticles = theObject.genParticleRefs();
  for(std::vector<reco::GenParticleRef>::const_iterator it = associatedGenParticles.begin(); it != associatedGenParticles.end(); ++it){
    if ( it->ref().isNonnull() && it->ref().isValid() ) {
      const reco::GenParticleRef& genParticle = (*it);
      if (abs(genParticle->pdgId()) == 11) {
        if(genParticle->mother(0)->pdgId() == genParticle->pdgId()) {
          motherCand = genParticle->mother(0)->mother(0);
          if(motherCand->mother(0)->pdgId() == motherCand->pdgId()) {grandMotherCand = motherCand->mother(0)->mother(0);} 
          else {grandMotherCand = motherCand->mother(0);}
        } else {
          motherCand = genParticle->mother(0);
          if(motherCand->mother(0)->pdgId() == motherCand->pdgId()) {grandMotherCand = motherCand->mother(0)->mother(0);} 
          else {grandMotherCand = motherCand->mother(0);}
        }
        if(_UseLeptonMotherId) {
          if(abs(motherCand->pdgId()) == _LeptonMotherId) {
            if(_UseLeptonGrandMotherId) {
              if(abs(grandMotherCand->pdgId()) == _LeptonGrandMotherId) {isGenMatched = true;theGenObject = genParticle->p4();}
            } else {isGenMatched = true;theGenObject = genParticle->p4();}
          }
        } else {isGenMatched = true;theGenObject = genParticle->p4();}
      }
    }
  }
  pair<bool, math::XYZTLorentzVector> GenMatchedInformation(isGenMatched,theGenObject);
  return GenMatchedInformation;
}

//-----Matching to generator level objects
//-----Muons
pair<bool, math::XYZTLorentzVector> HiMassTauAnalysis::matchToGen(const pat::Muon& theObject) {
  bool isGenMatched = false;
  math::XYZTLorentzVector theGenObject(0,0,0,0);
  associatedGenParticles = theObject.genParticleRefs();
  for(std::vector<reco::GenParticleRef>::const_iterator it = associatedGenParticles.begin(); it != associatedGenParticles.end(); ++it){
    if ( it->ref().isNonnull() && it->ref().isValid() ) {
      const reco::GenParticleRef& genParticle = (*it);
      if (abs(genParticle->pdgId()) == 13) {
        if(genParticle->mother(0)->pdgId() == genParticle->pdgId()) {
          motherCand = genParticle->mother(0)->mother(0);
          if(motherCand->mother(0)->pdgId() == motherCand->pdgId()) {grandMotherCand = motherCand->mother(0)->mother(0);}
          else {grandMotherCand = motherCand->mother(0);}
        } else {
          motherCand = genParticle->mother(0);
          if(motherCand->mother(0)->pdgId() == motherCand->pdgId()) {grandMotherCand = motherCand->mother(0)->mother(0);}
          else {grandMotherCand = motherCand->mother(0);}
        }
        if(_UseLeptonMotherId) {
          if(abs(motherCand->pdgId()) == _LeptonMotherId) {
            if(_UseLeptonGrandMotherId) {
              if(abs(grandMotherCand->pdgId()) == _LeptonGrandMotherId) {isGenMatched = true;theGenObject = genParticle->p4();}
            } else {isGenMatched = true;theGenObject = genParticle->p4();}
          }
        } else {isGenMatched = true;theGenObject = genParticle->p4();}
      }
    }
  }
  pair<bool, math::XYZTLorentzVector> GenMatchedInformation(isGenMatched,theGenObject);
  return GenMatchedInformation;
}

//-----Matching to generator level objects
//-----Taus
pair<bool, math::XYZTLorentzVector> HiMassTauAnalysis::matchToGen(const pat::Tau& theObject) {
  bool isGenMatched = false;
  math::XYZTLorentzVector theGenObject(0,0,0,0);
  for(GenParticleCollection::const_iterator genParticle = _genParticles->begin();genParticle != _genParticles->end();++genParticle) {
    if((abs(genParticle->pdgId()) == 15) && (genParticle->status() != 3)) {
      int neutrinos = 0;
      MChadtau = genParticle->p4();
      if(genParticle->mother(0)->pdgId() == genParticle->pdgId()) {
        motherCand = genParticle->mother(0)->mother(0);
        if(motherCand->mother(0)->pdgId() == motherCand->pdgId()) {grandMotherCand = motherCand->mother(0)->mother(0);}
        else {grandMotherCand = motherCand->mother(0);}
      } else {
        motherCand = genParticle->mother(0);
        if(motherCand->mother(0)->pdgId() == motherCand->pdgId()) {grandMotherCand = motherCand->mother(0)->mother(0);}
        else {grandMotherCand = motherCand->mother(0);}
      }
      for(int ii=0; ii<(int)(genParticle->numberOfDaughters()); ii++) {
        daughterCand = genParticle->daughter(ii);
        if( (abs(daughterCand->pdgId()) == 12) || (abs(daughterCand->pdgId()) == 14) || (abs(daughterCand->pdgId()) == 16) ) {
          neutrinos++;
          MChadtau = MChadtau - daughterCand->p4();
        }
      }
      if(neutrinos == 1) {
        if(reco::deltaR(MChadtau.eta(), MChadtau.phi(), theObject.eta(), theObject.phi()) < _TauToGenMatchingDeltaR) {
          if(_UseTauMotherId) {
            if(abs(motherCand->pdgId()) == _TauMotherId) {
              if(_UseTauGrandMotherId) {
                if(abs(grandMotherCand->pdgId()) == _TauGrandMotherId) {isGenMatched = true;theGenObject = MChadtau;}
              } else {
                isGenMatched = true;
                theGenObject = MChadtau;
              }
            }
          } else {isGenMatched = true;theGenObject = MChadtau;}
        }
      }
    }
  }
  pair<bool, math::XYZTLorentzVector> GenMatchedInformation(isGenMatched,theGenObject);
  return GenMatchedInformation;
}

//-----Calculate zeta variables
double HiMassTauAnalysis::CalculatePZeta(const pat::Tau& patTau, const pat::Muon& patMuon, bool smear, reco::Candidate::LorentzVector smearedLV, const pat::MET& patMET) {
  if(smear) {
    double zetaX = cos(patTau.phi()) + cos(smearedLV.phi());
    double zetaY = sin(patTau.phi()) + sin(smearedLV.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = patTau.px() + smearedLV.px();
    double visPy = patTau.py() + smearedLV.py();
    double px = visPx + patMET.px();
    double py = visPy + patMET.py();
    double pZeta = px*zetaX + py*zetaY;
    return pZeta;
  } else {
    double zetaX = cos(patTau.phi()) + cos(patMuon.phi());
    double zetaY = sin(patTau.phi()) + sin(patMuon.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = patTau.px() + patMuon.px();
    double visPy = patTau.py() + patMuon.py();
    double px = visPx + patMET.px();
    double py = visPy + patMET.py();
    double pZeta = px*zetaX + py*zetaY;
    return pZeta;
  }
}
double HiMassTauAnalysis::CalculatePZeta(const pat::Tau& patTau, const pat::Electron& patElectron, bool smear, reco::Candidate::LorentzVector smearedLV, const pat::MET& patMET) {
  if(smear) {
    double zetaX = cos(patTau.phi()) + cos(smearedLV.phi());
    double zetaY = sin(patTau.phi()) + sin(smearedLV.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = patTau.px() + smearedLV.px();
    double visPy = patTau.py() + smearedLV.py();
    double px = visPx + patMET.px();
    double py = visPy + patMET.py();
    double pZeta = px*zetaX + py*zetaY;
    return pZeta;
  } else {
    double zetaX = cos(patTau.phi()) + cos(patElectron.phi());
    double zetaY = sin(patTau.phi()) + sin(patElectron.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = patTau.px() + patElectron.px();
    double visPy = patTau.py() + patElectron.py();
    double px = visPx + patMET.px();
    double py = visPy + patMET.py();
    double pZeta = px*zetaX + py*zetaY;
    return pZeta;
  }
}
double HiMassTauAnalysis::CalculatePZeta(const pat::Electron& patElectron, bool smearE, reco::Candidate::LorentzVector smearedLVE, const pat::Muon& patMuon, bool smearM, reco::Candidate::LorentzVector smearedLVM, const pat::MET& patMET) {
  if(smearE) {
    if(smearM) {
      double zetaX = cos(smearedLVM.phi()) + cos(smearedLVE.phi());
      double zetaY = sin(smearedLVM.phi()) + sin(smearedLVE.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = smearedLVM.px() + smearedLVE.px();
      double visPy = smearedLVM.py() + smearedLVE.py();
      double px = visPx + patMET.px();
      double py = visPy + patMET.py();
      double pZeta = px*zetaX + py*zetaY;
      return pZeta;
    } else {
      double zetaX = cos(patMuon.phi()) + cos(smearedLVE.phi());
      double zetaY = sin(patMuon.phi()) + sin(smearedLVE.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = patMuon.px() + smearedLVE.px();
      double visPy = patMuon.py() + smearedLVE.py();
      double px = visPx + patMET.px();
      double py = visPy + patMET.py();
      double pZeta = px*zetaX + py*zetaY;
      return pZeta;
    }
  } else {
    if(_SmearTheMuon) {
      double zetaX = cos(smearedLVM.phi()) + cos(patElectron.phi());
      double zetaY = sin(smearedLVM.phi()) + sin(patElectron.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = smearedLVM.px() + patElectron.px();
      double visPy = smearedLVM.py() + patElectron.py();
      double px = visPx + patMET.px();
      double py = visPy + patMET.py();
      double pZeta = px*zetaX + py*zetaY;
      return pZeta;
    } else {
      double zetaX = cos(patMuon.phi()) + cos(patElectron.phi());
      double zetaY = sin(patMuon.phi()) + sin(patElectron.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = patMuon.px() + patElectron.px();
      double visPy = patMuon.py() + patElectron.py();
      double px = visPx + patMET.px();
      double py = visPy + patMET.py();
      double pZeta = px*zetaX + py*zetaY;
      return pZeta;
    }
  }
}
double HiMassTauAnalysis::CalculatePZeta(const pat::Muon& patMuon1, bool smear1, reco::Candidate::LorentzVector smearedLV1, const pat::Muon& patMuon2, bool smear2, reco::Candidate::LorentzVector smearedLV2, const pat::MET& patMET) {
  if(smear1) {
    double zetaX = cos(smearedLV1.phi()) + cos(smearedLV2.phi());
    double zetaY = sin(smearedLV1.phi()) + sin(smearedLV2.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = smearedLV1.px() + smearedLV2.px();
    double visPy = smearedLV1.py() + smearedLV2.py();
    double px = visPx + patMET.px();
    double py = visPy + patMET.py();
    double pZeta = px*zetaX + py*zetaY;
    return pZeta;
  } else {
    double zetaX = cos(patMuon1.phi()) + cos(patMuon2.phi());
    double zetaY = sin(patMuon1.phi()) + sin(patMuon2.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = patMuon1.px() + patMuon2.px();
    double visPy = patMuon1.py() + patMuon2.py();
    double px = visPx + patMET.px();
    double py = visPy + patMET.py();
    double pZeta = px*zetaX + py*zetaY;
    return pZeta;
  }
}
double HiMassTauAnalysis::CalculatePZeta(const pat::Electron& patElectron1,  bool smear1, reco::Candidate::LorentzVector smearedLV1, const pat::Electron& patElectron2,  bool smear2, reco::Candidate::LorentzVector smearedLV2, const pat::MET& patMET) {
  if(smear1) {
    double zetaX = cos(smearedLV1.phi()) + cos(smearedLV2.phi());
    double zetaY = sin(smearedLV1.phi()) + sin(smearedLV2.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = smearedLV1.px() + smearedLV2.px();
    double visPy = smearedLV1.py() + smearedLV2.py();
    double px = visPx + patMET.px();
    double py = visPy + patMET.py();
    double pZeta = px*zetaX + py*zetaY;
    return pZeta;
  } else {
    double zetaX = cos(patElectron1.phi()) + cos(patElectron2.phi());
    double zetaY = sin(patElectron1.phi()) + sin(patElectron2.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = patElectron1.px() + patElectron2.px();
    double visPy = patElectron1.py() + patElectron2.py();
    double px = visPx + patMET.px();
    double py = visPy + patMET.py();
    double pZeta = px*zetaX + py*zetaY;
    return pZeta;
  }
}
double HiMassTauAnalysis::CalculatePZeta(const pat::Tau& patTau1, const pat::Tau& patTau2, const pat::MET& patMET) {
  double zetaX = cos(patTau1.phi()) + cos(patTau2.phi());
  double zetaY = sin(patTau1.phi()) + sin(patTau2.phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = patTau1.px() + patTau2.px();
  double visPy = patTau1.py() + patTau2.py();
  double px = visPx + patMET.px();
  double py = visPy + patMET.py();
  double pZeta = px*zetaX + py*zetaY;
  return pZeta;
}
double HiMassTauAnalysis::CalculatePZetaVis(const pat::Tau& patTau, const pat::Muon& patMuon, bool smear, reco::Candidate::LorentzVector smearedLV, const pat::MET& patMET) {
  if(smear) {
    double zetaX = cos(patTau.phi()) + cos(smearedLV.phi());
    double zetaY = sin(patTau.phi()) + sin(smearedLV.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = patTau.px() + smearedLV.px();
    double visPy = patTau.py() + smearedLV.py();
    double pZetaVis = visPx*zetaX + visPy*zetaY;
    return pZetaVis;
  } else {
    double zetaX = cos(patTau.phi()) + cos(patMuon.phi());
    double zetaY = sin(patTau.phi()) + sin(patMuon.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = patTau.px() + patMuon.px();
    double visPy = patTau.py() + patMuon.py();
    double pZetaVis = visPx*zetaX + visPy*zetaY;
    return pZetaVis;
  }
}
double HiMassTauAnalysis::CalculatePZetaVis(const pat::Tau& patTau, const pat::Electron& patElectron, bool smear, reco::Candidate::LorentzVector smearedLV, const pat::MET& patMET) {
  if(smear) {
    double zetaX = cos(patTau.phi()) + cos(smearedLV.phi());
    double zetaY = sin(patTau.phi()) + sin(smearedLV.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = patTau.px() + smearedLV.px();
    double visPy = patTau.py() + smearedLV.py();
    double pZetaVis = visPx*zetaX + visPy*zetaY;
    return pZetaVis;
  } else {
    double zetaX = cos(patTau.phi()) + cos(patElectron.phi());
    double zetaY = sin(patTau.phi()) + sin(patElectron.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = patTau.px() + patElectron.px();
    double visPy = patTau.py() + patElectron.py();
    double pZetaVis = visPx*zetaX + visPy*zetaY;
    return pZetaVis;
  }
}
double HiMassTauAnalysis::CalculatePZetaVis(const pat::Electron& patElectron, bool smearE, reco::Candidate::LorentzVector smearedLVE, const pat::Muon& patMuon, bool smearM, reco::Candidate::LorentzVector smearedLVM, const pat::MET& patMET) {
  if(smearE) {
    if(smearM) {
      double zetaX = cos(smearedLVM.phi()) + cos(smearedLVE.phi());
      double zetaY = sin(smearedLVM.phi()) + sin(smearedLVE.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = smearedLVM.px() + smearedLVE.px();
      double visPy = smearedLVM.py() + smearedLVE.py();
      double pZetaVis = visPx*zetaX + visPy*zetaY;
      return pZetaVis;
    } else {
      double zetaX = cos(patMuon.phi()) + cos(smearedLVE.phi());
      double zetaY = sin(patMuon.phi()) + sin(smearedLVE.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = patMuon.px() + smearedLVE.px();
      double visPy = patMuon.py() + smearedLVE.py();
      double pZetaVis = visPx*zetaX + visPy*zetaY;
      return pZetaVis;
    }
  } else {
    if(smearM) {
      double zetaX = cos(smearedLVM.phi()) + cos(patElectron.phi());
      double zetaY = sin(smearedLVM.phi()) + sin(patElectron.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = smearedLVM.px() + patElectron.px();
      double visPy = smearedLVM.py() + patElectron.py();
      double pZetaVis = visPx*zetaX + visPy*zetaY;
      return pZetaVis;
    } else {
      double zetaX = cos(patMuon.phi()) + cos(patElectron.phi());
      double zetaY = sin(patMuon.phi()) + sin(patElectron.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = patMuon.px() + patElectron.px();
      double visPy = patMuon.py() + patElectron.py();
      double pZetaVis = visPx*zetaX + visPy*zetaY;
      return pZetaVis;
    }
  }
}
double HiMassTauAnalysis::CalculatePZetaVis(const pat::Muon& patMuon1, bool smear1, reco::Candidate::LorentzVector smearedLV1, const pat::Muon& patMuon2, bool smear2, reco::Candidate::LorentzVector smearedLV2, const pat::MET& patMET) {
  if(smear1) {
    double zetaX = cos(smearedLV1.phi()) + cos(smearedLV2.phi());
    double zetaY = sin(smearedLV1.phi()) + sin(smearedLV2.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = smearedLV1.px() + smearedLV2.px();
    double visPy = smearedLV1.py() + smearedLV2.py();
    double pZetaVis = visPx*zetaX + visPy*zetaY;
    return pZetaVis;
  } else {
    double zetaX = cos(patMuon1.phi()) + cos(patMuon2.phi());
    double zetaY = sin(patMuon1.phi()) + sin(patMuon2.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = patMuon1.px() + patMuon2.px();
    double visPy = patMuon1.py() + patMuon2.py();
    double pZetaVis = visPx*zetaX + visPy*zetaY;
    return pZetaVis;
  }
}
double HiMassTauAnalysis::CalculatePZetaVis(const pat::Electron& patElectron1, bool smear1, reco::Candidate::LorentzVector smearedLV1, const pat::Electron& patElectron2, bool smear2, reco::Candidate::LorentzVector smearedLV2, const pat::MET& patMET) {
  if(smear1) {
    double zetaX = cos(smearedLV1.phi()) + cos(smearedLV2.phi());
    double zetaY = sin(smearedLV1.phi()) + sin(smearedLV2.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = smearedLV1.px() + smearedLV2.px();
    double visPy = smearedLV1.py() + smearedLV2.py();
    double pZetaVis = visPx*zetaX + visPy*zetaY;
    return pZetaVis;
  } else {
    double zetaX = cos(patElectron1.phi()) + cos(patElectron2.phi());
    double zetaY = sin(patElectron1.phi()) + sin(patElectron2.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = patElectron1.px() + patElectron2.px();
    double visPy = patElectron1.py() + patElectron2.py();
    double pZetaVis = visPx*zetaX + visPy*zetaY;
    return pZetaVis;
  }
}
double HiMassTauAnalysis::CalculatePZetaVis(const pat::Tau& patTau1, const pat::Tau& patTau2, const pat::MET& patMET) {
  double zetaX = cos(patTau1.phi()) + cos(patTau2.phi());
  double zetaY = sin(patTau1.phi()) + sin(patTau2.phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = patTau1.px() + patTau2.px();
  double visPy = patTau1.py() + patTau2.py();
  double pZetaVis = visPx*zetaX + visPy*zetaY;
  return pZetaVis;
}

//-----Calculate mass reco variables
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::CalculateThe4Momentum(const pat::Tau& patTau, const pat::Electron& patElectron, bool smear, reco::Candidate::LorentzVector smearedLV, const pat::MET& patMET) {
  if(smear) {
    if(_UseVectorSumOfVisProductsAndMetMassReco) {
      double px = patTau.px() + smearedLV.px() + patMET.px();
      double py = patTau.py() + smearedLV.py() + patMET.py();
      double pz = patTau.pz() + smearedLV.pz();
      double e = patTau.energy() + smearedLV.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else if(_UseCollinerApproxMassReco) {
      double x1_numerator = (patTau.px() * smearedLV.py()) - (smearedLV.px() * patTau.py());
      double x1_denominator = (smearedLV.py() * (patTau.px() + patMET.px())) - (smearedLV.px() * (patTau.py() + patMET.py()));
      double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
      double x2_numerator = x1_numerator;
      double x2_denominator = (patTau.px() * (smearedLV.py() + patMET.py())) - (patTau.py() * (smearedLV.px() + patMET.px()));
      double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
      if ( (x1 > 0. && x1 < 1.) &&
           (x2 > 0. && x2 < 1.) ) {
        reco::Candidate::LorentzVector The_LorentzVect = (patTau.p4() / x1) + (smearedLV / x2);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else {
        double px = patTau.px() + smearedLV.px() + patMET.px();
        double py = patTau.py() + smearedLV.py() + patMET.py();
        double pz = patTau.pz() + smearedLV.pz();
        double e = patTau.energy() + smearedLV.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
        return MassRecoInformation;
      }
    } else {
      reco::Candidate::LorentzVector The_LorentzVect = patTau.p4() + smearedLV;
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    }
  } else {
    if(_UseVectorSumOfVisProductsAndMetMassReco) {
      double px = patTau.px() + patElectron.px() + patMET.px();
      double py = patTau.py() + patElectron.py() + patMET.py();
      double pz = patTau.pz() + patElectron.pz();
      double e = patTau.energy() + patElectron.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else if(_UseCollinerApproxMassReco) {
      double x1_numerator = (patTau.px() * patElectron.py()) - (patElectron.px() * patTau.py());
      double x1_denominator = (patElectron.py() * (patTau.px() + patMET.px())) - (patElectron.px() * (patTau.py() + patMET.py()));
      double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
      double x2_numerator = x1_numerator;
      double x2_denominator = (patTau.px() * (patElectron.py() + patMET.py())) - (patTau.py() * (patElectron.px() + patMET.px()));
      double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
      if ( (x1 > 0. && x1 < 1.) &&
           (x2 > 0. && x2 < 1.) ) {
        reco::Candidate::LorentzVector The_LorentzVect = (patTau.p4() / x1) + (patElectron.p4() / x2);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else {
        double px = patTau.px() + patElectron.px() + patMET.px();
        double py = patTau.py() + patElectron.py() + patMET.py();
        double pz = patTau.pz() + patElectron.pz();
        double e = patTau.energy() + patElectron.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
        return MassRecoInformation;
      }
    } else {
      reco::Candidate::LorentzVector The_LorentzVect = patTau.p4() + patElectron.p4();
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    }
  }
}
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::CalculateThe4Momentum(const pat::Tau& patTau, const pat::Muon& patMuon, bool smear, reco::Candidate::LorentzVector smearedLV, const pat::MET& patMET) {
  if(smear) {
    if(_UseVectorSumOfVisProductsAndMetMassReco) {
      double px = patTau.px() + smearedLV.px() + patMET.px();
      double py = patTau.py() + smearedLV.py() + patMET.py();
      double pz = patTau.pz() + smearedLV.pz();
      double e = patTau.energy() + smearedLV.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else if(_UseCollinerApproxMassReco) {
      double x1_numerator = (patTau.px() * smearedLV.py()) - (smearedLV.px() * patTau.py());
      double x1_denominator = (smearedLV.py() * (patTau.px() + patMET.px())) - (smearedLV.px() * (patTau.py() + patMET.py()));
      double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
      double x2_numerator = x1_numerator;
      double x2_denominator = (patTau.px() * (smearedLV.py() + patMET.py())) - (patTau.py() * (smearedLV.px() + patMET.px()));
      double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
      if ( (x1 > 0. && x1 < 1.) &&
           (x2 > 0. && x2 < 1.) ) {
        reco::Candidate::LorentzVector The_LorentzVect = (patTau.p4() / x1) + (smearedLV / x2);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else {
        double px = patTau.px() + smearedLV.px() + patMET.px();
        double py = patTau.py() + smearedLV.py() + patMET.py();
        double pz = patTau.pz() + smearedLV.pz();
        double e = patTau.energy() + smearedLV.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
        return MassRecoInformation;
      }
    } else {
      reco::Candidate::LorentzVector The_LorentzVect = patTau.p4() + smearedLV;
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    }
  } else {
    if(_UseVectorSumOfVisProductsAndMetMassReco) {
      double px = patTau.px() + patMuon.px() + patMET.px();
      double py = patTau.py() + patMuon.py() + patMET.py();
      double pz = patTau.pz() + patMuon.pz();
      double e = patTau.energy() + patMuon.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else if(_UseCollinerApproxMassReco) {
      double x1_numerator = (patTau.px() * patMuon.py()) - (patMuon.px() * patTau.py());
      double x1_denominator = (patMuon.py() * (patTau.px() + patMET.px())) - (patMuon.px() * (patTau.py() + patMET.py()));
      double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
      double x2_numerator = x1_numerator;
      double x2_denominator = (patTau.px() * (patMuon.py() + patMET.py())) - (patTau.py() * (patMuon.px() + patMET.px()));
      double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
      if ( (x1 > 0. && x1 < 1.) &&
           (x2 > 0. && x2 < 1.) ) {
        reco::Candidate::LorentzVector The_LorentzVect = (patTau.p4() / x1) + (patMuon.p4() / x2);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else {
        double px = patTau.px() + patMuon.px() + patMET.px();
        double py = patTau.py() + patMuon.py() + patMET.py();
        double pz = patTau.pz() + patMuon.pz();
        double e = patTau.energy() + patMuon.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
        return MassRecoInformation;
      }
    } else {
      reco::Candidate::LorentzVector The_LorentzVect = patTau.p4() + patMuon.p4();
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    }
  }
}
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::CalculateThe4Momentum(const pat::Electron& patElectron, bool smearE, reco::Candidate::LorentzVector smearedLVE, const pat::Muon& patMuon, bool smearM, reco::Candidate::LorentzVector smearedLVM, const pat::MET& patMET) {
  if(smearE) {
    if(smearM) {
      if(_UseVectorSumOfVisProductsAndMetMassReco) {
        double px = smearedLVE.px() + smearedLVM.px() + patMET.px();
        double py = smearedLVE.py() + smearedLVM.py() + patMET.py();
        double pz = smearedLVE.pz() + smearedLVM.pz();
        double e = smearedLVE.energy() + smearedLVM.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else if(_UseCollinerApproxMassReco) {
        double x1_numerator = (smearedLVE.px() * smearedLVM.py()) - (smearedLVM.px() * smearedLVE.py());
        double x1_denominator = (smearedLVM.py() * (smearedLVE.px() + patMET.px())) - (smearedLVM.px() * (smearedLVE.py() + patMET.py()));
        double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
        double x2_numerator = x1_numerator;
        double x2_denominator = (smearedLVE.px() * (smearedLVM.py() + patMET.py())) - (smearedLVE.py() * (smearedLVM.px() + patMET.px()));
        double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
        if ( (x1 > 0. && x1 < 1.) &&
             (x2 > 0. && x2 < 1.) ) {
          reco::Candidate::LorentzVector The_LorentzVect = (smearedLVE / x1) + (smearedLVM / x2);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
          return MassRecoInformation;
        } else {
          double px = smearedLVE.px() + smearedLVM.px() + patMET.px();
          double py = smearedLVE.py() + smearedLVM.py() + patMET.py();
          double pz = smearedLVE.pz() + smearedLVM.pz();
          double e = smearedLVE.energy() + smearedLVM.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
          reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
          return MassRecoInformation;
        }
      } else {
        reco::Candidate::LorentzVector The_LorentzVect = smearedLVE + smearedLVM;
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      }
    } else {
      if(_UseVectorSumOfVisProductsAndMetMassReco) {
        double px = smearedLVE.px() + patMuon.px() + patMET.px();
        double py = smearedLVE.py() + patMuon.py() + patMET.py();
        double pz = smearedLVE.pz() + patMuon.pz();
        double e = smearedLVE.energy() + patMuon.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else if(_UseCollinerApproxMassReco) {
        double x1_numerator = (smearedLVE.px() * patMuon.py()) - (patMuon.px() * smearedLVE.py());
        double x1_denominator = (patMuon.py() * (smearedLVE.px() + patMET.px())) - (patMuon.px() * (smearedLVE.py() + patMET.py()));
        double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
        double x2_numerator = x1_numerator;
        double x2_denominator = (smearedLVE.px() * (patMuon.py() + patMET.py())) - (smearedLVE.py() * (patMuon.px() + patMET.px()));
        double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
        if ( (x1 > 0. && x1 < 1.) &&
             (x2 > 0. && x2 < 1.) ) {
          reco::Candidate::LorentzVector The_LorentzVect = (smearedLVE / x1) + (patMuon.p4() / x2);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
          return MassRecoInformation;
        } else {
          double px = smearedLVE.px() + patMuon.px() + patMET.px();
          double py = smearedLVE.py() + patMuon.py() + patMET.py();
          double pz = smearedLVE.pz() + patMuon.pz();
          double e = smearedLVE.energy() + patMuon.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
          reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
          return MassRecoInformation;
        }
      } else {
        reco::Candidate::LorentzVector The_LorentzVect = smearedLVE + patMuon.p4();
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      }
    }
  } else {
    if(smearM) {
      if(_UseVectorSumOfVisProductsAndMetMassReco) {
        double px = patElectron.px() + smearedLVM.px() + patMET.px();
        double py = patElectron.py() + smearedLVM.py() + patMET.py();
        double pz = patElectron.pz() + smearedLVM.pz();
        double e = patElectron.energy() + smearedLVM.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else if(_UseCollinerApproxMassReco) {
        double x1_numerator = (patElectron.px() * smearedLVM.py()) - (smearedLVM.px() * patElectron.py());
        double x1_denominator = (smearedLVM.py() * (patElectron.px() + patMET.px())) - (smearedLVM.px() * (patElectron.py() + patMET.py()));
        double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
        double x2_numerator = x1_numerator;
        double x2_denominator = (patElectron.px() * (smearedLVM.py() + patMET.py())) - (patElectron.py() * (smearedLVM.px() + patMET.px()));
        double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
        if ( (x1 > 0. && x1 < 1.) &&
             (x2 > 0. && x2 < 1.) ) {
          reco::Candidate::LorentzVector The_LorentzVect = (patElectron.p4() / x1) + (smearedLVM / x2);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
          return MassRecoInformation;
        } else {
          double px = patElectron.px() + smearedLVM.px() + patMET.px();
          double py = patElectron.py() + smearedLVM.py() + patMET.py();
          double pz = patElectron.pz() + smearedLVM.pz();
          double e = patElectron.energy() + smearedLVM.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
          reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
          return MassRecoInformation;
        }
      } else {
        reco::Candidate::LorentzVector The_LorentzVect = patElectron.p4() + smearedLVM;
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      }
    } else {
      if(_UseVectorSumOfVisProductsAndMetMassReco) {
        double px = patElectron.px() + patMuon.px() + patMET.px();
        double py = patElectron.py() + patMuon.py() + patMET.py();
        double pz = patElectron.pz() + patMuon.pz();
        double e = patElectron.energy() + patMuon.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else if(_UseCollinerApproxMassReco) {
        double x1_numerator = (patElectron.px() * patMuon.py()) - (patMuon.px() * patElectron.py());
        double x1_denominator = (patMuon.py() * (patElectron.px() + patMET.px())) - (patMuon.px() * (patElectron.py() + patMET.py()));
        double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
        double x2_numerator = x1_numerator;
        double x2_denominator = (patElectron.px() * (patMuon.py() + patMET.py())) - (patElectron.py() * (patMuon.px() + patMET.px()));
        double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
        if ( (x1 > 0. && x1 < 1.) &&
             (x2 > 0. && x2 < 1.) ) {
          reco::Candidate::LorentzVector The_LorentzVect = (patElectron.p4() / x1) + (patMuon.p4() / x2);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
          return MassRecoInformation;
        } else {
          double px = patElectron.px() + patMuon.px() + patMET.px();
          double py = patElectron.py() + patMuon.py() + patMET.py();
          double pz = patElectron.pz() + patMuon.pz();
          double e = patElectron.energy() + patMuon.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
          reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
          return MassRecoInformation;
        }
      } else {
        reco::Candidate::LorentzVector The_LorentzVect = patElectron.p4() + patMuon.p4();
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      }
    }
  }
}
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::CalculateThe4Momentum(const pat::Muon& patMuon1, bool smear1, reco::Candidate::LorentzVector smearedLV1, const pat::Muon& patMuon2, bool smear2, reco::Candidate::LorentzVector smearedLV2, const pat::MET& patMET) {
  if(smear1) {
    if(_UseVectorSumOfVisProductsAndMetMassReco) {
      double px = smearedLV1.px() + smearedLV2.px() + patMET.px();
      double py = smearedLV1.py() + smearedLV2.py() + patMET.py();
      double pz = smearedLV1.pz() + smearedLV2.pz();
      double e = smearedLV1.energy() + smearedLV2.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else if(_UseCollinerApproxMassReco) {
      double x1_numerator = (smearedLV1.px() * smearedLV2.py()) - (smearedLV2.px() * smearedLV1.py());
      double x1_denominator = (smearedLV2.py() * (smearedLV1.px() + patMET.px())) - (smearedLV2.px() * (smearedLV1.py() + patMET.py()));
      double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
      double x2_numerator = x1_numerator;
      double x2_denominator = (smearedLV1.px() * (smearedLV2.py() + patMET.py())) - (smearedLV1.py() * (smearedLV2.px() + patMET.px()));
      double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
      if ( (x1 > 0. && x1 < 1.) &&
           (x2 > 0. && x2 < 1.) ) {
        reco::Candidate::LorentzVector The_LorentzVect = (smearedLV1 / x1) + (smearedLV2 / x2);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else {
        double px = smearedLV1.px() + smearedLV2.px() + patMET.px();
        double py = smearedLV1.py() + smearedLV2.py() + patMET.py();
        double pz = smearedLV1.pz() + smearedLV2.pz();
        double e = smearedLV1.energy() + smearedLV2.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
        return MassRecoInformation;
      }
    } else {
      reco::Candidate::LorentzVector The_LorentzVect = smearedLV1 + smearedLV2;
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    }
  } else {
    if(_UseVectorSumOfVisProductsAndMetMassReco) {
      double px = patMuon1.px() + patMuon2.px() + patMET.px();
      double py = patMuon1.py() + patMuon2.py() + patMET.py();
      double pz = patMuon1.pz() + patMuon2.pz();
      double e = patMuon1.energy() + patMuon2.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else if(_UseCollinerApproxMassReco) {
      double x1_numerator = (patMuon1.px() * patMuon2.py()) - (patMuon2.px() * patMuon1.py());
      double x1_denominator = (patMuon2.py() * (patMuon1.px() + patMET.px())) - (patMuon2.px() * (patMuon1.py() + patMET.py()));
      double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
      double x2_numerator = x1_numerator;
      double x2_denominator = (patMuon1.px() * (patMuon2.py() + patMET.py())) - (patMuon1.py() * (patMuon2.px() + patMET.px()));
      double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
      if ( (x1 > 0. && x1 < 1.) &&
           (x2 > 0. && x2 < 1.) ) {
        reco::Candidate::LorentzVector The_LorentzVect = (patMuon1.p4() / x1) + (patMuon2.p4() / x2);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else {
        double px = patMuon1.px() + patMuon2.px() + patMET.px();
        double py = patMuon1.py() + patMuon2.py() + patMET.py();
        double pz = patMuon1.pz() + patMuon2.pz();
        double e = patMuon1.energy() + patMuon2.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
        return MassRecoInformation;
      }
    } else {
      reco::Candidate::LorentzVector The_LorentzVect = patMuon1.p4() + patMuon2.p4();
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    }
  }
}
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::CalculateThe4Momentum(const pat::Electron& patElectron1, bool smear1, reco::Candidate::LorentzVector smearedLV1, const pat::Electron& patElectron2, bool smear2, reco::Candidate::LorentzVector smearedLV2, const pat::MET& patMET) {
  if(smear1) {
    if(_UseVectorSumOfVisProductsAndMetMassReco) {
      double px = smearedLV1.px() + smearedLV2.px() + patMET.px();
      double py = smearedLV1.py() + smearedLV2.py() + patMET.py();
      double pz = smearedLV1.pz() + smearedLV2.pz();
      double e = smearedLV1.energy() + smearedLV2.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else if(_UseCollinerApproxMassReco) {
      double x1_numerator = (smearedLV1.px() * smearedLV2.py()) - (smearedLV2.px() * smearedLV1.py());
      double x1_denominator = (smearedLV2.py() * (smearedLV1.px() + patMET.px())) - (smearedLV2.px() * (smearedLV1.py() + patMET.py()));
      double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
      double x2_numerator = x1_numerator;
      double x2_denominator = (smearedLV1.px() * (smearedLV2.py() + patMET.py())) - (smearedLV1.py() * (smearedLV2.px() + patMET.px()));
      double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
      if ( (x1 > 0. && x1 < 1.) &&
           (x2 > 0. && x2 < 1.) ) {
        reco::Candidate::LorentzVector The_LorentzVect = (smearedLV1 / x1) + (smearedLV2 / x2);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else {
        double px = smearedLV1.px() + smearedLV2.px() + patMET.px();
        double py = smearedLV1.py() + smearedLV2.py() + patMET.py();
        double pz = smearedLV1.pz() + smearedLV2.pz();
        double e = smearedLV1.energy() + smearedLV2.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
        return MassRecoInformation;
      }
    } else {
      reco::Candidate::LorentzVector The_LorentzVect = smearedLV1 + smearedLV2;
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    }
  } else {
    if(_UseVectorSumOfVisProductsAndMetMassReco) {
      double px = patElectron1.px() + patElectron2.px() + patMET.px();
      double py = patElectron1.py() + patElectron2.py() + patMET.py();
      double pz = patElectron1.pz() + patElectron2.pz();
      double e = patElectron1.energy() + patElectron2.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else if(_UseCollinerApproxMassReco) {
      double x1_numerator = (patElectron1.px() * patElectron2.py()) - (patElectron2.px() * patElectron1.py());
      double x1_denominator = (patElectron2.py() * (patElectron1.px() + patMET.px())) - (patElectron2.px() * (patElectron1.py() + patMET.py()));
      double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
      double x2_numerator = x1_numerator;
      double x2_denominator = (patElectron1.px() * (patElectron2.py() + patMET.py())) - (patElectron1.py() * (patElectron2.px() + patMET.px()));
      double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
      if ( (x1 > 0. && x1 < 1.) &&
           (x2 > 0. && x2 < 1.) ) {
        reco::Candidate::LorentzVector The_LorentzVect = (patElectron1.p4() / x1) + (patElectron2.p4() / x2);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else {
        double px = patElectron1.px() + patElectron2.px() + patMET.px();
        double py = patElectron1.py() + patElectron2.py() + patMET.py();
        double pz = patElectron1.pz() + patElectron2.pz();
        double e = patElectron1.energy() + patElectron2.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
        return MassRecoInformation;
      }
    } else {
      reco::Candidate::LorentzVector The_LorentzVect = patElectron1.p4() + patElectron2.p4();
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    }
  }
}
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::CalculateThe4Momentum(const pat::Tau& patTau1, const pat::Tau& patTau2, const pat::MET& patMET) {
  if(_UseVectorSumOfVisProductsAndMetMassReco) {
    double px = patTau1.px() + patTau2.px() + patMET.px();
    double py = patTau1.py() + patTau2.py() + patMET.py();
    double pz = patTau1.pz() + patTau2.pz();
    double e = patTau1.energy() + patTau2.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
    reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
    pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
    return MassRecoInformation;
  } else if(_UseCollinerApproxMassReco) {
    double x1_numerator = (patTau1.px() * patTau2.py()) - (patTau2.px() * patTau1.py());
    double x1_denominator = (patTau2.py() * (patTau1.px() + patMET.px())) - (patTau2.px() * (patTau1.py() + patMET.py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (patTau1.px() * (patTau2.py() + patMET.py())) - (patTau1.py() * (patTau2.px() + patMET.px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    if ( (x1 > 0. && x1 < 1.) &&
         (x2 > 0. && x2 < 1.) ) {
      reco::Candidate::LorentzVector The_LorentzVect = (patTau1.p4() / x1) + (patTau2.p4() / x2);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else {
      double px = patTau1.px() + patTau2.px() + patMET.px();
      double py = patTau1.py() + patTau2.py() + patMET.py();
      double pz = patTau1.pz() + patTau2.pz();
      double e = patTau1.energy() + patTau2.energy() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
      return MassRecoInformation;
    }
  } else {
    reco::Candidate::LorentzVector The_LorentzVect = patTau1.p4() + patTau2.p4();
    pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
    return MassRecoInformation;
  }
}

//-----Calculate lepton+met transverse mass
double HiMassTauAnalysis::CalculateLeptonMetMt(const pat::Muon& patMuon, bool smear, reco::Candidate::LorentzVector smearedLV, const pat::MET& patMET) {
  if(smear) {
    double px = smearedLV.px() + patMET.px();
    double py = smearedLV.py() + patMET.py();
    double et = smearedLV.Et() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
    double mt2 = et*et - (px*px + py*py);
    if ( mt2 < 0 ) { return -1.; }
    else { return sqrt(mt2); }
  } else {
    double px = patMuon.px() + patMET.px();
    double py = patMuon.py() + patMET.py();
    double et = patMuon.et() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
    double mt2 = et*et - (px*px + py*py);
    if ( mt2 < 0 ) { return -1.; }
    else { return sqrt(mt2); }
  }
}
double HiMassTauAnalysis::CalculateLeptonMetMt(const pat::Electron& patElectron, bool smear, reco::Candidate::LorentzVector smearedLV, const pat::MET& patMET) {
  if(smear) {
    double px = smearedLV.px() + patMET.px();
    double py = smearedLV.py() + patMET.py();
    double et = smearedLV.Et() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
    double mt2 = et*et - (px*px + py*py);
    if ( mt2 < 0 ) { return -1.; }
    else { return sqrt(mt2); }
  } else {
    double px = patElectron.px() + patMET.px();
    double py = patElectron.py() + patMET.py();
    double et = patElectron.et() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
    double mt2 = et*et - (px*px + py*py);
    if ( mt2 < 0 ) { return -1.; }
    else { return sqrt(mt2); }
  }
}
double HiMassTauAnalysis::CalculateLeptonMetMt(const pat::Tau& patTau, const pat::MET& patMET) {
  double px = patTau.px() + patMET.px();
  double py = patTau.py() + patMET.py();
  double et = patTau.et() + TMath::Sqrt((patMET.px() * patMET.px()) + (patMET.py() * patMET.py()));
  double mt2 = et*et - (px*px + py*py);
  if ( mt2 < 0 ) { return -1.; }
  else { return sqrt(mt2); }
}

//-----Calculate met correction
/*
double HiMassTauAnalysis::CorrectTheMet(const pat::Tau& patTau) {
  muonsrecoForMetCorrections = *(_recoMuonsForMetCorrections);
  vm_muCorrData = *(vm_muCorrData_h.product());
  muonsForMetCorrections = *(_patMuonsForMetCorrections);
  float px;
  float deltax;
  float corMetX;
  uint nMuons = muonsForMetCorrections.size();
  for(unsigned int iMu = 0; iMu<nMuons; iMu++) {
    if(muonsForMetCorrections[iMu].isGlobalMuon()) {
      uint nnMuons = muonsrecoForMetCorrections.size();
      for(unsigned int iiMu = 0; iiMu<nnMuons; iiMu++) {
        const reco::Muon *mu = &muonsrecoForMetCorrections[iiMu];
        if(muonsrecoForMetCorrections[iiMu].isGlobalMuon()) {
          px = muonsForMetCorrections[iMu].px();
          std::cout << "pt of pat muon  : " << px << std::endl;
          std::cout << "pt of reco muon : " << muonsrecoForMetCorrections[iiMu].px() << std::endl;
          reco::MuonMETCorrectionData muCorrData = (vm_muCorrData)[muonsrecoForMetCorrections.refAt(iiMu)];
          deltax = muCorrData.corrX();
          std::cout << "dx of muon: " << deltax << std::endl;
        }
      }
      corMetX = - px + deltax;
    }
  }
}
*/

//-----Calculate Tau Isolation Quantities
pair<int, double> HiMassTauAnalysis::CalculateTauTrackIsolation(const pat::Tau& patTau) {
  int nIsoTrks=0;
  double sumPtIsoTrks=0;
  if (patTau.isCaloTau()) {
    if (patTau.leadTrack().isNonnull()) {
      TrackRefVector TauIsolTracks = patTau.isolationTracks();
      for(TrackRefVector::const_iterator iTrk=TauIsolTracks.begin();iTrk!=TauIsolTracks.end();++iTrk) {
        if((reco::deltaR((**iTrk).eta(),(**iTrk).phi(),patTau.leadTrack()->eta(),patTau.leadTrack()->phi())<_RecoTauIsoDeltaRCone) && ((**iTrk).pt()>_RecoTauTrackIsoTrkThreshold)) {
          nIsoTrks++;
          sumPtIsoTrks = sumPtIsoTrks + (**iTrk).pt();
        }
      }
    }
  } else {
    if (patTau.leadPFChargedHadrCand().isNonnull()) { 
      PFCandidateRefVector TauIsolTracks = patTau.isolationPFChargedHadrCands();
      for(PFCandidateRefVector::const_iterator iTrk=TauIsolTracks.begin();iTrk!=TauIsolTracks.end();++iTrk) {
        if((reco::deltaR((**iTrk).eta(),(**iTrk).phi(),patTau.leadPFChargedHadrCand()->eta(),patTau.leadPFChargedHadrCand()->phi())<_RecoTauIsoDeltaRCone) && ((**iTrk).pt()>_RecoTauTrackIsoTrkThreshold)) {
          nIsoTrks++;
          sumPtIsoTrks = sumPtIsoTrks + (**iTrk).pt();
        }
      }
    }
  }
  pair<int, double> IsoTrksInformation(nIsoTrks,sumPtIsoTrks);
  return IsoTrksInformation;
}
pair<int, double> HiMassTauAnalysis::CalculateTauTrackIsolation(const pat::Tau& patTau, float deltaRCone, float trkMinPt) {
  int nIsoTrks=0;
  double sumPtIsoTrks=0;
  if (patTau.isCaloTau()) {
    if (patTau.leadTrack().isNonnull()) {
      TrackRefVector TauIsolTracks = patTau.isolationTracks();
      for(TrackRefVector::const_iterator iTrk=TauIsolTracks.begin();iTrk!=TauIsolTracks.end();++iTrk) {
        if((reco::deltaR((**iTrk).eta(),(**iTrk).phi(),patTau.leadTrack()->eta(),patTau.leadTrack()->phi())< deltaRCone) && ((**iTrk).pt()> trkMinPt)) {
          nIsoTrks++;
          sumPtIsoTrks = sumPtIsoTrks + (**iTrk).pt();
        }
      }
    }
  } else {
    if (patTau.leadPFChargedHadrCand().isNonnull()) { 
      PFCandidateRefVector TauIsolTracks = patTau.isolationPFChargedHadrCands();
      for(PFCandidateRefVector::const_iterator iTrk=TauIsolTracks.begin();iTrk!=TauIsolTracks.end();++iTrk) {
        if((reco::deltaR((**iTrk).eta(),(**iTrk).phi(),patTau.leadPFChargedHadrCand()->eta(),patTau.leadPFChargedHadrCand()->phi())< deltaRCone) &&((**iTrk).pt()> trkMinPt)) {
          nIsoTrks++;
          sumPtIsoTrks = sumPtIsoTrks + (**iTrk).pt();
        }
      }
    }
  }
  pair<int, double> IsoTrksInformation(nIsoTrks,sumPtIsoTrks);
  return IsoTrksInformation;
}

pair<int, double> HiMassTauAnalysis::CalculateTauEcalIsolation(const pat::Tau& patTau) {
  int nIsoGams=0;
  double sumPtIsoGams=0;
  if (_UseRecoTauEllipseForEcalIso) {
    if (patTau.isCaloTau()) {}
    else {
      if (patTau.leadPFChargedHadrCand().isNonnull()) { 
        PFCandidateRefVector TauIsolGammas = patTau.isolationPFGammaCands();
        for(PFCandidateRefVector::const_iterator iGam=TauIsolGammas.begin();iGam!=TauIsolGammas.end();++iGam) {
          if((reco::deltaR((**iGam).eta(),(**iGam).phi(),patTau.leadPFChargedHadrCand()->eta(),patTau.leadPFChargedHadrCand()->phi())<_RecoTauIsoDeltaRCone) && ((**iGam).pt()>_RecoTauGammaIsoGamThreshold)) {
            double deltaphi = TMath::Abs(normalizedPhi((**iGam).phi() - patTau.leadPFChargedHadrCand()->phi()));
            double deltaeta = (**iGam).eta() - patTau.leadPFChargedHadrCand()->eta();
            if( ( ((deltaphi * deltaphi) / (_RecoTauEcalIsoRphiForEllipse * _RecoTauEcalIsoRphiForEllipse)) + 
                  ((deltaeta * deltaeta) / (_RecoTauEcalIsoRetaForEllipse * _RecoTauEcalIsoRetaForEllipse)) ) > 1.0 ) {
              nIsoGams++;
              sumPtIsoGams = sumPtIsoGams + (**iGam).pt();
            }
          }
        }
        PFCandidateRefVector TauSigGammas = patTau.signalPFGammaCands();
        for(PFCandidateRefVector::const_iterator iGam=TauSigGammas.begin();iGam!=TauSigGammas.end();++iGam) {
          if((reco::deltaR((**iGam).eta(),(**iGam).phi(),patTau.leadPFChargedHadrCand()->eta(),patTau.leadPFChargedHadrCand()->phi())<_RecoTauIsoDeltaRCone) && ((**iGam).pt()>_RecoTauGammaIsoGamThreshold)) {
            double deltaphi = TMath::Abs(normalizedPhi((**iGam).phi() - patTau.leadPFChargedHadrCand()->phi()));
            double deltaeta = (**iGam).eta() - patTau.leadPFChargedHadrCand()->eta();
            if( ( ((deltaphi * deltaphi) / (_RecoTauEcalIsoRphiForEllipse * _RecoTauEcalIsoRphiForEllipse)) + 
                  ((deltaeta * deltaeta) / (_RecoTauEcalIsoRetaForEllipse * _RecoTauEcalIsoRetaForEllipse)) ) > 1.0 ) {
              nIsoGams++;
              sumPtIsoGams = sumPtIsoGams + (**iGam).pt();
            }
          }
        }
      }
    }
  } else {
    if (patTau.isCaloTau()) {}
    else {
      if (patTau.leadPFChargedHadrCand().isNonnull()) { 
        PFCandidateRefVector TauIsolGammas = patTau.isolationPFGammaCands();
        for(PFCandidateRefVector::const_iterator iGam=TauIsolGammas.begin();iGam!=TauIsolGammas.end();++iGam) {
          if((reco::deltaR((**iGam).eta(),(**iGam).phi(),patTau.leadPFChargedHadrCand()->eta(),patTau.leadPFChargedHadrCand()->phi())<_RecoTauIsoDeltaRCone) && ((**iGam).pt()>_RecoTauGammaIsoGamThreshold)) {
            nIsoGams++;
            sumPtIsoGams = sumPtIsoGams + (**iGam).pt();
          }
        }
      }
    }
  }
  pair<int, double> IsoGamsInformation(nIsoGams,sumPtIsoGams);
  return IsoGamsInformation;
}
pair<int, double> HiMassTauAnalysis::CalculateTauEcalIsolation(const pat::Tau& patTau, float deltaRCone, float gammaMinPt) {
  int nIsoGams=0;
  double sumPtIsoGams=0;
  if (_UseRecoTauEllipseForEcalIso) {
    if (patTau.isCaloTau()) {}
    else {
      if (patTau.leadPFChargedHadrCand().isNonnull()) { 
        PFCandidateRefVector TauIsolGammas = patTau.isolationPFGammaCands();
        for(PFCandidateRefVector::const_iterator iGam=TauIsolGammas.begin();iGam!=TauIsolGammas.end();++iGam) {
          if((reco::deltaR((**iGam).eta(),(**iGam).phi(),patTau.leadPFChargedHadrCand()->eta(),patTau.leadPFChargedHadrCand()->phi())<deltaRCone) && ((**iGam).pt()>gammaMinPt)) {
            double deltaphi = TMath::Abs(normalizedPhi((**iGam).phi() - patTau.leadPFChargedHadrCand()->phi()));
            double deltaeta = (**iGam).eta() - patTau.leadPFChargedHadrCand()->eta();
            if( ( ((deltaphi * deltaphi) / (_RecoTauEcalIsoRphiForEllipse * _RecoTauEcalIsoRphiForEllipse)) + 
                  ((deltaeta * deltaeta) / (_RecoTauEcalIsoRetaForEllipse * _RecoTauEcalIsoRetaForEllipse)) ) > 1.0 ) {
              nIsoGams++;
              sumPtIsoGams = sumPtIsoGams + (**iGam).pt();
            }
          }
        }
        PFCandidateRefVector TauSigGammas = patTau.signalPFGammaCands();
        for(PFCandidateRefVector::const_iterator iGam=TauSigGammas.begin();iGam!=TauSigGammas.end();++iGam) {
          if((reco::deltaR((**iGam).eta(),(**iGam).phi(),patTau.leadPFChargedHadrCand()->eta(),patTau.leadPFChargedHadrCand()->phi())<deltaRCone) && ((**iGam).pt()>gammaMinPt)) {
            double deltaphi = TMath::Abs(normalizedPhi((**iGam).phi() - patTau.leadPFChargedHadrCand()->phi()));
            double deltaeta = (**iGam).eta() - patTau.leadPFChargedHadrCand()->eta();
            if( ( ((deltaphi * deltaphi) / (_RecoTauEcalIsoRphiForEllipse * _RecoTauEcalIsoRphiForEllipse)) + 
                  ((deltaeta * deltaeta) / (_RecoTauEcalIsoRetaForEllipse * _RecoTauEcalIsoRetaForEllipse)) ) > 1.0 ) {
              nIsoGams++;
              sumPtIsoGams = sumPtIsoGams + (**iGam).pt();
            }
          }
        }
      }
    }
  } else {
    if (patTau.isCaloTau()) {}
    else {
      if (patTau.leadPFChargedHadrCand().isNonnull()) { 
        PFCandidateRefVector TauIsolGammas = patTau.isolationPFGammaCands();
        for(PFCandidateRefVector::const_iterator iGam=TauIsolGammas.begin();iGam!=TauIsolGammas.end();++iGam) {
          if((reco::deltaR((**iGam).eta(),(**iGam).phi(),patTau.leadPFChargedHadrCand()->eta(),patTau.leadPFChargedHadrCand()->phi())<deltaRCone) && ((**iGam).pt()>gammaMinPt)) {
            nIsoGams++;
            sumPtIsoGams = sumPtIsoGams + (**iGam).pt();
          }
        }
      }
    }
  }
  pair<int, double> IsoGamsInformation(nIsoGams,sumPtIsoGams);
  return IsoGamsInformation;
}

//-----Calculate Leading Gamma Information
int HiMassTauAnalysis::CalculateNumberSignalTauGammas(const pat::Tau& patTau) {
  int nSigGams=0;
  if (patTau.isCaloTau()) {
  } else {
    if (patTau.leadPFChargedHadrCand().isNonnull()) { 
      PFCandidateRefVector TauSigGammas = patTau.signalPFGammaCands();
      for(PFCandidateRefVector::const_iterator iGam=TauSigGammas.begin();iGam!=TauSigGammas.end();++iGam) {
        if((**iGam).pt()>_RecoTauSigGamThreshold) {
          nSigGams++;
        }
      }
    }
  }
  return nSigGams;
}

//-----Calculate Mass from Tau signal constituents
reco::Candidate::LorentzVector HiMassTauAnalysis::CalculateTauSignalTracksMass(const pat::Tau& patTau) {
  double px=0;
  double py=0;
  double pz=0;
  double e=0;
  if (patTau.isCaloTau()) {
    TrackRefVector TauSigTracks = patTau.signalTracks();
    for(TrackRefVector::const_iterator iTrk=TauSigTracks.begin();iTrk!=TauSigTracks.end();++iTrk) {
      px += (**iTrk).momentum().x(),
      py += (**iTrk).momentum().y(),
      pz += (**iTrk).momentum().z(),
      e += sqrt(pow((double)(**iTrk).momentum().r(),2)+pow(0.13957018,2));
    }
  } else {
    PFCandidateRefVector TauSigTracks = patTau.signalPFChargedHadrCands();
    for(PFCandidateRefVector::const_iterator iTrk=TauSigTracks.begin();iTrk!=TauSigTracks.end();++iTrk) {
      px += (**iTrk).momentum().x(),
      py += (**iTrk).momentum().y(),
      pz += (**iTrk).momentum().z(),
      e += sqrt(pow((double)(**iTrk).momentum().r(),2)+pow(0.13957018,2));
    }
  }
  reco::Candidate::LorentzVector TheSignalTracks_LorentzVect(px, py, pz, e);
  return TheSignalTracks_LorentzVect;
}

//-----Smear the light leptons (for studies of systematic uncertanties)
reco::Candidate::LorentzVector HiMassTauAnalysis::SmearLightLepton(const pat::Muon& patMuon) {
  if(matchToGen(patMuon).first) {
    math::XYZTLorentzVector unsmearedMomentum = matchToGen(patMuon).second;
    double smearedPt  = (unsmearedMomentum.pt() * CLHEP::RandGauss::shoot(_RelativeMuonPtOffset,_RelativeMuonPtSigma)) + unsmearedMomentum.pt();
    double smearedEta = CLHEP::RandGauss::shoot(_AbsoluteMuonEtaOffset,_AbsoluteMuonEtaSigma) + unsmearedMomentum.eta();
    double smearedPhi = CLHEP::RandGauss::shoot(_AbsoluteMuonPhiOffset,_AbsoluteMuonPhiSigma) + unsmearedMomentum.phi();
    math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector(smearedPt, smearedEta, smearedPhi, unsmearedMomentum.mass());
    reco::Candidate::LorentzVector smearedMomentum(smearedPtEtaPhiMVector.px(), smearedPtEtaPhiMVector.py(), smearedPtEtaPhiMVector.pz(), smearedPtEtaPhiMVector.energy());
    return smearedMomentum;
  } else {
    reco::Candidate::LorentzVector smearedMomentum = patMuon.p4();
    return smearedMomentum;
  }
}
reco::Candidate::LorentzVector HiMassTauAnalysis::SmearLightLepton(const pat::Electron& patElectron) {
  if(matchToGen(patElectron).first) {
    math::XYZTLorentzVector unsmearedMomentum = matchToGen(patElectron).second;
    double smearedPt  = (unsmearedMomentum.pt() * CLHEP::RandGauss::shoot(_RelativeElectronPtOffset,_RelativeElectronPtSigma)) + unsmearedMomentum.pt();
    double smearedEta = CLHEP::RandGauss::shoot(_AbsoluteElectronEtaOffset,_AbsoluteElectronEtaSigma) + unsmearedMomentum.eta();
    double smearedPhi = CLHEP::RandGauss::shoot(_AbsoluteElectronPhiOffset,_AbsoluteElectronPhiSigma) + unsmearedMomentum.phi();
    math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector(smearedPt, smearedEta, smearedPhi, unsmearedMomentum.mass());
    reco::Candidate::LorentzVector smearedMomentum(smearedPtEtaPhiMVector.px(), smearedPtEtaPhiMVector.py(), smearedPtEtaPhiMVector.pz(), smearedPtEtaPhiMVector.energy());
    return smearedMomentum;
  } else {
    reco::Candidate::LorentzVector smearedMomentum = patElectron.p4();
    return smearedMomentum;
  }
}

//-----Initialize information for the calculation of systematic uncertaintites
void HiMassTauAnalysis::InitializeInfoForPDFSystematicUncertaintites() {
  for (unsigned int i=0; i<pdfWeightTags_.size(); ++i) {
    std::cout << "\t" << pdfWeightTags_[i].label();
    pdfStart_Denominator_.push_back(-1);
    pdfStart_Numerator_.push_back(-1);
  }
}

bool HiMassTauAnalysis::isInTheCracks(float etaValue){
  return (fabs(etaValue) < 0.018 || 
	(fabs(etaValue)>0.423 && fabs(etaValue)<0.461) ||
	(fabs(etaValue)>0.770 && fabs(etaValue)<0.806) ||
	(fabs(etaValue)>1.127 && fabs(etaValue)<1.163) ||
	(fabs(etaValue)>1.460 && fabs(etaValue)<1.558));
}

// ---------------
void HiMassTauAnalysis::bookHistograms() {
  Service<TFileService> fs;

  _hEvents     = fs->make<TH1F>("Events", "Events", 2, 0., 2.);

  //--- book vertex histograms
  if (_FillRecoVertexHists) {
    _hVertexZposition = fs->make<TH1F>("VertexZposition", "VertexZposition", 50, -50., 50.);
    _hVertexNTracks   = fs->make<TH1F>("VertexNTracks", "VertexNTracks", 100, 0., 100.);
    _hNVertices       = fs->make<TH1F>("NVertices", "NVertices", 10, 0., 10.);
  }
 
  //--- book generator level histograms
  if (_FillGenTauHists) {
    _hNGenTau = fs->make<TH1F>("NGenTau", "NGenTau", 21, 0., 20.);
    _hGenTauEnergy = fs->make<TH1F>("GenTauEnergy", "GenTauEnergy", 200, 0., 500.);
    _hGenTauPt     = fs->make<TH1F>("GenTauPt", "GenTauPt", 200, 0., 500.);
    _hGenTauEta    = fs->make<TH1F>("GenTauEta", "GenTauEta", 72, -3.6, +3.6);
    _hGenTauPhi    = fs->make<TH1F>("GenTauPhi", "GenTauPhi", 36, -TMath::Pi(), +TMath::Pi());
  }
  
  //--- book reconstruction level histograms 
  if (_FillRecoTauHists) {
    _hNTau = fs->make<TH1F>("NTau", "NTau", 21, 0., 20.);
    _hTauJetEnergy = fs->make<TH1F>("TauJetEnergy", "TauJetEnergy", 200, 0., 500.);
    _hTauJetPt     = fs->make<TH1F>("TauJetPt", "TauJetPt", 200, 0., 500.);
    _hTauJetEta    = fs->make<TH1F>("TauJetEta", "TauJetEta", 72, -3.6, +3.6);
    _hTauJetPhi    = fs->make<TH1F>("TauJetPhi", "TauJetPhi", 36, -TMath::Pi(), +TMath::Pi());
    _hTauJetNumSignalTracks    = fs->make<TH1F>("TauJetNumSignalTracks", "TauJetNumSignalTracks", 10, 0, 10);
    _hTauJetNumSignalGammas    = fs->make<TH1F>("TauJetNumSignalGammas", "TauJetNumSignalGammas", 10, 0, 10);
    _hTauJetSeedTrackPt        = fs->make<TH1F>("TauJetSeedTrackPt", "TauJetSeedTrackPt", 200, 0., 500.);
    _hTauJetSeedTrackIpSignificance        = fs->make<TH1F>("TauJetSeedTrackIpSignificance", "TauJetSeedTrackIpSignificance", 100, 0., 100.);
    _hTauJetSeedTrackNhits = fs->make<TH1F>("TauJetSeedTrackNhits", "TauJetSeedTrackNhits", 40, 0., 40.);
    _hTauJetSeedTrackChi2  = fs->make<TH1F>("TauJetSeedTrackChi2", "TauJetSeedTrackChi2", 50, 0., 100.);
    _hTauJetCharge             = fs->make<TH1F>("TauJetCharge", "TauJetCharge", 5, 0., 5.);
    _hTauJetSignalTracksMass   = fs->make<TH1F>("TauJetSignalTracksMass", "TauJetSignalTracksMass", 50, 0., 5.);
    _hTauJetSignalTracksChargeFraction   = fs->make<TH1F>("TauJetSignalTracksChargeFraction", "TauJetSignalTracksChargeFraction", 30, 0., 1.5);
    _hTauJetNumIsoTracks       = fs->make<TH1F>("TauJetNumIsoTracks", "TauJetNumIsoTracks", 10, 0, 10);
    _hTauJetNumIsoGammas       = fs->make<TH1F>("TauJetNumIsoGammas", "TauJetNumIsoGammas", 10, 0, 10);
    _hTauJetNumIsoCands       = fs->make<TH1F>("TauJetNumIsoCands", "TauJetNumIsoCands", 10, 0, 10);
    _hTauJetSumPtIsoTracks = fs->make<TH1F>("TauJetSumPtIsoTracks", "TauJetSumPtIsoTracks", 100, 0, 50);
    _hTauJetSumPtIsoGammas = fs->make<TH1F>("TauJetSumPtIsoGammas", "TauJetSumPtIsoGammas", 100, 0, 50);
    _hTauJetSumPtIso = fs->make<TH1F>("TauJetSumPtIso", "TauJetSumPtIso", 100, 0, 50);
    _hTauJetGenTauDeltaPhi = fs->make<TH1F>("TauJetGenTauDeltaPhi", "TauJetGenTauDeltaPhi", 800, -0.2, 0.2);
    _hTauJetGenTauDeltaEta = fs->make<TH1F>("TauJetGenTauDeltaEta", "TauJetGenTauDeltaEta", 800, -0.2, 0.2);
    _hTauJetGenTauDeltaPt = fs->make<TH1F>("TauJetGenTauDeltaPt", "TauJetGenTauDeltaPt", 500, -5, 5);
  }

  //--- book reconstruction level histograms 
  if (_FillRecoMuonHists) {
    _hNMuon = fs->make<TH1F>("NMuon", "NMuon", 21, 0., 20.);
    _hMuonEnergy = fs->make<TH1F>("MuonEnergy", "MuonEnergy", 200, 0., 500.);
    _hMuonPt     = fs->make<TH1F>("MuonPt", "MuonPt", 200, 0., 500.);
    _hMuonEta    = fs->make<TH1F>("MuonEta", "MuonEta", 72, -3.6, +3.6);
    _hMuonPhi    = fs->make<TH1F>("MuonPhi", "MuonPhi", 36, -TMath::Pi(), +TMath::Pi());
    _hMuonTrackIso = fs->make<TH1F>("MuonTrackIso", "MuonTrackIso", 100, 0, 50);
    _hMuonEcalIso = fs->make<TH1F>("MuonEcalIso", "MuonEcalIso", 100, 0, 50);
    _hMuonIso = fs->make<TH1F>("MuonIso", "MuonIso", 100, 0, 50);
    _hMuonIp = fs->make<TH1F>("MuonIp", "MuonIp", 500, -1, +1);
    _hMuonIpSignificance        = fs->make<TH1F>("MuonIpSignificance", "MuonIpSignificance", 100, 0., 100.);
    _hMuonGenMuonDeltaPhi = fs->make<TH1F>("MuonGenMuonDeltaPhi", "MuonGenMuonDeltaPhi", 800, -0.2, 0.2);
    _hMuonGenMuonDeltaEta = fs->make<TH1F>("MuonGenMuonDeltaEta", "MuonGenMuonDeltaEta", 800, -0.2, 0.2);
    _hMuonGenMuonDeltaPt = fs->make<TH1F>("MuonGenMuonDeltaPt", "MuonGenMuonDeltaPt", 500, -5, 5);
  }
  
  if (_FillRecoElectronHists) {
    _hNElectron = fs->make<TH1F>("NElectron", "NElectron", 21, 0., 20.);
    _hElectronEnergy = fs->make<TH1F>("ElectronEnergy", "ElectronEnergy", 200, 0., 500.);
    _hElectronPt     = fs->make<TH1F>("ElectronPt", "ElectronPt", 200, 0., 500.);
    _hElectronEta    = fs->make<TH1F>("ElectronEta", "ElectronEta", 72, -3.6, +3.6);
    _hElectronPhi    = fs->make<TH1F>("ElectronPhi", "ElectronPhi", 36, -TMath::Pi(), +TMath::Pi());
    _hElectronTrackIso = fs->make<TH1F>("ElectronTrackIso", "ElectronTrackIso", 100, 0, 50);
    _hElectronEcalIso = fs->make<TH1F>("ElectronEcalIso", "ElectronEcalIso", 100, 0, 50);
    _hElectronIp = fs->make<TH1F>("ElectronIp", "ElectronIp", 500, -1, +1);
    _hElectronEoverP = fs->make<TH1F>("ElectronEoverP", "ElectronEoverP", 60, 0, +3);
    _hElectronHoverEm = fs->make<TH1F>("ElectronHoverEm", "ElectronHoverEm", 300, 0, +3);
    _hElectronClassification = fs->make<TH1F>("ElectronClassification", "ElectronClassification", 200, 0, 200);
    _hElectronGenElectronDeltaPhi = fs->make<TH1F>("ElectronGenElectronDeltaPhi", "ElectronGenElectronDeltaPhi", 800, -0.2, 0.2);
    _hElectronGenElectronDeltaEta = fs->make<TH1F>("ElectronGenElectronDeltaEta", "ElectronGenElectronDeltaEta", 800, -0.2, 0.2);
    _hElectronGenElectronDeltaPt = fs->make<TH1F>("ElectronGenElectronDeltaPt", "ElectronGenElectronDeltaPt", 500, -5, 5);
  }

  if (_FillRecoJetHists) {
    _hNJet = fs->make<TH1F>("NJet", "NJet", 21, 0., 20.);
    _hJetEnergy = fs->make<TH1F>("JetEnergy", "JetEnergy", 200, 0., 500.);
    _hJetPt = fs->make<TH1F>("JetPt", "JetPt", 200, 0., 500.);
    _hJetEta = fs->make<TH1F>("JetEta", "JetEta", 72, -3.6, +3.6);
    _hJetPhi = fs->make<TH1F>("JetPhi", "JetPhi", 36, -TMath::Pi(), +TMath::Pi());
    _hBJetDiscrByTrackCounting = fs->make<TH1F>("BJetDiscrByTrackCounting", "BJetDiscrByTrackCounting", 400, -20, 20);
    _hBJetDiscrBySimpleSecondaryV = fs->make<TH1F>("BJetDiscrBySimpleSecondaryV", "BJetDiscrBySimpleSecondaryV", 400, -20, 20);
    _hBJetDiscrByCombinedSecondaryV = fs->make<TH1F>("BJetDiscrByCombinedSecondaryV", "BJetDiscrByCombinedSecondaryV", 400, -20, 20);
  }

  if (_FillTopologyHists) {
    if( ((_AnalyzeMuonForLeg1) && (_AnalyzeTauForLeg2)) || ((_AnalyzeMuonForLeg2) && (_AnalyzeTauForLeg1)) ) {
      _hMuonPtVsTauPt = fs->make<TH2F>("MuonPtVsTauPt", "MuonPtVsTauPt", 100, 0, 500, 100, 0, 500);
      _hMuonTauDeltaR = fs->make<TH1F>("MuonTauDeltaR", "MuonTauDeltaR", 100, 0, 5.);
      _hMuonTauDeltaPtDivSumPt = fs->make<TH1F>("MuonTauDeltaPtDivSumPt", "MuonTauDeltaPtDivSumPt", 100, -5, 5.);
      _hMuonMetMt = fs->make<TH1F>("MuonMetMt", "MuonMetMt", 100, 0, 500);
      _hTauMetMt = fs->make<TH1F>("TauMetMt", "TauMetMt", 100, 0, 500);
      _hMuonTauOSLS = fs->make<TH1F>("MuonTauOSLS", "MuonTauOSLS", 20, -10, 10);
      _hMuonTauCosDphi = fs->make<TH1F>("MuonTauCosDphi", "MuonTauCosDphi", 220, -1.1, 1.1);
      _hMuonMetDeltaPhi = fs->make<TH1F>("MuonMetDeltaPhi", "MuonMetDeltaPhi", 72, 0, +TMath::Pi());
      _hTauMetDeltaPhi = fs->make<TH1F>("TauMetDeltaPhi", "TauMetDeltaPhi", 72, 0, +TMath::Pi());
      _hMuonMetDeltaPhiVsMuonTauCosDphi = fs->make<TH2F>("MuonMetDeltaPhiVsMuonTauCosDphi", "MuonMetDeltaPhiVsMuonTauCosDphi", 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
    }
    if( ((_AnalyzeElectronForLeg1) && (_AnalyzeTauForLeg2)) || ((_AnalyzeElectronForLeg2) && (_AnalyzeTauForLeg1)) ) {
      _hElectronPtVsTauPt = fs->make<TH2F>("ElectronPtVsTauPt", "ElectronPtVsTauPt", 100, 0, 500, 100, 0, 500);
      _hElectronTauDeltaR = fs->make<TH1F>("ElectronTauDeltaR", "ElectronTauDeltaR", 100, 0, 5.);
      _hElectronTauDeltaPtDivSumPt = fs->make<TH1F>("ElectronTauDeltaPtDivSumPt", "ElectronTauDeltaPtDivSumPt", 100, -5, 5.);
      _hElectronMetMt = fs->make<TH1F>("ElectronMetMt", "ElectronMetMt", 100, 0, 500);
      _hTauMetMt = fs->make<TH1F>("TauMetMt", "TauMetMt", 100, 0, 500);
      _hElectronTauOSLS = fs->make<TH1F>("ElectronTauOSLS", "ElectronTauOSLS", 20, -10, 10);
      _hElectronTauCosDphi = fs->make<TH1F>("ElectronTauCosDphi", "ElectronTauCosDphi", 220, -1.1, 1.1);
      _hElectronMetDeltaPhi = fs->make<TH1F>("ElectronMetDeltaPhi", "ElectronMetDeltaPhi", 72, 0, +TMath::Pi());
      _hTauMetDeltaPhi = fs->make<TH1F>("TauMetDeltaPhi", "TauMetDeltaPhi", 72, 0, +TMath::Pi());
      _hElectronMetDeltaPhiVsElectronTauCosDphi = fs->make<TH2F>("ElectronMetDeltaPhiVsElectronTauCosDphi", "ElectronMetDeltaPhiVsElectronTauCosDphi", 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
    }
    if( ((_AnalyzeMuonForLeg1) && (_AnalyzeElectronForLeg2)) || ((_AnalyzeMuonForLeg2) && (_AnalyzeElectronForLeg1)) ) {
      _hElectronPtVsMuonPt = fs->make<TH2F>("ElectronPtVsMuonPt", "ElectronPtVsMuonPt", 100, 0, 500, 100, 0, 500);
      _hElectronMuonDeltaR = fs->make<TH1F>("ElectronMuonDeltaR", "ElectronMuonDeltaR", 100, 0, 5.);
      _hElectronMuonDeltaPtDivSumPt = fs->make<TH1F>("ElectronMuonDeltaPtDivSumPt", "ElectronMuonDeltaPtDivSumPt", 100, -5, 5.);
      _hElectronMetMt = fs->make<TH1F>("ElectronMetMt", "ElectronMetMt", 100, 0, 500);
      _hMuonMetMt = fs->make<TH1F>("MuonMetMt", "MuonMetMt", 100, 0, 500);
      _hElectronMuonOSLS = fs->make<TH1F>("ElectronMuonOSLS", "ElectronMuonOSLS", 20, -10, 10);
      _hElectronMuonCosDphi = fs->make<TH1F>("ElectronMuonCosDphi", "ElectronMuonCosDphi", 220, -1.1, 1.1);
      _hElectronMetDeltaPhi = fs->make<TH1F>("ElectronMetDeltaPhi", "ElectronMetDeltaPhi", 72, 0, +TMath::Pi());
      _hMuonMetDeltaPhi = fs->make<TH1F>("MuonMetDeltaPhi", "MuonMetDeltaPhi", 72, 0, +TMath::Pi());
      _hElectronMetDeltaPhiVsElectronMuonCosDphi = fs->make<TH2F>("ElectronMetDeltaPhiVsElectronMuonCosDphi", "ElectronMetDeltaPhiVsElectronMuonCosDphi", 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
    }
    if( ((_AnalyzeTauForLeg1) && (_AnalyzeTauForLeg2)) ) {
      _hTau1PtVsTau2Pt = fs->make<TH2F>("Tau1PtVsTau2Pt", "Tau1PtVsTau2Pt", 100, 0, 500, 100, 0, 500);
      _hTau1Tau2DeltaR = fs->make<TH1F>("Tau1Tau2DeltaR", "Tau1Tau2DeltaR", 100, 0, 5.);
      _hTau1Tau2DeltaPtDivSumPt = fs->make<TH1F>("Tau1Tau2DeltaPtDivSumPt", "Tau1Tau2DeltaPtDivSumPt", 100, -5, 5.);
      _hTau1MetMt = fs->make<TH1F>("Tau1MetMt", "Tau1MetMt", 100, 0, 500);
      _hTau2MetMt = fs->make<TH1F>("Tau2MetMt", "Tau2MetMt", 100, 0, 500);
      _hTau1Tau2OSLS = fs->make<TH1F>("Tau1Tau2OSLS", "Tau1Tau2OSLS", 20, -10, 10);
      _hTau1Tau2CosDphi = fs->make<TH1F>("Tau1Tau2CosDphi", "Tau1Tau2CosDphi", 220, -1.1, 1.1);
      _hTau1MetDeltaPhi = fs->make<TH1F>("Tau1MetDeltaPhi", "Tau1MetDeltaPhi", 72, 0, +TMath::Pi());
      _hTau2MetDeltaPhi = fs->make<TH1F>("Tau2MetDeltaPhi", "Tau2MetDeltaPhi", 72, 0, +TMath::Pi());
      _hTau1MetDeltaPhiVsTau1Tau2CosDphi = fs->make<TH2F>("Tau1MetDeltaPhiVsTau1Tau2CosDphi", "Tau1MetDeltaPhiVsTau1Tau2CosDphi", 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
    }
    if( ((_AnalyzeMuonForLeg1) && (_AnalyzeMuonForLeg2)) ) {
      _hMuon1PtVsMuon2Pt = fs->make<TH2F>("Muon1PtVsMuon2Pt", "Muon1PtVsMuon2Pt", 100, 0, 500, 100, 0, 500);
      _hMuon1Muon2DeltaR = fs->make<TH1F>("Muon1Muon2DeltaR", "Muon1Muon2DeltaR", 100, 0, 5.);
      _hMuon1Muon2DeltaPtDivSumPt = fs->make<TH1F>("Muon1Muon2DeltaPtDivSumPt", "Muon1Muon2DeltaPtDivSumPt", 100, -5, 5.);
      _hMuon1MetMt = fs->make<TH1F>("Muon1MetMt", "Muon1MetMt", 100, 0, 500);
      _hMuon2MetMt = fs->make<TH1F>("Muon2MetMt", "Muon2MetMt", 100, 0, 500);
      _hMuon1Muon2OSLS = fs->make<TH1F>("Muon1Muon2OSLS", "Muon1Muon2OSLS", 20, -10, 10);
      _hMuon1Muon2CosDphi = fs->make<TH1F>("Muon1Muon2CosDphi", "Muon1Muon2CosDphi", 220, -1.1, 1.1);
      _hMuon1MetDeltaPhi = fs->make<TH1F>("Muon1MetDeltaPhi", "Muon1MetDeltaPhi", 72, 0, +TMath::Pi());
      _hMuon2MetDeltaPhi = fs->make<TH1F>("Muon2MetDeltaPhi", "Muon2MetDeltaPhi", 72, 0, +TMath::Pi());
      _hMuon1MetDeltaPhiVsMuon1Muon2CosDphi = fs->make<TH2F>("Muon1MetDeltaPhiVsMuon1Muon2CosDphi", "Muon1MetDeltaPhiVsMuon1Muon2CosDphi", 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
    }
    if( ((_AnalyzeElectronForLeg1) && (_AnalyzeElectronForLeg2)) ) {
      _hElectron1PtVsElectron2Pt = fs->make<TH2F>("Electron1PtVsElectron2Pt", "Electron1PtVsElectron2Pt", 100, 0, 500, 100, 0, 500);
      _hElectron1Electron2DeltaR = fs->make<TH1F>("Electron1Electron2DeltaR", "Electron1Electron2DeltaR", 100, 0, 5.);
      _hElectron1Electron2DeltaPtDivSumPt = fs->make<TH1F>("Electron1Electron2DeltaPtDivSumPt", "Electron1Electron2DeltaPtDivSumPt", 100, -5, 5.);
      _hElectron1MetMt = fs->make<TH1F>("Electron1MetMt", "Electron1MetMt", 100, 0, 500);
      _hElectron2MetMt = fs->make<TH1F>("Electron2MetMt", "Electron2MetMt", 100, 0, 500);
      _hElectron1Electron2OSLS = fs->make<TH1F>("Electron1Electron2OSLS", "Electron1Electron2OSLS", 20, -10, 10);
      _hElectron1Electron2CosDphi = fs->make<TH1F>("Electron1Electron2CosDphi", "Electron1Electron2CosDphi", 220, -1.1, 1.1);
      _hElectron1MetDeltaPhi = fs->make<TH1F>("Electron1MetDeltaPhi", "Electron1MetDeltaPhi", 72, 0, +TMath::Pi());
      _hElectron2MetDeltaPhi = fs->make<TH1F>("Electron2MetDeltaPhi", "Electron2MetDeltaPhi", 72, 0, +TMath::Pi());
      _hElectron1MetDeltaPhiVsElectron1Electron2CosDphi = fs->make<TH2F>("Electron1MetDeltaPhiVsElectron1Electron2CosDphi", "Electron1MetDeltaPhiVsElectron1Electron2CosDphi", 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
    }
/*
    _hTauJetSumPtIso_SeedOS = fs->make<TH1F>("TauJetSumPtIso_SeedOS", "TauJetSumPtIso_SeedOS", 100, 0, 50);
    _hTauJetSumPtIso_JetOS = fs->make<TH1F>("TauJetSumPtIso_JetOS", "TauJetSumPtIso_JetOS", 100, 0, 50);
    _hTauJetSumPtIso_SeedLS = fs->make<TH1F>("TauJetSumPtIso_SeedLS", "TauJetSumPtIso_SeedLS", 100, 0, 50);
    _hTauJetSumPtIso_JetLS = fs->make<TH1F>("TauJetSumPtIso_JetLS", "TauJetSumPtIso_JetLS", 100, 0, 50);
*/
    _hNotReconstructableMass = fs->make<TH1F>("NotReconstructableMass", "NotReconstructableMass", 150, 0, 1500);
    _hReconstructableMass = fs->make<TH1F>("ReconstructableMass", "ReconstructableMass", 150, 0, 1500);
    _hPZeta = fs->make<TH1F>("PZeta", "PZeta", 200, -100, 100);
    _hPZetaVis = fs->make<TH1F>("PZetaVis", "PZetaVis", 100, 0, 100);
    _hZeta2D = fs->make<TH2F>("Zeta2D", "Zeta2D", 100, 0, 100, 200, -100, 100);
    _hZeta1D = fs->make<TH1F>("Zeta1D", "Zeta1D", 50, -100, 100);
    _hMet  = fs->make<TH1F>("Met", "Met", 100, 0, 1000);
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void HiMassTauAnalysis::endJob() {
  printEfficiency();  

  if(_CalculatePdfSystematicUncertanties) {
    std::cout << "------------------------------PDFAnalysis------------------------------" << std::endl;
    std::cout << "PDF uncertainties will be determined for the following sets: ";
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "---------------------------------------------------------------------------\n";
    std::cout << "PDF weight systematics summary (Denominator)" << std::endl;
    std::cout << "---------------------------------------------------------------------------\n";
    std::cout << "Analyzed data (reference): " << _totalEvents << " [events]" << std::endl;
    if (_totalEvents==0) return;
    for (unsigned int i=0; i<pdfWeightTags_.size(); ++i) {
      unsigned int nmembers = weightedEvents_Denominator_.size()-pdfStart_Denominator_[i];
      if (i<pdfWeightTags_.size()-1) nmembers = pdfStart_Denominator_[i+1] - pdfStart_Denominator_[i];
      unsigned int npairs = (nmembers-1)/2;
      std::cout << "Results for PDF set " << pdfWeightTags_[i].label() << " ---->" << std::endl;
      double events_central = weightedEvents_Denominator_[pdfStart_Denominator_[i]]; 
      std::cout << "\tEstimate for central PDF member: " << int(events_central) << " [events]" << std::endl;
      std::cout << "\ti.e. " << std::setprecision(4) << 100*(events_central/_totalEvents-1.) << "% variation with respect to original PDF" << std::endl;
      if (npairs>0) {
        std::cout << "\tNumber of eigenvectors for uncertainty estimation: " << npairs << std::endl;
        double wplus = 0.; double wminus = 0.;
        for (unsigned int j=0; j<npairs; ++j) {
          double wa = weightedEvents_Denominator_[pdfStart_Denominator_[i]+2*j+1]/events_central-1.;
          double wb = weightedEvents_Denominator_[pdfStart_Denominator_[i]+2*j+2]/events_central-1.; 
          if (wa>wb) {
            if (wa<0.) wa = 0.;
            if (wb>0.) wb = 0.;
            wplus += wa*wa;
            wminus += wb*wb;
          } else {
            if (wb<0.) wb = 0.;
            if (wa>0.) wa = 0.;
            wplus += wb*wb;
            wminus += wa*wa;
          }
        }
        if (wplus>0) wplus = sqrt(wplus); if (wminus>0) wminus = sqrt(wminus);
        std::cout << "\tRelative uncertainty with respect to central member: +" << std::setprecision(4) << 100.*wplus << " / -" << std::setprecision(4) << 100.*wminus << " [%]" << std::endl;
      } else {std::cout << "\tNO eigenvectors for uncertainty estimation" << std::endl;}
    }
    std::cout << "End of PDF weight systematics summary (Denominator)" << std::endl;

    std::cout << "" << std::endl;
    std::cout << "---------------------------------------------------------------------------\n";
    std::cout << "PDF weight systematics summary (Numerator)" << std::endl;
    std::cout << "---------------------------------------------------------------------------\n";
    std::cout << "Analyzed data (reference): " << _totalEventsPassingCuts << " [events]" << std::endl;
    if (_totalEventsPassingCuts==0) return;
    for (unsigned int i=0; i<pdfWeightTags_.size(); ++i) {
      unsigned int nmembers = weightedEvents_Numerator_.size()-pdfStart_Numerator_[i];
      if (i<pdfWeightTags_.size()-1) nmembers = pdfStart_Numerator_[i+1] - pdfStart_Numerator_[i];
      unsigned int npairs = (nmembers-1)/2;
      std::cout << "Results for PDF set " << pdfWeightTags_[i].label() << " ---->" << std::endl;
      double events_central = weightedEvents_Numerator_[pdfStart_Numerator_[i]];
      std::cout << "\tEstimate for central PDF member: " << int(events_central) << " [events]" << std::endl;
      std::cout << "\ti.e. " << std::setprecision(4) << 100*(events_central/_totalEventsPassingCuts-1.) << "% variation with respect to original PDF" << std::endl;
      if (npairs>0) {
        std::cout << "\tNumber of eigenvectors for uncertainty estimation: " << npairs << std::endl;
        double wplus = 0.; double wminus = 0.;
        for (unsigned int j=0; j<npairs; ++j) {
          double wa = weightedEvents_Numerator_[pdfStart_Numerator_[i]+2*j+1]/events_central-1.;
          double wb = weightedEvents_Numerator_[pdfStart_Numerator_[i]+2*j+2]/events_central-1.;
          if (wa>wb) {
            if (wa<0.) wa = 0.;
            if (wb>0.) wb = 0.;
            wplus += wa*wa;
            wminus += wb*wb;
          } else {
            if (wb<0.) wb = 0.;
            if (wa>0.) wa = 0.;
            wplus += wb*wb;
            wminus += wa*wa;
          }
        }
        if (wplus>0) wplus = sqrt(wplus); if (wminus>0) wminus = sqrt(wminus);
        std::cout << "\tRelative uncertainty with respect to central member: +" << std::setprecision(4) << 100.*wplus << " / -" << std::setprecision(4) << 100.*wminus << " [%]" << std::endl;
      } else {std::cout << "\tNO eigenvectors for uncertainty estimation" << std::endl;}
    }
    std::cout << "End of PDF weight systematics summary (Numerator)" << std::endl;
  }

}

HiMassTauAnalysis::~HiMassTauAnalysis() { }

//define this as a plug-in
DEFINE_FWK_MODULE(HiMassTauAnalysis);

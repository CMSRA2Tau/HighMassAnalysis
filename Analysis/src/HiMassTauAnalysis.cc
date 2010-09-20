// Authors: Andres Florez, Alfredo Gurrola, Eduardo Luiggi, Chi Nhan Nguyen

#include "HighMassAnalysis/Analysis/interface/HiMassTauAnalysis.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

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
  _DoRecoTauDiscrBySignalTracksAndGammasMass = iConfig.getParameter<bool>("DoRecoTauDiscrBySignalTracksAndGammasMass");
  _RecoTauSignal3ProngAndGammasMassMinCutValue = iConfig.getParameter<double>("RecoTauSignal3ProngAndGammasMassMinCutValue");
  _RecoTauSignal3ProngAndGammasMassMaxCutValue = iConfig.getParameter<double>("RecoTauSignal3ProngAndGammasMassMaxCutValue");
  _RecoTauSignal1ProngAndGammasMassForPionMinCutValue = iConfig.getParameter<double>("RecoTauSignal1ProngAndGammasMassForPionMinCutValue");
  _RecoTauSignal1ProngAndGammasMassForPionMaxCutValue = iConfig.getParameter<double>("RecoTauSignal1ProngAndGammasMassForPionMaxCutValue");
  _RecoTauSignal1ProngAndGammasMassForKaonVetoMinCutValue = iConfig.getParameter<double>("RecoTauSignal1ProngAndGammasMassForKaonVetoMinCutValue");
  _RecoTauSignal1ProngAndGammasMassForKaonVetoMaxCutValue = iConfig.getParameter<double>("RecoTauSignal1ProngAndGammasMassForKaonVetoMaxCutValue");
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
  _SetTANC = iConfig.getParameter<bool>("SetTANC");

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
  _DoRecoMuonDiscrByPionVeto = iConfig.getParameter<bool>("DoRecoMuonDiscrByPionVeto");
  _RecoMuonCaloCompCoefficient = iConfig.getParameter<double>("RecoMuonCaloCompCoefficient");
  _RecoMuonSegmCompCoefficient = iConfig.getParameter<double>("RecoMuonSegmCompCoefficient");
  _RecoMuonAntiPionCut = iConfig.getParameter<double>("RecoMuonAntiPionCut");

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
  _DoRecoElectronDiscrByEcalDrivenSeed = iConfig.getParameter<bool>("DoRecoElectronDiscrByEcalDrivenSeed");
  _DoRecoElectronDiscrByTrackerDrivenSeed = iConfig.getParameter<bool>("DoRecoElectronDiscrByTrackerDrivenSeed");

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
  _DoTauDiscrByIsZeeCut = iConfig.getParameter<bool>("DoTauDiscrByIsZeeCut");

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
  _CalculateFSRSystematics = iConfig.getParameter<bool>("CalculateFSRSystematics");
  _CalculateISRGluonSystematics = iConfig.getParameter<bool>("CalculateISRGluonSystematics");
  _CalculateISRGammaSystematics = iConfig.getParameter<bool>("CalculateISRGammaSystematics");
  _SmearTheMuon = iConfig.getParameter<bool>("SmearTheMuon");
  _MuonPtScaleOffset = iConfig.getParameter<double>("MuonPtScaleOffset");
  _MuonPtSigmaOffset = iConfig.getParameter<double>("MuonPtSigmaOffset");
  _MuonEtaScaleOffset = iConfig.getParameter<double>("MuonEtaScaleOffset");
  _MuonEtaSigmaOffset = iConfig.getParameter<double>("MuonEtaSigmaOffset");
  _MuonPhiScaleOffset = iConfig.getParameter<double>("MuonPhiScaleOffset");
  _MuonPhiSigmaOffset = iConfig.getParameter<double>("MuonPhiSigmaOffset");
  _SmearTheElectron = iConfig.getParameter<bool>("SmearTheElectron");
  _ElectronPtScaleOffset = iConfig.getParameter<double>("ElectronPtScaleOffset");
  _ElectronPtSigmaOffset = iConfig.getParameter<double>("ElectronPtSigmaOffset");
  _ElectronEtaScaleOffset = iConfig.getParameter<double>("ElectronEtaScaleOffset");
  _ElectronEtaSigmaOffset = iConfig.getParameter<double>("ElectronEtaSigmaOffset");
  _ElectronPhiScaleOffset = iConfig.getParameter<double>("ElectronPhiScaleOffset");
  _ElectronPhiSigmaOffset = iConfig.getParameter<double>("ElectronPhiSigmaOffset");
  _SmearTheTau = iConfig.getParameter<bool>("SmearTheTau");
  _TauPtScaleOffset = iConfig.getParameter<double>("TauPtScaleOffset");
  _TauPtSigmaOffset = iConfig.getParameter<double>("TauPtSigmaOffset");
  _TauEtaScaleOffset = iConfig.getParameter<double>("TauEtaScaleOffset");
  _TauEtaSigmaOffset = iConfig.getParameter<double>("TauEtaSigmaOffset");
  _TauPhiScaleOffset = iConfig.getParameter<double>("TauPhiScaleOffset");
  _TauPhiSigmaOffset = iConfig.getParameter<double>("TauPhiSigmaOffset");
  _SmearTheJet = iConfig.getParameter<bool>("SmearTheJet");
  _JetEnergyScaleOffset = iConfig.getParameter<double>("JetEnergyScaleOffset");
  _SmearThePt = iConfig.getParameter<bool>("SmearThePt");
  _SmearTheEta = iConfig.getParameter<bool>("SmearTheEta");
  _SmearThePhi = iConfig.getParameter<bool>("SmearThePhi");

  if(_CalculateISRGluonSystematics) {
    // default pt bin edges
    std::vector<double> defPtEdges;
    defPtEdges.push_back(0.);
    defPtEdges.push_back(999999.);
    bosonPtBinEdges_ = iConfig.getUntrackedParameter<std::vector<double> > ("BosonPtBinEdges",defPtEdges);
    unsigned int ninputs_expected = bosonPtBinEdges_.size()-1;
    // default weights
    std::vector<double> defWeights;
    defWeights.push_back(1.);
    ptWeights_ = iConfig.getUntrackedParameter<std::vector<double> > ("PtWeights",defWeights);
    if (ptWeights_.size()==1 && ninputs_expected>1) {
      for (unsigned int i=1; i<ninputs_expected; i++){ ptWeights_.push_back(ptWeights_[0]);}
    }
  }

}

// ------------ method called once each job just before starting event loop  ------------
void  HiMassTauAnalysis::beginJob() {
  _totalEvents = 0;  
  _totalEventsPassingCuts = 0;  
  if(_CalculatePdfSystematicUncertanties) {InitializeInfoForPDFSystematicUncertaintites();}
  setMapSelectionAlgoIDs();
  initMapSelectionCounters();
  //  bookHistograms();  
  //if (_DoProduceNtuple) {
    //initializeVectors();
    //setupBranches();
  //}
}
// Set branches for the ntuple
void  HiMassTauAnalysis::setupBranches() {
   // For the taus
  _HMTTree = new TTree(_NtupleTreeName.c_str(), "HiMassDiTau Tree");
  _HMTTree->Branch("tauPdgId",&_tauPdgId);
  _HMTTree->Branch("tauMatched",&_tauMatched);
  _HMTTree->Branch("tauMotherPdgId",&_tauMotherPdgId);
  _HMTTree->Branch("tauE",&_tauE);
  _HMTTree->Branch("tauEt",&_tauEt);
  _HMTTree->Branch("tauPt",&_tauPt);
  _HMTTree->Branch("tauCharge",&_tauCharge);
  _HMTTree->Branch("tauEta",&_tauEta);
  _HMTTree->Branch("tauPhi",&_tauPhi);
  _HMTTree->Branch("tauNProngs",&_tauNProngs);
  _HMTTree->Branch("tauLTPt",&_tauLTPt);
  _HMTTree->Branch("tauLTChi2",&_tauLTChi2);
  _HMTTree->Branch("tauLTRecHitsSize",&_tauLTRecHitsSize);
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
  // For the electrons
  _HMTTree->Branch("eMatched",&_eMatched);
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
  _HMTTree->Branch("eUserEcalIso", &_eUserEcalIso);
  _HMTTree->Branch("eUserHcalIso", &_eUserHcalIso);
  _HMTTree->Branch("eUserTrkIso", &_eUserTrkIso);
  _HMTTree->Branch("eSCE1x5", &_eSCE1x5);
  _HMTTree->Branch("eSCE2x5", &_eSCE2x5);
  _HMTTree->Branch("eSCE5x5", &_eSCE5x5);
  _HMTTree->Branch("eIp", &_eIp);
  _HMTTree->Branch("eIpError", &_eIpError);
  _HMTTree->Branch("eIp_ctf", &_eIp_ctf);
  _HMTTree->Branch("eIpError_ctf",&_eIpError_ctf);
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

  //-------Jet Info
  _HMTTree->Branch("jetPt",  &_jetPt);
  _HMTTree->Branch("jetEt",  &_jetEt);
  _HMTTree->Branch("jetE",   &_jetE);
  _HMTTree->Branch("jetEta", &_jetEta); 
  _HMTTree->Branch("jetPhi", &_jetPhi);
  _HMTTree->Branch("jetEmFraction", &_jetEmFraction);
 
  //-----Gen resolution
  _HMTTree->Branch("egenPt",  &_egenPt);
  _HMTTree->Branch("egenE",   &_egenE);
  _HMTTree->Branch("egenPhi", &_egenPhi);
  _HMTTree->Branch("egenEta", &_egenEta);

  _HMTTree->Branch("taugenPt",  &_taugenPt);
  _HMTTree->Branch("taugenE",   &_taugenE);
  _HMTTree->Branch("taugenPhi", &_taugenPhi);
  _HMTTree->Branch("taugenEta", &_taugenEta);

  _HMTTree->Branch("isEE", &_isEEv);
  _HMTTree->Branch("isEB", &_isEBv);
  _HMTTree->Branch("isEBEEGap",&_isEBEEGap) ;  
  _HMTTree->Branch("isEBEtaGap",&_isEBEtaGap); 
  _HMTTree->Branch("isEBEtaGap",&_isEBPhiGap) ; 
  _HMTTree->Branch("isEBEtaGap",&_isEEDeeGap) ; 
  _HMTTree->Branch("isEBEtaGap",&_isEERingGap) ;

  _HMTTree->Branch("eEcalDriven",        &_eEcalDrivenv);
  _HMTTree->Branch("eTrkDriven",         &_eTrkDrivenv);
  _HMTTree->Branch("eMissingHits", &_eMissingHits);

  _HMTTree->Branch("PDF_weights",        &_PDF_weightsv);
  _HMTTree->Branch("ISR_gluon_weights",  &_ISR_gluon_weightsv);
  _HMTTree->Branch("ISR_gamma_weights",  &_ISR_gamma_weightsv);
  _HMTTree->Branch("FSR_weights",        &_FSR_weightsv);

  _HMTTree->Branch("bJetDiscrByTrackCounting",&_bJetDiscrByTrackCounting);
  _HMTTree->Branch("bJetDiscrBySimpleSecondaryV",&_bJetDiscrBySimpleSecondaryV);
  _HMTTree->Branch("bJetDiscrByCombinedSecondaryV",&_bJetDiscrByCombinedSecondaryV);

  _HMTTree->Branch("tauTancDiscOnePercent",&_tauTancDiscOnePercent);	      
  _HMTTree->Branch("tauTancDiscHalfPercent",&_tauTancDiscHalfPercent);	      
  _HMTTree->Branch("tauTancDiscQuarterPercent",&_tauTancDiscQuarterPercent);  
  _HMTTree->Branch("tauTancDiscTenthPercent",&_tauTancDiscTenthPercent);   

}
// Initialize the vectors for the ntuple
void HiMassTauAnalysis::initializeVectors(){

  _tauMatched              = NULL;
  _tauPdgId                = NULL;
  _tauMotherPdgId          = NULL;
  _tauE                    = NULL;
  _tauEt                   = NULL; 
  _tauPt                   = NULL; 
  _tauCharge               = NULL;
  _tauEta                  = NULL;
  _tauPhi                  = NULL;
  _tauNProngs              = NULL;
  _tauLTPt                 = NULL;
  _tauLTChi2               = NULL;
  _tauLTRecHitsSize        = NULL;
  _tauIsoTrackPtSum        = NULL;
  _tauIsoTrkPtSumDR1_0MinPt1_0    = NULL;
  _tauIsoTrkPtSumDR1_0MinPt0_5    = NULL;
  _tauIsoTrkPtSumDR0_75MinPt1_0   = NULL;
  _tauIsoTrkPtSumDR0_75MinPt0_5   = NULL;
  _tauIsoGammaEtSum               = NULL;
  _tauIsoGammaEtSumDR1_0MinPt1_5  = NULL;
  _tauIsoGammaEtSumDR1_0MinPt1_0  = NULL;
  _tauIsoGammaEtSumDR0_75MinPt1_5 = NULL;
  _tauIsoGammaEtSumDR0_75MinPt1_0 = NULL;
  _tauEmFraction                  = NULL;
  _tauHcalTotOverPLead            = NULL;
  _tauHcalMaxOverPLead            = NULL;
  _tauHcal3x3OverPLead            = NULL;
  _tauElectronPreId               = NULL;
  _tauModifiedEOverP              = NULL;
  _tauBremsRecoveryEOverPLead     = NULL;
  _tauDiscAgainstElec             = NULL;
  _tauLTCharge                    = NULL;
  _tauLTSignedIp                  = NULL;
  _tauIsInTheCraks                = NULL;
  
  _eventIsZee     = NULL;
  _zeeMass        = NULL;
  _zeePtAsymm     = NULL;
  _eMatched       = NULL;
  _ePdgId         = NULL;
  _eMotherPdgId   = NULL;
  _eE             = NULL;
  _eEt            = NULL;
  _ePt            = NULL;
  _eCharge        = NULL;
  _eEta           = NULL;
  _ePhi           = NULL;
  _eSigmaEtaEta   = NULL;
  _eSigmaIEtaIEta = NULL;
  _eEOverP        = NULL;
  _eHOverEm       = NULL;
  _eDeltaPhiIn    = NULL;
  _eDeltaEtaIn    = NULL;
  _eEcalIsoPat = NULL;
  _eHcalIsoPat = NULL;
  _eTrkIsoPat = NULL;
  _eIsoPat = NULL;
  _eUserEcalIso = NULL;
  _eUserTrkIso = NULL;
  _eUserHcalIso = NULL;  
  _eSCE1x5        = NULL;
  _eSCE2x5        = NULL;
  _eSCE5x5        = NULL;
  _eIp            = NULL;
  _eIpError       = NULL;
  _eIp_ctf        = NULL; 
  _eIpError_ctf   = NULL;
  _eClass         = NULL;
  _eMissingHits   = NULL;
  
  _mEt            = NULL;
  
  _eTauMass       = NULL;
  _eTauCosDPhi    = NULL;
  _eTauDelatR     = NULL;
  _eTauMetMass    = NULL;
  _eTauCollMass   = NULL;
  _eTauPZeta      = NULL;
  _eTauPZetaVis   = NULL;
    
  _diTauMass      = NULL;
  _diTauPt        = NULL;
  _diTauEt        = NULL;

  //------Jet info
  _jetPt  = NULL;
  _jetEt  = NULL;
  _jetE   = NULL;
  _jetPhi = NULL;
  _jetEta = NULL;
  _jetEmFraction = NULL;

  //-----Gen resolution
  _egenPt  = NULL;
  _egenE   = NULL;
  _egenEta = NULL;
  _egenPhi = NULL;

  _taugenPt  = NULL;
  _taugenE   = NULL;
  _taugenEta = NULL;
  _taugenPhi = NULL;
 
  _isEEv        = NULL;
  _isEBv        = NULL;
  _isEBEEGap    = NULL;  
  _isEBEtaGap   = NULL; 
  _isEBPhiGap   = NULL; 
  _isEEDeeGap   = NULL; 
  _isEERingGap  = NULL;

  _eEcalDrivenv        = NULL;
  _eTrkDrivenv         = NULL;
  _PDF_weightsv        = NULL;
  _ISR_gluon_weightsv  = NULL;
  _ISR_gamma_weightsv  = NULL;
  _FSR_weightsv        = NULL;

  _bJetDiscrByTrackCounting = NULL;	  
  _bJetDiscrBySimpleSecondaryV = NULL;   
  _bJetDiscrByCombinedSecondaryV = NULL; 

  _tauTancDiscOnePercent = NULL;
  _tauTancDiscHalfPercent = NULL;
  _tauTancDiscQuarterPercent = NULL;
  _tauTancDiscTenthPercent = NULL;


}
// clear the vectors vefore each event
void HiMassTauAnalysis::clearVectors(){

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
  _tauLTChi2->clear();
  _tauLTRecHitsSize->clear();
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
  _eUserEcalIso->clear();
  _eUserTrkIso->clear();
  _eUserHcalIso->clear();
 
  _eSCE1x5->clear();
  _eSCE2x5->clear();
  _eSCE5x5->clear();
  _eIp->clear();
  _eIpError->clear();
  _eIp_ctf->clear();
  _eIpError_ctf->clear();
  _eClass->clear();
  _eMissingHits->clear();
  
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

    //------Jet info
  _jetPt->clear();
  _jetEt->clear();
  _jetE->clear();
  _jetPhi->clear();
  _jetEta->clear();
  _jetEmFraction->clear();

  //-----Gen resolution
  _egenPt->clear();
  _egenE->clear();
  _egenEta->clear();
  _egenPhi->clear();

  _taugenPt->clear();
  _taugenE->clear();
  _taugenEta->clear();
  _taugenPhi->clear();

  _isEEv->clear();
  _isEBv->clear();
  _isEBEEGap->clear();
  _isEBEtaGap->clear();
  _isEBPhiGap->clear();
  _isEEDeeGap->clear();
  _isEERingGap->clear();

  _eEcalDrivenv->clear();
  _eTrkDrivenv ->clear();
  _PDF_weightsv->clear();
  _ISR_gluon_weightsv->clear();
  _ISR_gamma_weightsv->clear();
  _FSR_weightsv->clear();

  _tauTancDiscOnePercent->clear();
  _tauTancDiscHalfPercent->clear();
  _tauTancDiscQuarterPercent->clear();
  _tauTancDiscTenthPercent->clear();

  _bJetDiscrByTrackCounting->clear();	  
  _bJetDiscrBySimpleSecondaryV->clear();   
  _bJetDiscrByCombinedSecondaryV->clear(); 

}

// ------------ method called to for each event  ------------
void HiMassTauAnalysis::analyze(const Event& iEvent, const EventSetup& iSetup) {

  //------Number of events analyzed (denominator)
  _totalEvents++;

//  std::cout << "before pdf" << std::endl;

  //-----Get weights for the calculation of pdf systematic uncertainties for the denominator
  pdfWeightVector.clear();
  if(_CalculatePdfSystematicUncertanties) {
    for (unsigned int i=0; i<pdfWeightTags_.size(); ++i) {
//      std::cout << "pdf tag = " << pdfWeightTags_[i] << std::endl;
      edm::Handle<std::vector<double> > weightHandle;
      iEvent.getByLabel(pdfWeightTags_[i], weightHandle);
      std::vector<double> weights = (*weightHandle);
      pdfWeightVector = weights;
      unsigned int nmembers = weights.size();
      if (pdfStart_Denominator_[i]<0) { // if it's the first event
        pdfStart_Denominator_[i] = weightedEvents_Denominator_.size();
        for (unsigned int j=0; j<nmembers; ++j) {weightedEvents_Denominator_.push_back(0.);}
      }
      for (unsigned int j=0; j<nmembers; ++j) {
        weightedEvents_Denominator_[pdfStart_Denominator_[i]+j] += weights[j];
//        std::cout << "weight member " << j << " = " << weights[j] << std::endl;
      }
    }
  }else{
    pdfWeightVector.push_back(1);
  } 

  //------Grab the handle to the relevant collections
  getCollections(iEvent,iSetup);

  //------Calculate event weight for isr systematic uncertanties
  if(_CalculateISRGluonSystematics) {
    unsigned int gensize = _genParticles->size();
    // Set as default weight the asymptotic value at high pt (i.e. value of last bin)
    isrgluon_weight = ptWeights_[ptWeights_.size()-1];
    unsigned int nbins = bosonPtBinEdges_.size()-1;
    bool foundCorrectBoson = false;
    for(unsigned int i = 0; i<gensize; ++i) {
      if(foundCorrectBoson) continue;
      const reco::GenParticle& part = (*_genParticles)[i];
      int id = part.pdgId();
      if (id!=23 && abs(id)!=24 && abs(id)!=32) continue;
      int status = part.status();
      if (status!=3) continue;
      foundCorrectBoson = true;
      double ptofpart = part.pt();
      if (ptofpart>bosonPtBinEdges_[0] && ptofpart<bosonPtBinEdges_[nbins]) {
        bool foundCorrectBin = false;
        for (unsigned int j=1; j<=nbins; ++j) {
          if(foundCorrectBin) continue;
          if (ptofpart>bosonPtBinEdges_[j]) continue;
          isrgluon_weight = ptWeights_[j-1];
          foundCorrectBin = true;
        }
      }
    }
  } else {
    isrgluon_weight = 1;
  }

  if(_CalculateISRGammaSystematics) {
    isrgamma_weight = 1;
    // Find the boson at the hard scattering level
    const reco::GenParticle* boson = 0;
    int parton1Key = -1;
    int parton2Key = -1;
    unsigned int gensize = _genParticles->size();
    for (unsigned int i = 0; i<gensize; ++i) {
      const reco::GenParticle& part = (*_genParticles)[i];
      int status = abs(part.status());
      if (status!=3) continue;
      if (part.numberOfMothers()!=2) continue;
      int partId = abs(part.pdgId());
      if (status==3 && (partId==23||abs(partId)==24||abs(partId)==32)) {
	boson = &(*_genParticles)[i];
	parton1Key = part.motherRef(0).key();
	parton2Key = part.motherRef(1).key();
	break;
      }
    }
    
    // Consider only photons near the hard-scattering process
    const reco::GenParticle* photon = 0;
    if (boson) {
      for (unsigned int i = 0; i<gensize; ++i) {
	photon = 0;
	const reco::GenParticle& part = (*_genParticles)[i];
	int status = abs(part.status());
	if (status!=1) continue;
	int partId = abs(part.pdgId());
	if (partId!=22)  continue;
	if (part.numberOfMothers()!=1) continue;
	int keyM = part.motherRef(0).key();
	const reco::GenParticle* mother = &(*_genParticles)[keyM];
	if (mother->status()!=3) continue;
	int mId = mother->pdgId();
	if (abs(mId)>6 && mId!=2212) continue;
	for (unsigned int j=0; j<mother->numberOfDaughters(); ++j){ 
	  int keyD = mother->daughterRef(j).key();
	  if (keyD==parton1Key || keyD==parton2Key) {
	    photon = &part;
	    break;
	  }
	}
	if (photon) break;
      }  
    }
    
    if (boson && photon) {
      reco::Candidate::LorentzVector smom = boson->p4() + photon->p4();
      double s = smom.M2();
      double sqrts = smom.M();
      
      // Go to CM using the boost direction of the boson+photon system
      ROOT::Math::Boost cmboost(smom.BoostToCM());
      reco::Candidate::LorentzVector photonCM(cmboost(photon->p4()));
      double pcostheta = (  smom.x()*photonCM.x() + smom.y()*photonCM.y() + smom.z()*photonCM.z() ) / smom.P();
      
      // Determine kinematic invariants
      double t = - sqrts * (photonCM.t()-pcostheta);
      double MV = boson->mass();
      double u = MV*MV - s - t;
      isrgamma_weight = 1. - 2*t*u/(s*s+MV*MV*MV*MV);
      //printf(">>>>>>>>> s %f t %f u %f, MV %f, weight = %f\n", s, t, u, MV, (*weight));
    }
  } else {
    isrgamma_weight = 1;
  }

  //------Calculate event weight for fsr systematic uncertanties
  if(_CalculateFSRSystematics) {
    unsigned int gensize = _genParticles->size();
    fsr_weight = 1;
    for (unsigned int i = 0; i<gensize; ++i) {
      const reco::GenParticle& lepton = (*_genParticles)[i];
      if (lepton.status()!=3) continue;
      int leptonId = lepton.pdgId();
      if (abs(leptonId)!=11 && abs(leptonId)!=13 && abs(leptonId)!=15) continue;
      if (lepton.numberOfMothers()!=1) continue;
      const reco::Candidate * boson = lepton.mother();
      int bosonId = abs(boson->pdgId());
      if (bosonId!=23  && bosonId!=24 && bosonId!=32) continue;
      double bosonMass = boson->mass();
      double leptonMass = lepton.mass();
      double leptonEnergy = lepton.energy();
      double cosLeptonTheta = cos(lepton.theta());
      double sinLeptonTheta = sin(lepton.theta());
      double leptonPhi = lepton.phi();
      
      int trueKey = i;
      if (lepton.numberOfDaughters()==0) { continue; }
      else if (lepton.numberOfDaughters()==1) { 
	int otherleptonKey = lepton.daughterRef(0).key();
	const reco::GenParticle& otherlepton = (*_genParticles)[otherleptonKey];
	if (otherlepton.pdgId()!=leptonId) continue;
	if (otherlepton.numberOfDaughters()<=1) continue;
	trueKey = otherleptonKey;
      }
    
      const reco::GenParticle& trueLepton = (*_genParticles)[trueKey];
      unsigned int nDaughters = trueLepton.numberOfDaughters();

      for (unsigned int j = 0; j<nDaughters; ++j) {
	const reco::Candidate * photon = trueLepton.daughter(j);
	if (photon->pdgId()!=22) continue;
	double photonEnergy = photon->energy();
	double cosPhotonTheta = cos(photon->theta());
	double sinPhotonTheta = sin(photon->theta());
	double photonPhi = photon->phi();
	double costheta = sinLeptonTheta*sinPhotonTheta*cos(leptonPhi-photonPhi) + cosLeptonTheta*cosPhotonTheta;
	// Missing O(alpha) terms in soft-collinear approach
	// Only for W, from hep-ph/0303260
	if (bosonId==24) {
	  double betaLepton = sqrt(1-pow(leptonMass/leptonEnergy,2));
	  double delta = - 8*photonEnergy *(1-betaLepton*costheta) / pow(bosonMass,3) / (1-pow(leptonMass/bosonMass,2)) / (4-pow(leptonMass/bosonMass,2)) * leptonEnergy * (pow(leptonMass,2)/bosonMass+2*photonEnergy);
	  fsr_weight *= (1 + delta);
	}
	// Missing NLO QED orders in QED parton shower approach
	// Change coupling scale from 0 to kT to estimate this effect
	fsr_weight *= alphaRatio(photonEnergy*sqrt(1-pow(costheta,2)));
      }
    }
  } else {
    fsr_weight = 1;
  }

  // deltas for recalculation of MET (used when studying systematics)
  deltaForMEx = 0;
  deltaForMEy = 0;

  //-----Smearing momentum and position for systematic uncertanties
  smearedMuonMomentumVector.clear();
  if(((_AnalyzeMuonForLeg1) || (_AnalyzeMuonForLeg2))) {
    if(_SmearTheMuon) {
      for(pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end();++patMuon) {
        smearedMuonMomentumVector.push_back(SmearLightLepton(*patMuon));
        deltaForMEx = deltaForMEx + patMuon->px() - SmearLightLepton(*patMuon).px();
        deltaForMEy = deltaForMEy + patMuon->py() - SmearLightLepton(*patMuon).py();
      }
    } else {
      for(pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end();++patMuon) {
        smearedMuonMomentumVector.push_back(patMuon->p4());
      }
    }
  }
  smearedElectronMomentumVector.clear();
  if(((_AnalyzeElectronForLeg1) || (_AnalyzeElectronForLeg2))) {
    if(_SmearTheElectron) {
      for(pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end();++patElectron) {
        smearedElectronMomentumVector.push_back(SmearLightLepton(*patElectron));
        deltaForMEx = deltaForMEx + patElectron->px() - SmearLightLepton(*patElectron).px();
        deltaForMEy = deltaForMEy + patElectron->py() - SmearLightLepton(*patElectron).py();
      }
    } else {
      for(pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end();++patElectron) {
        smearedElectronMomentumVector.push_back(patElectron->p4());
      }
    }
  }
  smearedTauMomentumVector.clear();
  if(((_AnalyzeTauForLeg1) || (_AnalyzeTauForLeg2))) {
    if(_SmearTheTau) {
      for(pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end();++patTau) {
        smearedTauMomentumVector.push_back(SmearTau(*patTau));
        deltaForMEx = deltaForMEx + patTau->px() - SmearTau(*patTau).px();
        deltaForMEy = deltaForMEy + patTau->py() - SmearTau(*patTau).py();
      }
    } else {
      for(pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end();++patTau) {
        smearedTauMomentumVector.push_back(patTau->p4());
      }
    }
  }
  smearedJetMomentumVector.clear();
  if(_SmearTheJet) {
    for ( pat::JetCollection::const_iterator patJet = _patJets->begin(); patJet != _patJets->end(); ++patJet ) {
      smearedJetMomentumVector.push_back(SmearJet(*patJet));
      deltaForMEx = deltaForMEx + patJet->px() - SmearJet(*patJet).px();
      deltaForMEy = deltaForMEy + patJet->py() - SmearJet(*patJet).py();
    }
  } else {
    for ( pat::JetCollection::const_iterator patJet = _patJets->begin(); patJet != _patJets->end(); ++patJet ) {
      smearedJetMomentumVector.push_back(patJet->p4());
    }
  }

  if(_totalEvents == 1) {
    bookHistograms();
    if (_DoProduceNtuple) {
      initializeVectors();
      setupBranches();
    }
  }

  //------Number of events analyzed (denominator)
  for(unsigned int NpdfID = 0; NpdfID < pdfWeightVector.size();  NpdfID++){ 
    _hEvents[NpdfID]->Fill(0.0,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
  }

  //------Get the event flags (did the event pass the cuts?)
  getEventFlags(iEvent);
  if (_DoProduceNtuple){
    fillNtuple();
    clearVectors();
  }

  if (!passEventSelectionSequence()) return;  

  //------Number of events passing cuts (numerator)
  _totalEventsPassingCuts++;
  for(unsigned int NpdfID = 0; NpdfID < pdfWeightVector.size();  NpdfID++){   
    _hEvents[NpdfID]->Fill(1.0,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
  }

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

void HiMassTauAnalysis::getEventFlags(const Event& iEvent) {

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
  if(passRecoTriggerCuts(iEvent)) {nTriggersSatisfied++;}
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
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); 
	  patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if (!passRecoTauCuts((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1))) continue;
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
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); 
	  patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if (!passRecoTauCuts((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1))) continue;
      nGoodCandidatesLeg2++;
    }
  }
  if (nGoodCandidatesLeg2>=_RecoLeg2Nmin) _EventFlag[_mapSelectionAlgoID["RecoLeg2Nmin"]] = true;
  if (nGoodCandidatesLeg2<=_RecoLeg2Nmax) _EventFlag[_mapSelectionAlgoID["RecoLeg2Nmax"]] = true;

  // ------Number of Good Jets   
//  std::cout << "jet selections ..." << std::endl;
  int nGoodJets = 0;
  int theNumberOfJets = 0;
  for ( pat::JetCollection::const_iterator patJet = _patJets->begin(); 
	patJet != _patJets->end(); ++patJet ) {
    theNumberOfJets++;
    if (!passRecoJetCuts((*patJet),_SmearTheJet,smearedJetMomentumVector.at(theNumberOfJets-1))) continue;
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
      int theNumberOfTaus = 0;
      for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
        theNumberOfTaus++;
        if ((passRecoTauCuts((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1))) &&
            (passRecoMuonCuts((*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1))) && 
            (passTopologyCuts((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),
                              (*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),
                              (*(_patMETs->begin()))))) {
          nGoodCombinations++;
        }
      }
    }
  }
  if( ((_AnalyzeElectronForLeg1) && (_AnalyzeTauForLeg2)) || ((_AnalyzeElectronForLeg2) && (_AnalyzeTauForLeg1)) ) {
    int theNumberOfElectrons = 0;
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
      theNumberOfElectrons++;
      int theNumberOfTaus = 0;
      for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
        theNumberOfTaus++;
        if ((passRecoTauCuts((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1))) &&
            (passRecoElectronCuts((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1))) &&
            (passTopologyCuts((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),
                              (*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),
                              (*(_patMETs->begin()))))) {
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
    int theNumberOfTaus1 = 0;
    for ( pat::TauCollection::const_iterator patTau1 = _patTaus->begin();patTau1 != _patTaus->end(); ++patTau1 ) {
      theNumberOfTaus1++;
      int theNumberOfTaus2 = 0;
      for ( pat::TauCollection::const_iterator patTau2 = _patTaus->begin();patTau2 != _patTaus->end(); ++patTau2 ) {
        theNumberOfTaus2++;
        if ((passRecoTauCuts((*patTau1),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus1 - 1))) && 
            (passRecoTauCuts((*patTau2),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus2 - 1))) && 
            (passTopologyCuts((*patTau1),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus1 - 1),(*patTau2),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus2 - 1),(*(_patMETs->begin()))))) {
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

bool HiMassTauAnalysis::passRecoTriggerCuts(const Event& iEvent) {
  const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*_triggerResults);
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
bool HiMassTauAnalysis::passRecoTauCuts(const pat::Tau& patTau,bool smear,reco::Candidate::LorentzVector smearedLV) {
  // ----Matching to gen
  if(_MatchTauToGen) {
    if(_GenParticleSource.label() != "") {
      if(!(matchToGen(patTau).first)) {return false;}
    } else {return false;}
  }
  // ----Acceptance cuts
  if(smear) {
    if (fabs(smearedLV.eta())>_RecoTauEtaCut) {return false;}
    if (smearedLV.pt()<_RecoTauPtMinCut) {return false;}
    if (smearedLV.pt()>_RecoTauPtMaxCut) {return false;}
  } else {
    if (fabs(patTau.eta())>_RecoTauEtaCut) {return false;}
    if (patTau.pt()<_RecoTauPtMinCut) {return false;}
    if (patTau.pt()>_RecoTauPtMaxCut) {return false;}
  }
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
  // ----Signal constituents invariant mass cuts
  if (_DoRecoTauDiscrBySignalTracksAndGammasMass) {
    if (patTau.isCaloTau()) {
      if(patTau.signalTracks().size() == 1) {}
      if(patTau.signalTracks().size() == 3) {}
    } else {
      if(patTau.signalPFChargedHadrCands().size() == 1) {
        if(( (CalculateTauSignalTracksAndGammasMass(patTau).M() <= _RecoTauSignal1ProngAndGammasMassForPionMaxCutValue) && 
            (CalculateTauSignalTracksAndGammasMass(patTau).M() >= _RecoTauSignal1ProngAndGammasMassForPionMinCutValue) ) || 
          ( (CalculateTauSignalTracksAndGammasMass(patTau).M() <= _RecoTauSignal1ProngAndGammasMassForKaonVetoMaxCutValue) && 
            (CalculateTauSignalTracksAndGammasMass(patTau).M() >= _RecoTauSignal1ProngAndGammasMassForKaonVetoMinCutValue) )) { }
        else {return false;}
      }
      if(patTau.signalPFChargedHadrCands().size() == 3) {
        if(CalculateTauSignalTracksAndGammasMass(patTau).M() > _RecoTauSignal3ProngAndGammasMassMaxCutValue) {return false;}
        if(CalculateTauSignalTracksAndGammasMass(patTau).M() < _RecoTauSignal3ProngAndGammasMassMinCutValue) {return false;}
      }
    }
  }
  // ----Electron and Muon vetos
  if (_DoRecoTauDiscrAgainstElectron) { if ( (patTau.tauID(_RecoTauDiscrAgainstElectron.data()) < 0.5) ) {return false;} }
  if (_DoRecoTauDiscrByCrackCut) {
    if(smear) {
      if(isInTheCracks(smearedLV.eta())) {return false;}
    } else {
      if(isInTheCracks(patTau.eta())) {return false;}
    }
  }
  if (_DoRecoTauDiscrAgainstMuon) { if ( (patTau.tauID(_RecoTauDiscrAgainstMuon.data()) < 0.5) ) {return false;} }
  return true;
}

// ---------------Apply Muon Cuts
bool HiMassTauAnalysis::passRecoMuonCuts(const pat::Muon& patMuon,bool smear,reco::Candidate::LorentzVector smearedLV) {
  // ----Matching to gen
  if(_MatchLeptonToGen) {
    if(_GenParticleSource.label() != "") {
      if(!(matchToGen(patMuon).first)) {return false;}
    } else {return false;}
  }
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
  if (_DoRecoMuonDiscrByPionVeto) {
    if( ((_RecoMuonCaloCompCoefficient * muon::caloCompatibility(patMuon)) + (_RecoMuonSegmCompCoefficient * muon::segmentCompatibility(patMuon))) <= _RecoMuonAntiPionCut ) {return false;}
  }
  return true;
}

//--------------Apply Electron Cuts
bool HiMassTauAnalysis::passRecoElectronCuts(const pat::Electron& patElectron,bool smear,reco::Candidate::LorentzVector smearedLV) {
  // ----Matching to gen
  if(_MatchLeptonToGen) {
    if(_GenParticleSource.label() != "") {
      if(!(matchToGen(patElectron).first)) {return false;}
    } else {return false;}
  }
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
//  if (_DoRecoElectronDiscrByEcalDrivenSeed) { if(!(patElectron.ecalDrivenSeed())) {return false;} }
//  if (_DoRecoElectronDiscrByTrackerDrivenSeed) { if(!(patElectron.trackerDrivenSeed())) {return false;} }
  return true;
}

//--------------Apply Jet Cuts
bool HiMassTauAnalysis::passRecoJetCuts(const pat::Jet& patJet,bool smear,reco::Candidate::LorentzVector smearedLV) {
  // use corrected jet?
  if(_UseCorrectedJet) {
    // ----Acceptance cuts
    if(smear) {
      if (fabs(smearedLV.eta())>_RecoJetEtaMaxCut) {return false;}
      if (fabs(smearedLV.eta())<_RecoJetEtaMinCut) {return false;}
      if (smearedLV.pt()<_RecoJetPtCut) {return false;}
      // ----anti-overlap requirements
      if (_RemoveJetOverlapWithMuons) {
        int theNumberOfMuons = 0;
        for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
          theNumberOfMuons++;
          if( (passRecoMuonCuts((*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1))) ) {
            if(_SmearTheMuon) {
              if(reco::deltaR(smearedLV, smearedMuonMomentumVector.at(theNumberOfMuons-1)) < _JetMuonMatchingDeltaR) {return false;}
            } else {
              if(reco::deltaR(smearedLV, patMuon->p4()) < _JetMuonMatchingDeltaR) {return false;}
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
              if(reco::deltaR(smearedLV, smearedElectronMomentumVector.at(theNumberOfElectrons-1)) < _JetElectronMatchingDeltaR) {return false;}
            } else {
              if(reco::deltaR(smearedLV, patElectron->p4()) < _JetElectronMatchingDeltaR) {return false;}
            }
          }
        }
      }
      if (_RemoveJetOverlapWithTaus) {
        int theNumberOfTaus = 0;
        for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
          theNumberOfTaus++;
          if( (passRecoTauCuts((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1))) ) {
            if(_SmearTheTau) {
              if(reco::deltaR(smearedLV, smearedTauMomentumVector.at(theNumberOfTaus-1)) < _JetTauMatchingDeltaR) {return false;}
            } else {
              if(reco::deltaR(smearedLV, patTau->p4()) < _JetTauMatchingDeltaR) {return false;}
            }
          }
        }
      }
    } else {
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
        int theNumberOfTaus = 0;
        for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
          theNumberOfTaus++;
          if( (passRecoTauCuts((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1))) ) {
            if(_SmearTheTau) {
              if(reco::deltaR(patJet.p4(), smearedTauMomentumVector.at(theNumberOfTaus-1)) < _JetTauMatchingDeltaR) {return false;}
            } else {
              if(reco::deltaR(patJet.p4(), patTau->p4()) < _JetTauMatchingDeltaR) {return false;}
            }
          }
        }
      }
    }
  } else {
    // ----Acceptance cuts
    if(smear) {
      if (fabs(smearedLV.eta())>_RecoJetEtaMaxCut) {return false;}
      if (fabs(smearedLV.eta())<_RecoJetEtaMinCut) {return false;}
      if (smearedLV.pt()<_RecoJetPtCut) {return false;}
      // ----anti-overlap requirements
      if (_RemoveJetOverlapWithMuons) {
        int theNumberOfMuons = 0;
        for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
          theNumberOfMuons++;
          if( (passRecoMuonCuts((*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1))) ) {
            if(_SmearTheMuon) {
              if(reco::deltaR(smearedLV, smearedMuonMomentumVector.at(theNumberOfMuons-1)) < _JetMuonMatchingDeltaR) {return false;}
            } else {
              if(reco::deltaR(smearedLV, patMuon->p4()) < _JetMuonMatchingDeltaR) {return false;}
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
              if(reco::deltaR(smearedLV, smearedElectronMomentumVector.at(theNumberOfElectrons-1)) < _JetElectronMatchingDeltaR) {return false;}
            } else {
              if(reco::deltaR(smearedLV, patElectron->p4()) < _JetElectronMatchingDeltaR) {return false;}
            }
          }
        }
      }
      if (_RemoveJetOverlapWithTaus) {
        int theNumberOfTaus = 0;
        for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
          theNumberOfTaus++;
          if( (passRecoTauCuts((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1))) ) {
            if(_SmearTheTau) {
              if(reco::deltaR(smearedLV, smearedTauMomentumVector.at(theNumberOfTaus-1)) < _JetTauMatchingDeltaR) {return false;}
            } else {
              if(reco::deltaR(smearedLV, patTau->p4()) < _JetTauMatchingDeltaR) {return false;}
            }
          }
        }
      }
    } else {
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
        int theNumberOfTaus = 0;
        for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
          theNumberOfTaus++;
          if( (passRecoTauCuts((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1))) ) {
            if(_SmearTheTau) {
              if(reco::deltaR(patJet.correctedJet("raw","").p4(), smearedTauMomentumVector.at(theNumberOfTaus-1)) < _JetTauMatchingDeltaR) {return false;}
            } else {
              if(reco::deltaR(patJet.correctedJet("raw","").p4(), patTau->p4()) < _JetTauMatchingDeltaR) {return false;}
            }
          }
        }
      }
    }
  }
  return true;
}

// ---------------Apply Topology Cuts
bool HiMassTauAnalysis::passTopologyCuts(const pat::Tau& patTau, bool smearT, reco::Candidate::LorentzVector smearedLVT, const pat::Muon& patMuon, bool smearM, reco::Candidate::LorentzVector smearedLVM, const pat::MET& patMET) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if(smearT) {
    if(smearM) {
      if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(smearedLVM, smearedLVT) < _DiTauDeltaRCut) {return false;} }
    } else {
      if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(patMuon.p4(), smearedLVT) < _DiTauDeltaRCut) {return false;} }
    }
  } else {
    if(smearM) {
      if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(smearedLVM, patTau.p4()) < _DiTauDeltaRCut) {return false;} }
    } else {
      if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(patMuon.p4(), patTau.p4()) < _DiTauDeltaRCut) {return false;} }
    }
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
  if(smearT) {
    if(smearM) {
      if (_DoDiTauDiscrByCosDphi) {
        if(cos(TMath::Abs(normalizedPhi(smearedLVT.phi() - smearedLVM.phi()))) > _DiTauCosDphiMaxCut) {return false;}
        if(cos(TMath::Abs(normalizedPhi(smearedLVT.phi() - smearedLVM.phi()))) < _DiTauCosDphiMinCut) {return false;}
      }
    } else {
      if (_DoDiTauDiscrByCosDphi) {
        if(cos(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMuon.phi()))) > _DiTauCosDphiMaxCut) {return false;}
        if(cos(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMuon.phi()))) < _DiTauCosDphiMinCut) {return false;}
      }
    }
  } else {
    if(smearM) {
      if (_DoDiTauDiscrByCosDphi) {
        if(cos(TMath::Abs(normalizedPhi(patTau.phi() - smearedLVM.phi()))) > _DiTauCosDphiMaxCut) {return false;}
        if(cos(TMath::Abs(normalizedPhi(patTau.phi() - smearedLVM.phi()))) < _DiTauCosDphiMinCut) {return false;}
      }
    } else {
      if (_DoDiTauDiscrByCosDphi) {
        if(cos(TMath::Abs(normalizedPhi(patTau.phi() - patMuon.phi()))) > _DiTauCosDphiMaxCut) {return false;}
        if(cos(TMath::Abs(normalizedPhi(patTau.phi() - patMuon.phi()))) < _DiTauCosDphiMinCut) {return false;}
      }
    }
  }
  // ----Mass window requirement
  if (_DoDiscrByMassReco) { if( (CalculateThe4Momentum(patTau,smearT,smearedLVT,patMuon,smearM,smearedLVM,patMET).second.M() < _MassMinCut) || (CalculateThe4Momentum(patTau,smearT,smearedLVT,patMuon,smearM,smearedLVM,patMET).second.M() > _MassMaxCut) ) {return false;} }
  // ----Zeta requirement
  if (_DoDiTauDiscrByCDFzeta2D) {
    if( ((_PZetaCutCoefficient * CalculatePZeta(patTau,smearT,smearedLVT,patMuon,smearM,smearedLVM,patMET)) + 
         (_PZetaVisCutCoefficient * CalculatePZetaVis(patTau,smearT,smearedLVT,patMuon,smearM,smearedLVM,patMET))) < _CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {
    if((TMath::Sqrt((patMET.px() + deltaForMEx)*(patMET.px() + deltaForMEx) + (patMET.py() + deltaForMEy)*(patMET.py() + deltaForMEy))) < _RecoMetCut) {return false;}
    //    if (patMET.pt()<_RecoMetCut) {return false;}
  }
  if(smearT) {
    if(smearM) {
      if (_DoDiTauDiscrByDeltaPtDivSumPt) {
        if ( ((smearedLVT.pt() - smearedLVM.pt()) / (smearedLVT.pt() + smearedLVM.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
        if ( ((smearedLVT.pt() - smearedLVM.pt()) / (smearedLVT.pt() + smearedLVM.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
      }
      if (_DoDiscrByLeg1MetDphi) {
        if(_AnalyzeMuonForLeg1) {
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
        if(_AnalyzeTauForLeg1) {
          if(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
      }
      if (_DoDiscrByLeg2MetDphi) {
        if(_AnalyzeMuonForLeg2) {
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
        }
        if(_AnalyzeTauForLeg2) {
          if(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
        }
      }
    } else {
      if (_DoDiTauDiscrByDeltaPtDivSumPt) {
        if ( ((smearedLVT.pt() - patMuon.pt()) / (smearedLVT.pt() + patMuon.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
        if ( ((smearedLVT.pt() - patMuon.pt()) / (smearedLVT.pt() + patMuon.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
      }
      if (_DoDiscrByLeg1MetDphi) {
        if(_AnalyzeMuonForLeg1) {
          if(TMath::Abs(normalizedPhi(patMuon.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(patMuon.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
        if(_AnalyzeTauForLeg1) {
          if(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
      }
      if (_DoDiscrByLeg2MetDphi) {
        if(_AnalyzeMuonForLeg2) {
          if(TMath::Abs(normalizedPhi(patMuon.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(patMuon.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
        }
        if(_AnalyzeTauForLeg2) {
          if(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
        }
      }
    }
  } else {
    if(smearM) {
      if (_DoDiTauDiscrByDeltaPtDivSumPt) {
        if ( ((patTau.pt() - smearedLVM.pt()) / (patTau.pt() + smearedLVM.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
        if ( ((patTau.pt() - smearedLVM.pt()) / (patTau.pt() + smearedLVM.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
      }
      if (_DoDiscrByLeg1MetDphi) {
        if(_AnalyzeMuonForLeg1) {
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
        if(_AnalyzeTauForLeg1) {
          if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
      }
      if (_DoDiscrByLeg2MetDphi) {
        if(_AnalyzeMuonForLeg2) {
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
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
  }
  return true;
}
bool HiMassTauAnalysis::passTopologyCuts(const pat::Tau& patTau, bool smearT, reco::Candidate::LorentzVector smearedLVT, const pat::Electron& patElectron, bool smearM, reco::Candidate::LorentzVector smearedLVM, const pat::MET& patMET) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if(smearT) {
    if(smearM) {
      if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(smearedLVM, smearedLVT) < _DiTauDeltaRCut) {return false;} }
    } else {
      if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(patElectron.p4(), smearedLVT) < _DiTauDeltaRCut) {return false;} }
    }
  } else {
    if(smearM) {
      if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(smearedLVM, patTau.p4()) < _DiTauDeltaRCut) {return false;} }
    } else {
      if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(patElectron.p4(), patTau.p4()) < _DiTauDeltaRCut) {return false;} }
    }
  }

  // ---- apply Zee veto cut
  if(smearM) { if (_DoTauDiscrByIsZeeCut) { if(!(isZee(smearedLVM).first)) {return false;} }
  } else { if (_DoTauDiscrByIsZeeCut) { if(!(isZee(patElectron.p4()).first)) {return false;} } }

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
  if(smearT) {
    if(smearM) {
      if (_DoDiTauDiscrByCosDphi) {
        if(cos(TMath::Abs(normalizedPhi(smearedLVT.phi() - smearedLVM.phi()))) > _DiTauCosDphiMaxCut) {return false;}
        if(cos(TMath::Abs(normalizedPhi(smearedLVT.phi() - smearedLVM.phi()))) < _DiTauCosDphiMinCut) {return false;}
      }
    } else {
      if (_DoDiTauDiscrByCosDphi) {
        if(cos(TMath::Abs(normalizedPhi(smearedLVT.phi() - patElectron.phi()))) > _DiTauCosDphiMaxCut) {return false;}
        if(cos(TMath::Abs(normalizedPhi(smearedLVT.phi() - patElectron.phi()))) < _DiTauCosDphiMinCut) {return false;}
      }
    }
  } else {
    if(smearM) {
      if (_DoDiTauDiscrByCosDphi) {
        if(cos(TMath::Abs(normalizedPhi(patTau.phi() - smearedLVM.phi()))) > _DiTauCosDphiMaxCut) {return false;}
        if(cos(TMath::Abs(normalizedPhi(patTau.phi() - smearedLVM.phi()))) < _DiTauCosDphiMinCut) {return false;}
      }
    } else {
      if (_DoDiTauDiscrByCosDphi) {
        if(cos(TMath::Abs(normalizedPhi(patTau.phi() - patElectron.phi()))) > _DiTauCosDphiMaxCut) {return false;}
        if(cos(TMath::Abs(normalizedPhi(patTau.phi() - patElectron.phi()))) < _DiTauCosDphiMinCut) {return false;}
      }
    }
  }
  // ----Mass window requirement
  if (_DoDiscrByMassReco) { if( (CalculateThe4Momentum(patTau,smearT,smearedLVT,patElectron,smearM,smearedLVM,patMET).second.M() < _MassMinCut) || (CalculateThe4Momentum(patTau,smearT,smearedLVT,patElectron,smearM,smearedLVM,patMET).second.M() > _MassMaxCut) ) {return false;} }
  // ----Zeta requirement
  if (_DoDiTauDiscrByCDFzeta2D) {
    if( ((_PZetaCutCoefficient * CalculatePZeta(patTau,smearT,smearedLVT,patElectron,smearM,smearedLVM,patMET)) + 
         (_PZetaVisCutCoefficient * CalculatePZetaVis(patTau,smearT,smearedLVT,patElectron,smearM,smearedLVM,patMET))) < _CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {
    if((TMath::Sqrt((patMET.px() + deltaForMEx)*(patMET.px() + deltaForMEx) + (patMET.py() + deltaForMEy)*(patMET.py() + deltaForMEy))) < _RecoMetCut) {return false;}
    //    if (patMET.pt()<_RecoMetCut) {return false;}
  }
  if(smearT) {
    if(smearM) {
      if (_DoDiTauDiscrByDeltaPtDivSumPt) {
        if ( ((smearedLVT.pt() - smearedLVM.pt()) / (smearedLVT.pt() + smearedLVM.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
        if ( ((smearedLVT.pt() - smearedLVM.pt()) / (smearedLVT.pt() + smearedLVM.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
      }
      if (_DoDiscrByLeg1MetDphi) {
        if(_AnalyzeElectronForLeg1) {
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
        if(_AnalyzeTauForLeg1) {
          if(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
      }
      if (_DoDiscrByLeg2MetDphi) {
        if(_AnalyzeElectronForLeg2) {
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
        }
        if(_AnalyzeTauForLeg2) {
          if(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
        }
      }
    } else {
      if (_DoDiTauDiscrByDeltaPtDivSumPt) {
        if ( ((smearedLVT.pt() - patElectron.pt()) / (smearedLVT.pt() + patElectron.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
        if ( ((smearedLVT.pt() - patElectron.pt()) / (smearedLVT.pt() + patElectron.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
      }
      if (_DoDiscrByLeg1MetDphi) {
        if(_AnalyzeElectronForLeg1) {
          if(TMath::Abs(normalizedPhi(patElectron.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(patElectron.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
        if(_AnalyzeTauForLeg1) {
          if(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
      }
      if (_DoDiscrByLeg2MetDphi) {
        if(_AnalyzeElectronForLeg2) {
          if(TMath::Abs(normalizedPhi(patElectron.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(patElectron.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
        }
        if(_AnalyzeTauForLeg2) {
          if(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVT.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
        }
      }
    }
  } else {
    if(smearM) {
      if (_DoDiTauDiscrByDeltaPtDivSumPt) {
        if ( ((patTau.pt() - smearedLVM.pt()) / (patTau.pt() + smearedLVM.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
        if ( ((patTau.pt() - smearedLVM.pt()) / (patTau.pt() + smearedLVM.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
      }
      if (_DoDiscrByLeg1MetDphi) {
        if(_AnalyzeElectronForLeg1) {
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
        if(_AnalyzeTauForLeg1) {
          if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(patTau.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
        }
      }
      if (_DoDiscrByLeg2MetDphi) {
        if(_AnalyzeElectronForLeg2) {
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
          if(TMath::Abs(normalizedPhi(smearedLVM.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
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
  if (_DoDiscrByMet) {
    if(smearE) {
      if(smearM) {
        if (TMath::Sqrt(((patMET.px() + patElectron.px() + patMuon.px() - smearedLVM.px() - smearedLVE.px() + deltaForMEx) * 
                         (patMET.px() + patElectron.px() + patMuon.px() - smearedLVM.px() - smearedLVE.px() + deltaForMEx)) + 
                        ((patMET.py() + patElectron.py() + patMuon.py() - smearedLVM.py() - smearedLVE.py() + deltaForMEy) * 
                         (patMET.py() + patElectron.py() + patMuon.py() - smearedLVM.py() - smearedLVE.py() + deltaForMEy))) < _RecoMetCut) {return false;}
      } else {
        if (TMath::Sqrt(((patMET.px() + patElectron.px() + patMuon.px() - smearedLVM.px() - patMuon.px() + deltaForMEx) * 
                         (patMET.px() + patElectron.px() + patMuon.px() - smearedLVM.px() - patMuon.px() + deltaForMEx)) + 
                        ((patMET.py() + patElectron.py() + patMuon.py() - smearedLVM.py() - patMuon.py() + deltaForMEy) * 
                         (patMET.py() + patElectron.py() + patMuon.py() - smearedLVM.py() - patMuon.py() + deltaForMEy))) < _RecoMetCut) {return false;}
      }
    } else {
      if(smearM) {
        if (TMath::Sqrt(((patMET.px() + patElectron.px() + patMuon.px() - patElectron.px() - smearedLVM.px() + deltaForMEx) * 
                         (patMET.px() + patElectron.px() + patMuon.px() - patElectron.px() - smearedLVM.px() + deltaForMEx)) + 
                        ((patMET.py() + patElectron.py() + patMuon.py() - patElectron.py() - smearedLVM.py() + deltaForMEy) * 
                         (patMET.py() + patElectron.py() + patMuon.py() - patElectron.py() - smearedLVM.py() + deltaForMEy))) < _RecoMetCut) {return false;}
      } else {
        if((TMath::Sqrt((patMET.px() + deltaForMEx)*(patMET.px() + deltaForMEx) + (patMET.py() + deltaForMEy)*(patMET.py() + deltaForMEy))) < _RecoMetCut) {return false;}
//        if (patMET.pt()<_RecoMetCut) {return false;}
      }
    }
  }
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
  if (_DoDiscrByMet) {
    if(smear1) {
      if(smear2) {
        if (TMath::Sqrt(((patMET.px() + patMuon1.px() + patMuon2.px() - smearedLV1.px() - smearedLV2.px() + deltaForMEx) * 
                         (patMET.px() + patMuon1.px() + patMuon2.px() - smearedLV1.px() - smearedLV2.px() + deltaForMEx)) + 
                        ((patMET.py() + patMuon1.py() + patMuon2.py() - smearedLV1.py() - smearedLV2.py() + deltaForMEy) * 
                         (patMET.py() + patMuon1.py() + patMuon2.py() - smearedLV1.py() - smearedLV2.py() + deltaForMEy))) < _RecoMetCut) {return false;}
      } else {
        if (TMath::Sqrt(((patMET.px() + patMuon1.px() + patMuon2.px() - smearedLV1.px() - patMuon2.px() + deltaForMEx) * 
                         (patMET.px() + patMuon1.px() + patMuon2.px() - smearedLV1.px() - patMuon2.px() + deltaForMEx)) + 
                        ((patMET.py() + patMuon1.py() + patMuon2.py() - smearedLV1.py() - patMuon2.py() + deltaForMEy) * 
                         (patMET.py() + patMuon1.py() + patMuon2.py() - smearedLV1.py() - patMuon2.py() + deltaForMEy))) < _RecoMetCut) {return false;}
      }
    } else {
      if(smear2) {
        if (TMath::Sqrt(((patMET.px() + patMuon1.px() + patMuon2.px() - patMuon1.px() - smearedLV2.px() + deltaForMEx) * 
                         (patMET.px() + patMuon1.px() + patMuon2.px() - patMuon1.px() - smearedLV2.px() + deltaForMEx)) + 
                        ((patMET.py() + patMuon1.py() + patMuon2.py() - patMuon1.py() - smearedLV2.py() + deltaForMEy) * 
                         (patMET.py() + patMuon1.py() + patMuon2.py() - patMuon1.py() - smearedLV2.py() + deltaForMEy))) < _RecoMetCut) {return false;}
      } else {
        if((TMath::Sqrt((patMET.px() + deltaForMEx)*(patMET.px() + deltaForMEx) + (patMET.py() + deltaForMEy)*(patMET.py() + deltaForMEy))) < _RecoMetCut) {return false;}
//        if (patMET.pt()<_RecoMetCut) {return false;}
      }
    }
  }
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
  if (_DoDiscrByMet) {
    if(smear1) {
      if(smear2) {
        if (TMath::Sqrt(((patMET.px() + patElectron1.px() + patElectron2.px() - smearedLV1.px() - smearedLV2.px() + deltaForMEx) * 
                         (patMET.px() + patElectron1.px() + patElectron2.px() - smearedLV1.px() - smearedLV2.px() + deltaForMEx)) + 
                        ((patMET.py() + patElectron1.py() + patElectron2.py() - smearedLV1.py() - smearedLV2.py() + deltaForMEy) * 
                         (patMET.py() + patElectron1.py() + patElectron2.py() - smearedLV1.py() - smearedLV2.py() + deltaForMEy))) < _RecoMetCut) {return false;}
      } else {
        if (TMath::Sqrt(((patMET.px() + patElectron1.px() + patElectron2.px() - smearedLV1.px() - patElectron2.px() + deltaForMEx) * 
                         (patMET.px() + patElectron1.px() + patElectron2.px() - smearedLV1.px() - patElectron2.px() + deltaForMEx)) + 
                        ((patMET.py() + patElectron1.py() + patElectron2.py() - smearedLV1.py() - patElectron2.py() + deltaForMEy) * 
                         (patMET.py() + patElectron1.py() + patElectron2.py() - smearedLV1.py() - patElectron2.py() + deltaForMEy))) < _RecoMetCut) {return false;}
      }
    } else {
      if(smear2) {
        if (TMath::Sqrt(((patMET.px() + patElectron1.px() + patElectron2.px() - patElectron1.px() - smearedLV2.px() + deltaForMEx) * 
                         (patMET.px() + patElectron1.px() + patElectron2.px() - patElectron1.px() - smearedLV2.px() + deltaForMEx)) + 
                        ((patMET.py() + patElectron1.py() + patElectron2.py() - patElectron1.py() - smearedLV2.py() + deltaForMEy) * 
                         (patMET.py() + patElectron1.py() + patElectron2.py() - patElectron1.py() - smearedLV2.py() + deltaForMEy))) < _RecoMetCut) {return false;}
      } else {
        if((TMath::Sqrt((patMET.px() + deltaForMEx)*(patMET.px() + deltaForMEx) + (patMET.py() + deltaForMEy)*(patMET.py() + deltaForMEy))) < _RecoMetCut) {return false;}
//        if (patMET.pt()<_RecoMetCut) {return false;}
      }
    }
  }
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
bool HiMassTauAnalysis::passTopologyCuts(const pat::Tau& patTau1, bool smear1, reco::Candidate::LorentzVector smearedLV1, const pat::Tau& patTau2, bool smear2, reco::Candidate::LorentzVector smearedLV2, const pat::MET& patMET) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoDiTauDiscrByDeltaR) {
    if(smear1) {
      if(reco::deltaR(smearedLV1, smearedLV2) < _DiTauDeltaRCut) {return false;}
    } else {
      if(reco::deltaR(patTau1.p4(), patTau2.p4()) < _DiTauDeltaRCut) {return false;}
    }
  }
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
    if(smear1) {
      if(cos(TMath::Abs(normalizedPhi(smearedLV1.phi() - smearedLV2.phi()))) > _DiTauCosDphiMaxCut) {return false;}
      if(cos(TMath::Abs(normalizedPhi(smearedLV1.phi() - smearedLV2.phi()))) < _DiTauCosDphiMinCut) {return false;}
    } else {
      if(cos(TMath::Abs(normalizedPhi(patTau1.phi() - patTau2.phi()))) > _DiTauCosDphiMaxCut) {return false;}
      if(cos(TMath::Abs(normalizedPhi(patTau1.phi() - patTau2.phi()))) < _DiTauCosDphiMinCut) {return false;}
    }
  }
  // ----Mass window requirement
  if (_DoDiscrByMassReco) { if( (CalculateThe4Momentum(patTau1,smear1,smearedLV1,patTau2,smear2,smearedLV2,patMET).second.M() < _MassMinCut) || (CalculateThe4Momentum(patTau1,smear1,smearedLV1,patTau2,smear2,smearedLV2,patMET).second.M() > _MassMaxCut) ) {return false;} }
  // ----Zeta requirement
  if (_DoDiTauDiscrByCDFzeta2D) {
    if( ((_PZetaCutCoefficient * CalculatePZeta(patTau1,smear1,smearedLV1,patTau2,smear2,smearedLV2,patMET)) +
         (_PZetaVisCutCoefficient * CalculatePZetaVis(patTau1,smear1,smearedLV1,patTau2,smear2,smearedLV2,patMET))) < _CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {
    if(smear1) {
      if(smear2) {
        if (TMath::Sqrt(((patMET.px() + patTau1.px() + patTau2.px() - smearedLV1.px() - smearedLV2.px() + deltaForMEx) * 
                         (patMET.px() + patTau1.px() + patTau2.px() - smearedLV1.px() - smearedLV2.px() + deltaForMEx)) + 
                        ((patMET.py() + patTau1.py() + patTau2.py() - smearedLV1.py() - smearedLV2.py() + deltaForMEy) * 
                         (patMET.py() + patTau1.py() + patTau2.py() - smearedLV1.py() - smearedLV2.py() + deltaForMEy))) < _RecoMetCut) {return false;}
      } else {
        if (TMath::Sqrt(((patMET.px() + patTau1.px() + patTau2.px() - smearedLV1.px() - patTau2.px() + deltaForMEx) * 
                         (patMET.px() + patTau1.px() + patTau2.px() - smearedLV1.px() - patTau2.px() + deltaForMEx)) + 
                        ((patMET.py() + patTau1.py() + patTau2.py() - smearedLV1.py() - patTau2.py() + deltaForMEy) * 
                         (patMET.py() + patTau1.py() + patTau2.py() - smearedLV1.py() - patTau2.py() + deltaForMEy))) < _RecoMetCut) {return false;}
      }
    } else {
      if(smear2) {
        if (TMath::Sqrt(((patMET.px() + patTau1.px() + patTau2.px() - patTau1.px() - smearedLV2.px() + deltaForMEx) * 
                         (patMET.px() + patTau1.px() + patTau2.px() - patTau1.px() - smearedLV2.px() + deltaForMEx)) + 
                        ((patMET.py() + patTau1.py() + patTau2.py() - patTau1.py() - smearedLV2.py() + deltaForMEy) * 
                         (patMET.py() + patTau1.py() + patTau2.py() - patTau1.py() - smearedLV2.py() + deltaForMEy))) < _RecoMetCut) {return false;}
      } else {
        if((TMath::Sqrt((patMET.px() + deltaForMEx)*(patMET.px() + deltaForMEx) + (patMET.py() + deltaForMEy)*(patMET.py() + deltaForMEy))) < _RecoMetCut) {return false;}
//        if (patMET.pt()<_RecoMetCut) {return false;}
      }
    }
  }
  if (_DoDiTauDiscrByDeltaPtDivSumPt) {
    if(smear1) {
      if ( ((smearedLV1.pt() - smearedLV2.pt()) / (smearedLV1.pt() + smearedLV2.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
      if ( ((smearedLV1.pt() - smearedLV2.pt()) / (smearedLV1.pt() + smearedLV2.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
    } else {
      if ( ((patTau1.pt() - patTau2.pt()) / (patTau1.pt() + patTau2.pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
      if ( ((patTau1.pt() - patTau2.pt()) / (patTau1.pt() + patTau2.pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
    }
  }
  if (_DoDiscrByLeg1MetDphi) {
    if(_AnalyzeTauForLeg1) {
      if(smear1) {
        if(TMath::Abs(normalizedPhi(smearedLV1.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(smearedLV1.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
      } else {
        if(TMath::Abs(normalizedPhi(patTau1.phi() - patMET.phi())) > _Leg1MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patTau1.phi() - patMET.phi())) < _Leg1MetDphiMinCut) {return false;}
      }
    }
  }
  if (_DoDiscrByLeg2MetDphi) {
    if(_AnalyzeTauForLeg2) {
      if(smear1) {
        if(TMath::Abs(normalizedPhi(smearedLV2.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(smearedLV2.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
      } else {
        if(TMath::Abs(normalizedPhi(patTau2.phi() - patMET.phi())) > _Leg2MetDphiMaxCut) {return false;}
        if(TMath::Abs(normalizedPhi(patTau2.phi() - patMET.phi())) < _Leg2MetDphiMinCut) {return false;}
      }
    }
  }
  return true;
}

// ---------------Fill Ntuple
void HiMassTauAnalysis::fillNtuple() {

  // store vector of PDF weights
  _PDF_weightsv->push_back(pdfWeightVector);

  // store ISR weights
  _ISR_gluon_weightsv->push_back(isrgluon_weight);
  _ISR_gamma_weightsv->push_back(isrgamma_weight);

  // store FSR weights
  _FSR_weightsv->push_back(fsr_weight);

  // fill jet information used for jet based cuts (e.g. jet veto)
  int theNumberOfJets = 0;
  for ( pat::JetCollection::const_iterator patJet = _patJets->begin(); patJet != _patJets->end(); ++patJet ) {
    // use corrected jets or raw calorimeter jets
    _bJetDiscrByTrackCounting->push_back(patJet->bDiscriminator("trackCountingHighEffBJetTags"));	     
    _bJetDiscrBySimpleSecondaryV->push_back(patJet->bDiscriminator("simpleSecondaryVertexBJetTags"));	   
    _bJetDiscrByCombinedSecondaryV->push_back(patJet->bDiscriminator("combinedSecondaryVertexBJetTags")); 
    if(_UseCorrectedJet) {
      // smear the jets using EDAnalyzer function "SmearJet"
      if(_SmearTheJet) {
        _jetPt->push_back(smearedJetMomentumVector.at(theNumberOfJets).pt());
        _jetEt->push_back(smearedJetMomentumVector.at(theNumberOfJets).energy() * sin(smearedJetMomentumVector.at(theNumberOfJets).theta()));
        _jetE->push_back(smearedJetMomentumVector.at(theNumberOfJets).energy());
        _jetEta->push_back(smearedJetMomentumVector.at(theNumberOfJets).eta()); 
        _jetPhi->push_back(smearedJetMomentumVector.at(theNumberOfJets).phi());
        _jetEmFraction->push_back(patJet->correctedJet("raw", "").emEnergyFraction());
      } else {
        _jetPt->push_back(patJet->pt());
        _jetEt->push_back(patJet->et());
        _jetE->push_back(patJet->energy());
        _jetEta->push_back(patJet->eta()); 
        _jetPhi->push_back(patJet->phi());
        _jetEmFraction->push_back(patJet->correctedJet("raw", "").emEnergyFraction());
      }
    } else {
      // smear the jets using EDAnalyzer function "SmearJet"
      if(_SmearTheJet) {
        _jetPt->push_back(smearedJetMomentumVector.at(theNumberOfJets).pt());
        _jetEt->push_back(smearedJetMomentumVector.at(theNumberOfJets).energy() * sin(smearedJetMomentumVector.at(theNumberOfJets).theta()));
        _jetE->push_back(smearedJetMomentumVector.at(theNumberOfJets).energy());
        _jetEta->push_back(smearedJetMomentumVector.at(theNumberOfJets).eta()); 
        _jetPhi->push_back(smearedJetMomentumVector.at(theNumberOfJets).phi());
        _jetEmFraction->push_back(patJet->correctedJet("raw", "").emEnergyFraction());
      } else {
        _jetPt->push_back((float)patJet->correctedJet("raw","").pt());
        _jetEt->push_back(patJet->correctedJet("raw","").et());
        _jetE->push_back(patJet->correctedJet("raw","").energy());
        _jetEta->push_back(patJet->correctedJet("raw","").eta()); 
        _jetPhi->push_back(patJet->correctedJet("raw","").phi());
        _jetEmFraction->push_back(patJet->correctedJet("raw", "").emEnergyFraction());
      }
    }
    theNumberOfJets++;
  }

  // loop over electrons
  int theNumberOfElectrons = 0;
  for(pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron) {
    theNumberOfElectrons++;
    // loop over taus to create e-tau pairs
    int theNumberOfTaus = 0;
    for(pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau) {
      theNumberOfTaus++;
      // determine whether ntuple information will be filled based on systematics (smearing of resolution and scale for taus or electrons)      
      if(_SmearTheElectron){
        // only apply electron based smearing if the patElectron is matched to a "true" generator level electron
        if(matchToGen(*patElectron).first) {
          reco::Candidate::LorentzVector genElectron = matchToGen(*patElectron).second;
          double genElectronPt;
          double genElectronEta;
          double genElectronPhi;
          // determine whether we want to smear the resolution or change the scale of the pt, eta, or phi of the electron distribution
          if(_SmearThePt){
            genElectronPt  = ((genElectron.pt())*_ElectronPtScaleOffset + (patElectron->pt() - genElectron.pt())*_ElectronPtSigmaOffset);
          }else{genElectronPt  = patElectron->pt();}
          if(_SmearTheEta){
            genElectronEta = ((genElectron.eta())*_ElectronEtaScaleOffset + (patElectron->eta()- genElectron.eta())*_ElectronEtaSigmaOffset);
          }else{genElectronEta = patElectron->eta();}
          if(_SmearThePhi){
            genElectronPhi = ((genElectron.phi())*_ElectronPhiScaleOffset + (patElectron->phi()- genElectron.phi())*_ElectronPhiSigmaOffset);
          }else{genElectronPhi = patElectron->phi();}

          // create pt,eta,phi smeared root vector in order to recalculate the lorentz 4-momentum
          math::PtEtaPhiMLorentzVector genElectronLorentzVector(genElectronPt, genElectronEta, genElectronPhi, genElectron.mass());
          reco::Candidate::LorentzVector LorentzGenElectron(genElectronLorentzVector.px(), genElectronLorentzVector.py(), genElectronLorentzVector.pz(), genElectronLorentzVector.energy());
          // store the recalculated 4-momentum (will be used below)	  
          Electron = LorentzGenElectron;
        }
      }else{
        // if smearing is set to false, use default patElectron 4-momentum
	Electron = patElectron->p4();
      }

      // determine whether ntuple information will be filled based on systematics (smearing of resolution and scale for taus or electrons) 
      if(_SmearTheTau){
        // only apply tau based smearing if the patTau is matched to a "true" generator level tau
        if(matchToGen(*patTau).first) {
          reco::Candidate::LorentzVector genTau = matchToGen(*patTau).second;
          double genTauPt;
          double genTauEta;
          double genTauPhi;

          // determine whether we want to smear the resolution or change the scale of the pt, eta, or phi of the tau distribution
          if(_SmearThePt){
            genTauPt  = ((genTau.pt())*_TauPtScaleOffset + (patTau->pt() - genTau.pt())*_TauPtSigmaOffset);
          }else{genTauPt  = patTau->pt();}
          if(_SmearTheEta){
            genTauEta = ((genTau.eta())*_TauEtaScaleOffset + (patTau->eta()- genTau.eta())*_TauEtaSigmaOffset);
          }else{genTauEta = patTau->eta();}
          if(_SmearThePhi){
            genTauPhi = ((genTau.phi())*_TauPhiScaleOffset + (patTau->phi()- genTau.phi())*_TauPhiSigmaOffset);
          }else{genTauPhi = patTau->phi();}

          // create pt,eta,phi smeared root vector in order to recalculate the lorentz 4-momentum
          math::PtEtaPhiMLorentzVector genTauLorentzVector(genTauPt, genTauEta, genTauPhi, genTau.mass());
          reco::Candidate::LorentzVector LorentzGenTau(genTauLorentzVector.px(), genTauLorentzVector.py(), genTauLorentzVector.pz(), genTauLorentzVector.energy());

          // store the recalculated 4-momentum (will be used below)	  	  
          Tau = LorentzGenTau;
        }
      }else{
        // if smearing is set to false, use default patElectron 4-momentum
	Tau = patTau->p4();
      }

      // Recalculate MET for cases when electron or taus are smeared (using EDAnalyzer method and configuration parameters to recalculate the MET)
      // Recalculating MET for cases when electron is smeared ----> set SmearTheTau to False, SmearTheElectron to True, and specify correct smearing factors
      // Recalculating MET for cases when tau is smeared ----> set SmearTheTau to True, SmearTheElectron to False, and specify correct smearing factors
      // Recalculating MET for cases when jets are smeared ----> set SmearTheTau to False, SmearTheElectron to False, SmearTheJet to True, and specify correct smearing factors
      // SMEARING FACTORS: smear ONLY electron pt, tau pt, and jet energy scale factor

      double met_e = TMath::Sqrt(((*(_patMETs->begin())).px() + deltaForMEx)*((*(_patMETs->begin())).px() + deltaForMEx) + ((*(_patMETs->begin())).py() + deltaForMEy)*((*(_patMETs->begin())).py() + deltaForMEy));
      double met_px = ((*(_patMETs->begin())).px() + deltaForMEx);
      double met_py = ((*(_patMETs->begin())).py() + deltaForMEy);
      double met_pz = (*(_patMETs->begin())).pz();

      reco::Candidate::LorentzVector LorentzVectorMET(met_px,  met_py,  met_pz, met_e);
      MET = LorentzVectorMET;

      //----eEcalDriven? eTrkDriven?
//      const bool eEcalDriven = patElectron->ecalDrivenSeed();
//      const bool eTrkDriven  = patElectron->trackerDrivenSeed();
      
//      _eEcalDrivenv->push_back(eEcalDriven);
//      _eTrkDrivenv->push_back(eTrkDriven);
      
      // Get generator level information for pt,eta,phi and energy, matching the reco and gen objects using EDAnalyzer matching function 
      _taugenPt->push_back(matchToGen(*patTau).second.pt());
      _taugenE->push_back(matchToGen(*patTau).second.energy());
      _taugenEta->push_back(matchToGen(*patTau).second.eta());
      _taugenPhi->push_back(matchToGen(*patTau).second.phi());
      _egenPt->push_back(matchToGen(*patElectron).second.pt());
      _egenE->push_back(matchToGen(*patElectron).second.energy());
      _egenEta->push_back(matchToGen(*patElectron).second.eta());
      _egenPhi->push_back(matchToGen(*patElectron).second.phi());

      //  Get basic infomration of the tau candidate 
      _tauE->push_back(Tau.energy());    
      _tauEt->push_back((Tau.energy())*sin(Tau.theta()));									       
      _tauPt->push_back(Tau.pt());									       
      _tauCharge->push_back(patTau->charge());						       
      _tauEta->push_back(Tau.eta()); 						       
      _tauPhi->push_back(Tau.phi());
      // if patTau is a CaloTau, fill variables related to enery with dummy values (we use PFTaus)
      if (patTau->isCaloTau()){
        _tauNProngs->push_back(patTau->signalTracks().size());	       
        _tauEmFraction->push_back(-1.);									       
      	_tauHcalTotOverPLead->push_back(-1.);					    
      	_tauHcalMaxOverPLead->push_back(-1.); 									    
      	_tauHcal3x3OverPLead->push_back(-1.);									    
      	_tauElectronPreId->push_back(0);									    
      	_tauModifiedEOverP->push_back(-1.);						    
      	_tauBremsRecoveryEOverPLead->push_back(-1.);
        _tauDiscAgainstElec->push_back(patTau->tauID(_RecoTauDiscrAgainstElectron.data()));
        _tauIsInTheCraks->push_back(isInTheCracks(patTau->eta()));
        if(patTau->leadTrack().isNonnull()){
          _tauLTPt->push_back(patTau->leadTrack()->pt());
          _tauLTCharge->push_back(patTau->leadTrack()->charge());
          _tauLTChi2->push_back(patTau->leadTrack()->chi2());
          _tauLTRecHitsSize->push_back(patTau->leadTrack()->recHitsSize());
          _tauLTSignedIp->push_back(patTau->leadTracksignedSipt());
          _tauPdgId->push_back(getMatchedPdgId(patTau->leadTrack()->pt(), patTau->leadTrack()->eta(), patTau->leadTrack()->phi(), patTau->leadTrack()->charge()).first);
          _tauMotherPdgId->push_back(getMatchedPdgId(patTau->leadTrack()->pt(), patTau->leadTrack()->eta(), patTau->leadTrack()->phi(), patTau->leadTrack()->charge()).second);  
 
        }else{
          _tauLTPt->push_back(-1.);
          _tauLTCharge->push_back(-10);
          _tauLTChi2->push_back(100.);
          _tauLTRecHitsSize->push_back(0);
          _tauLTSignedIp->push_back(100);
          _tauPdgId->push_back(getMatchedPdgId(patTau->pt(), patTau->eta(), patTau->phi(), patTau->charge()).first);
          _tauMotherPdgId->push_back(getMatchedPdgId(patTau->pt(), patTau->eta(), patTau->phi(), patTau->charge()).second);
        }
        
      }else{
        // If there is a valid reference wrt a gen level object of a leadding PFCharged hadron candidate, then get all the track quility info
        // and the pdg particle and mother IDs     
        const PFCandidateRef& theLeadPFCand = patTau->leadPFChargedHadrCand();
        if(theLeadPFCand.isNonnull()){
	  _tauLTPt->push_back(theLeadPFCand->pt());
	  _tauLTCharge->push_back(theLeadPFCand->charge());
	  _tauLTChi2->push_back(theLeadPFCand->trackRef()->chi2());
	  _tauLTRecHitsSize->push_back(theLeadPFCand->trackRef()->recHitsSize());
	  _tauLTSignedIp->push_back(patTau->leadPFChargedHadrCandsignedSipt());
          _tauPdgId->push_back(getMatchedPdgId(theLeadPFCand->pt(), theLeadPFCand->eta(), theLeadPFCand->phi(), theLeadPFCand->charge()).first);
          _tauMotherPdgId->push_back(getMatchedPdgId(theLeadPFCand->pt(), theLeadPFCand->eta(), theLeadPFCand->phi(), theLeadPFCand->charge()).second); 
        }
        else {
          _tauLTPt->push_back(-1.);
          _tauLTCharge->push_back(-10);
          _tauLTChi2->push_back(100.);
          _tauLTRecHitsSize->push_back(0);
          _tauLTSignedIp->push_back(100);
          _tauPdgId->push_back(getMatchedPdgId(patTau->pt(), patTau->eta(), patTau->phi(), patTau->charge()).first);
          _tauMotherPdgId->push_back(getMatchedPdgId(patTau->pt(), patTau->eta(), patTau->phi(), patTau->charge()).second);
        } 									
        // Get all the energy related quantities of the tau (PFTau) and N prongs
        _tauNProngs->push_back(patTau->signalPFChargedHadrCands().size());				       
        _tauEmFraction->push_back(patTau->emFraction());							       
        _tauHcalTotOverPLead->push_back(patTau->hcalTotOverPLead());			       
        _tauHcalMaxOverPLead->push_back(patTau->hcalMaxOverPLead());					       
        _tauHcal3x3OverPLead->push_back(patTau->hcal3x3OverPLead());					       
        _tauElectronPreId->push_back(patTau->electronPreIDDecision());		       
        _tauModifiedEOverP->push_back(patTau->ecalStripSumEOverPLead());			       
        _tauBremsRecoveryEOverPLead->push_back(patTau->bremsRecoveryEOverPLead());       
        _tauDiscAgainstElec->push_back(patTau->tauID(_RecoTauDiscrAgainstElectron.data()));
	_tauIsInTheCraks->push_back(isInTheCracks(patTau->eta()));
      }

      // fill bool variable to determing if the tau candidate passed matching or not. 
      _tauMatched->push_back(int(matchToGen(*patTau).first));

      // Get tauIso using different cone sizes and/or thresholds 
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
      
      if(_SetTANC){
        _tauTancDiscOnePercent->push_back(patTau->tauID("byTaNCfrOnePercent"));
        _tauTancDiscHalfPercent->push_back(patTau->tauID("byTaNCfrHalfPercent"));
        _tauTancDiscQuarterPercent->push_back(patTau->tauID("byTaNCfrQuarterPercent"));
        _tauTancDiscTenthPercent->push_back(patTau->tauID("byTaNCfrTenthPercent"));
      }

      // ZeeVeto : look for 2 electrons that reconstruct the mass of the Z within 3 sigmas or that have a small Pt asymmetry.
      // if so, electron flaged as comming from Z->ee 
      _eventIsZee->push_back(int(isZee(Electron).first));
      _zeeMass->push_back(isZee(Electron).second.first);
      _zeePtAsymm->push_back(isZee(Electron).second.second);

      // fill bool variable to determing if the electron candidate passed matching or not.
      _eMatched->push_back(int(matchToGen(*patElectron).first));

      // Is the electron comming from the barrel or the endcap?
      _isEEv->push_back(patElectron->isEE());
      _isEBv->push_back(patElectron->isEB());
      // true if particle is in the crack between EB and EE
      _isEBEEGap   ->push_back(patElectron->isEBEEGap()) ;
      // true if particle is in EB, and inside the eta gaps between modules
      _isEBEtaGap  ->push_back(patElectron->isEBEtaGap());
      // true if particle is in EB, and inside the phi gaps between modules
      _isEBPhiGap  ->push_back(patElectron->isEBPhiGap());
      // true if particle is in EE, and inside the gaps between dees
      _isEEDeeGap  ->push_back(patElectron->isEEDeeGap());
      // true if particle is in EE, and inside the gaps between rings
      _isEERingGap ->push_back(patElectron->isEERingGap());

      // get Electron pdg Id 
      _ePdgId->push_back(getMatchedPdgId(patElectron->pt(), patElectron->eta(), patElectron->phi(), patElectron->charge()).first);
      _eMotherPdgId->push_back(getMatchedPdgId(patElectron->pt(), patElectron->eta(), patElectron->phi(), patElectron->charge()).second);
      // fill basic electron quantities
      _eE->push_back(Electron.energy());			       
      _ePt->push_back(Electron.pt());										       
      _eEt->push_back(patElectron->et());       
      _eCharge->push_back(patElectron->charge());	       
      _eEta->push_back(Electron.eta());										       
      _ePhi->push_back(Electron.phi());									       
      _eSigmaEtaEta->push_back(patElectron->scSigmaEtaEta());							       
      _eSigmaIEtaIEta->push_back(patElectron->scSigmaIEtaIEta());		       
      _eEOverP->push_back(patElectron->eSuperClusterOverP());					       
      _eHOverEm->push_back(patElectron->hadronicOverEm());		       
      _eDeltaPhiIn->push_back(patElectron->deltaEtaSuperClusterTrackAtVtx());	       
      _eDeltaEtaIn->push_back(patElectron->deltaPhiSuperClusterTrackAtVtx());
      // if set to tru use customized user Isolation, else use pat default values 

      _eEcalIsoPat->push_back(patElectron->ecalIso());
      _eHcalIsoPat->push_back(patElectron->hcalIso());										       
      _eTrkIsoPat->push_back(patElectron->trackIso());										       
     
      _eUserEcalIso->push_back(patElectron->userIsolation(pat::EcalIso));
      _eUserHcalIso->push_back(patElectron->userIsolation(pat::HcalIso));										     
      _eUserTrkIso->push_back(patElectron->userIsolation(pat::TrackIso));
      _eIsoPat->push_back(patElectron->caloIso());

      _eSCE1x5->push_back(patElectron->scE1x5());								       
      _eSCE2x5->push_back(patElectron->scE2x5Max());		       
      _eSCE5x5->push_back(patElectron->scE5x5()); 
      // fill the electrons impact parameter and its error
      const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
      if(patElectron->gsfTrack().isNonnull()) {
        _eIp->push_back(patElectron->gsfTrack()->dxy(thePrimaryEventVertex.position()));
        _eIpError->push_back(patElectron->gsfTrack()->d0Error()); 
      }
      else{
	_eIp->push_back(-100.);
	_eIpError->push_back(-100.);
      } 
      // get impact paramenter from of the CTF tracks
      if(patElectron->closestCtfTrackRef().isNonnull()){
	_eIp_ctf->push_back(patElectron->closestCtfTrackRef()->dxy(thePrimaryEventVertex.position()));
	_eIpError_ctf->push_back(patElectron->closestCtfTrackRef()->d0Error());
      } 
      else{ 
        _eIpError_ctf->push_back(-100.);
        _eIp_ctf->push_back(-100.); 
      }

      // Electron classification, is the electron classified as showering? narrow? big brem? golden? is in the cracks region? 
      _eClass->push_back(patElectron->classification());

      const Track* elTrack = (const reco::Track*)(patElectron->gsfTrack().get());
      // number of expected hits before the first track hit
      const HitPattern& pInner = elTrack->trackerExpectedHitsInner(); 
      _eMissingHits->push_back(pInner.numberOfHits());

      // Info for MET
      _mEt->push_back(MET.pt());
      // fill mass of the eTau visible producs												       
      reco::Candidate::LorentzVector The_LorentzVect = Tau + Electron;					       
      _eTauMass->push_back(The_LorentzVect.M());
      // CosDphi(e-Tau) between the electron and the tau 					       
      _eTauCosDPhi->push_back(cos(TMath::Abs(normalizedPhi(Electron.phi() - Tau.phi()))));	
      // DR between the electron and the tau       
      _eTauDelatR->push_back(reco::deltaR(Electron, Tau));
      
      //--------MET and Full Mass Reconstruction
      double full_px = Tau.px() + Electron.px() + MET.px();
      double full_py = Tau.py() + Electron.py() + MET.py();
      double full_pz = Tau.pz() + Electron.pz();
      double full_e  =  Tau.energy() + Electron.energy() + TMath::Sqrt((MET.px() * MET.px()) + (MET.py() * MET.py()));
      reco::Candidate::LorentzVector full_LorentzVector(full_px, full_py, full_pz, full_e);
      _eTauMetMass->push_back(full_LorentzVector.M());      

      // Colliniar mass approximation using the EDAnalyzer function
      _UseVectorSumOfVisProductsAndMetMassReco = false;
      _UseCollinerApproxMassReco = true;
      _eTauCollMass->push_back(CalculateThe4Momentum((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))).second.M());

      //=======Visible Pzeta
      double zetaX = cos(Tau.phi()) + cos(Electron.phi());
      double zetaY = sin(Tau.phi()) + sin(Electron.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = Tau.px() + Electron.px();
      double visPy = Tau.py() + Electron.py();
      double pZetaVis = visPx*zetaX + visPy*zetaY;
      _eTauPZetaVis->push_back(pZetaVis);
  
      //======Pzeta     
      double px = visPx + MET.px();
      double py = visPy + MET.py();
      double pZeta = px*zetaX + py*zetaY;
      _eTauPZeta->push_back(pZeta);  
    
    }										  
  }
  _HMTTree->Fill();


}

// ---------------Fill Histograms
void HiMassTauAnalysis::fillHistograms() {
  
  for(unsigned int NpdfID = 0; NpdfID < pdfWeightVector.size();  NpdfID++){

    // ------Vertices
    if (_FillRecoVertexHists) {
      int nVertices = 0;
      for(reco::VertexCollection::const_iterator primaryVertex = _primaryEventVertexCollection->begin();
	  primaryVertex != _primaryEventVertexCollection->end(); ++primaryVertex ) {
	if (!passRecoVertexCuts(*primaryVertex)) continue;
	_hVertexZposition[NpdfID]->Fill(primaryVertex->z(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	int pvntrk = 0;
	reco::Vertex::trackRef_iterator pvtrk;
	for(pvtrk=primaryVertex->tracks_begin();pvtrk!=primaryVertex->tracks_end();++pvtrk) {
	  if(primaryVertex->trackWeight(*pvtrk) > _RecoVertexTrackWeight) {pvntrk++;}
	}
	_hVertexNTracks[NpdfID]->Fill(pvntrk,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	nVertices++;
      }
      _hNVertices[NpdfID]->Fill(nVertices,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
    }
    
    // ------Generated Taus
    if ( (_FillGenTauHists) && (_GenParticleSource.label() != "") ) {
      reco::Candidate::LorentzVector theGenMotherObject(0,0,0,0);
      reco::Candidate::LorentzVector theGenGrandMotherObject(0,0,0,0);
      int nGenTaus = 0;
      for(GenParticleCollection::const_iterator genParticle = _genParticles->begin();genParticle != _genParticles->end();++genParticle) {
	if((abs(genParticle->pdgId()) == 15) && (genParticle->status() != 3)) {
	  int neutrinos = 0;
	  MChadtau = genParticle->p4();
	  if(genParticle->mother(0)->pdgId() == genParticle->pdgId()) {
	    motherCand = genParticle->mother(0)->mother(0);
	    theGenMotherObject = motherCand->p4();
	    if(motherCand->mother(0)->pdgId() == motherCand->pdgId()) {grandMotherCand = motherCand->mother(0)->mother(0);theGenGrandMotherObject = grandMotherCand->p4();}
	    else {grandMotherCand = motherCand->mother(0);theGenGrandMotherObject = grandMotherCand->p4();}
	  } else {
	    motherCand = genParticle->mother(0);
	    theGenMotherObject = motherCand->p4();
	    if(motherCand->mother(0)->pdgId() == motherCand->pdgId()) {grandMotherCand = motherCand->mother(0)->mother(0);theGenGrandMotherObject = grandMotherCand->p4();}
	    else {grandMotherCand = motherCand->mother(0);theGenGrandMotherObject = grandMotherCand->p4();}
	  }
	  for(int ii=0; ii<(int)(genParticle->numberOfDaughters()); ii++) {
	    daughterCand = genParticle->daughter(ii);
	    if( (abs(daughterCand->pdgId()) == 12) || (abs(daughterCand->pdgId()) == 14) || (abs(daughterCand->pdgId()) == 16) ) {
	      neutrinos++;
	      MChadtau = MChadtau - daughterCand->p4();
	    }
	  }
	  if(neutrinos == 1) {
	    if(_UseTauMotherId) {
	      if(abs(motherCand->pdgId()) == _TauMotherId) {
		if(_UseTauGrandMotherId) {
		  if(abs(grandMotherCand->pdgId()) == _TauGrandMotherId) {
		    _hGenTauEnergy[NpdfID]->Fill(MChadtau.energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		    _hGenTauPt[NpdfID]->Fill(MChadtau.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		    _hGenTauEta[NpdfID]->Fill(MChadtau.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		    _hGenTauPhi[NpdfID]->Fill(MChadtau.phi(),isrgluon_weight * isrgamma_weight * fsr_weight);
		    _hGenTauMotherEnergy[NpdfID]->Fill(theGenMotherObject.energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		    _hGenTauMotherPt[NpdfID]->Fill(theGenMotherObject.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		    _hGenTauMotherEta[NpdfID]->Fill(theGenMotherObject.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		    _hGenTauMotherPhi[NpdfID]->Fill(theGenMotherObject.phi(),isrgluon_weight * isrgamma_weight * fsr_weight);
		    _hGenTauGrandMotherEnergy[NpdfID]->Fill(theGenGrandMotherObject.energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		    _hGenTauGrandMotherPt[NpdfID]->Fill(theGenGrandMotherObject.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		    _hGenTauGrandMotherEta[NpdfID]->Fill(theGenGrandMotherObject.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		    _hGenTauGrandMotherPhi[NpdfID]->Fill(theGenGrandMotherObject.phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		    nGenTaus++;
		  }
		} else {
		  _hGenTauEnergy[NpdfID]->Fill(MChadtau.energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hGenTauPt[NpdfID]->Fill(MChadtau.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hGenTauEta[NpdfID]->Fill(MChadtau.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hGenTauPhi[NpdfID]->Fill(MChadtau.phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hGenTauMotherEnergy[NpdfID]->Fill(theGenMotherObject.energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hGenTauMotherPt[NpdfID]->Fill(theGenMotherObject.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hGenTauMotherEta[NpdfID]->Fill(theGenMotherObject.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hGenTauMotherPhi[NpdfID]->Fill(theGenMotherObject.phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hGenTauGrandMotherEnergy[NpdfID]->Fill(theGenGrandMotherObject.energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hGenTauGrandMotherPt[NpdfID]->Fill(theGenGrandMotherObject.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hGenTauGrandMotherEta[NpdfID]->Fill(theGenGrandMotherObject.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hGenTauGrandMotherPhi[NpdfID]->Fill(theGenGrandMotherObject.phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  nGenTaus++;
		}
	      }
	    } else {
	      _hGenTauEnergy[NpdfID]->Fill(MChadtau.energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hGenTauPt[NpdfID]->Fill(MChadtau.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hGenTauEta[NpdfID]->Fill(MChadtau.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hGenTauPhi[NpdfID]->Fill(MChadtau.phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hGenTauMotherEnergy[NpdfID]->Fill(theGenMotherObject.energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hGenTauMotherPt[NpdfID]->Fill(theGenMotherObject.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hGenTauMotherEta[NpdfID]->Fill(theGenMotherObject.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hGenTauMotherPhi[NpdfID]->Fill(theGenMotherObject.phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hGenTauGrandMotherEnergy[NpdfID]->Fill(theGenGrandMotherObject.energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hGenTauGrandMotherPt[NpdfID]->Fill(theGenGrandMotherObject.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hGenTauGrandMotherEta[NpdfID]->Fill(theGenGrandMotherObject.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hGenTauGrandMotherPhi[NpdfID]->Fill(theGenGrandMotherObject.phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      nGenTaus++;
	    }
	  }
	}
      }
      _hNGenTau[NpdfID]->Fill(nGenTaus,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
    }
    
    // ------Reco Tau Histograms
    if (_FillRecoTauHists) {
      int nTaus = 0;
      int theNumberOfTaus = 0;
      for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); 
	    patTau != _patTaus->end(); ++patTau ) {
	theNumberOfTaus++;
	if (!passRecoTauCuts((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1))) continue;
	if(_SmearTheTau) {
	  _hTauJetEnergy[NpdfID]->Fill(smearedTauMomentumVector.at(theNumberOfTaus-1).energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJetPt[NpdfID]->Fill(smearedTauMomentumVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJetEta[NpdfID]->Fill(smearedTauMomentumVector.at(theNumberOfTaus-1).eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJetPhi[NpdfID]->Fill(smearedTauMomentumVector.at(theNumberOfTaus-1).phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  if (patTau->isCaloTau()) {
	    _hTauJetNumSignalTracks[NpdfID]->Fill(patTau->signalTracks().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetNumSignalGammas[NpdfID]->Fill(0.0,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetCharge[NpdfID]->Fill(fabs(patTau->charge()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSignalTracksMass[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSignalTracksAndGammasMass[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSignalTracksChargeFraction[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).pt() / smearedTauMomentumVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    if(patTau->signalTracks().size() == 1) {_hTauJetSignalTracksMass1prong[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    if(patTau->signalTracks().size() == 3) {_hTauJetSignalTracksMass3prong[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    if(patTau->signalTracks().size() == 1) {_hTauJetSignalTracksAndGammasMass1prong[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    if(patTau->signalTracks().size() == 3) {_hTauJetSignalTracksAndGammasMass3prong[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    if( (patTau->leadTrack().isNonnull()) ) {
	      _hTauJetSeedTrackPt[NpdfID]->Fill(patTau->leadTrack()->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSeedTrackIpSignificance[NpdfID]->Fill(patTau->leadTracksignedSipt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      //            _hTauJetSeedTrackNhits[NpdfID]->Fill(patTau->leadTrack()->recHitsSize(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      //            _hTauJetSeedTrackChi2[NpdfID]->Fill(patTau->leadTrack()->chi2(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    }
	  } else {
	    _hTauJetNumSignalTracks[NpdfID]->Fill(patTau->signalPFChargedHadrCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetNumSignalGammas[NpdfID]->Fill(CalculateNumberSignalTauGammas(*patTau),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetCharge[NpdfID]->Fill(fabs(patTau->charge()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSignalTracksMass[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSignalTracksAndGammasMass[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSignalTracksChargeFraction[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).pt() / smearedTauMomentumVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    if(patTau->signalPFChargedHadrCands().size() == 1) {
	      _hTauJetSignalTracksMass1prong[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSignalTracksAndGammasMass1prong[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      if(CalculateNumberSignalTauGammas(*patTau) == 0) {_hTauJetMass1Prong0Gamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      if(CalculateNumberSignalTauGammas(*patTau) == 1) {_hTauJetMass1Prong1Gamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      if(CalculateNumberSignalTauGammas(*patTau) >= 2) {_hTauJetMass1Prong2orMoreGamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    }
	    if(patTau->signalPFChargedHadrCands().size() == 3) {
	      _hTauJetSignalTracksMass3prong[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSignalTracksAndGammasMass3prong[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      if(CalculateNumberSignalTauGammas(*patTau) == 0) {_hTauJetMass3Prong0Gamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      if(CalculateNumberSignalTauGammas(*patTau) == 1) {_hTauJetMass3Prong1Gamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      if(CalculateNumberSignalTauGammas(*patTau) >= 2) {_hTauJetMass3Prong2orMoreGamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    }
	    if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
	      _hTauJetSeedTrackPt[NpdfID]->Fill(patTau->leadPFChargedHadrCand()->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSeedTrackIpSignificance[NpdfID]->Fill(patTau->leadPFChargedHadrCandsignedSipt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      if( (patTau->leadPFChargedHadrCand()->trackRef().isNonnull()) ) {
		_hTauJetSeedTrackNhits[NpdfID]->Fill(patTau->leadPFChargedHadrCand()->trackRef()->recHitsSize(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hTauJetSeedTrackChi2[NpdfID]->Fill(patTau->leadPFChargedHadrCand()->trackRef()->chi2(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      }
	    }
	  }
	  if(_GenParticleSource.label() != "") {
	    if(matchToGen(*patTau).first) {
	      _hTauJetGenTauDeltaPhi[NpdfID]->Fill(normalizedPhi(smearedTauMomentumVector.at(theNumberOfTaus-1).phi() - matchToGen(*patTau).second.phi()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetGenTauDeltaEta[NpdfID]->Fill(smearedTauMomentumVector.at(theNumberOfTaus-1).eta() - matchToGen(*patTau).second.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetGenTauDeltaPt[NpdfID]->Fill((smearedTauMomentumVector.at(theNumberOfTaus-1).pt() - matchToGen(*patTau).second.pt()) / matchToGen(*patTau).second.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    }
	  }
	  if (_UseRecoTauDiscrByIsolationFlag) {
	    if (patTau->isCaloTau()) {
	      _hTauJetNumIsoTracks[NpdfID]->Fill(patTau->isolationTracks().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetNumIsoGammas[NpdfID]->Fill(0.0,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetNumIsoCands[NpdfID]->Fill(patTau->isolationTracks().size() + 0,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSumPtIsoTracks[NpdfID]->Fill(patTau->isolationTracksPtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSumPtIsoGammas[NpdfID]->Fill(patTau->isolationECALhitsEtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSumPtIso[NpdfID]->Fill(patTau->isolationTracksPtSum() + patTau->isolationECALhitsEtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    } else {
	      _hTauJetNumIsoTracks[NpdfID]->Fill(patTau->isolationPFChargedHadrCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetNumIsoGammas[NpdfID]->Fill(patTau->isolationPFGammaCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetNumIsoCands[NpdfID]->Fill(patTau->isolationPFChargedHadrCands().size() + patTau->isolationPFGammaCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSumPtIsoTracks[NpdfID]->Fill(patTau->isolationPFChargedHadrCandsPtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSumPtIsoGammas[NpdfID]->Fill(patTau->isolationPFGammaCandsEtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSumPtIso[NpdfID]->Fill(patTau->isolationPFChargedHadrCandsPtSum() + patTau->isolationPFGammaCandsEtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    }
	  } else {
	    _hTauJetNumIsoTracks[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetNumIsoGammas[NpdfID]->Fill(CalculateTauEcalIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetNumIsoCands[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).first + CalculateTauEcalIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIsoTracks[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIsoGammas[NpdfID]->Fill(CalculateTauEcalIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIso[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).second + CalculateTauEcalIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	} else {
	  _hTauJetEnergy[NpdfID]->Fill(patTau->energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJetPt[NpdfID]->Fill(patTau->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJetEta[NpdfID]->Fill(patTau->eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJetPhi[NpdfID]->Fill(patTau->phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  if (patTau->isCaloTau()) {
	    _hTauJetNumSignalTracks[NpdfID]->Fill(patTau->signalTracks().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetNumSignalGammas[NpdfID]->Fill(0.0,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetCharge[NpdfID]->Fill(fabs(patTau->charge()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSignalTracksMass[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSignalTracksAndGammasMass[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSignalTracksChargeFraction[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).pt() / patTau->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    if(patTau->signalTracks().size() == 1) {_hTauJetSignalTracksMass1prong[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    if(patTau->signalTracks().size() == 3) {_hTauJetSignalTracksMass3prong[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    if(patTau->signalTracks().size() == 1) {_hTauJetSignalTracksAndGammasMass1prong[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    if(patTau->signalTracks().size() == 3) {_hTauJetSignalTracksAndGammasMass3prong[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    if( (patTau->leadTrack().isNonnull()) ) {
	      _hTauJetSeedTrackPt[NpdfID]->Fill(patTau->leadTrack()->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSeedTrackIpSignificance[NpdfID]->Fill(patTau->leadTracksignedSipt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSeedTrackNhits[NpdfID]->Fill(patTau->leadTrack()->recHitsSize(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSeedTrackChi2[NpdfID]->Fill(patTau->leadTrack()->chi2(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    }
	  } else {
	    _hTauJetNumSignalTracks[NpdfID]->Fill(patTau->signalPFChargedHadrCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetNumSignalGammas[NpdfID]->Fill(CalculateNumberSignalTauGammas(*patTau),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetCharge[NpdfID]->Fill(fabs(patTau->charge()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSignalTracksMass[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSignalTracksAndGammasMass[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSignalTracksChargeFraction[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).pt() / patTau->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    if(patTau->signalPFChargedHadrCands().size() == 1) {
	      _hTauJetSignalTracksMass1prong[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSignalTracksAndGammasMass1prong[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      if(CalculateNumberSignalTauGammas(*patTau) == 0) {_hTauJetMass1Prong0Gamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      if(CalculateNumberSignalTauGammas(*patTau) == 1) {_hTauJetMass1Prong1Gamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      if(CalculateNumberSignalTauGammas(*patTau) >= 2) {_hTauJetMass1Prong2orMoreGamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    }
	    if(patTau->signalPFChargedHadrCands().size() == 3) {
	      _hTauJetSignalTracksMass3prong[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSignalTracksAndGammasMass3prong[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      if(CalculateNumberSignalTauGammas(*patTau) == 0) {_hTauJetMass3Prong0Gamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      if(CalculateNumberSignalTauGammas(*patTau) == 1) {_hTauJetMass3Prong1Gamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      if(CalculateNumberSignalTauGammas(*patTau) >= 2) {_hTauJetMass3Prong2orMoreGamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    }
	    if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
	      _hTauJetSeedTrackPt[NpdfID]->Fill(patTau->leadPFChargedHadrCand()->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSeedTrackIpSignificance[NpdfID]->Fill(patTau->leadPFChargedHadrCandsignedSipt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      if( (patTau->leadPFChargedHadrCand()->trackRef().isNonnull()) ) {
		_hTauJetSeedTrackNhits[NpdfID]->Fill(patTau->leadPFChargedHadrCand()->trackRef()->recHitsSize(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hTauJetSeedTrackChi2[NpdfID]->Fill(patTau->leadPFChargedHadrCand()->trackRef()->chi2(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      }
	    }
	  }
	  if(_GenParticleSource.label() != "") {
	    if(matchToGen(*patTau).first) {
	      _hTauJetGenTauDeltaPhi[NpdfID]->Fill(normalizedPhi(patTau->phi() - matchToGen(*patTau).second.phi()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetGenTauDeltaEta[NpdfID]->Fill(patTau->eta() - matchToGen(*patTau).second.eta(),isrgluon_weight * isrgamma_weight * fsr_weight);
	      _hTauJetGenTauDeltaPt[NpdfID]->Fill((patTau->pt() - matchToGen(*patTau).second.pt()) / matchToGen(*patTau).second.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    }
	  }
	  if (_UseRecoTauDiscrByIsolationFlag) {
	    if (patTau->isCaloTau()) {
	      _hTauJetNumIsoTracks[NpdfID]->Fill(patTau->isolationTracks().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetNumIsoGammas[NpdfID]->Fill(0.0,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetNumIsoCands[NpdfID]->Fill(patTau->isolationTracks().size() + 0,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSumPtIsoTracks[NpdfID]->Fill(patTau->isolationTracksPtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSumPtIsoGammas[NpdfID]->Fill(patTau->isolationECALhitsEtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSumPtIso[NpdfID]->Fill(patTau->isolationTracksPtSum() + patTau->isolationECALhitsEtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    } else {
	      _hTauJetNumIsoTracks[NpdfID]->Fill(patTau->isolationPFChargedHadrCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetNumIsoGammas[NpdfID]->Fill(patTau->isolationPFGammaCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetNumIsoCands[NpdfID]->Fill(patTau->isolationPFChargedHadrCands().size() + patTau->isolationPFGammaCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSumPtIsoTracks[NpdfID]->Fill(patTau->isolationPFChargedHadrCandsPtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSumPtIsoGammas[NpdfID]->Fill(patTau->isolationPFGammaCandsEtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauJetSumPtIso[NpdfID]->Fill(patTau->isolationPFChargedHadrCandsPtSum() + patTau->isolationPFGammaCandsEtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    }
	  } else {
	    _hTauJetNumIsoTracks[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetNumIsoGammas[NpdfID]->Fill(CalculateTauEcalIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetNumIsoCands[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).first + CalculateTauEcalIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIsoTracks[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIsoGammas[NpdfID]->Fill(CalculateTauEcalIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIso[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).second + CalculateTauEcalIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
	nTaus++;
      }
      _hNTau[NpdfID]->Fill(nTaus,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
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
	  _hMuonEnergy[NpdfID]->Fill(smearedMuonMomentumVector.at(theNumberOfMuons-1).energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hMuonPt[NpdfID]->Fill(smearedMuonMomentumVector.at(theNumberOfMuons-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hMuonEta[NpdfID]->Fill(smearedMuonMomentumVector.at(theNumberOfMuons-1).eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hMuonPhi[NpdfID]->Fill(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  if(_GenParticleSource.label() != "") {
	    if(matchToGen(*patMuon).first) {
	      _hMuonGenMuonDeltaPhi[NpdfID]->Fill(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi() - matchToGen(*patMuon).second.phi()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuonGenMuonDeltaEta[NpdfID]->Fill(smearedMuonMomentumVector.at(theNumberOfMuons-1).eta() - matchToGen(*patMuon).second.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuonGenMuonDeltaPt[NpdfID]->Fill((smearedMuonMomentumVector.at(theNumberOfMuons-1).pt() - matchToGen(*patMuon).second.pt()) / matchToGen(*patMuon).second.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    }
	  }
	} else {
	  _hMuonEnergy[NpdfID]->Fill(patMuon->energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hMuonPt[NpdfID]->Fill(patMuon->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hMuonEta[NpdfID]->Fill(patMuon->eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hMuonPhi[NpdfID]->Fill(patMuon->phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  if(_GenParticleSource.label() != "") {
	    if(matchToGen(*patMuon).first) {
	      _hMuonGenMuonDeltaPhi[NpdfID]->Fill(normalizedPhi(patMuon->phi() - matchToGen(*patMuon).second.phi()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuonGenMuonDeltaEta[NpdfID]->Fill(patMuon->eta() - matchToGen(*patMuon).second.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuonGenMuonDeltaPt[NpdfID]->Fill((patMuon->pt() - matchToGen(*patMuon).second.pt()) / matchToGen(*patMuon).second.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    }
	  }
	}
	_hMuonTrackIso[NpdfID]->Fill(patMuon->trackIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonTrackIsoTrkThreshold).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuonEcalIso[NpdfID]->Fill(patMuon->ecalIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonEcalIsoRecHitThreshold).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuonIso[NpdfID]->Fill(patMuon->trackIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonTrackIsoTrkThreshold).first +
				patMuon->ecalIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonEcalIsoRecHitThreshold).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	/*
	  _hMuonTrackIso[NpdfID]->Fill(patMuon->trackIso(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hMuonEcalIso[NpdfID]->Fill(patMuon->ecalIso(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hMuonIso[NpdfID]->Fill(patMuon->trackIso() + patMuon->ecalIso(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	*/
	_hMuonCaloCompatibility[NpdfID]->Fill(muon::caloCompatibility(*patMuon),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuonSegmentCompatibility[NpdfID]->Fill(muon::segmentCompatibility(*patMuon),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuonCaloCompatibilityVsSegmentCompatibility[NpdfID]->Fill(muon::caloCompatibility(*patMuon),muon::segmentCompatibility(*patMuon),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuonAntiPion[NpdfID]->Fill((_RecoMuonCaloCompCoefficient * muon::caloCompatibility(*patMuon)) + (_RecoMuonSegmCompCoefficient * muon::segmentCompatibility(*patMuon)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
	if ( patMuon->track().isNonnull() ) {
	  _hMuonIp[NpdfID]->Fill( patMuon->track()->dxy(thePrimaryEventVertex.position()) ,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  if(fabs(patMuon->track()->dxyError()) != 0) {
	    _hMuonIpSignificance[NpdfID]->Fill( fabs(patMuon->track()->dxy(thePrimaryEventVertex.position()) / patMuon->track()->dxyError()) ,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  } else {_hMuonIpSignificance[NpdfID]->Fill(-1,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	}
	nMuons++;
      }
      _hNMuon[NpdfID]->Fill(nMuons,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
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
	  _hElectronEnergy[NpdfID]->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hElectronPt[NpdfID]->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hElectronEta[NpdfID]->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hElectronPhi[NpdfID]->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  if(_GenParticleSource.label() != "") {
	    if(matchToGen(*patElectron).first) {
	      _hElectronGenElectronDeltaPhi[NpdfID]->Fill(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - matchToGen(*patElectron).second.phi()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectronGenElectronDeltaEta[NpdfID]->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).eta() - matchToGen(*patElectron).second.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectronGenElectronDeltaPt[NpdfID]->Fill((smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt() - matchToGen(*patElectron).second.pt()) / matchToGen(*patElectron).second.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    }
	  }
	} else {
	  _hElectronEnergy[NpdfID]->Fill(patElectron->energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hElectronPt[NpdfID]->Fill(patElectron->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hElectronEta[NpdfID]->Fill(patElectron->eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hElectronPhi[NpdfID]->Fill(patElectron->phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  if(_GenParticleSource.label() != "") {
	    if(matchToGen(*patElectron).first) {
	      _hElectronGenElectronDeltaPhi[NpdfID]->Fill(normalizedPhi(patElectron->phi() - matchToGen(*patElectron).second.phi()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectronGenElectronDeltaEta[NpdfID]->Fill(patElectron->eta() - matchToGen(*patElectron).second.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectronGenElectronDeltaPt[NpdfID]->Fill((patElectron->pt() - matchToGen(*patElectron).second.pt()) / matchToGen(*patElectron).second.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    }
	  }
	}
	_hElectronTrackIso[NpdfID]->Fill(patElectron->trackIso(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectronEcalIso[NpdfID]->Fill(patElectron->ecalIso(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	/*
	  _hElectronTrackIso[NpdfID]->Fill(patElectron->trackIsoDeposit()->depositAndCountWithin(_RecoElectronTrackIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoElectronTrackIsoTrkThreshold).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hElectronEcalIso[NpdfID]->Fill(patElectron->ecalIsoDeposit()->depositAndCountWithin(_RecoElectronEcalIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoElectronEcalIsoRecHitThreshold).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	*/
	const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
	if ( patElectron->track().isNonnull() ) { _hElectronIp[NpdfID]->Fill( patElectron->track()->dxy(thePrimaryEventVertex.position()) ,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID)); }
	//      if ( patElectron->gsfTrack().isNonnull() ) { _hElectronIp[NpdfID]->Fill( patElectron->gsfTrack()->dxy(thePrimaryEventVertex.position()) ,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID)); }
	_hElectronEoverP[NpdfID]->Fill(patElectron->eSuperClusterOverP(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectronHoverEm[NpdfID]->Fill(patElectron->hadronicOverEm(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectronClassification[NpdfID]->Fill(patElectron->classification(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	//_hElectronEcalDriven[NpdfID]->Fill(patElectron->ecalDrivenSeed(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	//_hElectronTrackerDriven[NpdfID]->Fill(patElectron->trackerDrivenSeed(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	nElectrons++;
      }
      _hNElectron[NpdfID]->Fill(nElectrons,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
    }
    
    // ------Reco Jet Histograms
    if (_FillRecoJetHists) {
      int nJets = 0;
      int theNumberOfJets = 0;
      for ( pat::JetCollection::const_iterator patJet = _patJets->begin(); 
	    patJet != _patJets->end(); ++patJet ) {
	theNumberOfJets++;
	if (!passRecoJetCuts((*patJet),_SmearTheJet,smearedJetMomentumVector.at(theNumberOfJets-1))) continue;
	if(_UseCorrectedJet) {
	  if(_SmearTheJet) {
	    _hJetEnergy[NpdfID]->Fill(smearedJetMomentumVector.at(theNumberOfJets-1).energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hJetPt[NpdfID]->Fill(smearedJetMomentumVector.at(theNumberOfJets-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hJetEta[NpdfID]->Fill(smearedJetMomentumVector.at(theNumberOfJets-1).eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hJetPhi[NpdfID]->Fill(smearedJetMomentumVector.at(theNumberOfJets-1).phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  } else {
	    _hJetEnergy[NpdfID]->Fill(patJet->energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hJetPt[NpdfID]->Fill(patJet->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hJetEta[NpdfID]->Fill(patJet->eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hJetPhi[NpdfID]->Fill(patJet->phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	  nJets++;
	} else {
	  if(_SmearTheJet) {
	    _hJetEnergy[NpdfID]->Fill(smearedJetMomentumVector.at(theNumberOfJets-1).energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hJetPt[NpdfID]->Fill(smearedJetMomentumVector.at(theNumberOfJets-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hJetEta[NpdfID]->Fill(smearedJetMomentumVector.at(theNumberOfJets-1).eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hJetPhi[NpdfID]->Fill(smearedJetMomentumVector.at(theNumberOfJets-1).phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  } else {
	    _hJetEnergy[NpdfID]->Fill(patJet->correctedJet("raw","").energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hJetPt[NpdfID]->Fill(patJet->correctedJet("raw","").pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hJetEta[NpdfID]->Fill(patJet->correctedJet("raw","").eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hJetPhi[NpdfID]->Fill(patJet->correctedJet("raw","").phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	  nJets++;
	}
	_hBJetDiscrByTrackCounting[NpdfID]->Fill(patJet->bDiscriminator("trackCountingHighEffBJetTags"),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBJetDiscrBySimpleSecondaryV[NpdfID]->Fill(patJet->bDiscriminator("simpleSecondaryVertexBJetTags"),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBJetDiscrByCombinedSecondaryV[NpdfID]->Fill(patJet->bDiscriminator("combinedSecondaryVertexBJetTags"),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
      }
      _hNJet[NpdfID]->Fill(nJets,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
    }
    
    // ------Topology Histograms
    if (_FillTopologyHists) {
      _hMet[NpdfID]->Fill(TMath::Sqrt((((*(_patMETs->begin())).px() + deltaForMEx)*((*(_patMETs->begin())).px() + deltaForMEx)) + (((*(_patMETs->begin())).py() + deltaForMEy)*((*(_patMETs->begin())).py() + deltaForMEy))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
      if( ((_AnalyzeMuonForLeg1) && (_AnalyzeTauForLeg2)) || ((_AnalyzeMuonForLeg2) && (_AnalyzeTauForLeg1)) ) {
	int theNumberOfMuons = 0;
	for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
	  theNumberOfMuons++;
	  int theNumberOfTaus = 0;
	  for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
	    theNumberOfTaus++;
	    if ((passRecoTauCuts((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1))) &&
		(passRecoMuonCuts((*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1))) &&
		(passTopologyCuts((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))))) {
	      if(_SmearTheTau) {
		if(_SmearTheMuon) {
		  _hMuonPtVsTauPt[NpdfID]->Fill(smearedMuonMomentumVector.at(theNumberOfMuons-1).pt(),smearedTauMomentumVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hMuonTauDeltaR[NpdfID]->Fill(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), smearedMuonMomentumVector.at(theNumberOfMuons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hMuonTauDeltaPtDivSumPt[NpdfID]->Fill((smearedTauMomentumVector.at(theNumberOfTaus-1).pt() - smearedMuonMomentumVector.at(theNumberOfMuons-1).pt()) / (smearedTauMomentumVector.at(theNumberOfTaus-1).pt() + smearedMuonMomentumVector.at(theNumberOfMuons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hMuonTauCosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi() - smearedTauMomentumVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  reco::Candidate::LorentzVector Temp_LorentzVect(
								  (*(_patMETs->begin())).px() + patMuon->px() - smearedMuonMomentumVector.at(theNumberOfMuons-1).px() + patTau->px() - smearedTauMomentumVector.at(theNumberOfTaus-1).px() + deltaForMEx,
								  (*(_patMETs->begin())).py() + patMuon->py() - smearedMuonMomentumVector.at(theNumberOfMuons-1).py() + patTau->py() - smearedTauMomentumVector.at(theNumberOfTaus-1).py() + deltaForMEy,
								  (*(_patMETs->begin())).pz(),
								  TMath::Sqrt( (((*(_patMETs->begin())).px() + deltaForMEx + patMuon->px() - smearedMuonMomentumVector.at(theNumberOfMuons-1).px() + patTau->px() - smearedTauMomentumVector.at(theNumberOfTaus-1).px()) *
										((*(_patMETs->begin())).px() + deltaForMEx + patMuon->px() - smearedMuonMomentumVector.at(theNumberOfMuons-1).px() + patTau->px() - smearedTauMomentumVector.at(theNumberOfTaus-1).px())) +
									       (((*(_patMETs->begin())).py() + deltaForMEy + patMuon->py() - smearedMuonMomentumVector.at(theNumberOfMuons-1).py() + patTau->py() - smearedTauMomentumVector.at(theNumberOfTaus-1).py()) *
										((*(_patMETs->begin())).py() + deltaForMEy + patMuon->py() - smearedMuonMomentumVector.at(theNumberOfMuons-1).py() + patTau->py() - smearedTauMomentumVector.at(theNumberOfTaus-1).py())) ));
		  _hMuonMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi() - Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hMuonMetDeltaPhiVsMuonTauCosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi() - Temp_LorentzVect.phi())), cos(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi() - smearedTauMomentumVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hTauMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauMomentumVector.at(theNumberOfTaus-1).phi() - Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		} else {
		  _hMuonPtVsTauPt[NpdfID]->Fill(patMuon->pt(),smearedTauMomentumVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hMuonTauDeltaR[NpdfID]->Fill(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), patMuon->p4()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hMuonTauDeltaPtDivSumPt[NpdfID]->Fill((smearedTauMomentumVector.at(theNumberOfTaus-1).pt() - patMuon->pt()) / (smearedTauMomentumVector.at(theNumberOfTaus-1).pt() + patMuon->pt()),isrgluon_weight * isrgamma_weight * fsr_weight);
		  _hMuonTauCosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(patMuon->phi() - smearedTauMomentumVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  reco::Candidate::LorentzVector Temp_LorentzVect(
								  (*(_patMETs->begin())).px() + patMuon->px() - patMuon->px() + patTau->px() - smearedTauMomentumVector.at(theNumberOfTaus-1).px() + deltaForMEx,
								  (*(_patMETs->begin())).py() + patMuon->py() - patMuon->py() + patTau->py() - smearedTauMomentumVector.at(theNumberOfTaus-1).py() + deltaForMEy,
								  (*(_patMETs->begin())).pz(),
								  TMath::Sqrt( (((*(_patMETs->begin())).px() + deltaForMEx + patMuon->px() - patMuon->px() + patTau->px() - smearedTauMomentumVector.at(theNumberOfTaus-1).px()) *
										((*(_patMETs->begin())).px() + deltaForMEx + patMuon->px() - patMuon->px() + patTau->px() - smearedTauMomentumVector.at(theNumberOfTaus-1).px())) +
									       (((*(_patMETs->begin())).py() + deltaForMEy + patMuon->py() - patMuon->py() + patTau->py() - smearedTauMomentumVector.at(theNumberOfTaus-1).py()) *
										((*(_patMETs->begin())).py() + deltaForMEy + patMuon->py() - patMuon->py() + patTau->py() - smearedTauMomentumVector.at(theNumberOfTaus-1).py())) ));
		  _hMuonMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patMuon->phi() - Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hMuonMetDeltaPhiVsMuonTauCosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patMuon->phi() - Temp_LorentzVect.phi())), cos(TMath::Abs(normalizedPhi(patMuon->phi() - smearedTauMomentumVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hTauMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauMomentumVector.at(theNumberOfTaus-1).phi() - Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		}
	      } else {
		if(_SmearTheMuon) {
		  _hMuonPtVsTauPt[NpdfID]->Fill(smearedMuonMomentumVector.at(theNumberOfMuons-1).pt(),patTau->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hMuonTauDeltaR[NpdfID]->Fill(reco::deltaR(patTau->p4(), smearedMuonMomentumVector.at(theNumberOfMuons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hMuonTauDeltaPtDivSumPt[NpdfID]->Fill((patTau->pt() - smearedMuonMomentumVector.at(theNumberOfMuons-1).pt()) / (patTau->pt() + smearedMuonMomentumVector.at(theNumberOfMuons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hMuonTauCosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi() - patTau->phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  reco::Candidate::LorentzVector Temp_LorentzVect(
								  (*(_patMETs->begin())).px() + patMuon->px() - smearedMuonMomentumVector.at(theNumberOfMuons-1).px() + deltaForMEx,
								  (*(_patMETs->begin())).py() + patMuon->py() - smearedMuonMomentumVector.at(theNumberOfMuons-1).py() + deltaForMEy,
								  (*(_patMETs->begin())).pz(),
								  TMath::Sqrt( (((*(_patMETs->begin())).px() + deltaForMEx + patMuon->px() - smearedMuonMomentumVector.at(theNumberOfMuons-1).px()) *
										((*(_patMETs->begin())).px() + deltaForMEx + patMuon->px() - smearedMuonMomentumVector.at(theNumberOfMuons-1).px())) +
									       (((*(_patMETs->begin())).py() + deltaForMEy + patMuon->py() - smearedMuonMomentumVector.at(theNumberOfMuons-1).py()) *
										((*(_patMETs->begin())).py() + deltaForMEy + patMuon->py() - smearedMuonMomentumVector.at(theNumberOfMuons-1).py())) ));
		  _hMuonMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi() - Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hMuonMetDeltaPhiVsMuonTauCosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi() - Temp_LorentzVect.phi())), cos(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons-1).phi() - patTau->phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hTauMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patTau->phi() - Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		} else {
		  _hMuonPtVsTauPt[NpdfID]->Fill(patMuon->pt(),patTau->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hMuonTauDeltaR[NpdfID]->Fill(reco::deltaR(patTau->p4(), patMuon->p4()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hMuonTauDeltaPtDivSumPt[NpdfID]->Fill((patTau->pt() - patMuon->pt()) / (patTau->pt() + patMuon->pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hMuonTauCosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(patMuon->phi() - patTau->phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hMuonMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patMuon->phi() - (*(_patMETs->begin())).phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hMuonMetDeltaPhiVsMuonTauCosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patMuon->phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(patMuon->phi() - patTau->phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hTauMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patTau->phi() - (*(_patMETs->begin())).phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		}
	      }
	      if(CalculateThe4Momentum((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))).first) {_hReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      else {_hNotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      _hMuonMetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauMetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      if(_UseTauSeedTrackForDiTauDiscrByOSLS) {
		if (patTau->isCaloTau()) {
		  if( (patTau->leadTrack().isNonnull()) ) {
		    _hMuonTauOSLS[NpdfID]->Fill(patMuon->charge() * patTau->leadTrack()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  }
		} else {
		  if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
		    _hMuonTauOSLS[NpdfID]->Fill(patMuon->charge() * patTau->leadPFChargedHadrCand()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  }
		}
	      } else {_hMuonTauOSLS[NpdfID]->Fill(patMuon->charge() * patTau->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      _hPZeta[NpdfID]->Fill(CalculatePZeta((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hPZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta2D[NpdfID]->Fill(CalculatePZetaVis((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))),CalculatePZeta((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta1D[NpdfID]->Fill((_PZetaCutCoefficient * CalculatePZeta((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin())))) + 
				     (_PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patMuon),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons-1),(*(_patMETs->begin())))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    }
	  }
	}
      }
      if( ((_AnalyzeElectronForLeg1) && (_AnalyzeTauForLeg2)) || ((_AnalyzeElectronForLeg2) && (_AnalyzeTauForLeg1)) ) {
	int theNumberOfElectrons = 0;
	for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
	  theNumberOfElectrons++;
	  int theNumberOfTaus = 0;
	  for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
	    theNumberOfTaus++;
	    if ((passRecoTauCuts((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1))) &&
		(passRecoElectronCuts((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1))) &&
		(passTopologyCuts((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))))) {
	      if(_SmearTheTau) {
		if(_SmearTheElectron) {
                  _hElectronIsZee[NpdfID]->Fill(isZee(smearedElectronMomentumVector.at(theNumberOfElectrons-1)).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronPtVsTauPt[NpdfID]->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt(),smearedTauMomentumVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronTauDeltaR[NpdfID]->Fill(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), smearedElectronMomentumVector.at(theNumberOfElectrons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronTauDeltaPtDivSumPt[NpdfID]->Fill((smearedTauMomentumVector.at(theNumberOfTaus-1).pt() - smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt()) / (smearedTauMomentumVector.at(theNumberOfTaus-1).pt() + smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronTauCosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - smearedTauMomentumVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  reco::Candidate::LorentzVector Temp_LorentzVect(
								  (*(_patMETs->begin())).px() + patElectron->px() - smearedElectronMomentumVector.at(theNumberOfElectrons-1).px() + patTau->px() - smearedTauMomentumVector.at(theNumberOfTaus-1).px() + deltaForMEx,
								  (*(_patMETs->begin())).py() + patElectron->py() - smearedElectronMomentumVector.at(theNumberOfElectrons-1).py() + patTau->py() - smearedTauMomentumVector.at(theNumberOfTaus-1).py() + deltaForMEy,
								  (*(_patMETs->begin())).pz(),
								  TMath::Sqrt( (((*(_patMETs->begin())).px() + deltaForMEx + patElectron->px() - smearedElectronMomentumVector.at(theNumberOfElectrons-1).px() + patTau->px() - smearedTauMomentumVector.at(theNumberOfTaus-1).px()) *
										((*(_patMETs->begin())).px() + deltaForMEx + patElectron->px() - smearedElectronMomentumVector.at(theNumberOfElectrons-1).px() + patTau->px() - smearedTauMomentumVector.at(theNumberOfTaus-1).px())) +
									       (((*(_patMETs->begin())).py() + deltaForMEy + patElectron->py() - smearedElectronMomentumVector.at(theNumberOfElectrons-1).py() + patTau->py() - smearedTauMomentumVector.at(theNumberOfTaus-1).py()) *
										((*(_patMETs->begin())).py() + deltaForMEy + patElectron->py() - smearedElectronMomentumVector.at(theNumberOfElectrons-1).py() + patTau->py() - smearedTauMomentumVector.at(theNumberOfTaus-1).py())) ));
		  _hElectronMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronMetDeltaPhiVsElectronTauCosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - Temp_LorentzVect.phi())), cos(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - smearedTauMomentumVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hTauMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauMomentumVector.at(theNumberOfTaus-1).phi() - Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		} else {
                  _hElectronIsZee[NpdfID]->Fill(isZee(patElectron->p4()).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronPtVsTauPt[NpdfID]->Fill(patElectron->pt(),smearedTauMomentumVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronTauDeltaR[NpdfID]->Fill(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), patElectron->p4()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronTauDeltaPtDivSumPt[NpdfID]->Fill((smearedTauMomentumVector.at(theNumberOfTaus-1).pt() - patElectron->pt()) / (smearedTauMomentumVector.at(theNumberOfTaus-1).pt() + patElectron->pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronTauCosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(patElectron->phi() - smearedTauMomentumVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  reco::Candidate::LorentzVector Temp_LorentzVect(
								  (*(_patMETs->begin())).px() + patElectron->px() - patElectron->px() + patTau->px() - smearedTauMomentumVector.at(theNumberOfTaus-1).px() + deltaForMEx,
								  (*(_patMETs->begin())).py() + patElectron->py() - patElectron->py() + patTau->py() - smearedTauMomentumVector.at(theNumberOfTaus-1).py() + deltaForMEy,
								  (*(_patMETs->begin())).pz(),
								  TMath::Sqrt( (((*(_patMETs->begin())).px() + deltaForMEx + patElectron->px() - patElectron->px() + patTau->px() - smearedTauMomentumVector.at(theNumberOfTaus-1).px()) *
										((*(_patMETs->begin())).px() + deltaForMEx + patElectron->px() - patElectron->px() + patTau->px() - smearedTauMomentumVector.at(theNumberOfTaus-1).px())) +
									       (((*(_patMETs->begin())).py() + deltaForMEy + patElectron->py() - patElectron->py() + patTau->py() - smearedTauMomentumVector.at(theNumberOfTaus-1).py()) *
										((*(_patMETs->begin())).py() + deltaForMEy + patElectron->py() - patElectron->py() + patTau->py() - smearedTauMomentumVector.at(theNumberOfTaus-1).py())) ));
		  _hElectronMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patElectron->phi() - Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronMetDeltaPhiVsElectronTauCosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patElectron->phi() - Temp_LorentzVect.phi())), cos(TMath::Abs(normalizedPhi(patElectron->phi() - smearedTauMomentumVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight);
		  _hTauMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauMomentumVector.at(theNumberOfTaus-1).phi() - Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		}
	      } else {
		if(_SmearTheElectron) {
                  _hElectronIsZee[NpdfID]->Fill(isZee(smearedElectronMomentumVector.at(theNumberOfElectrons-1)).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronPtVsTauPt[NpdfID]->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt(),patTau->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronTauDeltaR[NpdfID]->Fill(reco::deltaR(patTau->p4(), smearedElectronMomentumVector.at(theNumberOfElectrons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronTauDeltaPtDivSumPt[NpdfID]->Fill((patTau->pt() - smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt()) / (patTau->pt() + smearedElectronMomentumVector.at(theNumberOfElectrons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronTauCosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - patTau->phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  reco::Candidate::LorentzVector Temp_LorentzVect(
								  (*(_patMETs->begin())).px() + patElectron->px() - smearedElectronMomentumVector.at(theNumberOfElectrons-1).px() + deltaForMEx,
								  (*(_patMETs->begin())).py() + patElectron->py() - smearedElectronMomentumVector.at(theNumberOfElectrons-1).py() + deltaForMEy,
								  (*(_patMETs->begin())).pz(),
								  TMath::Sqrt( (((*(_patMETs->begin())).px() + deltaForMEx + patElectron->px() - smearedElectronMomentumVector.at(theNumberOfElectrons-1).px()) *
										((*(_patMETs->begin())).px() + deltaForMEx + patElectron->px() - smearedElectronMomentumVector.at(theNumberOfElectrons-1).px())) +
									       (((*(_patMETs->begin())).py() + deltaForMEy + patElectron->py() - smearedElectronMomentumVector.at(theNumberOfElectrons-1).py()) *
										((*(_patMETs->begin())).py() + deltaForMEy + patElectron->py() - smearedElectronMomentumVector.at(theNumberOfElectrons-1).py())) ));
		  _hElectronMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronMetDeltaPhiVsElectronTauCosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - Temp_LorentzVect.phi())), cos(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons-1).phi() - patTau->phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hTauMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patTau->phi() - Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		} else {
                  _hElectronIsZee[NpdfID]->Fill(isZee(patElectron->p4()).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronPtVsTauPt[NpdfID]->Fill(patElectron->pt(),patTau->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronTauDeltaR[NpdfID]->Fill(reco::deltaR(patTau->p4(), patElectron->p4()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronTauDeltaPtDivSumPt[NpdfID]->Fill((patTau->pt() - patElectron->pt()) / (patTau->pt() + patElectron->pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronTauCosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(patElectron->phi() - patTau->phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patElectron->phi() - (*(_patMETs->begin())).phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hElectronMetDeltaPhiVsElectronTauCosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patElectron->phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(patElectron->phi() - patTau->phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  _hTauMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patTau->phi() - (*(_patMETs->begin())).phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		}
	      }
	      if(CalculateThe4Momentum((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))).first) {_hReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      else {_hNotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      _hElectronMetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauMetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      if(_UseTauSeedTrackForDiTauDiscrByOSLS) {
		if (patTau->isCaloTau()) {
		  if( (patTau->leadTrack().isNonnull()) ) {
		    _hElectronTauOSLS[NpdfID]->Fill(patElectron->charge() * patTau->leadTrack()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  }
		} else {
		  if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
		    _hElectronTauOSLS[NpdfID]->Fill(patElectron->charge() * patTau->leadPFChargedHadrCand()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  }
		}
	      } else {_hElectronTauOSLS[NpdfID]->Fill(patElectron->charge() * patTau->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      _hPZeta[NpdfID]->Fill(CalculatePZeta((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hPZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta2D[NpdfID]->Fill(CalculatePZetaVis((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))),CalculatePZeta((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta1D[NpdfID]->Fill((_PZetaCutCoefficient * CalculatePZeta((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin())))) + (_PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus-1),(*patElectron),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons-1),(*(_patMETs->begin())))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
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
		_hMuon1PtVsMuon2Pt[NpdfID]->Fill(smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).pt(),smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hMuon1Muon2DeltaR[NpdfID]->Fill(reco::deltaR(smearedMuonMomentumVector.at(theNumberOfMuons1 - 1), smearedMuonMomentumVector.at(theNumberOfMuons2 - 1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hMuon1Muon2DeltaPtDivSumPt[NpdfID]->Fill((smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).pt() - smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).pt()) / (smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).pt() + smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hMuon1Muon2CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).phi() - smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		reco::Candidate::LorentzVector Temp_LorentzVect(
								(*(_patMETs->begin())).px() + patMuon1->px() - smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).px() + patMuon2->px() - smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).px() + deltaForMEx,
								(*(_patMETs->begin())).py() + patMuon1->py() - smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).py() + patMuon2->py() - smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).py() + deltaForMEy,
								(*(_patMETs->begin())).pz(),
								TMath::Sqrt( (((*(_patMETs->begin())).px() + deltaForMEx + patMuon1->px() - smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).px() + patMuon2->px() - smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).px()) *
									      ((*(_patMETs->begin())).px() + deltaForMEx + patMuon1->px() - smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).px() + patMuon2->px() - smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).px())) +
									     (((*(_patMETs->begin())).py() + deltaForMEy + patMuon1->py() - smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).py() + patMuon2->py() - smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).py()) *
									      ((*(_patMETs->begin())).py() + deltaForMEy + patMuon1->py() - smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).py() + patMuon2->py() - smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).py())) ));
		_hMuon1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).phi() - Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hMuon2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).phi() - Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hMuon1MetDeltaPhiVsMuon1Muon2CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).phi() - Temp_LorentzVect.phi())), cos(TMath::Abs(normalizedPhi(smearedMuonMomentumVector.at(theNumberOfMuons1 - 1).phi() - smearedMuonMomentumVector.at(theNumberOfMuons2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      } else {
		_hMuon1PtVsMuon2Pt[NpdfID]->Fill(patMuon1->pt(),patMuon2->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hMuon1Muon2DeltaR[NpdfID]->Fill(reco::deltaR(patMuon1->p4(), patMuon2->p4()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hMuon1Muon2DeltaPtDivSumPt[NpdfID]->Fill((patMuon1->pt() - patMuon2->pt()) / (patMuon1->pt() + patMuon2->pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hMuon1Muon2CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(patMuon1->phi() - patMuon2->phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hMuon1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patMuon1->phi() - (*(_patMETs->begin())).phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hMuon2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patMuon2->phi() - (*(_patMETs->begin())).phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hMuon1MetDeltaPhiVsMuon1Muon2CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patMuon1->phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(patMuon1->phi() - patMuon2->phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      }
	      if(CalculateThe4Momentum((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))).first) {_hReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      else {_hNotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      _hMuon1MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuon2MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuon1Muon2OSLS[NpdfID]->Fill(patMuon1->charge() * patMuon2->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hPZeta[NpdfID]->Fill(CalculatePZeta((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hPZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta2D[NpdfID]->Fill(CalculatePZetaVis((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))),CalculatePZeta((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta1D[NpdfID]->Fill((_PZetaCutCoefficient * CalculatePZeta((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin())))) + (_PZetaVisCutCoefficient * CalculatePZetaVis((*patMuon1),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons1 - 1),(*patMuon2),_SmearTheMuon,smearedMuonMomentumVector.at(theNumberOfMuons2 - 1),(*(_patMETs->begin())))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
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
		_hElectron1PtVsElectron2Pt[NpdfID]->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).pt(),smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hElectron1Electron2DeltaR[NpdfID]->Fill(reco::deltaR(smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1), smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hElectron1Electron2DeltaPtDivSumPt[NpdfID]->Fill((smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).pt() - smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).pt()) / (smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).pt() + smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hElectron1Electron2CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).phi() - smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		reco::Candidate::LorentzVector Temp_LorentzVect(
								(*(_patMETs->begin())).px() + patElectron1->px() - smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).px() + patElectron2->px() - smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).px() + deltaForMEx,
								(*(_patMETs->begin())).py() + patElectron1->py() - smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).py() + patElectron2->py() - smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).py() + deltaForMEy,
								(*(_patMETs->begin())).pz(),
								TMath::Sqrt( (((*(_patMETs->begin())).px() + deltaForMEx + patElectron1->px() - smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).px() + patElectron2->px() - smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).px()) *
									      ((*(_patMETs->begin())).px() + deltaForMEx + patElectron1->px() - smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).px() + patElectron2->px() - smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).px())) +
									     (((*(_patMETs->begin())).py() + deltaForMEy + patElectron1->py() - smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).py() + patElectron2->py() - smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).py()) *
									      ((*(_patMETs->begin())).py() + deltaForMEy + patElectron1->py() - smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).py() + patElectron2->py() - smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).py())) ));
		_hElectron1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).phi() - Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hElectron2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).phi() - Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hElectron1MetDeltaPhiVsElectron1Electron2CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).phi() - Temp_LorentzVect.phi())), cos(TMath::Abs(normalizedPhi(smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1).phi() - smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      } else {
		_hElectron1PtVsElectron2Pt[NpdfID]->Fill(patElectron1->pt(),patElectron2->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hElectron1Electron2DeltaR[NpdfID]->Fill(reco::deltaR(patElectron1->p4(), patElectron2->p4()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hElectron1Electron2DeltaPtDivSumPt[NpdfID]->Fill((patElectron1->pt() - patElectron2->pt()) / (patElectron1->pt() + patElectron2->pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hElectron1Electron2CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(patElectron1->phi() - patElectron2->phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hElectron1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patElectron1->phi() - (*(_patMETs->begin())).phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hElectron2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patElectron2->phi() - (*(_patMETs->begin())).phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hElectron1MetDeltaPhiVsElectron1Electron2CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patElectron1->phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(patElectron1->phi() - patElectron2->phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      }
	      if(CalculateThe4Momentum((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))).first) {_hReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      else {_hNotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      _hElectron1MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectron2MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectron1Electron2OSLS[NpdfID]->Fill(patElectron1->charge() * patElectron2->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hPZeta[NpdfID]->Fill(CalculatePZeta((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hPZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta2D[NpdfID]->Fill(CalculatePZetaVis((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))),CalculatePZeta((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta1D[NpdfID]->Fill((_PZetaCutCoefficient * CalculatePZeta((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin())))) + (_PZetaVisCutCoefficient * CalculatePZetaVis((*patElectron1),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1),(*patElectron2),_SmearTheElectron,smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1),(*(_patMETs->begin())))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    }
	  }
	}
      }
      if( ((_AnalyzeTauForLeg1) && (_AnalyzeTauForLeg2)) ) {
	int theNumberOfTaus1 = 0;
	for ( pat::TauCollection::const_iterator patTau1 = _patTaus->begin();patTau1 != _patTaus->end(); ++patTau1 ) {
	  theNumberOfTaus1++;
	  int theNumberOfTaus2 = 0;
	  for ( pat::TauCollection::const_iterator patTau2 = _patTaus->begin();patTau2 != _patTaus->end(); ++patTau2 ) {
	    theNumberOfTaus2++;
	    if ((passRecoTauCuts((*patTau1),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus1 - 1))) &&
		(passRecoTauCuts((*patTau2),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus2 - 1))) &&
		(passTopologyCuts((*patTau1),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus1 - 1),(*patTau2),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus2 - 1),(*(_patMETs->begin()))))) {
	      if(_SmearTheTau) {
		_hTau1PtVsTau2Pt[NpdfID]->Fill(smearedTauMomentumVector.at(theNumberOfTaus1 - 1).pt(),smearedTauMomentumVector.at(theNumberOfTaus2 - 1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hTau1Tau2DeltaR[NpdfID]->Fill(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus1 - 1), smearedTauMomentumVector.at(theNumberOfTaus2 - 1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hTau1Tau2DeltaPtDivSumPt[NpdfID]->Fill((smearedTauMomentumVector.at(theNumberOfTaus1 - 1).pt() - smearedTauMomentumVector.at(theNumberOfTaus2 - 1).pt()) / (smearedTauMomentumVector.at(theNumberOfTaus1 - 1).pt() + smearedTauMomentumVector.at(theNumberOfTaus2 - 1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		if(_UseTauSeedTrackForDiTauDiscrByOSLS) {
		  if (patTau1->isCaloTau()) {
		    if( (patTau1->leadTrack().isNonnull()) && (patTau2->leadTrack().isNonnull()) ) {
		      _hTau1Tau2OSLS[NpdfID]->Fill(patTau1->leadTrack()->charge() * patTau2->leadTrack()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		    }
		  } else {
		    if( (patTau1->leadPFChargedHadrCand().isNonnull()) && (patTau2->leadPFChargedHadrCand().isNonnull()) ) {
		      _hTau1Tau2OSLS[NpdfID]->Fill(patTau1->leadPFChargedHadrCand()->charge() * patTau2->leadPFChargedHadrCand()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		    }
		  }
		} else {_hTau1Tau2OSLS[NpdfID]->Fill(patTau1->charge() * patTau2->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
		reco::Candidate::LorentzVector Temp_LorentzVect(
								(*(_patMETs->begin())).px() + patTau1->px() - smearedTauMomentumVector.at(theNumberOfTaus1 - 1).px() + patTau2->px() - smearedTauMomentumVector.at(theNumberOfTaus2 - 1).px() + deltaForMEx,
								(*(_patMETs->begin())).py() + patTau1->py() - smearedTauMomentumVector.at(theNumberOfTaus1 - 1).py() + patTau2->py() - smearedTauMomentumVector.at(theNumberOfTaus2 - 1).py() + deltaForMEy,
								(*(_patMETs->begin())).pz(),
								TMath::Sqrt( (((*(_patMETs->begin())).px() + deltaForMEx + patTau1->px() - smearedTauMomentumVector.at(theNumberOfTaus1 - 1).px() + patTau2->px() - smearedTauMomentumVector.at(theNumberOfTaus2 - 1).px()) *
									      ((*(_patMETs->begin())).px() + deltaForMEx + patTau1->px() - smearedTauMomentumVector.at(theNumberOfTaus1 - 1).px() + patTau2->px() - smearedTauMomentumVector.at(theNumberOfTaus2 - 1).px())) +
									     (((*(_patMETs->begin())).py() + deltaForMEy + patTau1->py() - smearedTauMomentumVector.at(theNumberOfTaus1 - 1).py() + patTau2->py() - smearedTauMomentumVector.at(theNumberOfTaus2 - 1).py()) *
									      ((*(_patMETs->begin())).py() + deltaForMEy + patTau1->py() - smearedTauMomentumVector.at(theNumberOfTaus1 - 1).py() + patTau2->py() - smearedTauMomentumVector.at(theNumberOfTaus2 - 1).py())) ));
		_hTau1Tau2CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedTauMomentumVector.at(theNumberOfTaus1 - 1).phi() - smearedTauMomentumVector.at(theNumberOfTaus2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hTau1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauMomentumVector.at(theNumberOfTaus1 - 1).phi() -  Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hTau2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauMomentumVector.at(theNumberOfTaus2 - 1).phi() -  Temp_LorentzVect.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hTau1MetDeltaPhiVsTau1Tau2CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauMomentumVector.at(theNumberOfTaus1 - 1).phi() -  Temp_LorentzVect.phi())), cos(TMath::Abs(normalizedPhi(smearedTauMomentumVector.at(theNumberOfTaus1 - 1).phi() - smearedTauMomentumVector.at(theNumberOfTaus2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      } else {
		_hTau1PtVsTau2Pt[NpdfID]->Fill(patTau1->pt(),patTau2->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hTau1Tau2DeltaR[NpdfID]->Fill(reco::deltaR(patTau1->p4(), patTau2->p4()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hTau1Tau2DeltaPtDivSumPt[NpdfID]->Fill((patTau1->pt() - patTau2->pt()) / (patTau1->pt() + patTau2->pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		if(_UseTauSeedTrackForDiTauDiscrByOSLS) {
		  if (patTau1->isCaloTau()) {
		    if( (patTau1->leadTrack().isNonnull()) && (patTau2->leadTrack().isNonnull()) ) {
		      _hTau1Tau2OSLS[NpdfID]->Fill(patTau1->leadTrack()->charge() * patTau2->leadTrack()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		    }
		  } else {
		    if( (patTau1->leadPFChargedHadrCand().isNonnull()) && (patTau2->leadPFChargedHadrCand().isNonnull()) ) {
		      _hTau1Tau2OSLS[NpdfID]->Fill(patTau1->leadPFChargedHadrCand()->charge() * patTau2->leadPFChargedHadrCand()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		    }
		  }
		} else {_hTau1Tau2OSLS[NpdfID]->Fill(patTau1->charge() * patTau2->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
		_hTau1Tau2CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(patTau1->phi() - patTau2->phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hTau1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patTau1->phi() - (*(_patMETs->begin())).phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hTau2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patTau2->phi() - (*(_patMETs->begin())).phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		_hTau1MetDeltaPhiVsTau1Tau2CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(patTau1->phi() - (*(_patMETs->begin())).phi())), cos(TMath::Abs(normalizedPhi(patTau1->phi() - patTau2->phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      }
	      if(CalculateThe4Momentum((*patTau1),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus1 - 1),(*patTau2),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus2 - 1),(*(_patMETs->begin()))).first) {_hReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau1),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus1 - 1),(*patTau2),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus2 - 1),(*(_patMETs->begin()))).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      else {_hNotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau1),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus1 - 1),(*patTau2),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus2 - 1),(*(_patMETs->begin()))).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      _hTau1MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau1),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus1 - 1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTau2MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau2),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus2 - 1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hPZeta[NpdfID]->Fill(CalculatePZeta((*patTau1),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus1 - 1),(*patTau2),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus2 - 1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hPZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patTau1),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus1 - 1),(*patTau2),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus2 - 1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta2D[NpdfID]->Fill(CalculatePZetaVis((*patTau1),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus1 - 1),(*patTau2),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus2 - 1),(*(_patMETs->begin()))),CalculatePZeta((*patTau1),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus1 - 1),(*patTau2),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus2 - 1),(*(_patMETs->begin()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta1D[NpdfID]->Fill((_PZetaCutCoefficient * CalculatePZeta((*patTau1),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus1 - 1),(*patTau2),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus2 - 1),(*(_patMETs->begin())))) + (_PZetaVisCutCoefficient * CalculatePZetaVis((*patTau1),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus1 - 1),(*patTau2),_SmearTheTau,smearedTauMomentumVector.at(theNumberOfTaus2 - 1),(*(_patMETs->begin())))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    }
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
//  if (_DoProduceNtuple) { iEvent.getByLabel(_RecoDiTauSource, _patDiTaus); }
  iEvent.getByLabel(_RecoTauSource, _patTaus);
  iEvent.getByLabel(_RecoMuonSource, _patMuons);
  if(_GenParticleSource.label() != "") { iEvent.getByLabel(_GenParticleSource, _genParticles); }
  iEvent.getByLabel(_RecoElectronSource, _patElectrons);
  iEvent.getByLabel(_RecoJetSource, _patJets);
  iEvent.getByLabel(_RecoMetSource, _patMETs);
  iEvent.getByLabel(_RecoVertexSource, _primaryEventVertexCollection);
  iEvent.getByLabel(_RecoTriggerSource, _triggerResults);
/*
  iEvent.getByLabel("muons", _recoMuonsForMetCorrections);
  iEvent.getByLabel("muonMETValueMapProducer", "muCorrData", vm_muCorrData_h);
  iEvent.getByLabel("selectedLayer1Muons", _patMuonsForMetCorrections);
*/
}

pair<bool, pair<float, float> > HiMassTauAnalysis::isZee(reco::Candidate::LorentzVector theObject) {
  pair<bool, pair<float, float> > theOutPutPair;
  bool eventIsZee = false;
  bool massWindow = false;
  bool ptAsymmWindow = false;
  const float zeeMass = 90.1876;
  const float zeeWidht = 2.4952;
  float zeePtAsymmetry = -10.;
  unsigned int i = 0;
  pair<float, float> theMassPtAsymmPair;
  // if mass is within 3 sigmas of z or pt asymmetry is small set to true.				     
  for(pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron){
    reco::Candidate::LorentzVector Electron;
    
    if(_SmearTheElectron){Electron = smearedElectronMomentumVector.at(i); i++;}
    else{Electron = patElectron->p4();}
   
    if(reco::deltaR(theObject, Electron) < _DiTauDeltaRCut) continue;
    if(theObject == Electron)continue;
    reco::Candidate::LorentzVector The_LorentzVect = theObject + Electron;
    zeePtAsymmetry = (theObject.pt() - Electron.pt())/(theObject.pt() + Electron.pt());
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
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::matchToGen(const pat::Electron& theObject) {
  bool isGenMatched = false;
  reco::Candidate::LorentzVector theGenObject(0,0,0,0);
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
  pair<bool, reco::Candidate::LorentzVector> GenMatchedInformation(isGenMatched,theGenObject);
  return GenMatchedInformation;
}

//-----Matching to generator level objects
//-----Muons
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::matchToGen(const pat::Muon& theObject) {
  bool isGenMatched = false;
  reco::Candidate::LorentzVector theGenObject(0,0,0,0);
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
  pair<bool, reco::Candidate::LorentzVector> GenMatchedInformation(isGenMatched,theGenObject);
  return GenMatchedInformation;
}

//-----Matching to generator level objects
//-----Taus
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::matchToGen(const pat::Tau& theObject) {
  bool isGenMatched = false;
  reco::Candidate::LorentzVector theGenObject(0,0,0,0);
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
  pair<bool, reco::Candidate::LorentzVector> GenMatchedInformation(isGenMatched,theGenObject);
  return GenMatchedInformation;
}

//-----Calculate zeta variables
double HiMassTauAnalysis::CalculatePZeta(const pat::Tau& patTau, bool smearT, reco::Candidate::LorentzVector smearedLVT, const pat::Muon& patMuon, bool smearM, reco::Candidate::LorentzVector smearedLVM, const pat::MET& patMET) {
  if(smearT) {
    if(smearM) {
      double zetaX = cos(smearedLVT.phi()) + cos(smearedLVM.phi());
      double zetaY = sin(smearedLVT.phi()) + sin(smearedLVM.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = smearedLVT.px() + smearedLVM.px();
      double visPy = smearedLVT.py() + smearedLVM.py();
      double px = visPx + (patMET.px() + deltaForMEx + patMuon.px() - smearedLVM.px() + patTau.px() - smearedLVT.px());
      double py = visPy + (patMET.py() + deltaForMEy + patMuon.py() - smearedLVM.py() + patTau.py() - smearedLVT.py());
      double pZeta = px*zetaX + py*zetaY;
      return pZeta;
    } else {
      double zetaX = cos(smearedLVT.phi()) + cos(patMuon.phi());
      double zetaY = sin(smearedLVT.phi()) + sin(patMuon.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = smearedLVT.px() + patMuon.px();
      double visPy = smearedLVT.py() + patMuon.py();
      double px = visPx + (patMET.px() + deltaForMEx + patTau.px() - smearedLVT.px());
      double py = visPy + (patMET.py() + deltaForMEy + patTau.py() - smearedLVT.py());
      double pZeta = px*zetaX + py*zetaY;
      return pZeta;
    }
  } else {
    if(smearM) {
      double zetaX = cos(patTau.phi()) + cos(smearedLVM.phi());
      double zetaY = sin(patTau.phi()) + sin(smearedLVM.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = patTau.px() + smearedLVM.px();
      double visPy = patTau.py() + smearedLVM.py();
      double px = visPx + (patMET.px() + deltaForMEx + patMuon.px() - smearedLVM.px());
      double py = visPy + (patMET.py() + deltaForMEy + patMuon.py() - smearedLVM.py());
      double pZeta = px*zetaX + py*zetaY;
      return pZeta;
    } else {
      double zetaX = cos(patTau.phi()) + cos(patMuon.phi());
      double zetaY = sin(patTau.phi()) + sin(patMuon.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = patTau.px() + patMuon.px();
      double visPy = patTau.py() + patMuon.py();
      double px = visPx + patMET.px() + deltaForMEx;
      double py = visPy + patMET.py() + deltaForMEy;
      double pZeta = px*zetaX + py*zetaY;
      return pZeta;
    }
  }
}
double HiMassTauAnalysis::CalculatePZeta(const pat::Tau& patTau, bool smearT, reco::Candidate::LorentzVector smearedLVT, const pat::Electron& patElectron, bool smearM, reco::Candidate::LorentzVector smearedLVM, const pat::MET& patMET) {
  if(smearT) {
    if(smearM) {
      double zetaX = cos(smearedLVT.phi()) + cos(smearedLVM.phi());
      double zetaY = sin(smearedLVT.phi()) + sin(smearedLVM.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = smearedLVT.px() + smearedLVM.px();
      double visPy = smearedLVT.py() + smearedLVM.py();
      double px = visPx + (patMET.px() + deltaForMEx + patElectron.px() - smearedLVM.px() + patTau.px() - smearedLVT.px());
      double py = visPy + (patMET.py() + deltaForMEy + patElectron.py() - smearedLVM.py() + patTau.py() - smearedLVT.py());
      double pZeta = px*zetaX + py*zetaY;
      return pZeta;
    } else {
      double zetaX = cos(smearedLVT.phi()) + cos(patElectron.phi());
      double zetaY = sin(smearedLVT.phi()) + sin(patElectron.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = smearedLVT.px() + patElectron.px();
      double visPy = smearedLVT.py() + patElectron.py();
      double px = visPx + (patMET.px() + deltaForMEx + patTau.px() - smearedLVT.px());
      double py = visPy + (patMET.py() + deltaForMEy + patTau.py() - smearedLVT.py());
      double pZeta = px*zetaX + py*zetaY;
      return pZeta;
    }
  } else {
    if(smearM) {
      double zetaX = cos(patTau.phi()) + cos(smearedLVM.phi());
      double zetaY = sin(patTau.phi()) + sin(smearedLVM.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = patTau.px() + smearedLVM.px();
      double visPy = patTau.py() + smearedLVM.py();
      double px = visPx + (patMET.px() + deltaForMEx + patElectron.px() - smearedLVM.px());
      double py = visPy + (patMET.py() + deltaForMEy + patElectron.py() - smearedLVM.py());
      double pZeta = px*zetaX + py*zetaY;
      return pZeta;
    } else {
      double zetaX = cos(patTau.phi()) + cos(patElectron.phi());
      double zetaY = sin(patTau.phi()) + sin(patElectron.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = patTau.px() + patElectron.px();
      double visPy = patTau.py() + patElectron.py();
      double px = visPx + patMET.px() + deltaForMEx;
      double py = visPy + patMET.py() + deltaForMEy;
      double pZeta = px*zetaX + py*zetaY;
      return pZeta;
    }
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
      double px = visPx + (patMET.px() + deltaForMEx + patMuon.px() - smearedLVM.px() + patElectron.px() - smearedLVE.px());
      double py = visPy + (patMET.py() + deltaForMEy + patMuon.py() - smearedLVM.py() + patElectron.py() - smearedLVE.py());
      double pZeta = px*zetaX + py*zetaY;
      return pZeta;
    } else {
      double zetaX = cos(patMuon.phi()) + cos(smearedLVE.phi());
      double zetaY = sin(patMuon.phi()) + sin(smearedLVE.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = patMuon.px() + smearedLVE.px();
      double visPy = patMuon.py() + smearedLVE.py();
      double px = visPx + (patMET.px() + deltaForMEx + patElectron.px() - smearedLVE.px());
      double py = visPy + (patMET.py() + deltaForMEy + patElectron.py() - smearedLVE.py());
      double pZeta = px*zetaX + py*zetaY;
      return pZeta;
    }
  } else {
    if(smearM) {
      double zetaX = cos(smearedLVM.phi()) + cos(patElectron.phi());
      double zetaY = sin(smearedLVM.phi()) + sin(patElectron.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = smearedLVM.px() + patElectron.px();
      double visPy = smearedLVM.py() + patElectron.py();
      double px = visPx + (patMET.px() + deltaForMEx + patMuon.px() - smearedLVM.px());
      double py = visPy + (patMET.py() + deltaForMEy + patMuon.py() - smearedLVM.py());
      double pZeta = px*zetaX + py*zetaY;
      return pZeta;
    } else {
      double zetaX = cos(patMuon.phi()) + cos(patElectron.phi());
      double zetaY = sin(patMuon.phi()) + sin(patElectron.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = patMuon.px() + patElectron.px();
      double visPy = patMuon.py() + patElectron.py();
      double px = visPx + patMET.px() + deltaForMEx;
      double py = visPy + patMET.py() + deltaForMEy;
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
    double px = visPx + (patMET.px() + deltaForMEx + patMuon1.px() - smearedLV1.px() + patMuon2.px() - smearedLV2.px());
    double py = visPy + (patMET.py() + deltaForMEy + patMuon1.py() - smearedLV1.py() + patMuon2.py() - smearedLV2.py());
    double pZeta = px*zetaX + py*zetaY;
    return pZeta;
  } else {
    double zetaX = cos(patMuon1.phi()) + cos(patMuon2.phi());
    double zetaY = sin(patMuon1.phi()) + sin(patMuon2.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = patMuon1.px() + patMuon2.px();
    double visPy = patMuon1.py() + patMuon2.py();
    double px = visPx + patMET.px() + deltaForMEx;
    double py = visPy + patMET.py() + deltaForMEy;
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
    double px = visPx + (patMET.px() + deltaForMEx + patElectron1.px() - smearedLV1.px() + patElectron2.px() - smearedLV2.px());
    double py = visPy + (patMET.py() + deltaForMEy + patElectron1.py() - smearedLV1.py() + patElectron2.py() - smearedLV2.py());
    double pZeta = px*zetaX + py*zetaY;
    return pZeta;
  } else {
    double zetaX = cos(patElectron1.phi()) + cos(patElectron2.phi());
    double zetaY = sin(patElectron1.phi()) + sin(patElectron2.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = patElectron1.px() + patElectron2.px();
    double visPy = patElectron1.py() + patElectron2.py();
    double px = visPx + patMET.px() + deltaForMEx;
    double py = visPy + patMET.py() + deltaForMEy;
    double pZeta = px*zetaX + py*zetaY;
    return pZeta;
  }
}
double HiMassTauAnalysis::CalculatePZeta(const pat::Tau& patTau1,  bool smear1, reco::Candidate::LorentzVector smearedLV1, const pat::Tau& patTau2,  bool smear2, reco::Candidate::LorentzVector smearedLV2, const pat::MET& patMET) {
  if(smear1) {
    double zetaX = cos(smearedLV1.phi()) + cos(smearedLV2.phi());
    double zetaY = sin(smearedLV1.phi()) + sin(smearedLV2.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = smearedLV1.px() + smearedLV2.px();
    double visPy = smearedLV1.py() + smearedLV2.py();
    double px = visPx + (patMET.px() + deltaForMEx + patTau1.px() - smearedLV1.px() + patTau2.px() - smearedLV2.px());
    double py = visPy + (patMET.py() + deltaForMEy + patTau1.py() - smearedLV1.py() + patTau2.py() - smearedLV2.py());
    double pZeta = px*zetaX + py*zetaY;
    return pZeta;
  } else {
    double zetaX = cos(patTau1.phi()) + cos(patTau2.phi());
    double zetaY = sin(patTau1.phi()) + sin(patTau2.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = patTau1.px() + patTau2.px();
    double visPy = patTau1.py() + patTau2.py();
    double px = visPx + patMET.px() + deltaForMEx;
    double py = visPy + patMET.py() + deltaForMEy;
    double pZeta = px*zetaX + py*zetaY;
    return pZeta;
  }
}
double HiMassTauAnalysis::CalculatePZetaVis(const pat::Tau& patTau, bool smearT, reco::Candidate::LorentzVector smearedLVT, const pat::Muon& patMuon, bool smearM, reco::Candidate::LorentzVector smearedLVM, const pat::MET& patMET) {
  if(smearT) {
    if(smearM) {
      double zetaX = cos(smearedLVT.phi()) + cos(smearedLVM.phi());
      double zetaY = sin(smearedLVT.phi()) + sin(smearedLVM.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = smearedLVT.px() + smearedLVM.px();
      double visPy = smearedLVT.py() + smearedLVM.py();
      double pZetaVis = visPx*zetaX + visPy*zetaY;
      return pZetaVis;
    } else {
      double zetaX = cos(smearedLVT.phi()) + cos(patMuon.phi());
      double zetaY = sin(smearedLVT.phi()) + sin(patMuon.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = smearedLVT.px() + patMuon.px();
      double visPy = smearedLVT.py() + patMuon.py();
      double pZetaVis = visPx*zetaX + visPy*zetaY;
      return pZetaVis;
    }
  } else {
    if(smearM) {
      double zetaX = cos(patTau.phi()) + cos(smearedLVM.phi());
      double zetaY = sin(patTau.phi()) + sin(smearedLVM.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = patTau.px() + smearedLVM.px();
      double visPy = patTau.py() + smearedLVM.py();
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
}
double HiMassTauAnalysis::CalculatePZetaVis(const pat::Tau& patTau, bool smearT, reco::Candidate::LorentzVector smearedLVT, const pat::Electron& patElectron, bool smearM, reco::Candidate::LorentzVector smearedLVM, const pat::MET& patMET) {
  if(smearT) {
    if(smearM) {
      double zetaX = cos(smearedLVT.phi()) + cos(smearedLVM.phi());
      double zetaY = sin(smearedLVT.phi()) + sin(smearedLVM.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = smearedLVT.px() + smearedLVM.px();
      double visPy = smearedLVT.py() + smearedLVM.py();
      double pZetaVis = visPx*zetaX + visPy*zetaY;
      return pZetaVis;
    } else {
      double zetaX = cos(smearedLVT.phi()) + cos(patElectron.phi());
      double zetaY = sin(smearedLVT.phi()) + sin(patElectron.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = smearedLVT.px() + patElectron.px();
      double visPy = smearedLVT.py() + patElectron.py();
      double pZetaVis = visPx*zetaX + visPy*zetaY;
      return pZetaVis;
    }
  } else {
    if(smearM) {
      double zetaX = cos(patTau.phi()) + cos(smearedLVM.phi());
      double zetaY = sin(patTau.phi()) + sin(smearedLVM.phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
      double visPx = patTau.px() + smearedLVM.px();
      double visPy = patTau.py() + smearedLVM.py();
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
double HiMassTauAnalysis::CalculatePZetaVis(const pat::Tau& patTau1, bool smear1, reco::Candidate::LorentzVector smearedLV1, const pat::Tau& patTau2, bool smear2, reco::Candidate::LorentzVector smearedLV2, const pat::MET& patMET) {
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
    double zetaX = cos(patTau1.phi()) + cos(patTau2.phi());
    double zetaY = sin(patTau1.phi()) + sin(patTau2.phi());
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
    double visPx = patTau1.px() + patTau2.px();
    double visPy = patTau1.py() + patTau2.py();
    double pZetaVis = visPx*zetaX + visPy*zetaY;
    return pZetaVis;
  }
}

//-----Calculate mass reco variables
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::CalculateThe4Momentum(const pat::Tau& patTau, bool smearT, reco::Candidate::LorentzVector smearedLVT, const pat::Electron& patElectron, bool smearE, reco::Candidate::LorentzVector smearedLVE, const pat::MET& patMET) {
  if(smearT) {
    if(smearE) {
      double MEpx = patMET.px() + deltaForMEx + patElectron.px() - smearedLVE.px() + patTau.px() - smearedLVT.px();
      double MEpy = patMET.py() + deltaForMEy + patElectron.py() - smearedLVE.py() + patTau.py() - smearedLVT.py();
      if(_UseVectorSumOfVisProductsAndMetMassReco) {
        double px = smearedLVT.px() + smearedLVE.px() + MEpx;
        double py = smearedLVT.py() + smearedLVE.py() + MEpy;
        double pz = smearedLVT.pz() + smearedLVE.pz();
        double e = smearedLVT.energy() + smearedLVE.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else if(_UseCollinerApproxMassReco) {
        double x1_numerator = (smearedLVT.px() * smearedLVE.py()) - (smearedLVE.px() * smearedLVT.py());
        double x1_denominator = (smearedLVE.py() * (smearedLVT.px() + MEpx)) - (smearedLVE.px() * (smearedLVT.py() + MEpy));
        double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
        double x2_numerator = x1_numerator;
        double x2_denominator = (smearedLVT.px() * (smearedLVE.py() + MEpy)) - (smearedLVT.py() * (smearedLVE.px() + MEpx));
        double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
        if ( (x1 > 0. && x1 < 1.) &&
             (x2 > 0. && x2 < 1.) ) {
          reco::Candidate::LorentzVector The_LorentzVect = (smearedLVT / x1) + (smearedLVE / x2);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
          return MassRecoInformation;
        } else {
          double px = smearedLVT.px() + smearedLVE.px() + MEpx;
          double py = smearedLVT.py() + smearedLVE.py() + MEpy;
          double pz = smearedLVT.pz() + smearedLVE.pz();
          double e = smearedLVT.energy() + smearedLVE.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
          reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
          return MassRecoInformation;
        }
      } else {
        reco::Candidate::LorentzVector The_LorentzVect = smearedLVT + smearedLVE;
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      }
    } else {
      double MEpx = patMET.px() + deltaForMEx + patTau.px() - smearedLVT.px();
      double MEpy = patMET.py() + deltaForMEy + patTau.py() - smearedLVT.py();
      if(_UseVectorSumOfVisProductsAndMetMassReco) {
        double px = smearedLVT.px() + patElectron.px() + MEpx;
        double py = smearedLVT.py() + patElectron.py() + MEpy;
        double pz = smearedLVT.pz() + patElectron.pz();
        double e = smearedLVT.energy() + patElectron.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else if(_UseCollinerApproxMassReco) {
        double x1_numerator = (smearedLVT.px() * patElectron.py()) - (patElectron.px() * smearedLVT.py());
        double x1_denominator = (patElectron.py() * (smearedLVT.px() + MEpx)) - (patElectron.px() * (smearedLVT.py() + MEpy));
        double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
        double x2_numerator = x1_numerator;
        double x2_denominator = (smearedLVT.px() * (patElectron.py() + MEpy)) - (smearedLVT.py() * (patElectron.px() + MEpx));
        double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
        if ( (x1 > 0. && x1 < 1.) &&
             (x2 > 0. && x2 < 1.) ) {
          reco::Candidate::LorentzVector The_LorentzVect = (smearedLVT / x1) + (patElectron.p4() / x2);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
          return MassRecoInformation;
        } else {
          double px = smearedLVT.px() + patElectron.px() + MEpx;
          double py = smearedLVT.py() + patElectron.py() + MEpy;
          double pz = smearedLVT.pz() + patElectron.pz();
          double e = smearedLVT.energy() + patElectron.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
          reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
          return MassRecoInformation;
        }
      } else {
        reco::Candidate::LorentzVector The_LorentzVect = smearedLVT + patElectron.p4();
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      }
    }
  } else {
    if(smearE) {
      double MEpx = patMET.px() + deltaForMEx + patElectron.px() - smearedLVE.px();
      double MEpy = patMET.py() + deltaForMEy + patElectron.py() - smearedLVE.py();
      if(_UseVectorSumOfVisProductsAndMetMassReco) {
        double px = patTau.px() + smearedLVE.px() + MEpx;
        double py = patTau.py() + smearedLVE.py() + MEpy;
        double pz = patTau.pz() + smearedLVE.pz();
        double e = patTau.energy() + smearedLVE.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else if(_UseCollinerApproxMassReco) {
        double x1_numerator = (patTau.px() * smearedLVE.py()) - (smearedLVE.px() * patTau.py());
        double x1_denominator = (smearedLVE.py() * (patTau.px() + MEpx)) - (smearedLVE.px() * (patTau.py() + MEpy));
        double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
        double x2_numerator = x1_numerator;
        double x2_denominator = (patTau.px() * (smearedLVE.py() + MEpy)) - (patTau.py() * (smearedLVE.px() + MEpx));
        double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
        if ( (x1 > 0. && x1 < 1.) &&
             (x2 > 0. && x2 < 1.) ) {
          reco::Candidate::LorentzVector The_LorentzVect = (patTau.p4() / x1) + (smearedLVE / x2);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
          return MassRecoInformation;
        } else {
          double px = patTau.px() + smearedLVE.px() + MEpx;
          double py = patTau.py() + smearedLVE.py() + MEpy;
          double pz = patTau.pz() + smearedLVE.pz();
          double e = patTau.energy() + smearedLVE.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
          reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
          return MassRecoInformation;
        }
      } else {
        reco::Candidate::LorentzVector The_LorentzVect = patTau.p4() + smearedLVE;
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      }
    } else {
      if(_UseVectorSumOfVisProductsAndMetMassReco) {
        double px = patTau.px() + patElectron.px() + patMET.px() + deltaForMEx;
        double py = patTau.py() + patElectron.py() + patMET.py() + deltaForMEy;
        double pz = patTau.pz() + patElectron.pz();
        double e = patTau.energy() + patElectron.energy() + TMath::Sqrt(((patMET.px() + deltaForMEx) * (patMET.px() + deltaForMEx)) + ((patMET.py() + deltaForMEy) * (patMET.py() + deltaForMEy)));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else if(_UseCollinerApproxMassReco) {
        double x1_numerator = (patTau.px() * patElectron.py()) - (patElectron.px() * patTau.py());
        double x1_denominator = (patElectron.py() * (patTau.px() + patMET.px() + deltaForMEx)) - (patElectron.px() * (patTau.py() + patMET.py() + deltaForMEy));
        double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
        double x2_numerator = x1_numerator;
        double x2_denominator = (patTau.px() * (patElectron.py() + patMET.py() + deltaForMEy)) - (patTau.py() * (patElectron.px() + patMET.px() + deltaForMEx));
        double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
        if ( (x1 > 0. && x1 < 1.) &&
             (x2 > 0. && x2 < 1.) ) {
          reco::Candidate::LorentzVector The_LorentzVect = (patTau.p4() / x1) + (patElectron.p4() / x2);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
          return MassRecoInformation;
        } else {
          double px = patTau.px() + patElectron.px() + patMET.px() + deltaForMEx;
          double py = patTau.py() + patElectron.py() + patMET.py() + deltaForMEy;
          double pz = patTau.pz() + patElectron.pz();
          double e = patTau.energy() + patElectron.energy() + TMath::Sqrt(((patMET.px() + deltaForMEx) * (patMET.px() + deltaForMEx)) + ((patMET.py() + deltaForMEy) * (patMET.py() + deltaForMEy)));
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
}
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::CalculateThe4Momentum(const pat::Tau& patTau, bool smearT, reco::Candidate::LorentzVector smearedLVT, const pat::Muon& patMuon, bool smearE, reco::Candidate::LorentzVector smearedLVE, const pat::MET& patMET) {
  if(smearT) {
    if(smearE) {
      double MEpx = patMET.px() + deltaForMEx + patMuon.px() - smearedLVE.px() + patTau.px() - smearedLVT.px();
      double MEpy = patMET.py() + deltaForMEy + patMuon.py() - smearedLVE.py() + patTau.py() - smearedLVT.py();
      if(_UseVectorSumOfVisProductsAndMetMassReco) {
        double px = smearedLVT.px() + smearedLVE.px() + MEpx;
        double py = smearedLVT.py() + smearedLVE.py() + MEpy;
        double pz = smearedLVT.pz() + smearedLVE.pz();
        double e = smearedLVT.energy() + smearedLVE.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else if(_UseCollinerApproxMassReco) {
        double x1_numerator = (smearedLVT.px() * smearedLVE.py()) - (smearedLVE.px() * smearedLVT.py());
        double x1_denominator = (smearedLVE.py() * (smearedLVT.px() + MEpx)) - (smearedLVE.px() * (smearedLVT.py() + MEpy));
        double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
        double x2_numerator = x1_numerator;
        double x2_denominator = (smearedLVT.px() * (smearedLVE.py() + MEpy)) - (smearedLVT.py() * (smearedLVE.px() + MEpx));
        double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
        if ( (x1 > 0. && x1 < 1.) &&
             (x2 > 0. && x2 < 1.) ) {
          reco::Candidate::LorentzVector The_LorentzVect = (smearedLVT / x1) + (smearedLVE / x2);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
          return MassRecoInformation;
        } else {
          double px = smearedLVT.px() + smearedLVE.px() + MEpx;
          double py = smearedLVT.py() + smearedLVE.py() + MEpy;
          double pz = smearedLVT.pz() + smearedLVE.pz();
          double e = smearedLVT.energy() + smearedLVE.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
          reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
          return MassRecoInformation;
        }
      } else {
        reco::Candidate::LorentzVector The_LorentzVect = smearedLVT + smearedLVE;
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      }
    } else {
      double MEpx = patMET.px() + deltaForMEx + patTau.px() - smearedLVT.px();
      double MEpy = patMET.py() + deltaForMEy + patTau.py() - smearedLVT.py();
      if(_UseVectorSumOfVisProductsAndMetMassReco) {
        double px = smearedLVT.px() + patMuon.px() + MEpx;
        double py = smearedLVT.py() + patMuon.py() + MEpy;
        double pz = smearedLVT.pz() + patMuon.pz();
        double e = smearedLVT.energy() + patMuon.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else if(_UseCollinerApproxMassReco) {
        double x1_numerator = (smearedLVT.px() * patMuon.py()) - (patMuon.px() * smearedLVT.py());
        double x1_denominator = (patMuon.py() * (smearedLVT.px() + MEpx)) - (patMuon.px() * (smearedLVT.py() + MEpy));
        double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
        double x2_numerator = x1_numerator;
        double x2_denominator = (smearedLVT.px() * (patMuon.py() + MEpy)) - (smearedLVT.py() * (patMuon.px() + MEpx));
        double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
        if ( (x1 > 0. && x1 < 1.) &&
             (x2 > 0. && x2 < 1.) ) {
          reco::Candidate::LorentzVector The_LorentzVect = (smearedLVT / x1) + (patMuon.p4() / x2);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
          return MassRecoInformation;
        } else {
          double px = smearedLVT.px() + patMuon.px() + MEpx;
          double py = smearedLVT.py() + patMuon.py() + MEpy;
          double pz = smearedLVT.pz() + patMuon.pz();
          double e = smearedLVT.energy() + patMuon.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
          reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
          return MassRecoInformation;
        }
      } else {
        reco::Candidate::LorentzVector The_LorentzVect = smearedLVT + patMuon.p4();
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      }
    }
  } else {
    if(smearE) {
      double MEpx = patMET.px() + deltaForMEx + patMuon.px() - smearedLVE.px();
      double MEpy = patMET.py() + deltaForMEy + patMuon.py() - smearedLVE.py();
      if(_UseVectorSumOfVisProductsAndMetMassReco) {
        double px = patTau.px() + smearedLVE.px() + MEpx;
        double py = patTau.py() + smearedLVE.py() + MEpy;
        double pz = patTau.pz() + smearedLVE.pz();
        double e = patTau.energy() + smearedLVE.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else if(_UseCollinerApproxMassReco) {
        double x1_numerator = (patTau.px() * smearedLVE.py()) - (smearedLVE.px() * patTau.py());
        double x1_denominator = (smearedLVE.py() * (patTau.px() + MEpx)) - (smearedLVE.px() * (patTau.py() + MEpy));
        double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
        double x2_numerator = x1_numerator;
        double x2_denominator = (patTau.px() * (smearedLVE.py() + MEpy)) - (patTau.py() * (smearedLVE.px() + MEpx));
        double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
        if ( (x1 > 0. && x1 < 1.) &&
             (x2 > 0. && x2 < 1.) ) {
          reco::Candidate::LorentzVector The_LorentzVect = (patTau.p4() / x1) + (smearedLVE / x2);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
          return MassRecoInformation;
        } else {
          double px = patTau.px() + smearedLVE.px() + MEpx;
          double py = patTau.py() + smearedLVE.py() + MEpy;
          double pz = patTau.pz() + smearedLVE.pz();
          double e = patTau.energy() + smearedLVE.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
          reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
          return MassRecoInformation;
        }
      } else {
        reco::Candidate::LorentzVector The_LorentzVect = patTau.p4() + smearedLVE;
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      }
    } else {
      if(_UseVectorSumOfVisProductsAndMetMassReco) {
        double px = patTau.px() + patMuon.px() + patMET.px() + deltaForMEx;
        double py = patTau.py() + patMuon.py() + patMET.py() + deltaForMEy;
        double pz = patTau.pz() + patMuon.pz();
        double e = patTau.energy() + patMuon.energy() + TMath::Sqrt(((patMET.px() + deltaForMEx) * (patMET.px() + deltaForMEx)) + ((patMET.py() + deltaForMEy) * (patMET.py() + deltaForMEy)));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else if(_UseCollinerApproxMassReco) {
        double x1_numerator = (patTau.px() * patMuon.py()) - (patMuon.px() * patTau.py());
        double x1_denominator = (patMuon.py() * (patTau.px() + patMET.px() + deltaForMEx)) - (patMuon.px() * (patTau.py() + patMET.py() + deltaForMEy));
        double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
        double x2_numerator = x1_numerator;
        double x2_denominator = (patTau.px() * (patMuon.py() + patMET.py() + deltaForMEy)) - (patTau.py() * (patMuon.px() + patMET.px() + deltaForMEx));
        double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
        if ( (x1 > 0. && x1 < 1.) &&
             (x2 > 0. && x2 < 1.) ) {
          reco::Candidate::LorentzVector The_LorentzVect = (patTau.p4() / x1) + (patMuon.p4() / x2);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
          return MassRecoInformation;
        } else {
          double px = patTau.px() + patMuon.px() + patMET.px() + deltaForMEx;
          double py = patTau.py() + patMuon.py() + patMET.py() + deltaForMEy;
          double pz = patTau.pz() + patMuon.pz();
          double e = patTau.energy() + patMuon.energy() + TMath::Sqrt(((patMET.px() + deltaForMEx) * (patMET.px() + deltaForMEx)) + ((patMET.py() + deltaForMEy) * (patMET.py() + deltaForMEy)));
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
}
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::CalculateThe4Momentum(const pat::Electron& patElectron, bool smearE, reco::Candidate::LorentzVector smearedLVE, const pat::Muon& patMuon, bool smearM, reco::Candidate::LorentzVector smearedLVM, const pat::MET& patMET) {
  if(smearE) {
    if(smearM) {
      double MEpx = patMET.px() + deltaForMEx + patElectron.px() - smearedLVE.px() + patMuon.px() - smearedLVM.px();
      double MEpy = patMET.py() + deltaForMEy + patElectron.py() - smearedLVE.py() + patMuon.py() - smearedLVM.py();
      if(_UseVectorSumOfVisProductsAndMetMassReco) {
        double px = smearedLVE.px() + smearedLVM.px() + MEpx;
        double py = smearedLVE.py() + smearedLVM.py() + MEpy;
        double pz = smearedLVE.pz() + smearedLVM.pz();
        double e = smearedLVE.energy() + smearedLVM.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else if(_UseCollinerApproxMassReco) {
        double x1_numerator = (smearedLVE.px() * smearedLVM.py()) - (smearedLVM.px() * smearedLVE.py());
        double x1_denominator = (smearedLVM.py() * (smearedLVE.px() + MEpx)) - (smearedLVM.px() * (smearedLVE.py() + MEpy));
        double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
        double x2_numerator = x1_numerator;
        double x2_denominator = (smearedLVE.px() * (smearedLVM.py() + MEpy)) - (smearedLVE.py() * (smearedLVM.px() + MEpx));
        double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
        if ( (x1 > 0. && x1 < 1.) &&
             (x2 > 0. && x2 < 1.) ) {
          reco::Candidate::LorentzVector The_LorentzVect = (smearedLVE / x1) + (smearedLVM / x2);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
          return MassRecoInformation;
        } else {
          double px = smearedLVE.px() + smearedLVM.px() + MEpx;
          double py = smearedLVE.py() + smearedLVM.py() + MEpy;
          double pz = smearedLVE.pz() + smearedLVM.pz();
          double e = smearedLVE.energy() + smearedLVM.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
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
      double MEpx = patMET.px() + deltaForMEx + patElectron.px() - smearedLVE.px();
      double MEpy = patMET.py() + deltaForMEy + patElectron.py() - smearedLVE.py();
      if(_UseVectorSumOfVisProductsAndMetMassReco) {
        double px = smearedLVE.px() + patMuon.px() + MEpx;
        double py = smearedLVE.py() + patMuon.py() + MEpy;
        double pz = smearedLVE.pz() + patMuon.pz();
        double e = smearedLVE.energy() + patMuon.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else if(_UseCollinerApproxMassReco) {
        double x1_numerator = (smearedLVE.px() * patMuon.py()) - (patMuon.px() * smearedLVE.py());
        double x1_denominator = (patMuon.py() * (smearedLVE.px() + MEpx)) - (patMuon.px() * (smearedLVE.py() + MEpy));
        double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
        double x2_numerator = x1_numerator;
        double x2_denominator = (smearedLVE.px() * (patMuon.py() + MEpy)) - (smearedLVE.py() * (patMuon.px() + MEpx));
        double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
        if ( (x1 > 0. && x1 < 1.) &&
             (x2 > 0. && x2 < 1.) ) {
          reco::Candidate::LorentzVector The_LorentzVect = (smearedLVE / x1) + (patMuon.p4() / x2);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
          return MassRecoInformation;
        } else {
          double px = smearedLVE.px() + patMuon.px() + MEpx;
          double py = smearedLVE.py() + patMuon.py() + MEpy;
          double pz = smearedLVE.pz() + patMuon.pz();
          double e = smearedLVE.energy() + patMuon.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
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
      double MEpx = patMET.px() + deltaForMEx + patMuon.px() - smearedLVM.px();
      double MEpy = patMET.py() + deltaForMEy + patMuon.py() - smearedLVM.py();
      if(_UseVectorSumOfVisProductsAndMetMassReco) {
        double px = patElectron.px() + smearedLVM.px() + MEpx;
        double py = patElectron.py() + smearedLVM.py() + MEpy;
        double pz = patElectron.pz() + smearedLVM.pz();
        double e = patElectron.energy() + smearedLVM.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else if(_UseCollinerApproxMassReco) {
        double x1_numerator = (patElectron.px() * smearedLVM.py()) - (smearedLVM.px() * patElectron.py());
        double x1_denominator = (smearedLVM.py() * (patElectron.px() + MEpx)) - (smearedLVM.px() * (patElectron.py() + MEpy));
        double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
        double x2_numerator = x1_numerator;
        double x2_denominator = (patElectron.px() * (smearedLVM.py() + MEpy)) - (patElectron.py() * (smearedLVM.px() + MEpx));
        double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
        if ( (x1 > 0. && x1 < 1.) &&
             (x2 > 0. && x2 < 1.) ) {
          reco::Candidate::LorentzVector The_LorentzVect = (patElectron.p4() / x1) + (smearedLVM / x2);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
          return MassRecoInformation;
        } else {
          double px = patElectron.px() + smearedLVM.px() + MEpx;
          double py = patElectron.py() + smearedLVM.py() + MEpy;
          double pz = patElectron.pz() + smearedLVM.pz();
          double e = patElectron.energy() + smearedLVM.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
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
        double px = patElectron.px() + patMuon.px() + patMET.px() + deltaForMEx;
        double py = patElectron.py() + patMuon.py() + patMET.py() + deltaForMEy;
        double pz = patElectron.pz() + patMuon.pz();
        double e = patElectron.energy() + patMuon.energy() + TMath::Sqrt(((patMET.px() + deltaForMEx) * (patMET.px() + deltaForMEx)) + ((patMET.py() + deltaForMEy) * (patMET.py() + deltaForMEy)));
        reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else if(_UseCollinerApproxMassReco) {
        double x1_numerator = (patElectron.px() * patMuon.py()) - (patMuon.px() * patElectron.py());
        double x1_denominator = (patMuon.py() * (patElectron.px() + patMET.px() + deltaForMEx)) - (patMuon.px() * (patElectron.py() + patMET.py() + deltaForMEy));
        double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
        double x2_numerator = x1_numerator;
        double x2_denominator = (patElectron.px() * (patMuon.py() + patMET.py() + deltaForMEy)) - (patElectron.py() * (patMuon.px() + patMET.px() + deltaForMEx));
        double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
        if ( (x1 > 0. && x1 < 1.) &&
             (x2 > 0. && x2 < 1.) ) {
          reco::Candidate::LorentzVector The_LorentzVect = (patElectron.p4() / x1) + (patMuon.p4() / x2);
          pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
          return MassRecoInformation;
        } else {
          double px = patElectron.px() + patMuon.px() + patMET.px() + deltaForMEx;
          double py = patElectron.py() + patMuon.py() + patMET.py() + deltaForMEy;
          double pz = patElectron.pz() + patMuon.pz();
          double e = patElectron.energy() + patMuon.energy() + TMath::Sqrt(((patMET.px() + deltaForMEx) * (patMET.px() + deltaForMEx)) + ((patMET.py() + deltaForMEy) * (patMET.py() + deltaForMEy)));
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
    double MEpx = patMET.px() + deltaForMEx + patMuon1.px() - smearedLV1.px() + patMuon2.px() - smearedLV2.px();
    double MEpy = patMET.py() + deltaForMEy + patMuon1.py() - smearedLV1.py() + patMuon2.py() - smearedLV2.py();
    if(_UseVectorSumOfVisProductsAndMetMassReco) {
      double px = smearedLV1.px() + smearedLV2.px() + MEpx;
      double py = smearedLV1.py() + smearedLV2.py() + MEpy;
      double pz = smearedLV1.pz() + smearedLV2.pz();
      double e = smearedLV1.energy() + smearedLV2.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else if(_UseCollinerApproxMassReco) {
      double x1_numerator = (smearedLV1.px() * smearedLV2.py()) - (smearedLV2.px() * smearedLV1.py());
      double x1_denominator = (smearedLV2.py() * (smearedLV1.px() + MEpx)) - (smearedLV2.px() * (smearedLV1.py() + MEpy));
      double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
      double x2_numerator = x1_numerator;
      double x2_denominator = (smearedLV1.px() * (smearedLV2.py() + MEpy)) - (smearedLV1.py() * (smearedLV2.px() + MEpx));
      double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
      if ( (x1 > 0. && x1 < 1.) &&
           (x2 > 0. && x2 < 1.) ) {
        reco::Candidate::LorentzVector The_LorentzVect = (smearedLV1 / x1) + (smearedLV2 / x2);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else {
        double px = smearedLV1.px() + smearedLV2.px() + MEpx;
        double py = smearedLV1.py() + smearedLV2.py() + MEpy;
        double pz = smearedLV1.pz() + smearedLV2.pz();
        double e = smearedLV1.energy() + smearedLV2.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
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
      double px = patMuon1.px() + patMuon2.px() + patMET.px() + deltaForMEx;
      double py = patMuon1.py() + patMuon2.py() + patMET.py() + deltaForMEy;
      double pz = patMuon1.pz() + patMuon2.pz();
      double e = patMuon1.energy() + patMuon2.energy() + TMath::Sqrt(((patMET.px() + deltaForMEx) * (patMET.px() + deltaForMEx)) + ((patMET.py() + deltaForMEy) * (patMET.py() + deltaForMEy)));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else if(_UseCollinerApproxMassReco) {
      double x1_numerator = (patMuon1.px() * patMuon2.py()) - (patMuon2.px() * patMuon1.py());
      double x1_denominator = (patMuon2.py() * (patMuon1.px() + patMET.px() + deltaForMEx)) - (patMuon2.px() * (patMuon1.py() + patMET.py() + deltaForMEy));
      double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
      double x2_numerator = x1_numerator;
      double x2_denominator = (patMuon1.px() * (patMuon2.py() + patMET.py() + deltaForMEy)) - (patMuon1.py() * (patMuon2.px() + patMET.px() + deltaForMEx));
      double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
      if ( (x1 > 0. && x1 < 1.) &&
           (x2 > 0. && x2 < 1.) ) {
        reco::Candidate::LorentzVector The_LorentzVect = (patMuon1.p4() / x1) + (patMuon2.p4() / x2);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else {
        double px = patMuon1.px() + patMuon2.px() + patMET.px() + deltaForMEx;
        double py = patMuon1.py() + patMuon2.py() + patMET.py() + deltaForMEy;
        double pz = patMuon1.pz() + patMuon2.pz();
        double e = patMuon1.energy() + patMuon2.energy() + TMath::Sqrt(((patMET.px() + deltaForMEx) * (patMET.px() + deltaForMEx)) + ((patMET.py() + deltaForMEy) * (patMET.py() + deltaForMEy)));
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
    double MEpx = patMET.px() + deltaForMEx + patElectron1.px() - smearedLV1.px() + patElectron2.px() - smearedLV2.px();
    double MEpy = patMET.py() + deltaForMEy + patElectron1.py() - smearedLV1.py() + patElectron2.py() - smearedLV2.py();
    if(_UseVectorSumOfVisProductsAndMetMassReco) {
      double px = smearedLV1.px() + smearedLV2.px() + MEpx;
      double py = smearedLV1.py() + smearedLV2.py() + MEpy;
      double pz = smearedLV1.pz() + smearedLV2.pz();
      double e = smearedLV1.energy() + smearedLV2.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else if(_UseCollinerApproxMassReco) {
      double x1_numerator = (smearedLV1.px() * smearedLV2.py()) - (smearedLV2.px() * smearedLV1.py());
      double x1_denominator = (smearedLV2.py() * (smearedLV1.px() + MEpx)) - (smearedLV2.px() * (smearedLV1.py() + MEpy));
      double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
      double x2_numerator = x1_numerator;
      double x2_denominator = (smearedLV1.px() * (smearedLV2.py() + MEpy)) - (smearedLV1.py() * (smearedLV2.px() + MEpx));
      double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
      if ( (x1 > 0. && x1 < 1.) &&
           (x2 > 0. && x2 < 1.) ) {
        reco::Candidate::LorentzVector The_LorentzVect = (smearedLV1 / x1) + (smearedLV2 / x2);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else {
        double px = smearedLV1.px() + smearedLV2.px() + MEpx;
        double py = smearedLV1.py() + smearedLV2.py() + MEpy;
        double pz = smearedLV1.pz() + smearedLV2.pz();
        double e = smearedLV1.energy() + smearedLV2.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
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
      double px = patElectron1.px() + patElectron2.px() + patMET.px() + deltaForMEx;
      double py = patElectron1.py() + patElectron2.py() + patMET.py() + deltaForMEy;
      double pz = patElectron1.pz() + patElectron2.pz();
      double e = patElectron1.energy() + patElectron2.energy() + TMath::Sqrt(((patMET.px() + deltaForMEx) * (patMET.px() + deltaForMEx)) + ((patMET.py() + deltaForMEy) * (patMET.py() + deltaForMEy)));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else if(_UseCollinerApproxMassReco) {
      double x1_numerator = (patElectron1.px() * patElectron2.py()) - (patElectron2.px() * patElectron1.py());
      double x1_denominator = (patElectron2.py() * (patElectron1.px() + patMET.px() + deltaForMEx)) - (patElectron2.px() * (patElectron1.py() + patMET.py() + deltaForMEy));
      double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
      double x2_numerator = x1_numerator;
      double x2_denominator = (patElectron1.px() * (patElectron2.py() + patMET.py() + deltaForMEy)) - (patElectron1.py() * (patElectron2.px() + patMET.px() + deltaForMEx));
      double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
      if ( (x1 > 0. && x1 < 1.) &&
           (x2 > 0. && x2 < 1.) ) {
        reco::Candidate::LorentzVector The_LorentzVect = (patElectron1.p4() / x1) + (patElectron2.p4() / x2);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else {
        double px = patElectron1.px() + patElectron2.px() + patMET.px() + deltaForMEx;
        double py = patElectron1.py() + patElectron2.py() + patMET.py() + deltaForMEy;
        double pz = patElectron1.pz() + patElectron2.pz();
        double e = patElectron1.energy() + patElectron2.energy() + TMath::Sqrt(((patMET.px() + deltaForMEx) * (patMET.px() + deltaForMEx)) + ((patMET.py() + deltaForMEy) * (patMET.py() + deltaForMEy)));
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
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::CalculateThe4Momentum(const pat::Tau& patTau1, bool smear1, reco::Candidate::LorentzVector smearedLV1, const pat::Tau& patTau2, bool smear2, reco::Candidate::LorentzVector smearedLV2, const pat::MET& patMET) {
  if(smear1) {
    double MEpx = patMET.px() + deltaForMEx + patTau1.px() - smearedLV1.px() + patTau2.px() - smearedLV2.px();
    double MEpy = patMET.py() + deltaForMEy + patTau1.py() - smearedLV1.py() + patTau2.py() - smearedLV2.py();
    if(_UseVectorSumOfVisProductsAndMetMassReco) {
      double px = smearedLV1.px() + smearedLV2.px() + MEpx;
      double py = smearedLV1.py() + smearedLV2.py() + MEpy;
      double pz = smearedLV1.pz() + smearedLV2.pz();
      double e = smearedLV1.energy() + smearedLV2.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else if(_UseCollinerApproxMassReco) {
      double x1_numerator = (smearedLV1.px() * smearedLV2.py()) - (smearedLV2.px() * smearedLV1.py());
      double x1_denominator = (smearedLV2.py() * (smearedLV1.px() + MEpx)) - (smearedLV2.px() * (smearedLV1.py() + MEpy));
      double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
      double x2_numerator = x1_numerator;
      double x2_denominator = (smearedLV1.px() * (smearedLV2.py() + MEpy)) - (smearedLV1.py() * (smearedLV2.px() + MEpx));
      double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
      if ( (x1 > 0. && x1 < 1.) &&
           (x2 > 0. && x2 < 1.) ) {
        reco::Candidate::LorentzVector The_LorentzVect = (smearedLV1 / x1) + (smearedLV2 / x2);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else {
        double px = smearedLV1.px() + smearedLV2.px() + MEpx;
        double py = smearedLV1.py() + smearedLV2.py() + MEpy;
        double pz = smearedLV1.pz() + smearedLV2.pz();
        double e = smearedLV1.energy() + smearedLV2.energy() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
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
      double px = patTau1.px() + patTau2.px() + patMET.px() + deltaForMEx;
      double py = patTau1.py() + patTau2.py() + patMET.py() + deltaForMEy;
      double pz = patTau1.pz() + patTau2.pz();
      double e = patTau1.energy() + patTau2.energy() + TMath::Sqrt(((patMET.px() + deltaForMEx) * (patMET.px() + deltaForMEx)) + ((patMET.py() + deltaForMEy) * (patMET.py() + deltaForMEy)));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else if(_UseCollinerApproxMassReco) {
      double x1_numerator = (patTau1.px() * patTau2.py()) - (patTau2.px() * patTau1.py());
      double x1_denominator = (patTau2.py() * (patTau1.px() + patMET.px() + deltaForMEx)) - (patTau2.px() * (patTau1.py() + patMET.py() + deltaForMEy));
      double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
      double x2_numerator = x1_numerator;
      double x2_denominator = (patTau1.px() * (patTau2.py() + patMET.py() + deltaForMEy)) - (patTau1.py() * (patTau2.px() + patMET.px() + deltaForMEx));
      double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
      if ( (x1 > 0. && x1 < 1.) &&
           (x2 > 0. && x2 < 1.) ) {
        reco::Candidate::LorentzVector The_LorentzVect = (patTau1.p4() / x1) + (patTau2.p4() / x2);
        pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
        return MassRecoInformation;
      } else {
        double px = patTau1.px() + patTau2.px() + patMET.px() + deltaForMEx;
        double py = patTau1.py() + patTau2.py() + patMET.py() + deltaForMEy;
        double pz = patTau1.pz() + patTau2.pz();
        double e = patTau1.energy() + patTau2.energy() + TMath::Sqrt(((patMET.px() + deltaForMEx) * (patMET.px() + deltaForMEx)) + ((patMET.py() + deltaForMEy) * (patMET.py() + deltaForMEy)));
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
}

//-----Calculate lepton+met transverse mass
double HiMassTauAnalysis::CalculateLeptonMetMt(const pat::Muon& patMuon, bool smear, reco::Candidate::LorentzVector smearedLV, const pat::MET& patMET) {
  if(smear) {
    double MEpx = patMET.px() + deltaForMEx + patMuon.px() - smearedLV.px();
    double MEpy = patMET.py() + deltaForMEy + patMuon.py() - smearedLV.py();
    double px = smearedLV.px() + MEpx;
    double py = smearedLV.py() + MEpy;
    double et = smearedLV.Et() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
    double mt2 = et*et - (px*px + py*py);
    if ( mt2 < 0 ) { return -1.; }
    else { return sqrt(mt2); }
  } else {
    double px = patMuon.px() + patMET.px() + deltaForMEx;
    double py = patMuon.py() + patMET.py() + deltaForMEy;
    double et = patMuon.et() + TMath::Sqrt(((patMET.px() + deltaForMEx) * (patMET.px() + deltaForMEx)) + ((patMET.py() + deltaForMEy) * (patMET.py() + deltaForMEy)));
    double mt2 = et*et - (px*px + py*py);
    if ( mt2 < 0 ) { return -1.; }
    else { return sqrt(mt2); }
  }
}
double HiMassTauAnalysis::CalculateLeptonMetMt(const pat::Electron& patElectron, bool smear, reco::Candidate::LorentzVector smearedLV, const pat::MET& patMET) {
  if(smear) {
    double MEpx = patMET.px() + deltaForMEx + patElectron.px() - smearedLV.px();
    double MEpy = patMET.py() + deltaForMEy + patElectron.py() - smearedLV.py();
    double px = smearedLV.px() + MEpx;
    double py = smearedLV.py() + MEpy;
    double et = smearedLV.Et() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
    double mt2 = et*et - (px*px + py*py);
    if ( mt2 < 0 ) { return -1.; }
    else { return sqrt(mt2); }
  } else {
    double px = patElectron.px() + patMET.px() + deltaForMEx;
    double py = patElectron.py() + patMET.py() + deltaForMEy;
    double et = patElectron.et() + TMath::Sqrt(((patMET.px() + deltaForMEx) * (patMET.px() + deltaForMEx)) + ((patMET.py() + deltaForMEy) * (patMET.py() + deltaForMEy)));
    double mt2 = et*et - (px*px + py*py);
    if ( mt2 < 0 ) { return -1.; }
    else { return sqrt(mt2); }
  }
}
double HiMassTauAnalysis::CalculateLeptonMetMt(const pat::Tau& patTau, bool smear, reco::Candidate::LorentzVector smearedLV, const pat::MET& patMET) {
  if(smear) {
    double MEpx = patMET.px() + deltaForMEx + patTau.px() - smearedLV.px();
    double MEpy = patMET.py() + deltaForMEy + patTau.py() - smearedLV.py();
    double px = smearedLV.px() + MEpx;
    double py = smearedLV.py() + MEpy;
    double et = smearedLV.Et() + TMath::Sqrt((MEpx * MEpx) + (MEpy * MEpy));
    double mt2 = et*et - (px*px + py*py);
    if ( mt2 < 0 ) { return -1.; }
    else { return sqrt(mt2); }
  } else {
    double px = patTau.px() + patMET.px() + deltaForMEx;
    double py = patTau.py() + patMET.py() + deltaForMEy;
    double et = patTau.et() + TMath::Sqrt(((patMET.px() + deltaForMEx) * (patMET.px() + deltaForMEx)) + ((patMET.py() + deltaForMEy) * (patMET.py() + deltaForMEy)));
    double mt2 = et*et - (px*px + py*py);
    if ( mt2 < 0 ) { return -1.; }
    else { return sqrt(mt2); }
  }
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
  if (_UseRecoTauEllipseForEcalIso) {
    if (patTau.isCaloTau()) {
    } else {
      if (patTau.leadPFChargedHadrCand().isNonnull()) {
        PFCandidateRefVector TauSigGammas = patTau.signalPFGammaCands();
        for(PFCandidateRefVector::const_iterator iGam=TauSigGammas.begin();iGam!=TauSigGammas.end();++iGam) {
          if((**iGam).pt()>_RecoTauSigGamThreshold) {
            double deltaphi = TMath::Abs(normalizedPhi((**iGam).phi() - patTau.leadPFChargedHadrCand()->phi()));
            double deltaeta = (**iGam).eta() - patTau.leadPFChargedHadrCand()->eta();
            if( ( ((deltaphi * deltaphi) / (_RecoTauEcalIsoRphiForEllipse * _RecoTauEcalIsoRphiForEllipse)) + 
                  ((deltaeta * deltaeta) / (_RecoTauEcalIsoRetaForEllipse * _RecoTauEcalIsoRetaForEllipse)) ) <= 1.0 ) {
              nSigGams++;
            }
          }
        }
      }
    }
  } else {
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
  }
  return nSigGams;
}

//-----Calculate Mass from Tau signal track constituents
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
/*
      px += (**iTrk).momentum().x(),
      py += (**iTrk).momentum().y(),
      pz += (**iTrk).momentum().z(),
      e += sqrt(pow((double)(**iTrk).momentum().r(),2)+pow(0.13957018,2));
*/
      px += (**iTrk).px(),
      py += (**iTrk).py(),
      pz += (**iTrk).pz(),
      e += (**iTrk).energy();
    }
  }
  reco::Candidate::LorentzVector TheSignalTracks_LorentzVect(px, py, pz, e);
  return TheSignalTracks_LorentzVect;
}

//-----Calculate Mass from Tau signal track AND gamma constituents
reco::Candidate::LorentzVector HiMassTauAnalysis::CalculateTauSignalTracksAndGammasMass(const pat::Tau& patTau) {
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
/*
      px += (**iTrk).momentum().x(),
      py += (**iTrk).momentum().y(),
      pz += (**iTrk).momentum().z(),
      e += sqrt(pow((double)(**iTrk).momentum().r(),2)+pow(0.13957018,2));
*/
      px += (**iTrk).px(),
      py += (**iTrk).py(),
      pz += (**iTrk).pz(),
      e += (**iTrk).energy();
    }
  }
  if (_UseRecoTauEllipseForEcalIso) {
    if (patTau.isCaloTau()) {
    } else {
      if (patTau.leadPFChargedHadrCand().isNonnull()) {
        PFCandidateRefVector TauSigGammas = patTau.signalPFGammaCands();
        for(PFCandidateRefVector::const_iterator iGam=TauSigGammas.begin();iGam!=TauSigGammas.end();++iGam) {
          double deltaphi = TMath::Abs(normalizedPhi((**iGam).phi() - patTau.leadPFChargedHadrCand()->phi()));
          double deltaeta = (**iGam).eta() - patTau.leadPFChargedHadrCand()->eta();
          if( ( ((deltaphi * deltaphi) / (_RecoTauEcalIsoRphiForEllipse * _RecoTauEcalIsoRphiForEllipse)) +
                ((deltaeta * deltaeta) / (_RecoTauEcalIsoRetaForEllipse * _RecoTauEcalIsoRetaForEllipse)) ) <= 1.0 ) {
            if((**iGam).pt() > _RecoTauSigGamThreshold) {
              px += (**iGam).px(),
              py += (**iGam).py(),
              pz += (**iGam).pz(),
              e += (**iGam).energy();
            }
          }
        }
      }
    }
  } else {
    if (patTau.isCaloTau()) {
    } else {
      if (patTau.leadPFChargedHadrCand().isNonnull()) {
        PFCandidateRefVector TauSigGammas = patTau.signalPFGammaCands();
        for(PFCandidateRefVector::const_iterator iGam=TauSigGammas.begin();iGam!=TauSigGammas.end();++iGam) {
          if((**iGam).pt() > _RecoTauSigGamThreshold) {
            px += (**iGam).px(),
            py += (**iGam).py(),
            pz += (**iGam).pz(),
            e += (**iGam).energy();
          }
        }
      }
    }
  }
  reco::Candidate::LorentzVector TheSignalTracksAndGammas_LorentzVect(px, py, pz, e);
  return TheSignalTracksAndGammas_LorentzVect;
}

//-----Smear the light leptons (for studies of systematic uncertanties)
reco::Candidate::LorentzVector HiMassTauAnalysis::SmearLightLepton(const pat::Muon& patMuon) {
  double smearedPt;
  double smearedEta;
  double smearedPhi;
  if(_GenParticleSource.label() != "") {
    if(matchToGen(patMuon).first) {
      reco::Candidate::LorentzVector unsmearedMomentum = matchToGen(patMuon).second;
      if(_SmearThePt) {
        smearedPt = (unsmearedMomentum.pt() * _MuonPtScaleOffset) + ((patMuon.pt() -  unsmearedMomentum.pt()) * _MuonPtSigmaOffset);
      } else {smearedPt  = patMuon.pt();}
      if(_SmearTheEta) {
        smearedEta = (unsmearedMomentum.eta() * _MuonEtaScaleOffset) + ((patMuon.eta() -  unsmearedMomentum.eta()) * _MuonEtaSigmaOffset);
      } else {smearedEta  = patMuon.eta();}
      if(_SmearThePhi) {
        smearedPhi = (unsmearedMomentum.phi() * _MuonPhiScaleOffset) + ((patMuon.phi() -  unsmearedMomentum.phi()) * _MuonPhiSigmaOffset);
      } else {smearedPhi  = patMuon.phi();}
      math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector(smearedPt, smearedEta, smearedPhi, unsmearedMomentum.mass());
      reco::Candidate::LorentzVector smearedMomentum(smearedPtEtaPhiMVector.px(), smearedPtEtaPhiMVector.py(), smearedPtEtaPhiMVector.pz(), smearedPtEtaPhiMVector.energy());
      return smearedMomentum;
    } else {
      reco::Candidate::LorentzVector smearedMomentum = patMuon.p4();
      return smearedMomentum;
    }
  } else {
    reco::Candidate::LorentzVector smearedMomentum = patMuon.p4();
    return smearedMomentum;
  }
}
reco::Candidate::LorentzVector HiMassTauAnalysis::SmearLightLepton(const pat::Electron& patElectron) {
  double smearedPt;
  double smearedEta;
  double smearedPhi;
  if(_GenParticleSource.label() != "") {
    if(matchToGen(patElectron).first) {
      reco::Candidate::LorentzVector unsmearedMomentum = matchToGen(patElectron).second; 
      if(_SmearThePt) {
        smearedPt = (unsmearedMomentum.pt() * _ElectronPtScaleOffset) + ((patElectron.pt() -  unsmearedMomentum.pt()) * _ElectronPtSigmaOffset);
      } else {smearedPt  = patElectron.pt();}
      if(_SmearTheEta) {
        smearedEta = (unsmearedMomentum.eta() * _ElectronEtaScaleOffset) + ((patElectron.eta() -  unsmearedMomentum.eta()) * _ElectronEtaSigmaOffset);
      } else {smearedEta  = patElectron.eta();}
      if(_SmearThePhi) {
        smearedPhi = (unsmearedMomentum.phi() * _ElectronPhiScaleOffset) + ((patElectron.phi() -  unsmearedMomentum.phi()) * _ElectronPhiSigmaOffset);
      } else {smearedPhi  = patElectron.phi();}
      math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector(smearedPt, smearedEta, smearedPhi, unsmearedMomentum.mass());
      reco::Candidate::LorentzVector smearedMomentum(smearedPtEtaPhiMVector.px(), smearedPtEtaPhiMVector.py(), smearedPtEtaPhiMVector.pz(), smearedPtEtaPhiMVector.energy());
      return smearedMomentum;
    } else {
      reco::Candidate::LorentzVector smearedMomentum = patElectron.p4();
      return smearedMomentum;
    }
  } else {
    reco::Candidate::LorentzVector smearedMomentum = patElectron.p4();
    return smearedMomentum;
  }
}
reco::Candidate::LorentzVector HiMassTauAnalysis::SmearTau(const pat::Tau& patTau) {
  double smearedPt;
  double smearedEta;
  double smearedPhi;
  if(_GenParticleSource.label() != "") {
    if(matchToGen(patTau).first) {
      reco::Candidate::LorentzVector unsmearedMomentum = matchToGen(patTau).second; 
      if(_SmearThePt) {
        smearedPt = (unsmearedMomentum.pt() * _TauPtScaleOffset) + ((patTau.pt() -  unsmearedMomentum.pt()) * _TauPtSigmaOffset);
      } else {smearedPt  = patTau.pt();}
      if(_SmearTheEta) {
        smearedEta = (unsmearedMomentum.eta() * _TauEtaScaleOffset) + ((patTau.eta() -  unsmearedMomentum.eta()) * _TauEtaSigmaOffset);
      } else {smearedEta  = patTau.eta();}
      if(_SmearThePhi) {
        smearedPhi = (unsmearedMomentum.phi() * _TauPhiScaleOffset) + ((patTau.phi() -  unsmearedMomentum.phi()) * _TauPhiSigmaOffset);
      } else {smearedPhi  = patTau.phi();}
      math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector(smearedPt, smearedEta, smearedPhi, unsmearedMomentum.mass());
      reco::Candidate::LorentzVector smearedMomentum(smearedPtEtaPhiMVector.px(), smearedPtEtaPhiMVector.py(), smearedPtEtaPhiMVector.pz(), smearedPtEtaPhiMVector.energy());
      return smearedMomentum;
    } else {
      reco::Candidate::LorentzVector smearedMomentum = patTau.p4();
      return smearedMomentum;
    }
  } else {
    reco::Candidate::LorentzVector smearedMomentum = patTau.p4();
    return smearedMomentum;
  }
}
reco::Candidate::LorentzVector HiMassTauAnalysis::SmearJet(const pat::Jet& patJet) {
  bool isRealJet = true;
  reco::Candidate::LorentzVector tempJetVector;
  if(_UseCorrectedJet) {
    tempJetVector = patJet.p4();
  } else {
    tempJetVector = patJet.correctedJet("raw","").p4();
  }
  if(_GenParticleSource.label() != "") {
    for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
      if( (reco::deltaR(patMuon->p4(), tempJetVector) < _JetMuonMatchingDeltaR) && (matchToGen(*patMuon).first) ) {isRealJet = false;}
    }
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron ) {
      if( (reco::deltaR(patElectron->p4(), tempJetVector) < _JetElectronMatchingDeltaR) && (matchToGen(*patElectron).first) ) {isRealJet = false;}
    }
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
      if( (reco::deltaR(patTau->p4(), tempJetVector) < _JetTauMatchingDeltaR) && (matchToGen(*patTau).first) ) {isRealJet = false;}
    }
    if(isRealJet) {
      reco::Candidate::LorentzVector smearedMomentum = _JetEnergyScaleOffset * tempJetVector;
      return smearedMomentum;
    } else {
      reco::Candidate::LorentzVector smearedMomentum = tempJetVector;
      return smearedMomentum;
    }
  } else {
    reco::Candidate::LorentzVector smearedMomentum = tempJetVector;
    return smearedMomentum;
  }
}

//-----Initialize information for the calculation of pdf systematic uncertaintites
void HiMassTauAnalysis::InitializeInfoForPDFSystematicUncertaintites() {
  for (unsigned int i=0; i<pdfWeightTags_.size(); ++i) {
    std::cout << "\t" << pdfWeightTags_[i].label();
    pdfStart_Denominator_.push_back(-1);
    pdfStart_Numerator_.push_back(-1);
  }
}

double HiMassTauAnalysis::alphaRatio(double pt) {

      double pigaga = 0.;

      // Leptonic contribution (just one loop, precise at < 0.3% level)
      const double alphapi = 1/137.036/TMath::Pi();
      const double mass_e = 0.0005;
      const double mass_mu = 0.106;
      const double mass_tau = 1.777;
      const double mass_Z = 91.2;
      if (pt>mass_e) pigaga += alphapi * (2*log(pt/mass_e)/3.-5./9.);
      if (pt>mass_mu) pigaga += alphapi * (2*log(pt/mass_mu)/3.-5./9.);
      if (pt>mass_tau) pigaga += alphapi * (2*log(pt/mass_tau)/3.-5./9.);

      // Hadronic vaccum contribution
      // Using simple effective parametrization from Physics Letters B 513 (2001) 46.
      // Top contribution neglected
      double A = 0.; 
      double B = 0.; 
      double C = 0.; 
      if (pt<0.7) {
            A = 0.0; B = 0.0023092; C = 3.9925370;
      } else if (pt<2.0) {
            A = 0.0; B = 0.0022333; C = 4.2191779;
      } else if (pt<4.0) {
            A = 0.0; B = 0.0024402; C = 3.2496684;
      } else if (pt<10.0) {
            A = 0.0; B = 0.0027340; C = 2.0995092;
      } else if (pt<mass_Z) {
            A = 0.0010485; B = 0.0029431; C = 1.0;
      } else if (pt<10000.) {
            A = 0.0012234; B = 0.0029237; C = 1.0;
      } else {
            A = 0.0016894; B = 0.0028984; C = 1.0;
      }
      pigaga += A + B*log(1.+C*pt*pt);

      // Done
      return 1./(1.-pigaga);
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

  // Initialize TFileService
  Service<TFileService> fs;

  // Initialize stringstream used to name histograms for each PDF weight
  std::stringstream j;
  j.str("");

  // The loop below is used to create a different histogram for each event weighting factor (for systematic studies or renomarlization w/ respect to data).
  // If reweighting booleans are set to false, then event weight will be set to 1 and only 1 histogram per variable will be created.
  for(unsigned int NpdfCounter = 0; NpdfCounter < pdfWeightVector.size();  NpdfCounter++){
    j << NpdfCounter;

    //--- histogram containing the number of events analyzed and number passing specificied cuts
    _hEvents[NpdfCounter] = fs->make<TH1F>(("Events_"+j.str()).c_str(), ("Events_"+j.str()).c_str(), 2, 0., 2.);

    //--- book vertex histograms
    if (_FillRecoVertexHists) {
      _hVertexZposition[NpdfCounter] = fs->make<TH1F>(("VertexZposition_"+j.str()).c_str(), ("VertexZposition_"+j.str()).c_str(), 50, -50., 50.);
      _hVertexNTracks[NpdfCounter]   = fs->make<TH1F>(("VertexNTracks_"+j.str()).c_str(),   ("VertexNTracks_"+j.str()).c_str(),   100, 0., 100.);
      _hNVertices[NpdfCounter]       = fs->make<TH1F>(("NVertices_"+j.str()).c_str(),       ("NVertices_"+j.str()).c_str(),       10, 0., 10.);
    }
    
    //--- book generator level histograms
    if (_FillGenTauHists) {
      _hNGenTau[NpdfCounter]                 = fs->make<TH1F>(("NGenTau_"+j.str()).c_str(),           ("NGenTau_"+j.str()).c_str(),      21, 0., 20.);
      _hGenTauEnergy[NpdfCounter]            = fs->make<TH1F>(("GenTauEnergy_"+j.str()).c_str(),      ("GenTauEnergy_"+j.str()).c_str(), 200, 0., 500.);
      _hGenTauPt[NpdfCounter]                = fs->make<TH1F>(("GenTauPt_"+j.str()).c_str(),          ("GenTauPt_"+j.str()).c_str(),     200, 0., 500.);
      _hGenTauEta[NpdfCounter]               = fs->make<TH1F>(("GenTauEta_"+j.str()).c_str(),         ("GenTauEta_"+j.str()).c_str(), 72, -3.6, +3.6);
      _hGenTauPhi[NpdfCounter]               = fs->make<TH1F>(("GenTauPhi_"+j.str()).c_str(),         ("GenTauPhi_"+j.str()).c_str(), 36, -TMath::Pi(), +TMath::Pi());
      _hGenTauMotherEnergy[NpdfCounter]      = fs->make<TH1F>(("GenTauMotherEnergy_"+j.str()).c_str(),("GenTauMotherEnergy_"+j.str()).c_str(),200, 0., 500.);
      _hGenTauMotherPt[NpdfCounter]          = fs->make<TH1F>(("GenTauMotherPt_"+j.str()).c_str(),    ("GenTauMotherPt_"+j.str()).c_str(),    200, 0., 500.);
      _hGenTauMotherEta[NpdfCounter]         = fs->make<TH1F>(("GenTauMotherEta_"+j.str()).c_str(),   ("GenTauMotherEta_"+j.str()).c_str(),   72, -3.6, +3.6);
      _hGenTauMotherPhi[NpdfCounter]         = fs->make<TH1F>(("GenTauMotherPhi_"+j.str()).c_str(),   ("GenTauMotherPhi_"+j.str()).c_str(),   36, -TMath::Pi(), +TMath::Pi());
      _hGenTauGrandMotherEnergy[NpdfCounter] = fs->make<TH1F>(("GenTauGrandMotherEnergy_"+j.str()).c_str(), ("GenTauGrandMotherEnergy_"+j.str()).c_str(), 200, 0., 500.);
      _hGenTauGrandMotherPt[NpdfCounter]     = fs->make<TH1F>(("GenTauGrandMotherPt_"+j.str()).c_str(),("GenTauGrandMotherPt_"+j.str()).c_str(), 200, 0., 500.);
      _hGenTauGrandMotherEta[NpdfCounter]    = fs->make<TH1F>(("GenTauGrandMotherEta_"+j.str()).c_str(),("GenTauGrandMotherEta_"+j.str()).c_str(),72, -3.6, +3.6);
      _hGenTauGrandMotherPhi[NpdfCounter]    = fs->make<TH1F>(("GenTauGrandMotherPhi_"+j.str()).c_str(),("GenTauGrandMotherPhi_"+j.str()).c_str(),36, -TMath::Pi(), +TMath::Pi());
    }
    
    //--- book reconstruction level histograms 
    if (_FillRecoTauHists) {
      // Tau basic physics quantities
      _hNTau[NpdfCounter]                             = fs->make<TH1F>(("NTau_"+j.str()).c_str(),         ("NTau_"+j.str()).c_str(),         21, 0., 20.);
      _hTauJetEnergy[NpdfCounter]                     = fs->make<TH1F>(("TauJetEnergy_"+j.str()).c_str(), ("TauJetEnergy_"+j.str()).c_str(), 200, 0., 500.);
      _hTauJetPt[NpdfCounter]                         = fs->make<TH1F>(("TauJetPt_"+j.str()).c_str(),     ("TauJetPt_"+j.str()).c_str(),     200, 0., 500.);
      _hTauJetEta[NpdfCounter]                        = fs->make<TH1F>(("TauJetEta_"+j.str()).c_str(),    ("TauJetEta_"+j.str()).c_str(),    72, -3.6, +3.6);
      _hTauJetPhi[NpdfCounter]                        = fs->make<TH1F>(("TauJetPhi_"+j.str()).c_str(),    ("TauJetPhi_"+j.str()).c_str(),    36, -TMath::Pi(), +TMath::Pi());
      _hTauJetNumSignalTracks[NpdfCounter]            = fs->make<TH1F>(("TauJetNumSignalTracks_"+j.str()).c_str(),           ("TauJetNumSignalTracks_"+j.str()).c_str(), 10, 0, 10);
      _hTauJetNumSignalGammas[NpdfCounter]            = fs->make<TH1F>(("TauJetNumSignalGammas_"+j.str()).c_str(),           ("TauJetNumSignalGammas_"+j.str()).c_str(), 10, 0, 10);
      _hTauJetSeedTrackPt[NpdfCounter]                = fs->make<TH1F>(("TauJetSeedTrackPt_"+j.str()).c_str(),               ("TauJetSeedTrackPt_"+j.str()).c_str(),     200, 0., 500.);
      _hTauJetSeedTrackIpSignificance[NpdfCounter]    = fs->make<TH1F>(("TauJetSeedTrackIpSignificance_"+j.str()).c_str(),   ("TauJetSeedTrackIpSignificance_"+j.str()).c_str(), 100, 0., 100.);
      _hTauJetSeedTrackNhits[NpdfCounter]             = fs->make<TH1F>(("TauJetSeedTrackNhits_"+j.str()).c_str(),            ("TauJetSeedTrackNhits_"+j.str()).c_str(), 40, 0., 40.);
      _hTauJetSeedTrackChi2[NpdfCounter]              = fs->make<TH1F>(("TauJetSeedTrackChi2_"+j.str()).c_str(),             ("TauJetSeedTrackChi2_"+j.str()).c_str(), 50, 0., 100.);
      _hTauJetCharge[NpdfCounter]                     = fs->make<TH1F>(("TauJetCharge_"+j.str()).c_str(),                    ("TauJetCharge_"+j.str()).c_str(), 5, 0., 5.);
      _hTauJetSignalTracksMass[NpdfCounter]           = fs->make<TH1F>(("TauJetSignalTracksMass_"+j.str()).c_str(),          ("TauJetSignalTracksMass_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJetSignalTracksAndGammasMass[NpdfCounter]  = fs->make<TH1F>(("TauJetSignalTracksAndGammasMass_"+j.str()).c_str(), ("TauJetSignalTracksAndGammasMass_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJetSignalTracksChargeFraction[NpdfCounter] = fs->make<TH1F>(("TauJetSignalTracksChargeFraction_"+j.str()).c_str(),("TauJetSignalTracksChargeFraction_"+j.str()).c_str(), 30, 0., 1.5);
      _hTauJetNumIsoTracks[NpdfCounter]               = fs->make<TH1F>(("TauJetNumIsoTracks_"+j.str()).c_str(),              ("TauJetNumIsoTracks_"+j.str()).c_str(), 10, 0, 10);
      _hTauJetNumIsoGammas[NpdfCounter]               = fs->make<TH1F>(("TauJetNumIsoGammas_"+j.str()).c_str(),              ("TauJetNumIsoGammas_"+j.str()).c_str(), 10, 0, 10);
      _hTauJetNumIsoCands[NpdfCounter]                = fs->make<TH1F>(("TauJetNumIsoCands_"+j.str()).c_str(),               ("TauJetNumIsoCands_"+j.str()).c_str(), 10, 0, 10);
      _hTauJetSumPtIsoTracks[NpdfCounter]             = fs->make<TH1F>(("TauJetSumPtIsoTracks_"+j.str()).c_str(),            ("TauJetSumPtIsoTracks_"+j.str()).c_str(), 100, 0, 50);
      _hTauJetSumPtIsoGammas[NpdfCounter]             = fs->make<TH1F>(("TauJetSumPtIsoGammas_"+j.str()).c_str(),            ("TauJetSumPtIsoGammas_"+j.str()).c_str(), 100, 0, 50);
      _hTauJetSumPtIso[NpdfCounter]                   = fs->make<TH1F>(("TauJetSumPtIso_"+j.str()).c_str(),                  ("TauJetSumPtIso_"+j.str()).c_str(), 100, 0, 50);
      _hTauJetGenTauDeltaPhi[NpdfCounter]             = fs->make<TH1F>(("TauJetGenTauDeltaPhi_"+j.str()).c_str(),            ("TauJetGenTauDeltaPhi_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hTauJetGenTauDeltaEta[NpdfCounter]             = fs->make<TH1F>(("TauJetGenTauDeltaEta_"+j.str()).c_str(),            ("TauJetGenTauDeltaEta_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hTauJetGenTauDeltaPt[NpdfCounter]              = fs->make<TH1F>(("TauJetGenTauDeltaPt_"+j.str()).c_str(),             ("TauJetGenTauDeltaPt_"+j.str()).c_str(), 500, -5, 5);
      _hTauJetSignalTracksMass1prong[NpdfCounter]     = fs->make<TH1F>(("TauJetSignalTracksMass1prong_"+j.str()).c_str(),    ("TauJetSignalTracksMass1prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJetSignalTracksAndGammasMass1prong[NpdfCounter] = fs->make<TH1F>(("TauJetSignalTracksAndGammasMass1prong_"+j.str()).c_str(), ("TauJetSignalTracksAndGammasMass1prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJetSignalTracksMass3prong[NpdfCounter]          = fs->make<TH1F>(("TauJetSignalTracksMass3prong_"+j.str()).c_str(),          ("TauJetSignalTracksMass3prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJetSignalTracksAndGammasMass3prong[NpdfCounter] = fs->make<TH1F>(("TauJetSignalTracksAndGammasMass3prong_"+j.str()).c_str(), ("TauJetSignalTracksAndGammasMass3prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJetMass1Prong0Gamma[NpdfCounter]                = fs->make<TH1F>(("TauJetMass1Prong0Gamma_"+j.str()).c_str(),                ("TauJetMass1Prong0Gamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJetMass1Prong1Gamma[NpdfCounter]                = fs->make<TH1F>(("TauJetMass1Prong1Gamma_"+j.str()).c_str(),                ("TauJetMass1Prong1Gamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJetMass1Prong2orMoreGamma[NpdfCounter]          = fs->make<TH1F>(("TauJetMass1Prong2orMoreGamma_"+j.str()).c_str(),          ("TauJetMass1Prong2orMoreGamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJetMass3Prong0Gamma[NpdfCounter]                = fs->make<TH1F>(("TauJetMass3Prong0Gamma_"+j.str()).c_str(),                ("TauJetMass3Prong0Gamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJetMass3Prong1Gamma[NpdfCounter]                = fs->make<TH1F>(("TauJetMass3Prong1Gamma_"+j.str()).c_str(),                ("TauJetMass3Prong1Gamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJetMass3Prong2orMoreGamma[NpdfCounter]          = fs->make<TH1F>(("TauJetMass3Prong2orMoreGamma_"+j.str()).c_str(),          ("TauJetMass3Prong2orMoreGamma_"+j.str()).c_str(), 100, 0., 5.);
    }
    
    //--- book reconstruction level histograms 
    if (_FillRecoMuonHists) {
      _hNMuon[NpdfCounter]                       = fs->make<TH1F>(("NMuon_"+j.str()).c_str(),                   ("NMuon_"+j.str()).c_str(), 21, 0., 20.);
      _hMuonEnergy[NpdfCounter]                  = fs->make<TH1F>(("MuonEnergy_"+j.str()).c_str(),              ("MuonEnergy_"+j.str()).c_str(), 200, 0., 500.);
      _hMuonPt[NpdfCounter]                      = fs->make<TH1F>(("MuonPt_"+j.str()).c_str(),                  ("MuonPt_"+j.str()).c_str(),  200, 0., 500.);
      _hMuonEta[NpdfCounter]                     = fs->make<TH1F>(("MuonEta_"+j.str()).c_str(),                 ("MuonEta_"+j.str()).c_str(), 72, -3.6, +3.6);
      _hMuonPhi[NpdfCounter]                     = fs->make<TH1F>(("MuonPhi_"+j.str()).c_str(),                 ("MuonPhi_"+j.str()).c_str(), 36, -TMath::Pi(), +TMath::Pi());
      _hMuonTrackIso[NpdfCounter]                = fs->make<TH1F>(("MuonTrackIso_"+j.str()).c_str(),            ("MuonTrackIso_"+j.str()).c_str(), 100, 0, 50);
      _hMuonEcalIso[NpdfCounter]                 = fs->make<TH1F>(("MuonEcalIso_"+j.str()).c_str(),             ("MuonEcalIso_"+j.str()).c_str(), 100, 0, 50);
      _hMuonIso[NpdfCounter]                     = fs->make<TH1F>(("MuonIso_"+j.str()).c_str(),                 ("MuonIso_"+j.str()).c_str(), 100, 0, 50);
      _hMuonIp[NpdfCounter]                      = fs->make<TH1F>(("MuonIp_"+j.str()).c_str(),                  ("MuonIp_"+j.str()).c_str(), 500, -1, +1);
      _hMuonIpSignificance[NpdfCounter]          = fs->make<TH1F>(("MuonIpSignificance_"+j.str()).c_str(),      ("MuonIpSignificance_"+j.str()).c_str(), 100, 0., 100.);
      _hMuonGenMuonDeltaPhi[NpdfCounter]         = fs->make<TH1F>(("MuonGenMuonDeltaPhi_"+j.str()).c_str(),     ("MuonGenMuonDeltaPhi_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hMuonGenMuonDeltaEta[NpdfCounter]         = fs->make<TH1F>(("MuonGenMuonDeltaEta_"+j.str()).c_str(),     ("MuonGenMuonDeltaEta_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hMuonGenMuonDeltaPt[NpdfCounter]          = fs->make<TH1F>(("MuonGenMuonDeltaPt_"+j.str()).c_str(),      ("MuonGenMuonDeltaPt_"+j.str()).c_str(), 500, -5, 5);
      _hMuonCaloCompatibility[NpdfCounter]       = fs->make<TH1F>(("MuonCaloCompatibility_"+j.str()).c_str(),   ("MuonCaloCompatibility_"+j.str()).c_str(), 102, 0.0, 1.02);
      _hMuonSegmentCompatibility[NpdfCounter]    = fs->make<TH1F>(("MuonSegmentCompatibility_"+j.str()).c_str(),("MuonSegmentCompatibility_"+j.str()).c_str(), 102, 0.0, 1.02);
      _hMuonAntiPion[NpdfCounter]                = fs->make<TH1F>(("MuonAntiPion_"+j.str()).c_str(),            ("MuonAntiPion_"+j.str()).c_str(), 202, 0.0, 2.02);
      _hMuonCaloCompatibilityVsSegmentCompatibility[NpdfCounter] = fs->make<TH2F>(("MuonCaloCompatibilityVsSegmentCompatibility_"+j.str()).c_str(), ("MuonCaloCompatibilityVsSegmentCompatibility_"+j.str()).c_str(), 102, 0, 1.02, 102, 0, 1.02);
      
    }
    
    if (_FillRecoElectronHists) {
      _hNElectron[NpdfCounter]                   = fs->make<TH1F>(("NElectron_"+j.str()).c_str(),                  ("NElectron_"+j.str()).c_str(), 21, 0., 20.);
      _hElectronEnergy[NpdfCounter]              = fs->make<TH1F>(("ElectronEnergy_"+j.str()).c_str(),             ("ElectronEnergy_"+j.str()).c_str(), 200, 0., 500.);
      _hElectronPt[NpdfCounter]                  = fs->make<TH1F>(("ElectronPt_"+j.str()).c_str(),                 ("ElectronPt_"+j.str()).c_str(), 200, 0., 500.);
      _hElectronEta[NpdfCounter]                 = fs->make<TH1F>(("ElectronEta_"+j.str()).c_str(),                ("ElectronEta_"+j.str()).c_str(), 72, -3.6, +3.6);
      _hElectronPhi[NpdfCounter]                 = fs->make<TH1F>(("ElectronPhi_"+j.str()).c_str(),                ("ElectronPhi_"+j.str()).c_str(), 36, -TMath::Pi(), +TMath::Pi());
      _hElectronTrackIso[NpdfCounter]            = fs->make<TH1F>(("ElectronTrackIso_"+j.str()).c_str(),           ("ElectronTrackIso_"+j.str()).c_str(), 100, 0, 50);
      _hElectronEcalIso[NpdfCounter]             = fs->make<TH1F>(("ElectronEcalIso_"+j.str()).c_str(),            ("ElectronEcalIso_"+j.str()).c_str(), 100, 0, 50);
      _hElectronIp[NpdfCounter]                  = fs->make<TH1F>(("ElectronIp_"+j.str()).c_str(),                 ("ElectronIp_"+j.str()).c_str(), 500, -1, +1);
      _hElectronEoverP[NpdfCounter]              = fs->make<TH1F>(("ElectronEoverP_"+j.str()).c_str(),             ("ElectronEoverP_"+j.str()).c_str(), 60, 0, +3);
      _hElectronHoverEm[NpdfCounter]             = fs->make<TH1F>(("ElectronHoverEm_"+j.str()).c_str(),            ("ElectronHoverEm_"+j.str()).c_str(), 300, 0, +3);
      _hElectronClassification[NpdfCounter]      = fs->make<TH1F>(("ElectronClassification_"+j.str()).c_str(),     ("ElectronClassification_"+j.str()).c_str(), 200, 0, 200);
      _hElectronGenElectronDeltaPhi[NpdfCounter] = fs->make<TH1F>(("ElectronGenElectronDeltaPhi_"+j.str()).c_str(),("ElectronGenElectronDeltaPhi_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hElectronGenElectronDeltaEta[NpdfCounter] = fs->make<TH1F>(("ElectronGenElectronDeltaEta_"+j.str()).c_str(),("ElectronGenElectronDeltaEta_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hElectronGenElectronDeltaPt[NpdfCounter]  = fs->make<TH1F>(("ElectronGenElectronDeltaPt_"+j.str()).c_str(), ("ElectronGenElectronDeltaPt_"+j.str()).c_str(), 500, -5, 5);
      _hElectronEcalDriven[NpdfCounter]          = fs->make<TH1F>(("ElectronEcalDriven_"+j.str()).c_str(),         ("ElectronEcalDriven_"+j.str()).c_str(), 2, 0., 2.);
      _hElectronTrackerDriven[NpdfCounter]       = fs->make<TH1F>(("ElectronTrackerDriven_"+j.str()).c_str(),      ("ElectronTrackerDriven_"+j.str()).c_str(), 2, 0., 2.);
      _hElectronIsZee[NpdfCounter]               = fs->make<TH1F>(("ElectronIsZee_"+j.str()).c_str(),              ("ElectronIsZee_"+j.str()).c_str(), 2, 0., 2.);
    }
    
    if (_FillRecoJetHists) {
      _hNJet[NpdfCounter]        = fs->make<TH1F>(("NJet_"+j.str()).c_str(),      ("NJet_"+j.str()).c_str(), 21, 0., 20.);
      _hJetEnergy[NpdfCounter]   = fs->make<TH1F>(("JetEnergy_"+j.str()).c_str(), ("JetEnergy_"+j.str()).c_str(), 200, 0., 500.);
      _hJetPt[NpdfCounter]       = fs->make<TH1F>(("JetPt_"+j.str()).c_str(),     ("JetPt_"+j.str()).c_str(), 200, 0., 500.);
      _hJetEta[NpdfCounter]      = fs->make<TH1F>(("JetEta_"+j.str()).c_str(),    ("JetEta_"+j.str()).c_str(), 72, -3.6, +3.6);
      _hJetPhi[NpdfCounter]      = fs->make<TH1F>(("JetPhi_"+j.str()).c_str(),    ("JetPhi_"+j.str()).c_str(), 36, -TMath::Pi(), +TMath::Pi());
      _hBJetDiscrByTrackCounting[NpdfCounter]    = fs->make<TH1F>(("BJetDiscrByTrackCounting_"+j.str()).c_str(), ("BJetDiscrByTrackCounting_"+j.str()).c_str(), 400, -20, 20);
      _hBJetDiscrBySimpleSecondaryV[NpdfCounter] = fs->make<TH1F>(("BJetDiscrBySimpleSecondaryV_"+j.str()).c_str(), ("BJetDiscrBySimpleSecondaryV_"+j.str()).c_str(), 400, -20, 20);
      _hBJetDiscrByCombinedSecondaryV[NpdfCounter] = fs->make<TH1F>(("BJetDiscrByCombinedSecondaryV_"+j.str()).c_str(), ("BJetDiscrByCombinedSecondaryV_"+j.str()).c_str(), 400, -20, 20);
    }
    
    if (_FillTopologyHists) {
      if( ((_AnalyzeMuonForLeg1) && (_AnalyzeTauForLeg2)) || ((_AnalyzeMuonForLeg2) && (_AnalyzeTauForLeg1)) ) {
	_hMuonPtVsTauPt[NpdfCounter]                   = fs->make<TH2F>(("MuonPtVsTauPt_"+j.str()).c_str(),   ("MuonPtVsTauPt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hMuonTauDeltaR[NpdfCounter]                   = fs->make<TH1F>(("MuonTauDeltaR_"+j.str()).c_str(),   ("MuonTauDeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hMuonTauDeltaPtDivSumPt[NpdfCounter]          = fs->make<TH1F>(("MuonTauDeltaPtDivSumPt_"+j.str()).c_str(), ("MuonTauDeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hMuonMetMt[NpdfCounter]                       = fs->make<TH1F>(("MuonMetMt_"+j.str()).c_str(),       ("MuonMetMt_"+j.str()).c_str(), 100, 0, 500);
	_hTauMetMt[NpdfCounter]                        = fs->make<TH1F>(("TauMetMt_"+j.str()).c_str(),        ("TauMetMt_"+j.str()).c_str(), 100, 0, 500);
	_hMuonTauOSLS[NpdfCounter]                     = fs->make<TH1F>(("MuonTauOSLS_"+j.str()).c_str(),     ("MuonTauOSLS_"+j.str()).c_str(), 20, -10, 10);
	_hMuonTauCosDphi[NpdfCounter]                  = fs->make<TH1F>(("MuonTauCosDphi_"+j.str()).c_str(),  ("MuonTauCosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hMuonMetDeltaPhi[NpdfCounter]                 = fs->make<TH1F>(("MuonMetDeltaPhi_"+j.str()).c_str(), ("MuonMetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hTauMetDeltaPhi[NpdfCounter]                  = fs->make<TH1F>(("TauMetDeltaPhi_"+j.str()).c_str(), ("TauMetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hMuonMetDeltaPhiVsMuonTauCosDphi[NpdfCounter] = fs->make<TH2F>(("MuonMetDeltaPhiVsMuonTauCosDphi_"+j.str()).c_str(), ("MuonMetDeltaPhiVsMuonTauCosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
      }
      if( ((_AnalyzeElectronForLeg1) && (_AnalyzeTauForLeg2)) || ((_AnalyzeElectronForLeg2) && (_AnalyzeTauForLeg1)) ) {
	_hElectronPtVsTauPt[NpdfCounter]          = fs->make<TH2F>(("ElectronPtVsTauPt_"+j.str()).c_str(),          ("ElectronPtVsTauPt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hElectronTauDeltaR[NpdfCounter]          = fs->make<TH1F>(("ElectronTauDeltaR_"+j.str()).c_str(),          ("ElectronTauDeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hElectronTauDeltaPtDivSumPt[NpdfCounter] = fs->make<TH1F>(("ElectronTauDeltaPtDivSumPt_"+j.str()).c_str(), ("ElectronTauDeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hElectronMetMt[NpdfCounter]              = fs->make<TH1F>(("ElectronMetMt_"+j.str()).c_str(),              ("ElectronMetMt_"+j.str()).c_str(), 100, 0, 500);
	_hTauMetMt[NpdfCounter]                   = fs->make<TH1F>(("TauMetMt_"+j.str()).c_str(),                   ("TauMetMt_"+j.str()).c_str(), 100, 0, 500);
	_hElectronTauOSLS[NpdfCounter]            = fs->make<TH1F>(("ElectronTauOSLS_"+j.str()).c_str(),            ("ElectronTauOSLS_"+j.str()).c_str(), 20, -10, 10);
	_hElectronTauCosDphi[NpdfCounter]         = fs->make<TH1F>(("ElectronTauCosDphi_"+j.str()).c_str(),         ("ElectronTauCosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hElectronMetDeltaPhi[NpdfCounter]        = fs->make<TH1F>(("ElectronMetDeltaPhi_"+j.str()).c_str(),        ("ElectronMetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hTauMetDeltaPhi[NpdfCounter]             = fs->make<TH1F>(("TauMetDeltaPhi_"+j.str()).c_str(),             ("TauMetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hElectronMetDeltaPhiVsElectronTauCosDphi[NpdfCounter] = fs->make<TH2F>(("ElectronMetDeltaPhiVsElectronTauCosDphi_"+j.str()).c_str(), ("ElectronMetDeltaPhiVsElectronTauCosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
      }
      if( ((_AnalyzeMuonForLeg1) && (_AnalyzeElectronForLeg2)) || ((_AnalyzeMuonForLeg2) && (_AnalyzeElectronForLeg1)) ) {
	_hElectronPtVsMuonPt[NpdfCounter]          = fs->make<TH2F>(("ElectronPtVsMuonPt_"+j.str()).c_str(),          ("ElectronPtVsMuonPt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hElectronMuonDeltaR[NpdfCounter]          = fs->make<TH1F>(("ElectronMuonDeltaR_"+j.str()).c_str(),          ("ElectronMuonDeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hElectronMuonDeltaPtDivSumPt[NpdfCounter] = fs->make<TH1F>(("ElectronMuonDeltaPtDivSumPt_"+j.str()).c_str(), ("ElectronMuonDeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hElectronMetMt[NpdfCounter]               = fs->make<TH1F>(("ElectronMetMt_"+j.str()).c_str(),               ("ElectronMetMt_"+j.str()).c_str(), 100, 0, 500);
	_hMuonMetMt[NpdfCounter]                   = fs->make<TH1F>(("MuonMetMt_"+j.str()).c_str(),                   ("MuonMetMt_"+j.str()).c_str(), 100, 0, 500);
	_hElectronMuonOSLS[NpdfCounter]            = fs->make<TH1F>(("ElectronMuonOSLS_"+j.str()).c_str(),            ("ElectronMuonOSLS_"+j.str()).c_str(), 20, -10, 10);
	_hElectronMuonCosDphi[NpdfCounter]         = fs->make<TH1F>(("ElectronMuonCosDphi_"+j.str()).c_str(),         ("ElectronMuonCosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hElectronMetDeltaPhi[NpdfCounter]         = fs->make<TH1F>(("ElectronMetDeltaPhi_"+j.str()).c_str(),         ("ElectronMetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hMuonMetDeltaPhi[NpdfCounter]             = fs->make<TH1F>(("MuonMetDeltaPhi_"+j.str()).c_str(),             ("MuonMetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hElectronMetDeltaPhiVsElectronMuonCosDphi[NpdfCounter] = fs->make<TH2F>(("ElectronMetDeltaPhiVsElectronMuonCosDphi_"+j.str()).c_str(), ("ElectronMetDeltaPhiVsElectronMuonCosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
      }
      if( ((_AnalyzeTauForLeg1) && (_AnalyzeTauForLeg2)) ) {
	_hTau1PtVsTau2Pt[NpdfCounter]          = fs->make<TH2F>(("Tau1PtVsTau2Pt_"+j.str()).c_str(),          ("Tau1PtVsTau2Pt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hTau1Tau2DeltaR[NpdfCounter]          = fs->make<TH1F>(("Tau1Tau2DeltaR_"+j.str()).c_str(),          ("Tau1Tau2DeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hTau1Tau2DeltaPtDivSumPt[NpdfCounter] = fs->make<TH1F>(("Tau1Tau2DeltaPtDivSumPt_"+j.str()).c_str(), ("Tau1Tau2DeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hTau1MetMt[NpdfCounter]               = fs->make<TH1F>(("Tau1MetMt_"+j.str()).c_str(),               ("Tau1MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hTau2MetMt[NpdfCounter]               = fs->make<TH1F>(("Tau2MetMt_"+j.str()).c_str(),               ("Tau2MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hTau1Tau2OSLS[NpdfCounter]            = fs->make<TH1F>(("Tau1Tau2OSLS_"+j.str()).c_str(),            ("Tau1Tau2OSLS_"+j.str()).c_str(), 20, -10, 10);
	_hTau1Tau2CosDphi[NpdfCounter]         = fs->make<TH1F>(("Tau1Tau2CosDphi_"+j.str()).c_str(),         ("Tau1Tau2CosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hTau1MetDeltaPhi[NpdfCounter]         = fs->make<TH1F>(("Tau1MetDeltaPhi_"+j.str()).c_str(),         ("Tau1MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hTau2MetDeltaPhi[NpdfCounter]         = fs->make<TH1F>(("Tau2MetDeltaPhi_"+j.str()).c_str(),         ("Tau2MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hTau1MetDeltaPhiVsTau1Tau2CosDphi[NpdfCounter] = fs->make<TH2F>(("Tau1MetDeltaPhiVsTau1Tau2CosDphi_"+j.str()).c_str(), ("Tau1MetDeltaPhiVsTau1Tau2CosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
      }
      if( ((_AnalyzeMuonForLeg1) && (_AnalyzeMuonForLeg2)) ) {
	_hMuon1PtVsMuon2Pt[NpdfCounter]          = fs->make<TH2F>(("Muon1PtVsMuon2Pt_"+j.str()).c_str(),          ("Muon1PtVsMuon2Pt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hMuon1Muon2DeltaR[NpdfCounter]          = fs->make<TH1F>(("Muon1Muon2DeltaR_"+j.str()).c_str(),          ("Muon1Muon2DeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hMuon1Muon2DeltaPtDivSumPt[NpdfCounter] = fs->make<TH1F>(("Muon1Muon2DeltaPtDivSumPt_"+j.str()).c_str(), ("Muon1Muon2DeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hMuon1MetMt[NpdfCounter]                = fs->make<TH1F>(("Muon1MetMt_"+j.str()).c_str(),                ("Muon1MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hMuon2MetMt[NpdfCounter]                = fs->make<TH1F>(("Muon2MetMt_"+j.str()).c_str(),                ("Muon2MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hMuon1Muon2OSLS[NpdfCounter]            = fs->make<TH1F>(("Muon1Muon2OSLS_"+j.str()).c_str(),            ("Muon1Muon2OSLS_"+j.str()).c_str(), 20, -10, 10);
	_hMuon1Muon2CosDphi[NpdfCounter]         = fs->make<TH1F>(("Muon1Muon2CosDphi_"+j.str()).c_str(),         ("Muon1Muon2CosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hMuon1MetDeltaPhi[NpdfCounter]          = fs->make<TH1F>(("Muon1MetDeltaPhi_"+j.str()).c_str(),          ("Muon1MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hMuon2MetDeltaPhi[NpdfCounter]          = fs->make<TH1F>(("Muon2MetDeltaPhi_"+j.str()).c_str(),          ("Muon2MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hMuon1MetDeltaPhiVsMuon1Muon2CosDphi[NpdfCounter] = fs->make<TH2F>(("Muon1MetDeltaPhiVsMuon1Muon2CosDphi_"+j.str()).c_str(), ("Muon1MetDeltaPhiVsMuon1Muon2CosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
      }
      if( ((_AnalyzeElectronForLeg1) && (_AnalyzeElectronForLeg2)) ) {
	_hElectron1PtVsElectron2Pt[NpdfCounter]          = fs->make<TH2F>(("Electron1PtVsElectron2Pt_"+j.str()).c_str(), ("Electron1PtVsElectron2Pt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hElectron1Electron2DeltaR[NpdfCounter]          = fs->make<TH1F>(("Electron1Electron2DeltaR_"+j.str()).c_str(), ("Electron1Electron2DeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hElectron1Electron2DeltaPtDivSumPt[NpdfCounter] = fs->make<TH1F>(("Electron1Electron2DeltaPtDivSumPt_"+j.str()).c_str(), ("Electron1Electron2DeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hElectron1MetMt[NpdfCounter]            = fs->make<TH1F>(("Electron1MetMt_"+j.str()).c_str(),            ("Electron1MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hElectron2MetMt[NpdfCounter]            = fs->make<TH1F>(("Electron2MetMt_"+j.str()).c_str(),            ("Electron2MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hElectron1Electron2OSLS[NpdfCounter]    = fs->make<TH1F>(("Electron1Electron2OSLS_"+j.str()).c_str(),    ("Electron1Electron2OSLS_"+j.str()).c_str(), 20, -10, 10);
	_hElectron1Electron2CosDphi[NpdfCounter] = fs->make<TH1F>(("Electron1Electron2CosDphi_"+j.str()).c_str(), ("Electron1Electron2CosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hElectron1MetDeltaPhi[NpdfCounter]      = fs->make<TH1F>(("Electron1MetDeltaPhi_"+j.str()).c_str(),      ("Electron1MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hElectron2MetDeltaPhi[NpdfCounter]      = fs->make<TH1F>(("Electron2MetDeltaPhi_"+j.str()).c_str(),      ("Electron2MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hElectron1MetDeltaPhiVsElectron1Electron2CosDphi[NpdfCounter] = fs->make<TH2F>(("Electron1MetDeltaPhiVsElectron1Electron2CosDphi_"+j.str()).c_str(), ("Electron1MetDeltaPhiVsElectron1Electron2CosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
      }
      /*
	_hTauJetSumPtIso_SeedOS[NpdfCounter] = fs->make<TH1F>(("TauJetSumPtIso_SeedOS_"+j.str()).c_str(), ("TauJetSumPtIso_SeedOS_"+j.str()).c_str(), 100, 0, 50);
	_hTauJetSumPtIso_JetOS[NpdfCounter]  = fs->make<TH1F>(("TauJetSumPtIso_JetOS_"+j.str()).c_str(),  ("TauJetSumPtIso_JetOS_"+j.str()).c_str(), 100, 0, 50);
	_hTauJetSumPtIso_SeedLS[NpdfCounter] = fs->make<TH1F>(("TauJetSumPtIso_SeedLS_"+j.str()).c_str(), ("TauJetSumPtIso_SeedLS_"+j.str()).c_str(), 100, 0, 50);
	_hTauJetSumPtIso_JetLS[NpdfCounter]  = fs->make<TH1F>(("TauJetSumPtIso_JetLS_"+j.str()).c_str(),  ("TauJetSumPtIso_JetLS_"+j.str()).c_str(), 100, 0, 50);
      */
      _hNotReconstructableMass[NpdfCounter] = fs->make<TH1F>(("NotReconstructableMass_"+j.str()).c_str(), ("NotReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hReconstructableMass[NpdfCounter]    = fs->make<TH1F>(("ReconstructableMass_"+j.str()).c_str(),    ("ReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hPZeta[NpdfCounter]                  = fs->make<TH1F>(("PZeta_"+j.str()).c_str(),                  ("PZeta_"+j.str()).c_str(), 200, -100, 100);
      _hPZetaVis[NpdfCounter]               = fs->make<TH1F>(("PZetaVis_"+j.str()).c_str(),               ("PZetaVis_"+j.str()).c_str(), 100, 0, 100);
      _hZeta2D[NpdfCounter]                 = fs->make<TH2F>(("Zeta2D_"+j.str()).c_str(),                 ("Zeta2D_"+j.str()).c_str(), 100, 0, 100, 200, -100, 100);
      _hZeta1D[NpdfCounter]                 = fs->make<TH1F>(("Zeta1D_"+j.str()).c_str(),                 ("Zeta1D_"+j.str()).c_str(), 50, -100, 100);
      _hMet[NpdfCounter]                    = fs->make<TH1F>(("Met_"+j.str()).c_str(),                    ("Met_"+j.str()).c_str(), 100, 0, 1000);
    }
    j.str("");
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

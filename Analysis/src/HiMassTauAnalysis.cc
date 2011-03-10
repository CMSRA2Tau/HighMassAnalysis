// Authors: Andres Florez, Alfredo Gurrola, Eduardo Luiggi, Chi Nhan Nguyen

#include "HighMassAnalysis/Analysis/interface/HiMassTauAnalysis.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPEle.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPCutCodes.h"

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
  _DoRecoTauDiscrByH3x3OverP = iConfig.getParameter<bool>("DoRecoTauDiscrByH3x3OverP");
  _RecoTauH3x3OverP = iConfig.getParameter<double>("RecoTauH3x3OverP");
  _DoRecoTauDiscrByIsolation = iConfig.getParameter<bool>("DoRecoTauDiscrByIsolation");
  _UseRecoTauDiscrByIsolationFlag = iConfig.getParameter<bool>("UseRecoTauDiscrByIsolationFlag");
  _UseRecoTauIsoSumPtInsteadOfNiso = iConfig.getParameter<bool>("UseRecoTauIsoSumPtInsteadOfNiso");
  _UseRecoTauEllipseForEcalIso = iConfig.getParameter<bool>("UseRecoTauEllipseForEcalIso");
  _RecoTauEcalIsoRphiForEllipse = iConfig.getParameter<double>("RecoTauEcalIsoRphiForEllipse");
  _RecoTauEcalIsoRetaForEllipse = iConfig.getParameter<double>("RecoTauEcalIsoRetaForEllipse");
  _RecoTauTrackNisoMax = iConfig.getParameter<int>("RecoTauTrackNisoMax");
  _RecoTauEcalNisoMax = iConfig.getParameter<int>("RecoTauEcalNisoMax");
  _RecoTauTrackIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoTauTrackIsoSumPtMinCutValue");
  _RecoTauTrackIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoTauTrackIsoSumPtMaxCutValue");
  _RecoTauEcalIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoTauEcalIsoSumPtMinCutValue");
  _RecoTauEcalIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoTauEcalIsoSumPtMaxCutValue");
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
  _RecoMuonTrackIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoMuonTrackIsoSumPtMinCutValue");
  _RecoMuonTrackIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoMuonTrackIsoSumPtMaxCutValue");
  _RecoMuonEcalIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoMuonEcalIsoSumPtMinCutValue");
  _RecoMuonEcalIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoMuonEcalIsoSumPtMaxCutValue");
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
  _UseHeepInfo = iConfig.getParameter<bool>("UseHeepInfo");
  _RecoElectronEtaCut = iConfig.getParameter<double>("RecoElectronEtaCut");
  _RecoElectronPtMinCut = iConfig.getParameter<double>("RecoElectronPtMinCut");
  _RecoElectronPtMaxCut = iConfig.getParameter<double>("RecoElectronPtMaxCut");
  _DoRecoElectronDiscrByTrackIsolation = iConfig.getParameter<bool>("DoRecoElectronDiscrByTrackIsolation");
  _RecoElectronTrackIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoElectronTrackIsoSumPtMaxCutValue");
  _RecoElectronTrackIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoElectronTrackIsoSumPtMinCutValue");
  _RecoElectronTrackIsoDeltaRCone = iConfig.getParameter<double>("RecoElectronTrackIsoDeltaRCone");
  _RecoElectronTrackIsoTrkThreshold = iConfig.getParameter<double>("RecoElectronTrackIsoTrkThreshold");
  _DoRecoElectronDiscrByEcalIsolation = iConfig.getParameter<bool>("DoRecoElectronDiscrByEcalIsolation");
  _RecoElectronEcalIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoElectronEcalIsoSumPtMaxCutValue");
  _RecoElectronEcalIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoElectronEcalIsoSumPtMinCutValue");
  _RecoElectronEcalIsoDeltaRCone = iConfig.getParameter<double>("RecoElectronEcalIsoDeltaRCone");
  _RecoElectronEcalIsoRecHitThreshold = iConfig.getParameter<double>("RecoElectronEcalIsoRecHitThreshold");
  _DoRecoElectronDiscrByIp = iConfig.getParameter<bool>("DoRecoElectronDiscrByIp");
  _RecoElectronIpCut = iConfig.getParameter<double>("RecoElectronIpCut");
  _DoRecoElectronDiscrByEoverP = iConfig.getParameter<bool>("DoRecoElectronDiscrByEoverP");
  _RecoElectronEoverPMax = iConfig.getParameter<double>("RecoElectronEoverPMax");
  _RecoElectronEoverPMin = iConfig.getParameter<double>("RecoElectronEoverPMin");
  _DoRecoElectronDiscrByHoverEm = iConfig.getParameter<bool>("DoRecoElectronDiscrByHoverEm");
  _RecoElectronHoverEmCut = iConfig.getParameter<double>("RecoElectronHoverEmCut");
  _DoRecoElectronDiscrBySigmaIEtaIEta = iConfig.getParameter<bool>("DoRecoElectronDiscrBySigmaIEtaIEta");
  _RecoElectronSigmaIEtaIEta = iConfig.getParameter<double>("RecoElectronSigmaIEtaIEta");
  _DoRecoElectronDiscrByDEtaIn = iConfig.getParameter<bool>("DoRecoElectronDiscrByDEtaIn");
  _RecoElectronEEDEtaIn = iConfig.getParameter<double>("RecoElectronEEDEtaIn");
  _RecoElectronEBDEtaIn = iConfig.getParameter<double>("RecoElectronEBDEtaIn");
  _DoRecoElectronDiscrByDPhiIn = iConfig.getParameter<bool>("DoRecoElectronDiscrByDPhiIn");
  _RecoElectronEEDPhiIn = iConfig.getParameter<double>("RecoElectronEEDPhiIn");
  _RecoElectronEBDPhiIn = iConfig.getParameter<double>("RecoElectronEBDPhiIn");
  _DoRecoElectronDiscrBySCE2by5Over5by5 = iConfig.getParameter<bool>("DoRecoElectronDiscrBySCE2by5Over5by5");
  _RecoElectronEBscE1by5Over5by5 = iConfig.getParameter<double>("RecoElectronEBscE1by5Over5by5");
  _RecoElectronEBscE2by5Over5by5 = iConfig.getParameter<double>("RecoElectronEBscE2by5Over5by5");
  _DoRecoElectronDiscrByMissingHits = iConfig.getParameter<bool>("DoRecoElectronDiscrByMissingHits");
  _RecoElectronMissingHits = iConfig.getParameter<int>("RecoElectronMissingHits");
  _DoRecoElectronDiscrByEcalDrivenSeed = iConfig.getParameter<bool>("DoRecoElectronDiscrByEcalDrivenSeed");
  _DoRecoElectronDiscrByTrackerDrivenSeed = iConfig.getParameter<bool>("DoRecoElectronDiscrByTrackerDrivenSeed");

/*
  //-----Reco ParticleFlow Inputs
  _RecoParticleFlowSource = iConfig.getParameter<InputTag>("RecoParticleFlowSource");
  _UsePFlowBasedIsolationInsteadOfStandard = iConfig.getParameter<bool>("UsePFlowBasedIsolationInsteadOfStandard");
*/

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
  _ApplyJetBTagging = iConfig.getParameter<bool>("ApplyJetBTagging");
  _JetBTaggingTCHEcut = iConfig.getParameter<double>("JetBTaggingTCHEcut");
  _DoDiscrByFirstLeadingJet = iConfig.getParameter<bool>("DoDiscrByFirstLeadingJet");
  _RecoFirstLeadingJetPt = iConfig.getParameter<double>("RecoFirstLeadingJetPt");
  _RecoFirstLeadingJetEtaMinCut = iConfig.getParameter<double>("RecoFirstLeadingJetEtaMinCut");
  _RecoFirstLeadingJetEtaMaxCut = iConfig.getParameter<double>("RecoFirstLeadingJetEtaMaxCut");
  _DoDiscrBySecondLeadingJet = iConfig.getParameter<bool>("DoDiscrBySecondLeadingJet");
  _RecoSecondLeadingJetPt = iConfig.getParameter<double>("RecoSecondLeadingJetPt");
  _RecoSecondLeadingJetEtaMinCut = iConfig.getParameter<double>("RecoSecondLeadingJetEtaMinCut");
  _RecoSecondLeadingJetEtaMaxCut = iConfig.getParameter<double>("RecoSecondLeadingJetEtaMaxCut");

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
  _CalculateMetUsingOnlyLeg1AndLeg2 = iConfig.getParameter<bool>("CalculateMetUsingOnlyLeg1AndLeg2");
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
  _DoDiTauDiscrByDeltaPt = iConfig.getParameter<bool>("DoDiTauDiscrByDeltaPt");
  _DiTauDeltaPtMinCutValue = iConfig.getParameter<double>("DiTauDeltaPtMinCutValue");
  _DiTauDeltaPtMaxCutValue = iConfig.getParameter<double>("DiTauDeltaPtMaxCutValue");
  _DoDiscrByLeg1MetDphi = iConfig.getParameter<bool>("DoDiscrByLeg1MetDphi");
  _Leg1MetDphiMinCut = iConfig.getParameter<double>("Leg1MetDphiMinCut");
  _Leg1MetDphiMaxCut = iConfig.getParameter<double>("Leg1MetDphiMaxCut");
  _DoDiscrByLeg2MetDphi = iConfig.getParameter<bool>("DoDiscrByLeg2MetDphi");
  _Leg2MetDphiMinCut = iConfig.getParameter<double>("Leg2MetDphiMinCut");
  _Leg2MetDphiMaxCut = iConfig.getParameter<double>("Leg2MetDphiMaxCut");
  _DoTauDiscrByIsZeeCut = iConfig.getParameter<bool>("DoTauDiscrByIsZeeCut");

  //-----SUSY Specific Topology Inputs
  _DoSUSYDiscrByMHT = iConfig.getParameter<bool>("DoSUSYDiscrByMHT");
  _MhtCut = iConfig.getParameter<double>("MhtCut");
  _DoSUSYDiscrByR1 = iConfig.getParameter<bool>("DoSUSYDiscrByR1");
  _R1MinCut = iConfig.getParameter<double>("R1MinCut");
  _R1MaxCut = iConfig.getParameter<double>("R1MaxCut");
  _DoSUSYDiscrByR2 = iConfig.getParameter<bool>("DoSUSYDiscrByR2");
  _R2MinCut = iConfig.getParameter<double>("R2MinCut");
  _R2MaxCut = iConfig.getParameter<double>("R2MaxCut");
  _DoSUSYDiscrByAlpha = iConfig.getParameter<bool>("DoSUSYDiscrByAlpha");
  _AlphaMinCut = iConfig.getParameter<double>("AlphaMinCut");
  _AlphaMaxCut = iConfig.getParameter<double>("AlphaMaxCut");
  _DoSUSYDiscrByDphi1 = iConfig.getParameter<bool>("DoSUSYDiscrByDphi1");
  _Dphi1MinCut = iConfig.getParameter<double>("Dphi1MinCut");
  _Dphi1MaxCut = iConfig.getParameter<double>("Dphi1MaxCut");
  _DoSUSYDiscrByDphi2 = iConfig.getParameter<bool>("DoSUSYDiscrByDphi2");
  _Dphi2MinCut = iConfig.getParameter<double>("Dphi2MinCut");
  _Dphi2MaxCut = iConfig.getParameter<double>("Dphi2MaxCut");

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
  _RecoFirstLeadingJetNmin = iConfig.getParameter<int>("RecoFirstLeadingJetNmin");
  _RecoSecondLeadingJetNmin = iConfig.getParameter<int>("RecoSecondLeadingJetNmin");
  _SusyCombinationsNmin = iConfig.getParameter<int>("SusyCombinationsNmin");
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
  _CalculatePUSystematics = iConfig.getParameter<bool>("CalculatePUSystematics");
  _ApplyMuonTriggerScaleFactors = iConfig.getParameter<bool>("ApplyMuonTriggerScaleFactors");
  _ApplyElectronTriggerScaleFactors = iConfig.getParameter<bool>("ApplyElectronTriggerScaleFactors");
  _ApplyTauTriggerScaleFactors = iConfig.getParameter<bool>("ApplyTauTriggerScaleFactors");

  //---- Initialize everything needed to calculate systematics due to initial state gluon radiation
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

  //---- Initialize everything needed to apply Data/MC scale factors for trigger
  if(_ApplyMuonTriggerScaleFactors) {
    // default pt bin edges
    std::vector<double> defMuonTrigPtEdges;
    defMuonTrigPtEdges.push_back(0.);
    defMuonTrigPtEdges.push_back(999999.);
    MuonTrigPtBinEdges_ = iConfig.getUntrackedParameter<std::vector<double> > ("MuonTrigPtBinEdges",defMuonTrigPtEdges);
    unsigned int ninputs_expected = MuonTrigPtBinEdges_.size()-1;
    // default scale factors
    std::vector<double> defMuonTrigPtWeights;
    defMuonTrigPtWeights.push_back(1.);
    MuonTrigptWeights_ = iConfig.getUntrackedParameter<std::vector<double> > ("MuonTrigPtWeights",defMuonTrigPtWeights);
    if (MuonTrigptWeights_.size()==1 && ninputs_expected>1) {
      for (unsigned int i=1; i<ninputs_expected; i++){ MuonTrigptWeights_.push_back(MuonTrigptWeights_[0]);}
    }
    // default eta bin edges
    std::vector<double> defMuonTrigEtaEdges;
    defMuonTrigEtaEdges.push_back(-999999.);
    defMuonTrigEtaEdges.push_back(999999.);
    MuonTrigEtaBinEdges_ = iConfig.getUntrackedParameter<std::vector<double> > ("MuonTrigEtaBinEdges",defMuonTrigEtaEdges);
    ninputs_expected = MuonTrigEtaBinEdges_.size()-1;
    // default scale factors
    std::vector<double> defMuonTrigEtaWeights;
    defMuonTrigEtaWeights.push_back(1.);
    MuonTrigetaWeights_ = iConfig.getUntrackedParameter<std::vector<double> > ("MuonTrigEtaWeights",defMuonTrigEtaWeights);
    if (MuonTrigetaWeights_.size()==1 && ninputs_expected>1) {
      for (unsigned int i=1; i<ninputs_expected; i++){ MuonTrigetaWeights_.push_back(MuonTrigetaWeights_[0]);}
    }
  }
  if(_ApplyElectronTriggerScaleFactors) {
    // default pt bin edges
    std::vector<double> defElectronTrigPtEdges;
    defElectronTrigPtEdges.push_back(0.);
    defElectronTrigPtEdges.push_back(999999.);   
    ElectronTrigPtBinEdges_ = iConfig.getUntrackedParameter<std::vector<double> > ("ElectronTrigPtBinEdges",defElectronTrigPtEdges);
    unsigned int ninputs_expected = ElectronTrigPtBinEdges_.size()-1;
    // default scale factors
    std::vector<double> defElectronTrigPtWeights;
    defElectronTrigPtWeights.push_back(1.);
    ElectronTrigptWeights_ = iConfig.getUntrackedParameter<std::vector<double> > ("ElectronTrigPtWeights",defElectronTrigPtWeights);
    if (ElectronTrigptWeights_.size()==1 && ninputs_expected>1) {
      for (unsigned int i=1; i<ninputs_expected; i++){ ElectronTrigptWeights_.push_back(ElectronTrigptWeights_[0]);}
    }
    // default eta bin edges
    std::vector<double> defElectronTrigEtaEdges;
    defElectronTrigEtaEdges.push_back(-999999.);
    defElectronTrigEtaEdges.push_back(999999.);
    ElectronTrigEtaBinEdges_ = iConfig.getUntrackedParameter<std::vector<double> > ("ElectronTrigEtaBinEdges",defElectronTrigEtaEdges);
    ninputs_expected = ElectronTrigEtaBinEdges_.size()-1;
    // default scale factors
    std::vector<double> defElectronTrigEtaWeights;
    defElectronTrigEtaWeights.push_back(1.);
    ElectronTrigetaWeights_ = iConfig.getUntrackedParameter<std::vector<double> > ("ElectronTrigEtaWeights",defElectronTrigEtaWeights);
    if (ElectronTrigetaWeights_.size()==1 && ninputs_expected>1) {
      for (unsigned int i=1; i<ninputs_expected; i++){ ElectronTrigetaWeights_.push_back(ElectronTrigetaWeights_[0]);}
    }
  }
  if(_ApplyTauTriggerScaleFactors) {
    // default pt bin edges
    std::vector<double> defTauTrigPtEdges;
    defTauTrigPtEdges.push_back(0.);
    defTauTrigPtEdges.push_back(999999.);   
    TauTrigPtBinEdges_ = iConfig.getUntrackedParameter<std::vector<double> > ("TauTrigPtBinEdges",defTauTrigPtEdges);
    unsigned int ninputs_expected = TauTrigPtBinEdges_.size()-1;
    // default scale factors
    std::vector<double> defTauTrigPtWeights;
    defTauTrigPtWeights.push_back(1.);
    TauTrigptWeights_ = iConfig.getUntrackedParameter<std::vector<double> > ("TauTrigPtWeights",defTauTrigPtWeights);
    if (TauTrigptWeights_.size()==1 && ninputs_expected>1) {
      for (unsigned int i=1; i<ninputs_expected; i++){ TauTrigptWeights_.push_back(TauTrigptWeights_[0]);}
    }
    // default eta bin edges
    std::vector<double> defTauTrigEtaEdges;
    defTauTrigEtaEdges.push_back(-999999.);
    defTauTrigEtaEdges.push_back(999999.);
    TauTrigEtaBinEdges_ = iConfig.getUntrackedParameter<std::vector<double> > ("TauTrigEtaBinEdges",defTauTrigEtaEdges);
    ninputs_expected = TauTrigEtaBinEdges_.size()-1;
    // default scale factors
    std::vector<double> defTauTrigEtaWeights;
    defTauTrigEtaWeights.push_back(1.);
    TauTrigetaWeights_ = iConfig.getUntrackedParameter<std::vector<double> > ("TauTrigEtaWeights",defTauTrigEtaWeights);
    if (TauTrigetaWeights_.size()==1 && ninputs_expected>1) {
      for (unsigned int i=1; i<ninputs_expected; i++){ TauTrigetaWeights_.push_back(TauTrigetaWeights_[0]);}
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
}

// Set branches for the ntuple
void  HiMassTauAnalysis::setupBranches() {
  _HMTTree = new TTree(_NtupleTreeName.c_str(), "HiMassDiTau Tree");
  _HMTTree->Branch("tauTrkIsoPat",&_tauTrkIsoPat);
  _HMTTree->Branch("muTrkIsoPat", &_muTrkIsoPat);
  _HMTTree->Branch("mEt", &_mEt);
  _HMTTree->Branch("zeta", &_zeta);
  _HMTTree->Branch("muTauMetMass",&_muTauMetMass);
  _HMTTree->Branch("muMetMass",&_muMetMass);
  _HMTTree->Branch("OSLS",  &_OSLS);
  _HMTTree->Branch("jetPt",  &_jetPt);
  _HMTTree->Branch("jetEt",  &_jetEt);
  _HMTTree->Branch("jetEta", &_jetEta); 
  _HMTTree->Branch("jetPhi", &_jetPhi);
  _HMTTree->Branch("jetEmFraction", &_jetEmFraction);
  _HMTTree->Branch("bJetDiscrByTrackCounting",&_bJetDiscrByTrackCounting);
  _HMTTree->Branch("bJetDiscrBySimpleSecondaryV",&_bJetDiscrBySimpleSecondaryV);
  _HMTTree->Branch("bJetDiscrByCombinedSecondaryV",&_bJetDiscrByCombinedSecondaryV);

}

// Initialize the vectors for the ntuple
void HiMassTauAnalysis::initializeVectors(){

  _jetPt  = NULL;
  _jetEt  = NULL;
  _jetPhi = NULL;
  _jetEta = NULL;
  _jetEmFraction = NULL;
  _bJetDiscrByTrackCounting = NULL;	  
  _bJetDiscrBySimpleSecondaryV = NULL;   
  _bJetDiscrByCombinedSecondaryV = NULL; 
  _tauTrkIsoPat  = NULL;
  _muTrkIsoPat  = NULL;
  _mEt  = NULL;
  _zeta  = NULL;
  _muTauMetMass  = NULL;
  _muMetMass  = NULL;
  _OSLS  = NULL;

}

// clear the vectors vefore each event
void HiMassTauAnalysis::clearVectors(){

  _jetPt->clear();
  _jetEt->clear();
  _jetPhi->clear();
  _jetEta->clear();
  _jetEmFraction->clear();
  _bJetDiscrByTrackCounting->clear();
  _bJetDiscrBySimpleSecondaryV->clear();
  _bJetDiscrByCombinedSecondaryV->clear();
  _tauTrkIsoPat->clear();
  _muTrkIsoPat->clear();
  _mEt->clear();
  _zeta->clear();
  _muTauMetMass->clear();
  _muMetMass->clear();
  _OSLS->clear();

}

// ------------ method called to for each event  ------------
void HiMassTauAnalysis::analyze(const Event& iEvent, const EventSetup& iSetup) {

  //------Number of events analyzed (denominator)
  _totalEvents++;

  //-----Get weights for the calculation of pdf systematic uncertainties for the denominator
  pdfWeightVector.clear();
  if(_CalculatePdfSystematicUncertanties) {
    for (unsigned int i=0; i<pdfWeightTags_.size(); ++i) {
      //std::cout << "pdf tag = " << pdfWeightTags_[i] << std::endl;
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
        //std::cout << "weight member " << j << " = " << weights[j] << std::endl;
      }
    }
  } else{ pdfWeightVector.push_back(1); } 

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
  } else { isrgluon_weight = 1; }

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
  } else { isrgamma_weight = 1; }

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
  } else { fsr_weight = 1; }

  trig_weight = 1;

  // deltas for recalculation of MET (used when studying systematics)
  deltaForMEx = 0;
  deltaForMEy = 0;

  if(_CalculatePUSystematics) {
    float unirand = (float)(gRandom->Uniform(0,1));
    if(unirand<=0.225) {
      double rndphi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
      double rndmet = fabs(CLHEP::RandGauss::shoot(0, 3.5));
      deltaForMEx = deltaForMEx + (1.2 * rndmet * sqrt(1.0) * cos(rndphi));
      deltaForMEy = deltaForMEy + (1.2 * rndmet * sqrt(1.0) * sin(rndphi));
    } else if(unirand<=0.512) {
      double rndphi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
      double rndmet = fabs(CLHEP::RandGauss::shoot(0, 3.5));
      deltaForMEx = deltaForMEx + (1.2 * rndmet * sqrt(2.0) * cos(rndphi));
      deltaForMEy = deltaForMEy + (1.2 * rndmet * sqrt(2.0) * sin(rndphi));
    } else if(unirand<=0.780) {
      double rndphi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
      double rndmet = fabs(CLHEP::RandGauss::shoot(0, 3.5));
      deltaForMEx = deltaForMEx + (1.2 * rndmet * sqrt(3.0) * cos(rndphi));
      deltaForMEy = deltaForMEy + (1.2 * rndmet * sqrt(3.0) * sin(rndphi));
    } else if(unirand<=0.919) {
      double rndphi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
      double rndmet = fabs(CLHEP::RandGauss::shoot(0, 3.5));
      deltaForMEx = deltaForMEx + (1.2 * rndmet * sqrt(4.0) * cos(rndphi));
      deltaForMEy = deltaForMEy + (1.2 * rndmet * sqrt(4.0) * sin(rndphi));
    } else if(unirand<=0.976) {
      double rndphi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
      double rndmet = fabs(CLHEP::RandGauss::shoot(0, 3.5));
      deltaForMEx = deltaForMEx + (1.2 * rndmet * sqrt(5.0) * cos(rndphi));
      deltaForMEy = deltaForMEy + (1.2 * rndmet * sqrt(5.0) * sin(rndphi));
    } else if(unirand<=0.986) {
      double rndphi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
      double rndmet = fabs(CLHEP::RandGauss::shoot(0, 3.5));
      deltaForMEx = deltaForMEx + (1.2 * rndmet * sqrt(6.0) * cos(rndphi));
      deltaForMEy = deltaForMEy + (1.2 * rndmet * sqrt(6.0) * sin(rndphi));
    } else {
      double rndphi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
      double rndmet = fabs(CLHEP::RandGauss::shoot(0, 3.5));
      deltaForMEx = deltaForMEx + (1.2 * rndmet * sqrt(7.0) * cos(rndphi));
      deltaForMEy = deltaForMEy + (1.2 * rndmet * sqrt(7.0) * sin(rndphi));
    }
  }

  //-----Smearing momentum and position for systematic uncertanties and calculation of MET deltas
  smearedMuonMomentumVector.clear();
  smearedMuonPtEtaPhiMVector.clear();
  if(((_AnalyzeMuonForLeg1) || (_AnalyzeMuonForLeg2))) {
    if(_SmearTheMuon) {
      for(pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end();++patMuon) {
        smearedMuonMomentumVector.push_back(SmearLightLepton(*patMuon).first);
        smearedMuonPtEtaPhiMVector.push_back(SmearLightLepton(*patMuon).second);
        deltaForMEx = deltaForMEx + patMuon->px() - SmearLightLepton(*patMuon).first.px();
        deltaForMEy = deltaForMEy + patMuon->py() - SmearLightLepton(*patMuon).first.py();
      }
    } else {
      for(pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end();++patMuon) {
        smearedMuonMomentumVector.push_back(patMuon->p4());
        math::PtEtaPhiMLorentzVector thePtEtaPhiMVector(patMuon->pt(), patMuon->eta(), patMuon->phi(), patMuon->mass());
        smearedMuonPtEtaPhiMVector.push_back(thePtEtaPhiMVector);
      }
    }
  }
  smearedElectronMomentumVector.clear();
  smearedElectronPtEtaPhiMVector.clear();
  if(((_AnalyzeElectronForLeg1) || (_AnalyzeElectronForLeg2))) {
    if(_SmearTheElectron) {
      for(pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end();++patElectron) {
        smearedElectronMomentumVector.push_back(SmearLightLepton(*patElectron).first);
        smearedElectronPtEtaPhiMVector.push_back(SmearLightLepton(*patElectron).second);
        deltaForMEx = deltaForMEx + patElectron->px() - SmearLightLepton(*patElectron).first.px();
        deltaForMEy = deltaForMEy + patElectron->py() - SmearLightLepton(*patElectron).first.py();
      }
    } else {
      for(pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end();++patElectron) {
        smearedElectronMomentumVector.push_back(patElectron->p4());
        if(_UseHeepInfo) {
          heep::Ele theHeepElec(*patElectron);
          math::PtEtaPhiMLorentzVector thePtEtaPhiMVector(theHeepElec.et(), theHeepElec.scEta(), patElectron->phi(), patElectron->mass());
          smearedElectronPtEtaPhiMVector.push_back(thePtEtaPhiMVector);
        } else {
          math::PtEtaPhiMLorentzVector thePtEtaPhiMVector(patElectron->pt(), patElectron->eta(), patElectron->phi(), patElectron->mass());
          smearedElectronPtEtaPhiMVector.push_back(thePtEtaPhiMVector);
        }
      }
    }
  }
  smearedTauMomentumVector.clear();
  smearedTauPtEtaPhiMVector.clear();
  if(((_AnalyzeTauForLeg1) || (_AnalyzeTauForLeg2))) {
    if(_SmearTheTau) {
      for(pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end();++patTau) {
        smearedTauMomentumVector.push_back(SmearTau(*patTau).first);
        smearedTauPtEtaPhiMVector.push_back(SmearTau(*patTau).second);
        deltaForMEx = deltaForMEx + patTau->px() - SmearTau(*patTau).first.px();
        deltaForMEy = deltaForMEy + patTau->py() - SmearTau(*patTau).first.py();
      }
    } else {
      for(pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end();++patTau) {
        smearedTauMomentumVector.push_back(patTau->p4());
        math::PtEtaPhiMLorentzVector thePtEtaPhiMVector(patTau->pt(), patTau->eta(), patTau->phi(), patTau->mass());
        smearedTauPtEtaPhiMVector.push_back(thePtEtaPhiMVector);
      }
    }
  }
  smearedJetMomentumVector.clear();
  smearedJetPtEtaPhiMVector.clear();
  if(_SmearTheJet) {
    for ( pat::JetCollection::const_iterator patJet = _patJets->begin(); patJet != _patJets->end(); ++patJet ) {
      smearedJetMomentumVector.push_back(SmearJet(*patJet).first);
      smearedJetPtEtaPhiMVector.push_back(SmearJet(*patJet).second);
      deltaForMEx = deltaForMEx + patJet->px() - SmearJet(*patJet).first.px();
      deltaForMEy = deltaForMEy + patJet->py() - SmearJet(*patJet).first.py();
    }
  } else {
    for ( pat::JetCollection::const_iterator patJet = _patJets->begin(); patJet != _patJets->end(); ++patJet ) {
      if(_UseCorrectedJet) {
        smearedJetMomentumVector.push_back(patJet->p4());
        math::PtEtaPhiMLorentzVector thePtEtaPhiMVector(patJet->pt(), patJet->eta(), patJet->phi(), patJet->mass());
        smearedJetPtEtaPhiMVector.push_back(thePtEtaPhiMVector);
      } else {
        smearedJetMomentumVector.push_back(patJet->correctedJet("raw","").p4());
        math::PtEtaPhiMLorentzVector thePtEtaPhiMVector(patJet->correctedJet("raw","").pt(), patJet->correctedJet("raw","").eta(), patJet->correctedJet("raw","").phi(), patJet->correctedJet("raw","").mass());
        smearedJetPtEtaPhiMVector.push_back(thePtEtaPhiMVector);
      }
    }
  }
  double temppx = (*(_patMETs->begin())).px() + deltaForMEx;
  double temppy = (*(_patMETs->begin())).py() + deltaForMEy;
  double temppz = (*(_patMETs->begin())).pz();
  double temppt = TMath::Sqrt((temppx*temppx) + (temppy*temppy));
  reco::Candidate::LorentzVector theTempMETVector(temppx,temppy,temppz,temppt);
  theMETVector = theTempMETVector;

  if(_totalEvents == 1) {
    bookHistograms();
    if (_DoProduceNtuple) {
      initializeVectors();
      setupBranches();
    }
  }

  //------Number of events analyzed (denominator)
  for(unsigned int NpdfID = 0; NpdfID < pdfWeightVector.size();  NpdfID++){ 
//    _hEvents[NpdfID]->Fill(0.0,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
    _hEvents[NpdfID]->Fill(0.0);
  }

  //------Get the event flags (did the event pass the cuts?)
  getEventFlags(iEvent);

  if (!passEventSelectionSequence()) return;  

  isrgluon_weight = isrgluon_weight * trig_weight;

  //------Number of events passing cuts (numerator)
  _totalEventsPassingCuts++;
  for(unsigned int NpdfID = 0; NpdfID < pdfWeightVector.size();  NpdfID++){   
    _hEvents[NpdfID]->Fill(1.0,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
  }

  //-----Get weights for the calculation of pdf systematic uncertainties for the numerator
  if(_CalculatePdfSystematicUncertanties) {
    for (unsigned int i=0; i<pdfWeightTags_.size(); ++i) {
      //std::cout << "pdf tag = " << pdfWeightTags_[i] << std::endl;
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
        //std::cout << "weight member " << j << " = " << weights[j] << std::endl;
      }
    }
  }

  if (_DoProduceNtuple){
    fillNtuple();
    clearVectors();
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
  for (unsigned int i=0;i<_EventSelectionSequence.size();i++) { _mapSelectionAlgoID[_EventSelectionSequence[i]] = i; }
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
  int nTriggersSatisfied = 0;
  if(passRecoTriggerCuts(iEvent)) {nTriggersSatisfied++;}
  if (nTriggersSatisfied>=_RecoTriggersNmin) _EventFlag[_mapSelectionAlgoID["RecoTriggersNmin"]] = true;

  // ------Number of Good Vertices
  int nGoodVertices = 0;
  for ( reco::VertexCollection::const_iterator primaryVertex = _primaryEventVertexCollection->begin();
        primaryVertex != _primaryEventVertexCollection->end(); ++primaryVertex ) {
    if (!passRecoVertexCuts(*primaryVertex)) continue;
    nGoodVertices++;
  }
  if (nGoodVertices>=_RecoVertexNmin) _EventFlag[_mapSelectionAlgoID["RecoVertexNmin"]] = true;
  if (nGoodVertices<=_RecoVertexNmax) _EventFlag[_mapSelectionAlgoID["RecoVertexNmax"]] = true;

  //------Used only if attempting to recalculate MET using only leg1 and leg2 momenta
  double temppx = 0;
  double temppy = 0;
  std::vector<int> usedMuons;
  std::vector<int> usedElectrons;
  std::vector<int> usedTaus;
  usedMuons.clear();
  usedElectrons.clear();
  usedTaus.clear();
  if((_AnalyzeMuonForLeg1) || (_AnalyzeMuonForLeg2)) {
    int theNumberOfMuons = 0;
    double maxptfortrigscalefactor = 0;
    double etafortrigscalefactor = -6;
    for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
      usedMuons.push_back(0);
      theNumberOfMuons++;
      if (passRecoMuonCuts((*patMuon),theNumberOfMuons-1)) {
        if(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt() > maxptfortrigscalefactor) {
          maxptfortrigscalefactor = smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt();
          etafortrigscalefactor = smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).eta();
        }
      }
    }
    if(_ApplyMuonTriggerScaleFactors) {
      double trig_ptweight = MuonTrigptWeights_[MuonTrigptWeights_.size()-1];
      unsigned int nbins = MuonTrigPtBinEdges_.size()-1;
      if (maxptfortrigscalefactor>MuonTrigPtBinEdges_[0] && maxptfortrigscalefactor<MuonTrigPtBinEdges_[nbins]) {
        bool foundCorrectBin = false;
        for (unsigned int jtbin=1; jtbin<=nbins; ++jtbin) {
          if(foundCorrectBin) continue;
          if (maxptfortrigscalefactor>MuonTrigPtBinEdges_[jtbin]) continue;
          trig_ptweight = MuonTrigptWeights_[jtbin-1];
          foundCorrectBin = true;
        }  
      }  
      double trig_etaweight = MuonTrigetaWeights_[MuonTrigetaWeights_.size()-1];
      nbins = MuonTrigEtaBinEdges_.size()-1;
      if (etafortrigscalefactor>MuonTrigEtaBinEdges_[0] && etafortrigscalefactor<MuonTrigEtaBinEdges_[nbins]) {
        bool foundCorrectBin = false;
        for (unsigned int jtbin=1; jtbin<=nbins; ++jtbin) {
          if(foundCorrectBin) continue;
          if (etafortrigscalefactor>MuonTrigEtaBinEdges_[jtbin]) continue;
          trig_etaweight = MuonTrigetaWeights_[jtbin-1];
          foundCorrectBin = true;
        }
      }
      trig_weight = trig_weight * trig_ptweight * trig_etaweight;
    }
  }
  if((_AnalyzeElectronForLeg1) || (_AnalyzeElectronForLeg2)) {
    int theNumberOfElectrons = 0;
    double maxptfortrigscalefactor = 0;
    double etafortrigscalefactor = -6;
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron ) {
      usedElectrons.push_back(0);
      theNumberOfElectrons++;
      if (passRecoElectronCuts((*patElectron),theNumberOfElectrons-1)) {
        if(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt() > maxptfortrigscalefactor) {
          maxptfortrigscalefactor = smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt();
          etafortrigscalefactor = smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).eta();
        }
      }
    }
    if(_ApplyElectronTriggerScaleFactors) {
      double trig_ptweight = ElectronTrigptWeights_[ElectronTrigptWeights_.size()-1];
      unsigned int nbins = ElectronTrigPtBinEdges_.size()-1;
      if (maxptfortrigscalefactor>ElectronTrigPtBinEdges_[0] && maxptfortrigscalefactor<ElectronTrigPtBinEdges_[nbins]) {
        bool foundCorrectBin = false;
        for (unsigned int jtbin=1; jtbin<=nbins; ++jtbin) {
          if(foundCorrectBin) continue;
          if (maxptfortrigscalefactor>ElectronTrigPtBinEdges_[jtbin]) continue;
          trig_ptweight = ElectronTrigptWeights_[jtbin-1];
          foundCorrectBin = true;
        }
      }
      double trig_etaweight = ElectronTrigetaWeights_[ElectronTrigetaWeights_.size()-1];
      nbins = ElectronTrigEtaBinEdges_.size()-1;
      if (etafortrigscalefactor>ElectronTrigEtaBinEdges_[0] && etafortrigscalefactor<ElectronTrigEtaBinEdges_[nbins]) {
        bool foundCorrectBin = false;
        for (unsigned int jtbin=1; jtbin<=nbins; ++jtbin) {
          if(foundCorrectBin) continue;
          if (etafortrigscalefactor>ElectronTrigEtaBinEdges_[jtbin]) continue;
          trig_etaweight = ElectronTrigetaWeights_[jtbin-1];
          foundCorrectBin = true;
        }
      }  
      trig_weight = trig_weight * trig_ptweight * trig_etaweight;
    }
  }
  if((_AnalyzeTauForLeg1) || (_AnalyzeTauForLeg2)) {
    int theNumberOfTaus = 0;
    double maxptfortrigscalefactor = 0;
    double etafortrigscalefactor = -6;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
      usedTaus.push_back(0);
      theNumberOfTaus++;
      if (passRecoTauCuts((*patTau),theNumberOfTaus-1)) {
        if(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() > maxptfortrigscalefactor) {
          maxptfortrigscalefactor = smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt();
          etafortrigscalefactor = smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).eta();
        }
      }
    }
    if(_ApplyTauTriggerScaleFactors) {
      double trig_ptweight = TauTrigptWeights_[TauTrigptWeights_.size()-1];
      unsigned int nbins = TauTrigPtBinEdges_.size()-1;
      if (maxptfortrigscalefactor>TauTrigPtBinEdges_[0] && maxptfortrigscalefactor<TauTrigPtBinEdges_[nbins]) {
        bool foundCorrectBin = false;
        for (unsigned int jtbin=1; jtbin<=nbins; ++jtbin) {
          if(foundCorrectBin) continue;
          if (maxptfortrigscalefactor>TauTrigPtBinEdges_[jtbin]) continue;
          trig_ptweight = TauTrigptWeights_[jtbin-1];
          foundCorrectBin = true;
        }
      }
      double trig_etaweight = TauTrigetaWeights_[TauTrigetaWeights_.size()-1];
      nbins = TauTrigEtaBinEdges_.size()-1;
      if (etafortrigscalefactor>TauTrigEtaBinEdges_[0] && etafortrigscalefactor<TauTrigEtaBinEdges_[nbins]) {
        bool foundCorrectBin = false;
        for (unsigned int jtbin=1; jtbin<=nbins; ++jtbin) {
          if(foundCorrectBin) continue;
          if (etafortrigscalefactor>TauTrigEtaBinEdges_[jtbin]) continue;
          trig_etaweight = TauTrigetaWeights_[jtbin-1];
          foundCorrectBin = true;
        }
      }  
      trig_weight = trig_weight * trig_ptweight * trig_etaweight;
    }
  }

  //------Number of Good Candidates for leg1
  int nGoodCandidatesLeg1 = 0;
  if(_AnalyzeMuonForLeg1) {
    int theNumberOfMuons = 0;
    for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); 
  	  patMuon != _patMuons->end(); ++patMuon ) {
      theNumberOfMuons++;
      if (!passRecoMuonCuts((*patMuon),theNumberOfMuons-1)) continue;
      if(_CalculateMetUsingOnlyLeg1AndLeg2) {
        temppx += -smearedMuonMomentumVector.at(theNumberOfMuons-1).px();
        temppy += -smearedMuonMomentumVector.at(theNumberOfMuons-1).py();
        usedMuons[theNumberOfMuons-1] = 1;
      }
      nGoodCandidatesLeg1++;
    }
  }
  if(_AnalyzeElectronForLeg1) {
    int theNumberOfElectrons = 0;
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();
          patElectron != _patElectrons->end(); ++patElectron ) {
      theNumberOfElectrons++;
      if (!passRecoElectronCuts((*patElectron),theNumberOfElectrons-1)) continue;
      if(_CalculateMetUsingOnlyLeg1AndLeg2) {
        temppx += -smearedElectronMomentumVector.at(theNumberOfElectrons-1).px();
        temppy += -smearedElectronMomentumVector.at(theNumberOfElectrons-1).py();
        usedElectrons[theNumberOfElectrons-1] = 1;
      }
      nGoodCandidatesLeg1++;
    }
  }
  if(_AnalyzeTauForLeg1) {
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); 
	  patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if (!passRecoTauCuts((*patTau),theNumberOfTaus-1)) continue;
      if(_CalculateMetUsingOnlyLeg1AndLeg2) {
        temppx += -smearedTauMomentumVector.at(theNumberOfTaus-1).px();
        temppy += -smearedTauMomentumVector.at(theNumberOfTaus-1).py();
        usedTaus[theNumberOfTaus-1] = 1;
      }
      nGoodCandidatesLeg1++;
    }
  }
  if (nGoodCandidatesLeg1>=_RecoLeg1Nmin) _EventFlag[_mapSelectionAlgoID["RecoLeg1Nmin"]] = true;
  if (nGoodCandidatesLeg1<=_RecoLeg1Nmax) _EventFlag[_mapSelectionAlgoID["RecoLeg1Nmax"]] = true;

  //------Number of Good Candidates for leg2
  int nGoodCandidatesLeg2 = 0;
  if(_AnalyzeMuonForLeg2) {
    int theNumberOfMuons = 0;
    for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); 
  	  patMuon != _patMuons->end(); ++patMuon ) {
      theNumberOfMuons++;
      if (!passRecoMuonCuts((*patMuon),theNumberOfMuons-1)) continue;
      if((_CalculateMetUsingOnlyLeg1AndLeg2) && (usedMuons[theNumberOfMuons-1] == 0)) {
        temppx += -smearedMuonMomentumVector.at(theNumberOfMuons-1).px();
        temppy += -smearedMuonMomentumVector.at(theNumberOfMuons-1).py();
      }
      nGoodCandidatesLeg2++;
    }
  }
  if(_AnalyzeElectronForLeg2) {
    int theNumberOfElectrons = 0;
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();
          patElectron != _patElectrons->end(); ++patElectron ) {
      theNumberOfElectrons++;
      if (!passRecoElectronCuts((*patElectron),theNumberOfElectrons-1)) continue;
      if((_CalculateMetUsingOnlyLeg1AndLeg2) && (usedElectrons[theNumberOfElectrons-1] == 0)) {
        temppx += -smearedElectronMomentumVector.at(theNumberOfElectrons-1).px();
        temppy += -smearedElectronMomentumVector.at(theNumberOfElectrons-1).py();
      }
      nGoodCandidatesLeg2++;
    }
  }
  if(_AnalyzeTauForLeg2) {
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); 
	  patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if (!passRecoTauCuts((*patTau),theNumberOfTaus-1)) continue;
      if((_CalculateMetUsingOnlyLeg1AndLeg2) && (usedTaus[theNumberOfTaus-1] == 0)) {
        temppx += -smearedTauMomentumVector.at(theNumberOfTaus-1).px();
        temppy += -smearedTauMomentumVector.at(theNumberOfTaus-1).py();
      }
      nGoodCandidatesLeg2++;
    }
  }
  if (nGoodCandidatesLeg2>=_RecoLeg2Nmin) _EventFlag[_mapSelectionAlgoID["RecoLeg2Nmin"]] = true;
  if (nGoodCandidatesLeg2<=_RecoLeg2Nmax) _EventFlag[_mapSelectionAlgoID["RecoLeg2Nmax"]] = true;

  if(_CalculateMetUsingOnlyLeg1AndLeg2) {
    reco::Candidate::LorentzVector theTempMETVector(temppx,temppy,0.0,TMath::Sqrt((temppx*temppx) + (temppy*temppy)));
    theMETVector = theTempMETVector;
  }

  // ------Number of Good Jets   
  sumpxForMht = 0.0;
  sumpyForMht = 0.0;
  sumptForHt  = 0.0;
  int nGoodJets = 0;
  int theNumberOfJets = 0;
  for ( pat::JetCollection::const_iterator patJet = _patJets->begin(); 
	patJet != _patJets->end(); ++patJet ) {
    theNumberOfJets++;
    if (!passRecoJetCuts((*patJet),theNumberOfJets-1)) continue;
    sumpxForMht = sumpxForMht - smearedJetMomentumVector.at(theNumberOfJets-1).px();
    sumpyForMht = sumpyForMht - smearedJetMomentumVector.at(theNumberOfJets-1).py();
    sumptForHt  = sumptForHt  + smearedJetPtEtaPhiMVector.at(theNumberOfJets-1).pt();
    nGoodJets++;
  }
  if (nGoodJets>=_RecoJetNmin) _EventFlag[_mapSelectionAlgoID["RecoJetNmin"]] = true;
  if (nGoodJets<=_RecoJetNmax) _EventFlag[_mapSelectionAlgoID["RecoJetNmax"]] = true;  

  // ------Number of Good Leading Jets
  leadingjetpt = 0;
  theNumberOfJets = 0;
  theLeadingJetIndex = 0;
  for ( pat::JetCollection::const_iterator patJet = _patJets->begin();
        patJet != _patJets->end(); ++patJet ) {
    theNumberOfJets++;
    if(smearedJetPtEtaPhiMVector.at(theNumberOfJets-1).pt() > leadingjetpt) {leadingjetpt = smearedJetPtEtaPhiMVector.at(theNumberOfJets-1).pt(); theLeadingJetIndex = theNumberOfJets;}
  }
  int nGoodFirstLeadingJets = 0;
  theNumberOfJets = 0;
  for ( pat::JetCollection::const_iterator patJet = _patJets->begin();
        patJet != _patJets->end(); ++patJet ) {
    theNumberOfJets++;
    if (theNumberOfJets != theLeadingJetIndex) continue;
    if (!passRecoFirstLeadingJetCuts((*patJet),theNumberOfJets-1)) continue;
    nGoodFirstLeadingJets++;
  }
  if (nGoodFirstLeadingJets>=_RecoFirstLeadingJetNmin) _EventFlag[_mapSelectionAlgoID["RecoFirstLeadingJetNmin"]] = true;

  // ------Number of Good Second Leading Jets
  secondleadingjetpt = 0;
  theNumberOfJets = 0;
  theSecondLeadingJetIndex = 0;
  for ( pat::JetCollection::const_iterator patJet = _patJets->begin();
        patJet != _patJets->end(); ++patJet ) {
    theNumberOfJets++;
    if((smearedJetPtEtaPhiMVector.at(theNumberOfJets-1).pt() > secondleadingjetpt) && (theNumberOfJets != theLeadingJetIndex)) {secondleadingjetpt = smearedJetPtEtaPhiMVector.at(theNumberOfJets-1).pt(); theSecondLeadingJetIndex = theNumberOfJets;}
  }
  int nGoodSecondLeadingJets = 0;
  theNumberOfJets = 0;
  for ( pat::JetCollection::const_iterator patJet = _patJets->begin();
        patJet != _patJets->end(); ++patJet ) {
    theNumberOfJets++;
    if (theNumberOfJets != theSecondLeadingJetIndex) continue;
    if (!passRecoSecondLeadingJetCuts((*patJet),theNumberOfJets-1)) continue;
    nGoodSecondLeadingJets++;
  }
  if (nGoodSecondLeadingJets>=_RecoSecondLeadingJetNmin) _EventFlag[_mapSelectionAlgoID["RecoSecondLeadingJetNmin"]] = true;

  // ------Number of Good Susy Combinations (jet1+jet2+met combinations)
  int nGoodSusyCombinations = 0;
  int numberJets1 = 0;
  for ( pat::JetCollection::const_iterator patJet1 = _patJets->begin();patJet1 != _patJets->end(); ++patJet1 ) {
    numberJets1++;
    if( (numberJets1 == theLeadingJetIndex) && (passRecoFirstLeadingJetCuts((*patJet1),numberJets1 - 1)) ) {
      int numberJets2 = 0;
      for ( pat::JetCollection::const_iterator patJet2 = _patJets->begin();patJet2 != _patJets->end(); ++patJet2 ) {
        numberJets2++;
        if( (numberJets2 == theSecondLeadingJetIndex) && (passRecoSecondLeadingJetCuts((*patJet2),numberJets2 - 1)) ) {
          if(passSusyTopologyCuts(numberJets1 - 1,numberJets2 - 1)) {nGoodSusyCombinations++;}
        }
      }
    }
  }
  if (nGoodSusyCombinations>=_SusyCombinationsNmin) _EventFlag[_mapSelectionAlgoID["SusyCombinationsNmin"]] = true;

  // ------Number of Good Combinations (leg1+leg2+met combinations)
  int nGoodCombinations = 0;
  if( ((_AnalyzeMuonForLeg1) && (_AnalyzeTauForLeg2)) || ((_AnalyzeMuonForLeg2) && (_AnalyzeTauForLeg1)) ) {
    int theNumberOfMuons = 0;
    for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
      theNumberOfMuons++;
      int theNumberOfTaus = 0;
      for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
        theNumberOfTaus++;
        if ((passRecoTauCuts((*patTau),theNumberOfTaus-1)) &&
            (passRecoMuonCuts((*patMuon),theNumberOfMuons-1)) && 
            (passTopologyCuts((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1))) {
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
        if ((passRecoTauCuts((*patTau),theNumberOfTaus-1)) &&
            (passRecoElectronCuts((*patElectron),theNumberOfElectrons-1)) &&
            (passTopologyCuts((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1))) {
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
        if ((passRecoElectronCuts((*patElectron),theNumberOfElectrons-1)) && 
            (passRecoMuonCuts((*patMuon),theNumberOfMuons-1)) && 
            (passTopologyCuts((*patElectron),theNumberOfElectrons-1,(*patMuon),theNumberOfMuons-1))) {
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
        if ((passRecoTauCuts((*patTau1),theNumberOfTaus1 - 1)) && 
            (passRecoTauCuts((*patTau2),theNumberOfTaus2 - 1)) && 
            (passTopologyCuts((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1)) &&
            (theNumberOfTaus2 >= theNumberOfTaus1)) {
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
        if ((passRecoMuonCuts((*patMuon1),theNumberOfMuons1 - 1)) && 
            (passRecoMuonCuts((*patMuon2),theNumberOfMuons2 - 1)) && 
            (passTopologyCuts((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1)) &&
            (theNumberOfMuons2 >= theNumberOfMuons1)) {
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
        if ((passRecoElectronCuts((*patElectron1),theNumberOfElectrons1 - 1)) && 
            (passRecoElectronCuts((*patElectron2),theNumberOfElectrons2 - 1)) && 
            (passTopologyCuts((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1)) &&
            (theNumberOfElectrons2 >= theNumberOfElectrons1)) {
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
    if (_EventFlag[i]) { (_mapSelectionCounter[_EventSelectionSequence[i]])++; }
    if (cumulDecision) { (_mapSelectionCounterCumul[_EventSelectionSequence[i]])++; }
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
//      std::cout << "Please use one of the following triggers :" << std::endl;
      for(edm::TriggerNames::Strings::const_iterator triggerName = TheTriggerNames.triggerNames().begin();
          triggerName != TheTriggerNames.triggerNames().end(); ++triggerName ) {
        unsigned int index = TheTriggerNames.triggerIndex(*triggerName);
        if(index < TheTriggerNames.size()) {
//          std::string triggerDecision = (_triggerResults->accept(index)) ? "passed" : "failed";
//          std::cout << " Trigger Name = " << (*triggerName) << " " << triggerDecision << std::endl;
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
bool HiMassTauAnalysis::passRecoTauCuts(const pat::Tau& patTau,int nobj) {
  // ----Matching to gen
  if(_MatchTauToGen) {
    if(_GenParticleSource.label() != "") {
      if(!(matchToGen(patTau).first)) {return false;}
    } else {return false;}
  }
  // ----Acceptance cuts
  if (fabs(smearedTauPtEtaPhiMVector.at(nobj).eta())>_RecoTauEtaCut) {return false;}
  if (smearedTauPtEtaPhiMVector.at(nobj).pt()<_RecoTauPtMinCut) {return false;}
  if (smearedTauPtEtaPhiMVector.at(nobj).pt()>_RecoTauPtMaxCut) {return false;}
  // ----Lead track requirement
  if (_DoRecoTauDiscrByLeadTrack) {
    if (_UseRecoTauDiscrByLeadTrackFlag) { if ( (patTau.tauID(_RecoTauDiscrByLeadTrack.data()) < 0.5) ) {return false;} }
    else {
      if (patTau.isCaloTau()) { if( (!(patTau.leadTrack().isNonnull())) || (patTau.leadTrack()->pt() < _RecoTauLeadTrackThreshold) ) {return false;} }
      else { if( (!(patTau.leadPFChargedHadrCand().isNonnull())) || (patTau.leadPFChargedHadrCand()->pt() < _RecoTauLeadTrackThreshold) ) {return false;} }
    }
  }
  // ----"eVeto" - Lead track minimimum hits requirement && H3x3/LTpt
  if (_DoRecoTauDiscrByLeadTrackNhits) {
    if (patTau.isCaloTau()) {
      if( (!(patTau.leadTrack().isNonnull())) || ((int)(patTau.leadTrack()->recHitsSize()) < _RecoTauLeadTrackMinHits) ) {return false;}
    } else {
      if( (!(patTau.leadPFChargedHadrCand().isNonnull())) || (!(patTau.leadPFChargedHadrCand()->trackRef().isNonnull())) ) {return false;}
      if( (int)(patTau.leadPFChargedHadrCand()->trackRef()->recHitsSize()) < _RecoTauLeadTrackMinHits ) {return false;}
    }
  }
  if (_DoRecoTauDiscrByH3x3OverP) {
    if (patTau.isCaloTau()) {
      if( (!(patTau.leadTrack().isNonnull())) || ((patTau.hcal3x3OverPLead()) <= _RecoTauH3x3OverP) ) {return false;}
    } else {
      if( (!(patTau.leadPFChargedHadrCand().isNonnull())) || (patTau.hcal3x3OverPLead() <= _RecoTauH3x3OverP) ) {return false;}
    }
  }
  // ----Isolation requirement
  if (_DoRecoTauDiscrByIsolation) {
    if (_UseRecoTauDiscrByIsolationFlag) {if ( (patTau.tauID(_RecoTauDiscrByIsolation.data()) < 0.5) ) {return false;}}
    else {
      if(_UseRecoTauIsoSumPtInsteadOfNiso) {
        /*
        if( (CalculateTauTrackIsolation(patTau).second + CalculateTauEcalIsolation(patTau).second) >= _RecoTauIsoSumPtMaxCutValue ) {return false;}
        if( (CalculateTauTrackIsolation(patTau).second + CalculateTauEcalIsolation(patTau).second) < _RecoTauIsoSumPtMinCutValue ) {return false;}
        */
        if( (CalculateTauTrackIsolation(patTau).second) >= _RecoTauTrackIsoSumPtMaxCutValue ) {return false;}
        if( (CalculateTauTrackIsolation(patTau).second) < _RecoTauTrackIsoSumPtMinCutValue ) {return false;}
        if( (CalculateTauEcalIsolation(patTau).second) >= _RecoTauEcalIsoSumPtMaxCutValue ) {return false;}
        if( (CalculateTauEcalIsolation(patTau).second) < _RecoTauEcalIsoSumPtMinCutValue ) {return false;}
      } else {
        //if( (CalculateTauTrackIsolation(patTau).first + CalculateTauEcalIsolation(patTau).first) > _RecoTauNisoMax ) {return false;}
        if( (CalculateTauTrackIsolation(patTau).first) > _RecoTauTrackNisoMax ) {return false;}
        if( (CalculateTauEcalIsolation(patTau).first) > _RecoTauEcalNisoMax ) {return false;}
      }
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
/*
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
*/
  // ----Electron and Muon vetos
  if (_DoRecoTauDiscrAgainstElectron) { if ( (patTau.tauID(_RecoTauDiscrAgainstElectron.data()) < 0.5) ) {return false;} }
  if (_DoRecoTauDiscrByCrackCut) { if(isInTheCracks(smearedTauPtEtaPhiMVector.at(nobj).eta())) {return false;} }
  if (_DoRecoTauDiscrAgainstMuon) { if ( (patTau.tauID(_RecoTauDiscrAgainstMuon.data()) < 0.5) ) {return false;} }
//  if (_DoRecoTauDiscrAgainstMuon) { if ( (patTau.tauID(_RecoTauDiscrAgainstMuon.data()) > 0.5) ) {return false;} }
  return true;
}

// ---------------Apply Muon Cuts
bool HiMassTauAnalysis::passRecoMuonCuts(const pat::Muon& patMuon,int nobj) {
  // ----Matching to gen
  if(_MatchLeptonToGen) {
    if(_GenParticleSource.label() != "") {
      if(!(matchToGen(patMuon).first)) {return false;}
    } else {return false;}
  }
  // ----Require maching tracks in silicon tracker and muon chamber
  if (_DoRecoMuonDiscrByGlobal) {if (!patMuon.isGlobalMuon()) {return false;}}
  // ----Acceptance cuts
  if (fabs(smearedMuonPtEtaPhiMVector.at(nobj).eta())>_RecoMuonEtaCut) {return false;}
  if (smearedMuonPtEtaPhiMVector.at(nobj).pt()<_RecoMuonPtMinCut) {return false;}
  if (smearedMuonPtEtaPhiMVector.at(nobj).pt()>_RecoMuonPtMaxCut) {return false;}
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
    /*
    if( (patMuon.trackIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonTrackIsoTrkThreshold).first +
         patMuon.ecalIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonEcalIsoRecHitThreshold).first) 
         >= _RecoMuonIsoSumPtMaxCutValue) {return false;}
    if( (patMuon.trackIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonTrackIsoTrkThreshold).first +
         patMuon.ecalIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonEcalIsoRecHitThreshold).first) 
         < _RecoMuonIsoSumPtMinCutValue) {return false;}
    */
    if( (patMuon.trackIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonTrackIsoTrkThreshold).first) 
         >= _RecoMuonTrackIsoSumPtMaxCutValue) {return false;}
    if( (patMuon.trackIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonTrackIsoTrkThreshold).first) 
         < _RecoMuonTrackIsoSumPtMinCutValue) {return false;}
    if( (patMuon.ecalIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonEcalIsoRecHitThreshold).first) 
         >= _RecoMuonEcalIsoSumPtMaxCutValue) {return false;}
    if( (patMuon.ecalIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonEcalIsoRecHitThreshold).first) 
         < _RecoMuonEcalIsoSumPtMinCutValue) {return false;}
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
  if(_RecoTriggersNmin > 0) {
    bool matchedToTrigger = false;
    const trigger::TriggerObjectCollection & toc(handleTriggerEvent->getObjects());
//    std::vector<reco::Particle>  HLTMuMatched; 
//    int nMuHLT=0;
    for ( size_t ia = 0; ia < handleTriggerEvent->sizeFilters(); ++ ia) {
      std::string fullname = handleTriggerEvent->filterTag(ia).encode();
      std::string name;
      size_t p = fullname.find_first_of(':');
      if ( p != std::string::npos) {name = fullname.substr(0, p);}
      else {name = fullname;}
      if ( &toc !=0 ) {
        const trigger::Keys & k = handleTriggerEvent->filterKeys(ia);
        for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {
          if (name == "hltSingleMu9L3Filtered9"  ) { 
//            HLTMuMatched.push_back(toc[*ki].particle());
//            nMuHLT++;
//            std::cout << "name = " << name << std::endl;
//            std::cout << "pt = " << toc[*ki].particle().pt() << std::endl;
//            std::cout << "eta = " << toc[*ki].particle().eta() << std::endl;
            if( (reco::deltaR(smearedMuonPtEtaPhiMVector.at(nobj).eta(),smearedMuonPtEtaPhiMVector.at(nobj).phi(),toc[*ki].particle().eta(),toc[*ki].particle().phi()) <= 0.5) && 
                (toc[*ki].particle().pt() > 15.0) ) {matchedToTrigger = true;}
          }
        }    
      }
    }
//    std::cout << "number of trigs = " << nMuHLT << std::endl;
    if(!(matchedToTrigger)) {return false;}
  }
  return true;
}

//--------------Apply Electron Cuts
bool HiMassTauAnalysis::passRecoElectronCuts(const pat::Electron& patElectron,int nobj) {
  heep::Ele theHeepElec(patElectron);
  int cutResult = patElectron.userInt("HEEPId");
  const Track* elTrack = (const reco::Track*)(patElectron.gsfTrack().get());
  const HitPattern& pInner = elTrack->trackerExpectedHitsInner();
  // ----Matching to gen
  if(_MatchLeptonToGen) {
    if(_GenParticleSource.label() != "") {
      if(!(matchToGen(patElectron).first)) {return false;}
    } else {return false;}
  }
  // ----Acceptance cuts
  if (fabs(smearedElectronPtEtaPhiMVector.at(nobj).eta())>_RecoElectronEtaCut) {return false;}
  if (smearedElectronPtEtaPhiMVector.at(nobj).pt()<_RecoElectronPtMinCut) {return false;}
  if (smearedElectronPtEtaPhiMVector.at(nobj).pt()>_RecoElectronPtMaxCut) {return false;}
  // ----Isolation requirement
  if (_DoRecoElectronDiscrByTrackIsolation) {
//    if ( (patElectron.trackIsoDeposit()->depositAndCountWithin(_RecoElectronTrackIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoElectronTrackIsoTrkThreshold).first >= _RecoElectronTrackIsoSumPtCutValue) ) {return false;}
//    if(patElectron.trackIso() > _RecoElectronTrackIsoSumPtCutValue) {return false;}
    if(patElectron.dr04TkSumPt() > _RecoElectronTrackIsoSumPtMaxCutValue) {return false;}
    if(patElectron.dr04TkSumPt() < _RecoElectronTrackIsoSumPtMinCutValue) {return false;}
  }
  if (_DoRecoElectronDiscrByEcalIsolation) {
//    if ( (patElectron.ecalIsoDeposit()->depositAndCountWithin(_RecoElectronEcalIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoElectronEcalIsoRecHitThreshold).first >= _RecoElectronEcalIsoSumPtCutValue) ) {return false;}
//    if(patElectron.ecalIso() > _RecoElectronEcalIsoSumPtCutValue) {return false;}
    if(patElectron.dr04EcalRecHitSumEt() > _RecoElectronEcalIsoSumPtMaxCutValue) {return false;}
    if(patElectron.dr04EcalRecHitSumEt() < _RecoElectronEcalIsoSumPtMinCutValue) {return false;}
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
  if (_DoRecoElectronDiscrByEcalDrivenSeed) { if(!(patElectron.ecalDrivenSeed())) {return false;} }
  if (_DoRecoElectronDiscrByTrackerDrivenSeed) { if(!(patElectron.trackerDrivenSeed())) {return false;} }
  if (_DoRecoElectronDiscrByHoverEm) {
    if(_UseHeepInfo) {
      if(!(heep::CutCodes::passCuts(cutResult, "hadem"))) {return false;}
    } else { if(patElectron.hadronicOverEm() > _RecoElectronHoverEmCut) {return false;} }
  }
  if (_DoRecoElectronDiscrBySigmaIEtaIEta) {
    if(_UseHeepInfo) {
      if(!(heep::CutCodes::passCuts(cutResult, "sigmaIEtaIEta"))) {return false;}
    } else {
      if(patElectron.isEE()) { if(patElectron.scSigmaIEtaIEta() > _RecoElectronSigmaIEtaIEta) {return false;} }
    }
  }
  if (_DoRecoElectronDiscrByDEtaIn) {
    if(_UseHeepInfo) {
      if(!(heep::CutCodes::passCuts(cutResult, "dEtaIn"))) {return false;}
    } else {
      if(patElectron.isEE()) { if(fabs(patElectron.deltaEtaSuperClusterTrackAtVtx()) > _RecoElectronEEDEtaIn) {return false;} }
      if(patElectron.isEB()) { if(fabs(patElectron.deltaEtaSuperClusterTrackAtVtx()) > _RecoElectronEBDEtaIn) {return false;} }
    }
  }
  if (_DoRecoElectronDiscrByDPhiIn) {
    if(_UseHeepInfo) {
      if(!(heep::CutCodes::passCuts(cutResult, "dPhiIn"))) {return false;}
    } else {
      if(patElectron.isEE()) { if(fabs(patElectron.deltaPhiSuperClusterTrackAtVtx()) > _RecoElectronEEDPhiIn) {return false;} }
      if(patElectron.isEB()) { if(fabs(patElectron.deltaPhiSuperClusterTrackAtVtx()) > _RecoElectronEBDPhiIn) {return false;} }
    }
  }
  if (_DoRecoElectronDiscrBySCE2by5Over5by5) {
    if(_UseHeepInfo) {
      if(!(heep::CutCodes::passCuts(cutResult, "e2x5Over5x5"))) {return false;}
    } else {
      if(patElectron.isEB()) {
        if( ((patElectron.scE2x5Max() / patElectron.scE5x5()) > _RecoElectronEBscE2by5Over5by5) || 
            ((patElectron.scE1x5() / patElectron.scE5x5()) > _RecoElectronEBscE1by5Over5by5) ) { }
        else {return false;}
      }
    }
  }
  if (_DoRecoElectronDiscrByMissingHits) { if(pInner.numberOfHits() >= _RecoElectronMissingHits) {return false;} }
  return true;
}

//--------------Apply Jet Cuts
bool HiMassTauAnalysis::passRecoJetCuts(const pat::Jet& patJet,int nobj) {
  // ----Acceptance cuts
  if (fabs(smearedJetPtEtaPhiMVector.at(nobj).eta())>_RecoJetEtaMaxCut) {return false;}
  if (fabs(smearedJetPtEtaPhiMVector.at(nobj).eta())<_RecoJetEtaMinCut) {return false;}
  if (smearedJetPtEtaPhiMVector.at(nobj).pt()<_RecoJetPtCut) {return false;}
  // ----anti-overlap requirements
  if (_RemoveJetOverlapWithMuons) {
    int theNumberOfMuons = 0;
    for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
      theNumberOfMuons++;
      if( (passRecoMuonCuts((*patMuon),theNumberOfMuons-1)) ) {
        if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedMuonMomentumVector.at(theNumberOfMuons-1)) < _JetMuonMatchingDeltaR) {return false;}
      }
    }
  }
  if (_RemoveJetOverlapWithElectrons) {
    int theNumberOfElectrons = 0;
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron ) {
      theNumberOfElectrons++;
      if( (passRecoElectronCuts((*patElectron),theNumberOfElectrons-1)) ) {
        if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedElectronMomentumVector.at(theNumberOfElectrons-1)) < _JetElectronMatchingDeltaR) {return false;}
      }
    }
  }
  if (_RemoveJetOverlapWithTaus) {
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if( (passRecoTauCuts((*patTau),theNumberOfTaus-1)) ) {
        if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedTauMomentumVector.at(theNumberOfTaus-1)) < _JetTauMatchingDeltaR) {return false;}
      }
    }
  }
  if (_ApplyJetBTagging) {if(patJet.bDiscriminator("trackCountingHighEffBJetTags") <= _JetBTaggingTCHEcut) {return false;}}
  return true;
}

//--------------Apply Leading Jet Cuts
bool HiMassTauAnalysis::passRecoFirstLeadingJetCuts(const pat::Jet& patJet,int nobj) {
  if(_DoDiscrByFirstLeadingJet) {
    // ----Acceptance cuts
    if (fabs(smearedJetPtEtaPhiMVector.at(nobj).eta())>_RecoFirstLeadingJetEtaMaxCut) {return false;}
    if (fabs(smearedJetPtEtaPhiMVector.at(nobj).eta())<_RecoFirstLeadingJetEtaMinCut) {return false;}
    if (smearedJetPtEtaPhiMVector.at(nobj).pt()<_RecoFirstLeadingJetPt) {return false;}
  }
  return true;
}

//--------------Apply Second Leading Jet Cuts
bool HiMassTauAnalysis::passRecoSecondLeadingJetCuts(const pat::Jet& patJet,int nobj) {
  if(_DoDiscrBySecondLeadingJet) {
    // ----Acceptance cuts
    if (fabs(smearedJetPtEtaPhiMVector.at(nobj).eta())>_RecoSecondLeadingJetEtaMaxCut) {return false;}
    if (fabs(smearedJetPtEtaPhiMVector.at(nobj).eta())<_RecoSecondLeadingJetEtaMinCut) {return false;}
    if (smearedJetPtEtaPhiMVector.at(nobj).pt()<_RecoSecondLeadingJetPt) {return false;}
  }
  return true;
}

// ---------------Apply Topology Cuts
bool HiMassTauAnalysis::passTopologyCuts(const pat::Tau& patTau, int nobjT, const pat::Muon& patMuon, int nobjM) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(smearedTauMomentumVector.at(nobjT), smearedMuonMomentumVector.at(nobjM)) < _DiTauDeltaRCut) {return false;} }
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
  if (_DoDiTauDiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) > _DiTauCosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) < _DiTauCosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByMassReco) { if( (CalculateThe4Momentum(patTau,nobjT,patMuon,nobjM).second.M() < _MassMinCut) || (CalculateThe4Momentum(patTau,nobjT,patMuon,nobjM).second.M() > _MassMaxCut) ) {return false;} }
  // ----Zeta requirement
  if (_DoDiTauDiscrByCDFzeta2D) {
    if( ((_PZetaCutCoefficient * CalculatePZeta(patTau,nobjT,patMuon,nobjM)) + 
         (_PZetaVisCutCoefficient * CalculatePZetaVis(patTau,nobjT,patMuon,nobjM))) < _CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoDiTauDiscrByDeltaPtDivSumPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedMuonPtEtaPhiMVector.at(nobjM).pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedMuonPtEtaPhiMVector.at(nobjM).pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoDiTauDiscrByDeltaPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt())) < _DiTauDeltaPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt())) > _DiTauDeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByLeg1MetDphi) {
    if(_AnalyzeMuonForLeg1) {
      if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) > _Leg1MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) < _Leg1MetDphiMinCut) {return false;}
    }
    if(_AnalyzeTauForLeg1) {
      if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) > _Leg1MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) < _Leg1MetDphiMinCut) {return false;}
    }
  }
  if (_DoDiscrByLeg2MetDphi) {
    if(_AnalyzeMuonForLeg2) {
      if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) > _Leg2MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) < _Leg2MetDphiMinCut) {return false;}
    }
    if(_AnalyzeTauForLeg2) {
      if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) > _Leg2MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) < _Leg2MetDphiMinCut) {return false;}
    }
  }
  return true;
}
bool HiMassTauAnalysis::passTopologyCuts(const pat::Tau& patTau, int nobjT, const pat::Electron& patElectron, int nobjE) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(smearedTauMomentumVector.at(nobjT), smearedElectronMomentumVector.at(nobjE)) < _DiTauDeltaRCut) {return false;} }
  // ---- apply Zee veto cut
  if (_DoTauDiscrByIsZeeCut) { if((isZee(smearedElectronMomentumVector.at(nobjE)).first)) {return false;} }
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
  if (_DoDiTauDiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedElectronPtEtaPhiMVector.at(nobjE).phi()))) > _DiTauCosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedElectronPtEtaPhiMVector.at(nobjE).phi()))) < _DiTauCosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByMassReco) { if( (CalculateThe4Momentum(patTau,nobjT,patElectron,nobjE).second.M() < _MassMinCut) || (CalculateThe4Momentum(patTau,nobjT,patElectron,nobjE).second.M() > _MassMaxCut) ) {return false;} }
  // ----Zeta requirement
  if (_DoDiTauDiscrByCDFzeta2D) {
    if( ((_PZetaCutCoefficient * CalculatePZeta(patTau,nobjT,patElectron,nobjE)) + 
         (_PZetaVisCutCoefficient * CalculatePZetaVis(patTau,nobjT,patElectron,nobjE))) < _CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoDiTauDiscrByDeltaPtDivSumPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoDiTauDiscrByDeltaPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _DiTauDeltaPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _DiTauDeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByLeg1MetDphi) {
    if(_AnalyzeElectronForLeg1) {
      if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) > _Leg1MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) < _Leg1MetDphiMinCut) {return false;}
    }
    if(_AnalyzeTauForLeg1) {
      if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) > _Leg1MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) < _Leg1MetDphiMinCut) {return false;}
    }
  }
  if (_DoDiscrByLeg2MetDphi) {
    if(_AnalyzeElectronForLeg2) {
      if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) > _Leg2MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) < _Leg2MetDphiMinCut) {return false;}
    }
    if(_AnalyzeTauForLeg2) {
      if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) > _Leg2MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) < _Leg2MetDphiMinCut) {return false;}
    }
  }
  return true;
}
bool HiMassTauAnalysis::passTopologyCuts(const pat::Electron& patElectron, int nobjE, const pat::Muon& patMuon, int nobjM) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(smearedMuonMomentumVector.at(nobjM), smearedElectronMomentumVector.at(nobjE)) < _DiTauDeltaRCut) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_DiTauDiscrByOSLSType == "OS") {
    if((patElectron.charge() * patMuon.charge()) >= 0) {return false;}
  } else if (_DiTauDiscrByOSLSType == "LS") {
    if((patElectron.charge() * patMuon.charge()) <= 0) {return false;}
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoDiTauDiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) > _DiTauCosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) < _DiTauCosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByMassReco) { if( (CalculateThe4Momentum(patElectron,nobjE,patMuon,nobjM).second.M() < _MassMinCut) || (CalculateThe4Momentum(patElectron,nobjE,patMuon,nobjM).second.M() > _MassMaxCut) ) {return false;} }
  // ----Zeta requirement
  if (_DoDiTauDiscrByCDFzeta2D) {
    if( ((_PZetaCutCoefficient * CalculatePZeta(patElectron,nobjE,patMuon,nobjM)) +
         (_PZetaVisCutCoefficient * CalculatePZetaVis(patElectron,nobjE,patMuon,nobjM))) < _CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoDiTauDiscrByDeltaPtDivSumPt) {
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedMuonPtEtaPhiMVector.at(nobjM).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedMuonPtEtaPhiMVector.at(nobjM).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoDiTauDiscrByDeltaPt) {
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _DiTauDeltaPtMinCutValue ) {return false;}
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _DiTauDeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByLeg1MetDphi) {
    if(_AnalyzeElectronForLeg1) {
      if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) > _Leg1MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) < _Leg1MetDphiMinCut) {return false;}
    }
    if(_AnalyzeMuonForLeg1) {
      if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) > _Leg1MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) < _Leg1MetDphiMinCut) {return false;}
    }
  }
  if (_DoDiscrByLeg2MetDphi) {
    if(_AnalyzeElectronForLeg2) {
      if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) > _Leg2MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) < _Leg2MetDphiMinCut) {return false;}
    }
    if(_AnalyzeMuonForLeg2) {
      if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) > _Leg2MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) < _Leg2MetDphiMinCut) {return false;}
    }
  }
  return true;
}
bool HiMassTauAnalysis::passTopologyCuts(const pat::Muon& patMuon1, int nobj1, const pat::Muon& patMuon2, int nobj2) {
  // ----Separation cut between muon legs (remove double counting)
  if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(smearedMuonMomentumVector.at(nobj1), smearedMuonMomentumVector.at(nobj2)) < _DiTauDeltaRCut) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_DiTauDiscrByOSLSType == "OS") {
    if((patMuon1.charge() * patMuon2.charge()) >= 0) {return false;}
  } else if (_DiTauDiscrByOSLSType == "LS") {
    if((patMuon1.charge() * patMuon2.charge()) <= 0) {return false;}
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoDiTauDiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobj1).phi() - smearedMuonPtEtaPhiMVector.at(nobj2).phi()))) > _DiTauCosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobj1).phi() - smearedMuonPtEtaPhiMVector.at(nobj2).phi()))) < _DiTauCosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByMassReco) { if( (CalculateThe4Momentum(patMuon1,nobj1,patMuon2,nobj2).second.M() < _MassMinCut) || (CalculateThe4Momentum(patMuon1,nobj1,patMuon2,nobj2).second.M() > _MassMaxCut) ) {return false;} }
  // ----Zeta requirement
  if (_DoDiTauDiscrByCDFzeta2D) {
    if( ((_PZetaCutCoefficient * CalculatePZeta(patMuon1,nobj1,patMuon2,nobj2)) +
         (_PZetaVisCutCoefficient * CalculatePZetaVis(patMuon1,nobj1,patMuon2,nobj2))) < _CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoDiTauDiscrByDeltaPtDivSumPt) {
    if ( ((smearedMuonPtEtaPhiMVector.at(nobj1).pt() - smearedMuonPtEtaPhiMVector.at(nobj2).pt()) / (smearedMuonPtEtaPhiMVector.at(nobj1).pt() + smearedMuonPtEtaPhiMVector.at(nobj2).pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedMuonPtEtaPhiMVector.at(nobj1).pt() - smearedMuonPtEtaPhiMVector.at(nobj2).pt()) / (smearedMuonPtEtaPhiMVector.at(nobj1).pt() + smearedMuonPtEtaPhiMVector.at(nobj2).pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoDiTauDiscrByDeltaPt) {
    if ( ((smearedMuonPtEtaPhiMVector.at(nobj1).pt() - smearedMuonPtEtaPhiMVector.at(nobj2).pt())) < _DiTauDeltaPtMinCutValue ) {return false;}
    if ( ((smearedMuonPtEtaPhiMVector.at(nobj1).pt() - smearedMuonPtEtaPhiMVector.at(nobj2).pt())) > _DiTauDeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByLeg1MetDphi) {
    if(_AnalyzeMuonForLeg1) {
      if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobj1).phi() - theMETVector.phi())) > _Leg1MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobj1).phi() - theMETVector.phi())) < _Leg1MetDphiMinCut) {return false;}
    }
  }
  if (_DoDiscrByLeg2MetDphi) {
    if(_AnalyzeMuonForLeg2) {
      if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobj2).phi() - theMETVector.phi())) > _Leg2MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobj2).phi() - theMETVector.phi())) < _Leg2MetDphiMinCut) {return false;}
    }
  }
  return true;
}
bool HiMassTauAnalysis::passTopologyCuts(const pat::Electron& patElectron1, int nobj1, const pat::Electron& patElectron2, int nobj2) {
  // ----Separation cut between electron legs (remove double counting)
  if (_DoDiTauDiscrByDeltaR) { if(reco::deltaR(smearedElectronMomentumVector.at(nobj1), smearedElectronMomentumVector.at(nobj2)) < _DiTauDeltaRCut) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_DiTauDiscrByOSLSType == "OS") {
    if((patElectron1.charge() * patElectron2.charge()) >= 0) {return false;}
  } else if (_DiTauDiscrByOSLSType == "LS") {
    if((patElectron1.charge() * patElectron2.charge()) <= 0) {return false;}
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoDiTauDiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobj1).phi() - smearedElectronPtEtaPhiMVector.at(nobj2).phi()))) > _DiTauCosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobj1).phi() - smearedElectronPtEtaPhiMVector.at(nobj2).phi()))) < _DiTauCosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByMassReco) { if( (CalculateThe4Momentum(patElectron1,nobj1,patElectron2,nobj2).second.M() < _MassMinCut) || (CalculateThe4Momentum(patElectron1,nobj1,patElectron2,nobj2).second.M() > _MassMaxCut) ) {return false;} }
  // ----Zeta requirement
  if (_DoDiTauDiscrByCDFzeta2D) {
    if( ((_PZetaCutCoefficient * CalculatePZeta(patElectron1,nobj1,patElectron2,nobj2)) +
         (_PZetaVisCutCoefficient * CalculatePZetaVis(patElectron1,nobj1,patElectron2,nobj2))) < _CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoDiTauDiscrByDeltaPtDivSumPt) {
    if ( ((smearedElectronPtEtaPhiMVector.at(nobj1).pt() - smearedElectronPtEtaPhiMVector.at(nobj2).pt()) / (smearedElectronPtEtaPhiMVector.at(nobj1).pt() + smearedElectronPtEtaPhiMVector.at(nobj2).pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedElectronPtEtaPhiMVector.at(nobj1).pt() - smearedElectronPtEtaPhiMVector.at(nobj2).pt()) / (smearedElectronPtEtaPhiMVector.at(nobj1).pt() + smearedElectronPtEtaPhiMVector.at(nobj2).pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoDiTauDiscrByDeltaPt) {
    if ( ((smearedElectronPtEtaPhiMVector.at(nobj1).pt() - smearedElectronPtEtaPhiMVector.at(nobj2).pt())) < _DiTauDeltaPtMinCutValue ) {return false;}
    if ( ((smearedElectronPtEtaPhiMVector.at(nobj1).pt() - smearedElectronPtEtaPhiMVector.at(nobj2).pt())) > _DiTauDeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByLeg1MetDphi) {
    if(_AnalyzeElectronForLeg1) {
      if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobj1).phi() - theMETVector.phi())) > _Leg1MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobj1).phi() - theMETVector.phi())) < _Leg1MetDphiMinCut) {return false;}
    }
  }
  if (_DoDiscrByLeg2MetDphi) {
    if(_AnalyzeElectronForLeg2) {
      if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobj2).phi() - theMETVector.phi())) > _Leg2MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobj2).phi() - theMETVector.phi())) < _Leg2MetDphiMinCut) {return false;}
    }
  }
  return true;
}
bool HiMassTauAnalysis::passTopologyCuts(const pat::Tau& patTau1, int nobj1, const pat::Tau& patTau2, int nobj2) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoDiTauDiscrByDeltaR) {if(reco::deltaR(smearedTauMomentumVector.at(nobj1), smearedTauMomentumVector.at(nobj2)) < _DiTauDeltaRCut) {return false;}}
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
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobj1).phi() - smearedTauPtEtaPhiMVector.at(nobj2).phi()))) > _DiTauCosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobj1).phi() - smearedTauPtEtaPhiMVector.at(nobj2).phi()))) < _DiTauCosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByMassReco) { if( (CalculateThe4Momentum(patTau1,nobj1,patTau2,nobj2).second.M() < _MassMinCut) || (CalculateThe4Momentum(patTau1,nobj1,patTau2,nobj2).second.M() > _MassMaxCut) ) {return false;} }
  // ----Zeta requirement
  if (_DoDiTauDiscrByCDFzeta2D) {
    if( ((_PZetaCutCoefficient * CalculatePZeta(patTau1,nobj1,patTau2,nobj2)) +
         (_PZetaVisCutCoefficient * CalculatePZetaVis(patTau1,nobj1,patTau2,nobj2))) < _CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoDiTauDiscrByDeltaPtDivSumPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobj1).pt() - smearedTauPtEtaPhiMVector.at(nobj2).pt()) / (smearedTauPtEtaPhiMVector.at(nobj1).pt() + smearedTauPtEtaPhiMVector.at(nobj2).pt())) < _DiTauDeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobj1).pt() - smearedTauPtEtaPhiMVector.at(nobj2).pt()) / (smearedTauPtEtaPhiMVector.at(nobj1).pt() + smearedTauPtEtaPhiMVector.at(nobj2).pt())) > _DiTauDeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoDiTauDiscrByDeltaPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobj1).pt() - smearedTauPtEtaPhiMVector.at(nobj2).pt())) < _DiTauDeltaPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobj1).pt() - smearedTauPtEtaPhiMVector.at(nobj2).pt())) > _DiTauDeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByLeg1MetDphi) {
    if(_AnalyzeTauForLeg1) {
      if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobj1).phi() - theMETVector.phi())) > _Leg1MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobj1).phi() - theMETVector.phi())) < _Leg1MetDphiMinCut) {return false;}
    }
  }
  if (_DoDiscrByLeg2MetDphi) {
    if(_AnalyzeTauForLeg2) {
      if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobj2).phi() - theMETVector.phi())) > _Leg2MetDphiMaxCut) {return false;}
      if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobj2).phi() - theMETVector.phi())) < _Leg2MetDphiMinCut) {return false;}
    }
  }
  return true;
}
bool HiMassTauAnalysis::passSusyTopologyCuts(int nobj1, int nobj2) {
  if(_DoSUSYDiscrByMHT) { if(sqrt((sumpxForMht * sumpxForMht) + (sumpyForMht * sumpyForMht)) < _MhtCut) {return false;} }
  double dphi1;
  double dphi2;
  double r1;
  double r2;
  double alpha;
  if(_DoSUSYDiscrByR1) {
    dphi1 = TMath::Abs(normalizedPhi(smearedJetPtEtaPhiMVector.at(nobj1).phi() - theMETVector.phi()));
    dphi2 = TMath::Abs(normalizedPhi(smearedJetPtEtaPhiMVector.at(nobj2).phi() - theMETVector.phi()));
    r1 = sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0) );
    if(r1 < _R1MinCut) {return false;}
    if(r1 > _R1MaxCut) {return false;}
  }
  if(_DoSUSYDiscrByR2) {
    dphi1 = TMath::Abs(normalizedPhi(smearedJetPtEtaPhiMVector.at(nobj1).phi() - theMETVector.phi()));
    dphi2 = TMath::Abs(normalizedPhi(smearedJetPtEtaPhiMVector.at(nobj2).phi() - theMETVector.phi()));
    r2 = sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0) );
    if(r2 < _R2MinCut) {return false;}
    if(r2 > _R2MaxCut) {return false;}
  }
  if(_DoSUSYDiscrByAlpha) {
    double px = smearedJetMomentumVector.at(nobj1).px() + smearedJetMomentumVector.at(nobj2).px();
    double py = smearedJetMomentumVector.at(nobj1).py() + smearedJetMomentumVector.at(nobj2).py();
    double pz = smearedJetMomentumVector.at(nobj1).pz() + smearedJetMomentumVector.at(nobj2).pz();
    double e = smearedJetMomentumVector.at(nobj1).energy() + smearedJetMomentumVector.at(nobj2).energy();
    reco::Candidate::LorentzVector Susy_LorentzVect(px, py, pz, e);
    if(Susy_LorentzVect.M() > 0) {alpha = smearedJetPtEtaPhiMVector.at(nobj2).pt() / Susy_LorentzVect.M();}
    else {alpha = 999999999.0;}
    if(alpha < _AlphaMinCut) {return false;}
    if(alpha > _AlphaMaxCut) {return false;}
  }
  if(_DoSUSYDiscrByDphi1) {
    dphi1 = TMath::Abs(normalizedPhi(smearedJetPtEtaPhiMVector.at(nobj1).phi() - theMETVector.phi()));
    if(dphi1 < _Dphi1MinCut) {return false;}
    if(dphi1 > _Dphi1MaxCut) {return false;}
  }
  if(_DoSUSYDiscrByDphi2) {
    dphi2 = TMath::Abs(normalizedPhi(smearedJetPtEtaPhiMVector.at(nobj2).phi() - theMETVector.phi()));
    if(dphi2 < _Dphi2MinCut) {return false;}
    if(dphi2 > _Dphi2MaxCut) {return false;}
  }
  return true;
}

// ---------------Fill Ntuple
void HiMassTauAnalysis::fillNtuple() {

  // fill jet information used for jet based cuts (e.g. jet veto)
  int theNumberOfJets = 0;
  for ( pat::JetCollection::const_iterator patJet = _patJets->begin(); patJet != _patJets->end(); ++patJet ) {
    _bJetDiscrByTrackCounting->push_back(patJet->bDiscriminator("trackCountingHighEffBJetTags"));	     
    _bJetDiscrBySimpleSecondaryV->push_back(patJet->bDiscriminator("simpleSecondaryVertexBJetTags"));	   
    _bJetDiscrByCombinedSecondaryV->push_back(patJet->bDiscriminator("combinedSecondaryVertexBJetTags")); 
    _jetPt->push_back(smearedJetPtEtaPhiMVector.at(theNumberOfJets).pt());
    _jetEt->push_back(smearedJetMomentumVector.at(theNumberOfJets).energy() * sin(smearedJetMomentumVector.at(theNumberOfJets).theta()));
    _jetEta->push_back(smearedJetPtEtaPhiMVector.at(theNumberOfJets).eta()); 
    _jetPhi->push_back(smearedJetPtEtaPhiMVector.at(theNumberOfJets).phi());
    _jetEmFraction->push_back(patJet->correctedJet("raw", "").emEnergyFraction());
    theNumberOfJets++;
  }

  int theNumberOfMuons = 0;
  for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
    theNumberOfMuons++;
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if ((passRecoTauCuts((*patTau),theNumberOfTaus-1)) &&
	  (passRecoMuonCuts((*patMuon),theNumberOfMuons-1)) &&
	  (passTopologyCuts((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1))) {
        _mEt->push_back(theMETVector.pt());
        _muTauMetMass->push_back( CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M() );      
        _muMetMass->push_back(CalculateLeptonMetMt((*patMuon),theNumberOfMuons-1));
        _zeta->push_back( (_PZetaCutCoefficient * CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1)) +
                          (_PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1)) );  
        _muTrkIsoPat->push_back(patMuon->trackIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonTrackIsoTrkThreshold).first);
        if (patTau->isCaloTau()) {
          if (_UseRecoTauDiscrByIsolationFlag) {_tauTrkIsoPat->push_back(patTau->isolationTracksPtSum());}
          else {_tauTrkIsoPat->push_back(CalculateTauTrackIsolation(*patTau).second);}
        } else {
          if (_UseRecoTauDiscrByIsolationFlag) {_tauTrkIsoPat->push_back(patTau->isolationPFChargedHadrCandsPtSum());}
          else {_tauTrkIsoPat->push_back(CalculateTauTrackIsolation(*patTau).second);}
        }
	if(_UseTauSeedTrackForDiTauDiscrByOSLS) {
	  if (patTau->isCaloTau()) {
	    if( (patTau->leadTrack().isNonnull()) ) {_OSLS->push_back(patMuon->charge() * patTau->leadTrack()->charge());}
	  } else {
            if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {_OSLS->push_back(patMuon->charge() * patTau->leadPFChargedHadrCand()->charge());}
	  }
	} else {_OSLS->push_back(patMuon->charge() * patTau->charge());}
      }
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
      for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
	theNumberOfTaus++;
	if (!passRecoTauCuts((*patTau),theNumberOfTaus-1)) continue;
	_hTauJetEnergy[NpdfID]->Fill(smearedTauMomentumVector.at(theNumberOfTaus-1).energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJetPt[NpdfID]->Fill(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJetEta[NpdfID]->Fill(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJetPhi[NpdfID]->Fill(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	if (patTau->isCaloTau()) {
	  _hTauJetNumSignalTracks[NpdfID]->Fill(patTau->signalTracks().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJetNumSignalGammas[NpdfID]->Fill(0.0,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJetCharge[NpdfID]->Fill(fabs(patTau->charge()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJetSignalTracksMass[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJetSignalTracksAndGammasMass[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJetSignalTracksChargeFraction[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).pt() / smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  if(patTau->signalTracks().size() == 1) {_hTauJetSignalTracksMass1prong[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	  if(patTau->signalTracks().size() == 3) {_hTauJetSignalTracksMass3prong[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	  if(patTau->signalTracks().size() == 1) {_hTauJetSignalTracksAndGammasMass1prong[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	  if(patTau->signalTracks().size() == 3) {_hTauJetSignalTracksAndGammasMass3prong[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	  if( (patTau->leadTrack().isNonnull()) ) {
	    _hTauJetSeedTrackPt[NpdfID]->Fill(patTau->leadTrack()->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSeedTrackIpSignificance[NpdfID]->Fill(patTau->leadTracksignedSipt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSeedTrackNhits[NpdfID]->Fill(patTau->leadTrack()->recHitsSize(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSeedTrackChi2[NpdfID]->Fill(patTau->leadTrack()->chi2(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetH3x3OverP[NpdfID]->Fill(patTau->hcal3x3OverPLead(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
          if (_UseRecoTauDiscrByIsolationFlag) {
	    _hTauJetNumIsoTracks[NpdfID]->Fill(patTau->isolationTracks().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetNumIsoGammas[NpdfID]->Fill(0.0,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetNumIsoCands[NpdfID]->Fill(patTau->isolationTracks().size() + 0,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIsoTracks[NpdfID]->Fill(patTau->isolationTracksPtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIsoGammas[NpdfID]->Fill(patTau->isolationECALhitsEtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIso[NpdfID]->Fill(patTau->isolationTracksPtSum() + patTau->isolationECALhitsEtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          } else {
	    _hTauJetNumIsoTracks[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetNumIsoGammas[NpdfID]->Fill(CalculateTauEcalIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetNumIsoCands[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).first + CalculateTauEcalIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIsoTracks[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIsoGammas[NpdfID]->Fill(CalculateTauEcalIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIso[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).second + CalculateTauEcalIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          }
	} else {
	  _hTauJetNumSignalTracks[NpdfID]->Fill(patTau->signalPFChargedHadrCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJetNumSignalGammas[NpdfID]->Fill(CalculateNumberSignalTauGammas(*patTau),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJetCharge[NpdfID]->Fill(fabs(patTau->charge()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJetSignalTracksMass[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJetSignalTracksAndGammasMass[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJetSignalTracksChargeFraction[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).pt() / smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  if(patTau->signalPFChargedHadrCands().size() == 1) {
	    _hTauJetSignalTracksMass1prong[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSignalTracksAndGammasMass1prong[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSignalTracksAndPiZerosMass1prong[NpdfID]->Fill(CalculateTauSignalTracksAndPiZerosMass(*patTau).first.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetNumSignalPiZeros1prong[NpdfID]->Fill(CalculateTauSignalTracksAndPiZerosMass(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
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
	    _hTauJetH3x3OverP[NpdfID]->Fill(patTau->hcal3x3OverPLead(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
          if (_UseRecoTauDiscrByIsolationFlag) {
	    _hTauJetNumIsoTracks[NpdfID]->Fill(patTau->isolationPFChargedHadrCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetNumIsoGammas[NpdfID]->Fill(patTau->isolationPFGammaCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetNumIsoCands[NpdfID]->Fill(patTau->isolationPFChargedHadrCands().size() + patTau->isolationPFGammaCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIsoTracks[NpdfID]->Fill(patTau->isolationPFChargedHadrCandsPtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIsoGammas[NpdfID]->Fill(patTau->isolationPFGammaCandsEtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIso[NpdfID]->Fill(patTau->isolationPFChargedHadrCandsPtSum() + patTau->isolationPFGammaCandsEtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          } else {
	    _hTauJetNumIsoTracks[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetNumIsoGammas[NpdfID]->Fill(CalculateTauEcalIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetNumIsoCands[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).first + CalculateTauEcalIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIsoTracks[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIsoGammas[NpdfID]->Fill(CalculateTauEcalIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetSumPtIso[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).second + CalculateTauEcalIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          }
	}
	if(_GenParticleSource.label() != "") {
	  if(matchToGen(*patTau).first) {
	    _hTauJetGenTauDeltaPhi[NpdfID]->Fill(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi() - matchToGen(*patTau).second.phi()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetGenTauDeltaEta[NpdfID]->Fill(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).eta() - matchToGen(*patTau).second.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJetGenTauDeltaPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - matchToGen(*patTau).second.pt()) / matchToGen(*patTau).second.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
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
      double maxMuonPt = 0;
      for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); 
	    patMuon != _patMuons->end(); ++patMuon ) {
	theNumberOfMuons++;
	if (!passRecoMuonCuts((*patMuon),theNumberOfMuons-1)) continue;
         if(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt() > maxMuonPt){
           maxMuonPt = smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt();
           maxPtMuonVector = smearedMuonMomentumVector.at(theNumberOfMuons-1);
         }
        _hMuonEnergy[NpdfID]->Fill(smearedMuonMomentumVector.at(theNumberOfMuons-1).energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuonPt[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuonEta[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuonPhi[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	if(_GenParticleSource.label() != "") {
	  if(matchToGen(*patMuon).first) {
	    _hMuonGenMuonDeltaPhi[NpdfID]->Fill(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - matchToGen(*patMuon).second.phi()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuonGenMuonDeltaEta[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).eta() - matchToGen(*patMuon).second.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuonGenMuonDeltaPt[NpdfID]->Fill((smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt() - matchToGen(*patMuon).second.pt()) / matchToGen(*patMuon).second.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
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
      double maxElecEt = 0;
      for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();
	    patElectron != _patElectrons->end(); ++patElectron ) {
	theNumberOfElectrons++;
	if (!passRecoElectronCuts((*patElectron),theNumberOfElectrons-1)) continue;
        if(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt() > maxElecEt){
          maxElecEt = smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt();
          maxEtElectronVector = smearedElectronMomentumVector.at(theNumberOfElectrons-1); 
        }
	_hElectronEnergy[NpdfID]->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectronPt[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectronEta[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectronPhi[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	if(_GenParticleSource.label() != "") {
	  if(matchToGen(*patElectron).first) {
	    _hElectronGenElectronDeltaPhi[NpdfID]->Fill(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - matchToGen(*patElectron).second.phi()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectronGenElectronDeltaEta[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).eta() - matchToGen(*patElectron).second.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectronGenElectronDeltaPt[NpdfID]->Fill((smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt() - matchToGen(*patElectron).second.pt()) / matchToGen(*patElectron).second.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
	_hElectronTrackIso[NpdfID]->Fill(patElectron->dr04TkSumPt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectronEcalIso[NpdfID]->Fill(patElectron->dr04EcalRecHitSumEt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	/*
	  _hElectronTrackIso[NpdfID]->Fill(patElectron->trackIsoDeposit()->depositAndCountWithin(_RecoElectronTrackIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoElectronTrackIsoTrkThreshold).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hElectronEcalIso[NpdfID]->Fill(patElectron->ecalIsoDeposit()->depositAndCountWithin(_RecoElectronEcalIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoElectronEcalIsoRecHitThreshold).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	*/
	const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
	if ( patElectron->track().isNonnull() ) { _hElectronIp[NpdfID]->Fill( patElectron->track()->dxy(thePrimaryEventVertex.position()) ,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID)); }
	//if ( patElectron->gsfTrack().isNonnull() ) { _hElectronIp[NpdfID]->Fill( patElectron->gsfTrack()->dxy(thePrimaryEventVertex.position()) ,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID)); }
	_hElectronEoverP[NpdfID]->Fill(patElectron->eSuperClusterOverP(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        if(_UseHeepInfo) {
          heep::Ele theHeepElec(*patElectron);
	  _hElectronHoverEm[NpdfID]->Fill(theHeepElec.hOverE(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          if(patElectron->isEE()) {
	    _hElectronEESigmaIEtaIEta[NpdfID]->Fill(theHeepElec.sigmaIEtaIEta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectronEEDEta[NpdfID]->Fill(theHeepElec.dEtaIn(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectronEEDPhi[NpdfID]->Fill(theHeepElec.dPhiIn(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          }
          if(patElectron->isEB()) {
	    _hElectronEBSigmaIEtaIEta[NpdfID]->Fill(theHeepElec.sigmaIEtaIEta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectronEBDEta[NpdfID]->Fill(theHeepElec.dEtaIn(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectronEBDPhi[NpdfID]->Fill(theHeepElec.dPhiIn(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectronEB2by5Over5by5[NpdfID]->Fill(theHeepElec.scE2x5Max() / theHeepElec.scE5x5(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectronEB1by5Over5by5[NpdfID]->Fill(theHeepElec.scE1x5() / theHeepElec.scE5x5(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          }
        } else {
	  _hElectronHoverEm[NpdfID]->Fill(patElectron->hadronicOverEm(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          if(patElectron->isEE()) {
	    _hElectronEESigmaIEtaIEta[NpdfID]->Fill(patElectron->scSigmaIEtaIEta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectronEEDEta[NpdfID]->Fill(patElectron->deltaEtaSuperClusterTrackAtVtx(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectronEEDPhi[NpdfID]->Fill(patElectron->deltaPhiSuperClusterTrackAtVtx(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          }
          if(patElectron->isEB()) {
	    _hElectronEBSigmaIEtaIEta[NpdfID]->Fill(patElectron->scSigmaIEtaIEta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectronEBDEta[NpdfID]->Fill(patElectron->deltaEtaSuperClusterTrackAtVtx(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectronEBDPhi[NpdfID]->Fill(patElectron->deltaPhiSuperClusterTrackAtVtx(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectronEB2by5Over5by5[NpdfID]->Fill(patElectron->scE2x5Max() / patElectron->scE5x5(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectronEB1by5Over5by5[NpdfID]->Fill(patElectron->scE1x5() / patElectron->scE5x5(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          }
        }
        const Track* elTrack = (const reco::Track*)(patElectron->gsfTrack().get());
        const HitPattern& pInner = elTrack->trackerExpectedHitsInner();
	_hElectronMissingHits[NpdfID]->Fill(pInner.numberOfHits(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectronClassification[NpdfID]->Fill(patElectron->classification(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectronEcalDriven[NpdfID]->Fill(patElectron->ecalDrivenSeed(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectronTrackerDriven[NpdfID]->Fill(patElectron->trackerDrivenSeed(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
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
	if (!passRecoJetCuts((*patJet),theNumberOfJets-1)) continue;
	_hBJetDiscrByTrackCounting[NpdfID]->Fill(patJet->bDiscriminator("trackCountingHighEffBJetTags"),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBJetDiscrBySimpleSecondaryV[NpdfID]->Fill(patJet->bDiscriminator("simpleSecondaryVertexBJetTags"),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBJetDiscrByCombinedSecondaryV[NpdfID]->Fill(patJet->bDiscriminator("combinedSecondaryVertexBJetTags"),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        _hJetEnergy[NpdfID]->Fill(smearedJetMomentumVector.at(theNumberOfJets-1).energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hJetPt[NpdfID]->Fill(smearedJetPtEtaPhiMVector.at(theNumberOfJets-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hJetEta[NpdfID]->Fill(smearedJetPtEtaPhiMVector.at(theNumberOfJets-1).eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hJetPhi[NpdfID]->Fill(smearedJetPtEtaPhiMVector.at(theNumberOfJets-1).phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	nJets++;
      }
      _hNJet[NpdfID]->Fill(nJets,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
      _hMHT[NpdfID]->Fill(sqrt((sumpxForMht * sumpxForMht) + (sumpyForMht * sumpyForMht)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
      _hHT[NpdfID]->Fill(sumptForHt,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
      _hFirstLeadingJetPt[NpdfID]->Fill(leadingjetpt,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
      _hSecondLeadingJetPt[NpdfID]->Fill(secondleadingjetpt,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
    }
    
    // ------Topology Histograms
    double maximumLeptonPt = 0.0;
    double minimumLeptonEta = 9999.9;
    double maximumdR = 0.0;
    double minimumLeptonTrackIso = 9999.9;
    double minimumLeptonEcalIso = 9999.9;
    double maximumTauPt = 0.0;
    double maximumTauSeedPt = -1.0;
    double minimumTauEta = 9999.9;
    double minimumTauTrackIso = 9999.9;
    double minimumTauEcalIso = 9999.9;
    double minimumTauLeadNhits = 9999.9;
    double minimumTauH3x3OverP = 9999.9;
    double mimimumCosDphi = 1.0;
    double minimumOSLS = 9999.9;
    double maximumZeta = 0.0;
    double maximumLeptonMetMt = 0.0;
    if (_FillTopologyHists) {
      _hMet[NpdfID]->Fill(theMETVector.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
      if( ((_AnalyzeMuonForLeg1) && (_AnalyzeTauForLeg2)) || ((_AnalyzeMuonForLeg2) && (_AnalyzeTauForLeg1)) ) {
        double minimumIp = 999.9;
        double maximumSegComp = 0.0;
        double maximumCaloComp = 0.0;
        double maximumPionVeto = 0.0;
	int theNumberOfMuons = 0;
	for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
	  theNumberOfMuons++;
	  int theNumberOfTaus = 0;
	  for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
	    theNumberOfTaus++;
	    if ((passRecoTauCuts((*patTau),theNumberOfTaus-1)) &&
		(passRecoMuonCuts((*patMuon),theNumberOfMuons-1)) &&
		(passTopologyCuts((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1))) {
              //--- Best Lepton Pt
              if(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt() > maximumLeptonPt) {maximumLeptonPt = smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt();}
              //--- Best Lepton Eta
              if(fabs(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).eta()) < fabs(minimumLeptonEta)) {minimumLeptonEta = smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).eta();}
              //--- Best Lepton Track Iso
              if(patMuon->trackIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonTrackIsoTrkThreshold).first < minimumLeptonTrackIso) {minimumLeptonTrackIso = patMuon->trackIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonTrackIsoTrkThreshold).first;}
              //--- Best Lepton Ecal Iso
              if(patMuon->ecalIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonEcalIsoRecHitThreshold).first < minimumLeptonEcalIso) {minimumLeptonEcalIso = patMuon->ecalIsoDeposit()->depositAndCountWithin(_RecoMuonIsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuonEcalIsoRecHitThreshold).first;}
              //--- Best Pion Veto
	      if(muon::caloCompatibility(*patMuon) > maximumCaloComp) {maximumCaloComp = muon::caloCompatibility(*patMuon);}
	      if(muon::segmentCompatibility(*patMuon) > maximumSegComp) {maximumSegComp = muon::segmentCompatibility(*patMuon);}
              if(((_RecoMuonCaloCompCoefficient * muon::caloCompatibility(*patMuon)) + 
                  (_RecoMuonSegmCompCoefficient * muon::segmentCompatibility(*patMuon))) > maximumPionVeto) 
                {maximumPionVeto = (_RecoMuonCaloCompCoefficient * muon::caloCompatibility(*patMuon)) + 
                                   (_RecoMuonSegmCompCoefficient *muon::segmentCompatibility(*patMuon));}
              //--- Best Lepton Ip
	      const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
	      if ( patMuon->track().isNonnull() ) {
                if(fabs(patMuon->track()->dxy(thePrimaryEventVertex.position())) < fabs(minimumIp)) {minimumIp = patMuon->track()->dxy(thePrimaryEventVertex.position());}
              }
              //--- Best Tau Pt
              if(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() > maximumTauPt){
                if(reco::deltaR(maxPtMuonVector, smearedTauMomentumVector.at(theNumberOfTaus-1)) > _DiTauDeltaRCut)
                  {maximumTauPt = smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt();}
              }
              //--- Best Tau Eta
              if(fabs(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).eta()) < fabs(minimumTauEta)) {minimumTauEta = smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).eta();}
              //--- Best Tau Seed Pt, H3x3OverP, and lead track Nhits
              if (patTau->isCaloTau()) {
                if( (patTau->leadTrack().isNonnull()) ) {
                  if(patTau->leadTrack()->pt() > maximumTauSeedPt) {maximumTauSeedPt = patTau->leadTrack()->pt();}
                  if(patTau->hcal3x3OverPLead() < minimumTauH3x3OverP) {minimumTauH3x3OverP = patTau->hcal3x3OverPLead();}
                  if(patTau->leadTrack()->recHitsSize() < minimumTauLeadNhits) {minimumTauLeadNhits = patTau->leadTrack()->recHitsSize();}
                }
                if (_UseRecoTauDiscrByIsolationFlag) {
                  if(patTau->isolationTracksPtSum() < minimumTauTrackIso) {minimumTauTrackIso = patTau->isolationTracksPtSum();}
                  if(patTau->isolationECALhitsEtSum() < minimumTauEcalIso) {minimumTauEcalIso = patTau->isolationECALhitsEtSum();}
                } else {
                  if(CalculateTauTrackIsolation(*patTau).second < minimumTauTrackIso) {minimumTauTrackIso = CalculateTauTrackIsolation(*patTau).second;}
                  if(CalculateTauEcalIsolation(*patTau).second < minimumTauEcalIso) {minimumTauEcalIso = CalculateTauEcalIsolation(*patTau).second;}
                }
              } else {
                if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
                  if(patTau->leadPFChargedHadrCand()->pt() > maximumTauSeedPt) {maximumTauSeedPt = patTau->leadPFChargedHadrCand()->pt();}
                  if(patTau->hcal3x3OverPLead() < minimumTauH3x3OverP) {minimumTauH3x3OverP = patTau->hcal3x3OverPLead();}
	          if( (patTau->leadPFChargedHadrCand()->trackRef().isNonnull()) ) {
                    if(patTau->leadPFChargedHadrCand()->trackRef()->recHitsSize() < minimumTauLeadNhits) {minimumTauLeadNhits = patTau->leadPFChargedHadrCand()->trackRef()->recHitsSize();}
                  }
                }
                if (_UseRecoTauDiscrByIsolationFlag) {
                  if(patTau->isolationPFChargedHadrCandsPtSum() < minimumTauTrackIso) {minimumTauTrackIso = patTau->isolationPFChargedHadrCandsPtSum();}
                  if(patTau->isolationPFGammaCandsEtSum() < minimumTauEcalIso) {minimumTauEcalIso = patTau->isolationPFGammaCandsEtSum();}
                } else {
                  if(CalculateTauTrackIsolation(*patTau).second < minimumTauTrackIso) {minimumTauTrackIso = CalculateTauTrackIsolation(*patTau).second;}
                  if(CalculateTauEcalIsolation(*patTau).second < minimumTauEcalIso) {minimumTauEcalIso = CalculateTauEcalIsolation(*patTau).second;}
                }
              }
              //--- Best Lepton+Tau DeltaR
              if(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), smearedMuonMomentumVector.at(theNumberOfMuons-1)) > maximumdR) {maximumdR = reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), smearedMuonMomentumVector.at(theNumberOfMuons-1));}
              //--- Best Lepton+Tau CosDphi
              if(cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))) < mimimumCosDphi) {mimimumCosDphi = cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi())));}
	      _hMuonPtVsTauPt[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt(),smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuonTauDeltaR[NpdfID]->Fill(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), smearedMuonMomentumVector.at(theNumberOfMuons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuonTauDeltaPtDivSumPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt()) / (smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() + smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuonTauDeltaPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuonTauCosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuonMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuonMetDeltaPhiVsMuonTauCosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - theMETVector.phi())), cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      else {_hNotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
              //--- Best Lepton+Met Transverse mass
	      if(CalculateLeptonMetMt((*patMuon),theNumberOfMuons-1) > maximumLeptonMetMt) {maximumLeptonMetMt = CalculateLeptonMetMt((*patMuon),theNumberOfMuons-1);}
	      _hMuonMetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauMetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau),theNumberOfTaus-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              //--- Best Lepton+Tau OSLS
	      if(_UseTauSeedTrackForDiTauDiscrByOSLS) {
		if (patTau->isCaloTau()) {
		  if( (patTau->leadTrack().isNonnull()) ) {
		    _hMuonTauOSLS[NpdfID]->Fill(patMuon->charge() * patTau->leadTrack()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
                    if(patMuon->charge() * patTau->leadTrack()->charge() < minimumOSLS) {minimumOSLS = patMuon->charge() * patTau->leadTrack()->charge();}
		  }
		} else {
		  if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
		    _hMuonTauOSLS[NpdfID]->Fill(patMuon->charge() * patTau->leadPFChargedHadrCand()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
                    if(patMuon->charge() * patTau->leadPFChargedHadrCand()->charge() < minimumOSLS) {minimumOSLS = patMuon->charge() * patTau->leadPFChargedHadrCand()->charge();}
		  }
		}
	      } else {
                _hMuonTauOSLS[NpdfID]->Fill(patMuon->charge() * patTau->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
                if(patMuon->charge() * patTau->charge() < minimumOSLS) {minimumOSLS = patMuon->charge() * patTau->charge();}
              }
	      _hPZeta[NpdfID]->Fill(CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hPZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta2D[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta1D[NpdfID]->Fill((_PZetaCutCoefficient * CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1)) + 
				     (_PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              //--- Best Lepton+Tau+Met Zeta
	      if( (_PZetaCutCoefficient * CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1)) + 
                  (_PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1)) > maximumZeta ) {
                maximumZeta = (_PZetaCutCoefficient * CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1)) + 
                              (_PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1));
              }
	    }
	  }
	}
	_hBestMuonPt[NpdfID]->Fill(maximumLeptonPt,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestMuonEta[NpdfID]->Fill(minimumLeptonEta,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestMuonTrackIso[NpdfID]->Fill(minimumLeptonTrackIso,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestMuonEcalIso[NpdfID]->Fill(minimumLeptonEcalIso,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestMuonTauDeltaR[NpdfID]->Fill(maximumdR,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestMuonCaloCompatibility[NpdfID]->Fill(maximumCaloComp,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestMuonSegmentCompatibility[NpdfID]->Fill(maximumSegComp,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestMuonAntiPion[NpdfID]->Fill(maximumPionVeto,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestMuonIp[NpdfID]->Fill(minimumIp,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestTauJetPt[NpdfID]->Fill(maximumTauPt,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestTauJetEta[NpdfID]->Fill(minimumTauEta,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        _hBestTauJetSeedTrackPt[NpdfID]->Fill(maximumTauSeedPt,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        _hBestTauJetSumPtIsoTracks[NpdfID]->Fill(minimumTauTrackIso,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        _hBestTauJetSumPtIsoGammas[NpdfID]->Fill(minimumTauEcalIso,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestMuonTauCosDphi[NpdfID]->Fill(mimimumCosDphi,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestTauJetH3x3OverP[NpdfID]->Fill(minimumTauH3x3OverP,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestTauJetSeedTrackNhits[NpdfID]->Fill(minimumTauLeadNhits,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        _hBestMuonMetMt[NpdfID]->Fill(maximumLeptonMetMt,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        _hBestZeta1D[NpdfID]->Fill(maximumZeta,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        _hBestMuonTauOSLS[NpdfID]->Fill(minimumOSLS,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
      }
      if( ((_AnalyzeElectronForLeg1) && (_AnalyzeTauForLeg2)) || ((_AnalyzeElectronForLeg2) && (_AnalyzeTauForLeg1)) ) {
        double minimumHoverE = 9999.9;
        double minimumEESigmaIEtaIEta = 9999.9;
        double minimumEBSigmaIEtaIEta = 9999.9;
        double minimumEEDEta = 9999.9;
        double minimumEBDEta = 9999.9;
        double minimumEEDPhi = 9999.9;
        double minimumEBDPhi = 9999.9;
        double maximumEB2x2Over5x5 = 0.0;
        double maximumEB1x2Over5x5 = 0.0;
        double minimumMissingHits = 9999.9;
	int theNumberOfElectrons = 0;
	for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
	  theNumberOfElectrons++;
	  int theNumberOfTaus = 0;
	  for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
	    theNumberOfTaus++;
	    if ((passRecoTauCuts((*patTau),theNumberOfTaus-1)) &&
		(passRecoElectronCuts((*patElectron),theNumberOfElectrons-1)) &&
		(passTopologyCuts((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1))) {
              if(_UseHeepInfo) {
                heep::Ele theHeepElec(*patElectron);
                if(theHeepElec.hOverE() < minimumHoverE) {minimumHoverE = theHeepElec.hOverE();}
                if(patElectron->isEE()) {
                  if(theHeepElec.sigmaIEtaIEta() < minimumEESigmaIEtaIEta) {minimumEESigmaIEtaIEta = theHeepElec.sigmaIEtaIEta();}
                  if(fabs(theHeepElec.dEtaIn()) < fabs(minimumEEDEta)) {minimumEEDEta = theHeepElec.dEtaIn();}
                  if(fabs(theHeepElec.dPhiIn()) < fabs(minimumEEDPhi)) {minimumEEDPhi = theHeepElec.dPhiIn();}
                }
                if(patElectron->isEB()) {
                  if(theHeepElec.sigmaIEtaIEta() < minimumEBSigmaIEtaIEta) {minimumEBSigmaIEtaIEta = theHeepElec.sigmaIEtaIEta();}
                  if(fabs(theHeepElec.dEtaIn()) < fabs(minimumEBDEta)) {minimumEBDEta = theHeepElec.dEtaIn();}
                  if(fabs(theHeepElec.dPhiIn()) < fabs(minimumEBDPhi)) {minimumEBDPhi = theHeepElec.dPhiIn();}
                  if((theHeepElec.scE2x5Max() / theHeepElec.scE5x5()) > maximumEB2x2Over5x5) {maximumEB2x2Over5x5 = theHeepElec.scE2x5Max() / theHeepElec.scE5x5();}
                  if((theHeepElec.scE1x5() / theHeepElec.scE5x5()) > maximumEB1x2Over5x5) {maximumEB1x2Over5x5 = theHeepElec.scE1x5() / theHeepElec.scE5x5();}
                }
              } else {
                if(patElectron->hadronicOverEm() < minimumHoverE) {minimumHoverE = patElectron->hadronicOverEm();}
                if(patElectron->isEE()) {
                  if(patElectron->scSigmaIEtaIEta() < minimumEESigmaIEtaIEta) {minimumEESigmaIEtaIEta = patElectron->scSigmaIEtaIEta();}
                  if(fabs(patElectron->deltaEtaSuperClusterTrackAtVtx()) < fabs(minimumEEDEta)) {minimumEEDEta = patElectron->deltaEtaSuperClusterTrackAtVtx();}
                  if(fabs(patElectron->deltaPhiSuperClusterTrackAtVtx()) < fabs(minimumEEDPhi)) {minimumEEDPhi = patElectron->deltaPhiSuperClusterTrackAtVtx();}
                }
                if(patElectron->isEB()) {
                  if(patElectron->scSigmaIEtaIEta() < minimumEBSigmaIEtaIEta) {minimumEBSigmaIEtaIEta = patElectron->scSigmaIEtaIEta();}
                  if(fabs(patElectron->deltaEtaSuperClusterTrackAtVtx()) < fabs(minimumEBDEta)) {minimumEBDEta = patElectron->deltaEtaSuperClusterTrackAtVtx();}
                  if(fabs(patElectron->deltaPhiSuperClusterTrackAtVtx()) < fabs(minimumEBDPhi)) {minimumEBDPhi = patElectron->deltaPhiSuperClusterTrackAtVtx();}
                  if((patElectron->scE2x5Max() / patElectron->scE5x5()) > maximumEB2x2Over5x5) {maximumEB2x2Over5x5 = patElectron->scE2x5Max() / patElectron->scE5x5();}
                  if((patElectron->scE1x5() / patElectron->scE5x5()) > maximumEB1x2Over5x5) {maximumEB1x2Over5x5 = patElectron->scE1x5() / patElectron->scE5x5();}
                }
              }
              const Track* elTrack = (const reco::Track*)(patElectron->gsfTrack().get());
              const HitPattern& pInner = elTrack->trackerExpectedHitsInner();
              if(pInner.numberOfHits() < minimumMissingHits) {minimumMissingHits = pInner.numberOfHits();}
              if(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt() > maximumLeptonPt) {maximumLeptonPt = smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt();}
              if(fabs(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).eta()) < fabs(minimumLeptonEta)) {minimumLeptonEta = smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).eta();}
              if(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), smearedElectronMomentumVector.at(theNumberOfElectrons-1)) > maximumdR) {maximumdR = reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), smearedElectronMomentumVector.at(theNumberOfElectrons-1));}
              if(patElectron->dr04TkSumPt() < minimumLeptonTrackIso) {minimumLeptonTrackIso = patElectron->dr04TkSumPt();}
              if(patElectron->dr04EcalRecHitSumEt() < minimumLeptonEcalIso) {minimumLeptonEcalIso = patElectron->dr04EcalRecHitSumEt();}
              if(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() > maximumTauPt) {
                if(reco::deltaR(maxEtElectronVector, smearedTauMomentumVector.at(theNumberOfElectrons-1))>_DiTauDeltaRCut)
                  {maximumTauPt = smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt();}
              }
              if(fabs(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).eta()) < fabs(minimumTauEta)) {minimumTauEta = smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).eta();}
              if (patTau->isCaloTau()) {
                if( (patTau->leadTrack().isNonnull()) ) {
                  if(patTau->leadTrack()->pt() > maximumTauSeedPt) {maximumTauSeedPt = patTau->leadTrack()->pt();}
                  if(patTau->hcal3x3OverPLead() < minimumTauH3x3OverP) {minimumTauH3x3OverP = patTau->hcal3x3OverPLead();}
                  if(patTau->leadTrack()->recHitsSize() < minimumTauLeadNhits) {minimumTauLeadNhits = patTau->leadTrack()->recHitsSize();}
                }
                if (_UseRecoTauDiscrByIsolationFlag) {
                  if(patTau->isolationTracksPtSum() < minimumTauTrackIso) {minimumTauTrackIso = patTau->isolationTracksPtSum();}
                  if(patTau->isolationECALhitsEtSum() < minimumTauEcalIso) {minimumTauEcalIso = patTau->isolationECALhitsEtSum();}
                } else {
                  if(CalculateTauTrackIsolation(*patTau).second < minimumTauTrackIso) {minimumTauTrackIso = CalculateTauTrackIsolation(*patTau).second;}
                  if(CalculateTauEcalIsolation(*patTau).second < minimumTauEcalIso) {minimumTauEcalIso = CalculateTauEcalIsolation(*patTau).second;}
                }
              } else {
                if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
                  if(patTau->leadPFChargedHadrCand()->pt() > maximumTauSeedPt) {maximumTauSeedPt = patTau->leadPFChargedHadrCand()->pt();}
                  if(patTau->hcal3x3OverPLead() < minimumTauH3x3OverP) {minimumTauH3x3OverP = patTau->hcal3x3OverPLead();}
	          if( (patTau->leadPFChargedHadrCand()->trackRef().isNonnull()) ) {
                    if(patTau->leadPFChargedHadrCand()->trackRef()->recHitsSize() < minimumTauLeadNhits) {minimumTauLeadNhits = patTau->leadPFChargedHadrCand()->trackRef()->recHitsSize();}
                  }
                }
                if (_UseRecoTauDiscrByIsolationFlag) {
                  if(patTau->isolationPFChargedHadrCandsPtSum() < minimumTauTrackIso) {minimumTauTrackIso = patTau->isolationPFChargedHadrCandsPtSum();}
                  if(patTau->isolationPFGammaCandsEtSum() < minimumTauEcalIso) {minimumTauEcalIso = patTau->isolationPFGammaCandsEtSum();}
                } else {
                  if(CalculateTauTrackIsolation(*patTau).second < minimumTauTrackIso) {minimumTauTrackIso = CalculateTauTrackIsolation(*patTau).second;}
                  if(CalculateTauEcalIsolation(*patTau).second < minimumTauEcalIso) {minimumTauEcalIso = CalculateTauEcalIsolation(*patTau).second;}
                }
              }
              if(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))) < mimimumCosDphi) {mimimumCosDphi = cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi())));}
              _hElectronIsZee[NpdfID]->Fill(isZee(smearedElectronMomentumVector.at(theNumberOfElectrons-1)).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              _hElectronPtVsTauPt[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt(),smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectronTauDeltaR[NpdfID]->Fill(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), smearedElectronMomentumVector.at(theNumberOfElectrons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectronTauDeltaPtDivSumPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt()) / (smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() + smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectronTauDeltaPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectronTauCosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectronMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectronMetDeltaPhiVsElectronTauCosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - theMETVector.phi())), cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauMetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      else {_hNotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      if(CalculateLeptonMetMt((*patElectron),theNumberOfElectrons-1) > maximumLeptonMetMt) {maximumLeptonMetMt = CalculateLeptonMetMt((*patElectron),theNumberOfElectrons-1);}
	      _hElectronMetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTauMetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau),theNumberOfTaus-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      if(_UseTauSeedTrackForDiTauDiscrByOSLS) {
		if (patTau->isCaloTau()) {
		  if( (patTau->leadTrack().isNonnull()) ) {
		    _hElectronTauOSLS[NpdfID]->Fill(patElectron->charge() * patTau->leadTrack()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
                    if((patElectron->charge() * patTau->leadTrack()->charge()) < minimumOSLS) {minimumOSLS = patElectron->charge() * patTau->leadTrack()->charge();}
		  }
		} else {
		  if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
		    _hElectronTauOSLS[NpdfID]->Fill(patElectron->charge() * patTau->leadPFChargedHadrCand()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
                    if((patElectron->charge() * patTau->leadPFChargedHadrCand()->charge()) < minimumOSLS) {minimumOSLS = patElectron->charge() * patTau->leadPFChargedHadrCand()->charge();}
		  }
		}
	      } else {
                _hElectronTauOSLS[NpdfID]->Fill(patElectron->charge() * patTau->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
                if((patElectron->charge() * patTau->charge()) < minimumOSLS) {minimumOSLS = patElectron->charge() * patTau->charge();}
              }
	      _hPZeta[NpdfID]->Fill(CalculatePZeta((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hPZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta2D[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),CalculatePZeta((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta1D[NpdfID]->Fill((_PZetaCutCoefficient * CalculatePZeta((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1)) + (_PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              //--- Best Lepton+Tau+Met Zeta
	      if( (_PZetaCutCoefficient * CalculatePZeta((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1)) + 
                  (_PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1)) > maximumZeta ) {
                maximumZeta = (_PZetaCutCoefficient * CalculatePZeta((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1)) + 
                              (_PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1));
              }
	    }
	  }
	}
	_hBestElectronPt[NpdfID]->Fill(maximumLeptonPt,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestElectronEta[NpdfID]->Fill(minimumLeptonEta,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestElectronTrackIso[NpdfID]->Fill(minimumLeptonTrackIso,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestElectronEcalIso[NpdfID]->Fill(minimumLeptonEcalIso,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestElectronTauDeltaR[NpdfID]->Fill(maximumdR,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestTauJetPt[NpdfID]->Fill(maximumTauPt,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestTauJetEta[NpdfID]->Fill(minimumTauEta,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        _hBestTauJetSeedTrackPt[NpdfID]->Fill(maximumTauSeedPt,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        _hBestTauJetSumPtIsoTracks[NpdfID]->Fill(minimumTauTrackIso,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        _hBestTauJetSumPtIsoGammas[NpdfID]->Fill(minimumTauEcalIso,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestElectronHoverEm[NpdfID]->Fill(minimumHoverE,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestElectronEESigmaIEtaIEta[NpdfID]->Fill(minimumEESigmaIEtaIEta,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestElectronEEDEta[NpdfID]->Fill(minimumEEDEta,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestElectronEEDPhi[NpdfID]->Fill(minimumEEDPhi,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestElectronEBSigmaIEtaIEta[NpdfID]->Fill(minimumEBSigmaIEtaIEta,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestElectronEBDEta[NpdfID]->Fill(minimumEBDEta,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestElectronEBDPhi[NpdfID]->Fill(minimumEBDPhi,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestElectronEB2by5Over5by5[NpdfID]->Fill(maximumEB2x2Over5x5,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestElectronEB1by5Over5by5[NpdfID]->Fill(maximumEB1x2Over5x5,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestElectronMissingHits[NpdfID]->Fill(minimumMissingHits,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestElectronTauCosDphi[NpdfID]->Fill(mimimumCosDphi,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestTauJetH3x3OverP[NpdfID]->Fill(minimumTauH3x3OverP,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBestTauJetSeedTrackNhits[NpdfID]->Fill(minimumTauLeadNhits,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        _hBestElectronMetMt[NpdfID]->Fill(maximumLeptonMetMt,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        _hBestZeta1D[NpdfID]->Fill(maximumZeta,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        _hBestElectronTauOSLS[NpdfID]->Fill(minimumOSLS,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
      }
      if( ((_AnalyzeMuonForLeg1) && (_AnalyzeMuonForLeg2)) ) {
	int theNumberOfMuons1 = 0;
	for ( pat::MuonCollection::const_iterator patMuon1 = _patMuons->begin();patMuon1 != _patMuons->end(); ++patMuon1 ) {
	  theNumberOfMuons1++;
	  int theNumberOfMuons2 = 0;
	  for ( pat::MuonCollection::const_iterator patMuon2 = _patMuons->begin();patMuon2 != _patMuons->end(); ++patMuon2 ) {
	    theNumberOfMuons2++;
	    if ((passRecoMuonCuts((*patMuon1),theNumberOfMuons1 - 1)) && 
		(passRecoMuonCuts((*patMuon2),theNumberOfMuons2 - 1)) && 
		(passTopologyCuts((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1)) &&
                (theNumberOfMuons2 >= theNumberOfMuons1)) {
              _hMuon1PtVsMuon2Pt[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons1 - 1).pt(),smearedMuonPtEtaPhiMVector.at(theNumberOfMuons2 - 1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuon1Muon2DeltaR[NpdfID]->Fill(reco::deltaR(smearedMuonMomentumVector.at(theNumberOfMuons1 - 1), smearedMuonMomentumVector.at(theNumberOfMuons2 - 1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuon1Muon2DeltaPtDivSumPt[NpdfID]->Fill((smearedMuonPtEtaPhiMVector.at(theNumberOfMuons1 - 1).pt() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons2 - 1).pt()) / (smearedMuonPtEtaPhiMVector.at(theNumberOfMuons1 - 1).pt() + smearedMuonPtEtaPhiMVector.at(theNumberOfMuons2 - 1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuon1Muon2DeltaPt[NpdfID]->Fill((smearedMuonPtEtaPhiMVector.at(theNumberOfMuons1 - 1).pt() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons2 - 1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuon1Muon2CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons1 - 1).phi() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuon1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons1 - 1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuon2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons2 - 1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuon1MetDeltaPhiVsMuon1Muon2CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons1 - 1).phi() - theMETVector.phi())), cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons1 - 1).phi() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      if(CalculateThe4Momentum((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1).first) {_hReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      else {_hNotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      _hMuon1MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patMuon1),theNumberOfMuons1 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuon2MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patMuon2),theNumberOfMuons2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hMuon1Muon2OSLS[NpdfID]->Fill(patMuon1->charge() * patMuon2->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hPZeta[NpdfID]->Fill(CalculatePZeta((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hPZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta2D[NpdfID]->Fill(CalculatePZetaVis((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1),CalculatePZeta((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta1D[NpdfID]->Fill((_PZetaCutCoefficient * CalculatePZeta((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1)) + (_PZetaVisCutCoefficient * CalculatePZetaVis((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
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
	    if ((passRecoElectronCuts((*patElectron1),theNumberOfElectrons1 - 1)) && 
		(passRecoElectronCuts((*patElectron2),theNumberOfElectrons2 - 1)) && 
		(passTopologyCuts((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1)) &&
                (theNumberOfElectrons2 >= theNumberOfElectrons1)) {
	      _hElectron1PtVsElectron2Pt[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons1 - 1).pt(),smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons2 - 1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectron1Electron2DeltaR[NpdfID]->Fill(reco::deltaR(smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1), smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectron1Electron2DeltaPtDivSumPt[NpdfID]->Fill((smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons1 - 1).pt() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons2 - 1).pt()) / (smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons1 - 1).pt() + smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons2 - 1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectron1Electron2DeltaPt[NpdfID]->Fill((smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons1 - 1).pt() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons2 - 1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectron1Electron2CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons1 - 1).phi() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectron1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons1 - 1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectron2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons2 - 1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectron1MetDeltaPhiVsElectron1Electron2CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons1 - 1).phi() - theMETVector.phi())), cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons1 - 1).phi() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      if(CalculateThe4Momentum((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1).first) {_hReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      else {_hNotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      _hElectron1MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patElectron1),theNumberOfElectrons1 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectron2MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patElectron2),theNumberOfElectrons2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hElectron1Electron2OSLS[NpdfID]->Fill(patElectron1->charge() * patElectron2->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hPZeta[NpdfID]->Fill(CalculatePZeta((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hPZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta2D[NpdfID]->Fill(CalculatePZetaVis((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1),CalculatePZeta((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta1D[NpdfID]->Fill((_PZetaCutCoefficient * CalculatePZeta((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1)) + (_PZetaVisCutCoefficient * CalculatePZetaVis((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
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
	    if ((passRecoTauCuts((*patTau1),theNumberOfTaus1 - 1)) &&
		(passRecoTauCuts((*patTau2),theNumberOfTaus2 - 1)) &&
		(passTopologyCuts((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1)) &&
                (theNumberOfTaus2 >= theNumberOfTaus1)) {
	      _hTau1PtVsTau2Pt[NpdfID]->Fill(smearedTauPtEtaPhiMVector.at(theNumberOfTaus1 - 1).pt(),smearedTauPtEtaPhiMVector.at(theNumberOfTaus2 - 1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTau1Tau2DeltaR[NpdfID]->Fill(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus1 - 1), smearedTauMomentumVector.at(theNumberOfTaus2 - 1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTau1Tau2DeltaPtDivSumPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus1 - 1).pt() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus2 - 1).pt()) / (smearedTauPtEtaPhiMVector.at(theNumberOfTaus1 - 1).pt() + smearedTauPtEtaPhiMVector.at(theNumberOfTaus2 - 1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTau1Tau2DeltaPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus1 - 1).pt() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus2 - 1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      if(_UseTauSeedTrackForDiTauDiscrByOSLS) {
		if (patTau1->isCaloTau()) {
		  if( (patTau1->leadTrack().isNonnull()) && (patTau2->leadTrack().isNonnull()) ) {
		    _hTau1Tau2OSLS[NpdfID]->Fill(patTau1->leadTrack()->charge() * patTau2->leadTrack()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		  }
		} else {
		  if( (patTau1->leadPFChargedHadrCand().isNonnull()) && (patTau2->leadPFChargedHadrCand().isNonnull()) ) {
                    _hTau1Tau2OSLS[NpdfID]->Fill(patTau1->leadPFChargedHadrCand()->charge() * patTau2->leadPFChargedHadrCand()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
                    if((patTau1->leadPFChargedHadrCand()->charge() * patTau2->leadPFChargedHadrCand()->charge())<0) {
	              if(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).first) {_hReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	              else {_hNotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                    } else {
	              if(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).first) {_hReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	              else {_hNotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                    }
		  }
		}
	      } else {_hTau1Tau2OSLS[NpdfID]->Fill(patTau1->charge() * patTau2->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      _hTau1Tau2CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus1 - 1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTau1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus1 - 1).phi() -  theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTau2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus2 - 1).phi() -  theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTau1MetDeltaPhiVsTau1Tau2CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus1 - 1).phi() -  theMETVector.phi())), cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus1 - 1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      if(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).first) {_hReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      else {_hNotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	      _hTau1MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau1),theNumberOfTaus1 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hTau2MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau2),theNumberOfTaus2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hPZeta[NpdfID]->Fill(CalculatePZeta((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hPZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta2D[NpdfID]->Fill(CalculatePZetaVis((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1),CalculatePZeta((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	      _hZeta1D[NpdfID]->Fill((_PZetaCutCoefficient * CalculatePZeta((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1)) + (_PZetaVisCutCoefficient * CalculatePZetaVis((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    }
	  }
	}
      }
    }
    int numberJets1 = 0;
    for ( pat::JetCollection::const_iterator patJet1 = _patJets->begin();patJet1 != _patJets->end(); ++patJet1 ) {
      numberJets1++;
      if( (numberJets1 == theLeadingJetIndex) && (passRecoFirstLeadingJetCuts((*patJet1),numberJets1 - 1)) ) {
        int numberJets2 = 0;
        for ( pat::JetCollection::const_iterator patJet2 = _patJets->begin();patJet2 != _patJets->end(); ++patJet2 ) {
          numberJets2++;
          if( (numberJets2 == theSecondLeadingJetIndex) && (passRecoSecondLeadingJetCuts((*patJet2),numberJets2 - 1)) ) {
            if(passSusyTopologyCuts(numberJets1 - 1,numberJets2 - 1)) {
              double dphi1;
              double dphi2;
              double r1;
              double r2;
              double alpha;
              dphi1 = TMath::Abs(normalizedPhi(smearedJetPtEtaPhiMVector.at(numberJets1 - 1).phi() - theMETVector.phi()));
              dphi2 = TMath::Abs(normalizedPhi(smearedJetPtEtaPhiMVector.at(numberJets2 - 1).phi() - theMETVector.phi()));
              r1 = sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0) );
              r2 = sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0) );
              double px = smearedJetMomentumVector.at(numberJets1 - 1).px() + smearedJetMomentumVector.at(numberJets2 - 1).px();
              double py = smearedJetMomentumVector.at(numberJets1 - 1).py() + smearedJetMomentumVector.at(numberJets2 - 1).py();
              double pz = smearedJetMomentumVector.at(numberJets1 - 1).pz() + smearedJetMomentumVector.at(numberJets2 - 1).pz();
              double e = smearedJetMomentumVector.at(numberJets1 - 1).energy() + smearedJetMomentumVector.at(numberJets2 - 1).energy();
              reco::Candidate::LorentzVector Susy_LorentzVect(px, py, pz, e);
              if(Susy_LorentzVect.M() > 0) {alpha = smearedJetPtEtaPhiMVector.at(numberJets2 - 1).pt() / Susy_LorentzVect.M();}
              else {alpha = 999999999.0;}
              _hR1[NpdfID]->Fill(r1,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              _hR2[NpdfID]->Fill(r2,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              _hDphi1[NpdfID]->Fill(dphi1,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              _hDphi2[NpdfID]->Fill(dphi2,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              _hAlpha[NpdfID]->Fill(alpha,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
            }
          }
        }
      }
    }
  }
}

void HiMassTauAnalysis::getCollections(const Event& iEvent, const EventSetup& iSetup) {
  iEvent.getByLabel(_RecoTauSource, _patTaus);
  iEvent.getByLabel(_RecoMuonSource, _patMuons);
  if(_GenParticleSource.label() != "") { iEvent.getByLabel(_GenParticleSource, _genParticles); }
  iEvent.getByLabel(_RecoElectronSource, _patElectrons);
  iEvent.getByLabel(_RecoJetSource, _patJets);
  iEvent.getByLabel(_RecoMetSource, _patMETs);
  iEvent.getByLabel(_RecoVertexSource, _primaryEventVertexCollection);
  iEvent.getByLabel(_RecoTriggerSource, _triggerResults);
  iEvent.getByLabel( "hltTriggerSummaryAOD", handleTriggerEvent );
//  iEvent.getByLabel(_RecoParticleFlowSource, _pflow);
  iEvent.getByLabel("hpsPFTauProducer", _hpsTau);
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
    i++;
    if(reco::deltaR(theObject, smearedElectronMomentumVector.at(i - 1)) < _DiTauDeltaRCut) continue;
    if(theObject == smearedElectronMomentumVector.at(i - 1))continue;
    reco::Candidate::LorentzVector The_LorentzVect = theObject + smearedElectronMomentumVector.at(i - 1);
    zeePtAsymmetry = (theObject.pt() - smearedElectronPtEtaPhiMVector.at(i - 1).pt())/(theObject.pt() + smearedElectronPtEtaPhiMVector.at(i - 1).pt());
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

//-----Matching Electrons to generator level objects
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

//-----Matching Muons to generator level objects
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

//-----Matching Taus to generator level objects
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
double HiMassTauAnalysis::CalculatePZeta(const pat::Tau& patTau, int nobjT, const pat::Muon& patMuon, int nobjM) {
  double zetaX = cos(smearedTauPtEtaPhiMVector.at(nobjT).phi()) + cos(smearedMuonPtEtaPhiMVector.at(nobjM).phi());
  double zetaY = sin(smearedTauPtEtaPhiMVector.at(nobjT).phi()) + sin(smearedMuonPtEtaPhiMVector.at(nobjM).phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = smearedTauMomentumVector.at(nobjT).px() + smearedMuonMomentumVector.at(nobjM).px();
  double visPy = smearedTauMomentumVector.at(nobjT).py() + smearedMuonMomentumVector.at(nobjM).py();
  double px = visPx + theMETVector.px();
  double py = visPy + theMETVector.py();
  double pZeta = px*zetaX + py*zetaY;
  return pZeta;
}
double HiMassTauAnalysis::CalculatePZeta(const pat::Tau& patTau, int nobjT, const pat::Electron& patElectron, int nobjM) {
  double zetaX = cos(smearedTauPtEtaPhiMVector.at(nobjT).phi()) + cos(smearedElectronPtEtaPhiMVector.at(nobjM).phi());
  double zetaY = sin(smearedTauPtEtaPhiMVector.at(nobjT).phi()) + sin(smearedElectronPtEtaPhiMVector.at(nobjM).phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = smearedTauMomentumVector.at(nobjT).px() + smearedElectronMomentumVector.at(nobjM).px();
  double visPy = smearedTauMomentumVector.at(nobjT).py() + smearedElectronMomentumVector.at(nobjM).py();
  double px = visPx + theMETVector.px();
  double py = visPy + theMETVector.py();
  double pZeta = px*zetaX + py*zetaY;
  return pZeta;
}
double HiMassTauAnalysis::CalculatePZeta(const pat::Electron& patElectron, int nobjE, const pat::Muon& patMuon, int nobjM) {
  double zetaX = cos(smearedMuonPtEtaPhiMVector.at(nobjM).phi()) + cos(smearedElectronPtEtaPhiMVector.at(nobjE).phi());
  double zetaY = sin(smearedMuonPtEtaPhiMVector.at(nobjM).phi()) + sin(smearedElectronPtEtaPhiMVector.at(nobjE).phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = smearedMuonMomentumVector.at(nobjM).px() + smearedElectronMomentumVector.at(nobjE).px();
  double visPy = smearedMuonMomentumVector.at(nobjM).py() + smearedElectronMomentumVector.at(nobjE).py();
  double px = visPx + theMETVector.px();
  double py = visPy + theMETVector.py();
  double pZeta = px*zetaX + py*zetaY;
  return pZeta;
}
double HiMassTauAnalysis::CalculatePZeta(const pat::Muon& patMuon1, int nobj1, const pat::Muon& patMuon2, int nobj2) {
  double zetaX = cos(smearedMuonPtEtaPhiMVector.at(nobj1).phi()) + cos(smearedMuonPtEtaPhiMVector.at(nobj2).phi());
  double zetaY = sin(smearedMuonPtEtaPhiMVector.at(nobj1).phi()) + sin(smearedMuonPtEtaPhiMVector.at(nobj2).phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = smearedMuonMomentumVector.at(nobj1).px() + smearedMuonMomentumVector.at(nobj2).px();
  double visPy = smearedMuonMomentumVector.at(nobj1).py() + smearedMuonMomentumVector.at(nobj2).py();
  double px = visPx + theMETVector.px();
  double py = visPy + theMETVector.py();
  double pZeta = px*zetaX + py*zetaY;
  return pZeta;
}
double HiMassTauAnalysis::CalculatePZeta(const pat::Electron& patElectron1,  int nobj1, const pat::Electron& patElectron2,  int nobj2) {
  double zetaX = cos(smearedElectronPtEtaPhiMVector.at(nobj1).phi()) + cos(smearedElectronPtEtaPhiMVector.at(nobj2).phi());
  double zetaY = sin(smearedElectronPtEtaPhiMVector.at(nobj1).phi()) + sin(smearedElectronPtEtaPhiMVector.at(nobj2).phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = smearedElectronMomentumVector.at(nobj1).px() + smearedElectronMomentumVector.at(nobj2).px();
  double visPy = smearedElectronMomentumVector.at(nobj1).py() + smearedElectronMomentumVector.at(nobj2).py();
  double px = visPx + theMETVector.px();
  double py = visPy + theMETVector.py();
  double pZeta = px*zetaX + py*zetaY;
  return pZeta;
}
double HiMassTauAnalysis::CalculatePZeta(const pat::Tau& patTau1,  int nobj1, const pat::Tau& patTau2,  int nobj2) {
  double zetaX = cos(smearedTauPtEtaPhiMVector.at(nobj1).phi()) + cos(smearedTauPtEtaPhiMVector.at(nobj2).phi());
  double zetaY = sin(smearedTauPtEtaPhiMVector.at(nobj1).phi()) + sin(smearedTauPtEtaPhiMVector.at(nobj2).phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = smearedTauMomentumVector.at(nobj1).px() + smearedTauMomentumVector.at(nobj2).px();
  double visPy = smearedTauMomentumVector.at(nobj1).py() + smearedTauMomentumVector.at(nobj2).py();
  double px = visPx + theMETVector.px();
  double py = visPy + theMETVector.py();
  double pZeta = px*zetaX + py*zetaY;
  return pZeta;
}
double HiMassTauAnalysis::CalculatePZetaVis(const pat::Tau& patTau, int nobjT, const pat::Muon& patMuon, int nobjM) {
  double zetaX = cos(smearedTauPtEtaPhiMVector.at(nobjT).phi()) + cos(smearedMuonPtEtaPhiMVector.at(nobjM).phi());
  double zetaY = sin(smearedTauPtEtaPhiMVector.at(nobjT).phi()) + sin(smearedMuonPtEtaPhiMVector.at(nobjM).phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = smearedTauMomentumVector.at(nobjT).px() + smearedMuonMomentumVector.at(nobjM).px();
  double visPy = smearedTauMomentumVector.at(nobjT).py() + smearedMuonMomentumVector.at(nobjM).py();
  double pZetaVis = visPx*zetaX + visPy*zetaY;
  return pZetaVis;
}
double HiMassTauAnalysis::CalculatePZetaVis(const pat::Tau& patTau, int nobjT, const pat::Electron& patElectron, int nobjM) {
  double zetaX = cos(smearedTauPtEtaPhiMVector.at(nobjT).phi()) + cos(smearedElectronPtEtaPhiMVector.at(nobjM).phi());
  double zetaY = sin(smearedTauPtEtaPhiMVector.at(nobjT).phi()) + sin(smearedElectronPtEtaPhiMVector.at(nobjM).phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = smearedTauMomentumVector.at(nobjT).px() + smearedElectronMomentumVector.at(nobjM).px();
  double visPy = smearedTauMomentumVector.at(nobjT).py() + smearedElectronMomentumVector.at(nobjM).py();
  double pZetaVis = visPx*zetaX + visPy*zetaY;
  return pZetaVis;
}
double HiMassTauAnalysis::CalculatePZetaVis(const pat::Electron& patElectron, int nobjE, const pat::Muon& patMuon, int nobjM) {
  double zetaX = cos(smearedMuonPtEtaPhiMVector.at(nobjM).phi()) + cos(smearedElectronPtEtaPhiMVector.at(nobjE).phi());
  double zetaY = sin(smearedMuonPtEtaPhiMVector.at(nobjM).phi()) + sin(smearedElectronPtEtaPhiMVector.at(nobjE).phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = smearedMuonMomentumVector.at(nobjM).px() + smearedElectronMomentumVector.at(nobjE).px();
  double visPy = smearedMuonMomentumVector.at(nobjM).py() + smearedElectronMomentumVector.at(nobjE).py();
  double pZetaVis = visPx*zetaX + visPy*zetaY;
  return pZetaVis;
}
double HiMassTauAnalysis::CalculatePZetaVis(const pat::Muon& patMuon1, int nobj1, const pat::Muon& patMuon2, int nobj2) {
  double zetaX = cos(smearedMuonPtEtaPhiMVector.at(nobj1).phi()) + cos(smearedMuonPtEtaPhiMVector.at(nobj2).phi());
  double zetaY = sin(smearedMuonPtEtaPhiMVector.at(nobj1).phi()) + sin(smearedMuonPtEtaPhiMVector.at(nobj2).phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = smearedMuonMomentumVector.at(nobj1).px() + smearedMuonMomentumVector.at(nobj2).px();
  double visPy = smearedMuonMomentumVector.at(nobj1).py() + smearedMuonMomentumVector.at(nobj2).py();
  double pZetaVis = visPx*zetaX + visPy*zetaY;
  return pZetaVis;
}
double HiMassTauAnalysis::CalculatePZetaVis(const pat::Electron& patElectron1, int nobj1, const pat::Electron& patElectron2, int nobj2) {
  double zetaX = cos(smearedElectronPtEtaPhiMVector.at(nobj1).phi()) + cos(smearedElectronPtEtaPhiMVector.at(nobj2).phi());
  double zetaY = sin(smearedElectronPtEtaPhiMVector.at(nobj1).phi()) + sin(smearedElectronPtEtaPhiMVector.at(nobj2).phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = smearedElectronMomentumVector.at(nobj1).px() + smearedElectronMomentumVector.at(nobj2).px();
  double visPy = smearedElectronMomentumVector.at(nobj1).py() + smearedElectronMomentumVector.at(nobj2).py();
  double pZetaVis = visPx*zetaX + visPy*zetaY;
  return pZetaVis;
}
double HiMassTauAnalysis::CalculatePZetaVis(const pat::Tau& patTau1, int nobj1, const pat::Tau& patTau2, int nobj2) {
  double zetaX = cos(smearedTauPtEtaPhiMVector.at(nobj1).phi()) + cos(smearedTauPtEtaPhiMVector.at(nobj2).phi());
  double zetaY = sin(smearedTauPtEtaPhiMVector.at(nobj1).phi()) + sin(smearedTauPtEtaPhiMVector.at(nobj2).phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = smearedTauMomentumVector.at(nobj1).px() + smearedTauMomentumVector.at(nobj2).px();
  double visPy = smearedTauMomentumVector.at(nobj1).py() + smearedTauMomentumVector.at(nobj2).py();
  double pZetaVis = visPx*zetaX + visPy*zetaY;
  return pZetaVis;
}

//-----Calculate mass reco variables
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::CalculateThe4Momentum(const pat::Tau& patTau, int nobjT, const pat::Electron& patElectron, int nobjE) {
  if(_UseVectorSumOfVisProductsAndMetMassReco) {
    double px = smearedTauMomentumVector.at(nobjT).px() + smearedElectronMomentumVector.at(nobjE).px() + theMETVector.px();
    double py = smearedTauMomentumVector.at(nobjT).py() + smearedElectronMomentumVector.at(nobjE).py() + theMETVector.py();
    double pz = smearedTauMomentumVector.at(nobjT).pz() + smearedElectronMomentumVector.at(nobjE).pz();
    double e = smearedTauMomentumVector.at(nobjT).energy() + smearedElectronMomentumVector.at(nobjE).energy() + TMath::Sqrt((theMETVector.px() * theMETVector.px()) + (theMETVector.py() * theMETVector.py()));
    reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
    pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
    return MassRecoInformation;
  } else if(_UseCollinerApproxMassReco) {
    double x1_numerator = (smearedTauMomentumVector.at(nobjT).px() * smearedElectronMomentumVector.at(nobjE).py()) - (smearedElectronMomentumVector.at(nobjE).px() * smearedTauMomentumVector.at(nobjT).py());
    double x1_denominator = (smearedElectronMomentumVector.at(nobjE).py() * (smearedTauMomentumVector.at(nobjT).px() + theMETVector.px())) - (smearedElectronMomentumVector.at(nobjE).px() * (smearedTauMomentumVector.at(nobjT).py() + theMETVector.py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (smearedTauMomentumVector.at(nobjT).px() * (smearedElectronMomentumVector.at(nobjE).py() + theMETVector.py())) - (smearedTauMomentumVector.at(nobjT).py() * (smearedElectronMomentumVector.at(nobjE).px() + theMETVector.px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    if ( (x1 > 0. && x1 < 1.) &&
	 (x2 > 0. && x2 < 1.) ) {
      reco::Candidate::LorentzVector The_LorentzVect = (smearedTauMomentumVector.at(nobjT) / x1) + (smearedElectronMomentumVector.at(nobjE) / x2);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else {
      double px = smearedTauMomentumVector.at(nobjT).px() + smearedElectronMomentumVector.at(nobjE).px() + theMETVector.px();
      double py = smearedTauMomentumVector.at(nobjT).py() + smearedElectronMomentumVector.at(nobjE).py() + theMETVector.py();
      double pz = smearedTauMomentumVector.at(nobjT).pz() + smearedElectronMomentumVector.at(nobjE).pz();
      double e = smearedTauMomentumVector.at(nobjT).energy() + smearedElectronMomentumVector.at(nobjE).energy() + TMath::Sqrt((theMETVector.px() * theMETVector.px()) + (theMETVector.py() * theMETVector.py()));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
      return MassRecoInformation;
    }
  } else {
    reco::Candidate::LorentzVector The_LorentzVect = smearedTauMomentumVector.at(nobjT) + smearedElectronMomentumVector.at(nobjE);
    pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
    return MassRecoInformation;
  }
}
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::CalculateThe4Momentum(const pat::Tau& patTau, int nobjT, const pat::Muon& patMuon, int nobjE) {
  if(_UseVectorSumOfVisProductsAndMetMassReco) {
    double px = smearedTauMomentumVector.at(nobjT).px() + smearedMuonMomentumVector.at(nobjE).px() + theMETVector.px();
    double py = smearedTauMomentumVector.at(nobjT).py() + smearedMuonMomentumVector.at(nobjE).py() + theMETVector.py();
    double pz = smearedTauMomentumVector.at(nobjT).pz() + smearedMuonMomentumVector.at(nobjE).pz();
    double e = smearedTauMomentumVector.at(nobjT).energy() + smearedMuonMomentumVector.at(nobjE).energy() + TMath::Sqrt((theMETVector.px() * theMETVector.px()) + (theMETVector.py() * theMETVector.py()));
    reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
    pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
    return MassRecoInformation;
  } else if(_UseCollinerApproxMassReco) {
    double x1_numerator = (smearedTauMomentumVector.at(nobjT).px() * smearedMuonMomentumVector.at(nobjE).py()) - (smearedMuonMomentumVector.at(nobjE).px() * smearedTauMomentumVector.at(nobjT).py());
    double x1_denominator = (smearedMuonMomentumVector.at(nobjE).py() * (smearedTauMomentumVector.at(nobjT).px() + theMETVector.px())) - (smearedMuonMomentumVector.at(nobjE).px() * (smearedTauMomentumVector.at(nobjT).py() + theMETVector.py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (smearedTauMomentumVector.at(nobjT).px() * (smearedMuonMomentumVector.at(nobjE).py() + theMETVector.py())) - (smearedTauMomentumVector.at(nobjT).py() * (smearedMuonMomentumVector.at(nobjE).px() + theMETVector.px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    if ( (x1 > 0. && x1 < 1.) &&
	 (x2 > 0. && x2 < 1.) ) {
      reco::Candidate::LorentzVector The_LorentzVect = (smearedTauMomentumVector.at(nobjT) / x1) + (smearedMuonMomentumVector.at(nobjE) / x2);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else {
      double px = smearedTauMomentumVector.at(nobjT).px() + smearedMuonMomentumVector.at(nobjE).px() + theMETVector.px();
      double py = smearedTauMomentumVector.at(nobjT).py() + smearedMuonMomentumVector.at(nobjE).py() + theMETVector.py();
      double pz = smearedTauMomentumVector.at(nobjT).pz() + smearedMuonMomentumVector.at(nobjE).pz();
      double e = smearedTauMomentumVector.at(nobjT).energy() + smearedMuonMomentumVector.at(nobjE).energy() + TMath::Sqrt((theMETVector.px() * theMETVector.px()) + (theMETVector.py() * theMETVector.py()));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
      return MassRecoInformation;
    }
  } else {
    reco::Candidate::LorentzVector The_LorentzVect = smearedTauMomentumVector.at(nobjT) + smearedMuonMomentumVector.at(nobjE);
    pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
    return MassRecoInformation;
  }
}
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::CalculateThe4Momentum(const pat::Electron& patElectron, int nobjE, const pat::Muon& patMuon, int nobjM) {
  if(_UseVectorSumOfVisProductsAndMetMassReco) {
    double px = smearedElectronMomentumVector.at(nobjE).px() + smearedMuonMomentumVector.at(nobjM).px() + theMETVector.px();
    double py = smearedElectronMomentumVector.at(nobjE).py() + smearedMuonMomentumVector.at(nobjM).py() + theMETVector.py();
    double pz = smearedElectronMomentumVector.at(nobjE).pz() + smearedMuonMomentumVector.at(nobjM).pz();
    double e = smearedElectronMomentumVector.at(nobjE).energy() + smearedMuonMomentumVector.at(nobjM).energy() + TMath::Sqrt((theMETVector.px() * theMETVector.px()) + (theMETVector.py() * theMETVector.py()));
    reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
    pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
    return MassRecoInformation;
  } else if(_UseCollinerApproxMassReco) {
    double x1_numerator = (smearedElectronMomentumVector.at(nobjE).px() * smearedMuonMomentumVector.at(nobjM).py()) - (smearedMuonMomentumVector.at(nobjM).px() * smearedElectronMomentumVector.at(nobjE).py());
    double x1_denominator = (smearedMuonMomentumVector.at(nobjM).py() * (smearedElectronMomentumVector.at(nobjE).px() + theMETVector.px())) - (smearedMuonMomentumVector.at(nobjM).px() * (smearedElectronMomentumVector.at(nobjE).py() + theMETVector.py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (smearedElectronMomentumVector.at(nobjE).px() * (smearedMuonMomentumVector.at(nobjM).py() + theMETVector.py())) - (smearedElectronMomentumVector.at(nobjE).py() * (smearedMuonMomentumVector.at(nobjM).px() + theMETVector.px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    if ( (x1 > 0. && x1 < 1.) &&
	 (x2 > 0. && x2 < 1.) ) {
      reco::Candidate::LorentzVector The_LorentzVect = (smearedElectronMomentumVector.at(nobjE) / x1) + (smearedMuonMomentumVector.at(nobjM) / x2);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else {
      double px = smearedElectronMomentumVector.at(nobjE).px() + smearedMuonMomentumVector.at(nobjM).px() + theMETVector.px();
      double py = smearedElectronMomentumVector.at(nobjE).py() + smearedMuonMomentumVector.at(nobjM).py() + theMETVector.py();
      double pz = smearedElectronMomentumVector.at(nobjE).pz() + smearedMuonMomentumVector.at(nobjM).pz();
      double e = smearedElectronMomentumVector.at(nobjE).energy() + smearedMuonMomentumVector.at(nobjM).energy() + TMath::Sqrt((theMETVector.px() * theMETVector.px()) + (theMETVector.py() * theMETVector.py()));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
      return MassRecoInformation;
    }
  } else {
    reco::Candidate::LorentzVector The_LorentzVect = smearedElectronMomentumVector.at(nobjE) + smearedMuonMomentumVector.at(nobjM);
    pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
    return MassRecoInformation;
  }
}
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::CalculateThe4Momentum(const pat::Muon& patMuon1, int nobj1, const pat::Muon& patMuon2, int nobj2) {
  if(_UseVectorSumOfVisProductsAndMetMassReco) {
    double px = smearedMuonMomentumVector.at(nobj1).px() + smearedMuonMomentumVector.at(nobj2).px() + theMETVector.px();
    double py = smearedMuonMomentumVector.at(nobj1).py() + smearedMuonMomentumVector.at(nobj2).py() + theMETVector.py();
    double pz = smearedMuonMomentumVector.at(nobj1).pz() + smearedMuonMomentumVector.at(nobj2).pz();
    double e = smearedMuonMomentumVector.at(nobj1).energy() + smearedMuonMomentumVector.at(nobj2).energy() + TMath::Sqrt((theMETVector.px() * theMETVector.px()) + (theMETVector.py() * theMETVector.py()));
    reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
    pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
    return MassRecoInformation;
  } else if(_UseCollinerApproxMassReco) {
    double x1_numerator = (smearedMuonMomentumVector.at(nobj1).px() * smearedMuonMomentumVector.at(nobj2).py()) - (smearedMuonMomentumVector.at(nobj2).px() * smearedMuonMomentumVector.at(nobj1).py());
    double x1_denominator = (smearedMuonMomentumVector.at(nobj2).py() * (smearedMuonMomentumVector.at(nobj1).px() + theMETVector.px())) - (smearedMuonMomentumVector.at(nobj2).px() * (smearedMuonMomentumVector.at(nobj1).py() + theMETVector.py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (smearedMuonMomentumVector.at(nobj1).px() * (smearedMuonMomentumVector.at(nobj2).py() + theMETVector.py())) - (smearedMuonMomentumVector.at(nobj1).py() * (smearedMuonMomentumVector.at(nobj2).px() + theMETVector.px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    if ( (x1 > 0. && x1 < 1.) &&
	 (x2 > 0. && x2 < 1.) ) {
      reco::Candidate::LorentzVector The_LorentzVect = (smearedMuonMomentumVector.at(nobj1) / x1) + (smearedMuonMomentumVector.at(nobj2) / x2);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else {
      double px = smearedMuonMomentumVector.at(nobj1).px() + smearedMuonMomentumVector.at(nobj2).px() + theMETVector.px();
      double py = smearedMuonMomentumVector.at(nobj1).py() + smearedMuonMomentumVector.at(nobj2).py() + theMETVector.py();
      double pz = smearedMuonMomentumVector.at(nobj1).pz() + smearedMuonMomentumVector.at(nobj2).pz();
      double e = smearedMuonMomentumVector.at(nobj1).energy() + smearedMuonMomentumVector.at(nobj2).energy() + TMath::Sqrt((theMETVector.px() * theMETVector.px()) + (theMETVector.py() * theMETVector.py()));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
      return MassRecoInformation;
    }
  } else {
    reco::Candidate::LorentzVector The_LorentzVect = smearedMuonMomentumVector.at(nobj1) + smearedMuonMomentumVector.at(nobj2);
    pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
    return MassRecoInformation;
  }
}
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::CalculateThe4Momentum(const pat::Electron& patElectron1, int nobj1, const pat::Electron& patElectron2, int nobj2) {
  if(_UseVectorSumOfVisProductsAndMetMassReco) {
    double px = smearedElectronMomentumVector.at(nobj1).px() + smearedElectronMomentumVector.at(nobj2).px() + theMETVector.px();
    double py = smearedElectronMomentumVector.at(nobj1).py() + smearedElectronMomentumVector.at(nobj2).py() + theMETVector.py();
    double pz = smearedElectronMomentumVector.at(nobj1).pz() + smearedElectronMomentumVector.at(nobj2).pz();
    double e = smearedElectronMomentumVector.at(nobj1).energy() + smearedElectronMomentumVector.at(nobj2).energy() + TMath::Sqrt((theMETVector.px() * theMETVector.px()) + (theMETVector.py() * theMETVector.py()));
    reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
    pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
    return MassRecoInformation;
  } else if(_UseCollinerApproxMassReco) {
    double x1_numerator = (smearedElectronMomentumVector.at(nobj1).px() * smearedElectronMomentumVector.at(nobj2).py()) - (smearedElectronMomentumVector.at(nobj2).px() * smearedElectronMomentumVector.at(nobj1).py());
    double x1_denominator = (smearedElectronMomentumVector.at(nobj2).py() * (smearedElectronMomentumVector.at(nobj1).px() + theMETVector.px())) - (smearedElectronMomentumVector.at(nobj2).px() * (smearedElectronMomentumVector.at(nobj1).py() + theMETVector.py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (smearedElectronMomentumVector.at(nobj1).px() * (smearedElectronMomentumVector.at(nobj2).py() + theMETVector.py())) - (smearedElectronMomentumVector.at(nobj1).py() * (smearedElectronMomentumVector.at(nobj2).px() + theMETVector.px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    if ( (x1 > 0. && x1 < 1.) &&
	 (x2 > 0. && x2 < 1.) ) {
      reco::Candidate::LorentzVector The_LorentzVect = (smearedElectronMomentumVector.at(nobj1) / x1) + (smearedElectronMomentumVector.at(nobj2) / x2);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else {
      double px = smearedElectronMomentumVector.at(nobj1).px() + smearedElectronMomentumVector.at(nobj2).px() + theMETVector.px();
      double py = smearedElectronMomentumVector.at(nobj1).py() + smearedElectronMomentumVector.at(nobj2).py() + theMETVector.py();
      double pz = smearedElectronMomentumVector.at(nobj1).pz() + smearedElectronMomentumVector.at(nobj2).pz();
      double e = smearedElectronMomentumVector.at(nobj1).energy() + smearedElectronMomentumVector.at(nobj2).energy() + TMath::Sqrt((theMETVector.px() * theMETVector.px()) + (theMETVector.py() * theMETVector.py()));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
      return MassRecoInformation;
    }
  } else {
    reco::Candidate::LorentzVector The_LorentzVect = smearedElectronMomentumVector.at(nobj1) + smearedElectronMomentumVector.at(nobj2);
    pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
    return MassRecoInformation;
  }
}
pair<bool, reco::Candidate::LorentzVector> HiMassTauAnalysis::CalculateThe4Momentum(const pat::Tau& patTau1, int nobj1, const pat::Tau& patTau2, int nobj2) {
  if(_UseVectorSumOfVisProductsAndMetMassReco) {
    double px = smearedTauMomentumVector.at(nobj1).px() + smearedTauMomentumVector.at(nobj2).px() + theMETVector.px();
    double py = smearedTauMomentumVector.at(nobj1).py() + smearedTauMomentumVector.at(nobj2).py() + theMETVector.py();
    double pz = smearedTauMomentumVector.at(nobj1).pz() + smearedTauMomentumVector.at(nobj2).pz();
    double e = smearedTauMomentumVector.at(nobj1).energy() + smearedTauMomentumVector.at(nobj2).energy() + TMath::Sqrt((theMETVector.px() * theMETVector.px()) + (theMETVector.py() * theMETVector.py()));
    reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
    pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
    return MassRecoInformation;
  } else if(_UseCollinerApproxMassReco) {
    double x1_numerator = (smearedTauMomentumVector.at(nobj1).px() * smearedTauMomentumVector.at(nobj2).py()) - (smearedTauMomentumVector.at(nobj2).px() * smearedTauMomentumVector.at(nobj1).py());
    double x1_denominator = (smearedTauMomentumVector.at(nobj2).py() * (smearedTauMomentumVector.at(nobj1).px() + theMETVector.px())) - (smearedTauMomentumVector.at(nobj2).px() * (smearedTauMomentumVector.at(nobj1).py() + theMETVector.py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (smearedTauMomentumVector.at(nobj1).px() * (smearedTauMomentumVector.at(nobj2).py() + theMETVector.py())) - (smearedTauMomentumVector.at(nobj1).py() * (smearedTauMomentumVector.at(nobj2).px() + theMETVector.px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    if ( (x1 > 0. && x1 < 1.) &&
	 (x2 > 0. && x2 < 1.) ) {
      reco::Candidate::LorentzVector The_LorentzVect = (smearedTauMomentumVector.at(nobj1) / x1) + (smearedTauMomentumVector.at(nobj2) / x2);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
      return MassRecoInformation;
    } else {
      double px = smearedTauMomentumVector.at(nobj1).px() + smearedTauMomentumVector.at(nobj2).px() + theMETVector.px();
      double py = smearedTauMomentumVector.at(nobj1).py() + smearedTauMomentumVector.at(nobj2).py() + theMETVector.py();
      double pz = smearedTauMomentumVector.at(nobj1).pz() + smearedTauMomentumVector.at(nobj2).pz();
      double e = smearedTauMomentumVector.at(nobj1).energy() + smearedTauMomentumVector.at(nobj2).energy() + TMath::Sqrt((theMETVector.px() * theMETVector.px()) + (theMETVector.py() * theMETVector.py()));
      reco::Candidate::LorentzVector The_LorentzVect(px, py, pz, e);
      pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(false,The_LorentzVect);
      return MassRecoInformation;
    }
  } else {
    reco::Candidate::LorentzVector The_LorentzVect = smearedTauMomentumVector.at(nobj1) + smearedTauMomentumVector.at(nobj2);
    pair<bool, reco::Candidate::LorentzVector> MassRecoInformation(true,The_LorentzVect);
    return MassRecoInformation;
  }
}

//-----Calculate lepton+met transverse mass
double HiMassTauAnalysis::CalculateLeptonMetMt(const pat::Muon& patMuon, int nobj) {
  double px = smearedMuonMomentumVector.at(nobj).px() + theMETVector.px();
  double py = smearedMuonMomentumVector.at(nobj).py() + theMETVector.py();
  double et = smearedMuonMomentumVector.at(nobj).Et() + TMath::Sqrt((theMETVector.px() * theMETVector.px()) + (theMETVector.py() * theMETVector.py()));
  double mt2 = et*et - (px*px + py*py);
  if ( mt2 < 0 ) { return -1.; }
  else { return sqrt(mt2); }
}
double HiMassTauAnalysis::CalculateLeptonMetMt(const pat::Electron& patElectron, int nobj) {
  double px = smearedElectronMomentumVector.at(nobj).px() + theMETVector.px();
  double py = smearedElectronMomentumVector.at(nobj).py() + theMETVector.py();
  double et = smearedElectronMomentumVector.at(nobj).Et() + TMath::Sqrt((theMETVector.px() * theMETVector.px()) + (theMETVector.py() * theMETVector.py()));
  double mt2 = et*et - (px*px + py*py);
  if ( mt2 < 0 ) { return -1.; }
  else { return sqrt(mt2); }
}
double HiMassTauAnalysis::CalculateLeptonMetMt(const pat::Tau& patTau, int nobj) {
  double px = smearedTauMomentumVector.at(nobj).px() + theMETVector.px();
  double py = smearedTauMomentumVector.at(nobj).py() + theMETVector.py();
  double et = smearedTauMomentumVector.at(nobj).Et() + TMath::Sqrt((theMETVector.px() * theMETVector.px()) + (theMETVector.py() * theMETVector.py()));
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
  const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
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
        if((reco::deltaR((**iTrk).eta(),(**iTrk).phi(),patTau.leadPFChargedHadrCand()->eta(),patTau.leadPFChargedHadrCand()->phi())<_RecoTauIsoDeltaRCone) && ((**iTrk).pt()>_RecoTauTrackIsoTrkThreshold)
           && ((**iTrk).trackRef()->recHitsSize() > 8) && ((**iTrk).trackRef().isNonnull()) && (fabs((**iTrk).trackRef()->dxy(thePrimaryEventVertex.position())) < 0.03) 
           && (fabs((**iTrk).trackRef()->dz(thePrimaryEventVertex.position())) < 0.2) ) {
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
        if((reco::deltaR((**iTrk).eta(),(**iTrk).phi(),patTau.leadPFChargedHadrCand()->eta(),patTau.leadPFChargedHadrCand()->phi())< deltaRCone) &&((**iTrk).pt()> trkMinPt)
           && ((**iTrk).trackRef()->recHitsSize() > 8) ) {
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

//-----Calculate Mass from Tau signal track AND pi0 constituents
pair<reco::Candidate::LorentzVector,int> HiMassTauAnalysis::CalculateTauSignalTracksAndPiZerosMass(const pat::Tau& patTau) {
//reco::Candidate::LorentzVector HiMassTauAnalysis::CalculateTauSignalTracksAndPiZerosMass(const pat::Tau& patTau) {
//  double px=0;
//  double py=0;
//  double pz=0;
//  double e=0;
  int npi0 = 0;
  reco::Candidate::LorentzVector TheSignalTracksAndPiZeros_LorentzVect(0.,0.,0.,0.);
  if (patTau.isCaloTau()) {
  } else {
    if (patTau.leadPFChargedHadrCand().isNonnull()) {
      for ( reco::PFTauCollection::const_iterator hpsTau = _hpsTau->begin(); hpsTau != _hpsTau->end(); ++hpsTau ) {
        const std::vector<RecoTauPiZero> &signalpizeros = hpsTau->signalPiZeroCandidates();
        const std::vector<RecoTauPiZero> &isopizeros = hpsTau->isolationPiZeroCandidates();
        if((int)(signalpizeros.size()) > 0) {
          for(int nsigpi0s=0; nsigpi0s < (int)(signalpizeros.size()); nsigpi0s++) {
            const RecoTauPiZero &signalpizero = signalpizeros.at(nsigpi0s);
            if( reco::deltaR(patTau.leadPFChargedHadrCand()->p4(),signalpizero.p4()) < (5.0 / patTau.leadPFChargedHadrCand()->et()) ) {
//            std::cout << "tau #" << nsigpi0s+1 << " signal pi0 px  = " << signalpizero.px() << std::endl;
//            std::cout << "tau #" << nsigpi0s+1 << " signal pi0 energy = " << signalpizero.energy() << std::endl;
              TheSignalTracksAndPiZeros_LorentzVect += signalpizero.p4();
              npi0++;
            }
          }
        }
        if((int)(isopizeros.size()) > 0) {
          for(int nisopi0s=0; nisopi0s < (int)(isopizeros.size()); nisopi0s++) {
            const RecoTauPiZero &isopizero = isopizeros.at(nisopi0s);
            if( reco::deltaR(patTau.leadPFChargedHadrCand()->p4(),isopizero.p4()) < (5.0 / patTau.leadPFChargedHadrCand()->et()) ) {
//            std::cout << "tau #" << nisopi0s+1 << " signal pi0 px  = " << isopizero.px() << std::endl;
//            std::cout << "tau #" << nisopi0s+1 << " signal pi0 energy = " << isopizero.energy() << std::endl;
              TheSignalTracksAndPiZeros_LorentzVect += isopizero.p4();
              npi0++;
            }
          }
        }
      }
    }
  }
  pair<reco::Candidate::LorentzVector,int> TheSignalTracksAndPiZerosInfo(TheSignalTracksAndPiZeros_LorentzVect, npi0);
  return TheSignalTracksAndPiZerosInfo;
}

/*
//-----Calculate Lepton Isolation Quantities
pair<int, double> HiMassTauAnalysis::CalculateLeptonTrackIsolation(const pat::Muon& patMuon) {
  int nIsoTrks=0;
  double sumPtIsoTrks=0;

  for( unsigned ipfcand=0; ipfcand < _pflow->size(); ++ipfcand ) {
    const reco::PFCandidate& cand = (*_pflow)[ipfcand];
    if(cand.particleId()!=reco::PFCandidate::e)
  }

  for(PFCandidateRefVector::const_iterator iTrk=TauIsolTracks.begin();iTrk!=TauIsolTracks.end();++iTrk) {
    if((reco::deltaR((**iTrk).eta(),(**iTrk).phi(),patTau.leadPFChargedHadrCand()->eta(),patTau.leadPFChargedHadrCand()->phi())<_RecoTauIsoDeltaRCone) && ((**iTrk).pt()>_RecoTauTrackIsoTrkThreshold)) {
      nIsoTrks++;
      sumPtIsoTrks = sumPtIsoTrks + (**iTrk).pt();
    }
  }
  pair<int, double> LeptonIsoTrksInformation(nIsoTrks,sumPtIsoTrks);
  return LeptonIsoTrksInformation;
}
*/

//-----Smear the light leptons (for studies of systematic uncertanties)
pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> HiMassTauAnalysis::SmearLightLepton(const pat::Muon& patMuon) {
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
      reco::Candidate::LorentzVector smearedMomentum(smearedPtEtaPhiMVector.px(), smearedPtEtaPhiMVector.py(), smearedPtEtaPhiMVector.pz(), patMuon.energy());
      pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> theSmearedMomentumPair;
      theSmearedMomentumPair = make_pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector>(smearedMomentum, smearedPtEtaPhiMVector);
      return theSmearedMomentumPair;
    } else {
      math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector(patMuon.pt(), patMuon.eta(), patMuon.phi(), patMuon.mass());
      reco::Candidate::LorentzVector smearedMomentum = patMuon.p4();
      pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> theSmearedMomentumPair;
      theSmearedMomentumPair = make_pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector>(smearedMomentum, smearedPtEtaPhiMVector);
      return theSmearedMomentumPair;
    }
  } else {
    math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector(patMuon.pt(), patMuon.eta(), patMuon.phi(), patMuon.mass());
    reco::Candidate::LorentzVector smearedMomentum = patMuon.p4();
    pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> theSmearedMomentumPair;
    theSmearedMomentumPair = make_pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector>(smearedMomentum, smearedPtEtaPhiMVector);
    return theSmearedMomentumPair;
  }
}
pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> HiMassTauAnalysis::SmearLightLepton(const pat::Electron& patElectron) {
  double smearedPt;
  double smearedEt;
  double smearedEta;
  double smearedScEta;
  double smearedPhi;
  if(_GenParticleSource.label() != "") {
    if(matchToGen(patElectron).first) {
      reco::Candidate::LorentzVector unsmearedMomentum = matchToGen(patElectron).second; 
      if(_UseHeepInfo) {
        heep::Ele theHeepElec(patElectron);
        if(_SmearThePt) {
	  smearedEt = (unsmearedMomentum.energy() * sin(unsmearedMomentum.theta()) * _ElectronPtScaleOffset) + ((theHeepElec.et() -  (unsmearedMomentum.energy() * sin(unsmearedMomentum.theta()))) * _ElectronPtSigmaOffset);
	  smearedPt = (unsmearedMomentum.pt() * _ElectronPtScaleOffset) + ((patElectron.pt() -  unsmearedMomentum.pt()) * _ElectronPtSigmaOffset);
	} else {smearedEt  = theHeepElec.et(); smearedPt  = patElectron.pt();}
	if(_SmearTheEta) {
	  smearedScEta = (unsmearedMomentum.eta() * _ElectronEtaScaleOffset) + ((theHeepElec.scEta() -  unsmearedMomentum.eta()) * _ElectronEtaSigmaOffset);
	  smearedEta = (unsmearedMomentum.eta() * _ElectronEtaScaleOffset) + ((patElectron.eta() -  unsmearedMomentum.eta()) * _ElectronEtaSigmaOffset);
	} else {smearedScEta  = theHeepElec.scEta(); smearedEta  = patElectron.eta();}
	if(_SmearThePhi) {
	  smearedPhi = (unsmearedMomentum.phi() * _ElectronPhiScaleOffset) + ((patElectron.phi() -  unsmearedMomentum.phi()) * _ElectronPhiSigmaOffset);
	} else {smearedPhi  = patElectron.phi();}
	math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector(smearedPt, smearedEta, smearedPhi, unsmearedMomentum.mass());
	math::PtEtaPhiMLorentzVector smearedEtEtaPhiMVector(smearedEt, smearedScEta, smearedPhi, unsmearedMomentum.mass());
	reco::Candidate::LorentzVector smearedMomentum(smearedPtEtaPhiMVector.px(), smearedPtEtaPhiMVector.py(), smearedPtEtaPhiMVector.pz(), patElectron.energy());
	pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> theSmearedMomentumPair;
	theSmearedMomentumPair = make_pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector>(smearedMomentum, smearedEtEtaPhiMVector);
	return theSmearedMomentumPair;
      } else {
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
	reco::Candidate::LorentzVector smearedMomentum(smearedPtEtaPhiMVector.px(), smearedPtEtaPhiMVector.py(), smearedPtEtaPhiMVector.pz(), patElectron.energy());
	pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> theSmearedMomentumPair;
	theSmearedMomentumPair = make_pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector>(smearedMomentum, smearedPtEtaPhiMVector);
	return theSmearedMomentumPair;
      }
    } else {
      math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector(patElectron.pt(), patElectron.eta(), patElectron.phi(), patElectron.mass());
      reco::Candidate::LorentzVector smearedMomentum = patElectron.p4();
      pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> theSmearedMomentumPair;
      theSmearedMomentumPair = make_pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector>(smearedMomentum, smearedPtEtaPhiMVector);
      return theSmearedMomentumPair;
    }
  } else {
    math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector(patElectron.pt(), patElectron.eta(), patElectron.phi(), patElectron.mass());
    reco::Candidate::LorentzVector smearedMomentum = patElectron.p4();
    pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> theSmearedMomentumPair;
    theSmearedMomentumPair = make_pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector>(smearedMomentum, smearedPtEtaPhiMVector);
    return theSmearedMomentumPair;
  }
}
pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> HiMassTauAnalysis::SmearTau(const pat::Tau& patTau) {
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
      reco::Candidate::LorentzVector smearedMomentum(smearedPtEtaPhiMVector.px(), smearedPtEtaPhiMVector.py(), smearedPtEtaPhiMVector.pz(), patTau.energy());
      pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> theSmearedMomentumPair;
      theSmearedMomentumPair = make_pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector>(smearedMomentum, smearedPtEtaPhiMVector);
      return theSmearedMomentumPair;
    } else {
      math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector(patTau.pt(), patTau.eta(), patTau.phi(), patTau.mass());
      reco::Candidate::LorentzVector smearedMomentum = patTau.p4();
      pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> theSmearedMomentumPair;
      theSmearedMomentumPair = make_pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector>(smearedMomentum, smearedPtEtaPhiMVector);
      return theSmearedMomentumPair;
    }
  } else {
    math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector(patTau.pt(), patTau.eta(), patTau.phi(), patTau.mass());
    reco::Candidate::LorentzVector smearedMomentum = patTau.p4();
    pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> theSmearedMomentumPair;
    theSmearedMomentumPair = make_pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector>(smearedMomentum, smearedPtEtaPhiMVector);
    return theSmearedMomentumPair;
  }
}
pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> HiMassTauAnalysis::SmearJet(const pat::Jet& patJet) {
  bool isRealJet = true;
  reco::Candidate::LorentzVector tempJetVector;
  math::PtEtaPhiMLorentzVector tempPtEtaPhiMVector;
  if(_UseCorrectedJet) {
    tempJetVector = patJet.p4();
    math::PtEtaPhiMLorentzVector tempPtEtaPhiMVector(patJet.pt(), patJet.eta(), patJet.phi(), patJet.mass());
  } else {
    tempJetVector = patJet.correctedJet("raw","").p4();
    math::PtEtaPhiMLorentzVector tempPtEtaPhiMVector(patJet.correctedJet("raw","").pt(), patJet.correctedJet("raw","").eta(), patJet.correctedJet("raw","").phi(), patJet.correctedJet("raw","").mass());
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
      math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector(smearedMomentum.pt(), smearedMomentum.eta(), smearedMomentum.phi(), smearedMomentum.mass());
      pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> theSmearedMomentumPair;
      theSmearedMomentumPair = make_pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector>(smearedMomentum, smearedPtEtaPhiMVector);
      return theSmearedMomentumPair;
    } else {
      reco::Candidate::LorentzVector smearedMomentum = tempJetVector;
      math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector = tempPtEtaPhiMVector;
      pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> theSmearedMomentumPair;
      theSmearedMomentumPair = make_pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector>(smearedMomentum, smearedPtEtaPhiMVector);
      return theSmearedMomentumPair;
    }
  } else {
    reco::Candidate::LorentzVector smearedMomentum = tempJetVector;
    math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector = tempPtEtaPhiMVector;
    pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> theSmearedMomentumPair;
    theSmearedMomentumPair = make_pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector>(smearedMomentum, smearedPtEtaPhiMVector);
    return theSmearedMomentumPair;
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
      _hBestTauJetPt[NpdfCounter]                     = fs->make<TH1F>(("BestTauJetPt_"+j.str()).c_str(),     ("BestTauJetPt_"+j.str()).c_str(),     200, 0., 500.);
      _hBestTauJetEta[NpdfCounter]                    = fs->make<TH1F>(("BestTauJetEta_"+j.str()).c_str(),    ("BestTauJetEta_"+j.str()).c_str(),    72, -3.6, +3.6);
      _hTauJetPhi[NpdfCounter]                        = fs->make<TH1F>(("TauJetPhi_"+j.str()).c_str(),    ("TauJetPhi_"+j.str()).c_str(),    36, -TMath::Pi(), +TMath::Pi());
      _hTauJetNumSignalTracks[NpdfCounter]            = fs->make<TH1F>(("TauJetNumSignalTracks_"+j.str()).c_str(),           ("TauJetNumSignalTracks_"+j.str()).c_str(), 10, 0, 10);
      _hTauJetNumSignalGammas[NpdfCounter]            = fs->make<TH1F>(("TauJetNumSignalGammas_"+j.str()).c_str(),           ("TauJetNumSignalGammas_"+j.str()).c_str(), 10, 0, 10);
      _hTauJetSeedTrackPt[NpdfCounter]                = fs->make<TH1F>(("TauJetSeedTrackPt_"+j.str()).c_str(),               ("TauJetSeedTrackPt_"+j.str()).c_str(),     200, 0., 500.);
      _hBestTauJetSeedTrackPt[NpdfCounter]            = fs->make<TH1F>(("BestTauJetSeedTrackPt_"+j.str()).c_str(),           ("BestTauJetSeedTrackPt_"+j.str()).c_str(),     200, 0., 500.);
      _hTauJetSeedTrackIpSignificance[NpdfCounter]    = fs->make<TH1F>(("TauJetSeedTrackIpSignificance_"+j.str()).c_str(),   ("TauJetSeedTrackIpSignificance_"+j.str()).c_str(), 100, 0., 100.);
      _hTauJetSeedTrackNhits[NpdfCounter]             = fs->make<TH1F>(("TauJetSeedTrackNhits_"+j.str()).c_str(),            ("TauJetSeedTrackNhits_"+j.str()).c_str(), 40, 0., 40.);
      _hBestTauJetSeedTrackNhits[NpdfCounter]             = fs->make<TH1F>(("BestTauJetSeedTrackNhits_"+j.str()).c_str(),            ("BestTauJetSeedTrackNhits_"+j.str()).c_str(), 40, 0., 40.);
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
      _hBestTauJetSumPtIsoTracks[NpdfCounter]         = fs->make<TH1F>(("BestTauJetSumPtIsoTracks_"+j.str()).c_str(),        ("BestTauJetSumPtIsoTracks_"+j.str()).c_str(), 100, 0, 50);
      _hBestTauJetSumPtIsoGammas[NpdfCounter]         = fs->make<TH1F>(("BestTauJetSumPtIsoGammas_"+j.str()).c_str(),        ("BestTauJetSumPtIsoGammas_"+j.str()).c_str(), 100, 0, 50);
      _hTauJetSumPtIso[NpdfCounter]                   = fs->make<TH1F>(("TauJetSumPtIso_"+j.str()).c_str(),                  ("TauJetSumPtIso_"+j.str()).c_str(), 100, 0, 50);
      _hTauJetGenTauDeltaPhi[NpdfCounter]             = fs->make<TH1F>(("TauJetGenTauDeltaPhi_"+j.str()).c_str(),            ("TauJetGenTauDeltaPhi_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hTauJetGenTauDeltaEta[NpdfCounter]             = fs->make<TH1F>(("TauJetGenTauDeltaEta_"+j.str()).c_str(),            ("TauJetGenTauDeltaEta_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hTauJetGenTauDeltaPt[NpdfCounter]              = fs->make<TH1F>(("TauJetGenTauDeltaPt_"+j.str()).c_str(),             ("TauJetGenTauDeltaPt_"+j.str()).c_str(), 500, -5, 5);
      _hTauJetSignalTracksMass1prong[NpdfCounter]     = fs->make<TH1F>(("TauJetSignalTracksMass1prong_"+j.str()).c_str(),    ("TauJetSignalTracksMass1prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJetSignalTracksAndGammasMass1prong[NpdfCounter] = fs->make<TH1F>(("TauJetSignalTracksAndGammasMass1prong_"+j.str()).c_str(), ("TauJetSignalTracksAndGammasMass1prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJetSignalTracksAndPiZerosMass1prong[NpdfCounter] = fs->make<TH1F>(("TauJetSignalTracksAndPiZerosMass1prong_"+j.str()).c_str(), ("TauJetSignalTracksAndPiZerosMass1prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJetNumSignalPiZeros1prong[NpdfCounter] = fs->make<TH1F>(("TauJetNumSignalPiZeros1prong_"+j.str()).c_str(), ("TauJetNumSignalPiZeros1prong_"+j.str()).c_str(), 10, 0., 10.);
      _hTauJetSignalTracksMass3prong[NpdfCounter]          = fs->make<TH1F>(("TauJetSignalTracksMass3prong_"+j.str()).c_str(),          ("TauJetSignalTracksMass3prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJetSignalTracksAndGammasMass3prong[NpdfCounter] = fs->make<TH1F>(("TauJetSignalTracksAndGammasMass3prong_"+j.str()).c_str(), ("TauJetSignalTracksAndGammasMass3prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJetMass1Prong0Gamma[NpdfCounter]                = fs->make<TH1F>(("TauJetMass1Prong0Gamma_"+j.str()).c_str(),                ("TauJetMass1Prong0Gamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJetMass1Prong1Gamma[NpdfCounter]                = fs->make<TH1F>(("TauJetMass1Prong1Gamma_"+j.str()).c_str(),                ("TauJetMass1Prong1Gamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJetMass1Prong2orMoreGamma[NpdfCounter]          = fs->make<TH1F>(("TauJetMass1Prong2orMoreGamma_"+j.str()).c_str(),          ("TauJetMass1Prong2orMoreGamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJetMass3Prong0Gamma[NpdfCounter]                = fs->make<TH1F>(("TauJetMass3Prong0Gamma_"+j.str()).c_str(),                ("TauJetMass3Prong0Gamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJetMass3Prong1Gamma[NpdfCounter]                = fs->make<TH1F>(("TauJetMass3Prong1Gamma_"+j.str()).c_str(),                ("TauJetMass3Prong1Gamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJetMass3Prong2orMoreGamma[NpdfCounter]          = fs->make<TH1F>(("TauJetMass3Prong2orMoreGamma_"+j.str()).c_str(),          ("TauJetMass3Prong2orMoreGamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJetH3x3OverP[NpdfCounter]          = fs->make<TH1F>(("TauJetH3x3OverP_"+j.str()).c_str(),          ("TauJetH3x3OverP_"+j.str()).c_str(), 100, 0., 1.);
      _hBestTauJetH3x3OverP[NpdfCounter]          = fs->make<TH1F>(("BestTauJetH3x3OverP_"+j.str()).c_str(),          ("BestTauJetH3x3OverP_"+j.str()).c_str(), 100, 0., 1.);
    }
    
    //--- book reconstruction level histograms 
    if (_FillRecoMuonHists) {
      _hNMuon[NpdfCounter]                       = fs->make<TH1F>(("NMuon_"+j.str()).c_str(),                   ("NMuon_"+j.str()).c_str(), 21, 0., 20.);
      _hMuonEnergy[NpdfCounter]                  = fs->make<TH1F>(("MuonEnergy_"+j.str()).c_str(),              ("MuonEnergy_"+j.str()).c_str(), 200, 0., 500.);
      _hMuonPt[NpdfCounter]                      = fs->make<TH1F>(("MuonPt_"+j.str()).c_str(),                  ("MuonPt_"+j.str()).c_str(),  200, 0., 500.);
      _hMuonEta[NpdfCounter]                     = fs->make<TH1F>(("MuonEta_"+j.str()).c_str(),                 ("MuonEta_"+j.str()).c_str(), 72, -3.6, +3.6);
      _hBestMuonPt[NpdfCounter]                  = fs->make<TH1F>(("BestMuonPt_"+j.str()).c_str(),              ("BestMuonPt_"+j.str()).c_str(),  200, 0., 500.);
      _hBestMuonEta[NpdfCounter]                 = fs->make<TH1F>(("BestMuonEta_"+j.str()).c_str(),             ("BestMuonEta_"+j.str()).c_str(), 72, -3.6, +3.6);
      _hMuonPhi[NpdfCounter]                     = fs->make<TH1F>(("MuonPhi_"+j.str()).c_str(),                 ("MuonPhi_"+j.str()).c_str(), 36, -TMath::Pi(), +TMath::Pi());
      _hMuonTrackIso[NpdfCounter]                = fs->make<TH1F>(("MuonTrackIso_"+j.str()).c_str(),            ("MuonTrackIso_"+j.str()).c_str(), 100, 0, 50);
      _hMuonEcalIso[NpdfCounter]                 = fs->make<TH1F>(("MuonEcalIso_"+j.str()).c_str(),             ("MuonEcalIso_"+j.str()).c_str(), 100, 0, 50);
      _hBestMuonTrackIso[NpdfCounter]            = fs->make<TH1F>(("BestMuonTrackIso_"+j.str()).c_str(),        ("BestMuonTrackIso_"+j.str()).c_str(), 100, 0, 50);
      _hBestMuonEcalIso[NpdfCounter]             = fs->make<TH1F>(("BestMuonEcalIso_"+j.str()).c_str(),         ("BestMuonEcalIso_"+j.str()).c_str(), 100, 0, 50);
      _hMuonIso[NpdfCounter]                     = fs->make<TH1F>(("MuonIso_"+j.str()).c_str(),                 ("MuonIso_"+j.str()).c_str(), 100, 0, 50);
      _hMuonIp[NpdfCounter]                      = fs->make<TH1F>(("MuonIp_"+j.str()).c_str(),                  ("MuonIp_"+j.str()).c_str(), 500, -1, +1);
      _hBestMuonIp[NpdfCounter]                  = fs->make<TH1F>(("BestMuonIp_"+j.str()).c_str(),              ("BestMuonIp_"+j.str()).c_str(), 500, -1, +1);
      _hMuonIpSignificance[NpdfCounter]          = fs->make<TH1F>(("MuonIpSignificance_"+j.str()).c_str(),      ("MuonIpSignificance_"+j.str()).c_str(), 100, 0., 100.);
      _hMuonGenMuonDeltaPhi[NpdfCounter]         = fs->make<TH1F>(("MuonGenMuonDeltaPhi_"+j.str()).c_str(),     ("MuonGenMuonDeltaPhi_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hMuonGenMuonDeltaEta[NpdfCounter]         = fs->make<TH1F>(("MuonGenMuonDeltaEta_"+j.str()).c_str(),     ("MuonGenMuonDeltaEta_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hMuonGenMuonDeltaPt[NpdfCounter]          = fs->make<TH1F>(("MuonGenMuonDeltaPt_"+j.str()).c_str(),      ("MuonGenMuonDeltaPt_"+j.str()).c_str(), 500, -5, 5);
      _hMuonCaloCompatibility[NpdfCounter]       = fs->make<TH1F>(("MuonCaloCompatibility_"+j.str()).c_str(),   ("MuonCaloCompatibility_"+j.str()).c_str(), 102, 0.0, 1.02);
      _hMuonSegmentCompatibility[NpdfCounter]    = fs->make<TH1F>(("MuonSegmentCompatibility_"+j.str()).c_str(),("MuonSegmentCompatibility_"+j.str()).c_str(), 102, 0.0, 1.02);
      _hBestMuonCaloCompatibility[NpdfCounter]       = fs->make<TH1F>(("BestMuonCaloCompatibility_"+j.str()).c_str(),   ("BestMuonCaloCompatibility_"+j.str()).c_str(), 102, 0.0, 1.02);
      _hBestMuonSegmentCompatibility[NpdfCounter]    = fs->make<TH1F>(("BestMuonSegmentCompatibility_"+j.str()).c_str(),("BestMuonSegmentCompatibility_"+j.str()).c_str(), 102, 0.0, 1.02);
      _hMuonAntiPion[NpdfCounter]                = fs->make<TH1F>(("MuonAntiPion_"+j.str()).c_str(),            ("MuonAntiPion_"+j.str()).c_str(), 202, 0.0, 2.02);
      _hBestMuonAntiPion[NpdfCounter]            = fs->make<TH1F>(("BestMuonAntiPion_"+j.str()).c_str(),        ("BestMuonAntiPion_"+j.str()).c_str(), 202, 0.0, 2.02);
      _hMuonCaloCompatibilityVsSegmentCompatibility[NpdfCounter] = fs->make<TH2F>(("MuonCaloCompatibilityVsSegmentCompatibility_"+j.str()).c_str(), ("MuonCaloCompatibilityVsSegmentCompatibility_"+j.str()).c_str(), 102, 0, 1.02, 102, 0, 1.02);      
    }
    
    if (_FillRecoElectronHists) {
      _hNElectron[NpdfCounter]                   = fs->make<TH1F>(("NElectron_"+j.str()).c_str(),                  ("NElectron_"+j.str()).c_str(), 21, 0., 20.);
      _hElectronEnergy[NpdfCounter]              = fs->make<TH1F>(("ElectronEnergy_"+j.str()).c_str(),             ("ElectronEnergy_"+j.str()).c_str(), 200, 0., 500.);
      _hElectronPt[NpdfCounter]                  = fs->make<TH1F>(("ElectronPt_"+j.str()).c_str(),                 ("ElectronPt_"+j.str()).c_str(), 200, 0., 500.);
      _hElectronEta[NpdfCounter]                 = fs->make<TH1F>(("ElectronEta_"+j.str()).c_str(),                ("ElectronEta_"+j.str()).c_str(), 72, -3.6, +3.6);
      _hBestElectronPt[NpdfCounter]              = fs->make<TH1F>(("BestElectronPt_"+j.str()).c_str(),             ("BestElectronPt_"+j.str()).c_str(), 200, 0., 500.);
      _hBestElectronEta[NpdfCounter]             = fs->make<TH1F>(("BestElectronEta_"+j.str()).c_str(),            ("BestElectronEta_"+j.str()).c_str(), 72, -3.6, +3.6);
      _hElectronPhi[NpdfCounter]                 = fs->make<TH1F>(("ElectronPhi_"+j.str()).c_str(),                ("ElectronPhi_"+j.str()).c_str(), 36, -TMath::Pi(), +TMath::Pi());
      _hElectronTrackIso[NpdfCounter]            = fs->make<TH1F>(("ElectronTrackIso_"+j.str()).c_str(),           ("ElectronTrackIso_"+j.str()).c_str(), 100, 0, 50);
      _hElectronEcalIso[NpdfCounter]             = fs->make<TH1F>(("ElectronEcalIso_"+j.str()).c_str(),            ("ElectronEcalIso_"+j.str()).c_str(), 100, 0, 50);
      _hBestElectronTrackIso[NpdfCounter]        = fs->make<TH1F>(("BestElectronTrackIso_"+j.str()).c_str(),       ("BestElectronTrackIso_"+j.str()).c_str(), 100, 0, 50);
      _hBestElectronEcalIso[NpdfCounter]         = fs->make<TH1F>(("BestElectronEcalIso_"+j.str()).c_str(),        ("BestElectronEcalIso_"+j.str()).c_str(), 100, 0, 50);
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
      _hElectronHoverEm[NpdfCounter]             = fs->make<TH1F>(("ElectronHoverEm_"+j.str()).c_str(),            ("ElectronHoverEm_"+j.str()).c_str(), 100, 0., 0.5);
      _hElectronEESigmaIEtaIEta[NpdfCounter]     = fs->make<TH1F>(("ElectronEESigmaIEtaIEta_"+j.str()).c_str(),    ("ElectronEESigmaIEtaIEta_"+j.str()).c_str(), 100, 0., 0.05);
      _hElectronEEDEta[NpdfCounter]              = fs->make<TH1F>(("ElectronEEDEta_"+j.str()).c_str(),             ("ElectronEEDEta_"+j.str()).c_str(), 100, -0.05, 0.05);
      _hElectronEEDPhi[NpdfCounter]              = fs->make<TH1F>(("ElectronEEDPhi_"+j.str()).c_str(),             ("ElectronEEDPhi_"+j.str()).c_str(), 100, -0.2, 0.2);
      _hElectronEBSigmaIEtaIEta[NpdfCounter]     = fs->make<TH1F>(("ElectronEBSigmaIEtaIEta_"+j.str()).c_str(),    ("ElectronEBSigmaIEtaIEta_"+j.str()).c_str(), 100, 0., 0.05);
      _hElectronEBDEta[NpdfCounter]              = fs->make<TH1F>(("ElectronEBDEta_"+j.str()).c_str(),             ("ElectronEBDEta_"+j.str()).c_str(), 100, -0.05, 0.05);
      _hElectronEBDPhi[NpdfCounter]              = fs->make<TH1F>(("ElectronEBDPhi_"+j.str()).c_str(),             ("ElectronEBDPhi_"+j.str()).c_str(), 100, -0.2, 0.2);
      _hElectronEB2by5Over5by5[NpdfCounter]      = fs->make<TH1F>(("ElectronEB2by5Over5by5_"+j.str()).c_str(),     ("ElectronEB2by5Over5by5_"+j.str()).c_str(), 100, 0., 1.);
      _hElectronEB1by5Over5by5[NpdfCounter]      = fs->make<TH1F>(("ElectronEB1by5Over5by5_"+j.str()).c_str(),     ("ElectronEB1by5Over5by5_"+j.str()).c_str(), 100, 0., 1.);
      _hElectronMissingHits[NpdfCounter]         = fs->make<TH1F>(("ElectronMissingHits_"+j.str()).c_str(),        ("ElectronMissingHits_"+j.str()).c_str(), 10, 0., 10.);
      _hBestElectronHoverEm[NpdfCounter]             = fs->make<TH1F>(("BestElectronHoverEm_"+j.str()).c_str(),            ("BestElectronHoverEm_"+j.str()).c_str(), 100, 0., 0.5);
      _hBestElectronEESigmaIEtaIEta[NpdfCounter]     = fs->make<TH1F>(("BestElectronEESigmaIEtaIEta_"+j.str()).c_str(),    ("BestElectronEESigmaIEtaIEta_"+j.str()).c_str(), 100, 0., 0.05);
      _hBestElectronEEDEta[NpdfCounter]              = fs->make<TH1F>(("BestElectronEEDEta_"+j.str()).c_str(),             ("BestElectronEEDEta_"+j.str()).c_str(), 100, -0.05, 0.05);
      _hBestElectronEEDPhi[NpdfCounter]              = fs->make<TH1F>(("BestElectronEEDPhi_"+j.str()).c_str(),             ("BestElectronEEDPhi_"+j.str()).c_str(), 100, -0.2, 0.2);
      _hBestElectronEBSigmaIEtaIEta[NpdfCounter]     = fs->make<TH1F>(("BestElectronEBSigmaIEtaIEta_"+j.str()).c_str(),    ("BestElectronEBSigmaIEtaIEta_"+j.str()).c_str(), 100, 0., 0.05);
      _hBestElectronEBDEta[NpdfCounter]              = fs->make<TH1F>(("BestElectronEBDEta_"+j.str()).c_str(),             ("BestElectronEBDEta_"+j.str()).c_str(), 100, -0.05, 0.05);
      _hBestElectronEBDPhi[NpdfCounter]              = fs->make<TH1F>(("BestElectronEBDPhi_"+j.str()).c_str(),             ("BestElectronEBDPhi_"+j.str()).c_str(), 100, -0.2, 0.2);
      _hBestElectronEB2by5Over5by5[NpdfCounter]      = fs->make<TH1F>(("BestElectronEB2by5Over5by5_"+j.str()).c_str(),     ("BestElectronEB2by5Over5by5_"+j.str()).c_str(), 100, 0., 1.);
      _hBestElectronEB1by5Over5by5[NpdfCounter]      = fs->make<TH1F>(("BestElectronEB1by5Over5by5_"+j.str()).c_str(),     ("BestElectronEB1by5Over5by5_"+j.str()).c_str(), 100, 0., 1.);
      _hBestElectronMissingHits[NpdfCounter]         = fs->make<TH1F>(("BestElectronMissingHits_"+j.str()).c_str(),        ("BestElectronMissingHits_"+j.str()).c_str(), 10, 0., 10.);
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
      _hFirstLeadingJetPt[NpdfCounter]       = fs->make<TH1F>(("FirstLeadingJetPt_"+j.str()).c_str(),     ("FirstLeadingJetPt_"+j.str()).c_str(), 200, 0., 1000.);
      _hSecondLeadingJetPt[NpdfCounter]       = fs->make<TH1F>(("SecondLeadingJetPt_"+j.str()).c_str(),     ("SecondLeadingJetPt_"+j.str()).c_str(), 200, 0., 1000.);
      _hMHT[NpdfCounter]                    = fs->make<TH1F>(("MHT_"+j.str()).c_str(),                    ("MHT_"+j.str()).c_str(), 100, 0, 1000);
      _hHT[NpdfCounter]                    = fs->make<TH1F>(("HT_"+j.str()).c_str(),                    ("HT_"+j.str()).c_str(), 100, 0, 1000);
    }
    
    if (_FillTopologyHists) {
      if( ((_AnalyzeMuonForLeg1) && (_AnalyzeTauForLeg2)) || ((_AnalyzeMuonForLeg2) && (_AnalyzeTauForLeg1)) ) {
	_hMuonPtVsTauPt[NpdfCounter]                   = fs->make<TH2F>(("MuonPtVsTauPt_"+j.str()).c_str(),   ("MuonPtVsTauPt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hMuonTauDeltaR[NpdfCounter]                   = fs->make<TH1F>(("MuonTauDeltaR_"+j.str()).c_str(),   ("MuonTauDeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hBestMuonTauDeltaR[NpdfCounter]               = fs->make<TH1F>(("BestMuonTauDeltaR_"+j.str()).c_str(),   ("BestMuonTauDeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hMuonTauDeltaPtDivSumPt[NpdfCounter]          = fs->make<TH1F>(("MuonTauDeltaPtDivSumPt_"+j.str()).c_str(), ("MuonTauDeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hMuonTauDeltaPt[NpdfCounter]                  = fs->make<TH1F>(("MuonTauDeltaPt_"+j.str()).c_str(), ("MuonTauDeltaPt_"+j.str()).c_str(), 100, 0, 1000);
	_hMuonMetMt[NpdfCounter]                       = fs->make<TH1F>(("MuonMetMt_"+j.str()).c_str(),       ("MuonMetMt_"+j.str()).c_str(), 100, 0, 500);
	_hBestMuonMetMt[NpdfCounter]                   = fs->make<TH1F>(("BestMuonMetMt_"+j.str()).c_str(),       ("BestMuonMetMt_"+j.str()).c_str(), 100, 0, 500);
	_hTauMetMt[NpdfCounter]                        = fs->make<TH1F>(("TauMetMt_"+j.str()).c_str(),        ("TauMetMt_"+j.str()).c_str(), 100, 0, 500);
	_hMuonTauOSLS[NpdfCounter]                     = fs->make<TH1F>(("MuonTauOSLS_"+j.str()).c_str(),     ("MuonTauOSLS_"+j.str()).c_str(), 20, -10, 10);
	_hBestMuonTauOSLS[NpdfCounter]                 = fs->make<TH1F>(("BestMuonTauOSLS_"+j.str()).c_str(),     ("BestMuonTauOSLS_"+j.str()).c_str(), 20, -10, 10);
	_hMuonTauCosDphi[NpdfCounter]                  = fs->make<TH1F>(("MuonTauCosDphi_"+j.str()).c_str(),  ("MuonTauCosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hBestMuonTauCosDphi[NpdfCounter]                  = fs->make<TH1F>(("BestMuonTauCosDphi_"+j.str()).c_str(),  ("BestMuonTauCosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hMuonMetDeltaPhi[NpdfCounter]                 = fs->make<TH1F>(("MuonMetDeltaPhi_"+j.str()).c_str(), ("MuonMetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hTauMetDeltaPhi[NpdfCounter]                  = fs->make<TH1F>(("TauMetDeltaPhi_"+j.str()).c_str(), ("TauMetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hMuonMetDeltaPhiVsMuonTauCosDphi[NpdfCounter] = fs->make<TH2F>(("MuonMetDeltaPhiVsMuonTauCosDphi_"+j.str()).c_str(), ("MuonMetDeltaPhiVsMuonTauCosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
      }
      if( ((_AnalyzeElectronForLeg1) && (_AnalyzeTauForLeg2)) || ((_AnalyzeElectronForLeg2) && (_AnalyzeTauForLeg1)) ) {
	_hElectronPtVsTauPt[NpdfCounter]          = fs->make<TH2F>(("ElectronPtVsTauPt_"+j.str()).c_str(),          ("ElectronPtVsTauPt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hElectronTauDeltaR[NpdfCounter]          = fs->make<TH1F>(("ElectronTauDeltaR_"+j.str()).c_str(),          ("ElectronTauDeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hBestElectronTauDeltaR[NpdfCounter]      = fs->make<TH1F>(("BestElectronTauDeltaR_"+j.str()).c_str(),          ("BestElectronTauDeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hElectronTauDeltaPtDivSumPt[NpdfCounter] = fs->make<TH1F>(("ElectronTauDeltaPtDivSumPt_"+j.str()).c_str(), ("ElectronTauDeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hElectronTauDeltaPt[NpdfCounter]         = fs->make<TH1F>(("ElectronTauDeltaPt_"+j.str()).c_str(), ("ElectronTauDeltaPt_"+j.str()).c_str(), 100, 0, 1000);
	_hElectronMetMt[NpdfCounter]              = fs->make<TH1F>(("ElectronMetMt_"+j.str()).c_str(),              ("ElectronMetMt_"+j.str()).c_str(), 100, 0, 500);
	_hBestElectronMetMt[NpdfCounter]          = fs->make<TH1F>(("BestElectronMetMt_"+j.str()).c_str(),              ("BestElectronMetMt_"+j.str()).c_str(), 100, 0, 500);
	_hTauMetMt[NpdfCounter]                   = fs->make<TH1F>(("TauMetMt_"+j.str()).c_str(),                   ("TauMetMt_"+j.str()).c_str(), 100, 0, 500);
	_hElectronTauOSLS[NpdfCounter]            = fs->make<TH1F>(("ElectronTauOSLS_"+j.str()).c_str(),            ("ElectronTauOSLS_"+j.str()).c_str(), 20, -10, 10);
	_hBestElectronTauOSLS[NpdfCounter]        = fs->make<TH1F>(("BestElectronTauOSLS_"+j.str()).c_str(),            ("BestElectronTauOSLS_"+j.str()).c_str(), 20, -10, 10);
	_hElectronTauCosDphi[NpdfCounter]         = fs->make<TH1F>(("ElectronTauCosDphi_"+j.str()).c_str(),         ("ElectronTauCosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hBestElectronTauCosDphi[NpdfCounter]         = fs->make<TH1F>(("BestElectronTauCosDphi_"+j.str()).c_str(),         ("BestElectronTauCosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hElectronMetDeltaPhi[NpdfCounter]        = fs->make<TH1F>(("ElectronMetDeltaPhi_"+j.str()).c_str(),        ("ElectronMetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hTauMetDeltaPhi[NpdfCounter]             = fs->make<TH1F>(("TauMetDeltaPhi_"+j.str()).c_str(),             ("TauMetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hElectronMetDeltaPhiVsElectronTauCosDphi[NpdfCounter] = fs->make<TH2F>(("ElectronMetDeltaPhiVsElectronTauCosDphi_"+j.str()).c_str(), ("ElectronMetDeltaPhiVsElectronTauCosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
      }
      if( ((_AnalyzeMuonForLeg1) && (_AnalyzeElectronForLeg2)) || ((_AnalyzeMuonForLeg2) && (_AnalyzeElectronForLeg1)) ) {
	_hElectronPtVsMuonPt[NpdfCounter]          = fs->make<TH2F>(("ElectronPtVsMuonPt_"+j.str()).c_str(),          ("ElectronPtVsMuonPt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hElectronMuonDeltaR[NpdfCounter]          = fs->make<TH1F>(("ElectronMuonDeltaR_"+j.str()).c_str(),          ("ElectronMuonDeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hElectronMuonDeltaPtDivSumPt[NpdfCounter] = fs->make<TH1F>(("ElectronMuonDeltaPtDivSumPt_"+j.str()).c_str(), ("ElectronMuonDeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hElectronMuonDeltaPt[NpdfCounter]         = fs->make<TH1F>(("ElectronMuonDeltaPt_"+j.str()).c_str(), ("ElectronMuonDeltaPt_"+j.str()).c_str(), 100, 0, 1000);
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
	_hTau1Tau2DeltaPt[NpdfCounter]         = fs->make<TH1F>(("Tau1Tau2DeltaPt_"+j.str()).c_str(), ("Tau1Tau2DeltaPt_"+j.str()).c_str(), 100, 0, 1000);
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
	_hMuon1Muon2DeltaPt[NpdfCounter]         = fs->make<TH1F>(("Muon1Muon2DeltaPt_"+j.str()).c_str(), ("Muon1Muon2DeltaPt_"+j.str()).c_str(), 100, 0, 1000);
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
	_hElectron1Electron2DeltaPt[NpdfCounter]         = fs->make<TH1F>(("Electron1Electron2DeltaPt_"+j.str()).c_str(), ("Electron1Electron2DeltaPt_"+j.str()).c_str(), 100, 0, 1000);
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
      _hNotReconstructableMassOS[NpdfCounter] = fs->make<TH1F>(("NotReconstructableMassOS_"+j.str()).c_str(), ("NotReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hReconstructableMassOS[NpdfCounter]    = fs->make<TH1F>(("ReconstructableMassOS_"+j.str()).c_str(),    ("ReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hNotReconstructableMassLS[NpdfCounter] = fs->make<TH1F>(("NotReconstructableMassLS_"+j.str()).c_str(), ("NotReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hReconstructableMassLS[NpdfCounter]    = fs->make<TH1F>(("ReconstructableMassLS_"+j.str()).c_str(),    ("ReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hPZeta[NpdfCounter]                  = fs->make<TH1F>(("PZeta_"+j.str()).c_str(),                  ("PZeta_"+j.str()).c_str(), 200, -100, 100);
      _hPZetaVis[NpdfCounter]               = fs->make<TH1F>(("PZetaVis_"+j.str()).c_str(),               ("PZetaVis_"+j.str()).c_str(), 100, 0, 100);
      _hZeta2D[NpdfCounter]                 = fs->make<TH2F>(("Zeta2D_"+j.str()).c_str(),                 ("Zeta2D_"+j.str()).c_str(), 100, 0, 100, 200, -100, 100);
      _hZeta1D[NpdfCounter]                 = fs->make<TH1F>(("Zeta1D_"+j.str()).c_str(),                 ("Zeta1D_"+j.str()).c_str(), 50, -100, 100);
      _hBestZeta1D[NpdfCounter]             = fs->make<TH1F>(("BestZeta1D_"+j.str()).c_str(),             ("BestZeta1D_"+j.str()).c_str(), 50, -100, 100);
      _hMet[NpdfCounter]                    = fs->make<TH1F>(("Met_"+j.str()).c_str(),                    ("Met_"+j.str()).c_str(), 100, 0, 1000);
      _hR1[NpdfCounter] 		    = fs->make<TH1F>(("R1_"+j.str()).c_str(), ("R1_"+j.str()).c_str(), 60, 0, 6);
      _hR2[NpdfCounter] = fs->make<TH1F>(("R2_"+j.str()).c_str(), ("R2_"+j.str()).c_str(), 60, 0, 6);
      _hDphi1[NpdfCounter] = fs->make<TH1F>(("Dphi1_"+j.str()).c_str(), ("Dphi1_"+j.str()).c_str(), 36, -TMath::Pi(), +TMath::Pi());
      _hDphi2[NpdfCounter] = fs->make<TH1F>(("Dphi2_"+j.str()).c_str(), ("Dphi2_"+j.str()).c_str(), 36, -TMath::Pi(), +TMath::Pi());
      _hAlpha[NpdfCounter] = fs->make<TH1F>(("Alpha_"+j.str()).c_str(), ("Alpha_"+j.str()).c_str(), 50, 0, 2);
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

// Authors: Andres Florez, Alfredo Gurrola, Eduardo Luiggi, Chi Nhan Nguyen

#include "HighMassAnalysis/Analysis/interface/HiMassTauAnalysis.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPEle.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPCutCodes.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "boost/lexical_cast.hpp"
#include "boost/foreach.hpp"

#include <TMath.h>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace edm;
using namespace reco;

// constructors and destructor
//HiMassTauAnalysis::HiMassTauAnalysis(const ParameterSet& iConfig) {
HiMassTauAnalysis::HiMassTauAnalysis(const ParameterSet& iConfig) :
  Prefix  (iConfig.getParameter<std::string>  ("Prefix"  )),
  Suffix  (iConfig.getParameter<std::string>  ("Suffix"  )),
  scanFormat(iConfig.getParameter<std::string>("ScanFormat")),
  scanPars(iConfig.getParameter<std::vector<std::string> >("ScanParameters"))
{

  BOOST_FOREACH(const std::string& par, scanPars)
  //-----Generator level Inputs 
  _DoMSUGRApoint = iConfig.getParameter<bool>("DoMSUGRApoint");
  _GenParticleSource = iConfig.getUntrackedParameter<InputTag>("GenParticleSource");
  _SelectSusyScanPoint = iConfig.getParameter<bool>("SelectSusyScanPoint");
  _M0 = iConfig.getParameter<double>("M0");
  _M12 = iConfig.getParameter<double>("M12");

  _DoSMpoint = iConfig.getParameter<bool>("DoSMpoint");
  _mLSP = iConfig.getParameter<double>("mLSP");
  _mGL = iConfig.getParameter<double>("mGL");

  //-----Fill Histograms?
  _FillRecoVertexHists = true;
  _FillGenTauHists = true;
  _FillRecoTauHists = true;
  _FillRecoMuonHists = true;
  _FillRecoElectronHists = true;
  _FillRecoJetHists = true;
  _FillTopologyHists = true;

  //-----Reco Tau Inputs 
  _RecoTauSource = iConfig.getParameter<InputTag>("RecoTauSource");
  _RecoTau1EtaCut = iConfig.getParameter<double>("RecoTau1EtaCut");
  _RecoTau1PtMinCut = iConfig.getParameter<double>("RecoTau1PtMinCut");
  _RecoTau1PtMaxCut = iConfig.getParameter<double>("RecoTau1PtMaxCut");
  _DoRecoTau1DiscrByLeadTrack = iConfig.getParameter<bool>("DoRecoTau1DiscrByLeadTrack");
  _UseRecoTau1DiscrByLeadTrackFlag = iConfig.getParameter<bool>("UseRecoTau1DiscrByLeadTrackFlag");
  _RecoTau1DiscrByLeadTrack = iConfig.getUntrackedParameter<string>("RecoTau1DiscrByLeadTrack");
  _DoRecoTau1DiscrByLeadTrackNhits = iConfig.getParameter<bool>("DoRecoTau1DiscrByLeadTrackNhits");
  _RecoTau1LeadTrackMinHits = iConfig.getParameter<int>("RecoTau1LeadTrackMinHits");
  _DoRecoTau1DiscrByH3x3OverP = iConfig.getParameter<bool>("DoRecoTau1DiscrByH3x3OverP");
  _RecoTau1H3x3OverP = iConfig.getParameter<double>("RecoTau1H3x3OverP");
  _DoRecoTau1DiscrByIsolation = iConfig.getParameter<bool>("DoRecoTau1DiscrByIsolation");
  _UseRecoTau1DiscrByIsolationFlag = iConfig.getParameter<bool>("UseRecoTau1DiscrByIsolationFlag");
  _UseRecoTau1IsoSumPtInsteadOfNiso = iConfig.getParameter<bool>("UseRecoTau1IsoSumPtInsteadOfNiso");
  _UseRecoTau1EllipseForEcalIso = iConfig.getParameter<bool>("UseRecoTau1EllipseForEcalIso");
  _RecoTau1EcalIsoRphiForEllipse = iConfig.getParameter<double>("RecoTau1EcalIsoRphiForEllipse");
  _RecoTau1EcalIsoRetaForEllipse = iConfig.getParameter<double>("RecoTau1EcalIsoRetaForEllipse");
  _RecoTau1TrackNisoMax = iConfig.getParameter<int>("RecoTau1TrackNisoMax");
  _RecoTau1EcalNisoMax = iConfig.getParameter<int>("RecoTau1EcalNisoMax");
  _RecoTau1TrackIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoTau1TrackIsoSumPtMinCutValue");
  _RecoTau1TrackIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoTau1TrackIsoSumPtMaxCutValue");
  _RecoTau1EcalIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoTau1EcalIsoSumPtMinCutValue");
  _RecoTau1EcalIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoTau1EcalIsoSumPtMaxCutValue");
  _RecoTau1DiscrByIsolation = iConfig.getUntrackedParameter<string>("RecoTau1DiscrByIsolation");
  _RecoTau1DiscrByProngType = iConfig.getParameter<string>("RecoTau1DiscrByProngType");
  _RecoTau1LeadTrackThreshold = iConfig.getParameter<double>("RecoTau1LeadTrackThreshold");
  _RecoTau1SigGamThreshold = iConfig.getParameter<double>("RecoTau1SigGamThreshold");
  _RecoTau1IsoDeltaRCone = iConfig.getParameter<double>("RecoTau1IsoDeltaRCone");
  _RecoTau1TrackIsoTrkThreshold = iConfig.getParameter<double>("RecoTau1TrackIsoTrkThreshold");
  _RecoTau1GammaIsoGamThreshold = iConfig.getParameter<double>("RecoTau1GammaIsoGamThreshold");
  _DoRecoTau1DiscrAgainstElectron = iConfig.getParameter<bool>("DoRecoTau1DiscrAgainstElectron");
  _RecoTau1DiscrAgainstElectron = iConfig.getUntrackedParameter<string>("RecoTau1DiscrAgainstElectron");
  _DoRecoTau1DiscrByCrackCut = iConfig.getParameter<bool>("DoRecoTau1DiscrByCrackCut");
  _DoRecoTau1DiscrAgainstMuon = iConfig.getParameter<bool>("DoRecoTau1DiscrAgainstMuon");
  _RecoTau1DiscrAgainstMuon = iConfig.getUntrackedParameter<string>("RecoTau1DiscrAgainstMuon");
  _SelectTau1sThatAreMuons = iConfig.getParameter<bool>("SelectTau1sThatAreMuons");
  _SelectTau1sThatAreElectrons = iConfig.getParameter<bool>("SelectTau1sThatAreElectrons");

  _RecoTau2EtaCut = iConfig.getParameter<double>("RecoTau2EtaCut");
  _RecoTau2PtMinCut = iConfig.getParameter<double>("RecoTau2PtMinCut");
  _RecoTau2PtMaxCut = iConfig.getParameter<double>("RecoTau2PtMaxCut");
  _DoRecoTau2DiscrByLeadTrack = iConfig.getParameter<bool>("DoRecoTau2DiscrByLeadTrack");
  _UseRecoTau2DiscrByLeadTrackFlag = iConfig.getParameter<bool>("UseRecoTau2DiscrByLeadTrackFlag");
  _RecoTau2DiscrByLeadTrack = iConfig.getUntrackedParameter<string>("RecoTau2DiscrByLeadTrack");
  _DoRecoTau2DiscrByLeadTrackNhits = iConfig.getParameter<bool>("DoRecoTau2DiscrByLeadTrackNhits");
  _RecoTau2LeadTrackMinHits = iConfig.getParameter<int>("RecoTau2LeadTrackMinHits");
  _DoRecoTau2DiscrByH3x3OverP = iConfig.getParameter<bool>("DoRecoTau2DiscrByH3x3OverP");
  _RecoTau2H3x3OverP = iConfig.getParameter<double>("RecoTau2H3x3OverP");
  _DoRecoTau2DiscrByIsolation = iConfig.getParameter<bool>("DoRecoTau2DiscrByIsolation");
  _UseRecoTau2DiscrByIsolationFlag = iConfig.getParameter<bool>("UseRecoTau2DiscrByIsolationFlag");
  _UseRecoTau2IsoSumPtInsteadOfNiso = iConfig.getParameter<bool>("UseRecoTau2IsoSumPtInsteadOfNiso");
  _UseRecoTau2EllipseForEcalIso = iConfig.getParameter<bool>("UseRecoTau2EllipseForEcalIso");
  _RecoTau2EcalIsoRphiForEllipse = iConfig.getParameter<double>("RecoTau2EcalIsoRphiForEllipse");
  _RecoTau2EcalIsoRetaForEllipse = iConfig.getParameter<double>("RecoTau2EcalIsoRetaForEllipse");
  _RecoTau2TrackNisoMax = iConfig.getParameter<int>("RecoTau2TrackNisoMax");
  _RecoTau2EcalNisoMax = iConfig.getParameter<int>("RecoTau2EcalNisoMax");
  _RecoTau2TrackIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoTau2TrackIsoSumPtMinCutValue");
  _RecoTau2TrackIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoTau2TrackIsoSumPtMaxCutValue");
  _RecoTau2EcalIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoTau2EcalIsoSumPtMinCutValue");
  _RecoTau2EcalIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoTau2EcalIsoSumPtMaxCutValue");
  _RecoTau2DiscrByIsolation = iConfig.getUntrackedParameter<string>("RecoTau2DiscrByIsolation");
  _RecoTau2DiscrByProngType = iConfig.getParameter<string>("RecoTau2DiscrByProngType");
  _RecoTau2LeadTrackThreshold = iConfig.getParameter<double>("RecoTau2LeadTrackThreshold");
  _RecoTau2SigGamThreshold = iConfig.getParameter<double>("RecoTau2SigGamThreshold");
  _RecoTau2IsoDeltaRCone = iConfig.getParameter<double>("RecoTau2IsoDeltaRCone");
  _RecoTau2TrackIsoTrkThreshold = iConfig.getParameter<double>("RecoTau2TrackIsoTrkThreshold");
  _RecoTau2GammaIsoGamThreshold = iConfig.getParameter<double>("RecoTau2GammaIsoGamThreshold");
  _DoRecoTau2DiscrAgainstElectron = iConfig.getParameter<bool>("DoRecoTau2DiscrAgainstElectron");
  _RecoTau2DiscrAgainstElectron = iConfig.getUntrackedParameter<string>("RecoTau2DiscrAgainstElectron");
  _DoRecoTau2DiscrByCrackCut = iConfig.getParameter<bool>("DoRecoTau2DiscrByCrackCut");
  _DoRecoTau2DiscrAgainstMuon = iConfig.getParameter<bool>("DoRecoTau2DiscrAgainstMuon");
  _RecoTau2DiscrAgainstMuon = iConfig.getUntrackedParameter<string>("RecoTau2DiscrAgainstMuon");
  _SelectTau2sThatAreMuons = iConfig.getParameter<bool>("SelectTau2sThatAreMuons");
  _SelectTau2sThatAreElectrons = iConfig.getParameter<bool>("SelectTau2sThatAreElectrons");

  //-----Reco Muon Inputs
  _RecoMuonSource = iConfig.getParameter<InputTag>("RecoMuonSource");
  _RecoMuon1EtaCut = iConfig.getParameter<double>("RecoMuon1EtaCut");
  _RecoMuon1PtMinCut = iConfig.getParameter<double>("RecoMuon1PtMinCut");
  _RecoMuon1PtMaxCut = iConfig.getParameter<double>("RecoMuon1PtMaxCut");
  _DoRecoMuon1DiscrByGlobal = iConfig.getParameter<bool>("DoRecoMuon1DiscrByGlobal");
  _DoRecoMuon1DiscrByIsolation = iConfig.getParameter<bool>("DoRecoMuon1DiscrByIsolation");
  _RecoMuon1TrackIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoMuon1TrackIsoSumPtMinCutValue");
  _RecoMuon1TrackIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoMuon1TrackIsoSumPtMaxCutValue");
  _RecoMuon1EcalIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoMuon1EcalIsoSumPtMinCutValue");
  _RecoMuon1EcalIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoMuon1EcalIsoSumPtMaxCutValue");
  _RecoMuon1IsoDeltaRCone = iConfig.getParameter<double>("RecoMuon1IsoDeltaRCone");
  _RecoMuon1TrackIsoTrkThreshold = iConfig.getParameter<double>("RecoMuon1TrackIsoTrkThreshold");
  _RecoMuon1EcalIsoRecHitThreshold = iConfig.getParameter<double>("RecoMuon1EcalIsoRecHitThreshold");
  _DoRecoMuon1DiscrByIp = iConfig.getParameter<bool>("DoRecoMuon1DiscrByIp");
  _RecoMuon1IpCut = iConfig.getParameter<double>("RecoMuon1IpCut");
  _DoRecoMuon1DiscrByPionVeto = iConfig.getParameter<bool>("DoRecoMuon1DiscrByPionVeto");
  _RecoMuon1CaloCompCoefficient = iConfig.getParameter<double>("RecoMuon1CaloCompCoefficient");
  _RecoMuon1SegmCompCoefficient = iConfig.getParameter<double>("RecoMuon1SegmCompCoefficient");
  _RecoMuon1AntiPionCut = iConfig.getParameter<double>("RecoMuon1AntiPionCut");
  _RecoMuon2EtaCut = iConfig.getParameter<double>("RecoMuon2EtaCut");
  _RecoMuon2PtMinCut = iConfig.getParameter<double>("RecoMuon2PtMinCut");
  _RecoMuon2PtMaxCut = iConfig.getParameter<double>("RecoMuon2PtMaxCut");
  _DoRecoMuon2DiscrByGlobal = iConfig.getParameter<bool>("DoRecoMuon2DiscrByGlobal");
  _DoRecoMuon2DiscrByIsolation = iConfig.getParameter<bool>("DoRecoMuon2DiscrByIsolation");
  _RecoMuon2TrackIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoMuon2TrackIsoSumPtMinCutValue");
  _RecoMuon2TrackIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoMuon2TrackIsoSumPtMaxCutValue");
  _RecoMuon2EcalIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoMuon2EcalIsoSumPtMinCutValue");
  _RecoMuon2EcalIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoMuon2EcalIsoSumPtMaxCutValue");
  _RecoMuon2IsoDeltaRCone = iConfig.getParameter<double>("RecoMuon2IsoDeltaRCone");
  _RecoMuon2TrackIsoTrkThreshold = iConfig.getParameter<double>("RecoMuon2TrackIsoTrkThreshold");
  _RecoMuon2EcalIsoRecHitThreshold = iConfig.getParameter<double>("RecoMuon2EcalIsoRecHitThreshold");
  _DoRecoMuon2DiscrByIp = iConfig.getParameter<bool>("DoRecoMuon2DiscrByIp");
  _RecoMuon2IpCut = iConfig.getParameter<double>("RecoMuon2IpCut");
  _DoRecoMuon2DiscrByPionVeto = iConfig.getParameter<bool>("DoRecoMuon2DiscrByPionVeto");
  _RecoMuon2CaloCompCoefficient = iConfig.getParameter<double>("RecoMuon2CaloCompCoefficient");
  _RecoMuon2SegmCompCoefficient = iConfig.getParameter<double>("RecoMuon2SegmCompCoefficient");
  _RecoMuon2AntiPionCut = iConfig.getParameter<double>("RecoMuon2AntiPionCut");

  //-----Reco Electron Inputs
  _RecoElectronSource = iConfig.getParameter<InputTag>("RecoElectronSource");
  _UseHeepInfo = iConfig.getParameter<bool>("UseHeepInfo");
  _RecoElectron1EtaCut = iConfig.getParameter<double>("RecoElectron1EtaCut");
  _RecoElectron1PtMinCut = iConfig.getParameter<double>("RecoElectron1PtMinCut");
  _RecoElectron1PtMaxCut = iConfig.getParameter<double>("RecoElectron1PtMaxCut");
  _DoRecoElectron1DiscrByTrackIsolation = iConfig.getParameter<bool>("DoRecoElectron1DiscrByTrackIsolation");
  _RecoElectron1TrackIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoElectron1TrackIsoSumPtMaxCutValue");
  _RecoElectron1TrackIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoElectron1TrackIsoSumPtMinCutValue");
  _RecoElectron1TrackIsoDeltaRCone = iConfig.getParameter<double>("RecoElectron1TrackIsoDeltaRCone");
  _RecoElectron1TrackIsoTrkThreshold = iConfig.getParameter<double>("RecoElectron1TrackIsoTrkThreshold");
  _DoRecoElectron1DiscrByEcalIsolation = iConfig.getParameter<bool>("DoRecoElectron1DiscrByEcalIsolation");
  _RecoElectron1EcalIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoElectron1EcalIsoSumPtMaxCutValue");
  _RecoElectron1EcalIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoElectron1EcalIsoSumPtMinCutValue");
  _RecoElectron1EcalIsoDeltaRCone = iConfig.getParameter<double>("RecoElectron1EcalIsoDeltaRCone");
  _RecoElectron1EcalIsoRecHitThreshold = iConfig.getParameter<double>("RecoElectron1EcalIsoRecHitThreshold");
  _DoRecoElectron1DiscrByIp = iConfig.getParameter<bool>("DoRecoElectron1DiscrByIp");
//  _RecoElectron1IpCut = iConfig.getParameter<double>("RecoElectron1IpCut");
  _DoRecoElectron1DiscrByEoverP = iConfig.getParameter<bool>("DoRecoElectron1DiscrByEoverP");
//  _RecoElectron1EoverPMax = iConfig.getParameter<double>("RecoElectron1EoverPMax");
//  _RecoElectron1EoverPMin = iConfig.getParameter<double>("RecoElectron1EoverPMin");
  _DoRecoElectron1DiscrByHoverEm = iConfig.getParameter<bool>("DoRecoElectron1DiscrByHoverEm");
//  _RecoElectron1HoverEmCut = iConfig.getParameter<double>("RecoElectron1HoverEmCut");
  _DoRecoElectron1DiscrBySigmaIEtaIEta = iConfig.getParameter<bool>("DoRecoElectron1DiscrBySigmaIEtaIEta");
//  _RecoElectron1SigmaIEtaIEta = iConfig.getParameter<double>("RecoElectron1SigmaIEtaIEta");
  _DoRecoElectron1DiscrByDEtaIn = iConfig.getParameter<bool>("DoRecoElectron1DiscrByDEtaIn");
//  _RecoElectron1EEDEtaIn = iConfig.getParameter<double>("RecoElectron1EEDEtaIn");
//  _RecoElectron1EBDEtaIn = iConfig.getParameter<double>("RecoElectron1EBDEtaIn");
  _DoRecoElectron1DiscrByDPhiIn = iConfig.getParameter<bool>("DoRecoElectron1DiscrByDPhiIn");
//  _RecoElectron1EEDPhiIn = iConfig.getParameter<double>("RecoElectron1EEDPhiIn");
//  _RecoElectron1EBDPhiIn = iConfig.getParameter<double>("RecoElectron1EBDPhiIn");
  _DoRecoElectron1DiscrBySCE2by5Over5by5 = iConfig.getParameter<bool>("DoRecoElectron1DiscrBySCE2by5Over5by5");
//  _RecoElectron1EBscE1by5Over5by5 = iConfig.getParameter<double>("RecoElectron1EBscE1by5Over5by5");
//  _RecoElectron1EBscE2by5Over5by5 = iConfig.getParameter<double>("RecoElectron1EBscE2by5Over5by5");
  _DoRecoElectron1DiscrByMissingHits = iConfig.getParameter<bool>("DoRecoElectron1DiscrByMissingHits");
//  _RecoElectron1MissingHits = iConfig.getParameter<int>("RecoElectron1MissingHits");
//  _DoRecoElectron1DiscrByEcalDrivenSeed = iConfig.getParameter<bool>("DoRecoElectron1DiscrByEcalDrivenSeed");
//  _DoRecoElectron1DiscrByTrackerDrivenSeed = iConfig.getParameter<bool>("DoRecoElectron1DiscrByTrackerDrivenSeed");

  _RecoElectron2EtaCut = iConfig.getParameter<double>("RecoElectron2EtaCut");
  _RecoElectron2PtMinCut = iConfig.getParameter<double>("RecoElectron2PtMinCut");
  _RecoElectron2PtMaxCut = iConfig.getParameter<double>("RecoElectron2PtMaxCut");
  _DoRecoElectron2DiscrByTrackIsolation = iConfig.getParameter<bool>("DoRecoElectron2DiscrByTrackIsolation");
  _RecoElectron2TrackIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoElectron2TrackIsoSumPtMaxCutValue");
  _RecoElectron2TrackIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoElectron2TrackIsoSumPtMinCutValue");
  _RecoElectron2TrackIsoDeltaRCone = iConfig.getParameter<double>("RecoElectron2TrackIsoDeltaRCone");
  _RecoElectron2TrackIsoTrkThreshold = iConfig.getParameter<double>("RecoElectron2TrackIsoTrkThreshold");
  _DoRecoElectron2DiscrByEcalIsolation = iConfig.getParameter<bool>("DoRecoElectron2DiscrByEcalIsolation");
  _RecoElectron2EcalIsoSumPtMaxCutValue = iConfig.getParameter<double>("RecoElectron2EcalIsoSumPtMaxCutValue");
  _RecoElectron2EcalIsoSumPtMinCutValue = iConfig.getParameter<double>("RecoElectron2EcalIsoSumPtMinCutValue");
  _RecoElectron2EcalIsoDeltaRCone = iConfig.getParameter<double>("RecoElectron2EcalIsoDeltaRCone");
  _RecoElectron2EcalIsoRecHitThreshold = iConfig.getParameter<double>("RecoElectron2EcalIsoRecHitThreshold");
  _DoRecoElectron2DiscrByIp = iConfig.getParameter<bool>("DoRecoElectron2DiscrByIp");
//  _RecoElectron2IpCut = iConfig.getParameter<double>("RecoElectron2IpCut");
  _DoRecoElectron2DiscrByEoverP = iConfig.getParameter<bool>("DoRecoElectron2DiscrByEoverP");
//  _RecoElectron2EoverPMax = iConfig.getParameter<double>("RecoElectron2EoverPMax");
//  _RecoElectron2EoverPMin = iConfig.getParameter<double>("RecoElectron2EoverPMin");
  _DoRecoElectron2DiscrByHoverEm = iConfig.getParameter<bool>("DoRecoElectron2DiscrByHoverEm");
//  _RecoElectron2HoverEmCut = iConfig.getParameter<double>("RecoElectron2HoverEmCut");
  _DoRecoElectron2DiscrBySigmaIEtaIEta = iConfig.getParameter<bool>("DoRecoElectron2DiscrBySigmaIEtaIEta");
//  _RecoElectron2SigmaIEtaIEta = iConfig.getParameter<double>("RecoElectron2SigmaIEtaIEta");
  _DoRecoElectron2DiscrByDEtaIn = iConfig.getParameter<bool>("DoRecoElectron2DiscrByDEtaIn");
//  _RecoElectron2EEDEtaIn = iConfig.getParameter<double>("RecoElectron2EEDEtaIn");
//  _RecoElectron2EBDEtaIn = iConfig.getParameter<double>("RecoElectron2EBDEtaIn");
  _DoRecoElectron2DiscrByDPhiIn = iConfig.getParameter<bool>("DoRecoElectron2DiscrByDPhiIn");
//  _RecoElectron2EEDPhiIn = iConfig.getParameter<double>("RecoElectron2EEDPhiIn");
//  _RecoElectron2EBDPhiIn = iConfig.getParameter<double>("RecoElectron2EBDPhiIn");
  _DoRecoElectron2DiscrBySCE2by5Over5by5 = iConfig.getParameter<bool>("DoRecoElectron2DiscrBySCE2by5Over5by5");
//  _RecoElectron2EBscE1by5Over5by5 = iConfig.getParameter<double>("RecoElectron2EBscE1by5Over5by5");
//  _RecoElectron2EBscE2by5Over5by5 = iConfig.getParameter<double>("RecoElectron2EBscE2by5Over5by5");
  _DoRecoElectron2DiscrByMissingHits = iConfig.getParameter<bool>("DoRecoElectron2DiscrByMissingHits");
//  _RecoElectron2MissingHits = iConfig.getParameter<int>("RecoElectron2MissingHits");
//  _DoRecoElectron2DiscrByEcalDrivenSeed = iConfig.getParameter<bool>("DoRecoElectron2DiscrByEcalDrivenSeed");
//  _DoRecoElectron2DiscrByTrackerDrivenSeed = iConfig.getParameter<bool>("DoRecoElectron2DiscrByTrackerDrivenSeed");

  //-----Reco Jet Inputs
  _RecoJetSource = iConfig.getParameter<InputTag>("RecoJetSource");
  _RecoJetEtaMinCut = iConfig.getParameter<double>("RecoJetEtaMinCut");
  _RecoJetEtaMaxCut = iConfig.getParameter<double>("RecoJetEtaMaxCut");
  _RecoJetPtCut = iConfig.getParameter<double>("RecoJetPtCut");
  _UseCorrectedJet = iConfig.getParameter<bool>("UseCorrectedJet");
  _RemoveJetOverlapWithMuon1s = iConfig.getParameter<bool>("RemoveJetOverlapWithMuon1s");
  _JetMuon1MatchingDeltaR = iConfig.getParameter<double>("JetMuon1MatchingDeltaR");
  _RemoveJetOverlapWithElectron1s = iConfig.getParameter<bool>("RemoveJetOverlapWithElectron1s");
  _JetElectron1MatchingDeltaR = iConfig.getParameter<double>("JetElectron1MatchingDeltaR");
  _RemoveJetOverlapWithTau1s = iConfig.getParameter<bool>("RemoveJetOverlapWithTau1s");
  _JetTau1MatchingDeltaR = iConfig.getParameter<double>("JetTau1MatchingDeltaR");
  _RemoveJetOverlapWithMuon2s = iConfig.getParameter<bool>("RemoveJetOverlapWithMuon2s");
  _JetMuon2MatchingDeltaR = iConfig.getParameter<double>("JetMuon2MatchingDeltaR");
  _RemoveJetOverlapWithElectron2s = iConfig.getParameter<bool>("RemoveJetOverlapWithElectron2s");
  _JetElectron2MatchingDeltaR = iConfig.getParameter<double>("JetElectron2MatchingDeltaR");
  _RemoveJetOverlapWithTau2s = iConfig.getParameter<bool>("RemoveJetOverlapWithTau2s");
  _JetTau2MatchingDeltaR = iConfig.getParameter<double>("JetTau2MatchingDeltaR");
  _DoDiscrByFirstLeadingJet = iConfig.getParameter<bool>("DoDiscrByFirstLeadingJet");
  _RecoFirstLeadingJetPt = iConfig.getParameter<double>("RecoFirstLeadingJetPt");
  _RecoFirstLeadingJetEtaMinCut = iConfig.getParameter<double>("RecoFirstLeadingJetEtaMinCut");
  _RecoFirstLeadingJetEtaMaxCut = iConfig.getParameter<double>("RecoFirstLeadingJetEtaMaxCut");
  _RemoveFirstLeadingJetOverlapWithMuon1s = iConfig.getParameter<bool>("RemoveFirstLeadingJetOverlapWithMuon1s");
  _FirstLeadingJetMuon1MatchingDeltaR = iConfig.getParameter<double>("FirstLeadingJetMuon1MatchingDeltaR");
  _RemoveFirstLeadingJetOverlapWithElectron1s = iConfig.getParameter<bool>("RemoveFirstLeadingJetOverlapWithElectron1s");
  _FirstLeadingJetElectron1MatchingDeltaR = iConfig.getParameter<double>("FirstLeadingJetElectron1MatchingDeltaR");
  _RemoveFirstLeadingJetOverlapWithTau1s = iConfig.getParameter<bool>("RemoveFirstLeadingJetOverlapWithTau1s");
  _FirstLeadingJetTau1MatchingDeltaR = iConfig.getParameter<double>("FirstLeadingJetTau1MatchingDeltaR");
  _RemoveFirstLeadingJetOverlapWithMuon2s = iConfig.getParameter<bool>("RemoveFirstLeadingJetOverlapWithMuon2s");
  _FirstLeadingJetMuon2MatchingDeltaR = iConfig.getParameter<double>("FirstLeadingJetMuon2MatchingDeltaR");
  _RemoveFirstLeadingJetOverlapWithElectron2s = iConfig.getParameter<bool>("RemoveFirstLeadingJetOverlapWithElectron2s");
  _FirstLeadingJetElectron2MatchingDeltaR = iConfig.getParameter<double>("FirstLeadingJetElectron2MatchingDeltaR");
  _RemoveFirstLeadingJetOverlapWithTau2s = iConfig.getParameter<bool>("RemoveFirstLeadingJetOverlapWithTau2s");
  _FirstLeadingJetTau2MatchingDeltaR = iConfig.getParameter<double>("FirstLeadingJetTau2MatchingDeltaR");
  _DoDiscrBySecondLeadingJet = iConfig.getParameter<bool>("DoDiscrBySecondLeadingJet");
  _RecoSecondLeadingJetPt = iConfig.getParameter<double>("RecoSecondLeadingJetPt");
  _RecoSecondLeadingJetEtaMinCut = iConfig.getParameter<double>("RecoSecondLeadingJetEtaMinCut");
  _RecoSecondLeadingJetEtaMaxCut = iConfig.getParameter<double>("RecoSecondLeadingJetEtaMaxCut");
  _RemoveSecondLeadingJetOverlapWithMuon1s = iConfig.getParameter<bool>("RemoveSecondLeadingJetOverlapWithMuon1s");
  _SecondLeadingJetMuon1MatchingDeltaR = iConfig.getParameter<double>("SecondLeadingJetMuon1MatchingDeltaR");
  _RemoveSecondLeadingJetOverlapWithElectron1s = iConfig.getParameter<bool>("RemoveSecondLeadingJetOverlapWithElectron1s");
  _SecondLeadingJetElectron1MatchingDeltaR = iConfig.getParameter<double>("SecondLeadingJetElectron1MatchingDeltaR");
  _RemoveSecondLeadingJetOverlapWithTau1s = iConfig.getParameter<bool>("RemoveSecondLeadingJetOverlapWithTau1s");
  _SecondLeadingJetTau1MatchingDeltaR = iConfig.getParameter<double>("SecondLeadingJetTau1MatchingDeltaR");
  _RemoveSecondLeadingJetOverlapWithMuon2s = iConfig.getParameter<bool>("RemoveSecondLeadingJetOverlapWithMuon2s");
  _SecondLeadingJetMuon2MatchingDeltaR = iConfig.getParameter<double>("SecondLeadingJetMuon2MatchingDeltaR");
  _RemoveSecondLeadingJetOverlapWithElectron2s = iConfig.getParameter<bool>("RemoveSecondLeadingJetOverlapWithElectron2s");
  _SecondLeadingJetElectron2MatchingDeltaR = iConfig.getParameter<double>("SecondLeadingJetElectron2MatchingDeltaR");
  _RemoveSecondLeadingJetOverlapWithTau2s = iConfig.getParameter<bool>("RemoveSecondLeadingJetOverlapWithTau2s");
  _SecondLeadingJetTau2MatchingDeltaR = iConfig.getParameter<double>("SecondLeadingJetTau2MatchingDeltaR");

  //-----Reco b-Jet Inputs
  _RecoBJetEtaMinCut = iConfig.getParameter<double>("RecoBJetEtaMinCut");
  _RecoBJetEtaMaxCut = iConfig.getParameter<double>("RecoBJetEtaMaxCut");
  _RecoBJetPtCut = iConfig.getParameter<double>("RecoBJetPtCut");
  _RemoveBJetOverlapWithMuon1s = iConfig.getParameter<bool>("RemoveBJetOverlapWithMuon1s");
  _BJetMuon1MatchingDeltaR = iConfig.getParameter<double>("BJetMuon1MatchingDeltaR");
  _RemoveBJetOverlapWithElectron1s = iConfig.getParameter<bool>("RemoveBJetOverlapWithElectron1s");
  _BJetElectron1MatchingDeltaR = iConfig.getParameter<double>("BJetElectron1MatchingDeltaR");
  _RemoveBJetOverlapWithTau1s = iConfig.getParameter<bool>("RemoveBJetOverlapWithTau1s");
  _BJetTau1MatchingDeltaR = iConfig.getParameter<double>("BJetTau1MatchingDeltaR");
  _RemoveBJetOverlapWithMuon2s = iConfig.getParameter<bool>("RemoveBJetOverlapWithMuon2s");
  _BJetMuon2MatchingDeltaR = iConfig.getParameter<double>("BJetMuon2MatchingDeltaR");
  _RemoveBJetOverlapWithElectron2s = iConfig.getParameter<bool>("RemoveBJetOverlapWithElectron2s");
  _BJetElectron2MatchingDeltaR = iConfig.getParameter<double>("BJetElectron2MatchingDeltaR");
  _RemoveBJetOverlapWithTau2s = iConfig.getParameter<bool>("RemoveBJetOverlapWithTau2s");
  _BJetTau2MatchingDeltaR = iConfig.getParameter<double>("BJetTau2MatchingDeltaR");
  _ApplyJetBTagging = iConfig.getParameter<bool>("ApplyJetBTagging");
  _bTagger = iConfig.getUntrackedParameter<string>("bTagger");
  _JetBTaggingCut = iConfig.getParameter<double>("JetBTaggingCut");

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
  _DoDiscrByDiTauMassReco = iConfig.getParameter<bool>("DoDiscrByDiTauMassReco");
  _DiTauMassMinCut = iConfig.getParameter<double>("DiTauMassMinCut");
  _DiTauMassMaxCut = iConfig.getParameter<double>("DiTauMassMaxCut");
  _DoDiTauDiscrByCDFzeta2D = iConfig.getParameter<bool>("DoDiTauDiscrByCDFzeta2D");
  _DiTauPZetaCutCoefficient = iConfig.getParameter<double>("DiTauPZetaCutCoefficient");
  _DiTauPZetaVisCutCoefficient = iConfig.getParameter<double>("DiTauPZetaVisCutCoefficient");
  _DiTauCDFzeta2DCutValue = iConfig.getParameter<double>("DiTauCDFzeta2DCutValue");
  _DoDiTauDiscrByDeltaPtDivSumPt = iConfig.getParameter<bool>("DoDiTauDiscrByDeltaPtDivSumPt");
  _DiTauDeltaPtDivSumPtMinCutValue = iConfig.getParameter<double>("DiTauDeltaPtDivSumPtMinCutValue");
  _DiTauDeltaPtDivSumPtMaxCutValue = iConfig.getParameter<double>("DiTauDeltaPtDivSumPtMaxCutValue");
  _DoDiTauDiscrByDeltaPt = iConfig.getParameter<bool>("DoDiTauDiscrByDeltaPt");
  _DiTauDeltaPtMinCutValue = iConfig.getParameter<double>("DiTauDeltaPtMinCutValue");
  _DiTauDeltaPtMaxCutValue = iConfig.getParameter<double>("DiTauDeltaPtMaxCutValue");
  _DoDiMuonDiscrByDeltaR = iConfig.getParameter<bool>("DoDiMuonDiscrByDeltaR");
  _DiMuonDeltaRCut = iConfig.getParameter<double>("DiMuonDeltaRCut");
  _DiMuonDiscrByOSLSType = iConfig.getParameter<string>("DiMuonDiscrByOSLSType");
  _DoDiMuonDiscrByCosDphi = iConfig.getParameter<bool>("DoDiMuonDiscrByCosDphi");
  _DiMuonCosDphiMinCut = iConfig.getParameter<double>("DiMuonCosDphiMinCut");
  _DiMuonCosDphiMaxCut = iConfig.getParameter<double>("DiMuonCosDphiMaxCut");
  _DoDiscrByDiMuonMassReco = iConfig.getParameter<bool>("DoDiscrByDiMuonMassReco");
  _DiMuonMassMinCut = iConfig.getParameter<double>("DiMuonMassMinCut");
  _DiMuonMassMaxCut = iConfig.getParameter<double>("DiMuonMassMaxCut");
  _DoDiMuonDiscrByCDFzeta2D = iConfig.getParameter<bool>("DoDiMuonDiscrByCDFzeta2D");
  _DiMuonPZetaCutCoefficient = iConfig.getParameter<double>("DiMuonPZetaCutCoefficient");
  _DiMuonPZetaVisCutCoefficient = iConfig.getParameter<double>("DiMuonPZetaVisCutCoefficient");
  _DiMuonCDFzeta2DCutValue = iConfig.getParameter<double>("DiMuonCDFzeta2DCutValue");
  _DoDiMuonDiscrByDeltaPtDivSumPt = iConfig.getParameter<bool>("DoDiMuonDiscrByDeltaPtDivSumPt");
  _DiMuonDeltaPtDivSumPtMinCutValue = iConfig.getParameter<double>("DiMuonDeltaPtDivSumPtMinCutValue");
  _DiMuonDeltaPtDivSumPtMaxCutValue = iConfig.getParameter<double>("DiMuonDeltaPtDivSumPtMaxCutValue");
  _DoDiMuonDiscrByDeltaPt = iConfig.getParameter<bool>("DoDiMuonDiscrByDeltaPt");
  _DiMuonDeltaPtMinCutValue = iConfig.getParameter<double>("DiMuonDeltaPtMinCutValue");
  _DiMuonDeltaPtMaxCutValue = iConfig.getParameter<double>("DiMuonDeltaPtMaxCutValue");
  _DoDiElectronDiscrByDeltaR = iConfig.getParameter<bool>("DoDiElectronDiscrByDeltaR");
  _DiElectronDeltaRCut = iConfig.getParameter<double>("DiElectronDeltaRCut");
  _DiElectronDiscrByOSLSType = iConfig.getParameter<string>("DiElectronDiscrByOSLSType");
  _DoDiElectronDiscrByCosDphi = iConfig.getParameter<bool>("DoDiElectronDiscrByCosDphi");
  _DiElectronCosDphiMinCut = iConfig.getParameter<double>("DiElectronCosDphiMinCut");
  _DiElectronCosDphiMaxCut = iConfig.getParameter<double>("DiElectronCosDphiMaxCut");
  _DoDiscrByDiElectronMassReco = iConfig.getParameter<bool>("DoDiscrByDiElectronMassReco");
  _DiElectronMassMinCut = iConfig.getParameter<double>("DiElectronMassMinCut");
  _DiElectronMassMaxCut = iConfig.getParameter<double>("DiElectronMassMaxCut");
  _DoDiElectronDiscrByCDFzeta2D = iConfig.getParameter<bool>("DoDiElectronDiscrByCDFzeta2D");
  _DiElectronPZetaCutCoefficient = iConfig.getParameter<double>("DiElectronPZetaCutCoefficient");
  _DiElectronPZetaVisCutCoefficient = iConfig.getParameter<double>("DiElectronPZetaVisCutCoefficient");
  _DiElectronCDFzeta2DCutValue = iConfig.getParameter<double>("DiElectronCDFzeta2DCutValue");
  _DoDiElectronDiscrByDeltaPtDivSumPt = iConfig.getParameter<bool>("DoDiElectronDiscrByDeltaPtDivSumPt");
  _DiElectronDeltaPtDivSumPtMinCutValue = iConfig.getParameter<double>("DiElectronDeltaPtDivSumPtMinCutValue");
  _DiElectronDeltaPtDivSumPtMaxCutValue = iConfig.getParameter<double>("DiElectronDeltaPtDivSumPtMaxCutValue");
  _DoDiElectronDiscrByDeltaPt = iConfig.getParameter<bool>("DoDiElectronDiscrByDeltaPt");
  _DiElectronDeltaPtMinCutValue = iConfig.getParameter<double>("DiElectronDeltaPtMinCutValue");
  _DiElectronDeltaPtMaxCutValue = iConfig.getParameter<double>("DiElectronDeltaPtMaxCutValue");
  _DoMuon1Tau1DiscrByDeltaR = iConfig.getParameter<bool>("DoMuon1Tau1DiscrByDeltaR");
  _Muon1Tau1DeltaRCut = iConfig.getParameter<double>("Muon1Tau1DeltaRCut");
  _UseTauSeedTrackForMuon1Tau1DiscrByOSLS = iConfig.getParameter<bool>("UseTauSeedTrackForMuon1Tau1DiscrByOSLS");
  _Muon1Tau1DiscrByOSLSType = iConfig.getParameter<string>("Muon1Tau1DiscrByOSLSType");
  _DoMuon1Tau1DiscrByCosDphi = iConfig.getParameter<bool>("DoMuon1Tau1DiscrByCosDphi");
  _Muon1Tau1CosDphiMinCut = iConfig.getParameter<double>("Muon1Tau1CosDphiMinCut");
  _Muon1Tau1CosDphiMaxCut = iConfig.getParameter<double>("Muon1Tau1CosDphiMaxCut");
  _DoDiscrByMuon1Tau1MassReco = iConfig.getParameter<bool>("DoDiscrByMuon1Tau1MassReco");
  _Muon1Tau1MassMinCut = iConfig.getParameter<double>("Muon1Tau1MassMinCut");
  _Muon1Tau1MassMaxCut = iConfig.getParameter<double>("Muon1Tau1MassMaxCut");
  _DoMuon1Tau1DiscrByCDFzeta2D = iConfig.getParameter<bool>("DoMuon1Tau1DiscrByCDFzeta2D");
  _Muon1Tau1PZetaCutCoefficient = iConfig.getParameter<double>("Muon1Tau1PZetaCutCoefficient");
  _Muon1Tau1PZetaVisCutCoefficient = iConfig.getParameter<double>("Muon1Tau1PZetaVisCutCoefficient");
  _Muon1Tau1CDFzeta2DCutValue = iConfig.getParameter<double>("Muon1Tau1CDFzeta2DCutValue");
  _DoMuon1Tau1DiscrByDeltaPtDivSumPt = iConfig.getParameter<bool>("DoMuon1Tau1DiscrByDeltaPtDivSumPt");
  _Muon1Tau1DeltaPtDivSumPtMinCutValue = iConfig.getParameter<double>("Muon1Tau1DeltaPtDivSumPtMinCutValue");
  _Muon1Tau1DeltaPtDivSumPtMaxCutValue = iConfig.getParameter<double>("Muon1Tau1DeltaPtDivSumPtMaxCutValue");
  _DoMuon1Tau1DiscrByDeltaPt = iConfig.getParameter<bool>("DoMuon1Tau1DiscrByDeltaPt");
  _Muon1Tau1DeltaPtMinCutValue = iConfig.getParameter<double>("Muon1Tau1DeltaPtMinCutValue");
  _Muon1Tau1DeltaPtMaxCutValue = iConfig.getParameter<double>("Muon1Tau1DeltaPtMaxCutValue");
  _DoMuon1Tau2DiscrByDeltaR = iConfig.getParameter<bool>("DoMuon1Tau2DiscrByDeltaR");
  _Muon1Tau2DeltaRCut = iConfig.getParameter<double>("Muon1Tau2DeltaRCut");
  _UseTauSeedTrackForMuon1Tau2DiscrByOSLS = iConfig.getParameter<bool>("UseTauSeedTrackForMuon1Tau2DiscrByOSLS");
  _Muon1Tau2DiscrByOSLSType = iConfig.getParameter<string>("Muon1Tau2DiscrByOSLSType");
  _DoMuon1Tau2DiscrByCosDphi = iConfig.getParameter<bool>("DoMuon1Tau2DiscrByCosDphi");
  _Muon1Tau2CosDphiMinCut = iConfig.getParameter<double>("Muon1Tau2CosDphiMinCut");
  _Muon1Tau2CosDphiMaxCut = iConfig.getParameter<double>("Muon1Tau2CosDphiMaxCut");
  _DoDiscrByMuon1Tau2MassReco = iConfig.getParameter<bool>("DoDiscrByMuon1Tau2MassReco");
  _Muon1Tau2MassMinCut = iConfig.getParameter<double>("Muon1Tau2MassMinCut");
  _Muon1Tau2MassMaxCut = iConfig.getParameter<double>("Muon1Tau2MassMaxCut");
  _DoMuon1Tau2DiscrByCDFzeta2D = iConfig.getParameter<bool>("DoMuon1Tau2DiscrByCDFzeta2D");
  _Muon1Tau2PZetaCutCoefficient = iConfig.getParameter<double>("Muon1Tau2PZetaCutCoefficient");
  _Muon1Tau2PZetaVisCutCoefficient = iConfig.getParameter<double>("Muon1Tau2PZetaVisCutCoefficient");
  _Muon1Tau2CDFzeta2DCutValue = iConfig.getParameter<double>("Muon1Tau2CDFzeta2DCutValue");
  _DoMuon1Tau2DiscrByDeltaPtDivSumPt = iConfig.getParameter<bool>("DoMuon1Tau2DiscrByDeltaPtDivSumPt");
  _Muon1Tau2DeltaPtDivSumPtMinCutValue = iConfig.getParameter<double>("Muon1Tau2DeltaPtDivSumPtMinCutValue");
  _Muon1Tau2DeltaPtDivSumPtMaxCutValue = iConfig.getParameter<double>("Muon1Tau2DeltaPtDivSumPtMaxCutValue");
  _DoMuon1Tau2DiscrByDeltaPt = iConfig.getParameter<bool>("DoMuon1Tau2DiscrByDeltaPt");
  _Muon1Tau2DeltaPtMinCutValue = iConfig.getParameter<double>("Muon1Tau2DeltaPtMinCutValue");
  _Muon1Tau2DeltaPtMaxCutValue = iConfig.getParameter<double>("Muon1Tau2DeltaPtMaxCutValue");
  _DoMuon2Tau1DiscrByDeltaR = iConfig.getParameter<bool>("DoMuon2Tau1DiscrByDeltaR");
  _Muon2Tau1DeltaRCut = iConfig.getParameter<double>("Muon2Tau1DeltaRCut");
  _UseTauSeedTrackForMuon2Tau1DiscrByOSLS = iConfig.getParameter<bool>("UseTauSeedTrackForMuon2Tau1DiscrByOSLS");
  _Muon2Tau1DiscrByOSLSType = iConfig.getParameter<string>("Muon2Tau1DiscrByOSLSType");
  _DoMuon2Tau1DiscrByCosDphi = iConfig.getParameter<bool>("DoMuon2Tau1DiscrByCosDphi");
  _Muon2Tau1CosDphiMinCut = iConfig.getParameter<double>("Muon2Tau1CosDphiMinCut");
  _Muon2Tau1CosDphiMaxCut = iConfig.getParameter<double>("Muon2Tau1CosDphiMaxCut");
  _DoDiscrByMuon2Tau1MassReco = iConfig.getParameter<bool>("DoDiscrByMuon2Tau1MassReco");
  _Muon2Tau1MassMinCut = iConfig.getParameter<double>("Muon2Tau1MassMinCut");
  _Muon2Tau1MassMaxCut = iConfig.getParameter<double>("Muon2Tau1MassMaxCut");
  _DoMuon2Tau1DiscrByCDFzeta2D = iConfig.getParameter<bool>("DoMuon2Tau1DiscrByCDFzeta2D");
  _Muon2Tau1PZetaCutCoefficient = iConfig.getParameter<double>("Muon2Tau1PZetaCutCoefficient");
  _Muon2Tau1PZetaVisCutCoefficient = iConfig.getParameter<double>("Muon2Tau1PZetaVisCutCoefficient");
  _Muon2Tau1CDFzeta2DCutValue = iConfig.getParameter<double>("Muon2Tau1CDFzeta2DCutValue");
  _DoMuon2Tau1DiscrByDeltaPtDivSumPt = iConfig.getParameter<bool>("DoMuon2Tau1DiscrByDeltaPtDivSumPt");
  _Muon2Tau1DeltaPtDivSumPtMinCutValue = iConfig.getParameter<double>("Muon2Tau1DeltaPtDivSumPtMinCutValue");
  _Muon2Tau1DeltaPtDivSumPtMaxCutValue = iConfig.getParameter<double>("Muon2Tau1DeltaPtDivSumPtMaxCutValue");
  _DoMuon2Tau1DiscrByDeltaPt = iConfig.getParameter<bool>("DoMuon2Tau1DiscrByDeltaPt");
  _Muon2Tau1DeltaPtMinCutValue = iConfig.getParameter<double>("Muon2Tau1DeltaPtMinCutValue");
  _Muon2Tau1DeltaPtMaxCutValue = iConfig.getParameter<double>("Muon2Tau1DeltaPtMaxCutValue");
  _DoMuon2Tau2DiscrByDeltaR = iConfig.getParameter<bool>("DoMuon2Tau2DiscrByDeltaR");
  _Muon2Tau2DeltaRCut = iConfig.getParameter<double>("Muon2Tau2DeltaRCut");
  _UseTauSeedTrackForMuon2Tau2DiscrByOSLS = iConfig.getParameter<bool>("UseTauSeedTrackForMuon2Tau2DiscrByOSLS");
  _Muon2Tau2DiscrByOSLSType = iConfig.getParameter<string>("Muon2Tau2DiscrByOSLSType");
  _DoMuon2Tau2DiscrByCosDphi = iConfig.getParameter<bool>("DoMuon2Tau2DiscrByCosDphi");
  _Muon2Tau2CosDphiMinCut = iConfig.getParameter<double>("Muon2Tau2CosDphiMinCut");
  _Muon2Tau2CosDphiMaxCut = iConfig.getParameter<double>("Muon2Tau2CosDphiMaxCut");
  _DoDiscrByMuon2Tau2MassReco = iConfig.getParameter<bool>("DoDiscrByMuon2Tau2MassReco");
  _Muon2Tau2MassMinCut = iConfig.getParameter<double>("Muon2Tau2MassMinCut");
  _Muon2Tau2MassMaxCut = iConfig.getParameter<double>("Muon2Tau2MassMaxCut");
  _DoMuon2Tau2DiscrByCDFzeta2D = iConfig.getParameter<bool>("DoMuon2Tau2DiscrByCDFzeta2D");
  _Muon2Tau2PZetaCutCoefficient = iConfig.getParameter<double>("Muon2Tau2PZetaCutCoefficient");
  _Muon2Tau2PZetaVisCutCoefficient = iConfig.getParameter<double>("Muon2Tau2PZetaVisCutCoefficient");
  _Muon2Tau2CDFzeta2DCutValue = iConfig.getParameter<double>("Muon2Tau2CDFzeta2DCutValue");
  _DoMuon2Tau2DiscrByDeltaPtDivSumPt = iConfig.getParameter<bool>("DoMuon2Tau2DiscrByDeltaPtDivSumPt");
  _Muon2Tau2DeltaPtDivSumPtMinCutValue = iConfig.getParameter<double>("Muon2Tau2DeltaPtDivSumPtMinCutValue");
  _Muon2Tau2DeltaPtDivSumPtMaxCutValue = iConfig.getParameter<double>("Muon2Tau2DeltaPtDivSumPtMaxCutValue");
  _DoMuon2Tau2DiscrByDeltaPt = iConfig.getParameter<bool>("DoMuon2Tau2DiscrByDeltaPt");
  _Muon2Tau2DeltaPtMinCutValue = iConfig.getParameter<double>("Muon2Tau2DeltaPtMinCutValue");
  _Muon2Tau2DeltaPtMaxCutValue = iConfig.getParameter<double>("Muon2Tau2DeltaPtMaxCutValue");
  _DoElectron1Tau1DiscrByDeltaR = iConfig.getParameter<bool>("DoElectron1Tau1DiscrByDeltaR");
  _Electron1Tau1DeltaRCut = iConfig.getParameter<double>("Electron1Tau1DeltaRCut");
  _UseTauSeedTrackForElectron1Tau1DiscrByOSLS = iConfig.getParameter<bool>("UseTauSeedTrackForElectron1Tau1DiscrByOSLS");
  _Electron1Tau1DiscrByOSLSType = iConfig.getParameter<string>("Electron1Tau1DiscrByOSLSType");
  _DoElectron1Tau1DiscrByCosDphi = iConfig.getParameter<bool>("DoElectron1Tau1DiscrByCosDphi");
  _Electron1Tau1CosDphiMinCut = iConfig.getParameter<double>("Electron1Tau1CosDphiMinCut");
  _Electron1Tau1CosDphiMaxCut = iConfig.getParameter<double>("Electron1Tau1CosDphiMaxCut");
  _DoDiscrByElectron1Tau1MassReco = iConfig.getParameter<bool>("DoDiscrByElectron1Tau1MassReco");
  _Electron1Tau1MassMinCut = iConfig.getParameter<double>("Electron1Tau1MassMinCut");
  _Electron1Tau1MassMaxCut = iConfig.getParameter<double>("Electron1Tau1MassMaxCut");
  _DoElectron1Tau1DiscrByCDFzeta2D = iConfig.getParameter<bool>("DoElectron1Tau1DiscrByCDFzeta2D");
  _Electron1Tau1PZetaCutCoefficient = iConfig.getParameter<double>("Electron1Tau1PZetaCutCoefficient");
  _Electron1Tau1PZetaVisCutCoefficient = iConfig.getParameter<double>("Electron1Tau1PZetaVisCutCoefficient");
  _Electron1Tau1CDFzeta2DCutValue = iConfig.getParameter<double>("Electron1Tau1CDFzeta2DCutValue");
  _DoElectron1Tau1DiscrByDeltaPtDivSumPt = iConfig.getParameter<bool>("DoElectron1Tau1DiscrByDeltaPtDivSumPt");
  _Electron1Tau1DeltaPtDivSumPtMinCutValue = iConfig.getParameter<double>("Electron1Tau1DeltaPtDivSumPtMinCutValue");
  _Electron1Tau1DeltaPtDivSumPtMaxCutValue = iConfig.getParameter<double>("Electron1Tau1DeltaPtDivSumPtMaxCutValue");
  _DoElectron1Tau1DiscrByDeltaPt = iConfig.getParameter<bool>("DoElectron1Tau1DiscrByDeltaPt");
  _Electron1Tau1DeltaPtMinCutValue = iConfig.getParameter<double>("Electron1Tau1DeltaPtMinCutValue");
  _Electron1Tau1DeltaPtMaxCutValue = iConfig.getParameter<double>("Electron1Tau1DeltaPtMaxCutValue");
  _DoElectron1Tau2DiscrByDeltaR = iConfig.getParameter<bool>("DoElectron1Tau2DiscrByDeltaR");
  _Electron1Tau2DeltaRCut = iConfig.getParameter<double>("Electron1Tau2DeltaRCut");
  _UseTauSeedTrackForElectron1Tau2DiscrByOSLS = iConfig.getParameter<bool>("UseTauSeedTrackForElectron1Tau2DiscrByOSLS");
  _Electron1Tau2DiscrByOSLSType = iConfig.getParameter<string>("Electron1Tau2DiscrByOSLSType");
  _DoElectron1Tau2DiscrByCosDphi = iConfig.getParameter<bool>("DoElectron1Tau2DiscrByCosDphi");
  _Electron1Tau2CosDphiMinCut = iConfig.getParameter<double>("Electron1Tau2CosDphiMinCut");
  _Electron1Tau2CosDphiMaxCut = iConfig.getParameter<double>("Electron1Tau2CosDphiMaxCut");
  _DoDiscrByElectron1Tau2MassReco = iConfig.getParameter<bool>("DoDiscrByElectron1Tau2MassReco");
  _Electron1Tau2MassMinCut = iConfig.getParameter<double>("Electron1Tau2MassMinCut");
  _Electron1Tau2MassMaxCut = iConfig.getParameter<double>("Electron1Tau2MassMaxCut");
  _DoElectron1Tau2DiscrByCDFzeta2D = iConfig.getParameter<bool>("DoElectron1Tau2DiscrByCDFzeta2D");
  _Electron1Tau2PZetaCutCoefficient = iConfig.getParameter<double>("Electron1Tau2PZetaCutCoefficient");
  _Electron1Tau2PZetaVisCutCoefficient = iConfig.getParameter<double>("Electron1Tau2PZetaVisCutCoefficient");
  _Electron1Tau2CDFzeta2DCutValue = iConfig.getParameter<double>("Electron1Tau2CDFzeta2DCutValue");
  _DoElectron1Tau2DiscrByDeltaPtDivSumPt = iConfig.getParameter<bool>("DoElectron1Tau2DiscrByDeltaPtDivSumPt");
  _Electron1Tau2DeltaPtDivSumPtMinCutValue = iConfig.getParameter<double>("Electron1Tau2DeltaPtDivSumPtMinCutValue");
  _Electron1Tau2DeltaPtDivSumPtMaxCutValue = iConfig.getParameter<double>("Electron1Tau2DeltaPtDivSumPtMaxCutValue");
  _DoElectron1Tau2DiscrByDeltaPt = iConfig.getParameter<bool>("DoElectron1Tau2DiscrByDeltaPt");
  _Electron1Tau2DeltaPtMinCutValue = iConfig.getParameter<double>("Electron1Tau2DeltaPtMinCutValue");
  _Electron1Tau2DeltaPtMaxCutValue = iConfig.getParameter<double>("Electron1Tau2DeltaPtMaxCutValue");
  _DoElectron2Tau1DiscrByDeltaR = iConfig.getParameter<bool>("DoElectron2Tau1DiscrByDeltaR");
  _Electron2Tau1DeltaRCut = iConfig.getParameter<double>("Electron2Tau1DeltaRCut");
  _UseTauSeedTrackForElectron2Tau1DiscrByOSLS = iConfig.getParameter<bool>("UseTauSeedTrackForElectron2Tau1DiscrByOSLS");
  _Electron2Tau1DiscrByOSLSType = iConfig.getParameter<string>("Electron2Tau1DiscrByOSLSType");
  _DoElectron2Tau1DiscrByCosDphi = iConfig.getParameter<bool>("DoElectron2Tau1DiscrByCosDphi");
  _Electron2Tau1CosDphiMinCut = iConfig.getParameter<double>("Electron2Tau1CosDphiMinCut");
  _Electron2Tau1CosDphiMaxCut = iConfig.getParameter<double>("Electron2Tau1CosDphiMaxCut");
  _DoDiscrByElectron2Tau1MassReco = iConfig.getParameter<bool>("DoDiscrByElectron2Tau1MassReco");
  _Electron2Tau1MassMinCut = iConfig.getParameter<double>("Electron2Tau1MassMinCut");
  _Electron2Tau1MassMaxCut = iConfig.getParameter<double>("Electron2Tau1MassMaxCut");
  _DoElectron2Tau1DiscrByCDFzeta2D = iConfig.getParameter<bool>("DoElectron2Tau1DiscrByCDFzeta2D");
  _Electron2Tau1PZetaCutCoefficient = iConfig.getParameter<double>("Electron2Tau1PZetaCutCoefficient");
  _Electron2Tau1PZetaVisCutCoefficient = iConfig.getParameter<double>("Electron2Tau1PZetaVisCutCoefficient");
  _Electron2Tau1CDFzeta2DCutValue = iConfig.getParameter<double>("Electron2Tau1CDFzeta2DCutValue");
  _DoElectron2Tau1DiscrByDeltaPtDivSumPt = iConfig.getParameter<bool>("DoElectron2Tau1DiscrByDeltaPtDivSumPt");
  _Electron2Tau1DeltaPtDivSumPtMinCutValue = iConfig.getParameter<double>("Electron2Tau1DeltaPtDivSumPtMinCutValue");
  _Electron2Tau1DeltaPtDivSumPtMaxCutValue = iConfig.getParameter<double>("Electron2Tau1DeltaPtDivSumPtMaxCutValue");
  _DoElectron2Tau1DiscrByDeltaPt = iConfig.getParameter<bool>("DoElectron2Tau1DiscrByDeltaPt");
  _Electron2Tau1DeltaPtMinCutValue = iConfig.getParameter<double>("Electron2Tau1DeltaPtMinCutValue");
  _Electron2Tau1DeltaPtMaxCutValue = iConfig.getParameter<double>("Electron2Tau1DeltaPtMaxCutValue");
  _DoElectron2Tau2DiscrByDeltaR = iConfig.getParameter<bool>("DoElectron2Tau2DiscrByDeltaR");
  _Electron2Tau2DeltaRCut = iConfig.getParameter<double>("Electron2Tau2DeltaRCut");
  _UseTauSeedTrackForElectron2Tau2DiscrByOSLS = iConfig.getParameter<bool>("UseTauSeedTrackForElectron2Tau2DiscrByOSLS");
  _Electron2Tau2DiscrByOSLSType = iConfig.getParameter<string>("Electron2Tau2DiscrByOSLSType");
  _DoElectron2Tau2DiscrByCosDphi = iConfig.getParameter<bool>("DoElectron2Tau2DiscrByCosDphi");
  _Electron2Tau2CosDphiMinCut = iConfig.getParameter<double>("Electron2Tau2CosDphiMinCut");
  _Electron2Tau2CosDphiMaxCut = iConfig.getParameter<double>("Electron2Tau2CosDphiMaxCut");
  _DoDiscrByElectron2Tau2MassReco = iConfig.getParameter<bool>("DoDiscrByElectron2Tau2MassReco");
  _Electron2Tau2MassMinCut = iConfig.getParameter<double>("Electron2Tau2MassMinCut");
  _Electron2Tau2MassMaxCut = iConfig.getParameter<double>("Electron2Tau2MassMaxCut");
  _DoElectron2Tau2DiscrByCDFzeta2D = iConfig.getParameter<bool>("DoElectron2Tau2DiscrByCDFzeta2D");
  _Electron2Tau2PZetaCutCoefficient = iConfig.getParameter<double>("Electron2Tau2PZetaCutCoefficient");
  _Electron2Tau2PZetaVisCutCoefficient = iConfig.getParameter<double>("Electron2Tau2PZetaVisCutCoefficient");
  _Electron2Tau2CDFzeta2DCutValue = iConfig.getParameter<double>("Electron2Tau2CDFzeta2DCutValue");
  _DoElectron2Tau2DiscrByDeltaPtDivSumPt = iConfig.getParameter<bool>("DoElectron2Tau2DiscrByDeltaPtDivSumPt");
  _Electron2Tau2DeltaPtDivSumPtMinCutValue = iConfig.getParameter<double>("Electron2Tau2DeltaPtDivSumPtMinCutValue");
  _Electron2Tau2DeltaPtDivSumPtMaxCutValue = iConfig.getParameter<double>("Electron2Tau2DeltaPtDivSumPtMaxCutValue");
  _DoElectron2Tau2DiscrByDeltaPt = iConfig.getParameter<bool>("DoElectron2Tau2DiscrByDeltaPt");
  _Electron2Tau2DeltaPtMinCutValue = iConfig.getParameter<double>("Electron2Tau2DeltaPtMinCutValue");
  _Electron2Tau2DeltaPtMaxCutValue = iConfig.getParameter<double>("Electron2Tau2DeltaPtMaxCutValue");
  _DoElectron1Muon1DiscrByDeltaR = iConfig.getParameter<bool>("DoElectron1Muon1DiscrByDeltaR");
  _Electron1Muon1DeltaRCut = iConfig.getParameter<double>("Electron1Muon1DeltaRCut");
  _Electron1Muon1DiscrByOSLSType = iConfig.getParameter<string>("Electron1Muon1DiscrByOSLSType");
  _DoElectron1Muon1DiscrByCosDphi = iConfig.getParameter<bool>("DoElectron1Muon1DiscrByCosDphi");
  _Electron1Muon1CosDphiMinCut = iConfig.getParameter<double>("Electron1Muon1CosDphiMinCut");
  _Electron1Muon1CosDphiMaxCut = iConfig.getParameter<double>("Electron1Muon1CosDphiMaxCut");
  _DoDiscrByElectron1Muon1MassReco = iConfig.getParameter<bool>("DoDiscrByElectron1Muon1MassReco");
  _Electron1Muon1MassMinCut = iConfig.getParameter<double>("Electron1Muon1MassMinCut");
  _Electron1Muon1MassMaxCut = iConfig.getParameter<double>("Electron1Muon1MassMaxCut");
  _DoElectron1Muon1DiscrByCDFzeta2D = iConfig.getParameter<bool>("DoElectron1Muon1DiscrByCDFzeta2D");
  _Electron1Muon1PZetaCutCoefficient = iConfig.getParameter<double>("Electron1Muon1PZetaCutCoefficient");
  _Electron1Muon1PZetaVisCutCoefficient = iConfig.getParameter<double>("Electron1Muon1PZetaVisCutCoefficient");
  _Electron1Muon1CDFzeta2DCutValue = iConfig.getParameter<double>("Electron1Muon1CDFzeta2DCutValue");
  _DoElectron1Muon1DiscrByDeltaPtDivSumPt = iConfig.getParameter<bool>("DoElectron1Muon1DiscrByDeltaPtDivSumPt");
  _Electron1Muon1DeltaPtDivSumPtMinCutValue = iConfig.getParameter<double>("Electron1Muon1DeltaPtDivSumPtMinCutValue");
  _Electron1Muon1DeltaPtDivSumPtMaxCutValue = iConfig.getParameter<double>("Electron1Muon1DeltaPtDivSumPtMaxCutValue");
  _DoElectron1Muon1DiscrByDeltaPt = iConfig.getParameter<bool>("DoElectron1Muon1DiscrByDeltaPt");
  _Electron1Muon1DeltaPtMinCutValue = iConfig.getParameter<double>("Electron1Muon1DeltaPtMinCutValue");
  _Electron1Muon1DeltaPtMaxCutValue = iConfig.getParameter<double>("Electron1Muon1DeltaPtMaxCutValue");
  _DoElectron1Muon2DiscrByDeltaR = iConfig.getParameter<bool>("DoElectron1Muon2DiscrByDeltaR");
  _Electron1Muon2DeltaRCut = iConfig.getParameter<double>("Electron1Muon2DeltaRCut");
  _Electron1Muon2DiscrByOSLSType = iConfig.getParameter<string>("Electron1Muon2DiscrByOSLSType");
  _DoElectron1Muon2DiscrByCosDphi = iConfig.getParameter<bool>("DoElectron1Muon2DiscrByCosDphi");
  _Electron1Muon2CosDphiMinCut = iConfig.getParameter<double>("Electron1Muon2CosDphiMinCut");
  _Electron1Muon2CosDphiMaxCut = iConfig.getParameter<double>("Electron1Muon2CosDphiMaxCut");
  _DoDiscrByElectron1Muon2MassReco = iConfig.getParameter<bool>("DoDiscrByElectron1Muon2MassReco");
  _Electron1Muon2MassMinCut = iConfig.getParameter<double>("Electron1Muon2MassMinCut");
  _Electron1Muon2MassMaxCut = iConfig.getParameter<double>("Electron1Muon2MassMaxCut");
  _DoElectron1Muon2DiscrByCDFzeta2D = iConfig.getParameter<bool>("DoElectron1Muon2DiscrByCDFzeta2D");
  _Electron1Muon2PZetaCutCoefficient = iConfig.getParameter<double>("Electron1Muon2PZetaCutCoefficient");
  _Electron1Muon2PZetaVisCutCoefficient = iConfig.getParameter<double>("Electron1Muon2PZetaVisCutCoefficient");
  _Electron1Muon2CDFzeta2DCutValue = iConfig.getParameter<double>("Electron1Muon2CDFzeta2DCutValue");
  _DoElectron1Muon2DiscrByDeltaPtDivSumPt = iConfig.getParameter<bool>("DoElectron1Muon2DiscrByDeltaPtDivSumPt");
  _Electron1Muon2DeltaPtDivSumPtMinCutValue = iConfig.getParameter<double>("Electron1Muon2DeltaPtDivSumPtMinCutValue");
  _Electron1Muon2DeltaPtDivSumPtMaxCutValue = iConfig.getParameter<double>("Electron1Muon2DeltaPtDivSumPtMaxCutValue");
  _DoElectron1Muon2DiscrByDeltaPt = iConfig.getParameter<bool>("DoElectron1Muon2DiscrByDeltaPt");
  _Electron1Muon2DeltaPtMinCutValue = iConfig.getParameter<double>("Electron1Muon2DeltaPtMinCutValue");
  _Electron1Muon2DeltaPtMaxCutValue = iConfig.getParameter<double>("Electron1Muon2DeltaPtMaxCutValue");
  _DoElectron2Muon1DiscrByDeltaR = iConfig.getParameter<bool>("DoElectron2Muon1DiscrByDeltaR");
  _Electron2Muon1DeltaRCut = iConfig.getParameter<double>("Electron2Muon1DeltaRCut");
  _Electron2Muon1DiscrByOSLSType = iConfig.getParameter<string>("Electron2Muon1DiscrByOSLSType");
  _DoElectron2Muon1DiscrByCosDphi = iConfig.getParameter<bool>("DoElectron2Muon1DiscrByCosDphi");
  _Electron2Muon1CosDphiMinCut = iConfig.getParameter<double>("Electron2Muon1CosDphiMinCut");
  _Electron2Muon1CosDphiMaxCut = iConfig.getParameter<double>("Electron2Muon1CosDphiMaxCut");
  _DoDiscrByElectron2Muon1MassReco = iConfig.getParameter<bool>("DoDiscrByElectron2Muon1MassReco");
  _Electron2Muon1MassMinCut = iConfig.getParameter<double>("Electron2Muon1MassMinCut");
  _Electron2Muon1MassMaxCut = iConfig.getParameter<double>("Electron2Muon1MassMaxCut");
  _DoElectron2Muon1DiscrByCDFzeta2D = iConfig.getParameter<bool>("DoElectron2Muon1DiscrByCDFzeta2D");
  _Electron2Muon1PZetaCutCoefficient = iConfig.getParameter<double>("Electron2Muon1PZetaCutCoefficient");
  _Electron2Muon1PZetaVisCutCoefficient = iConfig.getParameter<double>("Electron2Muon1PZetaVisCutCoefficient");
  _Electron2Muon1CDFzeta2DCutValue = iConfig.getParameter<double>("Electron2Muon1CDFzeta2DCutValue");
  _DoElectron2Muon1DiscrByDeltaPtDivSumPt = iConfig.getParameter<bool>("DoElectron2Muon1DiscrByDeltaPtDivSumPt");
  _Electron2Muon1DeltaPtDivSumPtMinCutValue = iConfig.getParameter<double>("Electron2Muon1DeltaPtDivSumPtMinCutValue");
  _Electron2Muon1DeltaPtDivSumPtMaxCutValue = iConfig.getParameter<double>("Electron2Muon1DeltaPtDivSumPtMaxCutValue");
  _DoElectron2Muon1DiscrByDeltaPt = iConfig.getParameter<bool>("DoElectron2Muon1DiscrByDeltaPt");
  _Electron2Muon1DeltaPtMinCutValue = iConfig.getParameter<double>("Electron2Muon1DeltaPtMinCutValue");
  _Electron2Muon1DeltaPtMaxCutValue = iConfig.getParameter<double>("Electron2Muon1DeltaPtMaxCutValue");
  _DoElectron2Muon2DiscrByDeltaR = iConfig.getParameter<bool>("DoElectron2Muon2DiscrByDeltaR");
  _Electron2Muon2DeltaRCut = iConfig.getParameter<double>("Electron2Muon2DeltaRCut");
  _Electron2Muon2DiscrByOSLSType = iConfig.getParameter<string>("Electron2Muon2DiscrByOSLSType");
  _DoElectron2Muon2DiscrByCosDphi = iConfig.getParameter<bool>("DoElectron2Muon2DiscrByCosDphi");
  _Electron2Muon2CosDphiMinCut = iConfig.getParameter<double>("Electron2Muon2CosDphiMinCut");
  _Electron2Muon2CosDphiMaxCut = iConfig.getParameter<double>("Electron2Muon2CosDphiMaxCut");
  _DoDiscrByElectron2Muon2MassReco = iConfig.getParameter<bool>("DoDiscrByElectron2Muon2MassReco");
  _Electron2Muon2MassMinCut = iConfig.getParameter<double>("Electron2Muon2MassMinCut");
  _Electron2Muon2MassMaxCut = iConfig.getParameter<double>("Electron2Muon2MassMaxCut");
  _DoElectron2Muon2DiscrByCDFzeta2D = iConfig.getParameter<bool>("DoElectron2Muon2DiscrByCDFzeta2D");
  _Electron2Muon2PZetaCutCoefficient = iConfig.getParameter<double>("Electron2Muon2PZetaCutCoefficient");
  _Electron2Muon2PZetaVisCutCoefficient = iConfig.getParameter<double>("Electron2Muon2PZetaVisCutCoefficient");
  _Electron2Muon2CDFzeta2DCutValue = iConfig.getParameter<double>("Electron2Muon2CDFzeta2DCutValue");
  _DoElectron2Muon2DiscrByDeltaPtDivSumPt = iConfig.getParameter<bool>("DoElectron2Muon2DiscrByDeltaPtDivSumPt");
  _Electron2Muon2DeltaPtDivSumPtMinCutValue = iConfig.getParameter<double>("Electron2Muon2DeltaPtDivSumPtMinCutValue");
  _Electron2Muon2DeltaPtDivSumPtMaxCutValue = iConfig.getParameter<double>("Electron2Muon2DeltaPtDivSumPtMaxCutValue");
  _DoElectron2Muon2DiscrByDeltaPt = iConfig.getParameter<bool>("DoElectron2Muon2DiscrByDeltaPt");
  _Electron2Muon2DeltaPtMinCutValue = iConfig.getParameter<double>("Electron2Muon2DeltaPtMinCutValue");
  _Electron2Muon2DeltaPtMaxCutValue = iConfig.getParameter<double>("Electron2Muon2DeltaPtMaxCutValue");
  _UseVectorSumOfMuon1Tau1ProductsAndMetMassReco = iConfig.getParameter<bool>("UseVectorSumOfMuon1Tau1ProductsAndMetMassReco");
  _UseCollinerApproxMuon1Tau1MassReco = iConfig.getParameter<bool>("UseCollinerApproxMuon1Tau1MassReco");
  _UseVectorSumOfMuon1Tau2ProductsAndMetMassReco = iConfig.getParameter<bool>("UseVectorSumOfMuon1Tau2ProductsAndMetMassReco");
  _UseCollinerApproxMuon1Tau2MassReco = iConfig.getParameter<bool>("UseCollinerApproxMuon1Tau2MassReco");
  _UseVectorSumOfMuon2Tau1ProductsAndMetMassReco = iConfig.getParameter<bool>("UseVectorSumOfMuon2Tau1ProductsAndMetMassReco");
  _UseCollinerApproxMuon2Tau1MassReco = iConfig.getParameter<bool>("UseCollinerApproxMuon2Tau1MassReco");
  _UseVectorSumOfMuon2Tau2ProductsAndMetMassReco = iConfig.getParameter<bool>("UseVectorSumOfMuon2Tau2ProductsAndMetMassReco");
  _UseCollinerApproxMuon2Tau2MassReco = iConfig.getParameter<bool>("UseCollinerApproxMuon2Tau2MassReco");
  _UseVectorSumOfElectron1Tau1ProductsAndMetMassReco = iConfig.getParameter<bool>("UseVectorSumOfElectron1Tau1ProductsAndMetMassReco");
  _UseCollinerApproxElectron1Tau1MassReco = iConfig.getParameter<bool>("UseCollinerApproxElectron1Tau1MassReco");
  _UseVectorSumOfElectron1Tau2ProductsAndMetMassReco = iConfig.getParameter<bool>("UseVectorSumOfElectron1Tau2ProductsAndMetMassReco");
  _UseCollinerApproxElectron1Tau2MassReco = iConfig.getParameter<bool>("UseCollinerApproxElectron1Tau2MassReco");
  _UseVectorSumOfElectron2Tau1ProductsAndMetMassReco = iConfig.getParameter<bool>("UseVectorSumOfElectron2Tau1ProductsAndMetMassReco");
  _UseCollinerApproxElectron2Tau1MassReco = iConfig.getParameter<bool>("UseCollinerApproxElectron2Tau1MassReco");
  _UseVectorSumOfElectron2Tau2ProductsAndMetMassReco = iConfig.getParameter<bool>("UseVectorSumOfElectron2Tau2ProductsAndMetMassReco");
  _UseCollinerApproxElectron2Tau2MassReco = iConfig.getParameter<bool>("UseCollinerApproxElectron2Tau2MassReco");
  _UseVectorSumOfElectron1Muon1ProductsAndMetMassReco = iConfig.getParameter<bool>("UseVectorSumOfElectron1Muon1ProductsAndMetMassReco");
  _UseCollinerApproxElectron1Muon1MassReco = iConfig.getParameter<bool>("UseCollinerApproxElectron1Muon1MassReco");
  _UseVectorSumOfElectron1Muon2ProductsAndMetMassReco = iConfig.getParameter<bool>("UseVectorSumOfElectron1Muon2ProductsAndMetMassReco");
  _UseCollinerApproxElectron1Muon2MassReco = iConfig.getParameter<bool>("UseCollinerApproxElectron1Muon2MassReco");
  _UseVectorSumOfElectron2Muon1ProductsAndMetMassReco = iConfig.getParameter<bool>("UseVectorSumOfElectron2Muon1ProductsAndMetMassReco");
  _UseCollinerApproxElectron2Muon1MassReco = iConfig.getParameter<bool>("UseCollinerApproxElectron2Muon1MassReco");
  _UseVectorSumOfElectron2Muon2ProductsAndMetMassReco = iConfig.getParameter<bool>("UseVectorSumOfElectron2Muon2ProductsAndMetMassReco");
  _UseCollinerApproxElectron2Muon2MassReco = iConfig.getParameter<bool>("UseCollinerApproxElectron2Muon2MassReco");
  _UseVectorSumOfDiMuonProductsAndMetMassReco = iConfig.getParameter<bool>("UseVectorSumOfDiMuonProductsAndMetMassReco");
  _UseCollinerApproxDiMuonMassReco = iConfig.getParameter<bool>("UseCollinerApproxDiMuonMassReco");
  _UseVectorSumOfDiElectronProductsAndMetMassReco = iConfig.getParameter<bool>("UseVectorSumOfDiElectronProductsAndMetMassReco");
  _UseCollinerApproxDiElectronMassReco = iConfig.getParameter<bool>("UseCollinerApproxDiElectronMassReco");
  _UseVectorSumOfDiTauProductsAndMetMassReco = iConfig.getParameter<bool>("UseVectorSumOfDiTauProductsAndMetMassReco");
  _UseCollinerApproxDiTauMassReco = iConfig.getParameter<bool>("UseCollinerApproxDiTauMassReco");
  _DoDiscrByMuon1MetDphi = iConfig.getParameter<bool>("DoDiscrByMuon1MetDphi");
  _Muon1MetDphiMinCut = iConfig.getParameter<double>("Muon1MetDphiMinCut");
  _Muon1MetDphiMaxCut = iConfig.getParameter<double>("Muon1MetDphiMaxCut");
  _DoDiscrByMuon2MetDphi = iConfig.getParameter<bool>("DoDiscrByMuon2MetDphi");
  _Muon2MetDphiMinCut = iConfig.getParameter<double>("Muon2MetDphiMinCut");
  _Muon2MetDphiMaxCut = iConfig.getParameter<double>("Muon2MetDphiMaxCut");
  _DoDiscrByElectron1MetDphi = iConfig.getParameter<bool>("DoDiscrByElectron1MetDphi");
  _Electron1MetDphiMinCut = iConfig.getParameter<double>("Electron1MetDphiMinCut");
  _Electron1MetDphiMaxCut = iConfig.getParameter<double>("Electron1MetDphiMaxCut");
  _DoDiscrByElectron2MetDphi = iConfig.getParameter<bool>("DoDiscrByElectron2MetDphi");
  _Electron2MetDphiMinCut = iConfig.getParameter<double>("Electron2MetDphiMinCut");
  _Electron2MetDphiMaxCut = iConfig.getParameter<double>("Electron2MetDphiMaxCut");
  _DoDiscrByTau1MetDphi = iConfig.getParameter<bool>("DoDiscrByTau1MetDphi");
  _Tau1MetDphiMinCut = iConfig.getParameter<double>("Tau1MetDphiMinCut");
  _Tau1MetDphiMaxCut = iConfig.getParameter<double>("Tau1MetDphiMaxCut");
  _DoDiscrByTau2MetDphi = iConfig.getParameter<bool>("DoDiscrByTau2MetDphi");
  _Tau2MetDphiMinCut = iConfig.getParameter<double>("Tau2MetDphiMinCut");
  _Tau2MetDphiMaxCut = iConfig.getParameter<double>("Tau2MetDphiMaxCut");
  _DoDiscrByMuon1MetMt = iConfig.getParameter<bool>("DoDiscrByMuon1MetMt");
  _Muon1MetMtMinCut = iConfig.getParameter<double>("Muon1MetMtMinCut");
  _Muon1MetMtMaxCut = iConfig.getParameter<double>("Muon1MetMtMaxCut");
  _DoDiscrByMuon2MetMt = iConfig.getParameter<bool>("DoDiscrByMuon2MetMt");
  _Muon2MetMtMinCut = iConfig.getParameter<double>("Muon2MetMtMinCut");
  _Muon2MetMtMaxCut = iConfig.getParameter<double>("Muon2MetMtMaxCut");
  _DoDiscrByElectron1MetMt = iConfig.getParameter<bool>("DoDiscrByElectron1MetMt");
  _Electron1MetMtMinCut = iConfig.getParameter<double>("Electron1MetMtMinCut");
  _Electron1MetMtMaxCut = iConfig.getParameter<double>("Electron1MetMtMaxCut");
  _DoDiscrByElectron2MetMt = iConfig.getParameter<bool>("DoDiscrByElectron2MetMt");
  _Electron2MetMtMinCut = iConfig.getParameter<double>("Electron2MetMtMinCut");
  _Electron2MetMtMaxCut = iConfig.getParameter<double>("Electron2MetMtMaxCut");
  _DoDiscrByTau1MetMt = iConfig.getParameter<bool>("DoDiscrByTau1MetMt");
  _Tau1MetMtMinCut = iConfig.getParameter<double>("Tau1MetMtMinCut");
  _Tau1MetMtMaxCut = iConfig.getParameter<double>("Tau1MetMtMaxCut");
  _DoDiscrByTau2MetMt = iConfig.getParameter<bool>("DoDiscrByTau2MetMt");
  _Tau2MetMtMinCut = iConfig.getParameter<double>("Tau2MetMtMinCut");
  _Tau2MetMtMaxCut = iConfig.getParameter<double>("Tau2MetMtMaxCut");
  _DoMuon1DiscrByIsZllCut = iConfig.getParameter<bool>("DoMuon1DiscrByIsZllCut");
  _DoMuon2DiscrByIsZllCut = iConfig.getParameter<bool>("DoMuon2DiscrByIsZllCut");
  _DoElectron1DiscrByIsZllCut = iConfig.getParameter<bool>("DoElectron1DiscrByIsZllCut");
  _DoElectron2DiscrByIsZllCut = iConfig.getParameter<bool>("DoElectron2DiscrByIsZllCut");

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
//  _DoProduceNtuple = iConfig.getParameter<bool>("DoProduceNtuple");
//  _NtupleTreeName = (iConfig.getUntrackedParameter<std::string>("NtupleTreeName"));
  _DoProduceNtuple = false;

  //-----Event Sequence inputs
  _RecoTriggersNmin = iConfig.getParameter<int>("RecoTriggersNmin");
  _RecoVertexNmin = iConfig.getParameter<int>("RecoVertexNmin");
  _RecoVertexNmax = iConfig.getParameter<int>("RecoVertexNmax");
  _RecoMuon1Nmin = iConfig.getParameter<int>("RecoMuon1Nmin");
  _RecoMuon1Nmax = iConfig.getParameter<int>("RecoMuon1Nmax");
  _RecoMuon2Nmin = iConfig.getParameter<int>("RecoMuon2Nmin");
  _RecoMuon2Nmax = iConfig.getParameter<int>("RecoMuon2Nmax");
  _RecoElectron1Nmin = iConfig.getParameter<int>("RecoElectron1Nmin");
  _RecoElectron1Nmax = iConfig.getParameter<int>("RecoElectron1Nmax");
  _RecoElectron2Nmin = iConfig.getParameter<int>("RecoElectron2Nmin");
  _RecoElectron2Nmax = iConfig.getParameter<int>("RecoElectron2Nmax");
  _RecoTau1Nmin = iConfig.getParameter<int>("RecoTau1Nmin");
  _RecoTau1Nmax = iConfig.getParameter<int>("RecoTau1Nmax");
  _RecoTau2Nmin = iConfig.getParameter<int>("RecoTau2Nmin");
  _RecoTau2Nmax = iConfig.getParameter<int>("RecoTau2Nmax");
  _RecoJetNmin = iConfig.getParameter<int>("RecoJetNmin");
  _RecoJetNmax = iConfig.getParameter<int>("RecoJetNmax");
  _RecoFirstLeadingJetNmin = iConfig.getParameter<int>("RecoFirstLeadingJetNmin");
  _RecoSecondLeadingJetNmin = iConfig.getParameter<int>("RecoSecondLeadingJetNmin");
  _RecoBJetNmin = iConfig.getParameter<int>("RecoBJetNmin");
  _RecoBJetNmax = iConfig.getParameter<int>("RecoBJetNmax");
  _SusyCombinationsNmin = iConfig.getParameter<int>("SusyCombinationsNmin");
  _DiMuonCombinationsNmin = iConfig.getParameter<int>("DiMuonCombinationsNmin");
  _DiMuonCombinationsNmax = iConfig.getParameter<int>("DiMuonCombinationsNmax");
  _DiElectronCombinationsNmin = iConfig.getParameter<int>("DiElectronCombinationsNmin");
  _DiElectronCombinationsNmax = iConfig.getParameter<int>("DiElectronCombinationsNmax");
  _DiTauCombinationsNmin = iConfig.getParameter<int>("DiTauCombinationsNmin");
  _DiTauCombinationsNmax = iConfig.getParameter<int>("DiTauCombinationsNmax");
  _Muon1Tau1CombinationsNmin = iConfig.getParameter<int>("Muon1Tau1CombinationsNmin");
  _Muon1Tau1CombinationsNmax = iConfig.getParameter<int>("Muon1Tau1CombinationsNmax");
  _Muon1Tau2CombinationsNmin = iConfig.getParameter<int>("Muon1Tau2CombinationsNmin");
  _Muon1Tau2CombinationsNmax = iConfig.getParameter<int>("Muon1Tau2CombinationsNmax");
  _Muon2Tau1CombinationsNmin = iConfig.getParameter<int>("Muon2Tau1CombinationsNmin");
  _Muon2Tau1CombinationsNmax = iConfig.getParameter<int>("Muon2Tau1CombinationsNmax");
  _Muon2Tau2CombinationsNmin = iConfig.getParameter<int>("Muon2Tau2CombinationsNmin");
  _Muon2Tau2CombinationsNmax = iConfig.getParameter<int>("Muon2Tau2CombinationsNmax");
  _Electron1Tau1CombinationsNmin = iConfig.getParameter<int>("Electron1Tau1CombinationsNmin");
  _Electron1Tau1CombinationsNmax = iConfig.getParameter<int>("Electron1Tau1CombinationsNmax");
  _Electron1Tau2CombinationsNmin = iConfig.getParameter<int>("Electron1Tau2CombinationsNmin");
  _Electron1Tau2CombinationsNmax = iConfig.getParameter<int>("Electron1Tau2CombinationsNmax");
  _Electron2Tau1CombinationsNmin = iConfig.getParameter<int>("Electron2Tau1CombinationsNmin");
  _Electron2Tau1CombinationsNmax = iConfig.getParameter<int>("Electron2Tau1CombinationsNmax");
  _Electron2Tau2CombinationsNmin = iConfig.getParameter<int>("Electron2Tau2CombinationsNmin");
  _Electron2Tau2CombinationsNmax = iConfig.getParameter<int>("Electron2Tau2CombinationsNmax");
  _Electron1Muon1CombinationsNmin = iConfig.getParameter<int>("Electron1Muon1CombinationsNmin");
  _Electron1Muon1CombinationsNmax = iConfig.getParameter<int>("Electron1Muon1CombinationsNmax");
  _Electron1Muon2CombinationsNmin = iConfig.getParameter<int>("Electron1Muon2CombinationsNmin");
  _Electron1Muon2CombinationsNmax = iConfig.getParameter<int>("Electron1Muon2CombinationsNmax");
  _Electron2Muon1CombinationsNmin = iConfig.getParameter<int>("Electron2Muon1CombinationsNmin");
  _Electron2Muon1CombinationsNmax = iConfig.getParameter<int>("Electron2Muon1CombinationsNmax");
  _Electron2Muon2CombinationsNmin = iConfig.getParameter<int>("Electron2Muon2CombinationsNmin");
  _Electron2Muon2CombinationsNmax = iConfig.getParameter<int>("Electron2Muon2CombinationsNmax");
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
  _UseDataBaseForJEC = iConfig.getParameter<bool>("UseDataBaseForJEC");
  JES_UpOrDown = iConfig.getParameter<double>("JES_UpOrDown");
  _JetEnergyScaleOffset = iConfig.getParameter<double>("JetEnergyScaleOffset");
  _SmearThePt = iConfig.getParameter<bool>("SmearThePt");
  _SmearTheEta = iConfig.getParameter<bool>("SmearTheEta");
  _SmearThePhi = iConfig.getParameter<bool>("SmearThePhi");
  _CalculatePUSystematics = iConfig.getParameter<bool>("CalculatePUSystematics");
  _DataHistos = iConfig.getParameter<string>("DataHistos");
  _MCHistos = iConfig.getParameter<string>("MCHistos");
//  _ApplyMuonTriggerScaleFactors = iConfig.getParameter<bool>("ApplyMuonTriggerScaleFactors");
//  _ApplyElectronTriggerScaleFactors = iConfig.getParameter<bool>("ApplyElectronTriggerScaleFactors");
//  _ApplyTauTriggerScaleFactors = iConfig.getParameter<bool>("ApplyTauTriggerScaleFactors");

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

/*
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
*/

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
/*
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
*/
}

// Initialize the vectors for the ntuple
void HiMassTauAnalysis::initializeVectors(){
/*
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
*/
}

// clear the vectors vefore each event
void HiMassTauAnalysis::clearVectors(){
/*
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
*/
}

// ------------ method called to for each event  ------------
void HiMassTauAnalysis::analyze(const Event& iEvent, const EventSetup& iSetup) {

if(_DoMSUGRApoint)
{
  boost::smatch matches;
  Handle<LHEEventProduct> lhe;
  iEvent.getByLabel("source", lhe);

  if(lhe.isValid()) {
    const std::vector<std::string> comments(lhe->comments_begin(),lhe->comments_end());
    BOOST_FOREACH(const std::string& comment, comments)
      if (boost::regex_match(comment, matches, scanFormat)) break;
  }

/*
  float lspmass=0;
  float staumass=0;
  for(GenParticleCollection::const_iterator genParticle = _genParticles->begin();genParticle != _genParticles->end();++genParticle) {
    if((abs(genParticle->pdgId()) == 1000022)) {lspmass = genParticle->mass();}
    if((abs(genParticle->pdgId()) == 1000015)) {staumass = genParticle->mass();}
  }
//  if( ((staumass-lspmass) < 10.0) || ((staumass-lspmass) > 25.0) ) {return;}
*/

  double m0 = 0;
  double m12 = 0;
  for(unsigned i=0; i<scanPars.size(); ++i) {
    double xxx= atof((matches[i+1].str()).c_str());
//    if(i == 0) {std::cout << "m0=" << xxx << ", ";}
//    if(i == 1) {std::cout << "m12=" << xxx << ", ";}
//    if(i == 4) {std::cout << "signMu=" << xxx << ", dM=" << staumass-lspmass << std::endl;}
    if(i == 0) {m0 = xxx;}
    if(i == 1) {m12 = xxx;}
  }

  if(_SelectSusyScanPoint) {if((m0 != _M0) || (m12 != _M12)) {return;}}
}

if(_DoSMpoint){
  Handle<LHEEventProduct> product;
  iEvent.getByLabel("source", product);

//  comments_const_iterator c_begin = product->comments_begin();
//  comments_const_iterator c_end = product->comments_end();

  std::vector<std::string>::const_iterator c_begin = product->comments_begin();
  std::vector<std::string>::const_iterator c_end = product->comments_end();

  double mGL;
  double mLSP;
  double xCHI;

  for( std::vector<std::string>::const_iterator cit=c_begin; cit!=c_end; ++cit) {
    size_t found = (*cit).find("model");

    //# model T5zz_0.5_925.0_400.0to1000.0_450.0.lhe
    if( found != std::string::npos)   {
//         std::cout << *cit << std::endl;
         size_t foundLength = (*cit).size();
         found = (*cit).find("T5zz");
         std::string smaller = (*cit).substr(found+1,foundLength);
         found = smaller.find("_");
         smaller = smaller.substr(found+1,smaller.size());
//
         std::istringstream iss(smaller);
         iss >> xCHI;
         iss.clear();
//
         found = smaller.find("_");
         smaller = smaller.substr(found+1,smaller.size());
         iss.str(smaller);
         iss >> mGL;
         iss.clear();
//
         found = smaller.find("_");
         smaller = smaller.substr(found+1,smaller.size());
         iss.str(smaller);
         iss >> mLSP;
         iss.clear();
      
      }
   
  }

//  std::cout << "mGL: " << mGL << " mLSP: " << mLSP << std::endl;
//  std::cout << "_mGL: " << _mGL << " _mLSP: " << _mLSP << std::endl;
//  if(_mLSP == mLSP && _mGL == mGL)std::cout << "Victory!" << std::endl;
  if(_mLSP != mLSP || _mGL != mGL){return;}

}

  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl); 
  std::vector<JetCorrectionUncertainty*> jecUnc;
  jecUnc.clear();
  jecUnc.push_back(new JetCorrectionUncertainty((*JetCorParColl)["Uncertainty"]));

  //------Number of events analyzed (denominator)
  _totalEvents++;

  // Running over real data or MC?
  if(_totalEvents == 1) {isData = iEvent.isRealData();}

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

/*
   Handle<LHEEventProduct> product;
   iEvent.getByLabel("source", product);

   comments_const_iterator c_begin = product->comments_begin();
   comments_const_iterator c_end = product->comments_end();

   double m0, m12, tanb, A0, mu=1.0;
   double signMu;
   for( comments_const_iterator cit=c_begin; cit!=c_end; ++cit) {
      size_t found = (*cit).find("model");
      if( found != std::string::npos)   {    
//         std::cout << *cit << std::endl;  
         size_t foundLength = (*cit).size();
         found = (*cit).find("=");
         std::string smaller = (*cit).substr(found+1,foundLength);
         found = smaller.find("_");
         smaller = smaller.substr(found+1,smaller.size());
//
         std::istringstream iss(smaller);
         iss >> m0;
         iss.clear();
//
         found = smaller.find("_");
         smaller = smaller.substr(found+1,smaller.size());
         iss.str(smaller);
         iss >> m12;
         iss.clear();
//
         found = smaller.find("_");
         smaller = smaller.substr(found+1,smaller.size());
         iss.str(smaller);
         iss >> tanb;
         iss.clear();
      }
   }

   std::cout << "m0 = " << m0 << std::endl;
   std::cout << "m12 = " << m12 << std::endl;
   std::cout << "tanb = " << tanb << std::endl;
*/

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

/*
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
*/

  //This is where pileup weight is calculated **** Will Flanagan's code ****
  if((_CalculatePUSystematics) && (!isData)) {

    int nVertices = 0;
    for(reco::VertexCollection::const_iterator primaryVertex = _primaryEventVertexCollection->begin(); primaryVertex != _primaryEventVertexCollection->end(); ++primaryVertex ) {
      if (!passRecoVertexCuts(*primaryVertex)) continue;
      nVertices++;
    }

    int bin;
    double MCintegral; double MCvalue; double Dataintegral; double Datavalue;
    char * cstr1;
    char * cstr2;

    //Filenames must be c_strings below. Here is the conversion from strings
    cstr1 = new char [_MCHistos.size()+1];
    strcpy (cstr1, _MCHistos.c_str());
    cstr2 = new char [_DataHistos.size()+1];
    strcpy (cstr2, _DataHistos.c_str());

    //As you can see above cstr1 corresponds to MC and cstr2 corresponds to data.
    //Probability that file has N vertices is value / integral
    //The ratio of data probability to MC probability gives us our PU weight
    bin = getBin(cstr1,nVertices);
    MCvalue = getvalue(cstr1,bin);
    MCintegral = getintegral(cstr1,bin);

    Datavalue = getvalue(cstr2,bin);
    Dataintegral = getintegral(cstr2,bin);

    //Ratio of normalized histograms in given bin
    if((MCvalue * Dataintegral) != 0) {pu_weight = (Datavalue * MCintegral) / (MCvalue * Dataintegral);}
    else {pu_weight = 1.0;}

    //cout << "pu_weight: " << pu_weight << " nVertices: " << nVertices << std::endl;
    //cout << "MCvalue: " << MCvalue << " MCintegral: " << MCintegral << " Datavalue: " << Datavalue << " Dataintegral: " << Dataintegral << std::endl;
    if (isinf(pu_weight) || isnan(pu_weight)) {throw std::runtime_error("PU weight INF or NAN");}
  } else {pu_weight = 1.0;}

  //-----Smearing momentum and position for systematic uncertanties and calculation of MET deltas
  smearedMuonMomentumVector.clear();
  smearedMuonPtEtaPhiMVector.clear();
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
  smearedElectronMomentumVector.clear();
  smearedElectronPtEtaPhiMVector.clear();
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
  smearedTauMomentumVector.clear();
  smearedTauPtEtaPhiMVector.clear();
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
  smearedJetMomentumVector.clear();
  smearedJetPtEtaPhiMVector.clear();
  if(_SmearTheJet) {
    for ( pat::JetCollection::const_iterator patJet = _patJets->begin(); patJet != _patJets->end(); ++patJet ) {
      smearedJetMomentumVector.push_back(SmearJet((*patJet),jecUnc[0]).first);
      smearedJetPtEtaPhiMVector.push_back(SmearJet((*patJet),jecUnc[0]).second);
      deltaForMEx = deltaForMEx + patJet->px() - SmearJet((*patJet),jecUnc[0]).first.px();
      deltaForMEy = deltaForMEy + patJet->py() - SmearJet((*patJet),jecUnc[0]).first.py();
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
    _hEvents[NpdfID]->Fill(0.0);
  }

  //------Get the event flags (did the event pass the cuts?)
  getEventFlags(iEvent);

  if (!passEventSelectionSequence()) return;  

  isrgluon_weight = isrgluon_weight * trig_weight * pu_weight;

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

  sumpxForMht = 0.0;
  sumpyForMht = 0.0;
  sumptForHt  = 0.0;

/*
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
*/

  //------Number of Good Candidates
  int nGoodCandidatesMuon1 = 0;
  int nGoodCandidatesElectron1 = 0;
  int nGoodCandidatesTau1 = 0;
  int nGoodCandidatesMuon2 = 0;
  int nGoodCandidatesElectron2 = 0;
  int nGoodCandidatesTau2 = 0;
  int theNumberOfMuons = 0;
  for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
    theNumberOfMuons++;
    if (!passRecoMuon1Cuts((*patMuon),theNumberOfMuons-1)) continue;
    nGoodCandidatesMuon1++;
  }
  theNumberOfMuons = 0;
  for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
    theNumberOfMuons++;
    if (!passRecoMuon2Cuts((*patMuon),theNumberOfMuons-1)) continue;
    nGoodCandidatesMuon2++;
  }
  int theNumberOfElectrons = 0;
  for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron ) {
    theNumberOfElectrons++;
    if (!passRecoElectron1Cuts((*patElectron),theNumberOfElectrons-1)) continue;
    nGoodCandidatesElectron1++;
  }
  theNumberOfElectrons = 0;
  for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron ) {
    theNumberOfElectrons++;
    if (!passRecoElectron2Cuts((*patElectron),theNumberOfElectrons-1)) continue;
    nGoodCandidatesElectron2++;
  }
  int theNumberOfTaus = 0;
  for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
    theNumberOfTaus++;
    if (!passRecoTau1Cuts((*patTau),theNumberOfTaus-1)) continue;
    nGoodCandidatesTau1++;
  }
  theNumberOfTaus = 0;
  for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
    theNumberOfTaus++;
    if (!passRecoTau2Cuts((*patTau),theNumberOfTaus-1)) continue;
    nGoodCandidatesTau2++;
  }
  if (nGoodCandidatesMuon1>=_RecoMuon1Nmin) _EventFlag[_mapSelectionAlgoID["RecoMuon1Nmin"]] = true;
  if (nGoodCandidatesMuon1<=_RecoMuon1Nmax) _EventFlag[_mapSelectionAlgoID["RecoMuon1Nmax"]] = true;
  if (nGoodCandidatesMuon2>=_RecoMuon2Nmin) _EventFlag[_mapSelectionAlgoID["RecoMuon2Nmin"]] = true;
  if (nGoodCandidatesMuon2<=_RecoMuon2Nmax) _EventFlag[_mapSelectionAlgoID["RecoMuon2Nmax"]] = true;
  if (nGoodCandidatesElectron1>=_RecoElectron1Nmin) _EventFlag[_mapSelectionAlgoID["RecoElectron1Nmin"]] = true;
  if (nGoodCandidatesElectron1<=_RecoElectron1Nmax) _EventFlag[_mapSelectionAlgoID["RecoElectron1Nmax"]] = true;
  if (nGoodCandidatesElectron2>=_RecoElectron2Nmin) _EventFlag[_mapSelectionAlgoID["RecoElectron2Nmin"]] = true;
  if (nGoodCandidatesElectron2<=_RecoElectron2Nmax) _EventFlag[_mapSelectionAlgoID["RecoElectron2Nmax"]] = true;
  if (nGoodCandidatesTau1>=_RecoTau1Nmin) _EventFlag[_mapSelectionAlgoID["RecoTau1Nmin"]] = true;
  if (nGoodCandidatesTau1<=_RecoTau1Nmax) _EventFlag[_mapSelectionAlgoID["RecoTau1Nmax"]] = true;
  if (nGoodCandidatesTau2>=_RecoTau2Nmin) _EventFlag[_mapSelectionAlgoID["RecoTau2Nmin"]] = true;
  if (nGoodCandidatesTau2<=_RecoTau2Nmax) _EventFlag[_mapSelectionAlgoID["RecoTau2Nmax"]] = true;

  // ------Number of Good Jets   
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

  // ------Number of Good b-Jets
  int nGoodBJets = 0;
  theNumberOfJets = 0;
  for ( pat::JetCollection::const_iterator patJet = _patJets->begin();patJet != _patJets->end(); ++patJet ) {
    theNumberOfJets++;
    if (!passRecoBJetCuts((*patJet),theNumberOfJets-1)) continue;
    nGoodBJets++;
  }
  if (nGoodBJets>=_RecoBJetNmin) _EventFlag[_mapSelectionAlgoID["RecoBJetNmin"]] = true;
  if (nGoodBJets<=_RecoBJetNmax) _EventFlag[_mapSelectionAlgoID["RecoBJetNmax"]] = true;

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

  // ------Number of Good Combinations
  int nGoodMuon1Tau1Combinations = 0;
  theNumberOfMuons = 0;
  for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
    theNumberOfMuons++;
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if ((passRecoTau1Cuts((*patTau),theNumberOfTaus-1)) &&
          (passRecoMuon1Cuts((*patMuon),theNumberOfMuons-1)) && 
          (passMuon1Tau1TopologyCuts((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1))) {
        nGoodMuon1Tau1Combinations++;
      }
    }
  }
  int nGoodMuon1Tau2Combinations = 0;
  theNumberOfMuons = 0;
  for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
    theNumberOfMuons++;
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if ((passRecoTau2Cuts((*patTau),theNumberOfTaus-1)) &&
          (passRecoMuon1Cuts((*patMuon),theNumberOfMuons-1)) && 
          (passMuon1Tau2TopologyCuts((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1))) {
        nGoodMuon1Tau2Combinations++;
      }
    }
  }
  int nGoodMuon2Tau1Combinations = 0;
  theNumberOfMuons = 0;
  for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
    theNumberOfMuons++;
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if ((passRecoTau1Cuts((*patTau),theNumberOfTaus-1)) &&
          (passRecoMuon2Cuts((*patMuon),theNumberOfMuons-1)) && 
          (passMuon2Tau1TopologyCuts((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1))) {
        nGoodMuon2Tau1Combinations++;
      }
    }
  }
  int nGoodMuon2Tau2Combinations = 0;
  theNumberOfMuons = 0;
  for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
    theNumberOfMuons++;
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if ((passRecoTau2Cuts((*patTau),theNumberOfTaus-1)) &&
          (passRecoMuon2Cuts((*patMuon),theNumberOfMuons-1)) && 
          (passMuon2Tau2TopologyCuts((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1))) {
        nGoodMuon2Tau2Combinations++;
      }
    }
  }
  int nGoodElectron1Tau1Combinations = 0;
  theNumberOfElectrons = 0;
  for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
    theNumberOfElectrons++;
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if ((passRecoTau1Cuts((*patTau),theNumberOfTaus-1)) &&
          (passRecoElectron1Cuts((*patElectron),theNumberOfElectrons-1)) &&
          (passElectron1Tau1TopologyCuts((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1))) {
        nGoodElectron1Tau1Combinations++;
      }
    }
  }
  int nGoodElectron2Tau1Combinations = 0;
  theNumberOfElectrons = 0;
  for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
    theNumberOfElectrons++;
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if ((passRecoTau1Cuts((*patTau),theNumberOfTaus-1)) &&
          (passRecoElectron2Cuts((*patElectron),theNumberOfElectrons-1)) &&
          (passElectron2Tau1TopologyCuts((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1))) {
        nGoodElectron2Tau1Combinations++;
      }
    }
  }
  int nGoodElectron1Tau2Combinations = 0;
  theNumberOfElectrons = 0;
  for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
    theNumberOfElectrons++;
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if ((passRecoTau2Cuts((*patTau),theNumberOfTaus-1)) &&
          (passRecoElectron1Cuts((*patElectron),theNumberOfElectrons-1)) &&
          (passElectron1Tau2TopologyCuts((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1))) {
        nGoodElectron1Tau2Combinations++;
      }
    }
  }
  int nGoodElectron2Tau2Combinations = 0;
  theNumberOfElectrons = 0;
  for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
    theNumberOfElectrons++;
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if ((passRecoTau2Cuts((*patTau),theNumberOfTaus-1)) &&
          (passRecoElectron2Cuts((*patElectron),theNumberOfElectrons-1)) &&
          (passElectron2Tau2TopologyCuts((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1))) {
        nGoodElectron2Tau2Combinations++;
      }
    }
  }
  int nGoodElectron1Muon1Combinations = 0;
  theNumberOfMuons = 0;
  for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
    theNumberOfMuons++;
    int theNumberOfElectrons = 0;
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
      theNumberOfElectrons++;
      if ((passRecoElectron1Cuts((*patElectron),theNumberOfElectrons-1)) && 
          (passRecoMuon1Cuts((*patMuon),theNumberOfMuons-1)) && 
          (passElectron1Muon1TopologyCuts((*patElectron),theNumberOfElectrons-1,(*patMuon),theNumberOfMuons-1))) {
        nGoodElectron1Muon1Combinations++;
      }
    }
  }
  int nGoodElectron2Muon1Combinations = 0;
  theNumberOfMuons = 0;
  for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
    theNumberOfMuons++;
    int theNumberOfElectrons = 0;
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
      theNumberOfElectrons++;
      if ((passRecoElectron2Cuts((*patElectron),theNumberOfElectrons-1)) && 
          (passRecoMuon1Cuts((*patMuon),theNumberOfMuons-1)) && 
          (passElectron2Muon1TopologyCuts((*patElectron),theNumberOfElectrons-1,(*patMuon),theNumberOfMuons-1))) {
        nGoodElectron2Muon1Combinations++;
      }
    }
  }
  int nGoodElectron1Muon2Combinations = 0;
  theNumberOfMuons = 0;
  for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
    theNumberOfMuons++;
    int theNumberOfElectrons = 0;
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
      theNumberOfElectrons++;
      if ((passRecoElectron1Cuts((*patElectron),theNumberOfElectrons-1)) && 
          (passRecoMuon2Cuts((*patMuon),theNumberOfMuons-1)) && 
          (passElectron2Muon1TopologyCuts((*patElectron),theNumberOfElectrons-1,(*patMuon),theNumberOfMuons-1))) {
        nGoodElectron1Muon2Combinations++;
      }
    }
  }
  int nGoodElectron2Muon2Combinations = 0;
  theNumberOfMuons = 0;
  for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
    theNumberOfMuons++;
    int theNumberOfElectrons = 0;
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
      theNumberOfElectrons++;
      if ((passRecoElectron2Cuts((*patElectron),theNumberOfElectrons-1)) && 
          (passRecoMuon2Cuts((*patMuon),theNumberOfMuons-1)) && 
          (passElectron2Muon2TopologyCuts((*patElectron),theNumberOfElectrons-1,(*patMuon),theNumberOfMuons-1))) {
        nGoodElectron2Muon2Combinations++;
      }
    }
  }
  int nGoodDiTauCombinations = 0;
  int theNumberOfTaus1 = 0;
  for ( pat::TauCollection::const_iterator patTau1 = _patTaus->begin();patTau1 != _patTaus->end(); ++patTau1 ) {
    theNumberOfTaus1++;
    int theNumberOfTaus2 = 0;
    for ( pat::TauCollection::const_iterator patTau2 = _patTaus->begin();patTau2 != _patTaus->end(); ++patTau2 ) {
      theNumberOfTaus2++;
      if ((passRecoTau1Cuts((*patTau1),theNumberOfTaus1 - 1)) && 
          (passRecoTau2Cuts((*patTau2),theNumberOfTaus2 - 1)) && 
          (passDiTauTopologyCuts((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1)) ) {
        nGoodDiTauCombinations++;
      }
    }
  }
  int nGoodDiMuonCombinations = 0;
  int theNumberOfMuons1 = 0;
  for ( pat::MuonCollection::const_iterator patMuon1 = _patMuons->begin();patMuon1 != _patMuons->end(); ++patMuon1 ) {
    theNumberOfMuons1++;
    int theNumberOfMuons2 = 0;
    for ( pat::MuonCollection::const_iterator patMuon2 = _patMuons->begin();patMuon2 != _patMuons->end(); ++patMuon2 ) {
      theNumberOfMuons2++;
      if ((passRecoMuon1Cuts((*patMuon1),theNumberOfMuons1 - 1)) && 
          (passRecoMuon2Cuts((*patMuon2),theNumberOfMuons2 - 1)) && 
          (passDiMuonTopologyCuts((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1)) ) {
        nGoodDiMuonCombinations++;
      }
    }
  }
  int nGoodDiElectronCombinations = 0;
  int theNumberOfElectrons1 = 0;
  for ( pat::ElectronCollection::const_iterator patElectron1 = _patElectrons->begin();patElectron1 != _patElectrons->end(); ++patElectron1 ) {
    theNumberOfElectrons1++;
    int theNumberOfElectrons2 = 0;
    for ( pat::ElectronCollection::const_iterator patElectron2 = _patElectrons->begin();patElectron2 != _patElectrons->end(); ++patElectron2 ) {
      theNumberOfElectrons2++;
      if ((passRecoElectron1Cuts((*patElectron1),theNumberOfElectrons1 - 1)) && 
          (passRecoElectron2Cuts((*patElectron2),theNumberOfElectrons2 - 1)) && 
          (passDiElectronTopologyCuts((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1)) ) {
        nGoodDiElectronCombinations++;
      }
    }
  }
  if (nGoodDiMuonCombinations>=_DiMuonCombinationsNmin) _EventFlag[_mapSelectionAlgoID["DiMuonCombinationsNmin"]] = true;
  if (nGoodDiMuonCombinations<=_DiMuonCombinationsNmax) _EventFlag[_mapSelectionAlgoID["DiMuonCombinationsNmax"]] = true;
  if (nGoodDiElectronCombinations>=_DiElectronCombinationsNmin) _EventFlag[_mapSelectionAlgoID["DiElectronCombinationsNmin"]] = true;
  if (nGoodDiElectronCombinations<=_DiElectronCombinationsNmax) _EventFlag[_mapSelectionAlgoID["DiElectronCombinationsNmax"]] = true;
  if (nGoodDiTauCombinations>=_DiTauCombinationsNmin) _EventFlag[_mapSelectionAlgoID["DiTauCombinationsNmin"]] = true;
  if (nGoodDiTauCombinations<=_DiTauCombinationsNmax) _EventFlag[_mapSelectionAlgoID["DiTauCombinationsNmax"]] = true;
  if (nGoodMuon1Tau1Combinations>=_Muon1Tau1CombinationsNmin) _EventFlag[_mapSelectionAlgoID["Muon1Tau1CombinationsNmin"]] = true;
  if (nGoodMuon1Tau1Combinations<=_Muon1Tau1CombinationsNmax) _EventFlag[_mapSelectionAlgoID["Muon1Tau1CombinationsNmax"]] = true;
  if (nGoodMuon1Tau2Combinations>=_Muon1Tau2CombinationsNmin) _EventFlag[_mapSelectionAlgoID["Muon1Tau2CombinationsNmin"]] = true;
  if (nGoodMuon1Tau2Combinations<=_Muon1Tau2CombinationsNmax) _EventFlag[_mapSelectionAlgoID["Muon1Tau2CombinationsNmax"]] = true;
  if (nGoodMuon2Tau1Combinations>=_Muon2Tau1CombinationsNmin) _EventFlag[_mapSelectionAlgoID["Muon2Tau1CombinationsNmin"]] = true;
  if (nGoodMuon2Tau1Combinations<=_Muon2Tau1CombinationsNmax) _EventFlag[_mapSelectionAlgoID["Muon2Tau1CombinationsNmax"]] = true;
  if (nGoodMuon2Tau2Combinations>=_Muon2Tau2CombinationsNmin) _EventFlag[_mapSelectionAlgoID["Muon2Tau2CombinationsNmin"]] = true;
  if (nGoodMuon2Tau2Combinations<=_Muon2Tau2CombinationsNmax) _EventFlag[_mapSelectionAlgoID["Muon2Tau2CombinationsNmax"]] = true;
  if (nGoodElectron1Tau1Combinations>=_Electron1Tau1CombinationsNmin) _EventFlag[_mapSelectionAlgoID["Electron1Tau1CombinationsNmin"]] = true;
  if (nGoodElectron1Tau1Combinations<=_Electron1Tau1CombinationsNmax) _EventFlag[_mapSelectionAlgoID["Electron1Tau1CombinationsNmax"]] = true;
  if (nGoodElectron1Tau2Combinations>=_Electron1Tau2CombinationsNmin) _EventFlag[_mapSelectionAlgoID["Electron1Tau2CombinationsNmin"]] = true;
  if (nGoodElectron1Tau2Combinations<=_Electron1Tau2CombinationsNmax) _EventFlag[_mapSelectionAlgoID["Electron1Tau2CombinationsNmax"]] = true;
  if (nGoodElectron2Tau1Combinations>=_Electron2Tau1CombinationsNmin) _EventFlag[_mapSelectionAlgoID["Electron2Tau1CombinationsNmin"]] = true;
  if (nGoodElectron2Tau1Combinations<=_Electron2Tau1CombinationsNmax) _EventFlag[_mapSelectionAlgoID["Electron2Tau1CombinationsNmax"]] = true;
  if (nGoodElectron2Tau2Combinations>=_Electron2Tau2CombinationsNmin) _EventFlag[_mapSelectionAlgoID["Electron2Tau2CombinationsNmin"]] = true;
  if (nGoodElectron2Tau2Combinations<=_Electron2Tau2CombinationsNmax) _EventFlag[_mapSelectionAlgoID["Electron2Tau2CombinationsNmax"]] = true;
  if (nGoodElectron1Muon1Combinations>=_Electron1Muon1CombinationsNmin) _EventFlag[_mapSelectionAlgoID["Electron1Muon1CombinationsNmin"]] = true;
  if (nGoodElectron1Muon1Combinations<=_Electron1Muon1CombinationsNmax) _EventFlag[_mapSelectionAlgoID["Electron1Muon1CombinationsNmax"]] = true;
  if (nGoodElectron1Muon2Combinations>=_Electron1Muon2CombinationsNmin) _EventFlag[_mapSelectionAlgoID["Electron1Muon2CombinationsNmin"]] = true;
  if (nGoodElectron1Muon2Combinations<=_Electron1Muon2CombinationsNmax) _EventFlag[_mapSelectionAlgoID["Electron1Muon2CombinationsNmax"]] = true;
  if (nGoodElectron2Muon1Combinations>=_Electron2Muon1CombinationsNmin) _EventFlag[_mapSelectionAlgoID["Electron2Muon1CombinationsNmin"]] = true;
  if (nGoodElectron2Muon1Combinations<=_Electron2Muon1CombinationsNmax) _EventFlag[_mapSelectionAlgoID["Electron2Muon1CombinationsNmax"]] = true;
  if (nGoodElectron2Muon2Combinations>=_Electron2Muon2CombinationsNmin) _EventFlag[_mapSelectionAlgoID["Electron2Muon2CombinationsNmin"]] = true;
  if (nGoodElectron2Muon2Combinations<=_Electron2Muon2CombinationsNmax) _EventFlag[_mapSelectionAlgoID["Electron2Muon2CombinationsNmax"]] = true;

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
//      std::cout << "### HiMassTauAnalysis - CONFIGURATION ERROR:  Specified trigger " << (*TheTriggerPath) << " is not found/defined!!!!" << std::endl;
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
bool HiMassTauAnalysis::passRecoTau1Cuts(const pat::Tau& patTau,int nobj) {
  // ----Matching to gen
  if((_MatchTauToGen) && (!isData)) {
    if(_GenParticleSource.label() != "") {
      if(!(matchToGen(patTau).first)) {return false;}
    } else {return false;}
  }
  // ----Acceptance cuts
  if (fabs(smearedTauPtEtaPhiMVector.at(nobj).eta())>_RecoTau1EtaCut) {return false;}
  if (smearedTauPtEtaPhiMVector.at(nobj).pt()<_RecoTau1PtMinCut) {return false;}
  if (smearedTauPtEtaPhiMVector.at(nobj).pt()>_RecoTau1PtMaxCut) {return false;}
  // ----Lead track requirement
  if (_DoRecoTau1DiscrByLeadTrack) {
    if (_UseRecoTau1DiscrByLeadTrackFlag) { if ( (patTau.tauID(_RecoTau1DiscrByLeadTrack.data()) < 0.5) ) {return false;} }
    else {
      if (patTau.isCaloTau()) { if( (!(patTau.leadTrack().isNonnull())) || (patTau.leadTrack()->pt() < _RecoTau1LeadTrackThreshold) ) {return false;} }
      else { if( (!(patTau.leadPFChargedHadrCand().isNonnull())) || (patTau.leadPFChargedHadrCand()->pt() < _RecoTau1LeadTrackThreshold) ) {return false;} }
    }
  }
  // ----"eVeto" - Lead track minimimum hits requirement && H3x3/LTpt
  if (_DoRecoTau1DiscrByLeadTrackNhits) {
    if (patTau.isCaloTau()) {
      if( (!(patTau.leadTrack().isNonnull())) || ((int)(patTau.leadTrack()->numberOfValidHits()) < _RecoTau1LeadTrackMinHits) ) {return false;}
    } else {
      if( (!(patTau.leadPFChargedHadrCand().isNonnull())) || (!(patTau.leadPFChargedHadrCand()->trackRef().isNonnull())) ) {return false;}
      if( (int)(patTau.leadPFChargedHadrCand()->trackRef()->numberOfValidHits()) < _RecoTau1LeadTrackMinHits ) {return false;}
    }
  }
  if (_DoRecoTau1DiscrByH3x3OverP) {
    if (patTau.isCaloTau()) {
      if( (!(patTau.leadTrack().isNonnull())) || ((patTau.hcal3x3OverPLead()) <= _RecoTau1H3x3OverP) ) {return false;}
    } else {
      if( (!(patTau.leadPFChargedHadrCand().isNonnull())) || (patTau.hcal3x3OverPLead() <= _RecoTau1H3x3OverP) ) {return false;}
    }
  }
  // ----Isolation requirement
  if (_DoRecoTau1DiscrByIsolation) {
    if (_UseRecoTau1DiscrByIsolationFlag) {if ( (patTau.tauID(_RecoTau1DiscrByIsolation.data()) < 0.5) ) {return false;}}
    else {
      if(_UseRecoTau1IsoSumPtInsteadOfNiso) {
        _RecoTauIsoDeltaRCone = _RecoTau1IsoDeltaRCone;
        _RecoTauTrackIsoTrkThreshold = _RecoTau1TrackIsoTrkThreshold;
        _RecoTauGammaIsoGamThreshold = _RecoTau1GammaIsoGamThreshold;
        _UseRecoTauEllipseForEcalIso = _UseRecoTau1EllipseForEcalIso;
        _RecoTauEcalIsoRphiForEllipse = _RecoTau1EcalIsoRphiForEllipse;
        _RecoTauEcalIsoRetaForEllipse = _RecoTau1EcalIsoRetaForEllipse;
        if( (CalculateTauTrackIsolation(patTau).second) >= _RecoTau1TrackIsoSumPtMaxCutValue ) {return false;}
        if( (CalculateTauTrackIsolation(patTau).second) < _RecoTau1TrackIsoSumPtMinCutValue ) {return false;}
        if( (CalculateTauEcalIsolation(patTau).second) >= _RecoTau1EcalIsoSumPtMaxCutValue ) {return false;}
        if( (CalculateTauEcalIsolation(patTau).second) < _RecoTau1EcalIsoSumPtMinCutValue ) {return false;}
      } else {
        _RecoTauIsoDeltaRCone = _RecoTau1IsoDeltaRCone;
        _RecoTauTrackIsoTrkThreshold = _RecoTau1TrackIsoTrkThreshold;
        _RecoTauGammaIsoGamThreshold = _RecoTau1GammaIsoGamThreshold;
        _UseRecoTauEllipseForEcalIso = _UseRecoTau1EllipseForEcalIso;
        _RecoTauEcalIsoRphiForEllipse = _RecoTau1EcalIsoRphiForEllipse;
        _RecoTauEcalIsoRetaForEllipse = _RecoTau1EcalIsoRetaForEllipse;
        if( (CalculateTauTrackIsolation(patTau).first) > _RecoTau1TrackNisoMax ) {return false;}
        if( (CalculateTauEcalIsolation(patTau).first) > _RecoTau1EcalNisoMax ) {return false;}
      }
    }
  }
  // ----Require 1 or 3 prongs
  if (_RecoTau1DiscrByProngType == "1or3") {
    if (patTau.isCaloTau()) {
      if((patTau.signalTracks().size() == 1) ||(patTau.signalTracks().size() == 3)) {}
      else {return false;}
    } else {
      if((patTau.signalPFChargedHadrCands().size() == 1) ||(patTau.signalPFChargedHadrCands().size() == 3)) {}
      else {return false;}
    }
  } else if (_RecoTau1DiscrByProngType == "1") {
    if (patTau.isCaloTau()) {
      if(patTau.signalTracks().size() == 1) {}
      else {return false;}
    } else {
      if(patTau.signalPFChargedHadrCands().size() == 1) {}
      else {return false;}
    }
  } else if (_RecoTau1DiscrByProngType == "3") {
    if (patTau.isCaloTau()) {
      if(patTau.signalTracks().size() == 3) {}
      else {return false;}
    } else {
      if(patTau.signalPFChargedHadrCands().size() == 3) {}
      else {return false;}
    }
  } else if (_RecoTau1DiscrByProngType == "1or3hps") {
    if ( (patTau.tauID("byDecayModeFinding") < 0.5) ) {return false;}
  } else if (_RecoTau1DiscrByProngType == "1hps") {
    if (patTau.isCaloTau()) {
      if(patTau.signalTracks().size() == 1) {}
      else {return false;}
    } else {
      if( (patTau.signalPFChargedHadrCands().size() == 1) && (patTau.tauID("byDecayModeFinding") > 0.5) ) {}
      else {return false;}
    }
  } else {}

  // ----Electron and Muon vetos
  if ((_DoRecoTau1DiscrAgainstElectron) && (!(_SelectTau1sThatAreElectrons))) { if ( (patTau.tauID(_RecoTau1DiscrAgainstElectron.data()) < 0.5) ) {return false;} }
  if (_SelectTau1sThatAreElectrons) { if ( (patTau.tauID(_RecoTau1DiscrAgainstElectron.data()) > 0.5) ) {return false;} }
  if (_DoRecoTau1DiscrByCrackCut) { if(isInTheCracks(smearedTauPtEtaPhiMVector.at(nobj).eta())) {return false;} }
  if ((_DoRecoTau1DiscrAgainstMuon) && (!(_SelectTau1sThatAreMuons))) { if ( (patTau.tauID(_RecoTau1DiscrAgainstMuon.data()) < 0.5) ) {return false;} }
  if (_SelectTau1sThatAreMuons) { if ( (patTau.tauID(_RecoTau1DiscrAgainstMuon.data()) > 0.5) ) {return false;} }

  return true;
}
bool HiMassTauAnalysis::passRecoTau2Cuts(const pat::Tau& patTau,int nobj) {
  // ----Matching to gen
  if((_MatchTauToGen) && (!isData)) {
    if(_GenParticleSource.label() != "") {
      if(!(matchToGen(patTau).first)) {return false;}
    } else {return false;}
  }
  // ----Acceptance cuts
  if (fabs(smearedTauPtEtaPhiMVector.at(nobj).eta())>_RecoTau2EtaCut) {return false;}
  if (smearedTauPtEtaPhiMVector.at(nobj).pt()<_RecoTau2PtMinCut) {return false;}
  if (smearedTauPtEtaPhiMVector.at(nobj).pt()>_RecoTau2PtMaxCut) {return false;}
  // ----Lead track requirement
  if (_DoRecoTau2DiscrByLeadTrack) {
    if (_UseRecoTau2DiscrByLeadTrackFlag) { if ( (patTau.tauID(_RecoTau2DiscrByLeadTrack.data()) < 0.5) ) {return false;} }
    else {
      if (patTau.isCaloTau()) { if( (!(patTau.leadTrack().isNonnull())) || (patTau.leadTrack()->pt() < _RecoTau2LeadTrackThreshold) ) {return false;} }
      else { if( (!(patTau.leadPFChargedHadrCand().isNonnull())) || (patTau.leadPFChargedHadrCand()->pt() < _RecoTau2LeadTrackThreshold) ) {return false;} }
    }
  }
  // ----"eVeto" - Lead track minimimum hits requirement && H3x3/LTpt
  if (_DoRecoTau2DiscrByLeadTrackNhits) {
    if (patTau.isCaloTau()) {
      if( (!(patTau.leadTrack().isNonnull())) || ((int)(patTau.leadTrack()->numberOfValidHits()) < _RecoTau2LeadTrackMinHits) ) {return false;}
    } else {
      if( (!(patTau.leadPFChargedHadrCand().isNonnull())) || (!(patTau.leadPFChargedHadrCand()->trackRef().isNonnull())) ) {return false;}
      if( (int)(patTau.leadPFChargedHadrCand()->trackRef()->numberOfValidHits()) < _RecoTau2LeadTrackMinHits ) {return false;}
    }
  }
  if (_DoRecoTau2DiscrByH3x3OverP) {
    if (patTau.isCaloTau()) {
      if( (!(patTau.leadTrack().isNonnull())) || ((patTau.hcal3x3OverPLead()) <= _RecoTau2H3x3OverP) ) {return false;}
    } else {
      if( (!(patTau.leadPFChargedHadrCand().isNonnull())) || (patTau.hcal3x3OverPLead() <= _RecoTau2H3x3OverP) ) {return false;}
    }
  }
  // ----Isolation requirement
  if (_DoRecoTau2DiscrByIsolation) {
    if (_UseRecoTau2DiscrByIsolationFlag) {if ( (patTau.tauID(_RecoTau2DiscrByIsolation.data()) < 0.5) ) {return false;}}
    else {
      if(_UseRecoTau2IsoSumPtInsteadOfNiso) {
        _RecoTauIsoDeltaRCone = _RecoTau2IsoDeltaRCone;
        _RecoTauTrackIsoTrkThreshold = _RecoTau2TrackIsoTrkThreshold;
        _RecoTauGammaIsoGamThreshold = _RecoTau2GammaIsoGamThreshold;
        _UseRecoTauEllipseForEcalIso = _UseRecoTau2EllipseForEcalIso;
        _RecoTauEcalIsoRphiForEllipse = _RecoTau2EcalIsoRphiForEllipse;
        _RecoTauEcalIsoRetaForEllipse = _RecoTau2EcalIsoRetaForEllipse;
        if( (CalculateTauTrackIsolation(patTau).second) >= _RecoTau2TrackIsoSumPtMaxCutValue ) {return false;}
        if( (CalculateTauTrackIsolation(patTau).second) < _RecoTau2TrackIsoSumPtMinCutValue ) {return false;}
        if( (CalculateTauEcalIsolation(patTau).second) >= _RecoTau2EcalIsoSumPtMaxCutValue ) {return false;}
        if( (CalculateTauEcalIsolation(patTau).second) < _RecoTau2EcalIsoSumPtMinCutValue ) {return false;}
      } else {
        _RecoTauIsoDeltaRCone = _RecoTau2IsoDeltaRCone;
        _RecoTauTrackIsoTrkThreshold = _RecoTau2TrackIsoTrkThreshold;
        _RecoTauGammaIsoGamThreshold = _RecoTau2GammaIsoGamThreshold;
        _UseRecoTauEllipseForEcalIso = _UseRecoTau2EllipseForEcalIso;
        _RecoTauEcalIsoRphiForEllipse = _RecoTau2EcalIsoRphiForEllipse;
        _RecoTauEcalIsoRetaForEllipse = _RecoTau2EcalIsoRetaForEllipse;
        if( (CalculateTauTrackIsolation(patTau).first) > _RecoTau2TrackNisoMax ) {return false;}
        if( (CalculateTauEcalIsolation(patTau).first) > _RecoTau2EcalNisoMax ) {return false;}
      }
    }
  }
  // ----Require 1 or 3 prongs
  if (_RecoTau2DiscrByProngType == "1or3") {
    if (patTau.isCaloTau()) {
      if((patTau.signalTracks().size() == 1) ||(patTau.signalTracks().size() == 3)) {}
      else {return false;}
    } else {
      if((patTau.signalPFChargedHadrCands().size() == 1) ||(patTau.signalPFChargedHadrCands().size() == 3)) {}
      else {return false;}
    }
  } else if (_RecoTau2DiscrByProngType == "1") {
    if (patTau.isCaloTau()) {
      if(patTau.signalTracks().size() == 1) {}
      else {return false;}
    } else {
      if(patTau.signalPFChargedHadrCands().size() == 1) {}
      else {return false;}
    }
  } else if (_RecoTau2DiscrByProngType == "3") {
    if (patTau.isCaloTau()) {
      if(patTau.signalTracks().size() == 3) {}
      else {return false;}
    } else {
      if(patTau.signalPFChargedHadrCands().size() == 3) {}
      else {return false;}
    }
  } else if (_RecoTau2DiscrByProngType == "1or3hps") {
    if ( (patTau.tauID("byDecayModeFinding") < 0.5) ) {return false;}
  } else if (_RecoTau2DiscrByProngType == "1hps") {
    if (patTau.isCaloTau()) {
      if(patTau.signalTracks().size() == 1) {}
      else {return false;}
    } else {
      if( (patTau.signalPFChargedHadrCands().size() == 1) && (patTau.tauID("byDecayModeFinding") > 0.5) ) {}
      else {return false;}
    }
  } else {}

  // ----Electron and Muon vetos
  if ((_DoRecoTau2DiscrAgainstElectron) && (!(_SelectTau2sThatAreElectrons))) { if ( (patTau.tauID(_RecoTau2DiscrAgainstElectron.data()) < 0.5) ) {return false;} }
  if (_SelectTau2sThatAreElectrons) { if ( (patTau.tauID(_RecoTau2DiscrAgainstElectron.data()) > 0.5) ) {return false;} }
  if (_DoRecoTau2DiscrByCrackCut) { if(isInTheCracks(smearedTauPtEtaPhiMVector.at(nobj).eta())) {return false;} }
  if ((_DoRecoTau2DiscrAgainstMuon) && (!(_SelectTau2sThatAreMuons))) { if ( (patTau.tauID(_RecoTau2DiscrAgainstMuon.data()) < 0.5) ) {return false;} }
  if (_SelectTau2sThatAreMuons) { if ( (patTau.tauID(_RecoTau2DiscrAgainstMuon.data()) > 0.5) ) {return false;} }

  return true;
}

// ---------------Apply Muon Cuts
bool HiMassTauAnalysis::passRecoMuon1Cuts(const pat::Muon& patMuon,int nobj) {
  // ----Matching to gen
  if((_MatchLeptonToGen) && (!isData)) {
    if(_GenParticleSource.label() != "") {
      if(!(matchToGen(patMuon).first)) {return false;}
    } else {return false;}
  }
  // ----Require maching tracks in silicon tracker and muon chamber
  if (_DoRecoMuon1DiscrByGlobal) {if (!patMuon.isGlobalMuon()) {return false;}}
  // ----Acceptance cuts
  if (fabs(smearedMuonPtEtaPhiMVector.at(nobj).eta())>_RecoMuon1EtaCut) {return false;}
  if (smearedMuonPtEtaPhiMVector.at(nobj).pt()<_RecoMuon1PtMinCut) {return false;}
  if (smearedMuonPtEtaPhiMVector.at(nobj).pt()>_RecoMuon1PtMaxCut) {return false;}
  // ----Isolation requirement
  if (_DoRecoMuon1DiscrByIsolation) {
    if( (patMuon.trackIsoDeposit()->depositAndCountWithin(_RecoMuon1IsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuon1TrackIsoTrkThreshold).first) 
         >= _RecoMuon1TrackIsoSumPtMaxCutValue) {return false;}
    if( (patMuon.trackIsoDeposit()->depositAndCountWithin(_RecoMuon1IsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuon1TrackIsoTrkThreshold).first) 
         < _RecoMuon1TrackIsoSumPtMinCutValue) {return false;}
    if( (patMuon.ecalIsoDeposit()->depositAndCountWithin(_RecoMuon1IsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuon1EcalIsoRecHitThreshold).first) 
         >= _RecoMuon1EcalIsoSumPtMaxCutValue) {return false;}
    if( (patMuon.ecalIsoDeposit()->depositAndCountWithin(_RecoMuon1IsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuon1EcalIsoRecHitThreshold).first) 
         < _RecoMuon1EcalIsoSumPtMinCutValue) {return false;}
  }
  // ----Impact parameter requirement
  if (_DoRecoMuon1DiscrByIp) {
    const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
    if ( patMuon.track().isNonnull() ) {if ( (fabs(patMuon.track()->dxy(thePrimaryEventVertex.position())) >= _RecoMuon1IpCut) ) {return false;}}
    else {return false;}
  }
  if (_DoRecoMuon1DiscrByPionVeto) {
    if( ((_RecoMuon1CaloCompCoefficient * muon::caloCompatibility(patMuon)) + (_RecoMuon1SegmCompCoefficient * muon::segmentCompatibility(patMuon))) <= _RecoMuon1AntiPionCut ) {return false;}
  }
/*
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
*/
  return true;
}
bool HiMassTauAnalysis::passRecoMuon2Cuts(const pat::Muon& patMuon,int nobj) {
  // ----Matching to gen
  if((_MatchLeptonToGen) && (!isData)) {
    if(_GenParticleSource.label() != "") {
      if(!(matchToGen(patMuon).first)) {return false;}
    } else {return false;}
  }
  // ----Require maching tracks in silicon tracker and muon chamber
  if (_DoRecoMuon2DiscrByGlobal) {if (!patMuon.isGlobalMuon()) {return false;}}
  // ----Acceptance cuts
  if (fabs(smearedMuonPtEtaPhiMVector.at(nobj).eta())>_RecoMuon2EtaCut) {return false;}
  if (smearedMuonPtEtaPhiMVector.at(nobj).pt()<_RecoMuon2PtMinCut) {return false;}
  if (smearedMuonPtEtaPhiMVector.at(nobj).pt()>_RecoMuon2PtMaxCut) {return false;}
  // ----Isolation requirement
  if (_DoRecoMuon2DiscrByIsolation) {
    if( (patMuon.trackIsoDeposit()->depositAndCountWithin(_RecoMuon2IsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuon2TrackIsoTrkThreshold).first) 
         >= _RecoMuon2TrackIsoSumPtMaxCutValue) {return false;}
    if( (patMuon.trackIsoDeposit()->depositAndCountWithin(_RecoMuon2IsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuon2TrackIsoTrkThreshold).first) 
         < _RecoMuon2TrackIsoSumPtMinCutValue) {return false;}
    if( (patMuon.ecalIsoDeposit()->depositAndCountWithin(_RecoMuon2IsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuon2EcalIsoRecHitThreshold).first) 
         >= _RecoMuon2EcalIsoSumPtMaxCutValue) {return false;}
    if( (patMuon.ecalIsoDeposit()->depositAndCountWithin(_RecoMuon2IsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuon2EcalIsoRecHitThreshold).first) 
         < _RecoMuon2EcalIsoSumPtMinCutValue) {return false;}
  }
  // ----Impact parameter requirement
  if (_DoRecoMuon2DiscrByIp) {
    const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
    if ( patMuon.track().isNonnull() ) {if ( (fabs(patMuon.track()->dxy(thePrimaryEventVertex.position())) >= _RecoMuon2IpCut) ) {return false;}}
    else {return false;}
  }
  if (_DoRecoMuon2DiscrByPionVeto) {
    if( ((_RecoMuon2CaloCompCoefficient * muon::caloCompatibility(patMuon)) + (_RecoMuon2SegmCompCoefficient * muon::segmentCompatibility(patMuon))) <= _RecoMuon2AntiPionCut ) {return false;}
  }
/*
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
*/
  return true;
}

//--------------Apply Electron Cuts
bool HiMassTauAnalysis::passRecoElectron1Cuts(const pat::Electron& patElectron,int nobj) {
  heep::Ele theHeepElec(patElectron);
  int cutResult = patElectron.userInt("HEEPId");
  const Track* elTrack = (const reco::Track*)(patElectron.gsfTrack().get());
  const HitPattern& pInner = elTrack->trackerExpectedHitsInner();
  // ----Matching to gen
  if((_MatchLeptonToGen) && (!isData)) {
    if(_GenParticleSource.label() != "") {
      if(!(matchToGen(patElectron).first)) {return false;}
    } else {return false;}
  }
  // ----Acceptance cuts
  if (fabs(smearedElectronPtEtaPhiMVector.at(nobj).eta())>_RecoElectron1EtaCut) {return false;}
  if (smearedElectronPtEtaPhiMVector.at(nobj).pt()<_RecoElectron1PtMinCut) {return false;}
  if (smearedElectronPtEtaPhiMVector.at(nobj).pt()>_RecoElectron1PtMaxCut) {return false;}

/*
  std::cout << "Ip = " << patElectron.electronID("eidLooseOnlyIP") << std::endl;
  std::cout << "Conv. = " << patElectron.electronID("eidLooseOnlyConversionRejection") << std::endl;
  std::cout << "Iso = " << patElectron.electronID("eidLooseOnlyIso") << std::endl;
  std::cout << "ID = " << patElectron.electronID("eidLooseOnlyID") << std::endl;
  std::cout << "Loose = " << patElectron.electronID("eidLoose") << std::endl;
*/

  // ----Isolation requirement
  if ( (_DoRecoElectron1DiscrByTrackIsolation) || (_DoRecoElectron1DiscrByEcalIsolation) ) {
    if( (patElectron.electronID("eidLooseOnlyIso") == 15) ) { }
    else {return false;}
/*
    if (_DoRecoElectron1DiscrByTrackIsolation) {
      if(patElectron.dr04TkSumPt() > _RecoElectron1TrackIsoSumPtMaxCutValue) {return false;}
      if(patElectron.dr04TkSumPt() < _RecoElectron1TrackIsoSumPtMinCutValue) {return false;}
    }
    if (_DoRecoElectron1DiscrByEcalIsolation) {
      if(patElectron.dr04EcalRecHitSumEt() > _RecoElectron1EcalIsoSumPtMaxCutValue) {return false;}
      if(patElectron.dr04EcalRecHitSumEt() < _RecoElectron1EcalIsoSumPtMinCutValue) {return false;}
    }
*/
  }
  // ----Impact parameter requirement
  if (_DoRecoElectron1DiscrByIp) {
    if( (patElectron.electronID("eidLooseOnlyIP") == 15) ) { }
    else {return false;}
  }
/*
  if (_DoRecoElectron1DiscrByIp) {
    const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
//    if ( patElectron.gsfTrack().isNonnull() ) {if ( (fabs(patElectron.gsfTrack()->dxy(thePrimaryEventVertex.position())) >= _RecoElectron1IpCut) ) {return false;}}
    if ( patElectron.track().isNonnull() ) {if ( (fabs(patElectron.track()->dxy(thePrimaryEventVertex.position())) >= _RecoElectron1IpCut) ) {return false;}}
    else {return false;}
  }
*/
/*
  // ----E over p requirement
  if (_DoRecoElectron1DiscrByEoverP) { if((patElectron.eSuperClusterOverP()>_RecoElectron1EoverPMax) || (patElectron.eSuperClusterOverP()<_RecoElectron1EoverPMin)) {return false;} }
  // ----Electromagnetic energy fraction requirement
  if (_DoRecoElectron1DiscrByEcalDrivenSeed) { if(!(patElectron.ecalDrivenSeed())) {return false;} }
  if (_DoRecoElectron1DiscrByTrackerDrivenSeed) { if(!(patElectron.trackerDrivenSeed())) {return false;} }
  if (_DoRecoElectron1DiscrByHoverEm) {
    if(_UseHeepInfo) {
      if(!(heep::CutCodes::passCuts(cutResult, "hadem"))) {return false;}
    } else { if(patElectron.hadronicOverEm() > _RecoElectron1HoverEmCut) {return false;} }
  }
  if (_DoRecoElectron1DiscrBySigmaIEtaIEta) {
    if(_UseHeepInfo) {
      if(!(heep::CutCodes::passCuts(cutResult, "sigmaIEtaIEta"))) {return false;}
    } else {
      if(patElectron.isEE()) { if(patElectron.scSigmaIEtaIEta() > _RecoElectron1SigmaIEtaIEta) {return false;} }
    }
  }
  if (_DoRecoElectron1DiscrByDEtaIn) {
    if(_UseHeepInfo) {
      if(!(heep::CutCodes::passCuts(cutResult, "dEtaIn"))) {return false;}
    } else {
      if(patElectron.isEE()) { if(fabs(patElectron.deltaEtaSuperClusterTrackAtVtx()) > _RecoElectron1EEDEtaIn) {return false;} }
      if(patElectron.isEB()) { if(fabs(patElectron.deltaEtaSuperClusterTrackAtVtx()) > _RecoElectron1EBDEtaIn) {return false;} }
    }
  }
  if (_DoRecoElectron1DiscrByDPhiIn) {
    if(_UseHeepInfo) {
      if(!(heep::CutCodes::passCuts(cutResult, "dPhiIn"))) {return false;}
    } else {
      if(patElectron.isEE()) { if(fabs(patElectron.deltaPhiSuperClusterTrackAtVtx()) > _RecoElectron1EEDPhiIn) {return false;} }
      if(patElectron.isEB()) { if(fabs(patElectron.deltaPhiSuperClusterTrackAtVtx()) > _RecoElectron1EBDPhiIn) {return false;} }
    }
  }
  if (_DoRecoElectron1DiscrBySCE2by5Over5by5) {
    if(_UseHeepInfo) {
      if(!(heep::CutCodes::passCuts(cutResult, "e2x5Over5x5"))) {return false;}
    } else {
      if(patElectron.isEB()) {
        if( ((patElectron.scE2x5Max() / patElectron.scE5x5()) > _RecoElectron1EBscE2by5Over5by5) || 
            ((patElectron.scE1x5() / patElectron.scE5x5()) > _RecoElectron1EBscE1by5Over5by5) ) { }
        else {return false;}
      }
    }
  }
*/
  if (_DoRecoElectron1DiscrBySCE2by5Over5by5) {
    if( (patElectron.electronID("eidLooseOnlyID") == 15) ) { }
    else {return false;}
  }
//  if (_DoRecoElectron1DiscrByMissingHits) { if(pInner.numberOfHits() >= _RecoElectron1MissingHits) {return false;} }
  if (_DoRecoElectron1DiscrByMissingHits) {
    if( (patElectron.electronID("eidLooseOnlyConversionRejection") == 15) ) { }
    else {return false;}
  }

  return true;
}
bool HiMassTauAnalysis::passRecoElectron2Cuts(const pat::Electron& patElectron,int nobj) {
  heep::Ele theHeepElec(patElectron);
  int cutResult = patElectron.userInt("HEEPId");
  const Track* elTrack = (const reco::Track*)(patElectron.gsfTrack().get());
  const HitPattern& pInner = elTrack->trackerExpectedHitsInner();
  // ----Matching to gen
  if((_MatchLeptonToGen) && (!isData)) {
    if(_GenParticleSource.label() != "") {
      if(!(matchToGen(patElectron).first)) {return false;}
    } else {return false;}
  }
  // ----Acceptance cuts
  if (fabs(smearedElectronPtEtaPhiMVector.at(nobj).eta())>_RecoElectron2EtaCut) {return false;}
  if (smearedElectronPtEtaPhiMVector.at(nobj).pt()<_RecoElectron2PtMinCut) {return false;}
  if (smearedElectronPtEtaPhiMVector.at(nobj).pt()>_RecoElectron2PtMaxCut) {return false;}

/*
  std::cout << "Ip = " << patElectron.electronID("eidLooseOnlyIP") << std::endl;
  std::cout << "Conv. = " << patElectron.electronID("eidLooseOnlyConversionRejection") << std::endl;
  std::cout << "Iso = " << patElectron.electronID("eidLooseOnlyIso") << std::endl;
  std::cout << "ID = " << patElectron.electronID("eidLooseOnlyID") << std::endl;
  std::cout << "Loose = " << patElectron.electronID("eidLoose") << std::endl;
*/

  // ----Isolation requirement
  if ( (_DoRecoElectron2DiscrByTrackIsolation) || (_DoRecoElectron2DiscrByEcalIsolation) ) {
    if( (patElectron.electronID("eidLooseOnlyIso") == 15) ) { }
    else {return false;}
/*
    if (_DoRecoElectron2DiscrByTrackIsolation) {
      if(patElectron.dr04TkSumPt() > _RecoElectron2TrackIsoSumPtMaxCutValue) {return false;}
      if(patElectron.dr04TkSumPt() < _RecoElectron2TrackIsoSumPtMinCutValue) {return false;}
    }
    if (_DoRecoElectron2DiscrByEcalIsolation) {
      if(patElectron.dr04EcalRecHitSumEt() > _RecoElectron2EcalIsoSumPtMaxCutValue) {return false;}
      if(patElectron.dr04EcalRecHitSumEt() < _RecoElectron2EcalIsoSumPtMinCutValue) {return false;}
    }
*/
  }
  // ----Impact parameter requirement
  if (_DoRecoElectron2DiscrByIp) {
    if( (patElectron.electronID("eidLooseOnlyIP") == 15) ) { }
    else {return false;}
  }
/*
  if (_DoRecoElectron2DiscrByIp) {
    const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
//    if ( patElectron.gsfTrack().isNonnull() ) {if ( (fabs(patElectron.gsfTrack()->dxy(thePrimaryEventVertex.position())) >= _RecoElectron2IpCut) ) {return false;}}
    if ( patElectron.track().isNonnull() ) {if ( (fabs(patElectron.track()->dxy(thePrimaryEventVertex.position())) >= _RecoElectron2IpCut) ) {return false;}}
    else {return false;}
  }
*/
/*
  // ----E over p requirement
  if (_DoRecoElectron2DiscrByEoverP) { if((patElectron.eSuperClusterOverP()>_RecoElectron2EoverPMax) || (patElectron.eSuperClusterOverP()<_RecoElectron2EoverPMin)) {return false;} }
  // ----Electromagnetic energy fraction requirement
  if (_DoRecoElectron2DiscrByEcalDrivenSeed) { if(!(patElectron.ecalDrivenSeed())) {return false;} }
  if (_DoRecoElectron2DiscrByTrackerDrivenSeed) { if(!(patElectron.trackerDrivenSeed())) {return false;} }
  if (_DoRecoElectron2DiscrByHoverEm) {
    if(_UseHeepInfo) {
      if(!(heep::CutCodes::passCuts(cutResult, "hadem"))) {return false;}
    } else { if(patElectron.hadronicOverEm() > _RecoElectron2HoverEmCut) {return false;} }
  }
  if (_DoRecoElectron2DiscrBySigmaIEtaIEta) {
    if(_UseHeepInfo) {
      if(!(heep::CutCodes::passCuts(cutResult, "sigmaIEtaIEta"))) {return false;}
    } else {
      if(patElectron.isEE()) { if(patElectron.scSigmaIEtaIEta() > _RecoElectron2SigmaIEtaIEta) {return false;} }
    }
  }
  if (_DoRecoElectron2DiscrByDEtaIn) {
    if(_UseHeepInfo) {
      if(!(heep::CutCodes::passCuts(cutResult, "dEtaIn"))) {return false;}
    } else {
      if(patElectron.isEE()) { if(fabs(patElectron.deltaEtaSuperClusterTrackAtVtx()) > _RecoElectron2EEDEtaIn) {return false;} }
      if(patElectron.isEB()) { if(fabs(patElectron.deltaEtaSuperClusterTrackAtVtx()) > _RecoElectron2EBDEtaIn) {return false;} }
    }
  }
  if (_DoRecoElectron2DiscrByDPhiIn) {
    if(_UseHeepInfo) {
      if(!(heep::CutCodes::passCuts(cutResult, "dPhiIn"))) {return false;}
    } else {
      if(patElectron.isEE()) { if(fabs(patElectron.deltaPhiSuperClusterTrackAtVtx()) > _RecoElectron2EEDPhiIn) {return false;} }
      if(patElectron.isEB()) { if(fabs(patElectron.deltaPhiSuperClusterTrackAtVtx()) > _RecoElectron2EBDPhiIn) {return false;} }
    }
  }
  if (_DoRecoElectron2DiscrBySCE2by5Over5by5) {
    if(_UseHeepInfo) {
      if(!(heep::CutCodes::passCuts(cutResult, "e2x5Over5x5"))) {return false;}
    } else {
      if(patElectron.isEB()) {
        if( ((patElectron.scE2x5Max() / patElectron.scE5x5()) > _RecoElectron2EBscE2by5Over5by5) || 
            ((patElectron.scE1x5() / patElectron.scE5x5()) > _RecoElectron2EBscE1by5Over5by5) ) { }
        else {return false;}
      }
    }
  }
*/
  if (_DoRecoElectron2DiscrBySCE2by5Over5by5) {
    if( (patElectron.electronID("eidLooseOnlyID") == 15) ) { }
    else {return false;}
  }
//  if (_DoRecoElectron2DiscrByMissingHits) { if(pInner.numberOfHits() >= _RecoElectron2MissingHits) {return false;} }
  if (_DoRecoElectron2DiscrByMissingHits) {
    if( (patElectron.electronID("eidLooseOnlyConversionRejection") == 15) ) { }
    else {return false;}
  }

  return true;
}

//--------------Apply Jet Cuts
bool HiMassTauAnalysis::passRecoJetCuts(const pat::Jet& patJet,int nobj) {
  // ----Acceptance cuts
  if (fabs(smearedJetPtEtaPhiMVector.at(nobj).eta())>_RecoJetEtaMaxCut) {return false;}
  if (fabs(smearedJetPtEtaPhiMVector.at(nobj).eta())<_RecoJetEtaMinCut) {return false;}
  if (smearedJetPtEtaPhiMVector.at(nobj).pt()<_RecoJetPtCut) {return false;}
  // ----anti-overlap requirements
  if (_RemoveJetOverlapWithMuon1s) {
    int theNumberOfMuons = 0;
    for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
      theNumberOfMuons++;
      if( (passRecoMuon1Cuts((*patMuon),theNumberOfMuons-1)) ) {
        if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedMuonMomentumVector.at(theNumberOfMuons-1)) < _JetMuon1MatchingDeltaR) {return false;}
      }
    }
  }
  if (_RemoveJetOverlapWithMuon2s) {
    int theNumberOfMuons = 0;
    for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
      theNumberOfMuons++;
      if( (passRecoMuon2Cuts((*patMuon),theNumberOfMuons-1)) ) {
        if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedMuonMomentumVector.at(theNumberOfMuons-1)) < _JetMuon2MatchingDeltaR) {return false;}
      }
    }
  }
  if (_RemoveJetOverlapWithElectron1s) {
    int theNumberOfElectrons = 0;
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron ) {
      theNumberOfElectrons++;
      if( (passRecoElectron1Cuts((*patElectron),theNumberOfElectrons-1)) ) {
        if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedElectronMomentumVector.at(theNumberOfElectrons-1)) < _JetElectron1MatchingDeltaR) {return false;}
      }
    }
  }
  if (_RemoveJetOverlapWithElectron2s) {
    int theNumberOfElectrons = 0;
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron ) {
      theNumberOfElectrons++;
      if( (passRecoElectron2Cuts((*patElectron),theNumberOfElectrons-1)) ) {
        if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedElectronMomentumVector.at(theNumberOfElectrons-1)) < _JetElectron2MatchingDeltaR) {return false;}
      }
    }
  }
  if (_RemoveJetOverlapWithTau1s) {
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if( (passRecoTau1Cuts((*patTau),theNumberOfTaus-1)) ) {
        if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedTauMomentumVector.at(theNumberOfTaus-1)) < _JetTau1MatchingDeltaR) {return false;}
      }
    }
  }
  if (_RemoveJetOverlapWithTau2s) {
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if( (passRecoTau2Cuts((*patTau),theNumberOfTaus-1)) ) {
        if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedTauMomentumVector.at(theNumberOfTaus-1)) < _JetTau2MatchingDeltaR) {return false;}
      }
    }
  }
  return true;
}

//--------------Apply Jet Cuts
bool HiMassTauAnalysis::passRecoBJetCuts(const pat::Jet& patJet,int nobj) {
  // ----Acceptance cuts
  if (fabs(smearedJetPtEtaPhiMVector.at(nobj).eta())>_RecoBJetEtaMaxCut) {return false;}
  if (fabs(smearedJetPtEtaPhiMVector.at(nobj).eta())<_RecoBJetEtaMinCut) {return false;}
  if (smearedJetPtEtaPhiMVector.at(nobj).pt()<_RecoBJetPtCut) {return false;}
  // ----anti-overlap requirements
  if (_RemoveBJetOverlapWithMuon1s) {
    int theNumberOfMuons = 0;
    for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
      theNumberOfMuons++;
      if( (passRecoMuon1Cuts((*patMuon),theNumberOfMuons-1)) ) {
        if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedMuonMomentumVector.at(theNumberOfMuons-1)) < _BJetMuon1MatchingDeltaR) {return false;}
      }
    }
  }
  if (_RemoveBJetOverlapWithMuon2s) {
    int theNumberOfMuons = 0;
    for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
      theNumberOfMuons++;
      if( (passRecoMuon2Cuts((*patMuon),theNumberOfMuons-1)) ) {
        if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedMuonMomentumVector.at(theNumberOfMuons-1)) < _BJetMuon2MatchingDeltaR) {return false;}
      }
    }
  }
  if (_RemoveBJetOverlapWithElectron1s) {
    int theNumberOfElectrons = 0;
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron ) {
      theNumberOfElectrons++;
      if( (passRecoElectron1Cuts((*patElectron),theNumberOfElectrons-1)) ) {
        if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedElectronMomentumVector.at(theNumberOfElectrons-1)) < _BJetElectron1MatchingDeltaR) {return false;}
      }
    }
  }
  if (_RemoveBJetOverlapWithElectron2s) {
    int theNumberOfElectrons = 0;
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron ) {
      theNumberOfElectrons++;
      if( (passRecoElectron2Cuts((*patElectron),theNumberOfElectrons-1)) ) {
        if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedElectronMomentumVector.at(theNumberOfElectrons-1)) < _BJetElectron2MatchingDeltaR) {return false;}
      }
    }
  }
  if (_RemoveBJetOverlapWithTau1s) {
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if( (passRecoTau1Cuts((*patTau),theNumberOfTaus-1)) ) {
        if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedTauMomentumVector.at(theNumberOfTaus-1)) < _BJetTau1MatchingDeltaR) {return false;}
      }
    }
  }
  if (_RemoveBJetOverlapWithTau2s) {
    int theNumberOfTaus = 0;
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
      theNumberOfTaus++;
      if( (passRecoTau2Cuts((*patTau),theNumberOfTaus-1)) ) {
        if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedTauMomentumVector.at(theNumberOfTaus-1)) < _BJetTau2MatchingDeltaR) {return false;}
      }
    }
  }
  if (_ApplyJetBTagging) {if(patJet.bDiscriminator("trackCountingHighEffBJetTags") <= _JetBTaggingCut) {return false;}}
//  if (_ApplyJetBTagging) {if(patJet.bDiscriminator(_bTagger.data()) <= _JetBTaggingCut) {return false;}}
  return true;
}

//--------------Apply Leading Jet Cuts
bool HiMassTauAnalysis::passRecoFirstLeadingJetCuts(const pat::Jet& patJet,int nobj) {
  if(_DoDiscrByFirstLeadingJet) {
    // ----Acceptance cuts
    if (fabs(smearedJetPtEtaPhiMVector.at(nobj).eta())>_RecoFirstLeadingJetEtaMaxCut) {return false;}
    if (fabs(smearedJetPtEtaPhiMVector.at(nobj).eta())<_RecoFirstLeadingJetEtaMinCut) {return false;}
    if (smearedJetPtEtaPhiMVector.at(nobj).pt()<_RecoFirstLeadingJetPt) {return false;}
    if (_RemoveFirstLeadingJetOverlapWithMuon1s) {
      int theNumberOfMuons = 0;
      for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
        theNumberOfMuons++;
        if( (passRecoMuon1Cuts((*patMuon),theNumberOfMuons-1)) ) {
          if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedMuonMomentumVector.at(theNumberOfMuons-1)) < _FirstLeadingJetMuon1MatchingDeltaR) {return false;}
        }
      }
    }
    if (_RemoveFirstLeadingJetOverlapWithMuon2s) {
      int theNumberOfMuons = 0;
      for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
        theNumberOfMuons++;
        if( (passRecoMuon2Cuts((*patMuon),theNumberOfMuons-1)) ) {
          if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedMuonMomentumVector.at(theNumberOfMuons-1)) < _FirstLeadingJetMuon2MatchingDeltaR) {return false;}
        }
      }
    }
    if (_RemoveFirstLeadingJetOverlapWithElectron1s) {
      int theNumberOfElectrons = 0;
      for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron ) {
        theNumberOfElectrons++;
        if( (passRecoElectron1Cuts((*patElectron),theNumberOfElectrons-1)) ) {
          if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedElectronMomentumVector.at(theNumberOfElectrons-1)) < _FirstLeadingJetElectron1MatchingDeltaR) {return false;}
        }
      }
    }
    if (_RemoveFirstLeadingJetOverlapWithElectron2s) {
      int theNumberOfElectrons = 0;
      for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron ) {
        theNumberOfElectrons++;
        if( (passRecoElectron2Cuts((*patElectron),theNumberOfElectrons-1)) ) {
          if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedElectronMomentumVector.at(theNumberOfElectrons-1)) < _FirstLeadingJetElectron2MatchingDeltaR) {return false;}
        }
      }
    }
    if (_RemoveFirstLeadingJetOverlapWithTau1s) {
      int theNumberOfTaus = 0;
      for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
        theNumberOfTaus++;
        if( (passRecoTau1Cuts((*patTau),theNumberOfTaus-1)) ) {
          if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedTauMomentumVector.at(theNumberOfTaus-1)) < _FirstLeadingJetTau1MatchingDeltaR) {return false;}
        }
      }
    }
    if (_RemoveFirstLeadingJetOverlapWithTau2s) {
      int theNumberOfTaus = 0;
      for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
        theNumberOfTaus++;
        if( (passRecoTau2Cuts((*patTau),theNumberOfTaus-1)) ) {
          if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedTauMomentumVector.at(theNumberOfTaus-1)) < _FirstLeadingJetTau2MatchingDeltaR) {return false;}
        }
      }
    }
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
    if (_RemoveSecondLeadingJetOverlapWithMuon1s) {
      int theNumberOfMuons = 0;
      for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
        theNumberOfMuons++;
        if( (passRecoMuon1Cuts((*patMuon),theNumberOfMuons-1)) ) {
          if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedMuonMomentumVector.at(theNumberOfMuons-1)) < _SecondLeadingJetMuon1MatchingDeltaR) {return false;}
        }
      }
    }
    if (_RemoveSecondLeadingJetOverlapWithMuon2s) {
      int theNumberOfMuons = 0;
      for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
        theNumberOfMuons++;
        if( (passRecoMuon2Cuts((*patMuon),theNumberOfMuons-1)) ) {
          if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedMuonMomentumVector.at(theNumberOfMuons-1)) < _SecondLeadingJetMuon2MatchingDeltaR) {return false;}
        }
      }
    }
    if (_RemoveSecondLeadingJetOverlapWithElectron1s) {
      int theNumberOfElectrons = 0;
      for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron ) {
        theNumberOfElectrons++;
        if( (passRecoElectron1Cuts((*patElectron),theNumberOfElectrons-1)) ) {
          if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedElectronMomentumVector.at(theNumberOfElectrons-1)) < _SecondLeadingJetElectron1MatchingDeltaR) {return false;}
        }
      }
    }
    if (_RemoveSecondLeadingJetOverlapWithElectron2s) {
      int theNumberOfElectrons = 0;
      for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron ) {
        theNumberOfElectrons++;
        if( (passRecoElectron2Cuts((*patElectron),theNumberOfElectrons-1)) ) {
          if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedElectronMomentumVector.at(theNumberOfElectrons-1)) < _SecondLeadingJetElectron2MatchingDeltaR) {return false;}
        }
      }
    }
    if (_RemoveSecondLeadingJetOverlapWithTau1s) {
      int theNumberOfTaus = 0;
      for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
        theNumberOfTaus++;
        if( (passRecoTau1Cuts((*patTau),theNumberOfTaus-1)) ) {
          if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedTauMomentumVector.at(theNumberOfTaus-1)) < _SecondLeadingJetTau1MatchingDeltaR) {return false;}
        }
      }
    }
    if (_RemoveSecondLeadingJetOverlapWithTau2s) {
      int theNumberOfTaus = 0;
      for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
        theNumberOfTaus++;
        if( (passRecoTau2Cuts((*patTau),theNumberOfTaus-1)) ) {
          if(reco::deltaR(smearedJetMomentumVector.at(nobj), smearedTauMomentumVector.at(theNumberOfTaus-1)) < _SecondLeadingJetTau2MatchingDeltaR) {return false;}
        }
      }
    }
  }
  return true;
}

// ---------------Apply Topology Cuts
bool HiMassTauAnalysis::passMuon1Tau1TopologyCuts(const pat::Tau& patTau, int nobjT, const pat::Muon& patMuon, int nobjM) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoMuon1Tau1DiscrByDeltaR) { if(reco::deltaR(smearedTauMomentumVector.at(nobjT), smearedMuonMomentumVector.at(nobjM)) < _Muon1Tau1DeltaRCut) {return false;} }
  // ---- apply Zll veto cut
  if (_DoMuon1DiscrByIsZllCut) { if((isZmm(smearedMuonMomentumVector.at(nobjM)).first)) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_Muon1Tau1DiscrByOSLSType == "OS") {
    if(_UseTauSeedTrackForMuon1Tau1DiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patMuon.charge() * patTau.leadTrack()->charge()) >= 0) {return false;} }
      else { if((patMuon.charge() * patTau.leadPFChargedHadrCand()->charge()) >= 0) {return false;} }
    } else { if((patMuon.charge() * patTau.charge()) >= 0) {return false;} }
  } else if (_Muon1Tau1DiscrByOSLSType == "LS") {
    if(_UseTauSeedTrackForMuon1Tau1DiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patMuon.charge() * patTau.leadTrack()->charge()) <= 0) {return false;} }
      else { if((patMuon.charge() * patTau.leadPFChargedHadrCand()->charge()) <= 0) {return false;} }
    } else { if((patMuon.charge() * patTau.charge()) <= 0) {return false;} }
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoMuon1Tau1DiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) > _Muon1Tau1CosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) < _Muon1Tau1CosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByMuon1Tau1MassReco) {
    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon1Tau1ProductsAndMetMassReco;
    _UseCollinerApproxMassReco = _UseCollinerApproxMuon1Tau1MassReco;
    if( (CalculateThe4Momentum(patTau,nobjT,patMuon,nobjM).second.M() < _Muon1Tau1MassMinCut) || (CalculateThe4Momentum(patTau,nobjT,patMuon,nobjM).second.M() > _Muon1Tau1MassMaxCut) ) {return false;}
  }
  // ----Zeta requirement
  if (_DoMuon1Tau1DiscrByCDFzeta2D) {
    if( ((_Muon1Tau1PZetaCutCoefficient * CalculatePZeta(patTau,nobjT,patMuon,nobjM)) + 
         (_Muon1Tau1PZetaVisCutCoefficient * CalculatePZetaVis(patTau,nobjT,patMuon,nobjM))) < _Muon1Tau1CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoMuon1Tau1DiscrByDeltaPtDivSumPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedMuonPtEtaPhiMVector.at(nobjM).pt())) < _Muon1Tau1DeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedMuonPtEtaPhiMVector.at(nobjM).pt())) > _Muon1Tau1DeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoMuon1Tau1DiscrByDeltaPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt())) < _Muon1Tau1DeltaPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt())) > _Muon1Tau1DeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByMuon1MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) > _Muon1MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) < _Muon1MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByTau1MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) > _Tau1MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) < _Tau1MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByMuon1MetMt) {
    if( (CalculateLeptonMetMt(patMuon,nobjM)<_Muon1MetMtMinCut) || (CalculateLeptonMetMt(patMuon,nobjM)>_Muon1MetMtMaxCut) ) {return false;}
  }
  if (_DoDiscrByTau1MetMt) {
    if( (CalculateLeptonMetMt(patTau,nobjT)<_Tau1MetMtMinCut) || (CalculateLeptonMetMt(patTau,nobjT)>_Tau1MetMtMaxCut) ) {return false;}
  }
  return true;
}
bool HiMassTauAnalysis::passMuon1Tau2TopologyCuts(const pat::Tau& patTau, int nobjT, const pat::Muon& patMuon, int nobjM) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoMuon1Tau2DiscrByDeltaR) { if(reco::deltaR(smearedTauMomentumVector.at(nobjT), smearedMuonMomentumVector.at(nobjM)) < _Muon1Tau2DeltaRCut) {return false;} }
  // ---- apply Zll veto cut
  if (_DoMuon1DiscrByIsZllCut) { if((isZmm(smearedMuonMomentumVector.at(nobjM)).first)) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_Muon1Tau2DiscrByOSLSType == "OS") {
    if(_UseTauSeedTrackForMuon1Tau2DiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patMuon.charge() * patTau.leadTrack()->charge()) >= 0) {return false;} }
      else { if((patMuon.charge() * patTau.leadPFChargedHadrCand()->charge()) >= 0) {return false;} }
    } else { if((patMuon.charge() * patTau.charge()) >= 0) {return false;} }
  } else if (_Muon1Tau2DiscrByOSLSType == "LS") {
    if(_UseTauSeedTrackForMuon1Tau2DiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patMuon.charge() * patTau.leadTrack()->charge()) <= 0) {return false;} }
      else { if((patMuon.charge() * patTau.leadPFChargedHadrCand()->charge()) <= 0) {return false;} }
    } else { if((patMuon.charge() * patTau.charge()) <= 0) {return false;} }
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoMuon1Tau2DiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) > _Muon1Tau2CosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) < _Muon1Tau2CosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByMuon1Tau2MassReco) {
    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon1Tau2ProductsAndMetMassReco;
    _UseCollinerApproxMassReco = _UseCollinerApproxMuon1Tau2MassReco;
    if( (CalculateThe4Momentum(patTau,nobjT,patMuon,nobjM).second.M() < _Muon1Tau2MassMinCut) || (CalculateThe4Momentum(patTau,nobjT,patMuon,nobjM).second.M() > _Muon1Tau2MassMaxCut) ) {return false;}
  }
  // ----Zeta requirement
  if (_DoMuon1Tau2DiscrByCDFzeta2D) {
    if( ((_Muon1Tau2PZetaCutCoefficient * CalculatePZeta(patTau,nobjT,patMuon,nobjM)) + 
         (_Muon1Tau2PZetaVisCutCoefficient * CalculatePZetaVis(patTau,nobjT,patMuon,nobjM))) < _Muon1Tau2CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoMuon1Tau2DiscrByDeltaPtDivSumPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedMuonPtEtaPhiMVector.at(nobjM).pt())) < _Muon1Tau2DeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedMuonPtEtaPhiMVector.at(nobjM).pt())) > _Muon1Tau2DeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoMuon1Tau2DiscrByDeltaPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt())) < _Muon1Tau2DeltaPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt())) > _Muon1Tau2DeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByMuon1MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) > _Muon1MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) < _Muon1MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByTau2MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) > _Tau2MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) < _Tau2MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByMuon1MetMt) {
    if( (CalculateLeptonMetMt(patMuon,nobjM)<_Muon1MetMtMinCut) || (CalculateLeptonMetMt(patMuon,nobjM)>_Muon1MetMtMaxCut) ) {return false;}
  }
  if (_DoDiscrByTau2MetMt) {
    if( (CalculateLeptonMetMt(patTau,nobjT)<_Tau2MetMtMinCut) || (CalculateLeptonMetMt(patTau,nobjT)>_Tau2MetMtMaxCut) ) {return false;}
  }
  return true;
}
bool HiMassTauAnalysis::passMuon2Tau1TopologyCuts(const pat::Tau& patTau, int nobjT, const pat::Muon& patMuon, int nobjM) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoMuon2Tau1DiscrByDeltaR) { if(reco::deltaR(smearedTauMomentumVector.at(nobjT), smearedMuonMomentumVector.at(nobjM)) < _Muon2Tau1DeltaRCut) {return false;} }
  // ---- apply Zll veto cut
  if (_DoMuon2DiscrByIsZllCut) { if((isZmm(smearedMuonMomentumVector.at(nobjM)).first)) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_Muon2Tau1DiscrByOSLSType == "OS") {
    if(_UseTauSeedTrackForMuon2Tau1DiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patMuon.charge() * patTau.leadTrack()->charge()) >= 0) {return false;} }
      else { if((patMuon.charge() * patTau.leadPFChargedHadrCand()->charge()) >= 0) {return false;} }
    } else { if((patMuon.charge() * patTau.charge()) >= 0) {return false;} }
  } else if (_Muon2Tau1DiscrByOSLSType == "LS") {
    if(_UseTauSeedTrackForMuon2Tau1DiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patMuon.charge() * patTau.leadTrack()->charge()) <= 0) {return false;} }
      else { if((patMuon.charge() * patTau.leadPFChargedHadrCand()->charge()) <= 0) {return false;} }
    } else { if((patMuon.charge() * patTau.charge()) <= 0) {return false;} }
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoMuon2Tau1DiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) > _Muon2Tau1CosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) < _Muon2Tau1CosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByMuon2Tau1MassReco) {
    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon2Tau1ProductsAndMetMassReco;
    _UseCollinerApproxMassReco = _UseCollinerApproxMuon2Tau1MassReco;
    if( (CalculateThe4Momentum(patTau,nobjT,patMuon,nobjM).second.M() < _Muon2Tau1MassMinCut) || (CalculateThe4Momentum(patTau,nobjT,patMuon,nobjM).second.M() > _Muon2Tau1MassMaxCut) ) {return false;}
  }
  // ----Zeta requirement
  if (_DoMuon2Tau1DiscrByCDFzeta2D) {
    if( ((_Muon2Tau1PZetaCutCoefficient * CalculatePZeta(patTau,nobjT,patMuon,nobjM)) + 
         (_Muon2Tau1PZetaVisCutCoefficient * CalculatePZetaVis(patTau,nobjT,patMuon,nobjM))) < _Muon2Tau1CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoMuon2Tau1DiscrByDeltaPtDivSumPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedMuonPtEtaPhiMVector.at(nobjM).pt())) < _Muon2Tau1DeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedMuonPtEtaPhiMVector.at(nobjM).pt())) > _Muon2Tau1DeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoMuon2Tau1DiscrByDeltaPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt())) < _Muon2Tau1DeltaPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt())) > _Muon2Tau1DeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByMuon2MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) > _Muon2MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) < _Muon2MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByTau1MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) > _Tau1MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) < _Tau1MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByMuon2MetMt) {
    if( (CalculateLeptonMetMt(patMuon,nobjM)<_Muon2MetMtMinCut) || (CalculateLeptonMetMt(patMuon,nobjM)>_Muon2MetMtMaxCut) ) {return false;}
  }
  if (_DoDiscrByTau1MetMt) {
    if( (CalculateLeptonMetMt(patTau,nobjT)<_Tau1MetMtMinCut) || (CalculateLeptonMetMt(patTau,nobjT)>_Tau1MetMtMaxCut) ) {return false;}
  }
  return true;
}
bool HiMassTauAnalysis::passMuon2Tau2TopologyCuts(const pat::Tau& patTau, int nobjT, const pat::Muon& patMuon, int nobjM) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoMuon2Tau2DiscrByDeltaR) { if(reco::deltaR(smearedTauMomentumVector.at(nobjT), smearedMuonMomentumVector.at(nobjM)) < _Muon2Tau2DeltaRCut) {return false;} }
  // ---- apply Zll veto cut
  if (_DoMuon2DiscrByIsZllCut) { if((isZmm(smearedMuonMomentumVector.at(nobjM)).first)) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_Muon2Tau2DiscrByOSLSType == "OS") {
    if(_UseTauSeedTrackForMuon2Tau2DiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patMuon.charge() * patTau.leadTrack()->charge()) >= 0) {return false;} }
      else { if((patMuon.charge() * patTau.leadPFChargedHadrCand()->charge()) >= 0) {return false;} }
    } else { if((patMuon.charge() * patTau.charge()) >= 0) {return false;} }
  } else if (_Muon2Tau2DiscrByOSLSType == "LS") {
    if(_UseTauSeedTrackForMuon2Tau2DiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patMuon.charge() * patTau.leadTrack()->charge()) <= 0) {return false;} }
      else { if((patMuon.charge() * patTau.leadPFChargedHadrCand()->charge()) <= 0) {return false;} }
    } else { if((patMuon.charge() * patTau.charge()) <= 0) {return false;} }
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoMuon2Tau2DiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) > _Muon2Tau2CosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) < _Muon2Tau2CosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByMuon2Tau2MassReco) {
    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon2Tau2ProductsAndMetMassReco;
    _UseCollinerApproxMassReco = _UseCollinerApproxMuon2Tau2MassReco;
    if( (CalculateThe4Momentum(patTau,nobjT,patMuon,nobjM).second.M() < _Muon2Tau2MassMinCut) || (CalculateThe4Momentum(patTau,nobjT,patMuon,nobjM).second.M() > _Muon2Tau2MassMaxCut) ) {return false;}
  }
  // ----Zeta requirement
  if (_DoMuon2Tau2DiscrByCDFzeta2D) {
    if( ((_Muon2Tau2PZetaCutCoefficient * CalculatePZeta(patTau,nobjT,patMuon,nobjM)) + 
         (_Muon2Tau2PZetaVisCutCoefficient * CalculatePZetaVis(patTau,nobjT,patMuon,nobjM))) < _Muon2Tau2CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoMuon2Tau2DiscrByDeltaPtDivSumPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedMuonPtEtaPhiMVector.at(nobjM).pt())) < _Muon2Tau2DeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedMuonPtEtaPhiMVector.at(nobjM).pt())) > _Muon2Tau2DeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoMuon2Tau2DiscrByDeltaPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt())) < _Muon2Tau2DeltaPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedMuonPtEtaPhiMVector.at(nobjM).pt())) > _Muon2Tau2DeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByMuon2MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) > _Muon2MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) < _Muon2MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByTau2MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) > _Tau2MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) < _Tau2MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByMuon2MetMt) {
    if( (CalculateLeptonMetMt(patMuon,nobjM)<_Muon2MetMtMinCut) || (CalculateLeptonMetMt(patMuon,nobjM)>_Muon2MetMtMaxCut) ) {return false;}
  }
  if (_DoDiscrByTau2MetMt) {
    if( (CalculateLeptonMetMt(patTau,nobjT)<_Tau2MetMtMinCut) || (CalculateLeptonMetMt(patTau,nobjT)>_Tau2MetMtMaxCut) ) {return false;}
  }
  return true;
}
bool HiMassTauAnalysis::passElectron1Tau1TopologyCuts(const pat::Tau& patTau, int nobjT, const pat::Electron& patElectron, int nobjE) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoElectron1Tau1DiscrByDeltaR) { if(reco::deltaR(smearedTauMomentumVector.at(nobjT), smearedElectronMomentumVector.at(nobjE)) < _Electron1Tau1DeltaRCut) {return false;} }
  // ---- apply Zll veto cut
  if (_DoElectron1DiscrByIsZllCut) { if((isZee(smearedElectronMomentumVector.at(nobjE)).first)) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_Electron1Tau1DiscrByOSLSType == "OS") {
    if(_UseTauSeedTrackForElectron1Tau1DiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patElectron.charge() * patTau.leadTrack()->charge()) >= 0) {return false;} }
      else { if((patElectron.charge() * patTau.leadPFChargedHadrCand()->charge()) >= 0) {return false;} }
    } else { if((patElectron.charge() * patTau.charge()) >= 0) {return false;} }
  } else if (_Electron1Tau1DiscrByOSLSType == "LS") {
    if(_UseTauSeedTrackForElectron1Tau1DiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patElectron.charge() * patTau.leadTrack()->charge()) <= 0) {return false;} }
      else { if((patElectron.charge() * patTau.leadPFChargedHadrCand()->charge()) <= 0) {return false;} }
    } else { if((patElectron.charge() * patTau.charge()) <= 0) {return false;} }
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoElectron1Tau1DiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedElectronPtEtaPhiMVector.at(nobjE).phi()))) > _Electron1Tau1CosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedElectronPtEtaPhiMVector.at(nobjE).phi()))) < _Electron1Tau1CosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByElectron1Tau1MassReco) {
    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron1Tau1ProductsAndMetMassReco;
    _UseCollinerApproxMassReco = _UseCollinerApproxElectron1Tau1MassReco;
    if( (CalculateThe4Momentum(patTau,nobjT,patElectron,nobjE).second.M() < _Electron1Tau1MassMinCut) || (CalculateThe4Momentum(patTau,nobjT,patElectron,nobjE).second.M() > _Electron1Tau1MassMaxCut) ) {return false;}
  }
  // ----Zeta requirement
  if (_DoElectron1Tau1DiscrByCDFzeta2D) {
    if( ((_Electron1Tau1PZetaCutCoefficient * CalculatePZeta(patTau,nobjT,patElectron,nobjE)) + 
         (_Electron1Tau1PZetaVisCutCoefficient * CalculatePZetaVis(patTau,nobjT,patElectron,nobjE))) < _Electron1Tau1CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoElectron1Tau1DiscrByDeltaPtDivSumPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _Electron1Tau1DeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _Electron1Tau1DeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoElectron1Tau1DiscrByDeltaPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _Electron1Tau1DeltaPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _Electron1Tau1DeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByElectron1MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) > _Electron1MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) < _Electron1MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByTau1MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) > _Tau1MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) < _Tau1MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByElectron1MetMt) {
    if( (CalculateLeptonMetMt(patElectron,nobjE)<_Electron1MetMtMinCut) || (CalculateLeptonMetMt(patElectron,nobjE)>_Electron1MetMtMaxCut) ) {return false;}
  }
  if (_DoDiscrByTau1MetMt) {
    if( (CalculateLeptonMetMt(patTau,nobjT)<_Tau1MetMtMinCut) || (CalculateLeptonMetMt(patTau,nobjT)>_Tau1MetMtMaxCut) ) {return false;}
  }
  return true;
}
bool HiMassTauAnalysis::passElectron1Tau2TopologyCuts(const pat::Tau& patTau, int nobjT, const pat::Electron& patElectron, int nobjE) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoElectron1Tau2DiscrByDeltaR) { if(reco::deltaR(smearedTauMomentumVector.at(nobjT), smearedElectronMomentumVector.at(nobjE)) < _Electron1Tau2DeltaRCut) {return false;} }
  // ---- apply Zll veto cut
  if (_DoElectron1DiscrByIsZllCut) { if((isZee(smearedElectronMomentumVector.at(nobjE)).first)) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_Electron1Tau2DiscrByOSLSType == "OS") {
    if(_UseTauSeedTrackForElectron1Tau2DiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patElectron.charge() * patTau.leadTrack()->charge()) >= 0) {return false;} }
      else { if((patElectron.charge() * patTau.leadPFChargedHadrCand()->charge()) >= 0) {return false;} }
    } else { if((patElectron.charge() * patTau.charge()) >= 0) {return false;} }
  } else if (_Electron1Tau2DiscrByOSLSType == "LS") {
    if(_UseTauSeedTrackForElectron1Tau2DiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patElectron.charge() * patTau.leadTrack()->charge()) <= 0) {return false;} }
      else { if((patElectron.charge() * patTau.leadPFChargedHadrCand()->charge()) <= 0) {return false;} }
    } else { if((patElectron.charge() * patTau.charge()) <= 0) {return false;} }
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoElectron1Tau2DiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedElectronPtEtaPhiMVector.at(nobjE).phi()))) > _Electron1Tau2CosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedElectronPtEtaPhiMVector.at(nobjE).phi()))) < _Electron1Tau2CosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByElectron1Tau2MassReco) {
    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron1Tau2ProductsAndMetMassReco;
    _UseCollinerApproxMassReco = _UseCollinerApproxElectron1Tau2MassReco;
    if( (CalculateThe4Momentum(patTau,nobjT,patElectron,nobjE).second.M() < _Electron1Tau2MassMinCut) || (CalculateThe4Momentum(patTau,nobjT,patElectron,nobjE).second.M() > _Electron1Tau2MassMaxCut) ) {return false;}
  }
  // ----Zeta requirement
  if (_DoElectron1Tau2DiscrByCDFzeta2D) {
    if( ((_Electron1Tau2PZetaCutCoefficient * CalculatePZeta(patTau,nobjT,patElectron,nobjE)) + 
         (_Electron1Tau2PZetaVisCutCoefficient * CalculatePZetaVis(patTau,nobjT,patElectron,nobjE))) < _Electron1Tau2CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoElectron1Tau2DiscrByDeltaPtDivSumPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _Electron1Tau2DeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _Electron1Tau2DeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoElectron1Tau2DiscrByDeltaPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _Electron1Tau2DeltaPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _Electron1Tau2DeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByElectron1MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) > _Electron1MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) < _Electron1MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByTau2MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) > _Tau2MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) < _Tau2MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByElectron1MetMt) {
    if( (CalculateLeptonMetMt(patElectron,nobjE)<_Electron1MetMtMinCut) || (CalculateLeptonMetMt(patElectron,nobjE)>_Electron1MetMtMaxCut) ) {return false;}
  }
  if (_DoDiscrByTau2MetMt) {
    if( (CalculateLeptonMetMt(patTau,nobjT)<_Tau2MetMtMinCut) || (CalculateLeptonMetMt(patTau,nobjT)>_Tau2MetMtMaxCut) ) {return false;}
  }
  return true;
}
bool HiMassTauAnalysis::passElectron2Tau1TopologyCuts(const pat::Tau& patTau, int nobjT, const pat::Electron& patElectron, int nobjE) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoElectron2Tau1DiscrByDeltaR) { if(reco::deltaR(smearedTauMomentumVector.at(nobjT), smearedElectronMomentumVector.at(nobjE)) < _Electron2Tau1DeltaRCut) {return false;} }
  // ---- apply Zll veto cut
  if (_DoElectron2DiscrByIsZllCut) { if((isZee(smearedElectronMomentumVector.at(nobjE)).first)) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_Electron2Tau1DiscrByOSLSType == "OS") {
    if(_UseTauSeedTrackForElectron2Tau1DiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patElectron.charge() * patTau.leadTrack()->charge()) >= 0) {return false;} }
      else { if((patElectron.charge() * patTau.leadPFChargedHadrCand()->charge()) >= 0) {return false;} }
    } else { if((patElectron.charge() * patTau.charge()) >= 0) {return false;} }
  } else if (_Electron2Tau1DiscrByOSLSType == "LS") {
    if(_UseTauSeedTrackForElectron2Tau1DiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patElectron.charge() * patTau.leadTrack()->charge()) <= 0) {return false;} }
      else { if((patElectron.charge() * patTau.leadPFChargedHadrCand()->charge()) <= 0) {return false;} }
    } else { if((patElectron.charge() * patTau.charge()) <= 0) {return false;} }
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoElectron2Tau1DiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedElectronPtEtaPhiMVector.at(nobjE).phi()))) > _Electron2Tau1CosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedElectronPtEtaPhiMVector.at(nobjE).phi()))) < _Electron2Tau1CosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByElectron2Tau1MassReco) {
    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron2Tau1ProductsAndMetMassReco;
    _UseCollinerApproxMassReco = _UseCollinerApproxElectron2Tau1MassReco;
    if( (CalculateThe4Momentum(patTau,nobjT,patElectron,nobjE).second.M() < _Electron2Tau1MassMinCut) || (CalculateThe4Momentum(patTau,nobjT,patElectron,nobjE).second.M() > _Electron2Tau1MassMaxCut) ) {return false;}
  }
  // ----Zeta requirement
  if (_DoElectron2Tau1DiscrByCDFzeta2D) {
    if( ((_Electron2Tau1PZetaCutCoefficient * CalculatePZeta(patTau,nobjT,patElectron,nobjE)) + 
         (_Electron2Tau1PZetaVisCutCoefficient * CalculatePZetaVis(patTau,nobjT,patElectron,nobjE))) < _Electron2Tau1CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoElectron2Tau1DiscrByDeltaPtDivSumPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _Electron2Tau1DeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _Electron2Tau1DeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoElectron2Tau1DiscrByDeltaPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _Electron2Tau1DeltaPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _Electron2Tau1DeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByElectron2MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) > _Electron2MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) < _Electron2MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByTau1MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) > _Tau1MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) < _Tau1MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByElectron2MetMt) {
    if( (CalculateLeptonMetMt(patElectron,nobjE)<_Electron2MetMtMinCut) || (CalculateLeptonMetMt(patElectron,nobjE)>_Electron2MetMtMaxCut) ) {return false;}
  }
  if (_DoDiscrByTau1MetMt) {
    if( (CalculateLeptonMetMt(patTau,nobjT)<_Tau1MetMtMinCut) || (CalculateLeptonMetMt(patTau,nobjT)>_Tau1MetMtMaxCut) ) {return false;}
  }
  return true;
}
bool HiMassTauAnalysis::passElectron2Tau2TopologyCuts(const pat::Tau& patTau, int nobjT, const pat::Electron& patElectron, int nobjE) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoElectron2Tau2DiscrByDeltaR) { if(reco::deltaR(smearedTauMomentumVector.at(nobjT), smearedElectronMomentumVector.at(nobjE)) < _Electron2Tau2DeltaRCut) {return false;} }
  // ---- apply Zll veto cut
  if (_DoElectron2DiscrByIsZllCut) { if((isZee(smearedElectronMomentumVector.at(nobjE)).first)) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_Electron2Tau2DiscrByOSLSType == "OS") {
    if(_UseTauSeedTrackForElectron2Tau2DiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patElectron.charge() * patTau.leadTrack()->charge()) >= 0) {return false;} }
      else { if((patElectron.charge() * patTau.leadPFChargedHadrCand()->charge()) >= 0) {return false;} }
    } else { if((patElectron.charge() * patTau.charge()) >= 0) {return false;} }
  } else if (_Electron2Tau2DiscrByOSLSType == "LS") {
    if(_UseTauSeedTrackForElectron2Tau2DiscrByOSLS) {
      if (patTau.isCaloTau()) { if((patElectron.charge() * patTau.leadTrack()->charge()) <= 0) {return false;} }
      else { if((patElectron.charge() * patTau.leadPFChargedHadrCand()->charge()) <= 0) {return false;} }
    } else { if((patElectron.charge() * patTau.charge()) <= 0) {return false;} }
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoElectron2Tau2DiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedElectronPtEtaPhiMVector.at(nobjE).phi()))) > _Electron2Tau2CosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - smearedElectronPtEtaPhiMVector.at(nobjE).phi()))) < _Electron2Tau2CosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByElectron2Tau2MassReco) {
    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron2Tau2ProductsAndMetMassReco;
    _UseCollinerApproxMassReco = _UseCollinerApproxElectron2Tau2MassReco;
    if( (CalculateThe4Momentum(patTau,nobjT,patElectron,nobjE).second.M() < _Electron2Tau2MassMinCut) || (CalculateThe4Momentum(patTau,nobjT,patElectron,nobjE).second.M() > _Electron2Tau2MassMaxCut) ) {return false;}
  }
  // ----Zeta requirement
  if (_DoElectron2Tau2DiscrByCDFzeta2D) {
    if( ((_Electron2Tau2PZetaCutCoefficient * CalculatePZeta(patTau,nobjT,patElectron,nobjE)) + 
         (_Electron2Tau2PZetaVisCutCoefficient * CalculatePZetaVis(patTau,nobjT,patElectron,nobjE))) < _Electron2Tau2CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoElectron2Tau2DiscrByDeltaPtDivSumPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _Electron2Tau2DeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedTauPtEtaPhiMVector.at(nobjT).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _Electron2Tau2DeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoElectron2Tau2DiscrByDeltaPt) {
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _Electron2Tau2DeltaPtMinCutValue ) {return false;}
    if ( ((smearedTauPtEtaPhiMVector.at(nobjT).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _Electron2Tau2DeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByElectron2MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) > _Electron2MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) < _Electron2MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByTau2MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) > _Tau2MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobjT).phi() - theMETVector.phi())) < _Tau2MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByElectron2MetMt) {
    if( (CalculateLeptonMetMt(patElectron,nobjE)<_Electron2MetMtMinCut) || (CalculateLeptonMetMt(patElectron,nobjE)>_Electron2MetMtMaxCut) ) {return false;}
  }
  if (_DoDiscrByTau2MetMt) {
    if( (CalculateLeptonMetMt(patTau,nobjT)<_Tau2MetMtMinCut) || (CalculateLeptonMetMt(patTau,nobjT)>_Tau2MetMtMaxCut) ) {return false;}
  }
  return true;
}
bool HiMassTauAnalysis::passElectron1Muon1TopologyCuts(const pat::Electron& patElectron, int nobjE, const pat::Muon& patMuon, int nobjM) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoElectron1Muon1DiscrByDeltaR) { if(reco::deltaR(smearedMuonMomentumVector.at(nobjM), smearedElectronMomentumVector.at(nobjE)) < _Electron1Muon1DeltaRCut) {return false;} }
  // ---- apply Zll veto cut
  if (_DoElectron1DiscrByIsZllCut) { if((isZee(smearedElectronMomentumVector.at(nobjE)).first)) {return false;} }
  if (_DoMuon1DiscrByIsZllCut) { if((isZmm(smearedMuonMomentumVector.at(nobjM)).first)) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_Electron1Muon1DiscrByOSLSType == "OS") {
    if((patElectron.charge() * patMuon.charge()) >= 0) {return false;}
  } else if (_Electron1Muon1DiscrByOSLSType == "LS") {
    if((patElectron.charge() * patMuon.charge()) <= 0) {return false;}
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoElectron1Muon1DiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) > _Electron1Muon1CosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) < _Electron1Muon1CosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByElectron1Muon1MassReco) {
    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron1Muon1ProductsAndMetMassReco;
    _UseCollinerApproxMassReco = _UseCollinerApproxElectron1Muon1MassReco;
    if( (CalculateThe4Momentum(patElectron,nobjE,patMuon,nobjM).second.M() < _Electron1Muon1MassMinCut) || (CalculateThe4Momentum(patElectron,nobjE,patMuon,nobjM).second.M() > _Electron1Muon1MassMaxCut) ) {return false;}
  }
  // ----Zeta requirement
  if (_DoElectron1Muon1DiscrByCDFzeta2D) {
    if( ((_Electron1Muon1PZetaCutCoefficient * CalculatePZeta(patElectron,nobjE,patMuon,nobjM)) +
         (_Electron1Muon1PZetaVisCutCoefficient * CalculatePZetaVis(patElectron,nobjE,patMuon,nobjM))) < _Electron1Muon1CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoElectron1Muon1DiscrByDeltaPtDivSumPt) {
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedMuonPtEtaPhiMVector.at(nobjM).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _Electron1Muon1DeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedMuonPtEtaPhiMVector.at(nobjM).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _Electron1Muon1DeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoElectron1Muon1DiscrByDeltaPt) {
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _Electron1Muon1DeltaPtMinCutValue ) {return false;}
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _Electron1Muon1DeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByElectron1MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) > _Electron1MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) < _Electron1MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByMuon1MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) > _Muon1MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) < _Muon1MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByElectron1MetMt) {
    if( (CalculateLeptonMetMt(patElectron,nobjE)<_Electron1MetMtMinCut) || (CalculateLeptonMetMt(patElectron,nobjE)>_Electron1MetMtMaxCut) ) {return false;}
  }
  if (_DoDiscrByMuon1MetMt) {
    if( (CalculateLeptonMetMt(patMuon,nobjM)<_Muon1MetMtMinCut) || (CalculateLeptonMetMt(patMuon,nobjM)>_Muon1MetMtMaxCut) ) {return false;}
  }
  return true;
}
bool HiMassTauAnalysis::passElectron1Muon2TopologyCuts(const pat::Electron& patElectron, int nobjE, const pat::Muon& patMuon, int nobjM) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoElectron1Muon2DiscrByDeltaR) { if(reco::deltaR(smearedMuonMomentumVector.at(nobjM), smearedElectronMomentumVector.at(nobjE)) < _Electron1Muon2DeltaRCut) {return false;} }
  // ---- apply Zll veto cut
  if (_DoElectron1DiscrByIsZllCut) { if((isZee(smearedElectronMomentumVector.at(nobjE)).first)) {return false;} }
  if (_DoMuon2DiscrByIsZllCut) { if((isZmm(smearedMuonMomentumVector.at(nobjM)).first)) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_Electron1Muon2DiscrByOSLSType == "OS") {
    if((patElectron.charge() * patMuon.charge()) >= 0) {return false;}
  } else if (_Electron1Muon2DiscrByOSLSType == "LS") {
    if((patElectron.charge() * patMuon.charge()) <= 0) {return false;}
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoElectron1Muon2DiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) > _Electron1Muon2CosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) < _Electron1Muon2CosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByElectron1Muon2MassReco) {
    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron1Muon2ProductsAndMetMassReco;
    _UseCollinerApproxMassReco = _UseCollinerApproxElectron1Muon2MassReco;
    if( (CalculateThe4Momentum(patElectron,nobjE,patMuon,nobjM).second.M() < _Electron1Muon2MassMinCut) || (CalculateThe4Momentum(patElectron,nobjE,patMuon,nobjM).second.M() > _Electron1Muon2MassMaxCut) ) {return false;}
  }
  // ----Zeta requirement
  if (_DoElectron1Muon2DiscrByCDFzeta2D) {
    if( ((_Electron1Muon2PZetaCutCoefficient * CalculatePZeta(patElectron,nobjE,patMuon,nobjM)) +
         (_Electron1Muon2PZetaVisCutCoefficient * CalculatePZetaVis(patElectron,nobjE,patMuon,nobjM))) < _Electron1Muon2CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoElectron1Muon2DiscrByDeltaPtDivSumPt) {
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedMuonPtEtaPhiMVector.at(nobjM).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _Electron1Muon2DeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedMuonPtEtaPhiMVector.at(nobjM).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _Electron1Muon2DeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoElectron1Muon2DiscrByDeltaPt) {
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _Electron1Muon2DeltaPtMinCutValue ) {return false;}
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _Electron1Muon2DeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByElectron1MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) > _Electron1MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) < _Electron1MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByMuon2MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) > _Muon2MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) < _Muon2MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByElectron1MetMt) {
    if( (CalculateLeptonMetMt(patElectron,nobjE)<_Electron1MetMtMinCut) || (CalculateLeptonMetMt(patElectron,nobjE)>_Electron1MetMtMaxCut) ) {return false;}
  }
  if (_DoDiscrByMuon2MetMt) {
    if( (CalculateLeptonMetMt(patMuon,nobjM)<_Muon2MetMtMinCut) || (CalculateLeptonMetMt(patMuon,nobjM)>_Muon2MetMtMaxCut) ) {return false;}
  }
  return true;
}
bool HiMassTauAnalysis::passElectron2Muon1TopologyCuts(const pat::Electron& patElectron, int nobjE, const pat::Muon& patMuon, int nobjM) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoElectron2Muon1DiscrByDeltaR) { if(reco::deltaR(smearedMuonMomentumVector.at(nobjM), smearedElectronMomentumVector.at(nobjE)) < _Electron2Muon1DeltaRCut) {return false;} }
  // ---- apply Zll veto cut
  if (_DoElectron2DiscrByIsZllCut) { if((isZee(smearedElectronMomentumVector.at(nobjE)).first)) {return false;} }
  if (_DoMuon1DiscrByIsZllCut) { if((isZmm(smearedMuonMomentumVector.at(nobjM)).first)) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_Electron2Muon1DiscrByOSLSType == "OS") {
    if((patElectron.charge() * patMuon.charge()) >= 0) {return false;}
  } else if (_Electron2Muon1DiscrByOSLSType == "LS") {
    if((patElectron.charge() * patMuon.charge()) <= 0) {return false;}
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoElectron2Muon1DiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) > _Electron2Muon1CosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) < _Electron2Muon1CosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByElectron2Muon1MassReco) {
    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron2Muon1ProductsAndMetMassReco;
    _UseCollinerApproxMassReco = _UseCollinerApproxElectron2Muon1MassReco;
    if( (CalculateThe4Momentum(patElectron,nobjE,patMuon,nobjM).second.M() < _Electron2Muon1MassMinCut) || (CalculateThe4Momentum(patElectron,nobjE,patMuon,nobjM).second.M() > _Electron2Muon1MassMaxCut) ) {return false;}
  }
  // ----Zeta requirement
  if (_DoElectron2Muon1DiscrByCDFzeta2D) {
    if( ((_Electron2Muon1PZetaCutCoefficient * CalculatePZeta(patElectron,nobjE,patMuon,nobjM)) +
         (_Electron2Muon1PZetaVisCutCoefficient * CalculatePZetaVis(patElectron,nobjE,patMuon,nobjM))) < _Electron2Muon1CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoElectron2Muon1DiscrByDeltaPtDivSumPt) {
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedMuonPtEtaPhiMVector.at(nobjM).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _Electron2Muon1DeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedMuonPtEtaPhiMVector.at(nobjM).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _Electron2Muon1DeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoElectron2Muon1DiscrByDeltaPt) {
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _Electron2Muon1DeltaPtMinCutValue ) {return false;}
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _Electron2Muon1DeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByElectron2MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) > _Electron2MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) < _Electron2MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByMuon1MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) > _Muon1MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) < _Muon1MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByElectron2MetMt) {
    if( (CalculateLeptonMetMt(patElectron,nobjE)<_Electron2MetMtMinCut) || (CalculateLeptonMetMt(patElectron,nobjE)>_Electron2MetMtMaxCut) ) {return false;}
  }
  if (_DoDiscrByMuon1MetMt) {
    if( (CalculateLeptonMetMt(patMuon,nobjM)<_Muon1MetMtMinCut) || (CalculateLeptonMetMt(patMuon,nobjM)>_Muon1MetMtMaxCut) ) {return false;}
  }
  return true;
}
bool HiMassTauAnalysis::passElectron2Muon2TopologyCuts(const pat::Electron& patElectron, int nobjE, const pat::Muon& patMuon, int nobjM) {
  // ----Separation cut between lepton and tau jet (remove overlaps)
  if (_DoElectron2Muon2DiscrByDeltaR) { if(reco::deltaR(smearedMuonMomentumVector.at(nobjM), smearedElectronMomentumVector.at(nobjE)) < _Electron2Muon2DeltaRCut) {return false;} }
  // ---- apply Zll veto cut
  if (_DoElectron2DiscrByIsZllCut) { if((isZee(smearedElectronMomentumVector.at(nobjE)).first)) {return false;} }
  if (_DoMuon2DiscrByIsZllCut) { if((isZmm(smearedMuonMomentumVector.at(nobjM)).first)) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_Electron2Muon2DiscrByOSLSType == "OS") {
    if((patElectron.charge() * patMuon.charge()) >= 0) {return false;}
  } else if (_Electron2Muon2DiscrByOSLSType == "LS") {
    if((patElectron.charge() * patMuon.charge()) <= 0) {return false;}
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoElectron2Muon2DiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) > _Electron2Muon2CosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - smearedMuonPtEtaPhiMVector.at(nobjM).phi()))) < _Electron2Muon2CosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByElectron2Muon2MassReco) {
    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron2Muon2ProductsAndMetMassReco;
    _UseCollinerApproxMassReco = _UseCollinerApproxElectron2Muon2MassReco;
    if( (CalculateThe4Momentum(patElectron,nobjE,patMuon,nobjM).second.M() < _Electron2Muon2MassMinCut) || (CalculateThe4Momentum(patElectron,nobjE,patMuon,nobjM).second.M() > _Electron2Muon2MassMaxCut) ) {return false;}
  }
  // ----Zeta requirement
  if (_DoElectron2Muon2DiscrByCDFzeta2D) {
    if( ((_Electron2Muon2PZetaCutCoefficient * CalculatePZeta(patElectron,nobjE,patMuon,nobjM)) +
         (_Electron2Muon2PZetaVisCutCoefficient * CalculatePZetaVis(patElectron,nobjE,patMuon,nobjM))) < _Electron2Muon2CDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoElectron2Muon2DiscrByDeltaPtDivSumPt) {
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedMuonPtEtaPhiMVector.at(nobjM).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _Electron2Muon2DeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt()) / (smearedMuonPtEtaPhiMVector.at(nobjM).pt() + smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _Electron2Muon2DeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoElectron2Muon2DiscrByDeltaPt) {
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) < _Electron2Muon2DeltaPtMinCutValue ) {return false;}
    if ( ((smearedMuonPtEtaPhiMVector.at(nobjM).pt() - smearedElectronPtEtaPhiMVector.at(nobjE).pt())) > _Electron2Muon2DeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByElectron2MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) > _Electron2MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobjE).phi() - theMETVector.phi())) < _Electron2MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByMuon2MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) > _Muon2MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobjM).phi() - theMETVector.phi())) < _Muon2MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByElectron2MetMt) {
    if( (CalculateLeptonMetMt(patElectron,nobjE)<_Electron2MetMtMinCut) || (CalculateLeptonMetMt(patElectron,nobjE)>_Electron2MetMtMaxCut) ) {return false;}
  }
  if (_DoDiscrByMuon2MetMt) {
    if( (CalculateLeptonMetMt(patMuon,nobjM)<_Muon2MetMtMinCut) || (CalculateLeptonMetMt(patMuon,nobjM)>_Muon2MetMtMaxCut) ) {return false;}
  }
  return true;
}
bool HiMassTauAnalysis::passDiMuonTopologyCuts(const pat::Muon& patMuon1, int nobj1, const pat::Muon& patMuon2, int nobj2) {
  // ----Separation cut between muon legs (remove double counting)
  if (_DoDiMuonDiscrByDeltaR) { if(reco::deltaR(smearedMuonMomentumVector.at(nobj1), smearedMuonMomentumVector.at(nobj2)) < _DiMuonDeltaRCut) {return false;} }
  // ---- apply Zll veto cut
  if (_DoMuon1DiscrByIsZllCut) { if((isZmm(smearedMuonMomentumVector.at(nobj1)).first)) {return false;} }
  if (_DoMuon2DiscrByIsZllCut) { if((isZmm(smearedMuonMomentumVector.at(nobj2)).first)) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_DiMuonDiscrByOSLSType == "OS") {
    if((patMuon1.charge() * patMuon2.charge()) >= 0) {return false;}
  } else if (_DiMuonDiscrByOSLSType == "LS") {
    if((patMuon1.charge() * patMuon2.charge()) <= 0) {return false;}
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoDiMuonDiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobj1).phi() - smearedMuonPtEtaPhiMVector.at(nobj2).phi()))) > _DiMuonCosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobj1).phi() - smearedMuonPtEtaPhiMVector.at(nobj2).phi()))) < _DiMuonCosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByDiMuonMassReco) {
    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfDiMuonProductsAndMetMassReco;
    _UseCollinerApproxMassReco = _UseCollinerApproxDiMuonMassReco;
    if( (CalculateThe4Momentum(patMuon1,nobj1,patMuon2,nobj2).second.M() < _DiMuonMassMinCut) || (CalculateThe4Momentum(patMuon1,nobj1,patMuon2,nobj2).second.M() > _DiMuonMassMaxCut) ) {return false;}
  }
  // ----Zeta requirement
  if (_DoDiMuonDiscrByCDFzeta2D) {
    if( ((_DiMuonPZetaCutCoefficient * CalculatePZeta(patMuon1,nobj1,patMuon2,nobj2)) +
         (_DiMuonPZetaVisCutCoefficient * CalculatePZetaVis(patMuon1,nobj1,patMuon2,nobj2))) < _DiMuonCDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoDiMuonDiscrByDeltaPtDivSumPt) {
    if ( ((smearedMuonPtEtaPhiMVector.at(nobj1).pt() - smearedMuonPtEtaPhiMVector.at(nobj2).pt()) / (smearedMuonPtEtaPhiMVector.at(nobj1).pt() + smearedMuonPtEtaPhiMVector.at(nobj2).pt())) < _DiMuonDeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedMuonPtEtaPhiMVector.at(nobj1).pt() - smearedMuonPtEtaPhiMVector.at(nobj2).pt()) / (smearedMuonPtEtaPhiMVector.at(nobj1).pt() + smearedMuonPtEtaPhiMVector.at(nobj2).pt())) > _DiMuonDeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoDiMuonDiscrByDeltaPt) {
    if ( ((smearedMuonPtEtaPhiMVector.at(nobj1).pt() - smearedMuonPtEtaPhiMVector.at(nobj2).pt())) < _DiMuonDeltaPtMinCutValue ) {return false;}
    if ( ((smearedMuonPtEtaPhiMVector.at(nobj1).pt() - smearedMuonPtEtaPhiMVector.at(nobj2).pt())) > _DiMuonDeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByMuon1MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobj1).phi() - theMETVector.phi())) > _Muon1MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobj1).phi() - theMETVector.phi())) < _Muon1MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByMuon2MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobj2).phi() - theMETVector.phi())) > _Muon2MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(nobj2).phi() - theMETVector.phi())) < _Muon2MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByMuon1MetMt) {
    if( (CalculateLeptonMetMt(patMuon1,nobj1)<_Muon1MetMtMinCut) || (CalculateLeptonMetMt(patMuon1,nobj1)>_Muon1MetMtMaxCut) ) {return false;}
  }
  if (_DoDiscrByMuon2MetMt) {
    if( (CalculateLeptonMetMt(patMuon2,nobj2)<_Muon2MetMtMinCut) || (CalculateLeptonMetMt(patMuon2,nobj2)>_Muon2MetMtMaxCut) ) {return false;}
  }
  return true;
}
bool HiMassTauAnalysis::passDiElectronTopologyCuts(const pat::Electron& patElectron1, int nobj1, const pat::Electron& patElectron2, int nobj2) {
  // ----Separation cut between electron legs (remove double counting)
  if (_DoDiElectronDiscrByDeltaR) { if(reco::deltaR(smearedElectronMomentumVector.at(nobj1), smearedElectronMomentumVector.at(nobj2)) < _DiElectronDeltaRCut) {return false;} }
  // ---- apply Zll veto cut
  if (_DoElectron1DiscrByIsZllCut) { if((isZee(smearedMuonMomentumVector.at(nobj1)).first)) {return false;} }
  if (_DoElectron2DiscrByIsZllCut) { if((isZee(smearedMuonMomentumVector.at(nobj2)).first)) {return false;} }
  // ----Opposite sign - Like sign requirement
  if (_DiElectronDiscrByOSLSType == "OS") {
    if((patElectron1.charge() * patElectron2.charge()) >= 0) {return false;}
  } else if (_DiElectronDiscrByOSLSType == "LS") {
    if((patElectron1.charge() * patElectron2.charge()) <= 0) {return false;}
  } else {}
  // ----Require both legs to be almost back-to-back in phi
  if (_DoDiElectronDiscrByCosDphi) {
    if(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobj1).phi() - smearedElectronPtEtaPhiMVector.at(nobj2).phi()))) > _DiElectronCosDphiMaxCut) {return false;}
    if(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobj1).phi() - smearedElectronPtEtaPhiMVector.at(nobj2).phi()))) < _DiElectronCosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByDiElectronMassReco) {
    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfDiElectronProductsAndMetMassReco;
    _UseCollinerApproxMassReco = _UseCollinerApproxDiElectronMassReco;
    if( (CalculateThe4Momentum(patElectron1,nobj1,patElectron2,nobj2).second.M() < _DiElectronMassMinCut) || (CalculateThe4Momentum(patElectron1,nobj1,patElectron2,nobj2).second.M() > _DiElectronMassMaxCut) ) {return false;}
  }
  // ----Zeta requirement
  if (_DoDiElectronDiscrByCDFzeta2D) {
    if( ((_DiElectronPZetaCutCoefficient * CalculatePZeta(patElectron1,nobj1,patElectron2,nobj2)) +
         (_DiElectronPZetaVisCutCoefficient * CalculatePZetaVis(patElectron1,nobj1,patElectron2,nobj2))) < _DiElectronCDFzeta2DCutValue )
    {return false;}
  }
  // ----Missing transverse energy requirement
  if (_DoDiscrByMet) {if(theMETVector.pt() < _RecoMetCut) {return false;}}
  if (_DoDiElectronDiscrByDeltaPtDivSumPt) {
    if ( ((smearedElectronPtEtaPhiMVector.at(nobj1).pt() - smearedElectronPtEtaPhiMVector.at(nobj2).pt()) / (smearedElectronPtEtaPhiMVector.at(nobj1).pt() + smearedElectronPtEtaPhiMVector.at(nobj2).pt())) < _DiElectronDeltaPtDivSumPtMinCutValue ) {return false;}
    if ( ((smearedElectronPtEtaPhiMVector.at(nobj1).pt() - smearedElectronPtEtaPhiMVector.at(nobj2).pt()) / (smearedElectronPtEtaPhiMVector.at(nobj1).pt() + smearedElectronPtEtaPhiMVector.at(nobj2).pt())) > _DiElectronDeltaPtDivSumPtMaxCutValue ) {return false;}
  }
  if (_DoDiElectronDiscrByDeltaPt) {
    if ( ((smearedElectronPtEtaPhiMVector.at(nobj1).pt() - smearedElectronPtEtaPhiMVector.at(nobj2).pt())) < _DiElectronDeltaPtMinCutValue ) {return false;}
    if ( ((smearedElectronPtEtaPhiMVector.at(nobj1).pt() - smearedElectronPtEtaPhiMVector.at(nobj2).pt())) > _DiElectronDeltaPtMaxCutValue ) {return false;}
  }
  if (_DoDiscrByElectron1MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobj1).phi() - theMETVector.phi())) > _Electron1MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobj1).phi() - theMETVector.phi())) < _Electron1MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByElectron2MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobj2).phi() - theMETVector.phi())) > _Electron2MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(nobj2).phi() - theMETVector.phi())) < _Electron2MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByElectron1MetMt) {
    if( (CalculateLeptonMetMt(patElectron1,nobj1)<_Electron1MetMtMinCut) || (CalculateLeptonMetMt(patElectron1,nobj1)>_Electron1MetMtMaxCut) ) {return false;}
  }
  if (_DoDiscrByElectron2MetMt) {
    if( (CalculateLeptonMetMt(patElectron2,nobj2)<_Electron2MetMtMinCut) || (CalculateLeptonMetMt(patElectron2,nobj2)>_Electron2MetMtMaxCut) ) {return false;}
  }
  return true;
}
bool HiMassTauAnalysis::passDiTauTopologyCuts(const pat::Tau& patTau1, int nobj1, const pat::Tau& patTau2, int nobj2) {
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
  if (_DoDiscrByDiTauMassReco) {
    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfDiTauProductsAndMetMassReco;
    _UseCollinerApproxMassReco = _UseCollinerApproxDiTauMassReco;
    if( (CalculateThe4Momentum(patTau1,nobj1,patTau2,nobj2).second.M() < _DiTauMassMinCut) || (CalculateThe4Momentum(patTau1,nobj1,patTau2,nobj2).second.M() > _DiTauMassMaxCut) ) {return false;}
  }
  // ----Zeta requirement
  if (_DoDiTauDiscrByCDFzeta2D) {
    if( ((_DiTauPZetaCutCoefficient * CalculatePZeta(patTau1,nobj1,patTau2,nobj2)) +
         (_DiTauPZetaVisCutCoefficient * CalculatePZetaVis(patTau1,nobj1,patTau2,nobj2))) < _DiTauCDFzeta2DCutValue )
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
  if (_DoDiscrByTau1MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobj1).phi() - theMETVector.phi())) > _Tau1MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobj1).phi() - theMETVector.phi())) < _Tau1MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByTau2MetDphi) {
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobj2).phi() - theMETVector.phi())) > _Tau2MetDphiMaxCut) {return false;}
    if(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(nobj2).phi() - theMETVector.phi())) < _Tau2MetDphiMinCut) {return false;}
  }
  if (_DoDiscrByTau1MetMt) {
    if( (CalculateLeptonMetMt(patTau1,nobj1)<_Tau1MetMtMinCut) || (CalculateLeptonMetMt(patTau1,nobj1)>_Tau1MetMtMaxCut) ) {return false;}
  }
  if (_DoDiscrByTau2MetMt) {
    if( (CalculateLeptonMetMt(patTau2,nobj2)<_Tau2MetMtMinCut) || (CalculateLeptonMetMt(patTau2,nobj2)>_Tau2MetMtMaxCut) ) {return false;}
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
    dphi1 = normalizedPhi(smearedJetPtEtaPhiMVector.at(nobj1).phi() - theMETVector.phi());
    dphi2 = normalizedPhi(smearedJetPtEtaPhiMVector.at(nobj2).phi() - theMETVector.phi());
    r1 = sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0) );
    if(r1 < _R1MinCut) {return false;}
    if(r1 > _R1MaxCut) {return false;}
  }
  if(_DoSUSYDiscrByR2) {
    dphi1 = normalizedPhi(smearedJetPtEtaPhiMVector.at(nobj1).phi() - theMETVector.phi());
    dphi2 = normalizedPhi(smearedJetPtEtaPhiMVector.at(nobj2).phi() - theMETVector.phi());
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
  int numberJets1 = 0;
  for ( pat::JetCollection::const_iterator patJet1 = _patJets->begin();patJet1 != _patJets->end(); ++patJet1 ) {
    numberJets1++;
    if( (numberJets1 == theLeadingJetIndex) && (passRecoFirstLeadingJetCuts((*patJet1),numberJets1 - 1)) ) {
      int numberJets2 = 0;
      for ( pat::JetCollection::const_iterator patJet2 = _patJets->begin();patJet2 != _patJets->end(); ++patJet2 ) {
        numberJets2++;
        if( (numberJets2 == theSecondLeadingJetIndex) && (passRecoSecondLeadingJetCuts((*patJet2),numberJets2 - 1)) ) {
          int theNumberOfJets = 0;
          double sumpxForMht = 0.0;
          double sumpyForMht = 0.0;
          double phiForMht = 0.0;
          for ( pat::JetCollection::const_iterator patJet = _patJets->begin();
            patJet != _patJets->end(); ++patJet ) {
            theNumberOfJets++;
            if (!passRecoJetCuts((*patJet),theNumberOfJets-1)) continue;
            sumpxForMht = sumpxForMht - smearedJetMomentumVector.at(theNumberOfJets-1).px();
            sumpyForMht = sumpyForMht - smearedJetMomentumVector.at(theNumberOfJets-1).py();
          }
          if (sumpxForMht >= 0)
            phiForMht = atan(sumpyForMht/sumpxForMht);
          if (sumpxForMht < 0 && sumpyForMht >= 0)
            phiForMht = atan(sumpyForMht/sumpxForMht) + 3.1415;
          if (sumpxForMht < 0 && sumpyForMht < 0)
            phiForMht = atan(sumpyForMht/sumpxForMht) - 3.1415;
          double dphi1MHT = normalizedPhi(smearedJetPtEtaPhiMVector.at(numberJets1 - 1).phi() - phiForMht);
          double dphi2MHT = normalizedPhi(smearedJetPtEtaPhiMVector.at(numberJets2 - 1).phi() - phiForMht);
          if(_DoSUSYDiscrByDphi1) {
            if(dphi1MHT < _Dphi1MinCut) {return false;}              
            if(dphi1MHT > _Dphi1MaxCut) {return false;}
          }
          if(_DoSUSYDiscrByDphi2) {
            if(abs(dphi2MHT) < _Dphi2MinCut) {return false;}
            if(abs(dphi2MHT) > _Dphi2MaxCut) {return false;}
          }
        }
      }
    }
  }

  return true;
}

// ---------------Fill Ntuple
void HiMassTauAnalysis::fillNtuple() {

/*
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

*/

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
	if (!passRecoTau1Cuts((*patTau),theNumberOfTaus-1)) continue;
	_hTauJet1Energy[NpdfID]->Fill(smearedTauMomentumVector.at(theNumberOfTaus-1).energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJet1Pt[NpdfID]->Fill(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJet1Eta[NpdfID]->Fill(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJet1Phi[NpdfID]->Fill(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
        _UseRecoTauEllipseForEcalIso = _UseRecoTau1EllipseForEcalIso;
        _RecoTauEcalIsoRphiForEllipse = _RecoTau1EcalIsoRphiForEllipse;
        _RecoTauEcalIsoRetaForEllipse = _RecoTau1EcalIsoRetaForEllipse;
        _RecoTauSigGamThreshold = _RecoTau1SigGamThreshold;
        _hTauJet1NumSignalTracks[NpdfID]->Fill(patTau->signalPFChargedHadrCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJet1NumSignalGammas[NpdfID]->Fill(CalculateNumberSignalTauGammas(*patTau),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJet1Charge[NpdfID]->Fill(fabs(patTau->charge()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJet1SignalTracksMass[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJet1SignalTracksAndGammasMass[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJet1SignalTracksChargeFraction[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).pt() / smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	if(patTau->signalPFChargedHadrCands().size() == 1) {
          _UseRecoTauEllipseForEcalIso = _UseRecoTau1EllipseForEcalIso;
          _RecoTauEcalIsoRphiForEllipse = _RecoTau1EcalIsoRphiForEllipse;
          _RecoTauEcalIsoRetaForEllipse = _RecoTau1EcalIsoRetaForEllipse;
          _RecoTauSigGamThreshold = _RecoTau1SigGamThreshold;
	  _hTauJet1SignalTracksMass1prong[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet1SignalTracksAndGammasMass1prong[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet1SignalTracksAndPiZerosMass1prong[NpdfID]->Fill(CalculateTauSignalTracksAndPiZerosMass(*patTau).first.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet1NumSignalPiZeros1prong[NpdfID]->Fill(CalculateTauSignalTracksAndPiZerosMass(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  if(CalculateNumberSignalTauGammas(*patTau) == 0) {_hTauJet1Mass1Prong0Gamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	  if(CalculateNumberSignalTauGammas(*patTau) == 1) {_hTauJet1Mass1Prong1Gamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	  if(CalculateNumberSignalTauGammas(*patTau) >= 2) {_hTauJet1Mass1Prong2orMoreGamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	}
	if(patTau->signalPFChargedHadrCands().size() == 3) {
          _UseRecoTauEllipseForEcalIso = _UseRecoTau1EllipseForEcalIso;
          _RecoTauEcalIsoRphiForEllipse = _RecoTau1EcalIsoRphiForEllipse;
          _RecoTauEcalIsoRetaForEllipse = _RecoTau1EcalIsoRetaForEllipse;
          _RecoTauSigGamThreshold = _RecoTau1SigGamThreshold;
	  _hTauJet1SignalTracksMass3prong[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet1SignalTracksAndGammasMass3prong[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  if(CalculateNumberSignalTauGammas(*patTau) == 0) {_hTauJet1Mass3Prong0Gamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	  if(CalculateNumberSignalTauGammas(*patTau) == 1) {_hTauJet1Mass3Prong1Gamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	  if(CalculateNumberSignalTauGammas(*patTau) >= 2) {_hTauJet1Mass3Prong2orMoreGamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	}
	if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
	  _hTauJet1SeedTrackPt[NpdfID]->Fill(patTau->leadPFChargedHadrCand()->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet1SeedTrackIpSignificance[NpdfID]->Fill(patTau->leadPFChargedHadrCandsignedSipt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  if( (patTau->leadPFChargedHadrCand()->trackRef().isNonnull()) ) {
	    _hTauJet1SeedTrackNhits[NpdfID]->Fill(patTau->leadPFChargedHadrCand()->trackRef()->numberOfValidHits(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJet1SeedTrackChi2[NpdfID]->Fill(patTau->leadPFChargedHadrCand()->trackRef()->chi2(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          }
	  _hTauJet1H3x3OverP[NpdfID]->Fill(patTau->hcal3x3OverPLead(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	}
        if (_UseRecoTau1DiscrByIsolationFlag) {
	  _hTauJet1NumIsoTracks[NpdfID]->Fill(patTau->isolationPFChargedHadrCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet1NumIsoGammas[NpdfID]->Fill(patTau->isolationPFGammaCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet1NumIsoCands[NpdfID]->Fill(patTau->isolationPFChargedHadrCands().size() + patTau->isolationPFGammaCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet1SumPtIsoTracks[NpdfID]->Fill(patTau->isolationPFChargedHadrCandsPtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet1SumPtIsoGammas[NpdfID]->Fill(patTau->isolationPFGammaCandsEtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet1SumPtIso[NpdfID]->Fill(patTau->isolationPFChargedHadrCandsPtSum() + patTau->isolationPFGammaCandsEtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        } else {
          _RecoTauIsoDeltaRCone = _RecoTau1IsoDeltaRCone;
          _RecoTauTrackIsoTrkThreshold = _RecoTau1TrackIsoTrkThreshold;
          _RecoTauGammaIsoGamThreshold = _RecoTau1GammaIsoGamThreshold;
          _UseRecoTauEllipseForEcalIso = _UseRecoTau1EllipseForEcalIso;
          _RecoTauEcalIsoRphiForEllipse = _RecoTau1EcalIsoRphiForEllipse;
          _RecoTauEcalIsoRetaForEllipse = _RecoTau1EcalIsoRetaForEllipse;
	  _hTauJet1NumIsoTracks[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet1NumIsoGammas[NpdfID]->Fill(CalculateTauEcalIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet1NumIsoCands[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).first + CalculateTauEcalIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet1SumPtIsoTracks[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet1SumPtIsoGammas[NpdfID]->Fill(CalculateTauEcalIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet1SumPtIso[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).second + CalculateTauEcalIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          if (patTau->leadPFChargedHadrCand().isNonnull()) {
            for( unsigned ipfcand=0; ipfcand < _pflow->size(); ++ipfcand ) {
              const reco::PFCandidate& cand = (*_pflow)[ipfcand];
              if(cand.trackRef().isNonnull()) {
                if((cand.pt()>1.) && (cand.particleId()==1) && (cand.trackRef()->numberOfValidHits() > 8) && 
                   (cand.trackRef().isNonnull()) && (fabs(cand.trackRef()->dxy(thePrimaryEventVertex.position())) < 0.03) && 
                   (fabs(cand.trackRef()->dz(thePrimaryEventVertex.position())) < 0.2) ) {
   	          _hTauJet1NumberDensity[NpdfID]->Fill(reco::deltaR(cand.eta(),cand.phi(),patTau->leadPFChargedHadrCand()->eta(),patTau->leadPFChargedHadrCand()->phi()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
                }
              }
            }
          }
        }
	if((_GenParticleSource.label() != "") && (!isData)) {
	  if(matchToGen(*patTau).first) {
	    _hTauJet1GenTauDeltaPhi[NpdfID]->Fill(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi() - matchToGen(*patTau).second.phi()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJet1GenTauDeltaEta[NpdfID]->Fill(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).eta() - matchToGen(*patTau).second.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJet1GenTauDeltaPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - matchToGen(*patTau).second.pt()) / matchToGen(*patTau).second.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
        nTaus++;
      }
      _hNTau1[NpdfID]->Fill(nTaus,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
      nTaus = 0;
      theNumberOfTaus = 0;
      for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
	theNumberOfTaus++;
	if (!passRecoTau2Cuts((*patTau),theNumberOfTaus-1)) continue;
	_hTauJet2Energy[NpdfID]->Fill(smearedTauMomentumVector.at(theNumberOfTaus-1).energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJet2Pt[NpdfID]->Fill(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJet2Eta[NpdfID]->Fill(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJet2Phi[NpdfID]->Fill(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
        _UseRecoTauEllipseForEcalIso = _UseRecoTau2EllipseForEcalIso;
        _RecoTauEcalIsoRphiForEllipse = _RecoTau2EcalIsoRphiForEllipse;
        _RecoTauEcalIsoRetaForEllipse = _RecoTau2EcalIsoRetaForEllipse;
        _RecoTauSigGamThreshold = _RecoTau2SigGamThreshold;
        _hTauJet2NumSignalTracks[NpdfID]->Fill(patTau->signalPFChargedHadrCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJet2NumSignalGammas[NpdfID]->Fill(CalculateNumberSignalTauGammas(*patTau),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJet2Charge[NpdfID]->Fill(fabs(patTau->charge()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJet2SignalTracksMass[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJet2SignalTracksAndGammasMass[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hTauJet2SignalTracksChargeFraction[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).pt() / smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	if(patTau->signalPFChargedHadrCands().size() == 1) {
          _UseRecoTauEllipseForEcalIso = _UseRecoTau2EllipseForEcalIso;
          _RecoTauEcalIsoRphiForEllipse = _RecoTau2EcalIsoRphiForEllipse;
          _RecoTauEcalIsoRetaForEllipse = _RecoTau2EcalIsoRetaForEllipse;
          _RecoTauSigGamThreshold = _RecoTau2SigGamThreshold;
	  _hTauJet2SignalTracksMass1prong[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet2SignalTracksAndGammasMass1prong[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet2SignalTracksAndPiZerosMass1prong[NpdfID]->Fill(CalculateTauSignalTracksAndPiZerosMass(*patTau).first.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet2NumSignalPiZeros1prong[NpdfID]->Fill(CalculateTauSignalTracksAndPiZerosMass(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  if(CalculateNumberSignalTauGammas(*patTau) == 0) {_hTauJet2Mass1Prong0Gamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	  if(CalculateNumberSignalTauGammas(*patTau) == 1) {_hTauJet2Mass1Prong1Gamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	  if(CalculateNumberSignalTauGammas(*patTau) >= 2) {_hTauJet2Mass1Prong2orMoreGamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	}
	if(patTau->signalPFChargedHadrCands().size() == 3) {
          _UseRecoTauEllipseForEcalIso = _UseRecoTau2EllipseForEcalIso;
          _RecoTauEcalIsoRphiForEllipse = _RecoTau2EcalIsoRphiForEllipse;
          _RecoTauEcalIsoRetaForEllipse = _RecoTau2EcalIsoRetaForEllipse;
          _RecoTauSigGamThreshold = _RecoTau2SigGamThreshold;
	  _hTauJet2SignalTracksMass3prong[NpdfID]->Fill(CalculateTauSignalTracksMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet2SignalTracksAndGammasMass3prong[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  if(CalculateNumberSignalTauGammas(*patTau) == 0) {_hTauJet2Mass3Prong0Gamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	  if(CalculateNumberSignalTauGammas(*patTau) == 1) {_hTauJet2Mass3Prong1Gamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	  if(CalculateNumberSignalTauGammas(*patTau) >= 2) {_hTauJet2Mass3Prong2orMoreGamma[NpdfID]->Fill(CalculateTauSignalTracksAndGammasMass(*patTau).M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	}
	if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
	  _hTauJet2SeedTrackPt[NpdfID]->Fill(patTau->leadPFChargedHadrCand()->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet2SeedTrackIpSignificance[NpdfID]->Fill(patTau->leadPFChargedHadrCandsignedSipt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  if( (patTau->leadPFChargedHadrCand()->trackRef().isNonnull()) ) {
	    _hTauJet2SeedTrackNhits[NpdfID]->Fill(patTau->leadPFChargedHadrCand()->trackRef()->numberOfValidHits(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJet2SeedTrackChi2[NpdfID]->Fill(patTau->leadPFChargedHadrCand()->trackRef()->chi2(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          }
	  _hTauJet2H3x3OverP[NpdfID]->Fill(patTau->hcal3x3OverPLead(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	}
        if (_UseRecoTau2DiscrByIsolationFlag) {
	  _hTauJet2NumIsoTracks[NpdfID]->Fill(patTau->isolationPFChargedHadrCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet2NumIsoGammas[NpdfID]->Fill(patTau->isolationPFGammaCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet2NumIsoCands[NpdfID]->Fill(patTau->isolationPFChargedHadrCands().size() + patTau->isolationPFGammaCands().size(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet2SumPtIsoTracks[NpdfID]->Fill(patTau->isolationPFChargedHadrCandsPtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet2SumPtIsoGammas[NpdfID]->Fill(patTau->isolationPFGammaCandsEtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet2SumPtIso[NpdfID]->Fill(patTau->isolationPFChargedHadrCandsPtSum() + patTau->isolationPFGammaCandsEtSum(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        } else {
          _RecoTauIsoDeltaRCone = _RecoTau2IsoDeltaRCone;
          _RecoTauTrackIsoTrkThreshold = _RecoTau2TrackIsoTrkThreshold;
          _RecoTauGammaIsoGamThreshold = _RecoTau2GammaIsoGamThreshold;
          _UseRecoTauEllipseForEcalIso = _UseRecoTau2EllipseForEcalIso;
          _RecoTauEcalIsoRphiForEllipse = _RecoTau2EcalIsoRphiForEllipse;
          _RecoTauEcalIsoRetaForEllipse = _RecoTau2EcalIsoRetaForEllipse;
	  _hTauJet2NumIsoTracks[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet2NumIsoGammas[NpdfID]->Fill(CalculateTauEcalIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet2NumIsoCands[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).first + CalculateTauEcalIsolation(*patTau).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet2SumPtIsoTracks[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet2SumPtIsoGammas[NpdfID]->Fill(CalculateTauEcalIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  _hTauJet2SumPtIso[NpdfID]->Fill(CalculateTauTrackIsolation(*patTau).second + CalculateTauEcalIsolation(*patTau).second,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          if (patTau->leadPFChargedHadrCand().isNonnull()) {
            for( unsigned ipfcand=0; ipfcand < _pflow->size(); ++ipfcand ) {
              const reco::PFCandidate& cand = (*_pflow)[ipfcand];
              if(cand.trackRef().isNonnull()) {
                if((cand.pt()>1.) && (cand.particleId()==1) && (cand.trackRef()->numberOfValidHits() > 8) && 
                   (cand.trackRef().isNonnull()) && (fabs(cand.trackRef()->dxy(thePrimaryEventVertex.position())) < 0.03) && 
                   (fabs(cand.trackRef()->dz(thePrimaryEventVertex.position())) < 0.2) ) {
   	          _hTauJet2NumberDensity[NpdfID]->Fill(reco::deltaR(cand.eta(),cand.phi(),patTau->leadPFChargedHadrCand()->eta(),patTau->leadPFChargedHadrCand()->phi()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
                }
              }
            }
          }
        }
	if((_GenParticleSource.label() != "") && (!isData)) {
	  if(matchToGen(*patTau).first) {
	    _hTauJet2GenTauDeltaPhi[NpdfID]->Fill(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi() - matchToGen(*patTau).second.phi()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJet2GenTauDeltaEta[NpdfID]->Fill(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).eta() - matchToGen(*patTau).second.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTauJet2GenTauDeltaPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - matchToGen(*patTau).second.pt()) / matchToGen(*patTau).second.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
        nTaus++;
      }
      _hNTau2[NpdfID]->Fill(nTaus,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
    }

    // ------Reco Muon Histograms
    if (_FillRecoMuonHists) {
      int nMuons = 0;
      int theNumberOfMuons = 0;
      for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); 
	    patMuon != _patMuons->end(); ++patMuon ) {
	theNumberOfMuons++;
	if (!passRecoMuon1Cuts((*patMuon),theNumberOfMuons-1)) continue;
        _hMuon1Energy[NpdfID]->Fill(smearedMuonMomentumVector.at(theNumberOfMuons-1).energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon1Pt[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon1Eta[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon1Phi[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	if((_GenParticleSource.label() != "") && (!isData)) {
	  if(matchToGen(*patMuon).first) {
	    _hMuon1GenMuonDeltaPhi[NpdfID]->Fill(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - matchToGen(*patMuon).second.phi()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1GenMuonDeltaEta[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).eta() - matchToGen(*patMuon).second.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1GenMuonDeltaPt[NpdfID]->Fill((smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt() - matchToGen(*patMuon).second.pt()) / matchToGen(*patMuon).second.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
	_hMuon1TrackIso[NpdfID]->Fill(patMuon->trackIsoDeposit()->depositAndCountWithin(_RecoMuon1IsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuon1TrackIsoTrkThreshold).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon1EcalIso[NpdfID]->Fill(patMuon->ecalIsoDeposit()->depositAndCountWithin(_RecoMuon1IsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuon1EcalIsoRecHitThreshold).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon1Iso[NpdfID]->Fill(patMuon->trackIsoDeposit()->depositAndCountWithin(_RecoMuon1IsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuon1TrackIsoTrkThreshold).first +
				 patMuon->ecalIsoDeposit()->depositAndCountWithin(_RecoMuon1IsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuon1EcalIsoRecHitThreshold).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon1CaloCompatibility[NpdfID]->Fill(muon::caloCompatibility(*patMuon),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon1SegmentCompatibility[NpdfID]->Fill(muon::segmentCompatibility(*patMuon),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon1CaloCompatibilityVsSegmentCompatibility[NpdfID]->Fill(muon::caloCompatibility(*patMuon),muon::segmentCompatibility(*patMuon),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon1AntiPion[NpdfID]->Fill((_RecoMuon1CaloCompCoefficient * muon::caloCompatibility(*patMuon)) + (_RecoMuon1SegmCompCoefficient * muon::segmentCompatibility(*patMuon)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
	if ( patMuon->track().isNonnull() ) {
	  _hMuon1Ip[NpdfID]->Fill( patMuon->track()->dxy(thePrimaryEventVertex.position()) ,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  if(fabs(patMuon->track()->dxyError()) != 0) {
	    _hMuon1IpSignificance[NpdfID]->Fill( fabs(patMuon->track()->dxy(thePrimaryEventVertex.position()) / patMuon->track()->dxyError()) ,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  } else {_hMuon1IpSignificance[NpdfID]->Fill(-1,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	}
	nMuons++;
      }
      _hNMuon1[NpdfID]->Fill(nMuons,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
      nMuons = 0;
      theNumberOfMuons = 0;
      for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); 
	    patMuon != _patMuons->end(); ++patMuon ) {
	theNumberOfMuons++;
	if (!passRecoMuon2Cuts((*patMuon),theNumberOfMuons-1)) continue;
        _hMuon2Energy[NpdfID]->Fill(smearedMuonMomentumVector.at(theNumberOfMuons-1).energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon2Pt[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon2Eta[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon2Phi[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	if((_GenParticleSource.label() != "") && (!isData)) {
	  if(matchToGen(*patMuon).first) {
	    _hMuon2GenMuonDeltaPhi[NpdfID]->Fill(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - matchToGen(*patMuon).second.phi()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2GenMuonDeltaEta[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).eta() - matchToGen(*patMuon).second.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2GenMuonDeltaPt[NpdfID]->Fill((smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt() - matchToGen(*patMuon).second.pt()) / matchToGen(*patMuon).second.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
	_hMuon2TrackIso[NpdfID]->Fill(patMuon->trackIsoDeposit()->depositAndCountWithin(_RecoMuon2IsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuon2TrackIsoTrkThreshold).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon2EcalIso[NpdfID]->Fill(patMuon->ecalIsoDeposit()->depositAndCountWithin(_RecoMuon2IsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuon2EcalIsoRecHitThreshold).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon2Iso[NpdfID]->Fill(patMuon->trackIsoDeposit()->depositAndCountWithin(_RecoMuon2IsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuon2TrackIsoTrkThreshold).first +
				 patMuon->ecalIsoDeposit()->depositAndCountWithin(_RecoMuon2IsoDeltaRCone,reco::IsoDeposit::Vetos(),_RecoMuon2EcalIsoRecHitThreshold).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon2CaloCompatibility[NpdfID]->Fill(muon::caloCompatibility(*patMuon),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon2SegmentCompatibility[NpdfID]->Fill(muon::segmentCompatibility(*patMuon),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon2CaloCompatibilityVsSegmentCompatibility[NpdfID]->Fill(muon::caloCompatibility(*patMuon),muon::segmentCompatibility(*patMuon),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hMuon2AntiPion[NpdfID]->Fill((_RecoMuon2CaloCompCoefficient * muon::caloCompatibility(*patMuon)) + (_RecoMuon2SegmCompCoefficient * muon::segmentCompatibility(*patMuon)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
	if ( patMuon->track().isNonnull() ) {
	  _hMuon2Ip[NpdfID]->Fill( patMuon->track()->dxy(thePrimaryEventVertex.position()) ,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  if(fabs(patMuon->track()->dxyError()) != 0) {
	    _hMuon2IpSignificance[NpdfID]->Fill( fabs(patMuon->track()->dxy(thePrimaryEventVertex.position()) / patMuon->track()->dxyError()) ,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  } else {_hMuon2IpSignificance[NpdfID]->Fill(-1,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	}
	nMuons++;
      }
      _hNMuon2[NpdfID]->Fill(nMuons,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
    }

    //-----Fill Reco Electron Histograms
    if (_FillRecoElectronHists) {
      int nElectrons = 0;
      int theNumberOfElectrons = 0;
      for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();
	    patElectron != _patElectrons->end(); ++patElectron ) {
	theNumberOfElectrons++;
	if (!passRecoElectron1Cuts((*patElectron),theNumberOfElectrons-1)) continue;
	_hElectron1Energy[NpdfID]->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectron1Pt[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectron1Eta[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectron1Phi[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	if((_GenParticleSource.label() != "") && (!isData)) {
	  if(matchToGen(*patElectron).first) {
	    _hElectron1GenElectronDeltaPhi[NpdfID]->Fill(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - matchToGen(*patElectron).second.phi()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1GenElectronDeltaEta[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).eta() - matchToGen(*patElectron).second.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1GenElectronDeltaPt[NpdfID]->Fill((smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt() - matchToGen(*patElectron).second.pt()) / matchToGen(*patElectron).second.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
	_hElectron1TrackIso[NpdfID]->Fill(patElectron->dr04TkSumPt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectron1EcalIso[NpdfID]->Fill(patElectron->dr04EcalRecHitSumEt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
	if ( patElectron->track().isNonnull() ) { _hElectron1Ip[NpdfID]->Fill( patElectron->track()->dxy(thePrimaryEventVertex.position()) ,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID)); }
	//if ( patElectron->gsfTrack().isNonnull() ) { _hElectron1Ip[NpdfID]->Fill( patElectron->gsfTrack()->dxy(thePrimaryEventVertex.position()) ,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID)); }
	_hElectron1EoverP[NpdfID]->Fill(patElectron->eSuperClusterOverP(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        if(_UseHeepInfo) {
          heep::Ele theHeepElec(*patElectron);
	  _hElectron1HoverEm[NpdfID]->Fill(theHeepElec.hOverE(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          if(patElectron->isEE()) {
	    _hElectron1EESigmaIEtaIEta[NpdfID]->Fill(theHeepElec.sigmaIEtaIEta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1EEDEta[NpdfID]->Fill(theHeepElec.dEtaIn(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1EEDPhi[NpdfID]->Fill(theHeepElec.dPhiIn(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          }
          if(patElectron->isEB()) {
	    _hElectron1EBSigmaIEtaIEta[NpdfID]->Fill(theHeepElec.sigmaIEtaIEta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1EBDEta[NpdfID]->Fill(theHeepElec.dEtaIn(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1EBDPhi[NpdfID]->Fill(theHeepElec.dPhiIn(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1EB2by5Over5by5[NpdfID]->Fill(theHeepElec.scE2x5Max() / theHeepElec.scE5x5(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1EB1by5Over5by5[NpdfID]->Fill(theHeepElec.scE1x5() / theHeepElec.scE5x5(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          }
        } else {
	  _hElectron1HoverEm[NpdfID]->Fill(patElectron->hadronicOverEm(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          if(patElectron->isEE()) {
	    _hElectron1EESigmaIEtaIEta[NpdfID]->Fill(patElectron->scSigmaIEtaIEta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1EEDEta[NpdfID]->Fill(patElectron->deltaEtaSuperClusterTrackAtVtx(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1EEDPhi[NpdfID]->Fill(patElectron->deltaPhiSuperClusterTrackAtVtx(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          }
          if(patElectron->isEB()) {
	    _hElectron1EBSigmaIEtaIEta[NpdfID]->Fill(patElectron->scSigmaIEtaIEta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1EBDEta[NpdfID]->Fill(patElectron->deltaEtaSuperClusterTrackAtVtx(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1EBDPhi[NpdfID]->Fill(patElectron->deltaPhiSuperClusterTrackAtVtx(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1EB2by5Over5by5[NpdfID]->Fill(patElectron->scE2x5Max() / patElectron->scE5x5(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1EB1by5Over5by5[NpdfID]->Fill(patElectron->scE1x5() / patElectron->scE5x5(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          }
        }
        const Track* elTrack = (const reco::Track*)(patElectron->gsfTrack().get());
        const HitPattern& pInner = elTrack->trackerExpectedHitsInner();
	_hElectron1MissingHits[NpdfID]->Fill(pInner.numberOfHits(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectron1Classification[NpdfID]->Fill(patElectron->classification(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectron1EcalDriven[NpdfID]->Fill(patElectron->ecalDrivenSeed(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectron1TrackerDriven[NpdfID]->Fill(patElectron->trackerDrivenSeed(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	nElectrons++;
      }
      _hNElectron1[NpdfID]->Fill(nElectrons,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
      nElectrons = 0;
      theNumberOfElectrons = 0;
      for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();
	    patElectron != _patElectrons->end(); ++patElectron ) {
	theNumberOfElectrons++;
	if (!passRecoElectron2Cuts((*patElectron),theNumberOfElectrons-1)) continue;
	_hElectron2Energy[NpdfID]->Fill(smearedElectronMomentumVector.at(theNumberOfElectrons-1).energy(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectron2Pt[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectron2Eta[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectron2Phi[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	if((_GenParticleSource.label() != "") && (!isData)) {
	  if(matchToGen(*patElectron).first) {
	    _hElectron2GenElectronDeltaPhi[NpdfID]->Fill(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - matchToGen(*patElectron).second.phi()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2GenElectronDeltaEta[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).eta() - matchToGen(*patElectron).second.eta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2GenElectronDeltaPt[NpdfID]->Fill((smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt() - matchToGen(*patElectron).second.pt()) / matchToGen(*patElectron).second.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
	_hElectron2TrackIso[NpdfID]->Fill(patElectron->dr04TkSumPt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectron2EcalIso[NpdfID]->Fill(patElectron->dr04EcalRecHitSumEt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	const reco::Vertex& thePrimaryEventVertex = (*(_primaryEventVertexCollection)->begin());
	if ( patElectron->track().isNonnull() ) { _hElectron2Ip[NpdfID]->Fill( patElectron->track()->dxy(thePrimaryEventVertex.position()) ,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID)); }
	//if ( patElectron->gsfTrack().isNonnull() ) { _hElectron2Ip[NpdfID]->Fill( patElectron->gsfTrack()->dxy(thePrimaryEventVertex.position()) ,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID)); }
	_hElectron2EoverP[NpdfID]->Fill(patElectron->eSuperClusterOverP(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        if(_UseHeepInfo) {
          heep::Ele theHeepElec(*patElectron);
	  _hElectron2HoverEm[NpdfID]->Fill(theHeepElec.hOverE(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          if(patElectron->isEE()) {
	    _hElectron2EESigmaIEtaIEta[NpdfID]->Fill(theHeepElec.sigmaIEtaIEta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2EEDEta[NpdfID]->Fill(theHeepElec.dEtaIn(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2EEDPhi[NpdfID]->Fill(theHeepElec.dPhiIn(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          }
          if(patElectron->isEB()) {
	    _hElectron2EBSigmaIEtaIEta[NpdfID]->Fill(theHeepElec.sigmaIEtaIEta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2EBDEta[NpdfID]->Fill(theHeepElec.dEtaIn(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2EBDPhi[NpdfID]->Fill(theHeepElec.dPhiIn(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2EB2by5Over5by5[NpdfID]->Fill(theHeepElec.scE2x5Max() / theHeepElec.scE5x5(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2EB1by5Over5by5[NpdfID]->Fill(theHeepElec.scE1x5() / theHeepElec.scE5x5(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          }
        } else {
	  _hElectron2HoverEm[NpdfID]->Fill(patElectron->hadronicOverEm(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          if(patElectron->isEE()) {
	    _hElectron2EESigmaIEtaIEta[NpdfID]->Fill(patElectron->scSigmaIEtaIEta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2EEDEta[NpdfID]->Fill(patElectron->deltaEtaSuperClusterTrackAtVtx(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2EEDPhi[NpdfID]->Fill(patElectron->deltaPhiSuperClusterTrackAtVtx(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          }
          if(patElectron->isEB()) {
	    _hElectron2EBSigmaIEtaIEta[NpdfID]->Fill(patElectron->scSigmaIEtaIEta(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2EBDEta[NpdfID]->Fill(patElectron->deltaEtaSuperClusterTrackAtVtx(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2EBDPhi[NpdfID]->Fill(patElectron->deltaPhiSuperClusterTrackAtVtx(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2EB2by5Over5by5[NpdfID]->Fill(patElectron->scE2x5Max() / patElectron->scE5x5(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2EB1by5Over5by5[NpdfID]->Fill(patElectron->scE1x5() / patElectron->scE5x5(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
          }
        }
        const Track* elTrack = (const reco::Track*)(patElectron->gsfTrack().get());
        const HitPattern& pInner = elTrack->trackerExpectedHitsInner();
	_hElectron2MissingHits[NpdfID]->Fill(pInner.numberOfHits(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectron2Classification[NpdfID]->Fill(patElectron->classification(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectron2EcalDriven[NpdfID]->Fill(patElectron->ecalDrivenSeed(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hElectron2TrackerDriven[NpdfID]->Fill(patElectron->trackerDrivenSeed(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	nElectrons++;
      }
      _hNElectron2[NpdfID]->Fill(nElectrons,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
    }

    // ------Reco Jet Histograms
    if (_FillRecoJetHists) {
      int NbJets = 0;
      int theNumberOfJets = 0;
      for ( pat::JetCollection::const_iterator patJet = _patJets->begin(); 
	    patJet != _patJets->end(); ++patJet ) {
	theNumberOfJets++;
	if (!passRecoBJetCuts((*patJet),theNumberOfJets-1)) continue;
	_hBJetDiscrByTrackCounting[NpdfID]->Fill(patJet->bDiscriminator("trackCountingHighEffBJetTags"),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBJetDiscrBySimpleSecondaryV[NpdfID]->Fill(patJet->bDiscriminator("simpleSecondaryVertexBJetTags"),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	_hBJetDiscrByCombinedSecondaryV[NpdfID]->Fill(patJet->bDiscriminator("combinedSecondaryVertexBJetTags"),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
        NbJets++;
      }
      _hNBJet[NpdfID]->Fill(NbJets,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));

      int nJets = 0;
      theNumberOfJets = 0;
      for ( pat::JetCollection::const_iterator patJet = _patJets->begin(); 
	    patJet != _patJets->end(); ++patJet ) {
	theNumberOfJets++;
	if (!passRecoJetCuts((*patJet),theNumberOfJets-1)) continue;
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
    if (_FillTopologyHists) {
      _hMet[NpdfID]->Fill(theMETVector.pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
      if(!isData) {
        const GenMETCollection *genmetcol = genTrue.product();
        const GenMET *genMetTrue = &(genmetcol->front());
        _hMetResolution[NpdfID]->Fill((theMETVector.pt() - genMetTrue->pt()) / genMetTrue->pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
      }
      int theNumberOfMuons = 0;
      for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
	theNumberOfMuons++;
	int theNumberOfTaus = 0;
	for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
	  theNumberOfTaus++;
	  if ((passRecoTau1Cuts((*patTau),theNumberOfTaus-1)) &&
	      (passRecoMuon1Cuts((*patMuon),theNumberOfMuons-1)) &&
	      (passMuon1Tau1TopologyCuts((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1))) {
	    _hMuon1PtVsTau1Pt[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt(),smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau1DeltaR[NpdfID]->Fill(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), smearedMuonMomentumVector.at(theNumberOfMuons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau1DeltaPtDivSumPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt()) / (smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() + smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau1DeltaPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau1CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau1_Muon1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1MetDeltaPhiVsMuon1Tau1CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - theMETVector.phi())), cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau1_Tau1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
            _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon1Tau1ProductsAndMetMassReco;
            _UseCollinerApproxMassReco = _UseCollinerApproxMuon1Tau1MassReco;
	    if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon1Tau1ReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    else {_hMuon1Tau1NotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    _hMuon1Tau1_Muon1MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau1_Tau1MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau),theNumberOfTaus-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    if(_UseTauSeedTrackForMuon1Tau1DiscrByOSLS) {
	      if (patTau->isCaloTau()) {
		if( (patTau->leadTrack().isNonnull()) ) {
		  _hMuon1Tau1OSLS[NpdfID]->Fill(patMuon->charge() * patTau->leadTrack()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		}
              } else {
		if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
		  _hMuon1Tau1OSLS[NpdfID]->Fill(patMuon->charge() * patTau->leadPFChargedHadrCand()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
                  if((patTau->leadPFChargedHadrCand()->charge() * patMuon->charge())<0) {
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon1Tau1ProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxMuon1Tau1MassReco;
   	            if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon1Tau1ReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hMuon1Tau1NotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  } else {
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon1Tau1ProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxMuon1Tau1MassReco;
   	            if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon1Tau1ReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hMuon1Tau1NotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  }
		}
              }
	    } else {
              _hMuon1Tau1OSLS[NpdfID]->Fill(patMuon->charge() * patTau->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              if((patTau->charge() * patMuon->charge())<0) {
                _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon1Tau1ProductsAndMetMassReco;
                _UseCollinerApproxMassReco = _UseCollinerApproxMuon1Tau1MassReco;
                if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon1Tau1ReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	        else {_hMuon1Tau1NotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
              } else {
                _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon1Tau1ProductsAndMetMassReco;
                _UseCollinerApproxMassReco = _UseCollinerApproxMuon1Tau1MassReco;
   	        if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon1Tau1ReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	        else {_hMuon1Tau1NotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
              }
            }
	    _hMuon1Tau1PZeta[NpdfID]->Fill(CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau1PZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau1Zeta2D[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau1Zeta1D[NpdfID]->Fill((_Muon1Tau1PZetaCutCoefficient * CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1)) + 
				     (_Muon1Tau1PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
      }
      theNumberOfMuons = 0;
      for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
	theNumberOfMuons++;
	int theNumberOfTaus = 0;
	for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
	  theNumberOfTaus++;
	  if ((passRecoTau2Cuts((*patTau),theNumberOfTaus-1)) &&
	      (passRecoMuon1Cuts((*patMuon),theNumberOfMuons-1)) &&
	      (passMuon1Tau2TopologyCuts((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1))) {
	    _hMuon1PtVsTau2Pt[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt(),smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau2DeltaR[NpdfID]->Fill(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), smearedMuonMomentumVector.at(theNumberOfMuons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau2DeltaPtDivSumPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt()) / (smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() + smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau2DeltaPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau2CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau2_Muon1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1MetDeltaPhiVsMuon1Tau2CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - theMETVector.phi())), cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau2_Tau2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
            _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon1Tau2ProductsAndMetMassReco;
            _UseCollinerApproxMassReco = _UseCollinerApproxMuon1Tau2MassReco;
	    if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon1Tau2ReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    else {_hMuon1Tau2NotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    _hMuon1Tau2_Muon1MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau2_Tau2MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau),theNumberOfTaus-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    if(_UseTauSeedTrackForMuon1Tau2DiscrByOSLS) {
	      if (patTau->isCaloTau()) {
		if( (patTau->leadTrack().isNonnull()) ) {
		  _hMuon1Tau2OSLS[NpdfID]->Fill(patMuon->charge() * patTau->leadTrack()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		}
              } else {
		if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
		  _hMuon1Tau2OSLS[NpdfID]->Fill(patMuon->charge() * patTau->leadPFChargedHadrCand()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
                  if((patTau->leadPFChargedHadrCand()->charge() * patMuon->charge())<0) {
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon1Tau2ProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxMuon1Tau2MassReco;
   	            if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon1Tau2ReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hMuon1Tau2NotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  } else {
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon1Tau2ProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxMuon1Tau2MassReco;
   	            if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon1Tau2ReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hMuon1Tau2NotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  }
		}
              }
	    } else {
              _hMuon1Tau2OSLS[NpdfID]->Fill(patMuon->charge() * patTau->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              if((patTau->charge() * patMuon->charge())<0) {
                _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon1Tau2ProductsAndMetMassReco;
                _UseCollinerApproxMassReco = _UseCollinerApproxMuon1Tau2MassReco;
                if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon1Tau2ReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	        else {_hMuon1Tau2NotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
              } else {
                _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon1Tau2ProductsAndMetMassReco;
                _UseCollinerApproxMassReco = _UseCollinerApproxMuon1Tau2MassReco;
   	        if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon1Tau2ReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	        else {_hMuon1Tau2NotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
              }
            }
	    _hMuon1Tau2PZeta[NpdfID]->Fill(CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau2PZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau2Zeta2D[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Tau2Zeta1D[NpdfID]->Fill((_Muon1Tau2PZetaCutCoefficient * CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1)) + 
				     (_Muon1Tau2PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
      }
      theNumberOfMuons = 0;
      for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
	theNumberOfMuons++;
	int theNumberOfTaus = 0;
	for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
	  theNumberOfTaus++;
	  if ((passRecoTau1Cuts((*patTau),theNumberOfTaus-1)) &&
	      (passRecoMuon2Cuts((*patMuon),theNumberOfMuons-1)) &&
	      (passMuon2Tau1TopologyCuts((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1))) {
	    _hMuon2PtVsTau1Pt[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt(),smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau1DeltaR[NpdfID]->Fill(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), smearedMuonMomentumVector.at(theNumberOfMuons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau1DeltaPtDivSumPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt()) / (smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() + smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau1DeltaPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau1CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau1_Muon2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2MetDeltaPhiVsMuon2Tau1CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - theMETVector.phi())), cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau1_Tau1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
            _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon2Tau1ProductsAndMetMassReco;
            _UseCollinerApproxMassReco = _UseCollinerApproxMuon2Tau1MassReco;
	    if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon2Tau1ReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    else {_hMuon2Tau1NotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    _hMuon2Tau1_Muon2MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau1_Tau1MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau),theNumberOfTaus-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    if(_UseTauSeedTrackForMuon2Tau1DiscrByOSLS) {
	      if (patTau->isCaloTau()) {
		if( (patTau->leadTrack().isNonnull()) ) {
		  _hMuon2Tau1OSLS[NpdfID]->Fill(patMuon->charge() * patTau->leadTrack()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		}
              } else {
		if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
		  _hMuon2Tau1OSLS[NpdfID]->Fill(patMuon->charge() * patTau->leadPFChargedHadrCand()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
                  if((patTau->leadPFChargedHadrCand()->charge() * patMuon->charge())<0) {
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon2Tau1ProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxMuon2Tau1MassReco;
   	            if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon2Tau1ReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hMuon2Tau1NotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  } else {
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon2Tau1ProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxMuon2Tau1MassReco;
   	            if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon2Tau1ReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hMuon2Tau1NotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  }
		}
              }
	    } else {
              _hMuon2Tau1OSLS[NpdfID]->Fill(patMuon->charge() * patTau->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              if((patTau->charge() * patMuon->charge())<0) {
                _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon2Tau1ProductsAndMetMassReco;
                _UseCollinerApproxMassReco = _UseCollinerApproxMuon2Tau1MassReco;
                if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon2Tau1ReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	        else {_hMuon2Tau1NotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
              } else {
                _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon2Tau1ProductsAndMetMassReco;
                _UseCollinerApproxMassReco = _UseCollinerApproxMuon2Tau1MassReco;
   	        if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon2Tau1ReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	        else {_hMuon2Tau1NotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
              }
            }
	    _hMuon2Tau1PZeta[NpdfID]->Fill(CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau1PZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau1Zeta2D[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau1Zeta1D[NpdfID]->Fill((_Muon2Tau1PZetaCutCoefficient * CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1)) + 
				     (_Muon2Tau1PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
      }
      theNumberOfMuons = 0;
      for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin();patMuon != _patMuons->end(); ++patMuon ) {
	theNumberOfMuons++;
	int theNumberOfTaus = 0;
	for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
	  theNumberOfTaus++;
	  if ((passRecoTau2Cuts((*patTau),theNumberOfTaus-1)) &&
	      (passRecoMuon2Cuts((*patMuon),theNumberOfMuons-1)) &&
	      (passMuon2Tau2TopologyCuts((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1))) {
	    _hMuon2PtVsTau2Pt[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt(),smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau2DeltaR[NpdfID]->Fill(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), smearedMuonMomentumVector.at(theNumberOfMuons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau2DeltaPtDivSumPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt()) / (smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() + smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau2DeltaPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau2CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau2_Muon2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2MetDeltaPhiVsMuon2Tau2CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - theMETVector.phi())), cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau2_Tau2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
            _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon2Tau2ProductsAndMetMassReco;
            _UseCollinerApproxMassReco = _UseCollinerApproxMuon2Tau2MassReco;
	    if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon2Tau2ReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    else {_hMuon2Tau2NotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    _hMuon2Tau2_Muon2MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau2_Tau2MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau),theNumberOfTaus-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    if(_UseTauSeedTrackForMuon2Tau2DiscrByOSLS) {
	      if (patTau->isCaloTau()) {
		if( (patTau->leadTrack().isNonnull()) ) {
		  _hMuon2Tau2OSLS[NpdfID]->Fill(patMuon->charge() * patTau->leadTrack()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		}
              } else {
		if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
		  _hMuon2Tau2OSLS[NpdfID]->Fill(patMuon->charge() * patTau->leadPFChargedHadrCand()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
                  if((patTau->leadPFChargedHadrCand()->charge() * patMuon->charge())<0) {
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon2Tau2ProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxMuon2Tau2MassReco;
   	            if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon2Tau2ReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hMuon2Tau2NotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  } else {
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon2Tau2ProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxMuon2Tau2MassReco;
   	            if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon2Tau2ReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hMuon2Tau2NotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  }
		}
              }
	    } else {
              _hMuon2Tau2OSLS[NpdfID]->Fill(patMuon->charge() * patTau->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              if((patTau->charge() * patMuon->charge())<0) {
                _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon2Tau2ProductsAndMetMassReco;
                _UseCollinerApproxMassReco = _UseCollinerApproxMuon2Tau2MassReco;
                if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon2Tau2ReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	        else {_hMuon2Tau2NotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
              } else {
                _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfMuon2Tau2ProductsAndMetMassReco;
                _UseCollinerApproxMassReco = _UseCollinerApproxMuon2Tau2MassReco;
   	        if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).first) {_hMuon2Tau2ReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	        else {_hMuon2Tau2NotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
              }
            }
	    _hMuon2Tau2PZeta[NpdfID]->Fill(CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau2PZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau2Zeta2D[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon2Tau2Zeta1D[NpdfID]->Fill((_Muon2Tau2PZetaCutCoefficient * CalculatePZeta((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1)) + 
				     (_Muon2Tau2PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patMuon),theNumberOfMuons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
      }
      int theNumberOfElectrons = 0;
      for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
	theNumberOfElectrons++;
	int theNumberOfTaus = 0;
	for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
	  theNumberOfTaus++;
	  if ((passRecoTau1Cuts((*patTau),theNumberOfTaus-1)) &&
	      (passRecoElectron1Cuts((*patElectron),theNumberOfElectrons-1)) &&
	      (passElectron1Tau1TopologyCuts((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1))) {
            _hElectron1Tau1_Electron1IsZee[NpdfID]->Fill(isZee(smearedElectronMomentumVector.at(theNumberOfElectrons-1)).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
            _hElectron1PtVsTau1Pt[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt(),smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau1DeltaR[NpdfID]->Fill(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), smearedElectronMomentumVector.at(theNumberOfElectrons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau1DeltaPtDivSumPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt()) / (smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() + smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau1DeltaPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau1CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau1_Electron1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1MetDeltaPhiVsElectron1Tau1CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - theMETVector.phi())), cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau1_Tau1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
            _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron1Tau1ProductsAndMetMassReco;
            _UseCollinerApproxMassReco = _UseCollinerApproxElectron1Tau1MassReco;
	    if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron1Tau1ReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    else {_hElectron1Tau1NotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    _hElectron1Tau1_Electron1MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau1_Tau1MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau),theNumberOfTaus-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    if(_UseTauSeedTrackForElectron1Tau1DiscrByOSLS) {
	      if (patTau->isCaloTau()) {
		if( (patTau->leadTrack().isNonnull()) ) {
		  _hElectron1Tau1OSLS[NpdfID]->Fill(patElectron->charge() * patTau->leadTrack()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		}
	      } else {
		if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
		  _hElectron1Tau1OSLS[NpdfID]->Fill(patElectron->charge() * patTau->leadPFChargedHadrCand()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
                  if((patTau->leadPFChargedHadrCand()->charge() * patElectron->charge())<0) {
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron1Tau1ProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxElectron1Tau1MassReco;
	            if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron1Tau1ReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hElectron1Tau1NotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  } else {
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron1Tau1ProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxElectron1Tau1MassReco;
	            if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron1Tau1ReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hElectron1Tau1NotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  }
		}
              }
	    } else {
              _hElectron1Tau1OSLS[NpdfID]->Fill(patElectron->charge() * patTau->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              if((patTau->charge() * patElectron->charge())<0) {
                _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron1Tau1ProductsAndMetMassReco;
                _UseCollinerApproxMassReco = _UseCollinerApproxElectron1Tau1MassReco;
                if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron1Tau1ReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	        else {_hElectron1Tau1NotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
              } else {
                _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron1Tau1ProductsAndMetMassReco;
                _UseCollinerApproxMassReco = _UseCollinerApproxElectron1Tau1MassReco;
	        if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron1Tau1ReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	        else {_hElectron1Tau1NotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
              }
            }
	    _hElectron1Tau1PZeta[NpdfID]->Fill(CalculatePZeta((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau1PZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau1Zeta2D[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),CalculatePZeta((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau1Zeta1D[NpdfID]->Fill((_Electron1Tau1PZetaCutCoefficient * CalculatePZeta((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1)) + (_Electron1Tau1PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
      }
      theNumberOfElectrons = 0;
      for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
	theNumberOfElectrons++;
	int theNumberOfTaus = 0;
	for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
	  theNumberOfTaus++;
	  if ((passRecoTau2Cuts((*patTau),theNumberOfTaus-1)) &&
	      (passRecoElectron1Cuts((*patElectron),theNumberOfElectrons-1)) &&
	      (passElectron1Tau2TopologyCuts((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1))) {
            _hElectron1Tau2_Electron1IsZee[NpdfID]->Fill(isZee(smearedElectronMomentumVector.at(theNumberOfElectrons-1)).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
            _hElectron1PtVsTau2Pt[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt(),smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau2DeltaR[NpdfID]->Fill(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), smearedElectronMomentumVector.at(theNumberOfElectrons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau2DeltaPtDivSumPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt()) / (smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() + smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau2DeltaPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau2CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau2_Electron1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1MetDeltaPhiVsElectron1Tau2CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - theMETVector.phi())), cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau2_Tau2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
            _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron1Tau2ProductsAndMetMassReco;
            _UseCollinerApproxMassReco = _UseCollinerApproxElectron1Tau2MassReco;
	    if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron1Tau2ReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    else {_hElectron1Tau2NotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    _hElectron1Tau2_Electron1MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau2_Tau2MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau),theNumberOfTaus-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    if(_UseTauSeedTrackForElectron1Tau2DiscrByOSLS) {
	      if (patTau->isCaloTau()) {
		if( (patTau->leadTrack().isNonnull()) ) {
		  _hElectron1Tau2OSLS[NpdfID]->Fill(patElectron->charge() * patTau->leadTrack()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		}
	      } else {
		if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
		  _hElectron1Tau2OSLS[NpdfID]->Fill(patElectron->charge() * patTau->leadPFChargedHadrCand()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
                  if((patTau->leadPFChargedHadrCand()->charge() * patElectron->charge())<0) {
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron1Tau2ProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxElectron1Tau2MassReco;
	            if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron1Tau2ReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hElectron1Tau2NotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  } else {
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron1Tau2ProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxElectron1Tau2MassReco;
	            if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron1Tau2ReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hElectron1Tau2NotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  }
		}
              }
	    } else {
              _hElectron1Tau2OSLS[NpdfID]->Fill(patElectron->charge() * patTau->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              if((patTau->charge() * patElectron->charge())<0) {
                _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron1Tau2ProductsAndMetMassReco;
                _UseCollinerApproxMassReco = _UseCollinerApproxElectron1Tau2MassReco;
                if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron1Tau2ReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	        else {_hElectron1Tau2NotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
              } else {
                _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron1Tau2ProductsAndMetMassReco;
                _UseCollinerApproxMassReco = _UseCollinerApproxElectron1Tau2MassReco;
	        if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron1Tau2ReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	        else {_hElectron1Tau2NotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
              }
            }
	    _hElectron1Tau2PZeta[NpdfID]->Fill(CalculatePZeta((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau2PZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau2Zeta2D[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),CalculatePZeta((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Tau2Zeta1D[NpdfID]->Fill((_Electron1Tau2PZetaCutCoefficient * CalculatePZeta((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1)) + (_Electron1Tau2PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
      }
      theNumberOfElectrons = 0;
      for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
	theNumberOfElectrons++;
	int theNumberOfTaus = 0;
	for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
	  theNumberOfTaus++;
	  if ((passRecoTau1Cuts((*patTau),theNumberOfTaus-1)) &&
	      (passRecoElectron2Cuts((*patElectron),theNumberOfElectrons-1)) &&
	      (passElectron2Tau1TopologyCuts((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1))) {
            _hElectron2Tau1_Electron2IsZee[NpdfID]->Fill(isZee(smearedElectronMomentumVector.at(theNumberOfElectrons-1)).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
            _hElectron2PtVsTau1Pt[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt(),smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau1DeltaR[NpdfID]->Fill(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), smearedElectronMomentumVector.at(theNumberOfElectrons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau1DeltaPtDivSumPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt()) / (smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() + smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau1DeltaPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau1CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau1_Electron2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2MetDeltaPhiVsElectron2Tau1CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - theMETVector.phi())), cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau1_Tau1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
            _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron2Tau1ProductsAndMetMassReco;
            _UseCollinerApproxMassReco = _UseCollinerApproxElectron2Tau1MassReco;
	    if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron2Tau1ReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    else {_hElectron2Tau1NotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    _hElectron2Tau1_Electron2MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau1_Tau1MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau),theNumberOfTaus-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    if(_UseTauSeedTrackForElectron2Tau1DiscrByOSLS) {
	      if (patTau->isCaloTau()) {
		if( (patTau->leadTrack().isNonnull()) ) {
		  _hElectron2Tau1OSLS[NpdfID]->Fill(patElectron->charge() * patTau->leadTrack()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		}
	      } else {
		if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
		  _hElectron2Tau1OSLS[NpdfID]->Fill(patElectron->charge() * patTau->leadPFChargedHadrCand()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
                  if((patTau->leadPFChargedHadrCand()->charge() * patElectron->charge())<0) {
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron2Tau1ProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxElectron2Tau1MassReco;
	            if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron2Tau1ReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hElectron2Tau1NotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  } else {
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron2Tau1ProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxElectron2Tau1MassReco;
	            if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron2Tau1ReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hElectron2Tau1NotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  }
		}
              }
	    } else {
              _hElectron2Tau1OSLS[NpdfID]->Fill(patElectron->charge() * patTau->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              if((patTau->charge() * patElectron->charge())<0) {
                _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron2Tau1ProductsAndMetMassReco;
                _UseCollinerApproxMassReco = _UseCollinerApproxElectron2Tau1MassReco;
                if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron2Tau1ReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	        else {_hElectron2Tau1NotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
              } else {
                _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron2Tau1ProductsAndMetMassReco;
                _UseCollinerApproxMassReco = _UseCollinerApproxElectron2Tau1MassReco;
	        if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron2Tau1ReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	        else {_hElectron2Tau1NotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
              }
            }
	    _hElectron2Tau1PZeta[NpdfID]->Fill(CalculatePZeta((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau1PZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau1Zeta2D[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),CalculatePZeta((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau1Zeta1D[NpdfID]->Fill((_Electron2Tau1PZetaCutCoefficient * CalculatePZeta((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1)) + (_Electron2Tau1PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
      }
      theNumberOfElectrons = 0;
      for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin();patElectron != _patElectrons->end(); ++patElectron ) {
	theNumberOfElectrons++;
	int theNumberOfTaus = 0;
	for ( pat::TauCollection::const_iterator patTau = _patTaus->begin();patTau != _patTaus->end(); ++patTau ) {
	  theNumberOfTaus++;
	  if ((passRecoTau2Cuts((*patTau),theNumberOfTaus-1)) &&
	      (passRecoElectron2Cuts((*patElectron),theNumberOfElectrons-1)) &&
	      (passElectron2Tau2TopologyCuts((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1))) {
            _hElectron2Tau2_Electron2IsZee[NpdfID]->Fill(isZee(smearedElectronMomentumVector.at(theNumberOfElectrons-1)).first,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
            _hElectron2PtVsTau2Pt[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt(),smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau2DeltaR[NpdfID]->Fill(reco::deltaR(smearedTauMomentumVector.at(theNumberOfTaus-1), smearedElectronMomentumVector.at(theNumberOfElectrons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau2DeltaPtDivSumPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt()) / (smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() + smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau2DeltaPt[NpdfID]->Fill((smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).pt() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau2CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau2_Electron2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2MetDeltaPhiVsElectron2Tau2CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - theMETVector.phi())), cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons-1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau2_Tau2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus-1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
            _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron2Tau2ProductsAndMetMassReco;
            _UseCollinerApproxMassReco = _UseCollinerApproxElectron2Tau2MassReco;
	    if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron2Tau2ReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    else {_hElectron2Tau2NotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    _hElectron2Tau2_Electron2MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau2_Tau2MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau),theNumberOfTaus-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    if(_UseTauSeedTrackForElectron2Tau2DiscrByOSLS) {
	      if (patTau->isCaloTau()) {
		if( (patTau->leadTrack().isNonnull()) ) {
		  _hElectron2Tau2OSLS[NpdfID]->Fill(patElectron->charge() * patTau->leadTrack()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
		}
	      } else {
		if( (patTau->leadPFChargedHadrCand().isNonnull()) ) {
		  _hElectron2Tau2OSLS[NpdfID]->Fill(patElectron->charge() * patTau->leadPFChargedHadrCand()->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
                  if((patTau->leadPFChargedHadrCand()->charge() * patElectron->charge())<0) {
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron2Tau2ProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxElectron2Tau2MassReco;
	            if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron2Tau2ReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hElectron2Tau2NotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  } else {
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron2Tau2ProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxElectron2Tau2MassReco;
	            if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron2Tau2ReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hElectron2Tau2NotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  }
		}
              }
	    } else {
              _hElectron2Tau2OSLS[NpdfID]->Fill(patElectron->charge() * patTau->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              if((patTau->charge() * patElectron->charge())<0) {
                _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron2Tau2ProductsAndMetMassReco;
                _UseCollinerApproxMassReco = _UseCollinerApproxElectron2Tau2MassReco;
                if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron2Tau2ReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	        else {_hElectron2Tau2NotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
              } else {
                _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfElectron2Tau2ProductsAndMetMassReco;
                _UseCollinerApproxMassReco = _UseCollinerApproxElectron2Tau2MassReco;
	        if(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).first) {_hElectron2Tau2ReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	        else {_hElectron2Tau2NotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
              }
            }
	    _hElectron2Tau2PZeta[NpdfID]->Fill(CalculatePZeta((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau2PZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau2Zeta2D[NpdfID]->Fill(CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),CalculatePZeta((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron2Tau2Zeta1D[NpdfID]->Fill((_Electron2Tau2PZetaCutCoefficient * CalculatePZeta((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1)) + (_Electron2Tau2PZetaVisCutCoefficient * CalculatePZetaVis((*patTau),theNumberOfTaus-1,(*patElectron),theNumberOfElectrons-1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
      }
      int theNumberOfMuons1 = 0;
      for ( pat::MuonCollection::const_iterator patMuon1 = _patMuons->begin();patMuon1 != _patMuons->end(); ++patMuon1 ) {
	theNumberOfMuons1++;
	int theNumberOfMuons2 = 0;
	for ( pat::MuonCollection::const_iterator patMuon2 = _patMuons->begin();patMuon2 != _patMuons->end(); ++patMuon2 ) {
	  theNumberOfMuons2++;
	  if ((passRecoMuon1Cuts((*patMuon1),theNumberOfMuons1 - 1)) && 
	      (passRecoMuon2Cuts((*patMuon2),theNumberOfMuons2 - 1)) && 
	      (passDiMuonTopologyCuts((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1)) &&
              (theNumberOfMuons2 > theNumberOfMuons1)) {
            _hMuon1PtVsMuon2Pt[NpdfID]->Fill(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons1 - 1).pt(),smearedMuonPtEtaPhiMVector.at(theNumberOfMuons2 - 1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Muon2DeltaR[NpdfID]->Fill(reco::deltaR(smearedMuonMomentumVector.at(theNumberOfMuons1 - 1), smearedMuonMomentumVector.at(theNumberOfMuons2 - 1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Muon2DeltaPtDivSumPt[NpdfID]->Fill((smearedMuonPtEtaPhiMVector.at(theNumberOfMuons1 - 1).pt() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons2 - 1).pt()) / (smearedMuonPtEtaPhiMVector.at(theNumberOfMuons1 - 1).pt() + smearedMuonPtEtaPhiMVector.at(theNumberOfMuons2 - 1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Muon2DeltaPt[NpdfID]->Fill((smearedMuonPtEtaPhiMVector.at(theNumberOfMuons1 - 1).pt() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons2 - 1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Muon2CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons1 - 1).phi() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiMuon_Muon1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons1 - 1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiMuon_Muon2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons2 - 1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1MetDeltaPhiVsMuon1Muon2CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons1 - 1).phi() - theMETVector.phi())), cos(TMath::Abs(normalizedPhi(smearedMuonPtEtaPhiMVector.at(theNumberOfMuons1 - 1).phi() - smearedMuonPtEtaPhiMVector.at(theNumberOfMuons2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
            _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfDiMuonProductsAndMetMassReco;
            _UseCollinerApproxMassReco = _UseCollinerApproxDiMuonMassReco;
	    if(CalculateThe4Momentum((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1).first) {_hDiMuonReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    else {_hDiMuonNotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    _hDiMuon_Muon1MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patMuon1),theNumberOfMuons1 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiMuon_Muon2MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patMuon2),theNumberOfMuons2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hMuon1Muon2OSLS[NpdfID]->Fill(patMuon1->charge() * patMuon2->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiMuonPZeta[NpdfID]->Fill(CalculatePZeta((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiMuonPZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiMuonZeta2D[NpdfID]->Fill(CalculatePZetaVis((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1),CalculatePZeta((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiMuonZeta1D[NpdfID]->Fill((_DiMuonPZetaCutCoefficient * CalculatePZeta((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1)) + (_DiMuonPZetaVisCutCoefficient * CalculatePZetaVis((*patMuon1),theNumberOfMuons1 - 1,(*patMuon2),theNumberOfMuons2 - 1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
      }
      int theNumberOfElectrons1 = 0;
      for ( pat::ElectronCollection::const_iterator patElectron1 = _patElectrons->begin();patElectron1 != _patElectrons->end(); ++patElectron1 ) {
	theNumberOfElectrons1++;
	int theNumberOfElectrons2 = 0;
	for ( pat::ElectronCollection::const_iterator patElectron2 = _patElectrons->begin();patElectron2 != _patElectrons->end(); ++patElectron2 ) {
	  theNumberOfElectrons2++;
	  if ((passRecoElectron1Cuts((*patElectron1),theNumberOfElectrons1 - 1)) && 
	      (passRecoElectron2Cuts((*patElectron2),theNumberOfElectrons2 - 1)) && 
	      (passDiElectronTopologyCuts((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1)) &&
              (theNumberOfElectrons2 > theNumberOfElectrons1)) {
	    _hElectron1PtVsElectron2Pt[NpdfID]->Fill(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons1 - 1).pt(),smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons2 - 1).pt(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Electron2DeltaR[NpdfID]->Fill(reco::deltaR(smearedElectronMomentumVector.at(theNumberOfElectrons1 - 1), smearedElectronMomentumVector.at(theNumberOfElectrons2 - 1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Electron2DeltaPtDivSumPt[NpdfID]->Fill((smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons1 - 1).pt() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons2 - 1).pt()) / (smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons1 - 1).pt() + smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons2 - 1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Electron2DeltaPt[NpdfID]->Fill((smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons1 - 1).pt() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons2 - 1).pt()),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Electron2CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons1 - 1).phi() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiElectron_Electron1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons1 - 1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiElectron_Electron2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons2 - 1).phi() - theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1MetDeltaPhiVsElectron1Electron2CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons1 - 1).phi() - theMETVector.phi())), cos(TMath::Abs(normalizedPhi(smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons1 - 1).phi() - smearedElectronPtEtaPhiMVector.at(theNumberOfElectrons2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
            _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfDiElectronProductsAndMetMassReco;
            _UseCollinerApproxMassReco = _UseCollinerApproxDiElectronMassReco;
	    if(CalculateThe4Momentum((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1).first) {_hDiElectronReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    else {_hDiElectronNotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    _hDiElectron_Electron1MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patElectron1),theNumberOfElectrons1 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiElectron_Electron2MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patElectron2),theNumberOfElectrons2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hElectron1Electron2OSLS[NpdfID]->Fill(patElectron1->charge() * patElectron2->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiElectronPZeta[NpdfID]->Fill(CalculatePZeta((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiElectronPZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiElectronZeta2D[NpdfID]->Fill(CalculatePZetaVis((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1),CalculatePZeta((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiElectronZeta1D[NpdfID]->Fill((_DiElectronPZetaCutCoefficient * CalculatePZeta((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1)) + (_DiElectronPZetaVisCutCoefficient * CalculatePZetaVis((*patElectron1),theNumberOfElectrons1 - 1,(*patElectron2),theNumberOfElectrons2 - 1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	  }
	}
      }
      int theNumberOfTaus1 = 0;
      for ( pat::TauCollection::const_iterator patTau1 = _patTaus->begin();patTau1 != _patTaus->end(); ++patTau1 ) {
	theNumberOfTaus1++;
	int theNumberOfTaus2 = 0;
	for ( pat::TauCollection::const_iterator patTau2 = _patTaus->begin();patTau2 != _patTaus->end(); ++patTau2 ) {
	  theNumberOfTaus2++;
	  if ((passRecoTau1Cuts((*patTau1),theNumberOfTaus1 - 1)) &&
	      (passRecoTau2Cuts((*patTau2),theNumberOfTaus2 - 1)) &&
	      (passDiTauTopologyCuts((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1)) &&
              (theNumberOfTaus2 > theNumberOfTaus1)) {
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
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfDiTauProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxDiTauMassReco;
	            if(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).first) {_hDiTauReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hDiTauNotReconstructableMassOS[NpdfID]->Fill(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  } else {
                    _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfDiTauProductsAndMetMassReco;
                    _UseCollinerApproxMassReco = _UseCollinerApproxDiTauMassReco;
	            if(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).first) {_hDiTauReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	            else {_hDiTauNotReconstructableMassLS[NpdfID]->Fill(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
                  }
		}
              }
	    } else {_hTau1Tau2OSLS[NpdfID]->Fill(patTau1->charge() * patTau2->charge(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    _hTau1Tau2CosDphi[NpdfID]->Fill(cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus1 - 1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiTau_Tau1MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus1 - 1).phi() -  theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiTau_Tau2MetDeltaPhi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus2 - 1).phi() -  theMETVector.phi())),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hTau1MetDeltaPhiVsTau1Tau2CosDphi[NpdfID]->Fill(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus1 - 1).phi() -  theMETVector.phi())), cos(TMath::Abs(normalizedPhi(smearedTauPtEtaPhiMVector.at(theNumberOfTaus1 - 1).phi() - smearedTauPtEtaPhiMVector.at(theNumberOfTaus2 - 1).phi()))),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
            _UseVectorSumOfVisProductsAndMetMassReco = _UseVectorSumOfDiTauProductsAndMetMassReco;
            _UseCollinerApproxMassReco = _UseCollinerApproxDiTauMassReco;
	    if(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).first) {_hDiTauReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    else {_hDiTauNotReconstructableMass[NpdfID]->Fill(CalculateThe4Momentum((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1).second.M(),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));}
	    _hDiTau_Tau1MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau1),theNumberOfTaus1 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiTau_Tau2MetMt[NpdfID]->Fill(CalculateLeptonMetMt((*patTau2),theNumberOfTaus2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiTauPZeta[NpdfID]->Fill(CalculatePZeta((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiTauPZetaVis[NpdfID]->Fill(CalculatePZetaVis((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiTauZeta2D[NpdfID]->Fill(CalculatePZetaVis((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1),CalculatePZeta((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
	    _hDiTauZeta1D[NpdfID]->Fill((_DiTauPZetaCutCoefficient * CalculatePZeta((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1)) + (_DiTauPZetaVisCutCoefficient * CalculatePZetaVis((*patTau1),theNumberOfTaus1 - 1,(*patTau2),theNumberOfTaus2 - 1)),isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
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
              int theNumberOfJets = 0;
              double sumpxForMht = 0.0;
              double sumpyForMht = 0.0;
              double phiForMht = 0.0;
              for ( pat::JetCollection::const_iterator patJet = _patJets->begin();
                patJet != _patJets->end(); ++patJet ) {
                theNumberOfJets++;
                if (!passRecoJetCuts((*patJet),theNumberOfJets-1)) continue;
                sumpxForMht = sumpxForMht - smearedJetMomentumVector.at(theNumberOfJets-1).px();
                sumpyForMht = sumpyForMht - smearedJetMomentumVector.at(theNumberOfJets-1).py();
              }
              if (sumpxForMht >= 0)
                phiForMht = atan(sumpyForMht/sumpxForMht);
              if (sumpxForMht < 0 && sumpyForMht >= 0)
                phiForMht = atan(sumpyForMht/sumpxForMht) + 3.1415;
              if (sumpxForMht < 0 && sumpyForMht < 0)
                phiForMht = atan(sumpyForMht/sumpxForMht) - 3.1415;
              double dphi1MHT = normalizedPhi(smearedJetPtEtaPhiMVector.at(numberJets1 - 1).phi() - phiForMht);
              double dphi2MHT = normalizedPhi(smearedJetPtEtaPhiMVector.at(numberJets2 - 1).phi() - phiForMht);
              double r1;
              double r2;
              double alpha;
              dphi1 = normalizedPhi(smearedJetPtEtaPhiMVector.at(numberJets1 - 1).phi() - theMETVector.phi());
              dphi2 = normalizedPhi(smearedJetPtEtaPhiMVector.at(numberJets2 - 1).phi() - theMETVector.phi());
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
              _hDphi1MHT[NpdfID]->Fill(dphi1MHT,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              _hDphi2MHT[NpdfID]->Fill(dphi2MHT,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              _hDphi1[NpdfID]->Fill(dphi1,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              _hDphi2[NpdfID]->Fill(dphi2,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
              _hDphi1VsDphi2[NpdfID]->Fill(dphi1,dphi2,isrgluon_weight * isrgamma_weight * fsr_weight * pdfWeightVector.at(NpdfID));
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
  iEvent.getByLabel("particleFlow", _pflow);
  iEvent.getByLabel("hpsPFTauProducer", _hpsTau);
  if(!isData) {iEvent.getByLabel("genMetTrue", genTrue);}
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
    if(reco::deltaR(theObject, smearedElectronMomentumVector.at(i - 1)) < _DiElectronDeltaRCut) continue;
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

pair<bool, pair<float, float> > HiMassTauAnalysis::isZmm(reco::Candidate::LorentzVector theObject) {
  pair<bool, pair<float, float> > theOutPutPair;
  bool eventIsZmm = false;
  bool massWindow = false;
  bool ptAsymmWindow = false;
  const float zmmMass = 90.1876;
  const float zmmWidht = 2.4952;
  float zmmPtAsymmetry = -10.;
  unsigned int i = 0;
  pair<float, float> theMassPtAsymmPair;
  // if mass is within 3 sigmas of z or pt asymmetry is small set to true.				     
  for(pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon){
    i++;
    if(reco::deltaR(theObject, smearedMuonMomentumVector.at(i - 1)) < _DiMuonDeltaRCut) continue;
    if(theObject == smearedMuonMomentumVector.at(i - 1))continue;
    reco::Candidate::LorentzVector The_LorentzVect = theObject + smearedMuonMomentumVector.at(i - 1);
    zmmPtAsymmetry = (theObject.pt() - smearedMuonPtEtaPhiMVector.at(i - 1).pt())/(theObject.pt() + smearedMuonPtEtaPhiMVector.at(i - 1).pt());
    theMassPtAsymmPair = make_pair<float, float>(The_LorentzVect.M(), zmmPtAsymmetry);
    
    if(The_LorentzVect.M() < (zmmMass + 3.*zmmWidht) && The_LorentzVect.M() > (zmmMass - 3.*zmmWidht))massWindow = true;
    if(fabs(zmmPtAsymmetry) < 0.20) ptAsymmWindow = true;
    if(massWindow || ptAsymmWindow){
      eventIsZmm = true;
      break;
    }
  }
  theOutPutPair = make_pair<bool, pair<float, float> >(eventIsZmm, theMassPtAsymmPair);
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
//    if ( it->ref().isNonnull() && it->ref().isValid() ) {
    if ( it->isNonnull() && it->isAvailable() ) {
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
//    if ( it->ref().isNonnull() && it->ref().isValid() ) {
    if ( it->isNonnull() && it->isAvailable() ) {
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
        if((**iTrk).trackRef().isNonnull()) {
          if((reco::deltaR((**iTrk).eta(),(**iTrk).phi(),patTau.leadPFChargedHadrCand()->eta(),patTau.leadPFChargedHadrCand()->phi())<_RecoTauIsoDeltaRCone) && ((**iTrk).pt()>_RecoTauTrackIsoTrkThreshold)
             && ((**iTrk).trackRef()->numberOfValidHits() > 8) && ((**iTrk).trackRef().isNonnull()) && (fabs((**iTrk).trackRef()->dxy(thePrimaryEventVertex.position())) < 0.03) 
             && (fabs((**iTrk).trackRef()->dz(thePrimaryEventVertex.position())) < 0.2) ) {
            nIsoTrks++;
            sumPtIsoTrks = sumPtIsoTrks + (**iTrk).pt();
          }
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
           && ((**iTrk).trackRef()->numberOfValidHits() > 8) ) {
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

//-----Smear the light leptons (for studies of systematic uncertanties)
pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> HiMassTauAnalysis::SmearLightLepton(const pat::Muon& patMuon) {
  double smearedPt;
  double smearedEta;
  double smearedPhi;
  if( (_GenParticleSource.label() != "") && (!isData) ) {
    if(matchToGen(patMuon).first) {
      reco::Candidate::LorentzVector unsmearedMomentum = matchToGen(patMuon).second;
      if(_SmearTheMuon) {
        smearedPt = (unsmearedMomentum.pt() * _MuonPtScaleOffset) + ((patMuon.pt() -  unsmearedMomentum.pt()) * _MuonPtSigmaOffset);
        smearedEta = (unsmearedMomentum.eta() * _MuonEtaScaleOffset) + ((patMuon.eta() -  unsmearedMomentum.eta()) * _MuonEtaSigmaOffset);
        smearedPhi = (unsmearedMomentum.phi() * _MuonPhiScaleOffset) + ((patMuon.phi() -  unsmearedMomentum.phi()) * _MuonPhiSigmaOffset);
      } else {
        smearedPt  = patMuon.pt();
        smearedEta  = patMuon.eta();
        smearedPhi  = patMuon.phi();
      }
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
  if( (_GenParticleSource.label() != "") && (!isData) ) {
    if(matchToGen(patElectron).first) {
      reco::Candidate::LorentzVector unsmearedMomentum = matchToGen(patElectron).second; 
      if(_UseHeepInfo) {
        heep::Ele theHeepElec(patElectron);
        if(_SmearTheElectron) {
	  smearedEt = (unsmearedMomentum.energy() * sin(unsmearedMomentum.theta()) * _ElectronPtScaleOffset) + ((theHeepElec.et() -  (unsmearedMomentum.energy() * sin(unsmearedMomentum.theta()))) * _ElectronPtSigmaOffset);
	  smearedPt = (unsmearedMomentum.pt() * _ElectronPtScaleOffset) + ((patElectron.pt() -  unsmearedMomentum.pt()) * _ElectronPtSigmaOffset);
	  smearedScEta = (unsmearedMomentum.eta() * _ElectronEtaScaleOffset) + ((theHeepElec.scEta() -  unsmearedMomentum.eta()) * _ElectronEtaSigmaOffset);
	  smearedEta = (unsmearedMomentum.eta() * _ElectronEtaScaleOffset) + ((patElectron.eta() -  unsmearedMomentum.eta()) * _ElectronEtaSigmaOffset);
	  smearedPhi = (unsmearedMomentum.phi() * _ElectronPhiScaleOffset) + ((patElectron.phi() -  unsmearedMomentum.phi()) * _ElectronPhiSigmaOffset);
        } else {
          smearedEt  = theHeepElec.et();
          smearedPt  = patElectron.pt();
          smearedScEta  = theHeepElec.scEta();
          smearedEta  = patElectron.eta();
          smearedPhi  = patElectron.phi();
        }
	math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector(smearedPt, smearedEta, smearedPhi, unsmearedMomentum.mass());
	math::PtEtaPhiMLorentzVector smearedEtEtaPhiMVector(smearedEt, smearedScEta, smearedPhi, unsmearedMomentum.mass());
	reco::Candidate::LorentzVector smearedMomentum(smearedPtEtaPhiMVector.px(), smearedPtEtaPhiMVector.py(), smearedPtEtaPhiMVector.pz(), patElectron.energy());
	pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> theSmearedMomentumPair;
	theSmearedMomentumPair = make_pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector>(smearedMomentum, smearedEtEtaPhiMVector);
	return theSmearedMomentumPair;
      } else {
        heep::Ele theHeepElec(patElectron);
        if(_SmearTheElectron) {
	  smearedPt = (unsmearedMomentum.pt() * _ElectronPtScaleOffset) + ((patElectron.pt() -  unsmearedMomentum.pt()) * _ElectronPtSigmaOffset);
	  smearedEta = (unsmearedMomentum.eta() * _ElectronEtaScaleOffset) + ((patElectron.eta() -  unsmearedMomentum.eta()) * _ElectronEtaSigmaOffset);
	  smearedPhi = (unsmearedMomentum.phi() * _ElectronPhiScaleOffset) + ((patElectron.phi() -  unsmearedMomentum.phi()) * _ElectronPhiSigmaOffset);
        } else {
          smearedEt  = theHeepElec.et();
          smearedPt  = patElectron.pt();
          smearedScEta  = theHeepElec.scEta();
          smearedEta  = patElectron.eta();
          smearedPhi  = patElectron.phi();
        }
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
  if( (_GenParticleSource.label() != "") && (!isData) ) {
    if(matchToGen(patTau).first) {
      reco::Candidate::LorentzVector unsmearedMomentum = matchToGen(patTau).second; 
      if(_SmearTheTau) {
        smearedPt = (unsmearedMomentum.pt() * _TauPtScaleOffset) + ((patTau.pt() -  unsmearedMomentum.pt()) * _TauPtSigmaOffset);
        smearedEta = (unsmearedMomentum.eta() * _TauEtaScaleOffset) + ((patTau.eta() -  unsmearedMomentum.eta()) * _TauEtaSigmaOffset);
        smearedPhi = (unsmearedMomentum.phi() * _TauPhiScaleOffset) + ((patTau.phi() -  unsmearedMomentum.phi()) * _TauPhiSigmaOffset);
      } else {
        smearedPt  = patTau.pt();
        smearedEta  = patTau.eta();
        smearedPhi  = patTau.phi();
      }
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
pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> HiMassTauAnalysis::SmearJet(const pat::Jet& patJet,JetCorrectionUncertainty* jesUncertainty) {
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
  if( (_GenParticleSource.label() != "") && (!isData) ) {
    for ( pat::MuonCollection::const_iterator patMuon = _patMuons->begin(); patMuon != _patMuons->end(); ++patMuon ) {
      if( ((reco::deltaR(patMuon->p4(), tempJetVector) < _JetMuon1MatchingDeltaR) || 
           (reco::deltaR(patMuon->p4(), tempJetVector) < _JetMuon2MatchingDeltaR)) 
          && (matchToGen(*patMuon).first) ) {isRealJet = false;}
    }
    for ( pat::ElectronCollection::const_iterator patElectron = _patElectrons->begin(); patElectron != _patElectrons->end(); ++patElectron ) {
      if( ((reco::deltaR(patElectron->p4(), tempJetVector) < _JetElectron1MatchingDeltaR) ||
           (reco::deltaR(patElectron->p4(), tempJetVector) < _JetElectron2MatchingDeltaR))
          && (matchToGen(*patElectron).first) ) {isRealJet = false;}
    }
    for ( pat::TauCollection::const_iterator patTau = _patTaus->begin(); patTau != _patTaus->end(); ++patTau ) {
      if( ((reco::deltaR(patTau->p4(), tempJetVector) < _JetTau1MatchingDeltaR) ||
           (reco::deltaR(patTau->p4(), tempJetVector) < _JetTau2MatchingDeltaR))
          && (matchToGen(*patTau).first) ) {isRealJet = false;}
    }
    if(isRealJet) {
      if(_UseDataBaseForJEC) {
        jesUncertainty->setJetEta(tempJetVector.eta());
        jesUncertainty->setJetPt(tempJetVector.pt());
        double xxx = jesUncertainty->getUncertainty(true);
        _JetEnergyScaleOffset = 1. + (JES_UpOrDown * xxx);
        reco::Candidate::LorentzVector smearedMomentum = _JetEnergyScaleOffset * tempJetVector;
        math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector(smearedMomentum.pt(), smearedMomentum.eta(), smearedMomentum.phi(), smearedMomentum.mass());
        pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector> theSmearedMomentumPair;
        theSmearedMomentumPair = make_pair<reco::Candidate::LorentzVector,math::PtEtaPhiMLorentzVector>(smearedMomentum, smearedPtEtaPhiMVector);
        return theSmearedMomentumPair;
      } else {
        reco::Candidate::LorentzVector smearedMomentum = _JetEnergyScaleOffset * tempJetVector;
        math::PtEtaPhiMLorentzVector smearedPtEtaPhiMVector(smearedMomentum.pt(), smearedMomentum.eta(), smearedMomentum.phi(), smearedMomentum.mass());
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
      _hNVertices[NpdfCounter]       = fs->make<TH1F>(("NVertices_"+j.str()).c_str(),       ("NVertices_"+j.str()).c_str(),       40, 0., 40.);
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
      _hNTau1[NpdfCounter]                             = fs->make<TH1F>(("NTau1_"+j.str()).c_str(),         ("NTau1_"+j.str()).c_str(),         21, 0., 20.);
      _hTauJet1Energy[NpdfCounter]                     = fs->make<TH1F>(("TauJet1Energy_"+j.str()).c_str(), ("TauJet1Energy_"+j.str()).c_str(), 200, 0., 500.);
      _hTauJet1Pt[NpdfCounter]                         = fs->make<TH1F>(("TauJet1Pt_"+j.str()).c_str(),     ("TauJet1Pt_"+j.str()).c_str(),     200, 0., 500.);
      _hTauJet1Eta[NpdfCounter]                        = fs->make<TH1F>(("TauJet1Eta_"+j.str()).c_str(),    ("TauJet1Eta_"+j.str()).c_str(),    72, -3.6, +3.6);
      _hTauJet1Phi[NpdfCounter]                        = fs->make<TH1F>(("TauJet1Phi_"+j.str()).c_str(),    ("TauJet1Phi_"+j.str()).c_str(),    36, -TMath::Pi(), +TMath::Pi());
      _hTauJet1NumSignalTracks[NpdfCounter]            = fs->make<TH1F>(("TauJet1NumSignalTracks_"+j.str()).c_str(),           ("TauJet1NumSignalTracks_"+j.str()).c_str(), 10, 0, 10);
      _hTauJet1NumSignalGammas[NpdfCounter]            = fs->make<TH1F>(("TauJet1NumSignalGammas_"+j.str()).c_str(),           ("TauJet1NumSignalGammas_"+j.str()).c_str(), 10, 0, 10);
      _hTauJet1SeedTrackPt[NpdfCounter]                = fs->make<TH1F>(("TauJet1SeedTrackPt_"+j.str()).c_str(),               ("TauJet1SeedTrackPt_"+j.str()).c_str(),     200, 0., 500.);
      _hTauJet1SeedTrackIpSignificance[NpdfCounter]    = fs->make<TH1F>(("TauJet1SeedTrackIpSignificance_"+j.str()).c_str(),   ("TauJet1SeedTrackIpSignificance_"+j.str()).c_str(), 100, 0., 100.);
      _hTauJet1SeedTrackNhits[NpdfCounter]             = fs->make<TH1F>(("TauJet1SeedTrackNhits_"+j.str()).c_str(),            ("TauJet1SeedTrackNhits_"+j.str()).c_str(), 40, 0., 40.);
      _hTauJet1SeedTrackChi2[NpdfCounter]              = fs->make<TH1F>(("TauJet1SeedTrackChi2_"+j.str()).c_str(),             ("TauJet1SeedTrackChi2_"+j.str()).c_str(), 50, 0., 100.);
      _hTauJet1Charge[NpdfCounter]                     = fs->make<TH1F>(("TauJet1Charge_"+j.str()).c_str(),                    ("TauJet1Charge_"+j.str()).c_str(), 5, 0., 5.);
      _hTauJet1SignalTracksMass[NpdfCounter]           = fs->make<TH1F>(("TauJet1SignalTracksMass_"+j.str()).c_str(),          ("TauJet1SignalTracksMass_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJet1SignalTracksAndGammasMass[NpdfCounter]  = fs->make<TH1F>(("TauJet1SignalTracksAndGammasMass_"+j.str()).c_str(), ("TauJet1SignalTracksAndGammasMass_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJet1SignalTracksChargeFraction[NpdfCounter] = fs->make<TH1F>(("TauJet1SignalTracksChargeFraction_"+j.str()).c_str(),("TauJet1SignalTracksChargeFraction_"+j.str()).c_str(), 30, 0., 1.5);
      _hTauJet1NumIsoTracks[NpdfCounter]               = fs->make<TH1F>(("TauJet1NumIsoTracks_"+j.str()).c_str(),              ("TauJet1NumIsoTracks_"+j.str()).c_str(), 10, 0, 10);
      _hTauJet1NumIsoGammas[NpdfCounter]               = fs->make<TH1F>(("TauJet1NumIsoGammas_"+j.str()).c_str(),              ("TauJet1NumIsoGammas_"+j.str()).c_str(), 10, 0, 10);
      _hTauJet1NumIsoCands[NpdfCounter]                = fs->make<TH1F>(("TauJet1NumIsoCands_"+j.str()).c_str(),               ("TauJet1NumIsoCands_"+j.str()).c_str(), 10, 0, 10);
      _hTauJet1SumPtIsoTracks[NpdfCounter]             = fs->make<TH1F>(("TauJet1SumPtIsoTracks_"+j.str()).c_str(),            ("TauJet1SumPtIsoTracks_"+j.str()).c_str(), 100, 0, 50);
      _hTauJet1SumPtIsoGammas[NpdfCounter]             = fs->make<TH1F>(("TauJet1SumPtIsoGammas_"+j.str()).c_str(),            ("TauJet1SumPtIsoGammas_"+j.str()).c_str(), 100, 0, 50);
      _hTauJet1SumPtIso[NpdfCounter]                   = fs->make<TH1F>(("TauJet1SumPtIso_"+j.str()).c_str(),                  ("TauJet1SumPtIso_"+j.str()).c_str(), 100, 0, 50);
      _hTauJet1NumberDensity[NpdfCounter]              = fs->make<TH1F>(("TauJet1NumberDensity_"+j.str()).c_str(),                  ("TauJet1NumberDensity_"+j.str()).c_str(), 500, 0, 10);
      _hTauJet1GenTauDeltaPhi[NpdfCounter]             = fs->make<TH1F>(("TauJet1GenTauDeltaPhi_"+j.str()).c_str(),            ("TauJet1GenTauDeltaPhi_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hTauJet1GenTauDeltaEta[NpdfCounter]             = fs->make<TH1F>(("TauJet1GenTauDeltaEta_"+j.str()).c_str(),            ("TauJet1GenTauDeltaEta_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hTauJet1GenTauDeltaPt[NpdfCounter]              = fs->make<TH1F>(("TauJet1GenTauDeltaPt_"+j.str()).c_str(),             ("TauJet1GenTauDeltaPt_"+j.str()).c_str(), 500, -5, 5);
      _hTauJet1SignalTracksMass1prong[NpdfCounter]     = fs->make<TH1F>(("TauJet1SignalTracksMass1prong_"+j.str()).c_str(),    ("TauJet1SignalTracksMass1prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJet1SignalTracksAndGammasMass1prong[NpdfCounter] = fs->make<TH1F>(("TauJet1SignalTracksAndGammasMass1prong_"+j.str()).c_str(), ("TauJet1SignalTracksAndGammasMass1prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJet1SignalTracksAndPiZerosMass1prong[NpdfCounter] = fs->make<TH1F>(("TauJet1SignalTracksAndPiZerosMass1prong_"+j.str()).c_str(), ("TauJet1SignalTracksAndPiZerosMass1prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJet1NumSignalPiZeros1prong[NpdfCounter] = fs->make<TH1F>(("TauJet1NumSignalPiZeros1prong_"+j.str()).c_str(), ("TauJet1NumSignalPiZeros1prong_"+j.str()).c_str(), 10, 0., 10.);
      _hTauJet1SignalTracksMass3prong[NpdfCounter]          = fs->make<TH1F>(("TauJet1SignalTracksMass3prong_"+j.str()).c_str(),          ("TauJet1SignalTracksMass3prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJet1SignalTracksAndGammasMass3prong[NpdfCounter] = fs->make<TH1F>(("TauJet1SignalTracksAndGammasMass3prong_"+j.str()).c_str(), ("TauJet1SignalTracksAndGammasMass3prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJet1Mass1Prong0Gamma[NpdfCounter]                = fs->make<TH1F>(("TauJet1Mass1Prong0Gamma_"+j.str()).c_str(),                ("TauJet1Mass1Prong0Gamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJet1Mass1Prong1Gamma[NpdfCounter]                = fs->make<TH1F>(("TauJet1Mass1Prong1Gamma_"+j.str()).c_str(),                ("TauJet1Mass1Prong1Gamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJet1Mass1Prong2orMoreGamma[NpdfCounter]          = fs->make<TH1F>(("TauJet1Mass1Prong2orMoreGamma_"+j.str()).c_str(),          ("TauJet1Mass1Prong2orMoreGamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJet1Mass3Prong0Gamma[NpdfCounter]                = fs->make<TH1F>(("TauJet1Mass3Prong0Gamma_"+j.str()).c_str(),                ("TauJet1Mass3Prong0Gamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJet1Mass3Prong1Gamma[NpdfCounter]                = fs->make<TH1F>(("TauJet1Mass3Prong1Gamma_"+j.str()).c_str(),                ("TauJet1Mass3Prong1Gamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJet1Mass3Prong2orMoreGamma[NpdfCounter]          = fs->make<TH1F>(("TauJet1Mass3Prong2orMoreGamma_"+j.str()).c_str(),          ("TauJet1Mass3Prong2orMoreGamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJet1H3x3OverP[NpdfCounter]          = fs->make<TH1F>(("TauJet1H3x3OverP_"+j.str()).c_str(),          ("TauJet1H3x3OverP_"+j.str()).c_str(), 100, 0., 1.);
      _hNTau2[NpdfCounter]                             = fs->make<TH1F>(("NTau2_"+j.str()).c_str(),         ("NTau2_"+j.str()).c_str(),         21, 0., 20.);
      _hTauJet2Energy[NpdfCounter]                     = fs->make<TH1F>(("TauJet2Energy_"+j.str()).c_str(), ("TauJet2Energy_"+j.str()).c_str(), 200, 0., 500.);
      _hTauJet2Pt[NpdfCounter]                         = fs->make<TH1F>(("TauJet2Pt_"+j.str()).c_str(),     ("TauJet2Pt_"+j.str()).c_str(),     200, 0., 500.);
      _hTauJet2Eta[NpdfCounter]                        = fs->make<TH1F>(("TauJet2Eta_"+j.str()).c_str(),    ("TauJet2Eta_"+j.str()).c_str(),    72, -3.6, +3.6);
      _hTauJet2Phi[NpdfCounter]                        = fs->make<TH1F>(("TauJet2Phi_"+j.str()).c_str(),    ("TauJet2Phi_"+j.str()).c_str(),    36, -TMath::Pi(), +TMath::Pi());
      _hTauJet2NumSignalTracks[NpdfCounter]            = fs->make<TH1F>(("TauJet2NumSignalTracks_"+j.str()).c_str(),           ("TauJet2NumSignalTracks_"+j.str()).c_str(), 10, 0, 10);
      _hTauJet2NumSignalGammas[NpdfCounter]            = fs->make<TH1F>(("TauJet2NumSignalGammas_"+j.str()).c_str(),           ("TauJet2NumSignalGammas_"+j.str()).c_str(), 10, 0, 10);
      _hTauJet2SeedTrackPt[NpdfCounter]                = fs->make<TH1F>(("TauJet2SeedTrackPt_"+j.str()).c_str(),               ("TauJet2SeedTrackPt_"+j.str()).c_str(),     200, 0., 500.);
      _hTauJet2SeedTrackIpSignificance[NpdfCounter]    = fs->make<TH1F>(("TauJet2SeedTrackIpSignificance_"+j.str()).c_str(),   ("TauJet2SeedTrackIpSignificance_"+j.str()).c_str(), 100, 0., 100.);
      _hTauJet2SeedTrackNhits[NpdfCounter]             = fs->make<TH1F>(("TauJet2SeedTrackNhits_"+j.str()).c_str(),            ("TauJet2SeedTrackNhits_"+j.str()).c_str(), 40, 0., 40.);
      _hTauJet2SeedTrackChi2[NpdfCounter]              = fs->make<TH1F>(("TauJet2SeedTrackChi2_"+j.str()).c_str(),             ("TauJet2SeedTrackChi2_"+j.str()).c_str(), 50, 0., 100.);
      _hTauJet2Charge[NpdfCounter]                     = fs->make<TH1F>(("TauJet2Charge_"+j.str()).c_str(),                    ("TauJet2Charge_"+j.str()).c_str(), 5, 0., 5.);
      _hTauJet2SignalTracksMass[NpdfCounter]           = fs->make<TH1F>(("TauJet2SignalTracksMass_"+j.str()).c_str(),          ("TauJet2SignalTracksMass_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJet2SignalTracksAndGammasMass[NpdfCounter]  = fs->make<TH1F>(("TauJet2SignalTracksAndGammasMass_"+j.str()).c_str(), ("TauJet2SignalTracksAndGammasMass_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJet2SignalTracksChargeFraction[NpdfCounter] = fs->make<TH1F>(("TauJet2SignalTracksChargeFraction_"+j.str()).c_str(),("TauJet2SignalTracksChargeFraction_"+j.str()).c_str(), 30, 0., 1.5);
      _hTauJet2NumIsoTracks[NpdfCounter]               = fs->make<TH1F>(("TauJet2NumIsoTracks_"+j.str()).c_str(),              ("TauJet2NumIsoTracks_"+j.str()).c_str(), 10, 0, 10);
      _hTauJet2NumIsoGammas[NpdfCounter]               = fs->make<TH1F>(("TauJet2NumIsoGammas_"+j.str()).c_str(),              ("TauJet2NumIsoGammas_"+j.str()).c_str(), 10, 0, 10);
      _hTauJet2NumIsoCands[NpdfCounter]                = fs->make<TH1F>(("TauJet2NumIsoCands_"+j.str()).c_str(),               ("TauJet2NumIsoCands_"+j.str()).c_str(), 10, 0, 10);
      _hTauJet2SumPtIsoTracks[NpdfCounter]             = fs->make<TH1F>(("TauJet2SumPtIsoTracks_"+j.str()).c_str(),            ("TauJet2SumPtIsoTracks_"+j.str()).c_str(), 100, 0, 50);
      _hTauJet2SumPtIsoGammas[NpdfCounter]             = fs->make<TH1F>(("TauJet2SumPtIsoGammas_"+j.str()).c_str(),            ("TauJet2SumPtIsoGammas_"+j.str()).c_str(), 100, 0, 50);
      _hTauJet2SumPtIso[NpdfCounter]                   = fs->make<TH1F>(("TauJet2SumPtIso_"+j.str()).c_str(),                  ("TauJet2SumPtIso_"+j.str()).c_str(), 100, 0, 50);
      _hTauJet2NumberDensity[NpdfCounter]              = fs->make<TH1F>(("TauJet2NumberDensity_"+j.str()).c_str(),                  ("TauJet2NumberDensity_"+j.str()).c_str(), 500, 0, 10);
      _hTauJet2GenTauDeltaPhi[NpdfCounter]             = fs->make<TH1F>(("TauJet2GenTauDeltaPhi_"+j.str()).c_str(),            ("TauJet2GenTauDeltaPhi_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hTauJet2GenTauDeltaEta[NpdfCounter]             = fs->make<TH1F>(("TauJet2GenTauDeltaEta_"+j.str()).c_str(),            ("TauJet2GenTauDeltaEta_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hTauJet2GenTauDeltaPt[NpdfCounter]              = fs->make<TH1F>(("TauJet2GenTauDeltaPt_"+j.str()).c_str(),             ("TauJet2GenTauDeltaPt_"+j.str()).c_str(), 500, -5, 5);
      _hTauJet2SignalTracksMass1prong[NpdfCounter]     = fs->make<TH1F>(("TauJet2SignalTracksMass1prong_"+j.str()).c_str(),    ("TauJet2SignalTracksMass1prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJet2SignalTracksAndGammasMass1prong[NpdfCounter] = fs->make<TH1F>(("TauJet2SignalTracksAndGammasMass1prong_"+j.str()).c_str(), ("TauJet2SignalTracksAndGammasMass1prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJet2SignalTracksAndPiZerosMass1prong[NpdfCounter] = fs->make<TH1F>(("TauJet2SignalTracksAndPiZerosMass1prong_"+j.str()).c_str(), ("TauJet2SignalTracksAndPiZerosMass1prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJet2NumSignalPiZeros1prong[NpdfCounter] = fs->make<TH1F>(("TauJet2NumSignalPiZeros1prong_"+j.str()).c_str(), ("TauJet2NumSignalPiZeros1prong_"+j.str()).c_str(), 10, 0., 10.);
      _hTauJet2SignalTracksMass3prong[NpdfCounter]          = fs->make<TH1F>(("TauJet2SignalTracksMass3prong_"+j.str()).c_str(),          ("TauJet2SignalTracksMass3prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJet2SignalTracksAndGammasMass3prong[NpdfCounter] = fs->make<TH1F>(("TauJet2SignalTracksAndGammasMass3prong_"+j.str()).c_str(), ("TauJet2SignalTracksAndGammasMass3prong_"+j.str()).c_str(), 50, 0., 5.);
      _hTauJet2Mass1Prong0Gamma[NpdfCounter]                = fs->make<TH1F>(("TauJet2Mass1Prong0Gamma_"+j.str()).c_str(),                ("TauJet2Mass1Prong0Gamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJet2Mass1Prong1Gamma[NpdfCounter]                = fs->make<TH1F>(("TauJet2Mass1Prong1Gamma_"+j.str()).c_str(),                ("TauJet2Mass1Prong1Gamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJet2Mass1Prong2orMoreGamma[NpdfCounter]          = fs->make<TH1F>(("TauJet2Mass1Prong2orMoreGamma_"+j.str()).c_str(),          ("TauJet2Mass1Prong2orMoreGamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJet2Mass3Prong0Gamma[NpdfCounter]                = fs->make<TH1F>(("TauJet2Mass3Prong0Gamma_"+j.str()).c_str(),                ("TauJet2Mass3Prong0Gamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJet2Mass3Prong1Gamma[NpdfCounter]                = fs->make<TH1F>(("TauJet2Mass3Prong1Gamma_"+j.str()).c_str(),                ("TauJet2Mass3Prong1Gamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJet2Mass3Prong2orMoreGamma[NpdfCounter]          = fs->make<TH1F>(("TauJet2Mass3Prong2orMoreGamma_"+j.str()).c_str(),          ("TauJet2Mass3Prong2orMoreGamma_"+j.str()).c_str(), 100, 0., 5.);
      _hTauJet2H3x3OverP[NpdfCounter]          = fs->make<TH1F>(("TauJet2H3x3OverP_"+j.str()).c_str(),          ("TauJet2H3x3OverP_"+j.str()).c_str(), 100, 0., 1.);
    }
    
    //--- book reconstruction level histograms 
    if (_FillRecoMuonHists) {
      _hNMuon1[NpdfCounter]                       = fs->make<TH1F>(("NMuon1_"+j.str()).c_str(),                   ("NMuon1_"+j.str()).c_str(), 21, 0., 20.);
      _hMuon1Energy[NpdfCounter]                  = fs->make<TH1F>(("Muon1Energy_"+j.str()).c_str(),              ("Muon1Energy_"+j.str()).c_str(), 200, 0., 500.);
      _hMuon1Pt[NpdfCounter]                      = fs->make<TH1F>(("Muon1Pt_"+j.str()).c_str(),                  ("Muon1Pt_"+j.str()).c_str(),  200, 0., 500.);
      _hMuon1Eta[NpdfCounter]                     = fs->make<TH1F>(("Muon1Eta_"+j.str()).c_str(),                 ("Muon1Eta_"+j.str()).c_str(), 72, -3.6, +3.6);
      _hMuon1Phi[NpdfCounter]                     = fs->make<TH1F>(("Muon1Phi_"+j.str()).c_str(),                 ("Muon1Phi_"+j.str()).c_str(), 36, -TMath::Pi(), +TMath::Pi());
      _hMuon1TrackIso[NpdfCounter]                = fs->make<TH1F>(("Muon1TrackIso_"+j.str()).c_str(),            ("Muon1TrackIso_"+j.str()).c_str(), 100, 0, 50);
      _hMuon1EcalIso[NpdfCounter]                 = fs->make<TH1F>(("Muon1EcalIso_"+j.str()).c_str(),             ("Muon1EcalIso_"+j.str()).c_str(), 100, 0, 50);
      _hMuon1Iso[NpdfCounter]                     = fs->make<TH1F>(("Muon1Iso_"+j.str()).c_str(),                 ("Muon1Iso_"+j.str()).c_str(), 100, 0, 50);
      _hMuon1Ip[NpdfCounter]                      = fs->make<TH1F>(("Muon1Ip_"+j.str()).c_str(),                  ("Muon1Ip_"+j.str()).c_str(), 500, -1, +1);
      _hMuon1IpSignificance[NpdfCounter]          = fs->make<TH1F>(("Muon1IpSignificance_"+j.str()).c_str(),      ("Muon1IpSignificance_"+j.str()).c_str(), 100, 0., 100.);
      _hMuon1GenMuonDeltaPhi[NpdfCounter]         = fs->make<TH1F>(("Muon1GenMuonDeltaPhi_"+j.str()).c_str(),     ("Muon1GenMuonDeltaPhi_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hMuon1GenMuonDeltaEta[NpdfCounter]         = fs->make<TH1F>(("Muon1GenMuonDeltaEta_"+j.str()).c_str(),     ("Muon1GenMuonDeltaEta_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hMuon1GenMuonDeltaPt[NpdfCounter]          = fs->make<TH1F>(("Muon1GenMuonDeltaPt_"+j.str()).c_str(),      ("Muon1GenMuonDeltaPt_"+j.str()).c_str(), 500, -5, 5);
      _hMuon1CaloCompatibility[NpdfCounter]       = fs->make<TH1F>(("Muon1CaloCompatibility_"+j.str()).c_str(),   ("Muon1CaloCompatibility_"+j.str()).c_str(), 102, 0.0, 1.02);
      _hMuon1SegmentCompatibility[NpdfCounter]    = fs->make<TH1F>(("Muon1SegmentCompatibility_"+j.str()).c_str(),("Muon1SegmentCompatibility_"+j.str()).c_str(), 102, 0.0, 1.02);
      _hMuon1AntiPion[NpdfCounter]                = fs->make<TH1F>(("Muon1AntiPion_"+j.str()).c_str(),            ("Muon1AntiPion_"+j.str()).c_str(), 202, 0.0, 2.02);
      _hMuon1CaloCompatibilityVsSegmentCompatibility[NpdfCounter] = fs->make<TH2F>(("Muon1CaloCompatibilityVsSegmentCompatibility_"+j.str()).c_str(), ("Muon1CaloCompatibilityVsSegmentCompatibility_"+j.str()).c_str(), 102, 0, 1.02, 102, 0, 1.02);      
      _hNMuon2[NpdfCounter]                       = fs->make<TH1F>(("NMuon2_"+j.str()).c_str(),                   ("NMuon2_"+j.str()).c_str(), 21, 0., 20.);
      _hMuon2Energy[NpdfCounter]                  = fs->make<TH1F>(("Muon2Energy_"+j.str()).c_str(),              ("Muon2Energy_"+j.str()).c_str(), 200, 0., 500.);
      _hMuon2Pt[NpdfCounter]                      = fs->make<TH1F>(("Muon2Pt_"+j.str()).c_str(),                  ("Muon2Pt_"+j.str()).c_str(),  200, 0., 500.);
      _hMuon2Eta[NpdfCounter]                     = fs->make<TH1F>(("Muon2Eta_"+j.str()).c_str(),                 ("Muon2Eta_"+j.str()).c_str(), 72, -3.6, +3.6);
      _hMuon2Phi[NpdfCounter]                     = fs->make<TH1F>(("Muon2Phi_"+j.str()).c_str(),                 ("Muon2Phi_"+j.str()).c_str(), 36, -TMath::Pi(), +TMath::Pi());
      _hMuon2TrackIso[NpdfCounter]                = fs->make<TH1F>(("Muon2TrackIso_"+j.str()).c_str(),            ("Muon2TrackIso_"+j.str()).c_str(), 100, 0, 50);
      _hMuon2EcalIso[NpdfCounter]                 = fs->make<TH1F>(("Muon2EcalIso_"+j.str()).c_str(),             ("Muon2EcalIso_"+j.str()).c_str(), 100, 0, 50);
      _hMuon2Iso[NpdfCounter]                     = fs->make<TH1F>(("Muon2Iso_"+j.str()).c_str(),                 ("Muon2Iso_"+j.str()).c_str(), 100, 0, 50);
      _hMuon2Ip[NpdfCounter]                      = fs->make<TH1F>(("Muon2Ip_"+j.str()).c_str(),                  ("Muon2Ip_"+j.str()).c_str(), 500, -1, +1);
      _hMuon2IpSignificance[NpdfCounter]          = fs->make<TH1F>(("Muon2IpSignificance_"+j.str()).c_str(),      ("Muon2IpSignificance_"+j.str()).c_str(), 100, 0., 100.);
      _hMuon2GenMuonDeltaPhi[NpdfCounter]         = fs->make<TH1F>(("Muon2GenMuonDeltaPhi_"+j.str()).c_str(),     ("Muon2GenMuonDeltaPhi_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hMuon2GenMuonDeltaEta[NpdfCounter]         = fs->make<TH1F>(("Muon2GenMuonDeltaEta_"+j.str()).c_str(),     ("Muon2GenMuonDeltaEta_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hMuon2GenMuonDeltaPt[NpdfCounter]          = fs->make<TH1F>(("Muon2GenMuonDeltaPt_"+j.str()).c_str(),      ("Muon2GenMuonDeltaPt_"+j.str()).c_str(), 500, -5, 5);
      _hMuon2CaloCompatibility[NpdfCounter]       = fs->make<TH1F>(("Muon2CaloCompatibility_"+j.str()).c_str(),   ("Muon2CaloCompatibility_"+j.str()).c_str(), 102, 0.0, 1.02);
      _hMuon2SegmentCompatibility[NpdfCounter]    = fs->make<TH1F>(("Muon2SegmentCompatibility_"+j.str()).c_str(),("Muon2SegmentCompatibility_"+j.str()).c_str(), 102, 0.0, 1.02);
      _hMuon2AntiPion[NpdfCounter]                = fs->make<TH1F>(("Muon2AntiPion_"+j.str()).c_str(),            ("Muon2AntiPion_"+j.str()).c_str(), 202, 0.0, 2.02);
      _hMuon2CaloCompatibilityVsSegmentCompatibility[NpdfCounter] = fs->make<TH2F>(("Muon2CaloCompatibilityVsSegmentCompatibility_"+j.str()).c_str(), ("Muon2CaloCompatibilityVsSegmentCompatibility_"+j.str()).c_str(), 102, 0, 1.02, 102, 0, 1.02);      
    }
    
    if (_FillRecoElectronHists) {
      _hNElectron1[NpdfCounter]                   = fs->make<TH1F>(("NElectron1_"+j.str()).c_str(),                  ("NElectron1_"+j.str()).c_str(), 21, 0., 20.);
      _hElectron1Energy[NpdfCounter]              = fs->make<TH1F>(("Electron1Energy_"+j.str()).c_str(),             ("Electron1Energy_"+j.str()).c_str(), 200, 0., 500.);
      _hElectron1Pt[NpdfCounter]                  = fs->make<TH1F>(("Electron1Pt_"+j.str()).c_str(),                 ("Electron1Pt_"+j.str()).c_str(), 200, 0., 500.);
      _hElectron1Eta[NpdfCounter]                 = fs->make<TH1F>(("Electron1Eta_"+j.str()).c_str(),                ("Electron1Eta_"+j.str()).c_str(), 72, -3.6, +3.6);
      _hElectron1Phi[NpdfCounter]                 = fs->make<TH1F>(("Electron1Phi_"+j.str()).c_str(),                ("Electron1Phi_"+j.str()).c_str(), 36, -TMath::Pi(), +TMath::Pi());
      _hElectron1TrackIso[NpdfCounter]            = fs->make<TH1F>(("Electron1TrackIso_"+j.str()).c_str(),           ("Electron1TrackIso_"+j.str()).c_str(), 100, 0, 50);
      _hElectron1EcalIso[NpdfCounter]             = fs->make<TH1F>(("Electron1EcalIso_"+j.str()).c_str(),            ("Electron1EcalIso_"+j.str()).c_str(), 100, 0, 50);
      _hElectron1Ip[NpdfCounter]                  = fs->make<TH1F>(("Electron1Ip_"+j.str()).c_str(),                 ("Electron1Ip_"+j.str()).c_str(), 500, -1, +1);
      _hElectron1EoverP[NpdfCounter]              = fs->make<TH1F>(("Electron1EoverP_"+j.str()).c_str(),             ("Electron1EoverP_"+j.str()).c_str(), 60, 0, +3);
      _hElectron1HoverEm[NpdfCounter]             = fs->make<TH1F>(("Electron1HoverEm_"+j.str()).c_str(),            ("Electron1HoverEm_"+j.str()).c_str(), 300, 0, +3);
      _hElectron1Classification[NpdfCounter]      = fs->make<TH1F>(("Electron1Classification_"+j.str()).c_str(),     ("Electron1Classification_"+j.str()).c_str(), 200, 0, 200);
      _hElectron1GenElectronDeltaPhi[NpdfCounter] = fs->make<TH1F>(("Electron1GenElectronDeltaPhi_"+j.str()).c_str(),("Electron1GenElectronDeltaPhi_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hElectron1GenElectronDeltaEta[NpdfCounter] = fs->make<TH1F>(("Electron1GenElectronDeltaEta_"+j.str()).c_str(),("Electron1GenElectronDeltaEta_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hElectron1GenElectronDeltaPt[NpdfCounter]  = fs->make<TH1F>(("Electron1GenElectronDeltaPt_"+j.str()).c_str(), ("Electron1GenElectronDeltaPt_"+j.str()).c_str(), 500, -5, 5);
      _hElectron1EcalDriven[NpdfCounter]          = fs->make<TH1F>(("Electron1EcalDriven_"+j.str()).c_str(),         ("Electron1EcalDriven_"+j.str()).c_str(), 2, 0., 2.);
      _hElectron1TrackerDriven[NpdfCounter]       = fs->make<TH1F>(("Electron1TrackerDriven_"+j.str()).c_str(),      ("Electron1TrackerDriven_"+j.str()).c_str(), 2, 0., 2.);
      _hElectron1HoverEm[NpdfCounter]             = fs->make<TH1F>(("Electron1HoverEm_"+j.str()).c_str(),            ("Electron1HoverEm_"+j.str()).c_str(), 100, 0., 0.5);
      _hElectron1EESigmaIEtaIEta[NpdfCounter]     = fs->make<TH1F>(("Electron1EESigmaIEtaIEta_"+j.str()).c_str(),    ("Electron1EESigmaIEtaIEta_"+j.str()).c_str(), 100, 0., 0.05);
      _hElectron1EEDEta[NpdfCounter]              = fs->make<TH1F>(("Electron1EEDEta_"+j.str()).c_str(),             ("Electron1EEDEta_"+j.str()).c_str(), 100, -0.05, 0.05);
      _hElectron1EEDPhi[NpdfCounter]              = fs->make<TH1F>(("Electron1EEDPhi_"+j.str()).c_str(),             ("Electron1EEDPhi_"+j.str()).c_str(), 100, -0.2, 0.2);
      _hElectron1EBSigmaIEtaIEta[NpdfCounter]     = fs->make<TH1F>(("Electron1EBSigmaIEtaIEta_"+j.str()).c_str(),    ("Electron1EBSigmaIEtaIEta_"+j.str()).c_str(), 100, 0., 0.05);
      _hElectron1EBDEta[NpdfCounter]              = fs->make<TH1F>(("Electron1EBDEta_"+j.str()).c_str(),             ("Electron1EBDEta_"+j.str()).c_str(), 100, -0.05, 0.05);
      _hElectron1EBDPhi[NpdfCounter]              = fs->make<TH1F>(("Electron1EBDPhi_"+j.str()).c_str(),             ("Electron1EBDPhi_"+j.str()).c_str(), 100, -0.2, 0.2);
      _hElectron1EB2by5Over5by5[NpdfCounter]      = fs->make<TH1F>(("Electron1EB2by5Over5by5_"+j.str()).c_str(),     ("Electron1EB2by5Over5by5_"+j.str()).c_str(), 100, 0., 1.);
      _hElectron1EB1by5Over5by5[NpdfCounter]      = fs->make<TH1F>(("Electron1EB1by5Over5by5_"+j.str()).c_str(),     ("Electron1EB1by5Over5by5_"+j.str()).c_str(), 100, 0., 1.);
      _hElectron1MissingHits[NpdfCounter]         = fs->make<TH1F>(("Electron1MissingHits_"+j.str()).c_str(),        ("Electron1MissingHits_"+j.str()).c_str(), 10, 0., 10.);
      _hNElectron2[NpdfCounter]                   = fs->make<TH1F>(("NElectron2_"+j.str()).c_str(),                  ("NElectron2_"+j.str()).c_str(), 21, 0., 20.);
      _hElectron2Energy[NpdfCounter]              = fs->make<TH1F>(("Electron2Energy_"+j.str()).c_str(),             ("Electron2Energy_"+j.str()).c_str(), 200, 0., 500.);
      _hElectron2Pt[NpdfCounter]                  = fs->make<TH1F>(("Electron2Pt_"+j.str()).c_str(),                 ("Electron2Pt_"+j.str()).c_str(), 200, 0., 500.);
      _hElectron2Eta[NpdfCounter]                 = fs->make<TH1F>(("Electron2Eta_"+j.str()).c_str(),                ("Electron2Eta_"+j.str()).c_str(), 72, -3.6, +3.6);
      _hElectron2Phi[NpdfCounter]                 = fs->make<TH1F>(("Electron2Phi_"+j.str()).c_str(),                ("Electron2Phi_"+j.str()).c_str(), 36, -TMath::Pi(), +TMath::Pi());
      _hElectron2TrackIso[NpdfCounter]            = fs->make<TH1F>(("Electron2TrackIso_"+j.str()).c_str(),           ("Electron2TrackIso_"+j.str()).c_str(), 100, 0, 50);
      _hElectron2EcalIso[NpdfCounter]             = fs->make<TH1F>(("Electron2EcalIso_"+j.str()).c_str(),            ("Electron2EcalIso_"+j.str()).c_str(), 100, 0, 50);
      _hElectron2Ip[NpdfCounter]                  = fs->make<TH1F>(("Electron2Ip_"+j.str()).c_str(),                 ("Electron2Ip_"+j.str()).c_str(), 500, -1, +1);
      _hElectron2EoverP[NpdfCounter]              = fs->make<TH1F>(("Electron2EoverP_"+j.str()).c_str(),             ("Electron2EoverP_"+j.str()).c_str(), 60, 0, +3);
      _hElectron2HoverEm[NpdfCounter]             = fs->make<TH1F>(("Electron2HoverEm_"+j.str()).c_str(),            ("Electron2HoverEm_"+j.str()).c_str(), 300, 0, +3);
      _hElectron2Classification[NpdfCounter]      = fs->make<TH1F>(("Electron2Classification_"+j.str()).c_str(),     ("Electron2Classification_"+j.str()).c_str(), 200, 0, 200);
      _hElectron2GenElectronDeltaPhi[NpdfCounter] = fs->make<TH1F>(("Electron2GenElectronDeltaPhi_"+j.str()).c_str(),("Electron2GenElectronDeltaPhi_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hElectron2GenElectronDeltaEta[NpdfCounter] = fs->make<TH1F>(("Electron2GenElectronDeltaEta_"+j.str()).c_str(),("Electron2GenElectronDeltaEta_"+j.str()).c_str(), 800, -0.2, 0.2);
      _hElectron2GenElectronDeltaPt[NpdfCounter]  = fs->make<TH1F>(("Electron2GenElectronDeltaPt_"+j.str()).c_str(), ("Electron2GenElectronDeltaPt_"+j.str()).c_str(), 500, -5, 5);
      _hElectron2EcalDriven[NpdfCounter]          = fs->make<TH1F>(("Electron2EcalDriven_"+j.str()).c_str(),         ("Electron2EcalDriven_"+j.str()).c_str(), 2, 0., 2.);
      _hElectron2TrackerDriven[NpdfCounter]       = fs->make<TH1F>(("Electron2TrackerDriven_"+j.str()).c_str(),      ("Electron2TrackerDriven_"+j.str()).c_str(), 2, 0., 2.);
//      _hElectron2IsZee[NpdfCounter]               = fs->make<TH1F>(("Electron2IsZee_"+j.str()).c_str(),              ("Electron2IsZee_"+j.str()).c_str(), 2, 0., 2.);
      _hElectron2HoverEm[NpdfCounter]             = fs->make<TH1F>(("Electron2HoverEm_"+j.str()).c_str(),            ("Electron2HoverEm_"+j.str()).c_str(), 100, 0., 0.5);
      _hElectron2EESigmaIEtaIEta[NpdfCounter]     = fs->make<TH1F>(("Electron2EESigmaIEtaIEta_"+j.str()).c_str(),    ("Electron2EESigmaIEtaIEta_"+j.str()).c_str(), 100, 0., 0.05);
      _hElectron2EEDEta[NpdfCounter]              = fs->make<TH1F>(("Electron2EEDEta_"+j.str()).c_str(),             ("Electron2EEDEta_"+j.str()).c_str(), 100, -0.05, 0.05);
      _hElectron2EEDPhi[NpdfCounter]              = fs->make<TH1F>(("Electron2EEDPhi_"+j.str()).c_str(),             ("Electron2EEDPhi_"+j.str()).c_str(), 100, -0.2, 0.2);
      _hElectron2EBSigmaIEtaIEta[NpdfCounter]     = fs->make<TH1F>(("Electron2EBSigmaIEtaIEta_"+j.str()).c_str(),    ("Electron2EBSigmaIEtaIEta_"+j.str()).c_str(), 100, 0., 0.05);
      _hElectron2EBDEta[NpdfCounter]              = fs->make<TH1F>(("Electron2EBDEta_"+j.str()).c_str(),             ("Electron2EBDEta_"+j.str()).c_str(), 100, -0.05, 0.05);
      _hElectron2EBDPhi[NpdfCounter]              = fs->make<TH1F>(("Electron2EBDPhi_"+j.str()).c_str(),             ("Electron2EBDPhi_"+j.str()).c_str(), 100, -0.2, 0.2);
      _hElectron2EB2by5Over5by5[NpdfCounter]      = fs->make<TH1F>(("Electron2EB2by5Over5by5_"+j.str()).c_str(),     ("Electron2EB2by5Over5by5_"+j.str()).c_str(), 100, 0., 1.);
      _hElectron2EB1by5Over5by5[NpdfCounter]      = fs->make<TH1F>(("Electron2EB1by5Over5by5_"+j.str()).c_str(),     ("Electron2EB1by5Over5by5_"+j.str()).c_str(), 100, 0., 1.);
      _hElectron2MissingHits[NpdfCounter]         = fs->make<TH1F>(("Electron2MissingHits_"+j.str()).c_str(),        ("Electron2MissingHits_"+j.str()).c_str(), 10, 0., 10.);
    }
    
    if (_FillRecoJetHists) {
      _hNJet[NpdfCounter]        = fs->make<TH1F>(("NJet_"+j.str()).c_str(),      ("NJet_"+j.str()).c_str(), 21, 0., 20.);
      _hNBJet[NpdfCounter]        = fs->make<TH1F>(("NBJet_"+j.str()).c_str(),      ("NBJet_"+j.str()).c_str(), 21, 0., 20.);
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
	_hMuon1PtVsTau1Pt[NpdfCounter]                   = fs->make<TH2F>(("Muon1PtVsTau1Pt_"+j.str()).c_str(),   ("Muon1PtVsTau1Pt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hMuon1Tau1DeltaR[NpdfCounter]                   = fs->make<TH1F>(("Muon1Tau1DeltaR_"+j.str()).c_str(),   ("Muon1Tau1DeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hMuon1Tau1DeltaPtDivSumPt[NpdfCounter]          = fs->make<TH1F>(("Muon1Tau1DeltaPtDivSumPt_"+j.str()).c_str(), ("Muon1Tau1DeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hMuon1Tau1DeltaPt[NpdfCounter]                  = fs->make<TH1F>(("Muon1Tau1DeltaPt_"+j.str()).c_str(), ("Muon1Tau1DeltaPt_"+j.str()).c_str(), 100, 0, 1000);
	_hMuon1Tau1_Muon1MetMt[NpdfCounter]                       = fs->make<TH1F>(("Muon1Tau1_Muon1MetMt_"+j.str()).c_str(),       ("Muon1Tau1_Muon1MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hMuon1Tau1_Tau1MetMt[NpdfCounter]                        = fs->make<TH1F>(("Muon1Tau1_Tau1MetMt_"+j.str()).c_str(),        ("Muon1Tau1_Tau1MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hMuon1Tau1OSLS[NpdfCounter]                     = fs->make<TH1F>(("Muon1Tau1OSLS_"+j.str()).c_str(),     ("Muon1Tau1OSLS_"+j.str()).c_str(), 20, -10, 10);
	_hMuon1Tau1CosDphi[NpdfCounter]                  = fs->make<TH1F>(("Muon1Tau1CosDphi_"+j.str()).c_str(),  ("Muon1Tau1CosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hMuon1Tau1_Muon1MetDeltaPhi[NpdfCounter]                 = fs->make<TH1F>(("Muon1Tau1_Muon1MetDeltaPhi_"+j.str()).c_str(), ("Muon1Tau1_Muon1MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hMuon1Tau1_Tau1MetDeltaPhi[NpdfCounter]                  = fs->make<TH1F>(("Muon1Tau1_Tau1MetDeltaPhi_"+j.str()).c_str(), ("Muon1Tau1_Tau1MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hMuon1MetDeltaPhiVsMuon1Tau1CosDphi[NpdfCounter] = fs->make<TH2F>(("Muon1MetDeltaPhiVsMuon1Tau1CosDphi_"+j.str()).c_str(), ("Muon1MetDeltaPhiVsMuon1Tau1CosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
	_hMuon1PtVsTau2Pt[NpdfCounter]                   = fs->make<TH2F>(("Muon1PtVsTau2Pt_"+j.str()).c_str(),   ("Muon1PtVsTau2Pt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hMuon1Tau2DeltaR[NpdfCounter]                   = fs->make<TH1F>(("Muon1Tau2DeltaR_"+j.str()).c_str(),   ("Muon1Tau2DeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hMuon1Tau2DeltaPtDivSumPt[NpdfCounter]          = fs->make<TH1F>(("Muon1Tau2DeltaPtDivSumPt_"+j.str()).c_str(), ("Muon1Tau2DeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hMuon1Tau2DeltaPt[NpdfCounter]                  = fs->make<TH1F>(("Muon1Tau2DeltaPt_"+j.str()).c_str(), ("Muon1Tau2DeltaPt_"+j.str()).c_str(), 100, 0, 1000);
	_hMuon1Tau2_Muon1MetMt[NpdfCounter]                       = fs->make<TH1F>(("Muon1Tau2_Muon1MetMt_"+j.str()).c_str(),       ("Muon1Tau2_Muon1MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hMuon1Tau2_Tau2MetMt[NpdfCounter]                        = fs->make<TH1F>(("Muon1Tau2_Tau2MetMt_"+j.str()).c_str(),        ("Muon1Tau2_Tau2MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hMuon1Tau2OSLS[NpdfCounter]                     = fs->make<TH1F>(("Muon1Tau2OSLS_"+j.str()).c_str(),     ("Muon1Tau2OSLS_"+j.str()).c_str(), 20, -10, 10);
	_hMuon1Tau2CosDphi[NpdfCounter]                  = fs->make<TH1F>(("Muon1Tau2CosDphi_"+j.str()).c_str(),  ("Muon1Tau2CosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hMuon1Tau2_Muon1MetDeltaPhi[NpdfCounter]                 = fs->make<TH1F>(("Muon1Tau2_Muon1MetDeltaPhi_"+j.str()).c_str(), ("Muon1Tau2_Muon1MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hMuon1Tau2_Tau2MetDeltaPhi[NpdfCounter]                  = fs->make<TH1F>(("Muon1Tau2_Tau2MetDeltaPhi_"+j.str()).c_str(), ("Muon1Tau2_Tau2MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hMuon1MetDeltaPhiVsMuon1Tau2CosDphi[NpdfCounter] = fs->make<TH2F>(("Muon1MetDeltaPhiVsMuon1Tau2CosDphi_"+j.str()).c_str(), ("Muon1MetDeltaPhiVsMuon1Tau2CosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
	_hMuon2PtVsTau1Pt[NpdfCounter]                   = fs->make<TH2F>(("Muon2PtVsTau1Pt_"+j.str()).c_str(),   ("Muon2PtVsTau1Pt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hMuon2Tau1DeltaR[NpdfCounter]                   = fs->make<TH1F>(("Muon2Tau1DeltaR_"+j.str()).c_str(),   ("Muon2Tau1DeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hMuon2Tau1DeltaPtDivSumPt[NpdfCounter]          = fs->make<TH1F>(("Muon2Tau1DeltaPtDivSumPt_"+j.str()).c_str(), ("Muon2Tau1DeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hMuon2Tau1DeltaPt[NpdfCounter]                  = fs->make<TH1F>(("Muon2Tau1DeltaPt_"+j.str()).c_str(), ("Muon2Tau1DeltaPt_"+j.str()).c_str(), 100, 0, 1000);
	_hMuon2Tau1_Muon2MetMt[NpdfCounter]                       = fs->make<TH1F>(("Muon2Tau1_Muon2MetMt_"+j.str()).c_str(),       ("Muon2Tau1_Muon2MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hMuon2Tau1_Tau1MetMt[NpdfCounter]                        = fs->make<TH1F>(("Muon2Tau1_Tau1MetMt_"+j.str()).c_str(),        ("Muon2Tau1_Tau1MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hMuon2Tau1OSLS[NpdfCounter]                     = fs->make<TH1F>(("Muon2Tau1OSLS_"+j.str()).c_str(),     ("Muon2Tau1OSLS_"+j.str()).c_str(), 20, -10, 10);
	_hMuon2Tau1CosDphi[NpdfCounter]                  = fs->make<TH1F>(("Muon2Tau1CosDphi_"+j.str()).c_str(),  ("Muon2Tau1CosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hMuon2Tau1_Muon2MetDeltaPhi[NpdfCounter]                 = fs->make<TH1F>(("Muon2Tau1_Muon2MetDeltaPhi_"+j.str()).c_str(), ("Muon2Tau1_Muon2MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hMuon2Tau1_Tau1MetDeltaPhi[NpdfCounter]                  = fs->make<TH1F>(("Muon2Tau1_Tau1MetDeltaPhi_"+j.str()).c_str(), ("Muon2Tau1_Tau1MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hMuon2MetDeltaPhiVsMuon2Tau1CosDphi[NpdfCounter] = fs->make<TH2F>(("Muon2MetDeltaPhiVsMuon2Tau1CosDphi_"+j.str()).c_str(), ("Muon2MetDeltaPhiVsMuon2Tau1CosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
	_hMuon2PtVsTau2Pt[NpdfCounter]                   = fs->make<TH2F>(("Muon2PtVsTau2Pt_"+j.str()).c_str(),   ("Muon2PtVsTau2Pt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hMuon2Tau2DeltaR[NpdfCounter]                   = fs->make<TH1F>(("Muon2Tau2DeltaR_"+j.str()).c_str(),   ("Muon2Tau2DeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hMuon2Tau2DeltaPtDivSumPt[NpdfCounter]          = fs->make<TH1F>(("Muon2Tau2DeltaPtDivSumPt_"+j.str()).c_str(), ("Muon2Tau2DeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hMuon2Tau2DeltaPt[NpdfCounter]                  = fs->make<TH1F>(("Muon2Tau2DeltaPt_"+j.str()).c_str(), ("Muon2Tau2DeltaPt_"+j.str()).c_str(), 100, 0, 1000);
	_hMuon2Tau2_Muon2MetMt[NpdfCounter]                       = fs->make<TH1F>(("Muon2Tau2_Muon2MetMt_"+j.str()).c_str(),       ("Muon2Tau2_Muon2MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hMuon2Tau2_Tau2MetMt[NpdfCounter]                        = fs->make<TH1F>(("Muon2Tau2_Tau2MetMt_"+j.str()).c_str(),        ("Muon2Tau2_Tau2MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hMuon2Tau2OSLS[NpdfCounter]                     = fs->make<TH1F>(("Muon2Tau2OSLS_"+j.str()).c_str(),     ("Muon2Tau2OSLS_"+j.str()).c_str(), 20, -10, 10);
	_hMuon2Tau2CosDphi[NpdfCounter]                  = fs->make<TH1F>(("Muon2Tau2CosDphi_"+j.str()).c_str(),  ("Muon2Tau2CosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hMuon2Tau2_Muon2MetDeltaPhi[NpdfCounter]                 = fs->make<TH1F>(("Muon2Tau2_Muon2MetDeltaPhi_"+j.str()).c_str(), ("Muon2Tau2_Muon2MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hMuon2Tau2_Tau2MetDeltaPhi[NpdfCounter]                  = fs->make<TH1F>(("Muon2Tau2_Tau2MetDeltaPhi_"+j.str()).c_str(), ("Muon2Tau2_Tau2MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hMuon2MetDeltaPhiVsMuon2Tau2CosDphi[NpdfCounter] = fs->make<TH2F>(("Muon2MetDeltaPhiVsMuon2Tau2CosDphi_"+j.str()).c_str(), ("Muon2MetDeltaPhiVsMuon2Tau2CosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
        _hElectron1Tau1_Electron1IsZee[NpdfCounter]               = fs->make<TH1F>(("Electron1Tau1_Electron1IsZee_"+j.str()).c_str(),              ("Electron1Tau1_Electron1IsZee_"+j.str()).c_str(), 2, 0., 2.);
        _hElectron1Tau2_Electron1IsZee[NpdfCounter]               = fs->make<TH1F>(("Electron1Tau2_Electron1IsZee_"+j.str()).c_str(),              ("Electron1Tau2_Electron1IsZee_"+j.str()).c_str(), 2, 0., 2.);
        _hElectron2Tau1_Electron2IsZee[NpdfCounter]               = fs->make<TH1F>(("Electron2Tau1_Electron2IsZee_"+j.str()).c_str(),              ("Electron2Tau1_Electron2IsZee_"+j.str()).c_str(), 2, 0., 2.);
        _hElectron2Tau2_Electron2IsZee[NpdfCounter]               = fs->make<TH1F>(("Electron2Tau2_Electron2IsZee_"+j.str()).c_str(),              ("Electron2Tau2_Electron2IsZee_"+j.str()).c_str(), 2, 0., 2.);
	_hElectron1PtVsTau1Pt[NpdfCounter]          = fs->make<TH2F>(("Electron1PtVsTau1Pt_"+j.str()).c_str(),          ("Electron1PtVsTau1Pt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hElectron1Tau1DeltaR[NpdfCounter]          = fs->make<TH1F>(("Electron1Tau1DeltaR_"+j.str()).c_str(),          ("Electron1Tau1DeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hElectron1Tau1DeltaPtDivSumPt[NpdfCounter] = fs->make<TH1F>(("Electron1Tau1DeltaPtDivSumPt_"+j.str()).c_str(), ("Electron1Tau1DeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hElectron1Tau1DeltaPt[NpdfCounter]         = fs->make<TH1F>(("Electron1Tau1DeltaPt_"+j.str()).c_str(), ("Electron1Tau1DeltaPt_"+j.str()).c_str(), 100, 0, 1000);
	_hElectron1Tau1_Electron1MetMt[NpdfCounter]              = fs->make<TH1F>(("Electron1Tau1_Electron1MetMt_"+j.str()).c_str(),              ("Electron1Tau1_Electron1MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hElectron1Tau1_Tau1MetMt[NpdfCounter]                   = fs->make<TH1F>(("Electron1Tau1_Tau1MetMt_"+j.str()).c_str(),                   ("Electron1Tau1_Tau1MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hElectron1Tau1OSLS[NpdfCounter]            = fs->make<TH1F>(("Electron1Tau1OSLS_"+j.str()).c_str(),            ("Electron1Tau1OSLS_"+j.str()).c_str(), 20, -10, 10);
	_hElectron1Tau1CosDphi[NpdfCounter]         = fs->make<TH1F>(("Electron1Tau1CosDphi_"+j.str()).c_str(),         ("Electron1Tau1CosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hElectron1Tau1_Electron1MetDeltaPhi[NpdfCounter]        = fs->make<TH1F>(("Electron1Tau1_Electron1MetDeltaPhi_"+j.str()).c_str(),        ("Electron1Tau1_Electron1MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hElectron1Tau1_Tau1MetDeltaPhi[NpdfCounter]             = fs->make<TH1F>(("Electron1Tau1_Tau1MetDeltaPhi_"+j.str()).c_str(),             ("Electron1Tau1_Tau1MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hElectron1MetDeltaPhiVsElectron1Tau1CosDphi[NpdfCounter] = fs->make<TH2F>(("Electron1MetDeltaPhiVsElectron1Tau1CosDphi_"+j.str()).c_str(), ("Electron1MetDeltaPhiVsElectron1Tau1CosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
	_hElectron1PtVsTau2Pt[NpdfCounter]          = fs->make<TH2F>(("Electron1PtVsTau2Pt_"+j.str()).c_str(),          ("Electron1PtVsTau2Pt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hElectron1Tau2DeltaR[NpdfCounter]          = fs->make<TH1F>(("Electron1Tau2DeltaR_"+j.str()).c_str(),          ("Electron1Tau2DeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hElectron1Tau2DeltaPtDivSumPt[NpdfCounter] = fs->make<TH1F>(("Electron1Tau2DeltaPtDivSumPt_"+j.str()).c_str(), ("Electron1Tau2DeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hElectron1Tau2DeltaPt[NpdfCounter]         = fs->make<TH1F>(("Electron1Tau2DeltaPt_"+j.str()).c_str(), ("Electron1Tau2DeltaPt_"+j.str()).c_str(), 100, 0, 1000);
	_hElectron1Tau2_Electron1MetMt[NpdfCounter]              = fs->make<TH1F>(("Electron1Tau2_Electron1MetMt_"+j.str()).c_str(),              ("Electron1Tau2_Electron1MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hElectron1Tau2_Tau2MetMt[NpdfCounter]                   = fs->make<TH1F>(("Electron1Tau2_Tau2MetMt_"+j.str()).c_str(),                   ("Electron1Tau2_Tau2MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hElectron1Tau2OSLS[NpdfCounter]            = fs->make<TH1F>(("Electron1Tau2OSLS_"+j.str()).c_str(),            ("Electron1Tau2OSLS_"+j.str()).c_str(), 20, -10, 10);
	_hElectron1Tau2CosDphi[NpdfCounter]         = fs->make<TH1F>(("Electron1Tau2CosDphi_"+j.str()).c_str(),         ("Electron1Tau2CosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hElectron1Tau2_Electron1MetDeltaPhi[NpdfCounter]        = fs->make<TH1F>(("Electron1Tau2_Electron1MetDeltaPhi_"+j.str()).c_str(),        ("Electron1Tau2_Electron1MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hElectron1Tau2_Tau2MetDeltaPhi[NpdfCounter]             = fs->make<TH1F>(("Electron1Tau2_Tau2MetDeltaPhi_"+j.str()).c_str(),             ("Electron1Tau2_Tau2MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hElectron1MetDeltaPhiVsElectron1Tau2CosDphi[NpdfCounter] = fs->make<TH2F>(("Electron1MetDeltaPhiVsElectron1Tau2CosDphi_"+j.str()).c_str(), ("Electron1MetDeltaPhiVsElectron1Tau2CosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
	_hElectron2PtVsTau1Pt[NpdfCounter]          = fs->make<TH2F>(("Electron2PtVsTau1Pt_"+j.str()).c_str(),          ("Electron2PtVsTau1Pt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hElectron2Tau1DeltaR[NpdfCounter]          = fs->make<TH1F>(("Electron2Tau1DeltaR_"+j.str()).c_str(),          ("Electron2Tau1DeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hElectron2Tau1DeltaPtDivSumPt[NpdfCounter] = fs->make<TH1F>(("Electron2Tau1DeltaPtDivSumPt_"+j.str()).c_str(), ("Electron2Tau1DeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hElectron2Tau1DeltaPt[NpdfCounter]         = fs->make<TH1F>(("Electron2Tau1DeltaPt_"+j.str()).c_str(), ("Electron2Tau1DeltaPt_"+j.str()).c_str(), 100, 0, 1000);
	_hElectron2Tau1_Electron2MetMt[NpdfCounter]              = fs->make<TH1F>(("Electron2Tau1_Electron2MetMt_"+j.str()).c_str(),              ("Electron2Tau1_Electron2MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hElectron2Tau1_Tau1MetMt[NpdfCounter]                   = fs->make<TH1F>(("Electron2Tau1_Tau1MetMt_"+j.str()).c_str(),                   ("Electron2Tau1_Tau1MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hElectron2Tau1OSLS[NpdfCounter]            = fs->make<TH1F>(("Electron2Tau1OSLS_"+j.str()).c_str(),            ("Electron2Tau1OSLS_"+j.str()).c_str(), 20, -10, 10);
	_hElectron2Tau1CosDphi[NpdfCounter]         = fs->make<TH1F>(("Electron2Tau1CosDphi_"+j.str()).c_str(),         ("Electron2Tau1CosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hElectron2Tau1_Electron2MetDeltaPhi[NpdfCounter]        = fs->make<TH1F>(("Electron2Tau1_Electron2MetDeltaPhi_"+j.str()).c_str(),        ("Electron2Tau1_Electron2MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hElectron2Tau1_Tau1MetDeltaPhi[NpdfCounter]             = fs->make<TH1F>(("Electron2Tau1_Tau1MetDeltaPhi_"+j.str()).c_str(),             ("Electron2Tau1_Tau1MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hElectron2MetDeltaPhiVsElectron2Tau1CosDphi[NpdfCounter] = fs->make<TH2F>(("Electron2MetDeltaPhiVsElectron2Tau1CosDphi_"+j.str()).c_str(), ("Electron2MetDeltaPhiVsElectron2Tau1CosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
	_hElectron2PtVsTau2Pt[NpdfCounter]          = fs->make<TH2F>(("Electron2PtVsTau2Pt_"+j.str()).c_str(),          ("Electron2PtVsTau2Pt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hElectron2Tau2DeltaR[NpdfCounter]          = fs->make<TH1F>(("Electron2Tau2DeltaR_"+j.str()).c_str(),          ("Electron2Tau2DeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hElectron2Tau2DeltaPtDivSumPt[NpdfCounter] = fs->make<TH1F>(("Electron2Tau2DeltaPtDivSumPt_"+j.str()).c_str(), ("Electron2Tau2DeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hElectron2Tau2DeltaPt[NpdfCounter]         = fs->make<TH1F>(("Electron2Tau2DeltaPt_"+j.str()).c_str(), ("Electron2Tau2DeltaPt_"+j.str()).c_str(), 100, 0, 1000);
	_hElectron2Tau2_Electron2MetMt[NpdfCounter]              = fs->make<TH1F>(("Electron2Tau2_Electron2MetMt_"+j.str()).c_str(),              ("Electron2Tau2_Electron2MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hElectron2Tau2_Tau2MetMt[NpdfCounter]                   = fs->make<TH1F>(("Electron2Tau2_Tau2MetMt_"+j.str()).c_str(),                   ("Electron2Tau2_Tau2MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hElectron2Tau2OSLS[NpdfCounter]            = fs->make<TH1F>(("Electron2Tau2OSLS_"+j.str()).c_str(),            ("Electron2Tau2OSLS_"+j.str()).c_str(), 20, -10, 10);
	_hElectron2Tau2CosDphi[NpdfCounter]         = fs->make<TH1F>(("Electron2Tau2CosDphi_"+j.str()).c_str(),         ("Electron2Tau2CosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hElectron2Tau2_Electron2MetDeltaPhi[NpdfCounter]        = fs->make<TH1F>(("Electron2Tau2_Electron2MetDeltaPhi_"+j.str()).c_str(),        ("Electron2Tau2_Electron2MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hElectron2Tau2_Tau2MetDeltaPhi[NpdfCounter]             = fs->make<TH1F>(("Electron2Tau2_Tau2MetDeltaPhi_"+j.str()).c_str(),             ("Electron2Tau2_Tau2MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hElectron2MetDeltaPhiVsElectron2Tau2CosDphi[NpdfCounter] = fs->make<TH2F>(("Electron2MetDeltaPhiVsElectron2Tau2CosDphi_"+j.str()).c_str(), ("Electron2MetDeltaPhiVsElectron2Tau2CosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
/*
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
*/
	_hTau1PtVsTau2Pt[NpdfCounter]          = fs->make<TH2F>(("Tau1PtVsTau2Pt_"+j.str()).c_str(),          ("Tau1PtVsTau2Pt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hTau1Tau2DeltaR[NpdfCounter]          = fs->make<TH1F>(("Tau1Tau2DeltaR_"+j.str()).c_str(),          ("Tau1Tau2DeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hTau1Tau2DeltaPtDivSumPt[NpdfCounter] = fs->make<TH1F>(("Tau1Tau2DeltaPtDivSumPt_"+j.str()).c_str(), ("Tau1Tau2DeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hTau1Tau2DeltaPt[NpdfCounter]         = fs->make<TH1F>(("Tau1Tau2DeltaPt_"+j.str()).c_str(), ("Tau1Tau2DeltaPt_"+j.str()).c_str(), 100, 0, 1000);
	_hDiTau_Tau1MetMt[NpdfCounter]               = fs->make<TH1F>(("DiTau_Tau1MetMt_"+j.str()).c_str(),               ("DiTau_Tau1MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hDiTau_Tau2MetMt[NpdfCounter]               = fs->make<TH1F>(("DiTau_Tau2MetMt_"+j.str()).c_str(),               ("DiTau_Tau2MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hTau1Tau2OSLS[NpdfCounter]            = fs->make<TH1F>(("Tau1Tau2OSLS_"+j.str()).c_str(),            ("Tau1Tau2OSLS_"+j.str()).c_str(), 20, -10, 10);
	_hTau1Tau2CosDphi[NpdfCounter]         = fs->make<TH1F>(("Tau1Tau2CosDphi_"+j.str()).c_str(),         ("Tau1Tau2CosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hDiTau_Tau1MetDeltaPhi[NpdfCounter]         = fs->make<TH1F>(("DiTau_Tau1MetDeltaPhi_"+j.str()).c_str(),         ("DiTau_Tau1MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hDiTau_Tau2MetDeltaPhi[NpdfCounter]         = fs->make<TH1F>(("DiTau_Tau2MetDeltaPhi_"+j.str()).c_str(),         ("DiTau_Tau2MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hTau1MetDeltaPhiVsTau1Tau2CosDphi[NpdfCounter] = fs->make<TH2F>(("Tau1MetDeltaPhiVsTau1Tau2CosDphi_"+j.str()).c_str(), ("Tau1MetDeltaPhiVsTau1Tau2CosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
	_hMuon1PtVsMuon2Pt[NpdfCounter]          = fs->make<TH2F>(("Muon1PtVsMuon2Pt_"+j.str()).c_str(),          ("Muon1PtVsMuon2Pt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hMuon1Muon2DeltaR[NpdfCounter]          = fs->make<TH1F>(("Muon1Muon2DeltaR_"+j.str()).c_str(),          ("Muon1Muon2DeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hMuon1Muon2DeltaPtDivSumPt[NpdfCounter] = fs->make<TH1F>(("Muon1Muon2DeltaPtDivSumPt_"+j.str()).c_str(), ("Muon1Muon2DeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hMuon1Muon2DeltaPt[NpdfCounter]         = fs->make<TH1F>(("Muon1Muon2DeltaPt_"+j.str()).c_str(), ("Muon1Muon2DeltaPt_"+j.str()).c_str(), 100, 0, 1000);
	_hDiMuon_Muon1MetMt[NpdfCounter]                = fs->make<TH1F>(("DiMuon_Muon1MetMt_"+j.str()).c_str(),                ("DiMuon_Muon1MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hDiMuon_Muon2MetMt[NpdfCounter]                = fs->make<TH1F>(("DiMuon_Muon2MetMt_"+j.str()).c_str(),                ("DiMuon_Muon2MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hMuon1Muon2OSLS[NpdfCounter]            = fs->make<TH1F>(("Muon1Muon2OSLS_"+j.str()).c_str(),            ("Muon1Muon2OSLS_"+j.str()).c_str(), 20, -10, 10);
	_hMuon1Muon2CosDphi[NpdfCounter]         = fs->make<TH1F>(("Muon1Muon2CosDphi_"+j.str()).c_str(),         ("Muon1Muon2CosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hDiMuon_Muon1MetDeltaPhi[NpdfCounter]          = fs->make<TH1F>(("DiMuon_Muon1MetDeltaPhi_"+j.str()).c_str(),          ("DiMuon_Muon1MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hDiMuon_Muon2MetDeltaPhi[NpdfCounter]          = fs->make<TH1F>(("DiMuon_Muon2MetDeltaPhi_"+j.str()).c_str(),          ("DiMuon_Muon2MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hMuon1MetDeltaPhiVsMuon1Muon2CosDphi[NpdfCounter] = fs->make<TH2F>(("Muon1MetDeltaPhiVsMuon1Muon2CosDphi_"+j.str()).c_str(), ("Muon1MetDeltaPhiVsMuon1Muon2CosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
	_hElectron1PtVsElectron2Pt[NpdfCounter]          = fs->make<TH2F>(("Electron1PtVsElectron2Pt_"+j.str()).c_str(), ("Electron1PtVsElectron2Pt_"+j.str()).c_str(), 100, 0, 500, 100, 0, 500);
	_hElectron1Electron2DeltaR[NpdfCounter]          = fs->make<TH1F>(("Electron1Electron2DeltaR_"+j.str()).c_str(), ("Electron1Electron2DeltaR_"+j.str()).c_str(), 100, 0, 5.);
	_hElectron1Electron2DeltaPtDivSumPt[NpdfCounter] = fs->make<TH1F>(("Electron1Electron2DeltaPtDivSumPt_"+j.str()).c_str(), ("Electron1Electron2DeltaPtDivSumPt_"+j.str()).c_str(), 100, -5, 5.);
	_hElectron1Electron2DeltaPt[NpdfCounter]         = fs->make<TH1F>(("Electron1Electron2DeltaPt_"+j.str()).c_str(), ("Electron1Electron2DeltaPt_"+j.str()).c_str(), 100, 0, 1000);
	_hDiElectron_Electron1MetMt[NpdfCounter]            = fs->make<TH1F>(("DiElectron_Electron1MetMt_"+j.str()).c_str(),            ("DiElectron_Electron1MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hDiElectron_Electron2MetMt[NpdfCounter]            = fs->make<TH1F>(("DiElectron_Electron2MetMt_"+j.str()).c_str(),            ("DiElectron_Electron2MetMt_"+j.str()).c_str(), 100, 0, 500);
	_hElectron1Electron2OSLS[NpdfCounter]    = fs->make<TH1F>(("Electron1Electron2OSLS_"+j.str()).c_str(),    ("Electron1Electron2OSLS_"+j.str()).c_str(), 20, -10, 10);
	_hElectron1Electron2CosDphi[NpdfCounter] = fs->make<TH1F>(("Electron1Electron2CosDphi_"+j.str()).c_str(), ("Electron1Electron2CosDphi_"+j.str()).c_str(), 220, -1.1, 1.1);
	_hDiElectron_Electron1MetDeltaPhi[NpdfCounter]      = fs->make<TH1F>(("DiElectron_Electron1MetDeltaPhi_"+j.str()).c_str(),      ("DiElectron_Electron1MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hDiElectron_Electron2MetDeltaPhi[NpdfCounter]      = fs->make<TH1F>(("DiElectron_Electron2MetDeltaPhi_"+j.str()).c_str(),      ("DiElectron_Electron2MetDeltaPhi_"+j.str()).c_str(), 72, 0, +TMath::Pi());
	_hElectron1MetDeltaPhiVsElectron1Electron2CosDphi[NpdfCounter] = fs->make<TH2F>(("Electron1MetDeltaPhiVsElectron1Electron2CosDphi_"+j.str()).c_str(), ("Electron1MetDeltaPhiVsElectron1Electron2CosDphi_"+j.str()).c_str(), 72, 0, +TMath::Pi(), 220, -1.1, 1.1);
      /*
	_hTauJetSumPtIso_SeedOS[NpdfCounter] = fs->make<TH1F>(("TauJetSumPtIso_SeedOS_"+j.str()).c_str(), ("TauJetSumPtIso_SeedOS_"+j.str()).c_str(), 100, 0, 50);
	_hTauJetSumPtIso_JetOS[NpdfCounter]  = fs->make<TH1F>(("TauJetSumPtIso_JetOS_"+j.str()).c_str(),  ("TauJetSumPtIso_JetOS_"+j.str()).c_str(), 100, 0, 50);
	_hTauJetSumPtIso_SeedLS[NpdfCounter] = fs->make<TH1F>(("TauJetSumPtIso_SeedLS_"+j.str()).c_str(), ("TauJetSumPtIso_SeedLS_"+j.str()).c_str(), 100, 0, 50);
	_hTauJetSumPtIso_JetLS[NpdfCounter]  = fs->make<TH1F>(("TauJetSumPtIso_JetLS_"+j.str()).c_str(),  ("TauJetSumPtIso_JetLS_"+j.str()).c_str(), 100, 0, 50);
      */
      _hDiTauNotReconstructableMass[NpdfCounter] = fs->make<TH1F>(("DiTauNotReconstructableMass_"+j.str()).c_str(), ("DiTauNotReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hDiTauReconstructableMass[NpdfCounter]    = fs->make<TH1F>(("DiTauReconstructableMass_"+j.str()).c_str(),    ("DiTauReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hDiTauNotReconstructableMassOS[NpdfCounter] = fs->make<TH1F>(("DiTauNotReconstructableMassOS_"+j.str()).c_str(), ("DiTauNotReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hDiTauReconstructableMassOS[NpdfCounter]    = fs->make<TH1F>(("DiTauReconstructableMassOS_"+j.str()).c_str(),    ("DiTauReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hDiTauNotReconstructableMassLS[NpdfCounter] = fs->make<TH1F>(("DiTauNotReconstructableMassLS_"+j.str()).c_str(), ("DiTauNotReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hDiTauReconstructableMassLS[NpdfCounter]    = fs->make<TH1F>(("DiTauReconstructableMassLS_"+j.str()).c_str(),    ("DiTauReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hDiMuonNotReconstructableMass[NpdfCounter] = fs->make<TH1F>(("DiMuonNotReconstructableMass_"+j.str()).c_str(), ("DiMuonNotReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hDiMuonReconstructableMass[NpdfCounter]    = fs->make<TH1F>(("DiMuonReconstructableMass_"+j.str()).c_str(),    ("DiMuonReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hDiMuonNotReconstructableMassOS[NpdfCounter] = fs->make<TH1F>(("DiMuonNotReconstructableMassOS_"+j.str()).c_str(), ("DiMuonNotReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hDiMuonReconstructableMassOS[NpdfCounter]    = fs->make<TH1F>(("DiMuonReconstructableMassOS_"+j.str()).c_str(),    ("DiMuonReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hDiMuonNotReconstructableMassLS[NpdfCounter] = fs->make<TH1F>(("DiMuonNotReconstructableMassLS_"+j.str()).c_str(), ("DiMuonNotReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hDiMuonReconstructableMassLS[NpdfCounter]    = fs->make<TH1F>(("DiMuonReconstructableMassLS_"+j.str()).c_str(),    ("DiMuonReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hDiElectronNotReconstructableMass[NpdfCounter] = fs->make<TH1F>(("DiElectronNotReconstructableMass_"+j.str()).c_str(), ("DiElectronNotReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hDiElectronReconstructableMass[NpdfCounter]    = fs->make<TH1F>(("DiElectronReconstructableMass_"+j.str()).c_str(),    ("DiElectronReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hDiElectronNotReconstructableMassOS[NpdfCounter] = fs->make<TH1F>(("DiElectronNotReconstructableMassOS_"+j.str()).c_str(), ("DiElectronNotReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hDiElectronReconstructableMassOS[NpdfCounter]    = fs->make<TH1F>(("DiElectronReconstructableMassOS_"+j.str()).c_str(),    ("DiElectronReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hDiElectronNotReconstructableMassLS[NpdfCounter] = fs->make<TH1F>(("DiElectronNotReconstructableMassLS_"+j.str()).c_str(), ("DiElectronNotReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hDiElectronReconstructableMassLS[NpdfCounter]    = fs->make<TH1F>(("DiElectronReconstructableMassLS_"+j.str()).c_str(),    ("DiElectronReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon1Tau1NotReconstructableMass[NpdfCounter] = fs->make<TH1F>(("Muon1Tau1NotReconstructableMass_"+j.str()).c_str(), ("Muon1Tau1NotReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon1Tau1ReconstructableMass[NpdfCounter]    = fs->make<TH1F>(("Muon1Tau1ReconstructableMass_"+j.str()).c_str(),    ("Muon1Tau1ReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon1Tau1NotReconstructableMassOS[NpdfCounter] = fs->make<TH1F>(("Muon1Tau1NotReconstructableMassOS_"+j.str()).c_str(), ("Muon1Tau1NotReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon1Tau1ReconstructableMassOS[NpdfCounter]    = fs->make<TH1F>(("Muon1Tau1ReconstructableMassOS_"+j.str()).c_str(),    ("Muon1Tau1ReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon1Tau1NotReconstructableMassLS[NpdfCounter] = fs->make<TH1F>(("Muon1Tau1NotReconstructableMassLS_"+j.str()).c_str(), ("Muon1Tau1NotReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon1Tau1ReconstructableMassLS[NpdfCounter]    = fs->make<TH1F>(("Muon1Tau1ReconstructableMassLS_"+j.str()).c_str(),    ("Muon1Tau1ReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon1Tau2NotReconstructableMass[NpdfCounter] = fs->make<TH1F>(("Muon1Tau2NotReconstructableMass_"+j.str()).c_str(), ("Muon1Tau2NotReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon1Tau2ReconstructableMass[NpdfCounter]    = fs->make<TH1F>(("Muon1Tau2ReconstructableMass_"+j.str()).c_str(),    ("Muon1Tau2ReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon1Tau2NotReconstructableMassOS[NpdfCounter] = fs->make<TH1F>(("Muon1Tau2NotReconstructableMassOS_"+j.str()).c_str(), ("Muon1Tau2NotReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon1Tau2ReconstructableMassOS[NpdfCounter]    = fs->make<TH1F>(("Muon1Tau2ReconstructableMassOS_"+j.str()).c_str(),    ("Muon1Tau2ReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon1Tau2NotReconstructableMassLS[NpdfCounter] = fs->make<TH1F>(("Muon1Tau2NotReconstructableMassLS_"+j.str()).c_str(), ("Muon1Tau2NotReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon1Tau2ReconstructableMassLS[NpdfCounter]    = fs->make<TH1F>(("Muon1Tau2ReconstructableMassLS_"+j.str()).c_str(),    ("Muon1Tau2ReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon2Tau1NotReconstructableMass[NpdfCounter] = fs->make<TH1F>(("Muon2Tau1NotReconstructableMass_"+j.str()).c_str(), ("Muon2Tau1NotReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon2Tau1ReconstructableMass[NpdfCounter]    = fs->make<TH1F>(("Muon2Tau1ReconstructableMass_"+j.str()).c_str(),    ("Muon2Tau1ReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon2Tau1NotReconstructableMassOS[NpdfCounter] = fs->make<TH1F>(("Muon2Tau1NotReconstructableMassOS_"+j.str()).c_str(), ("Muon2Tau1NotReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon2Tau1ReconstructableMassOS[NpdfCounter]    = fs->make<TH1F>(("Muon2Tau1ReconstructableMassOS_"+j.str()).c_str(),    ("Muon2Tau1ReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon2Tau1NotReconstructableMassLS[NpdfCounter] = fs->make<TH1F>(("Muon2Tau1NotReconstructableMassLS_"+j.str()).c_str(), ("Muon2Tau1NotReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon2Tau1ReconstructableMassLS[NpdfCounter]    = fs->make<TH1F>(("Muon2Tau1ReconstructableMassLS_"+j.str()).c_str(),    ("Muon2Tau1ReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon2Tau2NotReconstructableMass[NpdfCounter] = fs->make<TH1F>(("Muon2Tau2NotReconstructableMass_"+j.str()).c_str(), ("Muon2Tau2NotReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon2Tau2ReconstructableMass[NpdfCounter]    = fs->make<TH1F>(("Muon2Tau2ReconstructableMass_"+j.str()).c_str(),    ("Muon2Tau2ReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon2Tau2NotReconstructableMassOS[NpdfCounter] = fs->make<TH1F>(("Muon2Tau2NotReconstructableMassOS_"+j.str()).c_str(), ("Muon2Tau2NotReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon2Tau2ReconstructableMassOS[NpdfCounter]    = fs->make<TH1F>(("Muon2Tau2ReconstructableMassOS_"+j.str()).c_str(),    ("Muon2Tau2ReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon2Tau2NotReconstructableMassLS[NpdfCounter] = fs->make<TH1F>(("Muon2Tau2NotReconstructableMassLS_"+j.str()).c_str(), ("Muon2Tau2NotReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hMuon2Tau2ReconstructableMassLS[NpdfCounter]    = fs->make<TH1F>(("Muon2Tau2ReconstructableMassLS_"+j.str()).c_str(),    ("Muon2Tau2ReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron1Tau1NotReconstructableMass[NpdfCounter] = fs->make<TH1F>(("Electron1Tau1NotReconstructableMass_"+j.str()).c_str(), ("Electron1Tau1NotReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron1Tau1ReconstructableMass[NpdfCounter]    = fs->make<TH1F>(("Electron1Tau1ReconstructableMass_"+j.str()).c_str(),    ("Electron1Tau1ReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron1Tau1NotReconstructableMassOS[NpdfCounter] = fs->make<TH1F>(("Electron1Tau1NotReconstructableMassOS_"+j.str()).c_str(), ("Electron1Tau1NotReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron1Tau1ReconstructableMassOS[NpdfCounter]    = fs->make<TH1F>(("Electron1Tau1ReconstructableMassOS_"+j.str()).c_str(),    ("Electron1Tau1ReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron1Tau1NotReconstructableMassLS[NpdfCounter] = fs->make<TH1F>(("Electron1Tau1NotReconstructableMassLS_"+j.str()).c_str(), ("Electron1Tau1NotReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron1Tau1ReconstructableMassLS[NpdfCounter]    = fs->make<TH1F>(("Electron1Tau1ReconstructableMassLS_"+j.str()).c_str(),    ("Electron1Tau1ReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron1Tau2NotReconstructableMass[NpdfCounter] = fs->make<TH1F>(("Electron1Tau2NotReconstructableMass_"+j.str()).c_str(), ("Electron1Tau2NotReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron1Tau2ReconstructableMass[NpdfCounter]    = fs->make<TH1F>(("Electron1Tau2ReconstructableMass_"+j.str()).c_str(),    ("Electron1Tau2ReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron1Tau2NotReconstructableMassOS[NpdfCounter] = fs->make<TH1F>(("Electron1Tau2NotReconstructableMassOS_"+j.str()).c_str(), ("Electron1Tau2NotReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron1Tau2ReconstructableMassOS[NpdfCounter]    = fs->make<TH1F>(("Electron1Tau2ReconstructableMassOS_"+j.str()).c_str(),    ("Electron1Tau2ReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron1Tau2NotReconstructableMassLS[NpdfCounter] = fs->make<TH1F>(("Electron1Tau2NotReconstructableMassLS_"+j.str()).c_str(), ("Electron1Tau2NotReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron1Tau2ReconstructableMassLS[NpdfCounter]    = fs->make<TH1F>(("Electron1Tau2ReconstructableMassLS_"+j.str()).c_str(),    ("Electron1Tau2ReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron2Tau1NotReconstructableMass[NpdfCounter] = fs->make<TH1F>(("Electron2Tau1NotReconstructableMass_"+j.str()).c_str(), ("Electron2Tau1NotReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron2Tau1ReconstructableMass[NpdfCounter]    = fs->make<TH1F>(("Electron2Tau1ReconstructableMass_"+j.str()).c_str(),    ("Electron2Tau1ReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron2Tau1NotReconstructableMassOS[NpdfCounter] = fs->make<TH1F>(("Electron2Tau1NotReconstructableMassOS_"+j.str()).c_str(), ("Electron2Tau1NotReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron2Tau1ReconstructableMassOS[NpdfCounter]    = fs->make<TH1F>(("Electron2Tau1ReconstructableMassOS_"+j.str()).c_str(),    ("Electron2Tau1ReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron2Tau1NotReconstructableMassLS[NpdfCounter] = fs->make<TH1F>(("Electron2Tau1NotReconstructableMassLS_"+j.str()).c_str(), ("Electron2Tau1NotReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron2Tau1ReconstructableMassLS[NpdfCounter]    = fs->make<TH1F>(("Electron2Tau1ReconstructableMassLS_"+j.str()).c_str(),    ("Electron2Tau1ReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron2Tau2NotReconstructableMass[NpdfCounter] = fs->make<TH1F>(("Electron2Tau2NotReconstructableMass_"+j.str()).c_str(), ("Electron2Tau2NotReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron2Tau2ReconstructableMass[NpdfCounter]    = fs->make<TH1F>(("Electron2Tau2ReconstructableMass_"+j.str()).c_str(),    ("Electron2Tau2ReconstructableMass_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron2Tau2NotReconstructableMassOS[NpdfCounter] = fs->make<TH1F>(("Electron2Tau2NotReconstructableMassOS_"+j.str()).c_str(), ("Electron2Tau2NotReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron2Tau2ReconstructableMassOS[NpdfCounter]    = fs->make<TH1F>(("Electron2Tau2ReconstructableMassOS_"+j.str()).c_str(),    ("Electron2Tau2ReconstructableMassOS_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron2Tau2NotReconstructableMassLS[NpdfCounter] = fs->make<TH1F>(("Electron2Tau2NotReconstructableMassLS_"+j.str()).c_str(), ("Electron2Tau2NotReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hElectron2Tau2ReconstructableMassLS[NpdfCounter]    = fs->make<TH1F>(("Electron2Tau2ReconstructableMassLS_"+j.str()).c_str(),    ("Electron2Tau2ReconstructableMassLS_"+j.str()).c_str(), 150, 0, 1500);
      _hDiTauPZeta[NpdfCounter]                  = fs->make<TH1F>(("DiTauPZeta_"+j.str()).c_str(),                  ("DiTauPZeta_"+j.str()).c_str(), 200, -100, 100);
      _hDiTauPZetaVis[NpdfCounter]               = fs->make<TH1F>(("DiTauPZetaVis_"+j.str()).c_str(),               ("DiTauPZetaVis_"+j.str()).c_str(), 100, 0, 100);
      _hDiTauZeta2D[NpdfCounter]                 = fs->make<TH2F>(("DiTauZeta2D_"+j.str()).c_str(),                 ("DiTauZeta2D_"+j.str()).c_str(), 100, 0, 100, 200, -100, 100);
      _hDiTauZeta1D[NpdfCounter]                 = fs->make<TH1F>(("DiTauZeta1D_"+j.str()).c_str(),                 ("DiTauZeta1D_"+j.str()).c_str(), 50, -100, 100);
      _hDiMuonPZeta[NpdfCounter]                  = fs->make<TH1F>(("DiMuonPZeta_"+j.str()).c_str(),                  ("DiMuonPZeta_"+j.str()).c_str(), 200, -100, 100);
      _hDiMuonPZetaVis[NpdfCounter]               = fs->make<TH1F>(("DiMuonPZetaVis_"+j.str()).c_str(),               ("DiMuonPZetaVis_"+j.str()).c_str(), 100, 0, 100);
      _hDiMuonZeta2D[NpdfCounter]                 = fs->make<TH2F>(("DiMuonZeta2D_"+j.str()).c_str(),                 ("DiMuonZeta2D_"+j.str()).c_str(), 100, 0, 100, 200, -100, 100);
      _hDiMuonZeta1D[NpdfCounter]                 = fs->make<TH1F>(("DiMuonZeta1D_"+j.str()).c_str(),                 ("DiMuonZeta1D_"+j.str()).c_str(), 50, -100, 100);
      _hDiElectronPZeta[NpdfCounter]                  = fs->make<TH1F>(("DiElectronPZeta_"+j.str()).c_str(),                  ("DiElectronPZeta_"+j.str()).c_str(), 200, -100, 100);
      _hDiElectronPZetaVis[NpdfCounter]               = fs->make<TH1F>(("DiElectronPZetaVis_"+j.str()).c_str(),               ("DiElectronPZetaVis_"+j.str()).c_str(), 100, 0, 100);
      _hDiElectronZeta2D[NpdfCounter]                 = fs->make<TH2F>(("DiElectronZeta2D_"+j.str()).c_str(),                 ("DiElectronZeta2D_"+j.str()).c_str(), 100, 0, 100, 200, -100, 100);
      _hDiElectronZeta1D[NpdfCounter]                 = fs->make<TH1F>(("DiElectronZeta1D_"+j.str()).c_str(),                 ("DiElectronZeta1D_"+j.str()).c_str(), 50, -100, 100);
      _hMuon1Tau1PZeta[NpdfCounter]                  = fs->make<TH1F>(("Muon1Tau1PZeta_"+j.str()).c_str(),                  ("Muon1Tau1PZeta_"+j.str()).c_str(), 200, -100, 100);
      _hMuon1Tau1PZetaVis[NpdfCounter]               = fs->make<TH1F>(("Muon1Tau1PZetaVis_"+j.str()).c_str(),               ("Muon1Tau1PZetaVis_"+j.str()).c_str(), 100, 0, 100);
      _hMuon1Tau1Zeta2D[NpdfCounter]                 = fs->make<TH2F>(("Muon1Tau1Zeta2D_"+j.str()).c_str(),                 ("Muon1Tau1Zeta2D_"+j.str()).c_str(), 100, 0, 100, 200, -100, 100);
      _hMuon1Tau1Zeta1D[NpdfCounter]                 = fs->make<TH1F>(("Muon1Tau1Zeta1D_"+j.str()).c_str(),                 ("Muon1Tau1Zeta1D_"+j.str()).c_str(), 50, -100, 100);
      _hMuon1Tau2PZeta[NpdfCounter]                  = fs->make<TH1F>(("Muon1Tau2PZeta_"+j.str()).c_str(),                  ("Muon1Tau2PZeta_"+j.str()).c_str(), 200, -100, 100);
      _hMuon1Tau2PZetaVis[NpdfCounter]               = fs->make<TH1F>(("Muon1Tau2PZetaVis_"+j.str()).c_str(),               ("Muon1Tau2PZetaVis_"+j.str()).c_str(), 100, 0, 100);
      _hMuon1Tau2Zeta2D[NpdfCounter]                 = fs->make<TH2F>(("Muon1Tau2Zeta2D_"+j.str()).c_str(),                 ("Muon1Tau2Zeta2D_"+j.str()).c_str(), 100, 0, 100, 200, -100, 100);
      _hMuon1Tau2Zeta1D[NpdfCounter]                 = fs->make<TH1F>(("Muon1Tau2Zeta1D_"+j.str()).c_str(),                 ("Muon1Tau2Zeta1D_"+j.str()).c_str(), 50, -100, 100);
      _hMuon2Tau1PZeta[NpdfCounter]                  = fs->make<TH1F>(("Muon2Tau1PZeta_"+j.str()).c_str(),                  ("Muon2Tau1PZeta_"+j.str()).c_str(), 200, -100, 100);
      _hMuon2Tau1PZetaVis[NpdfCounter]               = fs->make<TH1F>(("Muon2Tau1PZetaVis_"+j.str()).c_str(),               ("Muon2Tau1PZetaVis_"+j.str()).c_str(), 100, 0, 100);
      _hMuon2Tau1Zeta2D[NpdfCounter]                 = fs->make<TH2F>(("Muon2Tau1Zeta2D_"+j.str()).c_str(),                 ("Muon2Tau1Zeta2D_"+j.str()).c_str(), 100, 0, 100, 200, -100, 100);
      _hMuon2Tau1Zeta1D[NpdfCounter]                 = fs->make<TH1F>(("Muon2Tau1Zeta1D_"+j.str()).c_str(),                 ("Muon2Tau1Zeta1D_"+j.str()).c_str(), 50, -100, 100);
      _hMuon2Tau2PZeta[NpdfCounter]                  = fs->make<TH1F>(("Muon2Tau2PZeta_"+j.str()).c_str(),                  ("Muon2Tau2PZeta_"+j.str()).c_str(), 200, -100, 100);
      _hMuon2Tau2PZetaVis[NpdfCounter]               = fs->make<TH1F>(("Muon2Tau2PZetaVis_"+j.str()).c_str(),               ("Muon2Tau2PZetaVis_"+j.str()).c_str(), 100, 0, 100);
      _hMuon2Tau2Zeta2D[NpdfCounter]                 = fs->make<TH2F>(("Muon2Tau2Zeta2D_"+j.str()).c_str(),                 ("Muon2Tau2Zeta2D_"+j.str()).c_str(), 100, 0, 100, 200, -100, 100);
      _hMuon2Tau2Zeta1D[NpdfCounter]                 = fs->make<TH1F>(("Muon2Tau2Zeta1D_"+j.str()).c_str(),                 ("Muon2Tau2Zeta1D_"+j.str()).c_str(), 50, -100, 100);
      _hElectron1Tau1PZeta[NpdfCounter]                  = fs->make<TH1F>(("Electron1Tau1PZeta_"+j.str()).c_str(),                  ("Electron1Tau1PZeta_"+j.str()).c_str(), 200, -100, 100);
      _hElectron1Tau1PZetaVis[NpdfCounter]               = fs->make<TH1F>(("Electron1Tau1PZetaVis_"+j.str()).c_str(),               ("Electron1Tau1PZetaVis_"+j.str()).c_str(), 100, 0, 100);
      _hElectron1Tau1Zeta2D[NpdfCounter]                 = fs->make<TH2F>(("Electron1Tau1Zeta2D_"+j.str()).c_str(),                 ("Electron1Tau1Zeta2D_"+j.str()).c_str(), 100, 0, 100, 200, -100, 100);
      _hElectron1Tau1Zeta1D[NpdfCounter]                 = fs->make<TH1F>(("Electron1Tau1Zeta1D_"+j.str()).c_str(),                 ("Electron1Tau1Zeta1D_"+j.str()).c_str(), 50, -100, 100);
      _hElectron1Tau2PZeta[NpdfCounter]                  = fs->make<TH1F>(("Electron1Tau2PZeta_"+j.str()).c_str(),                  ("Electron1Tau2PZeta_"+j.str()).c_str(), 200, -100, 100);
      _hElectron1Tau2PZetaVis[NpdfCounter]               = fs->make<TH1F>(("Electron1Tau2PZetaVis_"+j.str()).c_str(),               ("Electron1Tau2PZetaVis_"+j.str()).c_str(), 100, 0, 100);
      _hElectron1Tau2Zeta2D[NpdfCounter]                 = fs->make<TH2F>(("Electron1Tau2Zeta2D_"+j.str()).c_str(),                 ("Electron1Tau2Zeta2D_"+j.str()).c_str(), 100, 0, 100, 200, -100, 100);
      _hElectron1Tau2Zeta1D[NpdfCounter]                 = fs->make<TH1F>(("Electron1Tau2Zeta1D_"+j.str()).c_str(),                 ("Electron1Tau2Zeta1D_"+j.str()).c_str(), 50, -100, 100);
      _hElectron2Tau1PZeta[NpdfCounter]                  = fs->make<TH1F>(("Electron2Tau1PZeta_"+j.str()).c_str(),                  ("Electron2Tau1PZeta_"+j.str()).c_str(), 200, -100, 100);
      _hElectron2Tau1PZetaVis[NpdfCounter]               = fs->make<TH1F>(("Electron2Tau1PZetaVis_"+j.str()).c_str(),               ("Electron2Tau1PZetaVis_"+j.str()).c_str(), 100, 0, 100);
      _hElectron2Tau1Zeta2D[NpdfCounter]                 = fs->make<TH2F>(("Electron2Tau1Zeta2D_"+j.str()).c_str(),                 ("Electron2Tau1Zeta2D_"+j.str()).c_str(), 100, 0, 100, 200, -100, 100);
      _hElectron2Tau1Zeta1D[NpdfCounter]                 = fs->make<TH1F>(("Electron2Tau1Zeta1D_"+j.str()).c_str(),                 ("Electron2Tau1Zeta1D_"+j.str()).c_str(), 50, -100, 100);
      _hElectron2Tau2PZeta[NpdfCounter]                  = fs->make<TH1F>(("Electron2Tau2PZeta_"+j.str()).c_str(),                  ("Electron2Tau2PZeta_"+j.str()).c_str(), 200, -100, 100);
      _hElectron2Tau2PZetaVis[NpdfCounter]               = fs->make<TH1F>(("Electron2Tau2PZetaVis_"+j.str()).c_str(),               ("Electron2Tau2PZetaVis_"+j.str()).c_str(), 100, 0, 100);
      _hElectron2Tau2Zeta2D[NpdfCounter]                 = fs->make<TH2F>(("Electron2Tau2Zeta2D_"+j.str()).c_str(),                 ("Electron2Tau2Zeta2D_"+j.str()).c_str(), 100, 0, 100, 200, -100, 100);
      _hElectron2Tau2Zeta1D[NpdfCounter]                 = fs->make<TH1F>(("Electron2Tau2Zeta1D_"+j.str()).c_str(),                 ("Electron2Tau2Zeta1D_"+j.str()).c_str(), 50, -100, 100);
      _hMet[NpdfCounter]                    = fs->make<TH1F>(("Met_"+j.str()).c_str(),                    ("Met_"+j.str()).c_str(), 100, 0, 1000);
      _hMetResolution[NpdfCounter]          = fs->make<TH1F>(("MetResolution_"+j.str()).c_str(),                    ("MetResolution_"+j.str()).c_str(), 500, -5, 5);
      _hR1[NpdfCounter] 		    = fs->make<TH1F>(("R1_"+j.str()).c_str(), ("R1_"+j.str()).c_str(), 60, 0, 6);
      _hR2[NpdfCounter] = fs->make<TH1F>(("R2_"+j.str()).c_str(), ("R2_"+j.str()).c_str(), 60, 0, 6);
      _hDphi1MHT[NpdfCounter] = fs->make<TH1F>(("Dphi1_"+j.str()).c_str(), ("Dphi1_"+j.str()).c_str(), 72, -2.0 * TMath::Pi(), +2.0 * TMath::Pi());
      _hDphi2MHT[NpdfCounter] = fs->make<TH1F>(("Dphi2MHT_"+j.str()).c_str(), ("Dphi2MHT_"+j.str()).c_str(), 72, -2.0 * TMath::Pi(), +2.0 * TMath::Pi());
      _hDphi1[NpdfCounter] = fs->make<TH1F>(("Dphi1MHT_"+j.str()).c_str(), ("Dphi1MHT_"+j.str()).c_str(), 72, -2.0 * TMath::Pi(), +2.0 * TMath::Pi());
      _hDphi2[NpdfCounter] = fs->make<TH1F>(("Dphi2_"+j.str()).c_str(), ("Dphi2_"+j.str()).c_str(), 72, -2.0 * TMath::Pi(), +2.0 * TMath::Pi());
      _hDphi1VsDphi2[NpdfCounter] = fs->make<TH2F>(("Dphi1VsDphi2_"+j.str()).c_str(), ("Dphi1VsDphi2_"+j.str()).c_str(), 72, -2.0 * TMath::Pi(), +2.0 * TMath::Pi(), 72, -2.0 * TMath::Pi(), +2.0 * TMath::Pi());
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

//////////////////////////////////////////////////////////////////////////////
// The following three functions are for PU reweighting
// The PU histograms are NVertices_0 in our current framework 11/23/11
// getBin retrieves the correct bin for our N_vertex value
// getValue determines the number of events with that N_vertex
// getintegral determines the total number of events in the N_vertex histogram
// so that we can normalize to obtain the probability that each file has
// X vertices... Will 11/23/11
//////////////////////////////////////////////////////////////////////////////

int HiMassTauAnalysis::getBin(char * cstr, int nVertices) {
    TFile *file1 = new TFile (cstr);
    TH1 *hist = dynamic_cast<TH1*>(file1->Get("analyzeHiMassTau/NVertices_0"));
    if (!hist) {throw std::runtime_error("failed to extract histogram");}

    int result = hist->GetBin(nVertices+1);
    file1->Close();
    return result;
}

double HiMassTauAnalysis::getvalue(char * cstr, int bin) {
    TFile *file1 = new TFile (cstr);
    TH1 *hist = dynamic_cast<TH1*>(file1->Get("analyzeHiMassTau/NVertices_0"));
    if (!hist) {throw std::runtime_error("failed to extract histogram");}

    double result = hist->GetBinContent(bin);
    file1->Close();
    return result;
    
}

double HiMassTauAnalysis::getintegral(char * cstr, int bin) {
    TFile *file1 = new TFile (cstr);
    TH1 *hist = dynamic_cast<TH1*>(file1->Get("analyzeHiMassTau/NVertices_0"));
    if (!hist) {throw std::runtime_error("failed to extract histogram");}

    double result = hist->Integral();
    file1->Close();
    return result;
}

HiMassTauAnalysis::~HiMassTauAnalysis() { }

//define this as a plug-in
DEFINE_FWK_MODULE(HiMassTauAnalysis);


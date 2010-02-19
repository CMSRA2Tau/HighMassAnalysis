import FWCore.ParameterSet.Config as cms
import copy

# Standard pat sequences
from PhysicsTools.PatAlgos.patSequences_cff import * 

# --------------------Modifications for electrons--------------------

from RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi import *
#elecIdCutBasedRobust = eidCutBasedExt.copy()
#elecIdCutBasedRobust.src = cms.InputTag("gsfElectrons")
#elecIdCutBasedRobust.electronQuality = 'robust'
#elecIdCutBasedRobust.robustEleIDCuts = cms.PSet(
#   barrel = cms.vdouble(0.015, 0.012, 0.02, 0.0025), # [0.015, 0.0092, 0.020, 0.0025]
#   endcap = cms.vdouble(0.018, 0.025, 0.02, 0.0040)  # [0.018, 0.025, 0.020, 0.0040]
#)
elecIdCutBasedLoose = eidCutBasedExt.copy()
elecIdCutBasedLoose.electronQuality = 'loose'
elecIdCutBasedTight = eidCutBasedExt.copy()
elecIdCutBasedTight.electronQuality = 'tight'
electronIdCutBased = cms.Sequence( 
                                  elecIdCutBasedLoose
                                  *elecIdCutBasedTight )

from RecoEgamma.EgammaIsolationAlgos.eleTrackExtractorBlocks_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleEcalExtractorBlocks_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleHcalExtractorBlocks_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositTk_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositEcalFromHits_cff import *
from PhysicsTools.PatAlgos.recoLayer0.electronIsolation_cff import *
#from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositHcalFromTowers_cff import *
EleIsoTrackExtractorBlock.DR_Max = cms.double(1.)
eleIsoDepositTk.src = cms.InputTag("gsfElectrons")
eleIsoDepositTk.ExtractorPSet = cms.PSet(EleIsoTrackExtractorBlock)
EleIsoEcalFromHitsExtractorBlock.extRadius = cms.double(1.)
eleIsoDepositEcalFromHits.src = cms.InputTag("gsfElectrons")
eleIsoDepositEcalFromHits.ExtractorPSet = cms.PSet(EleIsoEcalFromHitsExtractorBlock)
EleIsoHcalFromTowersExtractorBlock.extRadius = cms.double(1.)
eleIsoDepositHcalFromTowers.src = cms.InputTag("gsfElectrons")
eleIsoDepositHcalFromTowers.ExtractorPSet = cms.PSet(EleIsoHcalFromTowersExtractorBlock)
electronIsoDeposits = cms.Sequence( eleIsoDepositTk
                                   *eleIsoDepositEcalFromHits 
                                   *eleIsoDepositHcalFromTowers )
recoElectronIsolation = cms.Sequence( electronIsoDeposits )

from PhysicsTools.PatAlgos.recoLayer0.electronId_cff import *
from PhysicsTools.PatAlgos.recoLayer0.electronIsolation_cff import *
#from PhysicsTools.PatAlgos.recoLayer0.aodReco_cff import *
#from PhysicsTools.PatAlgos.triggerLayer0.trigMatchSequences_cff import *
from PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi import *
from PhysicsTools.PatAlgos.cleaningLayer1.electronCleaner_cfi import *
#patTrigMatchElectron = cms.Sequence( electronTrigMatchHLT1Electron )
#patTrigMatch._seq = patTrigMatch._seq * patHLT1Electron * patTrigMatchElectron
allLayer1Electrons.isolation.tracker.src = cms.InputTag("eleIsoDepositTk")
allLayer1Electrons.isolation.tracker.deltaR = cms.double(0.6)
allLayer1Electrons.isolation.tracker.vetos = cms.vstring(
    '0.015',         # inner radius veto cone (was 0.015)
    'Threshold(1.0)' # threshold on individual track pt
)
allLayer1Electrons.isolation.tracker.skipDefaultVeto = cms.bool(True)
allLayer1Electrons.isolation.ecal.src = cms.InputTag("eleIsoDepositEcalFromHits")
allLayer1Electrons.isolation.ecal.deltaR = cms.double(0.6)
allLayer1Electrons.isolation.ecal.vetos = cms.vstring(
    'EcalBarrel:0.045', 
    'EcalBarrel:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)',
    'EcalEndcaps:0.07',                         #0.07
    'EcalEndcaps:RectangularEtaPhiVeto(-0.05,0.05,-0.5,0.5)',
    'EcalBarrel:ThresholdFromTransverse(0.12)', # default is 3*sigmaBarrel; sigmaBarrel = 0.04
    'EcalEndcaps:ThresholdFromTransverse(0.45)' # default is 3*sigmaEndCap; sigmaEndCap = 0.15
)
allLayer1Electrons.isolation.ecal.skipDefaultVeto = cms.bool(True)
allLayer1Electrons.isolation.hcal.src = cms.InputTag("eleIsoDepositHcalFromTowers")
allLayer1Electrons.isolation.hcal.deltaR = cms.double(0.6)
allLayer1Electrons.isoDeposits = cms.PSet(
   tracker         = allLayer1Electrons.isolation.tracker.src,
   ecal            = allLayer1Electrons.isolation.ecal.src,
   hcal            = allLayer1Electrons.isolation.hcal.src,
)
allLayer1Electrons.addElectronID = cms.bool(True)
allLayer1Electrons.electronIDSources = cms.PSet(
   #robust = cms.InputTag("elecIdCutBasedRobust"),
   loose  = cms.InputTag("elecIdCutBasedLoose"),
   tight  = cms.InputTag("elecIdCutBasedTight")        
)
#allLayer1Electrons.addTrigMatch = cms.bool(True)
#allLayer1Electrons.trigPrimMatch = cms.VInputTag(
#    cms.InputTag("electronTrigMatchHLT1Electron")
#)
allLayer1Electrons.addGenMatch = cms.bool(True)
allLayer1Electrons.genParticleMatch = cms.InputTag("electronMatch")
cleanLayer1Electrons.checkOverlaps = cms.PSet()

# --------------------Modifications for muons--------------------

#from PhysicsTools.PatAlgos.recoLayer0.muonIsolation_cff import *
#from PhysicsTools.PatAlgos.triggerLayer0.trigMatchSequences_cff import *
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import *
from PhysicsTools.PatAlgos.cleaningLayer1.muonCleaner_cfi import *
#patTrigMatch._seq = patTrigMatch._seq * patTrigMatchHLT1MuonIso
allLayer1Muons.isolation.tracker.deltaR = cms.double(0.6)
allLayer1Muons.isolation.ecal.deltaR = cms.double(0.6)
allLayer1Muons.isolation.hcal.deltaR = cms.double(0.6)
allLayer1Muons.isolation.user.deltaR = cms.double(0.6)
allLayer1Muons.isoDeposits = cms.PSet(
   tracker         = allLayer1Muons.isolation.tracker.src,
   ecal            = allLayer1Muons.isolation.ecal.src,
   hcal            = allLayer1Muons.isolation.hcal.src,
)
#allLayer1Muons.addTrigMatch = cms.bool(True)
#allLayer1Muons.trigPrimMatch = cms.VInputTag(
#    cms.InputTag("muonTrigMatchHLT1MuonNonIso"),
#    cms.InputTag("muonTrigMatchHLT1MuonIso")
#)
allLayer1Muons.addGenMatch = cms.bool(True)

# --------------------Modifications for taus--------------------

from PhysicsTools.PatAlgos.producersLayer1.tauProducer_cfi import *
from PhysicsTools.PatAlgos.cleaningLayer1.tauCleaner_cfi import * 
allLayer1Taus.addDecayMode = cms.bool(True)
cleanLayer1Taus.preselection = cms.string('')
cleanLayer1Taus.checkOverlaps = cms.PSet()

# --------------------Modifications for jets--------------------

cleanLayer1Jets.checkOverlaps = cms.PSet()

# --------------------Modifications for met--------------------

#from PhysicsTools.PatAlgos.recoLayer0.jetMETCorrections_cff import *
#from PhysicsTools.PatAlgos.producersLayer1.metProducer_cff import *
#metJESCorIC5CaloJetMuons = corMetGlobalMuons.clone(uncorMETInputTag = cms.InputTag('metJESCorIC5CaloJet'))
## apply tau-jet specific corrections
#from JetMETCorrections.Type1MET.TauMetCorrections_cff import * 
#tauMetCorr.InputTausLabel = cms.string('fixedConePFTauProducer')
#tauMetCorr.InputCaloJetsLabel = metJESCorIC5CaloJet.inputUncorJetsLabel
#tauMetCorr.jetPTthreshold = metJESCorIC5CaloJet.jetPTthreshold
#tauMetCorr.jetEMfracLimit = metJESCorIC5CaloJet.jetEMfracLimit
#tauMetCorr.correctorLabel = metJESCorIC5CaloJet.corrector
#tauMetCorr.InputMETLabel = cms.string('metJESCorIC5CaloJetMuons')
#tauMetCorr.seedTrackPt = cms.double(5.0)
##patMETCorrections._seq = patMETCorrections._seq * MetTauCorrections
## change input for pat::MET production
#layer1METs.metSource  = cms.InputTag("tauMetCorr")
##layer1METs.addTrigMatch  = cms.bool(False)
#makeLayer1METs.replace( patMETCorrections, metJESCorIC5CaloJet * metJESCorIC5CaloJetMuons * MetTauCorrections )

# --------------------MC matching------------------------------

from PhysicsTools.PatAlgos.mcMatchLayer0.tauMatch_cfi import *
#from PhysicsTools.PatAlgos.triggerLayer0.patTrigMatcher_cfi import *
tauMatchFixedConeHighEff = tauMatch.copy()
tauMatchFixedConeHighEff.src = cms.InputTag("fixedConeHighEffPFTauProducer")
tauGenJetMatchFixedConeHighEff = tauGenJetMatch.copy()
tauGenJetMatchFixedConeHighEff.src = cms.InputTag("fixedConeHighEffPFTauProducer")
#tauTrigMatchHLT1TauFixedConeHighEff = tauTrigMatchHLT1Tau.copy()
#tauTrigMatchHLT1TauFixedConeHighEff.src = cms.InputTag("fixedConeHighEffPFTauProducer")
tauMatchFixedCone = tauMatch.copy()
tauMatchFixedCone.src = cms.InputTag("fixedConePFTauProducer")
tauGenJetMatchFixedCone = tauGenJetMatch.copy()
tauGenJetMatchFixedCone.src = cms.InputTag("fixedConePFTauProducer")
#tauTrigMatchHLT1TauFixedCone = tauTrigMatchHLT1Tau.copy()
#tauTrigMatchHLT1TauFixedCone.src = cms.InputTag("fixedConePFTauProducer")
tauMatchShrinkingCone = tauMatch.copy()
tauMatchShrinkingCone.src = cms.InputTag("shrinkingConePFTauProducer")
tauGenJetMatchShrinkingCone = tauGenJetMatch.copy()
tauGenJetMatchShrinkingCone.src = cms.InputTag("shrinkingConePFTauProducer")
#tauTrigMatchHLT1TauShrinkingCone = tauTrigMatchHLT1Tau.copy()
#tauTrigMatchHLT1TauShrinkingCone.src = cms.InputTag("shrinkingConePFTauProducer")
tauMatchShrinkingTightCone = tauMatch.copy()
tauMatchShrinkingTightCone.src = cms.InputTag("shrinkingTightConePFTauProducer")
tauGenJetMatchShrinkingTightCone = tauGenJetMatch.copy()
tauGenJetMatchShrinkingTightCone.src = cms.InputTag("shrinkingTightConePFTauProducer")
#tauTrigMatchHLT1TauShrinkingTightCone = tauTrigMatchHLT1Tau.copy()
#tauTrigMatchHLT1TauShrinkingTightCone.src = cms.InputTag("shrinkingTightConePFTauProducer")


import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.patEventContent_cff import *
from Configuration.EventContent.EventContent_cff import *

#--------------------------------------------------------------------------------
# per default, drop everything that is not specifically kept
#--------------------------------------------------------------------------------

patTupleEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *')
)

#--------------------------------------------------------------------------------
# keep PAT layer 1 objects
#--------------------------------------------------------------------------------

patTupleEventContent.outputCommands.extend(patEventContentNoCleaning)
patTupleEventContent.outputCommands.extend(patExtraAodEventContent)
patTupleEventContent.outputCommands.extend(patTriggerEventContent)

#--------------------------------------------------------------------------------
# keep collections of additional collections
#--------------------------------------------------------------------------------

patTupleEventContent.outputCommands.extend(
    [
    'keep *_TriggerResults_*_*',
    'keep *_hcalnoise_*_*',
    'keep recoPFCandidate*_particleFlow__*',
    'keep recoTracks_generalTracks_*_*',
    'keep recoTrackExtras_generalTracks_*_*',
    'keep recoMuons_muons_*_*',
    'keep *_muonMETValueMapProducer_*_*',
    'keep booledmValueMap_muid*_*_*',
    'keep *_selectedLayer1CaloTaus_*_*',
    'keep *_selectedLayer1FixedConeHighEffPFTaus_*_*',
    'keep *_selectedLayer1FixedConePFTaus_*_*',
    'keep *_selectedLayer1ShrinkingConeHighEffPFTaus_*_*',
    'keep *_selectedLayer1ShrinkingConePFTaus_*_*',
    'keep *_selectedLayer1HPSPFTaus_*_*',
    'keep patJets_patJetsAK5PF_*_*',
    'keep *_layer1PFMETs_*_*',
    'keep *_layer1METsPF_*_*',
    'keep *_patMETs_*_*',
    'keep *_patMETsPF_*_*',
    'keep *_genMetCalo__hiMassTau',
    'keep *_genMetTrue__hiMassTau',
    'keep *_generator_*_*',      
    'keep recoGsfElectronCores_*_*_*',
    'keep recoSuperClusters_*_*_*',
    'keep DcsStatuss_scalersRawToDigi_*_*',
    'keep *_hpsPFTauDiscriminationAgainstElectron_*_*',
    'keep *_hpsPFTauDiscriminationAgainstMuon_*_*',
    'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*',
    'keep *_hpsPFTauDiscriminationByLooseIsolation_*_*',
    'keep *_hpsPFTauDiscriminationByMediumIsolation_*_*',
    'keep *_hpsPFTauDiscriminationByTightIsolation_*_*',
    'keep *_ak5PFJetsLegacyHPSPiZeros_*_*',
    'keep *_ak5PFJetsLegacyTaNCPiZeros_*_*',
    'keep *_ak5PFJetsRecoTauPiZeros_*_*',
    'keep *_hpsPFTauProducer_*_*'
    ]
)

#--------------------------------------------------------------------------------
# keep collections of selected PAT layer 1 particles
#--------------------------------------------------------------------------------

patTupleEventContent.outputCommands.extend(
    [ 'keep *_selectedPatElectrons_*_*',
      'keep *_heepPatElectrons_*_*',
      'keep *_selectedPatMuons_*_*',
      'keep *_selectedPatJets_*_*' ]
)

#--------------------------------------------------------------------------------
# keep generator level tau-jets
#--------------------------------------------------------------------------------

patTupleEventContent.outputCommands.extend(
    [ 'keep *_tauGenJets*_*_*', ]
)

#--------------------------------------------------------------------------------
# 
#--------------------------------------------------------------------------------
patTupleEventContent.outputCommands.extend(HLTriggerRECO.outputCommands)

patTupleEventContent.outputCommands.extend(
    [ 'drop CaloTowers*_*_*_*',
      'drop patPhotons_*_*_*',
      'drop *_selectedLayer1Taus_*_*' ]
)

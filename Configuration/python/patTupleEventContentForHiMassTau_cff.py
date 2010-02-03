import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.patEventContent_cff import *

#--------------------------------------------------------------------------------
# per default, drop everything that is not specifically kept
#--------------------------------------------------------------------------------

patTupleEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *')
)

#--------------------------------------------------------------------------------
# keep PAT layer 1 objects
#--------------------------------------------------------------------------------

#patTupleEventContent.outputCommands.extend(patEventContent)
patTupleEventContent.outputCommands.extend(patEventContentNoLayer1Cleaning)
patTupleEventContent.outputCommands.extend(patExtraAodEventContent)
patTupleEventContent.outputCommands.extend(patTriggerEventContent)


#--------------------------------------------------------------------------------
# keep collections of additional collections
#--------------------------------------------------------------------------------

patTupleEventContent.outputCommands.extend(
    [
    'keep recoPFCandidate*_particleFlow__*',
    'keep recoTracks_generalTracks_*_*',
    'keep recoTrackExtras_generalTracks_*_*',
    'keep recoMuons_muons_*_*',
    'keep *_muonMETValueMapProducer_*_*',
    'keep *_selectedLayer1CaloTaus_*_*',
    'keep *_selectedLayer1FixedConeHighEffPFTaus_*_*',
    'keep *_selectedLayer1FixedConePFTaus_*_*',
    'keep *_selectedLayer1ShrinkingConeHighEffPFTaus_*_*',
    'keep *_selectedLayer1ShrinkingConePFTaus_*_*',
    'keep *_allLayer1JetsPF_*_*',
    'keep *_layer1PFMETs_*_*',
    'keep *_layer1METsPF_*_*',
    'keep *_genMetCalo__hiMassTau',
    'keep *_genMetTrue__hiMassTau',
#    'keep *_cteq65PdfWeights_*_*',
#    'keep *_MRST2006nnloPdfWeights_*_*',
#    'keep *_MRST2007lomodPdfWeights_*_*',
    'keep *_generator_*_*',
    ]
)

#--------------------------------------------------------------------------------
# keep collections of selected PAT layer 1 particles
#--------------------------------------------------------------------------------

patTupleEventContent.outputCommands.extend(
    [ 'keep *_selectedLayer1Electrons_*_*',
      'keep *_selectedLayer1Muons_*_*',
#      'keep *_selectedLayer1Taus_*_*',
      'keep *_selectedLayer1Jets_*_*' ]
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

patTupleEventContent.outputCommands.extend(
    [ 'drop CaloTowers*_*_*_*',
      'drop patPhotons_*_*_*',
      'drop *_selectedLayer1Taus_*_*' ]
)


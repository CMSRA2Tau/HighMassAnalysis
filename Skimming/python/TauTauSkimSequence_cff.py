import FWCore.ParameterSet.Config as cms

from HighMassAnalysis.Skimming.tauSelector_cfi import *

TauTauPairs = cms.EDProducer("DeltaRMinCandCombiner",
  decay = cms.string('selectedPFTaus@+ selectedPFTaus@-'),
  checkCharge = cms.bool(False),
  cut = cms.string( ''),
  name = cms.string('TauTauCandidates'),
  deltaRMin = cms.double(0.7)
)

selectedTauTauPairs = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('TauTauPairs'),
  minNumber = cms.uint32(1)                      
)

TauTauSkimSequence = cms.Sequence(
    selectedPFTaus
  * TauTauPairs
  * selectedTauTauPairs
)

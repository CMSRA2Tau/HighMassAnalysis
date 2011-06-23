import FWCore.ParameterSet.Config as cms

from HighMassAnalysis.Skimming.tauSelector_cfi import *

TauTauPairs = cms.EDProducer("DeltaRMinCandCombiner",
  decay = cms.string('selectedPFTaus@+ selectedPFTaus@-'),
  checkCharge = cms.bool(False),
  cut = cms.string( 'daughter(0).pt > 15 && daughter(1).pt > 15 && abs(daughter(0).eta) < 2.1 && abs(daughter(1).eta) < 2.1'),
  name = cms.string('TauTauCandidates'),
  deltaRMin = cms.double(0.3)
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

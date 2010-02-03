import FWCore.ParameterSet.Config as cms

from HighMassAnalysis.Skimming.tauSelector_cfi import *
from HighMassAnalysis.Skimming.muonSelector_cfi import *

muTauPairs = cms.EDProducer("DeltaRMinCandCombiner",
  decay = cms.string('selectedPFTaus@+ selectedMuons@-'),
  checkCharge = cms.bool(False),
  cut = cms.string( ''),
  name = cms.string('muTauCandidates'),
  deltaRMin = cms.double(0.7)
)

selectedMuTauPairs = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('muTauPairs'),
  minNumber = cms.uint32(1)                      
)

muTauSkimSequence = cms.Sequence(
  ( selectedPFTaus + selectedMuons )
  * muTauPairs
  * selectedMuTauPairs
)

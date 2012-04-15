import FWCore.ParameterSet.Config as cms

from HighMassAnalysis.Skimming.tauSelector_cfi import *
from HighMassAnalysis.Skimming.elecSelector_cfi import *

elecTauPairs = cms.EDProducer("DeltaRMinCandCombiner",
  decay = cms.string('selectedLooseHPSPatTau@+ selectedElectrons@-'),
  checkCharge = cms.bool(False),
  cut = cms.string( ''),
  name = cms.string('etauCandidates'),
  deltaRMin = cms.double(0.3)
)

selectedElecTauPairs = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('elecTauPairs'),
  minNumber = cms.uint32(1)                      
)

elecTauSkimSequence = cms.Sequence(
  ( selectedLooseHPSPatTau + selectedElectrons )
  * elecTauPairs
  * selectedElecTauPairs
)

import FWCore.ParameterSet.Config as cms

from HighMassAnalysis.Skimming.tauSelector_cfi import *
from HighMassAnalysis.Skimming.muonSelector_cfi import *
from HighMassAnalysis.Skimming.elecSelector_cfi import *

selectedTausInclusive = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('selectedLooseHPSPatTau'),
  minNumber = cms.uint32(1)                      
)

InclusiveTauSkimSequence = cms.Sequence(
    selectedLooseHPSPatTau
  * selectedTausInclusive
)

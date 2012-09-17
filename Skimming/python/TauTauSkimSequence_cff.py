import FWCore.ParameterSet.Config as cms

from HighMassAnalysis.Skimming.tauSelector_cfi import *

TauTauPairs = cms.EDProducer("DeltaRMinCandCombiner",
  decay = cms.string('selectedLooseHPSPatTau@+ selectedLooseHPSPatTau@-'),
  checkCharge = cms.bool(False),
  cut = cms.string( ''),
  name = cms.string('TauTauCandidates'),
  deltaRMin = cms.double(0.3)
)

selectedTauTauPairs = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('TauTauPairs'),
  minNumber = cms.uint32(1)                      
)

TauTauSkimSequence = cms.Sequence(
    selectedLooseHPSPatTau
  * TauTauPairs
  * selectedTauTauPairs
)

MuTauTauSkimSequence = cms.Sequence(
   selectedLooseHPSPatTau
  * TauTauPairs
  * selectedTauTauPairs
  * selectedMuons 
)

ElecTauTauSkimSequence = cms.Sequence(
   selectedLooseHPSPatTau
  * TauTauPairs
  * selectedTauTauPairs
  * selectedElectrons 
)

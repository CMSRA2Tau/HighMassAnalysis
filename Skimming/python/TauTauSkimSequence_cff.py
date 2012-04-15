import FWCore.ParameterSet.Config as cms

from HighMassAnalysis.Skimming.tauSelector_cfi import *

TauTauPairs = cms.EDProducer("DeltaRMinCandCombiner",
  decay = cms.string('selectedHPSPatTau@+ selectedLooseHPSPatTau@-'),
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
    selectedHPSTaus
  * selectedLooseHPSTaus
  * TauTauPairs
  * selectedTauTauPairs
)

MuTauTauSkimSequence = cms.Sequence(
   selectedLooseHPSPatTau
  * selectedHPSPatTau
  * TauTauPairs
  * selectedTauTauPairs
  * selectedMuons 
)

ElecTauTauSkimSequence = cms.Sequence(
   selectedLooseHPSPatTau
  * selectedHPSPatTau
  * TauTauPairs
  * selectedTauTauPairs
  * selectedElectrons 
)

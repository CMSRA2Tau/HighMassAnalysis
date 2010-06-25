import FWCore.ParameterSet.Config as cms

from HighMassAnalysis.Skimming.elecSelector_cfi import *
from HighMassAnalysis.Skimming.muonSelector_cfi import *
from HighMassAnalysis.Skimming.tauSelector_cfi import *

muTauPairs = cms.EDProducer("DeltaRMinCandCombiner",
  decay = cms.string('selectedPFTaus@+ selectedMuons@-'),
  checkCharge = cms.bool(False),
  cut = cms.string( ''),
  name = cms.string('muTauCandidates'),
  deltaRMin = cms.double(0.25)
)

selectedMuTauPairs = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('muTauPairs'),
  minNumber = cms.uint32(1)
)

selectedMuTauPairs2 = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('muTauPairs'),
  minNumber = cms.uint32(2)
)

elecTauPairs = cms.EDProducer("DeltaRMinCandCombiner",
  decay = cms.string('selectedPFTaus@+ selectedElectrons@-'),
  checkCharge = cms.bool(False),
  cut = cms.string( ''),
  name = cms.string('etauCandidates'),
  deltaRMin = cms.double(0.25)
)

selectedElecTauPairs = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('elecTauPairs'),
  minNumber = cms.uint32(1)
)

selectedElecTauPairs2 = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('elecTauPairs'),
  minNumber = cms.uint32(2)
)

eMuPairs = cms.EDProducer("DeltaRMinCandCombiner",
  decay = cms.string('selectedElectrons@+ selectedMuons@-'),
  checkCharge = cms.bool(False),
  cut = cms.string( ''),
  name = cms.string('eMuCandidates'),
  deltaRMin = cms.double(0.25)
)

selectedEMuPairs = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('eMuPairs'),
  minNumber = cms.uint32(1)                      
)

muMuPairs = cms.EDProducer("DeltaRMinCandCombiner",
  decay = cms.string('selectedMuons@+ selectedMuons@-'),
  checkCharge = cms.bool(False),
  cut = cms.string( ''),
  name = cms.string('muMuCandidates'),
  deltaRMin = cms.double(0.25)
)

selectedmuMuPairs = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('muMuPairs'),
  minNumber = cms.uint32(1)                      
)

elecElecPairs = cms.EDProducer("DeltaRMinCandCombiner",
  decay = cms.string('selectedElectrons@+ selectedElectrons@-'),
  checkCharge = cms.bool(False),
  cut = cms.string( ''),
  name = cms.string('elecElecCandidates'),
  deltaRMin = cms.double(0.25)
)

selectedelecElecPairs = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('elecElecPairs'),
  minNumber = cms.uint32(1)                      
)

muElecTauSkimSequence = cms.Sequence(
  ( selectedPFTaus + selectedElectrons + selectedMuons )
  * eMuPairs
  * selectedEMuPairs
  * muTauPairs
  * selectedMuTauPairs
  * elecTauPairs
  * selectedElecTauPairs
)

muMuTauSkimSequence = cms.Sequence(
  ( selectedPFTaus + selectedMuons )
  * muMuPairs
  * selectedmuMuPairs
  * muTauPairs
  * selectedMuTauPairs2
)

elecElecTauSkimSequence = cms.Sequence(
  ( selectedPFTaus + selectedElectrons )
  * elecElecPairs
  * selectedelecElecPairs
  * elecTauPairs
  * selectedElecTauPairs2
)

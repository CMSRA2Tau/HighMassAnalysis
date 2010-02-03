import FWCore.ParameterSet.Config as cms

from HighMassAnalysis.Skimming.elecSelector_cfi import *
from HighMassAnalysis.Skimming.muonSelector_cfi import *

eMuPairs = cms.EDProducer("DeltaRMinCandCombiner",
  decay = cms.string('selectedElectrons@+ selectedMuons@-'),
  checkCharge = cms.bool(False),
  cut = cms.string( ''),
  name = cms.string('eMuCandidates'),
  deltaRMin = cms.double(0.3)
)

selectedEMuPairs = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('eMuPairs'),
  minNumber = cms.uint32(1)                      
)

muElecSkimSequence = cms.Sequence(
  ( selectedElectrons + selectedMuons )
  * eMuPairs
  * selectedEMuPairs
)

muMuPairs = cms.EDProducer("DeltaRMinCandCombiner",
  decay = cms.string('selectedMuons@+ selectedMuons@-'),
  checkCharge = cms.bool(False),
  cut = cms.string( ''),
  name = cms.string('muMuCandidates'),
  deltaRMin = cms.double(0.3)
)

selectedmuMuPairs = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('muMuPairs'),
  minNumber = cms.uint32(1)                      
)

muMuSkimSequence = cms.Sequence(
  selectedMuons
  * muMuPairs
  * selectedmuMuPairs
)

elecElecPairs = cms.EDProducer("DeltaRMinCandCombiner",
  decay = cms.string('selectedElectrons@+ selectedElectrons@-'),
  checkCharge = cms.bool(False),
  cut = cms.string( ''),
  name = cms.string('elecElecCandidates'),
  deltaRMin = cms.double(0.3)
)

selectedelecElecPairs = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('elecElecPairs'),
  minNumber = cms.uint32(1)                      
)

elecElecSkimSequence = cms.Sequence(
  selectedElectrons
  * elecElecPairs
  * selectedelecElecPairs
)

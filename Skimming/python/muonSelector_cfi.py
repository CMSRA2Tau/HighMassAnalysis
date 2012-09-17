import FWCore.ParameterSet.Config as cms

selectedMuons = cms.EDFilter("MuonSelector",
  src = cms.InputTag('muons'),
  cut = cms.string("isGlobalMuon & pt > 8 & abs(eta) < 2.1"),
  filter = cms.bool(True)
)

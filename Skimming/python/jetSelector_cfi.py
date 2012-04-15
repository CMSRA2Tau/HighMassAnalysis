import FWCore.ParameterSet.Config as cms

selectedPatJets = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("patJets"),
    cut = cms.string(''),
    filter = cms.bool(True)
)

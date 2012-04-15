import FWCore.ParameterSet.Config as cms

WFilter = cms.EDFilter("MCGenEventFilter",
   GenSource = cms.InputTag("genParticles"),
   ParticleID = cms.untracked.int32(24),
   ParticleStatus = cms.untracked.int32(3),
   DaughterParticleID = cms.untracked.int32(0)
)

WtoMuNuFilter = cms.EDFilter("MCGenEventFilter",
   GenSource = cms.InputTag("genParticles"),
   ParticleID = cms.untracked.int32(24),
   ParticleStatus = cms.untracked.int32(3),
   DaughterParticleID = cms.untracked.int32(13)
)

WtoElecNuFilter = cms.EDFilter("MCGenEventFilter",
   GenSource = cms.InputTag("genParticles"),
   ParticleID = cms.untracked.int32(24),
   ParticleStatus = cms.untracked.int32(3),
   DaughterParticleID = cms.untracked.int32(11)
)

HiggsFilter = cms.EDFilter("MCGenEventFilter",
   GenSource = cms.InputTag("genParticles"),
   ParticleID = cms.untracked.int32(25),
   ParticleStatus = cms.untracked.int32(3),
   DaughterParticleID = cms.untracked.int32(0)
)

HiggsToTauTauFilter = cms.EDFilter("MCGenEventFilter",
   GenSource = cms.InputTag("genParticles"),
   ParticleID = cms.untracked.int32(25),
   ParticleStatus = cms.untracked.int32(3),
   DaughterParticleID = cms.untracked.int32(15)
)

WHToMuTauTauFilter = cms.Sequence(
  WFilter *
  WtoMuNuFilter *
  HiggsFilter *
  HiggsToTauTauFilter
)

WHToElecTauTauFilter = cms.Sequence(
  WFilter *
  WtoElecNuFilter *
  HiggsFilter *
  HiggsToTauTauFilter
)


import FWCore.ParameterSet.Config as cms

#################################################################################

genParticlesFromZPrimes = cms.EDProducer("GenParticlePruner",
  src = cms.InputTag("genParticles"),
  select = cms.vstring(
    "drop * ", # this is the default
    "keep+ pdgId = {Z'0}",
    "drop pdgId = {Z'0}"
  )
)

genElectronsFromZPrimes = cms.EDProducer("GenParticlePruner",
  src = cms.InputTag("genParticlesFromZPrimes"),
  select = cms.vstring(
    "drop * ", # this is the default
    "keep pdgId = {e+}",
    "keep pdgId = {e-}"
  )
)

genMuonsFromZPrimes = cms.EDProducer("GenParticlePruner",
  src = cms.InputTag("genParticlesFromZPrimes"),
  select = cms.vstring(
    "drop * ", # this is the default
    "keep pdgId = {mu+}",
    "keep pdgId = {mu-}"
  )
)

genTausFromZPrimes = cms.EDProducer("GenParticlePruner",
  src = cms.InputTag("genParticlesFromZPrimes"),
  select = cms.vstring(
    "drop * ", # this is the default
    "keep pdgId = {tau-}",
    "keep pdgId = {tau+}"
  )
)

produceGenDecayProductsFromZPrimes = cms.Sequence(genParticlesFromZPrimes * genElectronsFromZPrimes * genMuonsFromZPrimes * genTausFromZPrimes)



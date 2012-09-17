import FWCore.ParameterSet.Config as cms

selectedHPSPatTau = cms.EDFilter("PATTauSelector",
  src = cms.InputTag('selectedPatTaus'),
  cut = cms.string('pt > 15 & abs(eta) < 2.1 & tauID("decayModeFinding") > 0.5'),
#  cut = cms.string('pt > 15 & abs(eta) < 2.1 & tauID("decayModeFinding") > 0.5 & tauID("byLooseIsolationDeltaBetaCorr") > 0.5'),
  filter = cms.bool(True)
)

selectedLooseHPSPatTau = cms.EDFilter("PATTauSelector",
  src = cms.InputTag('selectedPatTaus'),
  cut = cms.string('pt > 15 & abs(eta) < 2.1 & tauID("decayModeFinding") > 0.5'),
  filter = cms.bool(True)
)

selectedVeryLooseHPSPatTau = cms.EDFilter("PATTauSelector",
  src = cms.InputTag('selectedPatTaus'),
  cut = cms.string('pt > 15 & abs(eta) < 2.1'),
  filter = cms.bool(True)
)


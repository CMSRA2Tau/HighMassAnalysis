import FWCore.ParameterSet.Config as cms
import copy

from Configuration.EventContent.EventContent_cff import *

skimEventContent = copy.deepcopy(AODSIMEventContent)

#tauAnalysisEventContent.outputCommands.extend(RecoEcalRECO.outputCommands)
#tauAnalysisEventContent.outputCommands.extend(RecoParticleFlowRECO.outputCommands)

# Keep track extra collection
allTrackExtra = cms.PSet(
  outputCommands = cms.untracked.vstring('keep recoTrackExtras_generalTracks_*_*'
  )
)

#--------------------------------------------------------------------------------
# keep all ECAL + HCAL recHits for computation of IsoDeposits
#--------------------------------------------------------------------------------
allEcalRecHits = cms.PSet(
  outputCommands = cms.untracked.vstring('keep *_ecalRecHit_EcalRecHitsEB_*',
                                         'keep *_ecalRecHit_EcalRecHitsEE_*',
                                         'keep *_ecalRecHit_EcalRecHitsES_*')
)
skimEventContent.outputCommands.extend(allEcalRecHits.outputCommands)

allHcalRecHits = cms.PSet(
  outputCommands = cms.untracked.vstring('keep *_hbhereco_*_*',
                                         'keep *_hfreco_*_*',
                                         'keep *_horeco_*_*')
)
skimEventContent.outputCommands.extend(allHcalRecHits.outputCommands)
skimEventContent.outputCommands.extend(allTrackExtra.outputCommands)

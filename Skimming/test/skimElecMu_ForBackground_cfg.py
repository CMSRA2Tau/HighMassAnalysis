import FWCore.ParameterSet.Config as cms

process = cms.Process("eMuSkim")

process.load("FWCore/MessageService/MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = 'INFO'
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

from HighMassAnalysis.Skimming.EventContent_cff import *

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
#process.GlobalTag.globaltag = 'STARTUP3XY_V9::All'
process.GlobalTag.globaltag = 'START3X_V16::All'

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)
#process.load("HighMassAnalysis.Skimming.FILESTOREAD")

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
        'dcache:/pnfs/cms/WAX/resilient/luiggi/Zprime/GEN-SIM-RECO/ZprimeTauTau500_GEN-SIM-RECO_0.root',
  )
)

# Generator cuts
process.load("HighMassAnalysis.Skimming.genLevelSequence_cff")

# Skim sequence
process.load("HighMassAnalysis.Skimming.leptonLeptonSkimSequence_cff")

process.eMuSkimPath = cms.Path(
#    process.genLevelElecMuSequence * 
    process.muElecSkimSequence
)

eMuEventSelection = cms.untracked.PSet(
  SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('eMuSkimPath')
  )
)

process.elecMuSkimOutputModule = cms.OutputModule("PoolOutputModule",                                 
  skimEventContent,                                               
  eMuEventSelection,
  fileName = cms.untracked.string("eMuSkim.root")
)

process.o = cms.EndPath( process.elecMuSkimOutputModule )


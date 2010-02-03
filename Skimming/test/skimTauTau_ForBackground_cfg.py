import FWCore.ParameterSet.Config as cms

process = cms.Process("leptonTauSkim")

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

process.PFTauProducerSequence = cms.Sequence(
    process.ic5PFJetTracksAssociatorAtVertex
 *  process.pfRecoTauTagInfoProducer
 *  process.produceAndDiscriminateFixedConeHighEffPFTaus
)

# Generator cuts
process.load("HighMassAnalysis.Skimming.genLevelSequence_cff")

#tau-tau Skim sequence & selection
process.load("HighMassAnalysis.Skimming.TauTauSkimSequence_cff")

process.TauTauSkimPath = cms.Path(
#    process.genLevelTauTauSequence *
    process.PFTauProducerSequence
 *  process.TauTauSkimSequence
)

TauTauEventSelection = cms.untracked.PSet(
  SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('TauTauSkimPath')
  )
)

process.TauTauSkimOutputModule = cms.OutputModule("PoolOutputModule",                                 
  skimEventContent,                                               
  TauTauEventSelection,
  fileName = cms.untracked.string("TauTauSkim.root")
)

process.o = cms.EndPath( process.TauTauSkimOutputModule )


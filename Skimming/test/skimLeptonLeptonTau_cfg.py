import FWCore.ParameterSet.Config as cms

process = cms.Process("leptonLeptonTauSkim")

process.load("FWCore/MessageService/MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = 'INFO'
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

from HighMassAnalysis.Skimming.EventContent_cff import *

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_31X_V3::All' # need to change to the correct one
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
 duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
  fileNames = cms.untracked.vstring(
	'/store/user/lpctau/HighMassTau/ZprimeTauTau400_7TeV_341_GEN-SIM-RECO/zprimeTauTau400_7TeV_341_GEN-SIM-RECO_0.root',
	'/store/user/lpctau/HighMassTau/ZprimeTauTau400_7TeV_341_GEN-SIM-RECO/zprimeTauTau400_7TeV_341_GEN-SIM-RECO_1.root'
  )
)

process.PFTauProducerSequence = cms.Sequence(
    process.ic5PFJetTracksAssociatorAtVertex
 *  process.pfRecoTauTagInfoProducer
 *  process.produceAndDiscriminateFixedConeHighEffPFTaus
)

process.load("HighMassAnalysis.Skimming.leptonLeptonTauSkimSequence_cff")

#e-mu-tau Skim sequence & selection
process.muElecTauSkimPath = cms.Path(
    process.PFTauProducerSequence
 *  process.muElecTauSkimSequence
)

#mu-mu-tau Skim sequence & selection
process.muMuTauSkimPath = cms.Path(
    process.PFTauProducerSequence
 *  process.muMuTauSkimSequence
)

#e-e-tau Skim sequence & selection
process.elecElecTauSkimPath = cms.Path(
    process.PFTauProducerSequence
 *  process.elecElecTauSkimSequence
)

leptonLeptonTauEventSelection = cms.untracked.PSet(
  SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('muElecTauSkimPath','muMuTauSkimPath','elecElecTauSkimPath')
  )
)

process.leptonLeptonTauSkimOutputModule = cms.OutputModule("PoolOutputModule",                                 
  skimEventContent,                                               
  leptonLeptonTauEventSelection,
  fileName = cms.untracked.string("leptonLeptonTauSkim.root")
)

process.o = cms.EndPath( process.leptonLeptonTauSkimOutputModule )


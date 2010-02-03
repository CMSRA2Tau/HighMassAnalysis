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
#process.GlobalTag.globaltag = 'START3X_V16::All'  # for 341 signal
process.GlobalTag.globaltag = 'MC_31X_V3::All'
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

#process.load("HighMassAnalysis.Skimming.FILESTOREAD")

process.source = cms.Source("PoolSource",
 duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
  fileNames = cms.untracked.vstring(
       #'/store/mc/Summer08/WJets-madgraph/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0016/FECF78F2-BBE9-DD11-8EAA-003048D15DCA.root'
       #'/store/user/eluiggi/mcOut/recoFiles/ZprimeTauTau500_GEN-SIM-RECO_0.root'
       #'/store/mc/Summer09/Ztautau/GEN-SIM-RECO/MC_31X_V3_7TeV-v1/0013/FA162037-EE9B-DE11-A028-001E0B5A645A.root'
	'/store/user/lpctau/HighMassTau/ZprimeTauTau400_7TeV_341_GEN-SIM-RECO/zprimeTauTau400_7TeV_341_GEN-SIM-RECO_0.root',
	'/store/user/lpctau/HighMassTau/ZprimeTauTau400_7TeV_341_GEN-SIM-RECO/zprimeTauTau400_7TeV_341_GEN-SIM-RECO_1.root'
  )
)

process.PFTauProducerSequence = cms.Sequence(
    process.ic5PFJetTracksAssociatorAtVertex
 *  process.pfRecoTauTagInfoProducer
 *  process.produceAndDiscriminateFixedConeHighEffPFTaus
)

# Generator cuts
process.load("HighMassAnalysis.Skimming.genLevelSequence_cff")

#mu-tau Skim sequence & selection
process.load("HighMassAnalysis.Skimming.muTauSkimSequence_cff")
process.muTauSkimPath = cms.Path(
    process.genLevelMuTauSequence *
    process.PFTauProducerSequence
 *  process.muTauSkimSequence
)

muTauEventSelection = cms.untracked.PSet(
  SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('muTauSkimPath')
  )
)

process.muTauSkimOutputModule = cms.OutputModule("PoolOutputModule",                                 
  skimEventContent,                                               
  muTauEventSelection,
  fileName = cms.untracked.string("zprimeMuTauSkim.root")
  #fileName = cms.untracked.string("muTau_outputFILENAME")
)


# e-tau Skim sequence & selection
process.load("HighMassAnalysis.Skimming.elecTauSkimSequence_cff")
process.elecTauSkimPath = cms.Path(
    process.genLevelElecTauSequence *
    process.PFTauProducerSequence
 *  process.elecTauSkimSequence
)


elecTauEventSelection = cms.untracked.PSet(
  SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('elecTauSkimPath')
  )
)

process.elecTauSkimOutputModule = cms.OutputModule("PoolOutputModule",                                 
  skimEventContent,                                               
  elecTauEventSelection,
  fileName = cms.untracked.string("zprimeETauSkim.root")
  #fileName = cms.untracked.string("eTau_outputFILENAME")
)

#tau-tau Skim sequence & selection
process.load("HighMassAnalysis.Skimming.TauTauSkimSequence_cff")
process.TauTauSkimPath = cms.Path(
    process.genLevelTauTauSequence *
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
  fileName = cms.untracked.string("zprimeTauTauSkim.root")
  #fileName = cms.untracked.string("tauTau_outputFILENAME")
)

# lepton-lepton skim sequence & selection
process.load("HighMassAnalysis.Skimming.leptonLeptonSkimSequence_cff")

#e-mu Skim sequence & selection

process.eMuSkimPath = cms.Path(
    process.genLevelElecMuSequence * 
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
  fileName = cms.untracked.string("zprimeEMuSkim.root")
  #fileName = cms.untracked.string("eMu_outputFILENAME")
)

#mu-mu Skim sequence & selection

process.muMuSkimPath = cms.Path(
   process.genLevelMuMuSequence * 
   process.muMuSkimSequence
)

muMuEventSelection = cms.untracked.PSet(
  SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('muMuSkimPath')
  )
)

process.muMuSkimOutputModule = cms.OutputModule("PoolOutputModule",                                 
  skimEventContent,                                               
  muMuEventSelection,
  fileName = cms.untracked.string("zprimeMuMuSkim.root")
  #fileName = cms.untracked.string("muMu_outputFILENAME")
)

#e-e Skim sequence & selection

process.elecElecSkimPath = cms.Path(
    process.genLevelElecElecSequence * 
    process.elecElecSkimSequence
)

elecElecEventSelection = cms.untracked.PSet(
  SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('elecElecSkimPath')
  )
)

process.elecElecSkimOutputModule = cms.OutputModule("PoolOutputModule",                                 
  skimEventContent,                                               
  elecElecEventSelection,
  fileName = cms.untracked.string("zprimeEESkim.root")
  #fileName = cms.untracked.string("eE_outputFILENAME")
)



process.o = cms.EndPath( 
 			process.muTauSkimOutputModule +  
		        process.elecTauSkimOutputModule + 
		        process.TauTauSkimOutputModule + 
		        process.elecMuSkimOutputModule + 
		        process.muMuSkimOutputModule + 
		        process.elecElecSkimOutputModule
		       )


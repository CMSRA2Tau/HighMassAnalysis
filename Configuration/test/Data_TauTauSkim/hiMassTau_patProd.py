import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
                  
options.register ('data',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Run this on real data")
    
options.register ('signal',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Is this the signal?")
    
options.register ('channel',
                  'mutau',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Release to use: 42x or 52x")

options.parseArguments()

data    = options.data
signal  = options.signal
channel = options.channel

process = cms.Process("PATTuple")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
process.load("Configuration.StandardSequences.MagneticField_cff")

#if(data):
#  process.GlobalTag.globaltag = 'GR_R_42_V12::All'
#else:
#  process.GlobalTag.globaltag = 'START42_V12::All'

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.source = cms.Source("PoolSource", 
     fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_5_2_3_patch3/RelValZTT/GEN-SIM-RECO/START52_V9_special_120410-v1/0122/6E6C7970-0283-E111-A6CB-003048FFD728.root'
      )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

## Output Module Configuration (expects a path 'p')
from HighMassAnalysis.Configuration.patTupleEventContentForHiMassTau_cff import *
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('skimPat.root'),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('drop *', *patTupleEventContent ),
    fastCloning = cms.untracked.bool(False)
)
process.outpath = cms.EndPath(process.out)

process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                     applyfilter = cms.untracked.bool(True),
                                     debugOn = cms.untracked.bool(False),
                                     numtrack = cms.untracked.uint32(10),
                                     thresh = cms.untracked.double(0.2)
                                     )

# trigger + Skim sequence
process.load("HighMassAnalysis.Skimming.triggerReq_cfi")
process.load("HighMassAnalysis.Skimming.genLevelSequence_cff")

if(channel == "emu"):
  process.load("HighMassAnalysis.Skimming.leptonLeptonSkimSequence_cff")
  process.theSkim = cms.Sequence(
    process.muElecSkimSequence
  )
  process.hltFilter = cms.Sequence(
    process.emuHLTFilter
  )
  process.genLevelSelection = cms.Sequence(
    process.genLevelElecMuSequence
  )
if(channel == "etau"):
  process.load("HighMassAnalysis.Skimming.elecTauSkimSequence_cff")
  process.theSkim = cms.Sequence(
    process.elecTauSkimSequence
  )
  process.hltFilter = cms.Sequence(
     process.etauHLTFilter
  )
  process.genLevelSelection = cms.Sequence(
    process.genLevelElecTauSequence
  )
if(channel == "mutau"):
  process.load("HighMassAnalysis.Skimming.muTauSkimSequence_cff")
  process.theSkim = cms.Sequence(
    process.muTauSkimSequence
  )
  process.hltFilter = cms.Sequence(
     process.mutauHLTFilter
  )
  process.genLevelSelection = cms.Sequence(
    process.genLevelMuTauSequence
  )
if(channel == "tautau"):
  process.load("HighMassAnalysis.Skimming.TauTauSkimSequence_cff")
  process.theSkim = cms.Sequence(
    process.TauTauSkimSequence
  )
  process.hltFilter = cms.Sequence(
     process.tautauHLTFilter
  )
  process.genLevelSelection = cms.Sequence(
    process.genLevelTauTauSequence
  )
if(channel == "mutautau"):
  process.load("HighMassAnalysis.Skimming.TauTauSkimSequence_cff")
  process.load("HighMassAnalysis.Skimming.WHTauTauGenLevel_cfi")
  process.theSkim = cms.Sequence(
    process.MuTauTauSkimSequence
  )
  process.hltFilter = cms.Sequence(
     process.mutauHLTFilter
  )
  process.genLevelSelection = cms.Sequence(
    process.WHToMuTauTauFilter +
    process.genLevelTauTauSequence
  )
if(channel == "electautau"):
  process.load("HighMassAnalysis.Skimming.TauTauSkimSequence_cff")
  process.load("HighMassAnalysis.Skimming.WHTauTauGenLevel_cfi")
  process.theSkim = cms.Sequence(
    process.ElecTauTauSkimSequence
  )
  process.hltFilter = cms.Sequence(
     process.etauHLTFilter
  )
  process.genLevelSelection = cms.Sequence(
    process.WHToElecTauTauFilter +
    process.genLevelTauTauSequence
  )
if(channel == "susy"):
  process.theSkim = cms.Sequence(  )
  process.hltFilter = cms.Sequence(
     process.SusyHLTFilter
  )
  process.genLevelSelection = cms.Sequence(  )

# Standard pat sequences
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# resolutions
process.load("TopQuarkAnalysis.TopObjectResolutions.stringResolutions_etEtaPhi_Summer11_cff")

# PF Met Corrections
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")

if(data):
  from PhysicsTools.PatAlgos.tools.coreTools import *
  removeMCMatching(process, ['All'])

# --------------------Modifications for taus--------------------

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

# --------------------Adding PF Isolation to Leptons--------------------

if(data):
  process.patMuons.addResolutions = cms.bool(False)
else:
  process.patMuons.addResolutions = cms.bool(True)
  process.patMuons.resolutions = cms.PSet( default = cms.string("muonResolution") )


from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
process.eleIsoSequence = setupPFElectronIso(process, 'selectedPatElectrons')
process.muIsoSequence = setupPFMuonIso(process, 'selectedPatMuons')

# --------------------Modifications for electrons--------------------

if(data):
  process.patElectrons.addResolutions = cms.bool(False)
else:
  process.patElectrons.addResolutions = cms.bool(True)
  process.patElectrons.resolutions = cms.PSet( default = cms.string("elecResolution") )


from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
process.heepPatElectrons = cms.EDProducer("HEEPAttStatusToPAT",
                                          eleLabel = cms.InputTag("selectedPatElectrons"),
                                          barrelCuts = cms.PSet(heepBarrelCuts),
                                          endcapCuts = cms.PSet(heepEndcapCuts),
                                          applyRhoCorrToEleIsol = cms.bool(False), 
                                          eleIsolEffectiveAreas = cms.PSet (
                                              trackerBarrel = cms.double(0.),
                                              trackerEndcap = cms.double(0.),
                                              ecalBarrel = cms.double(0.101),
                                              ecalEndcap = cms.double(0.046),
                                              hcalBarrel = cms.double(0.021),
                                              hcalEndcap = cms.double(0.040)
                                              ),
                                          eleRhoCorrLabel = cms.InputTag("kt6PFJets","rho"),
                                          )

# --------------------Modifications for jets--------------------

process.ak5PFJetsForPatJets = process.ak5PFJets.clone( 
                                doAreaFastjet = cms.bool(False)
                              )
    
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(process,cms.InputTag('ak5PFJetsForPatJets'),
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute'])),
                 doType1MET   = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True 
                 )

if(data):
  process.patJets.addResolutions = cms.bool(False)
else:
  process.patJets.addResolutions = cms.bool(True)
  process.patJets.resolutions = cms.PSet( default = cms.string("udscResolutionPF") )
  
# --------------------Modifications for MET--------------------

if(data):
  process.ak5PFJets.doAreaFastjet = cms.bool(True)
  process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
  process.metAnalysisSequence=cms.Sequence(process.producePFMETCorrections)
  process.patPFType1Type0METs = process.patMETs.clone(
                                      metSource = cms.InputTag('pfType1CorrectedMet')
                                      addMuonCorrections = cms.bool(False)
                                      addGenMET    = cms.bool(False)
                                )
  process.patPFType1Type2Type0METs = process.patMETs.clone(
                                      metSource = cms.InputTag('pfType1p2CorrectedMet')
                                      addMuonCorrections = cms.bool(False)
                                      addGenMET    = cms.bool(False)
                                )
else:
  process.ak5PFJets.doAreaFastjet = cms.bool(True)
  process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")
  process.metAnalysisSequence=cms.Sequence(process.producePFMETCorrections)
  process.patPFType1Type0METs = process.patMETs.clone(
                                      metSource = cms.InputTag('pfType1CorrectedMet')
                                      addMuonCorrections = cms.bool(False)
                                      addGenMET    = cms.bool(True)
                                )
  process.patPFType1Type2Type0METs = process.patMETs.clone(
                                      metSource = cms.InputTag('pfType1p2CorrectedMet')
                                      addMuonCorrections = cms.bool(False)
                                      addGenMET    = cms.bool(True)
                                      addResolutions = cms.bool(True)
                                      resolutions = cms.PSet( default = cms.string("metResolutionPF") )
                                )


# Let it run
if(data):
  process.p = cms.Path(
    process.scrapingVeto + 
    process.PFTau +
    process.ak5PFJetsForPatJets +
    process.metAnalysisSequence +
    process.patDefaultSequence + 
    process.eleIsoSequence +
    process.muIsoSequence +
    process.heepPatElectrons +
    process.patPFType1Type0METs +
    process.patPFType1Type2Type0METs +
    process.hltFilter +       
    process.theSkim	       
  )
else:
  if(signal):
    process.p = cms.Path(
      process.genLevelSelection +
      process.PFTau +
      process.ak5PFJetsForPatJets +
      process.metAnalysisSequence + 
      process.patDefaultSequence + 
      process.eleIsoSequence +
      process.muIsoSequence +
      process.heepPatElectrons +
      process.patPFType1Type0METs +
      process.patPFType1Type2Type0METs
    )
  else:
    process.p = cms.Path(
      process.PFTau +
      process.ak5PFJetsForPatJets +
      process.metAnalysisSequence +
      process.patDefaultSequence + 
      process.eleIsoSequence +
      process.muIsoSequence +
      process.heepPatElectrons +
      process.patPFType1Type0METs +
      process.patPFType1Type2Type0METs +
      process.theSkim	       
    )


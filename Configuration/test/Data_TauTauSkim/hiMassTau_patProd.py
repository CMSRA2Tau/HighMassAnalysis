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
                  "Desired channel")

options.parseArguments()

data    = options.data
signal  = options.signal
channel = options.channel
process = cms.Process("PATTuple")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

## Geometry and Detector Conditions (needed for a few patTuple production steps)
if data:
  process.load("Configuration.Geometry.GeometryIdeal_cff")
  process.load('Configuration.StandardSequences.Services_cff')
  process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# Jul13 re-reco dataset (A and B)
#  process.GlobalTag.globaltag = 'FT_53_V6_AN3::All'
# Aug06 recover/re-reco dataset
#  process.GlobalTag.globaltag = 'FT_53_V6C_AN3::All'
# 2012C Aug24 re-reco dataset
  process.GlobalTag.globaltag = 'FT_53_V10_AN3::All'
# 2012C promptreco dataset & Dec11
#  process.GlobalTag.globaltag = 'FT_P_V42C_AN3::All'
  process.load("Configuration.StandardSequences.MagneticField_cff")
  process.source = cms.Source("PoolSource", 
       fileNames = cms.untracked.vstring(
          '/store/relval/CMSSW_5_2_2/Jet/RECO/GR_R_52_V4_RelVal_jet2011B-v2/0252/96518387-A174-E111-95A6-001A928116E8.root'
        )
  )
else:
  process.load("Configuration.Geometry.GeometryIdeal_cff")
  process.load('Configuration.StandardSequences.Services_cff')
  process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
  process.GlobalTag.globaltag = 'START53_V7G::All'
  process.load("Configuration.StandardSequences.MagneticField_cff")
  process.source = cms.Source("PoolSource", 
       fileNames = cms.untracked.vstring(
          '/store/relval/CMSSW_5_3_2-START53_V6/RelValZTT/GEN-SIM-RECO/v1/0000/4AF6A8D2-7EB9-E111-9F26-003048678F78.root'
        )
  )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

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

if(channel == "mumu"):
  process.load("HighMassAnalysis.Skimming.leptonLeptonSkimSequence_cff")
  process.theSkim = cms.Sequence(
    process.muMuSkimSequence 
  )
  process.hltFilter = cms.Sequence(
    process.mumuHLTFilter
  )
  process.genLevelSelection = cms.Sequence( )
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
if(channel == "vbf"):
  process.theVBFSkim = cms.EDFilter("VBFSkim")
  process.theSkim = cms.Sequence( process.theVBFSkim )
  process.hltFilter = cms.Sequence( )
  process.genLevelSelection = cms.Sequence(  )

# Standard pat sequences
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# resolutions
#process.load("TopQuarkAnalysis.TopObjectResolutions.stringResolutions_etEtaPhi_Summer11_cff")

# PF Met Corrections
#process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
#process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")

if(data):
  from PhysicsTools.PatAlgos.tools.coreTools import *
  removeMCMatching(process, ['All'])

# --------------------Modifications for taus--------------------

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

# --------------------Adding PF Isolation to Leptons--------------------

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
process.eleIsoSequence = setupPFElectronIso(process, 'patElectrons')
process.muIsoSequence = setupPFMuonIso(process, 'patMuons')

# --------------------Modifications for electrons--------------------

from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
process.heepPatElectrons = cms.EDProducer("HEEPAttStatusToPAT",
                                          eleLabel = cms.InputTag("patElectrons"),
                                          barrelCuts = cms.PSet(heepBarrelCuts),
                                          endcapCuts = cms.PSet(heepEndcapCuts),
                                          applyRhoCorrToEleIsol = cms.bool(True), 
                                          eleIsolEffectiveAreas = cms.PSet (heepEffectiveAreas),
                                          eleRhoCorrLabel = cms.InputTag("kt6PFJetsForIsolation","rho"),
                                          verticesLabel = cms.InputTag("offlinePrimaryVerticesWithBS"),
                                          )

from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

# --------------------Modifications for jets--------------------

from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(process, cms.InputTag('ak5PFJets'),
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])),
                 doType1MET   = True,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True,
)
process.selectedPatJets.cut = cms.string("pt > 15 && abs(eta) < 5")
  
# --------------------Modifications for MET--------------------

process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfMETcorrType0'),
    cms.InputTag('pfJetMETcorr', 'type1')        
)


# Let it run
if(data):
  process.p = cms.Path(
    process.scrapingVeto + 
    process.type0PFMEtCorrection +
    process.recoTauClassicHPSSequence +
    process.patDefaultSequence + 
    process.kt6PFJetsForIsolation + 
    process.eleIsoSequence +
    process.muIsoSequence +
    process.heepPatElectrons +
    process.hltFilter +       
    process.theSkim
  )
else:
  if(signal):
    process.p = cms.Path(
      process.genLevelSelection +
      process.type0PFMEtCorrection +
      process.recoTauClassicHPSSequence +
      process.patDefaultSequence + 
      process.kt6PFJetsForIsolation + 
      process.eleIsoSequence +
      process.muIsoSequence +
      process.heepPatElectrons + 
      process.theSkim
    )
  else:
    process.p = cms.Path(
      process.type0PFMEtCorrection +
      process.recoTauClassicHPSSequence +
      process.patDefaultSequence + 
      process.kt6PFJetsForIsolation + 
      process.eleIsoSequence +
      process.muIsoSequence +
      process.heepPatElectrons +
      process.theSkim
    )


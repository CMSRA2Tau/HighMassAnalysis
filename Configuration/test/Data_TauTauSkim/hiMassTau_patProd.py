import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
                  
options.register ('data',
                  True,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Run this on real data")
    
options.register ('signal',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Is this the signal?")
    
options.register ('channel',
                  'susy',
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
  process.GlobalTag.globaltag = 'FT53_V21A_AN6::All'
  process.load("Configuration.StandardSequences.MagneticField_cff")
  process.source = cms.Source("PoolSource", 
       fileNames = cms.untracked.vstring(
#          '/store/relval/CMSSW_5_2_2/Jet/RECO/GR_R_52_V4_RelVal_jet2011B-v2/0252/96518387-A174-E111-95A6-001A928116E8.root'
#        '/store/data/Run2012A/SingleMu/AOD/22Jan2013-v1/20000/002F5062-346F-E211-BF00-1CC1DE04DF20.root'
'/store/data/Run2012D/SingleMu/AOD/22Jan2013-v1/10000/E47A49D5-A5A7-E211-9AA0-00259073E382.root'
        )
  )
else:
  process.load("Configuration.Geometry.GeometryIdeal_cff")
  process.load('Configuration.StandardSequences.Services_cff')
  process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
  process.GlobalTag.globaltag = 'START53_V27::All'
  process.load("Configuration.StandardSequences.MagneticField_cff")
  process.source = cms.Source("PoolSource", 
       fileNames = cms.untracked.vstring(
#          '/store/relval/CMSSW_5_3_2-START53_V6/RelValZTT/GEN-SIM-RECO/v1/0000/4AF6A8D2-7EB9-E111-9F26-003048678F78.root'
'/store/mc/Summer12_DR53X/WZJetsTo2Q2Nu_TuneZ2star_8TeV-madgraph-tauloa/AODSIM/PU_S10_START53_V7A-v1/0000/08878E54-90DD-E111-91F5-003048678B14.root'
        )
  )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

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

# filter out scraping events
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                     applyfilter = cms.untracked.bool(True),
                                     debugOn = cms.untracked.bool(False),
                                     numtrack = cms.untracked.uint32(10),
                                     thresh = cms.untracked.double(0.2)
                                     )

# filter out anomalous MET from detector noise, cosmic rays, and beam halo particles 
process.load("RecoMET.METFilters.metFilters_cff")

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
  process.hltFilter = cms.Sequence(  )
#     process.SusyHLTFilter
  process.genLevelSelection = cms.Sequence(  )
if(channel == "vbf"):
  process.theVBFSkim = cms.EDFilter("VBFSkim")
  process.theSkim = cms.Sequence( process.theVBFSkim )
  process.hltFilter = cms.Sequence( )
  process.genLevelSelection = cms.Sequence(  )
if(channel == "inclusiveTau"):
  process.load("HighMassAnalysis.Skimming.InclusiveTauSkimSequence_cff")
  process.theSkim = cms.Sequence(
    process.InclusiveTauSkimSequence
  )
  process.hltFilter = cms.Sequence( )
  process.genLevelSelection = cms.Sequence(  )

# Standard pat sequences
process.load("PhysicsTools.PatAlgos.patSequences_cff")

if(data):
  from PhysicsTools.PatAlgos.tools.coreTools import *
  removeMCMatching(process, ['All'])

# --------------------Modifications for taus--------------------

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

# --------------------Adding PF Isolation to Leptons--------------------

#from PhysicsTools.PatAlgos.tools.pfTools import *
#usePFIso( process )

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.muIsoSequence = setupPFMuonIso(process, 'muons')

process.patElectrons.isoDeposits = cms.PSet(
        pfChargedHadrons = cms.InputTag("elPFIsoDepositChargedPFIso"),
        pfChargedAll = cms.InputTag("elPFIsoDepositChargedAllPFIso"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoDepositPUPFIso"),
        pfNeutralHadrons = cms.InputTag("elPFIsoDepositNeutralPFIso"),
        pfPhotons = cms.InputTag("elPFIsoDepositGammaPFIso")
)
process.patElectrons.isolationValues = cms.PSet(
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFIdPFIso"),
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFIdPFIso"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFIdPFIso"),
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFIdPFIso"),
        pfPhotons = cms.InputTag("elPFIsoValueGamma03PFIdPFIso")
)
process.patElectrons.isolationValuesNoPFId = cms.PSet(
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03NoPFIdPFIso"),
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03NoPFIdPFIso"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03NoPFIdPFIso"),
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03NoPFIdPFIso"),
        pfPhotons = cms.InputTag("elPFIsoValueGamma03NoPFIdPFIso")
)

process.patMuons.isoDeposits = cms.PSet(
        pfChargedHadrons = cms.InputTag("muPFIsoDepositChargedPFIso"),
        pfChargedAll = cms.InputTag("muPFIsoDepositChargedAllPFIso"),
        pfPUChargedHadrons = cms.InputTag("muPFIsoDepositPUPFIso"),
        pfNeutralHadrons = cms.InputTag("muPFIsoDepositNeutralPFIso"),
        pfPhotons = cms.InputTag("muPFIsoDepositGammaPFIso")
)
process.patMuons.isolationValues = cms.PSet(
        pfChargedHadrons = cms.InputTag("muPFIsoValueCharged04PFIso"),
        pfChargedAll = cms.InputTag("muPFIsoValueChargedAll04PFIso"),
        pfPUChargedHadrons = cms.InputTag("muPFIsoValuePU04PFIso"),
        pfNeutralHadrons = cms.InputTag("muPFIsoValueNeutral04PFIso"),
        pfPhotons = cms.InputTag("muPFIsoValueGamma04PFIso")
)

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

process.patConversions = cms.EDProducer("PATConversionProducer",
    # this should be your last selected electron collection name since currently index is used to match with electron later.
    electronSource = cms.InputTag("heepPatElectrons")
)

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
process.load("CMGTools.External.pujetidsequence_cff")
  
# --------------------Modifications for MET--------------------

#process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
#process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")
#if data:
#  process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data
#else:
#  process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc
#
#process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
#process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
#process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
#    cms.InputTag('pfMETcorrType0'),
#    cms.InputTag('pfMEtSysShiftCorr'),
#    cms.InputTag('pfJetMETcorr', 'type1')        
#)

process.load("JetMETCorrections.Type1MET.correctionTermsCaloMet_cff")
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff")


if data:
  process.corrPfMetType1.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
else:
  process.corrPfMetType1.jetCorrLabel = cms.string("ak5PFL1FastL2L3")


process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType0PFCandidate_cff")
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType0RecoTrack_cff")
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetShiftXY_cff")

if data:
  process.corrPfMetShiftXY.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data
else:
  process.corrPfMetShiftXY.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc


process.load("JetMETCorrections.Type1MET.correctedMet_cff")

from PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi import patMETs
process.patPfMetT0pcT1Txy = patMETs.clone(
    metSource = cms.InputTag('pfMetT0pcT1Txy'),
    addMuonCorrections = cms.bool(False),
    addGenMET    = cms.bool(False)
)

# Let it run
if(data):
  process.p = cms.Path(
    process.scrapingVeto + 
    process.metFilters + 

#    process.type0PFMEtCorrection + 
#    process.pfMEtSysShiftCorrSequence + 

    process.correctionTermsPfMetType1Type2 +
    process.correctionTermsPfMetType0RecoTrack +
    process.correctionTermsPfMetType0PFCandidate +
    process.correctionTermsPfMetShiftXY +
#    process.correctionTermsCaloMet +
#    process.caloMetT1 + 
#    process.caloMetT1T2 + 
    process.pfMetT0rt +
    process.pfMetT0rtT1 +
    process.pfMetT0pc +
    process.pfMetT0pcT1 +
    process.pfMetT0rtTxy +
    process.pfMetT0rtT1Txy +
    process.pfMetT0pcTxy +
    process.pfMetT0pcT1Txy +
    process.pfMetT1 +
    process.pfMetT1Txy + 

    process.pfParticleSelectionSequence +
    process.eleIsoSequence +
    process.muIsoSequence +
    process.recoTauClassicHPSSequence + 
    process.kt6PFJetsForIsolation + 
    process.patDefaultSequence + 

    process.patPfMetT0pcT1Txy + 

    process.puJetIdSqeuence + 
    process.heepPatElectrons + 
    process.patConversions + 
    process.hltFilter +       
    process.theSkim
  )
else:
  if(signal):
    process.p = cms.Path(
      process.genLevelSelection +
      process.scrapingVeto + 
      process.metFilters + 

#      process.type0PFMEtCorrection + 
#      process.pfMEtSysShiftCorrSequence + 

      process.correctionTermsPfMetType1Type2 +
      process.correctionTermsPfMetType0RecoTrack +
      process.correctionTermsPfMetType0PFCandidate +
      process.correctionTermsPfMetShiftXY +
#      process.correctionTermsCaloMet +
#      process.caloMetT1 +
#      process.caloMetT1T2 +
      process.pfMetT0rt +
      process.pfMetT0rtT1 +
      process.pfMetT0pc +
      process.pfMetT0pcT1 +
      process.pfMetT0rtTxy +
      process.pfMetT0rtT1Txy +
      process.pfMetT0pcTxy +
      process.pfMetT0pcT1Txy +
      process.pfMetT1 +
      process.pfMetT1Txy + 

      process.pfParticleSelectionSequence +
      process.eleIsoSequence +
      process.muIsoSequence +
      process.recoTauClassicHPSSequence +
      process.kt6PFJetsForIsolation + 
      process.patDefaultSequence + 

      process.patPfMetT0pcT1Txy +

      process.puJetIdSqeuence + 
      process.heepPatElectrons + 
      process.patConversions + 
      process.theSkim
    )
  else:
    process.p = cms.Path(
      process.scrapingVeto + 
      process.metFilters + 

#      process.type0PFMEtCorrection +
#      process.pfMEtSysShiftCorrSequence + 

      process.correctionTermsPfMetType1Type2 +
      process.correctionTermsPfMetType0RecoTrack +
      process.correctionTermsPfMetType0PFCandidate +
      process.correctionTermsPfMetShiftXY +
#      process.correctionTermsCaloMet +
#      process.caloMetT1 +
#      process.caloMetT1T2 +
      process.pfMetT0rt +
      process.pfMetT0rtT1 +
      process.pfMetT0pc +
      process.pfMetT0pcT1 +
      process.pfMetT0rtTxy +
      process.pfMetT0rtT1Txy +
      process.pfMetT0pcTxy +
      process.pfMetT0pcT1Txy +
      process.pfMetT1 +
      process.pfMetT1Txy +

      process.pfParticleSelectionSequence +
      process.eleIsoSequence +
      process.muIsoSequence +
      process.recoTauClassicHPSSequence +
      process.kt6PFJetsForIsolation + 
      process.patDefaultSequence + 

      process.patPfMetT0pcT1Txy +

      process.puJetIdSqeuence + 
      process.heepPatElectrons +
      process.patConversions + 
      process.theSkim
    )

import FWCore.ParameterSet.Config as cms
import copy

### This files controls whether we are working with data/MC, signalMC (for skim), and channel.
### Make sure to set it to the desired process.
from HighMassAnalysis.Configuration.hiMassSetup_cfi import *

process = cms.Process('hiMassTau')

# import of standard configurations for RECOnstruction
process.load('Configuration/StandardSequences/Services_cff')

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = 'INFO'
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# standard sequences
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_noesprefer_cff')
if(data):
#  process.GlobalTag.globaltag = 'GR_P_V16::All'
  process.GlobalTag.globaltag = 'GR_R_42_V12::All'
else:
#  process.GlobalTag.globaltag = 'START311_V1::All'
  process.GlobalTag.globaltag = 'START42_V12::All'

# import particle data table - needed for print-out of generator level information
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# import sequence for PAT-tuple production
process.load("HighMassAnalysis.Configuration.producePatTuple_cff")

# import event-content definition of products to be stored in patTuple
#from HighMassAnalysis.Configuration.patTupleEventContentForHiMassTau_cff import *
from HighMassAnalysis.Configuration.patTupleEventContentForHiMassTau_cff import patTupleEventContent

#process.saveThePatTuple = cms.OutputModule("PoolOutputModule",
#    patTupleEventContent,                                               
#    fastCloning = cms.untracked.bool(False),
#    fileName = cms.untracked.string('skimPat.root')
#)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('skimPat.root'),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('drop *', *patTupleEventContent ),
    fastCloning = cms.untracked.bool(False)
)

process.outpath = cms.EndPath(process.out)

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32( 100 )
)

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
#        '/store/relval/CMSSW_4_2_0/RelValZTT/GEN-SIM-RECO/START42_V9-v1/0059/607DAB3A-0F5F-E011-9BFB-003048678FA6.root',
#        '/store/relval/CMSSW_4_2_0/RelValZTT/GEN-SIM-RECO/START42_V9-v1/0054/38039DE1-735E-E011-B7CA-0018F3D09616.root',
#        '/store/relval/CMSSW_4_2_0/RelValZTT/GEN-SIM-RECO/START42_V9-v1/0054/107DB9B4-7D5E-E011-91E9-001A92810AEA.root'
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/FEF9BFD2-8C7E-E011-9DD4-0026189438F7.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/FEBD919F-8E7E-E011-9306-003048678A6A.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/FE9386E5-8C7E-E011-8249-002618943947.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/FE76039F-877E-E011-89B1-0026189438BA.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/FCFD375D-8D7E-E011-AA64-003048678D52.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/FC6B2925-897E-E011-A5B0-00261894396E.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/FC39BCA1-887E-E011-80FF-0026189438CE.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/FABE5404-887E-E011-931F-00248C55CC7F.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/FA5D7982-8C7E-E011-AD5B-0018F3D096EA.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/F4C5D1E2-8F7E-E011-9448-003048678F0C.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/F4830D13-897E-E011-B42A-00261894393C.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/F4179404-8C7E-E011-9255-0026189438C0.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/F28B9311-8B7E-E011-9A65-003048678F74.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/F26CE071-8E7E-E011-992F-00304867C0EA.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/F24025BE-887E-E011-B703-002618943831.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/F0A8ABA4-887E-E011-B510-002618943898.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/F05FE181-877E-E011-A0FB-00261894397B.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/ECFC2E34-8A7E-E011-AEE5-0026189438B8.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/ECF50A9C-897E-E011-9C96-001A92971AA4.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/ECEE933B-8B7E-E011-8B62-0026189437ED.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/ECED21FE-8B7E-E011-8B77-003048678E2A.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/EC8AAED8-8C7E-E011-AB2C-00248C55CC97.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/EC3A4259-877E-E011-8B6A-00261894393D.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/EAEB4106-8B7E-E011-A778-001A928116D4.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/EA3824D5-8C7E-E011-949F-0026189438E0.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/E83DA30A-8F7E-E011-BAD3-003048678F92.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/E6A58E11-917E-E011-899F-00304867920A.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/E4E42950-8E7E-E011-B912-001A92810AF2.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/E4CF8E66-8C7E-E011-AE31-00304867C1BA.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/E4118CCC-8C7E-E011-BEF8-003048679248.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/E2EDF24F-8C7E-E011-B0FE-002618943810.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/E294D986-8B7E-E011-8E3C-0030486790B0.root',
        '/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0003/E206F63E-8B7E-E011-B6C1-001A928116DE.root',
)
    #skipBadFiles = cms.untracked.bool(True) 
)

#process.load("CondCore.DBCommon.CondDBCommon_cfi")
#process.jec = cms.ESSource("PoolDBESSource",
#      DBParameters = cms.PSet(
#        messageLevel = cms.untracked.int32(0)
#        ),
#      timetype = cms.string('runnumber'),
#      toGet = cms.VPSet(
#          cms.PSet(
#                record = cms.string('JetCorrectionsRecord'),
#                tag    = cms.string('JetCorrectorParametersCollection_Jec10V3_AK5PF'),
#                label  = cms.untracked.string('AK5PF')),
#          cms.PSet(
#                record = cms.string('JetCorrectionsRecord'),
#                tag    = cms.string('JetCorrectorParametersCollection_Jec10V3_AK5Calo'),
#                label  = cms.untracked.string('AK5Calo')),
#      ),
#      connect = cms.string('sqlite_file:Jec10V3.db')
#)
#process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                     applyfilter = cms.untracked.bool(True),
                                     debugOn = cms.untracked.bool(False),
                                     numtrack = cms.untracked.uint32(10),
                                     thresh = cms.untracked.double(0.2)
                                     )

# include particle flow based MET
from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')

# include particle flow based jets
from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process,
                cms.InputTag('ak5PFJets'),         # Jet collection;must be already in the event when patLayer0 sequence is executed
                'AK5', 'PF',
                doJTA=True,            # Run Jet-Track association & JetCharge
                doBTagging=True,       # Run b-tagging
                jetCorrLabel=('AK5PF', ['L1FastJet', 'L2Relative','L3Absolute']),
                doType1MET=False,
                doL1Cleaning=False,
                doL1Counters=False,
                genJetCollection = cms.InputTag(""),
                doJetID = False
                )

from PhysicsTools.PatAlgos.tools.electronTools import *
#addElectronUserIsolation(process,["Tracker"])
#addElectronUserIsolation(process,["Ecal"])
#addElectronUserIsolation(process,["Hcal"])

from PhysicsTools.PatAlgos.tools.muonTools import *
addMuonUserIsolation(process)

# trigger + Skim sequence
process.load("HighMassAnalysis.Skimming.triggerReq_cfi")

if(channel == "emu"):
  process.load("HighMassAnalysis.Skimming.leptonLeptonSkimSequence_cff")
  process.theSkim = cms.Sequence(
    process.muElecSkimSequence
  )
  process.hltFilter = cms.Sequence(
    process.emuHLTFilter
  )
if(channel == "etau"):
  process.load("HighMassAnalysis.Skimming.elecTauSkimSequence_cff")
  process.theSkim = cms.Sequence(
    process.elecTauSkimSequence
  )
  process.hltFilter = cms.Sequence(
     process.etauHLTFilter
  )
if(channel == "mutau"):
  process.load("HighMassAnalysis.Skimming.muTauSkimSequence_cff")
  process.theSkim = cms.Sequence(
    process.muTauSkimSequence
  )
  process.hltFilter = cms.Sequence(
     process.mutauHLTFilter
  )
if(channel == "tautau"):
  process.load("HighMassAnalysis.Skimming.TauTauSkimSequence_cff")
  process.theSkim = cms.Sequence(
    process.TauTauSkimSequence
  )
  process.hltFilter = cms.Sequence(
     process.tautauHLTFilter
  )
if(channel == "susy"):
  process.theSkim = cms.Sequence(  )
  process.hltFilter = cms.Sequence(
     process.SusyHLTFilter
  )

if(data):
  process.p = cms.Path(
    process.scrapingVeto + 
    process.hltFilter +       
    process.produceAndDiscriminatePFTausCustomized +
    process.theSkim +	       
    process.producePatTuple
  )
else:
  if(signal):
    process.p = cms.Path(
      process.produceAndDiscriminatePFTausCustomized +
      process.producePatTuple
    )
  else:
    process.p = cms.Path(
      process.produceAndDiscriminatePFTausCustomized +
      process.theSkim +	       
      process.producePatTuple
    )


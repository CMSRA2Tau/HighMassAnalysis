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
  process.GlobalTag.globaltag = 'GR10_P_V11::All'
else:
  process.GlobalTag.globaltag = 'START38_V12::All'

# import particle data table - needed for print-out of generator level information
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# import sequence for PAT-tuple production
process.load("HighMassAnalysis.Configuration.producePatTuple_cff")

# import event-content definition of products to be stored in patTuple
from HighMassAnalysis.Configuration.patTupleEventContentForHiMassTau_cff import *

process.savePatTuple = cms.OutputModule("PoolOutputModule",
    patTupleEventContent,                                               
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('skimPat.root')
)

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32( -1 )
)

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/FE79D77E-1BEE-DF11-B185-002618943932.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/EADAAF18-1BEE-DF11-8C80-0018F3D0969C.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/E0EA85B8-1CEE-DF11-BEC1-002618943944.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/DE0F3721-1CEE-DF11-830E-001A92971B06.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/C83CE2DC-1BEE-DF11-BE47-003048678B76.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/C4AED029-1DEE-DF11-A5BA-0018F3D09626.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/C41B33DF-1AEE-DF11-AA19-001A92811748.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/BAFA850A-1DEE-DF11-B3BD-00248C65A3EC.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/A2B6E85D-1DEE-DF11-8CDB-0026189437FA.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/A02A2523-1DEE-DF11-8066-002618943932.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/96939707-1CEE-DF11-8F9C-001A928116B0.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/94432DD8-1BEE-DF11-A735-003048678D86.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/92AB1A0B-1DEE-DF11-B2B0-003048679070.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/9008D6F5-1BEE-DF11-B61A-003048678FB2.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/821A9F6A-1BEE-DF11-9060-003048679070.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/688FAB07-1CEE-DF11-97AF-001BFCDBD1BE.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/489BC0DE-1BEE-DF11-BCF8-002618943810.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/402BBF29-1DEE-DF11-8DF0-0026189437FD.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/38929044-1DEE-DF11-9547-0018F3D096A6.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/3683BA3F-1DEE-DF11-B9E0-001A92810AD8.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/288AF2B1-1AEE-DF11-9AAC-0018F3D0965E.root',
#        '/store/data/Run2010A/JetMETTau/RECO/Nov4ReReco_v1/0153/2656870A-1DEE-DF11-8259-0018F3D0969C.root',

#        '/store/data/Run2010A/Mu/RECO/Nov4ReReco_v1/0011/82BC70AF-50EC-DF11-AB25-00237DA1FD7C.root',
#        '/store/data/Run2010A/Mu/RECO/Nov4ReReco_v1/0010/FEDA8FD8-36EC-DF11-95FA-00237DA1985A.root',
#        '/store/data/Run2010A/Mu/RECO/Nov4ReReco_v1/0010/FEBE31FC-F1EB-DF11-ADC6-00237DA0F456.root',
#        '/store/data/Run2010A/Mu/RECO/Nov4ReReco_v1/0010/FE84B981-0AEC-DF11-B31B-001E0B5F27BC.root',
#        '/store/data/Run2010A/Mu/RECO/Nov4ReReco_v1/0010/FC305B3D-26EC-DF11-A581-1CC1DE0437C8.root',
#        '/store/data/Run2010A/Mu/RECO/Nov4ReReco_v1/0010/FAE7513C-26EC-DF11-9423-1CC1DE0590E8.root',
#        '/store/data/Run2010A/Mu/RECO/Nov4ReReco_v1/0010/FACE7FA0-25EC-DF11-B161-00237DA1A548.root',
#        '/store/data/Run2010A/Mu/RECO/Nov4ReReco_v1/0010/FA6ECDEE-25EC-DF11-8394-001F29C464E4.root',
#        '/store/data/Run2010A/Mu/RECO/Nov4ReReco_v1/0010/F8AE46C1-0AEC-DF11-8BA7-1CC1DE04FF48.root',
#        '/store/data/Run2010A/Mu/RECO/Nov4ReReco_v1/0010/F8781E40-26EC-DF11-A041-001F2969CD10.root',
#        '/store/data/Run2010A/Mu/RECO/Nov4ReReco_v1/0010/F82A9E86-26EC-DF11-AA69-1CC1DE05B0C8.root',
#        '/store/data/Run2010A/Mu/RECO/Nov4ReReco_v1/0010/F6FE1928-ECEB-DF11-933D-00237DA1ED4A.root',
#        '/store/data/Run2010A/Mu/RECO/Nov4ReReco_v1/0010/F6FC6DF2-24EC-DF11-85C9-0017A4770424.root',
#        '/store/data/Run2010A/Mu/RECO/Nov4ReReco_v1/0010/F64C1269-0AEC-DF11-99DE-00237DA1EDE0.root',

        '/store/mc/Fall10/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola/GEN-SIM-RECO/START38_V12-v1/0004/E6D9BDE2-86C8-DF11-9827-00215E221782.root',
        '/store/mc/Fall10/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola/GEN-SIM-RECO/START38_V12-v1/0004/D6116835-8FC8-DF11-BE6C-00215E21D64E.root',
        '/store/mc/Fall10/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola/GEN-SIM-RECO/START38_V12-v1/0004/D60E54A5-7CC8-DF11-9A02-E41F131816A0.root',
        '/store/mc/Fall10/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola/GEN-SIM-RECO/START38_V12-v1/0004/D05136F1-D4C8-DF11-9F4B-E41F13181AB4.root',
        '/store/mc/Fall10/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola/GEN-SIM-RECO/START38_V12-v1/0004/AE387075-C5C8-DF11-8784-00215E21DCF6.root',
        '/store/mc/Fall10/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola/GEN-SIM-RECO/START38_V12-v1/0004/A8AED05A-7EC8-DF11-AE69-00215E21D972.root',
        '/store/mc/Fall10/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola/GEN-SIM-RECO/START38_V12-v1/0004/9C43BDD4-89C8-DF11-839B-00215E221812.root',
        '/store/mc/Fall10/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola/GEN-SIM-RECO/START38_V12-v1/0004/6E1C9934-8FC8-DF11-8781-00215E21D64E.root',
    )
    #skipBadFiles = cms.untracked.bool(True) 
)

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
addJetCollection(process,cms.InputTag('ak5PFJets'),
                  'AK5', 'PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5','PF'),
                 doType1MET   = False,
                 doL1Cleaning = False,                 
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID          = False
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

if(data):
  process.p = cms.Path(
    process.scrapingVeto + 
    process.hltFilter +       
    process.theSkim +	       
    process.producePatTuple +
    process.savePatTuple	       
  )
else:
  if(signal):
    process.p = cms.Path(
      process.producePatTuple +
      process.savePatTuple	       
    )
  else:
    process.p = cms.Path(
      process.theSkim +	       
      process.producePatTuple +
      process.savePatTuple	       
    )

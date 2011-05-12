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
  process.GlobalTag.globaltag = 'GR_P_V16::All'
else:
  process.GlobalTag.globaltag = 'START311_V1::All'

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
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0008/8890F217-334F-E011-BABA-0022198F5B49.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0005/145708DA-C74E-E011-BBAF-E0CB4E55368D.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0003/FE917076-4E4E-E011-94C7-90E6BA442F0C.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0003/FE19F8A6-4A4E-E011-9B20-0030487E3026.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0003/FAAFA997-434E-E011-81E4-E0CB4E29C4FF.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0003/FAADF2DE-464E-E011-BCA4-90E6BA442F2A.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0003/F69DC08F-4A4E-E011-B44E-E0CB4E1A117B.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0003/F2770B3C-414E-E011-9AB1-90E6BA442F0A.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0003/F07175A8-4A4E-E011-B9C0-003048D28EB2.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0003/ECAABCB8-4A4E-E011-85C3-E0CB4E29C4D4.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0003/EC89FA70-574E-E011-997E-E0CB4E19F9B2.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0003/EC5B5088-404E-E011-9938-485B39800B90.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0003/EA8AB777-4E4E-E011-B891-E0CB4E19F9A6.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0003/E8FEF73B-414E-E011-9C9A-E0CB4E29C4ED.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0003/E8478397-434E-E011-A14C-90E6BA442F00.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0003/E83672B4-484E-E011-AA1D-E0CB4E1A118B.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0003/E454EE96-434E-E011-9747-E0CB4E1A1152.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0003/E44015E0-464E-E011-9707-E0CB4E19F96B.root',
        
'/store/mc/Spring11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0003/DAB69276-4E4E-E011-9D64-E0CB4E1A117E.root',    
)
    #skipBadFiles = cms.untracked.bool(True) 
)

process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
          cms.PSet(
                record = cms.string('JetCorrectionsRecord'),
                tag    = cms.string('JetCorrectorParametersCollection_Jec10V3_AK5PF'),
                label  = cms.untracked.string('AK5PF')),
          cms.PSet(
                record = cms.string('JetCorrectionsRecord'),
                tag    = cms.string('JetCorrectorParametersCollection_Jec10V3_AK5Calo'),
                label  = cms.untracked.string('AK5Calo')),
      ),
      connect = cms.string('sqlite_file:Jec10V3.db')
)
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

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
  process.hltFilter = cms.Sequence(
     process.SusyHLTFilter
  )

if(data):
  process.p = cms.Path(
    process.scrapingVeto + 
    process.hltFilter +       
    process.theSkim +	       
    process.producePatTuple
#    process.saveThePatTuple	       
  )
else:
  if(signal):
    process.p = cms.Path(
      process.producePatTuple
#      process.saveThePatTuple	       
    )
  else:
    process.p = cms.Path(
      process.theSkim +	       
      process.producePatTuple
#      process.saveThePatTuple	       
    )


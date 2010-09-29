import FWCore.ParameterSet.Config as cms
import copy

signal = 0

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
process.GlobalTag.globaltag = 'START36_V10::All'
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_noesprefer_cff')
#process.GlobalTag.globaltag = 'STARTUP3XY_V9::All'

# import particle data table - needed for print-out of generator level information
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# import sequence for PAT-tuple production
process.load("HighMassAnalysis.Configuration.producePatTuple_cff")

# import event-content definition of products to be stored in patTuple
from HighMassAnalysis.Configuration.patTupleEventContentForHiMassTau_cff import *

process.savePatTuple = cms.OutputModule("PoolOutputModule",
    patTupleEventContent,                                               
    fileName = cms.untracked.string('zTauTauElecTauSkimPat.root')
)

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32( -1 )
)

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
       '/store/mc/Summer10/QCD_EMEnriched_Pt20to30/GEN-SIM-RECO/START36_V9_S09-v2/0009/EE61B95C-8088-DF11-9003-002618943930.root'
    )
    #skipBadFiles = cms.untracked.bool(True) 
)
process.source.inputCommands = cms.untracked.vstring(
	"keep *", 
	"drop *_MEtoEDMConverter_*_*", 
	"drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__HLT"
)

# Skimming step
process.load("HighMassAnalysis.Skimming.genLevelSequence_cff")
# Skim sequence
process.load("HighMassAnalysis.Skimming.elecTauSkimSequence_cff")

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
addElectronUserIsolation(process,["Tracker"])
addElectronUserIsolation(process,["Ecal"])
addElectronUserIsolation(process,["Hcal"])

from PhysicsTools.PatAlgos.tools.muonTools import *
addMuonUserIsolation(process)

if(signal):
	process.p = cms.Path( 
		      process.producePatTuple +
                      process.savePatTuple
	)
else:
	process.p = cms.Path( 
		      process.elecTauSkimSequence +
		      process.producePatTuple +
                      process.savePatTuple
	)
	

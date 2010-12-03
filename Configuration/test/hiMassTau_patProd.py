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
  #process.GlobalTag.globaltag = 'GR10_P_V7::All'
else:
  process.GlobalTag.globaltag = 'START36_V10::All'

# import particle data table - needed for print-out of generator level information
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# import sequence for PAT-tuple production
process.load("HighMassAnalysis.Configuration.producePatTuple_cff")

# import event-content definition of products to be stored in patTuple
from HighMassAnalysis.Configuration.patTupleEventContentForHiMassTau_cff import *

process.savePatTuple = cms.OutputModule("PoolOutputModule",
    patTupleEventContent,                                               
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('eMuSkimPat.root')
)

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32( -1 )
)

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
      '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/711/F47056D9-8DE7-DF11-9355-003048F024FA.root',
    )
    #skipBadFiles = cms.untracked.bool(True) 
)

process.source.inputCommands = cms.untracked.vstring(
       "keep *", 
       "drop *_MEtoEDMConverter_*_*", 
       "drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__HLT"
)

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

process.PFTauProducerSequence = cms.Sequence(
      process.ic5PFJetTracksAssociatorAtVertex *
      process.pfRecoTauTagInfoProducer *
      process.shrinkingConePFTauProducer*
      process.shrinkingConePFTauDiscriminationByLeadingTrackFinding*
      process.shrinkingConePFTauDiscriminationByLeadingTrackPtCut*
      process.shrinkingConePFTauDiscriminationByLeadingPionPtCut*
      process.shrinkingConePFTauDiscriminationByIsolation*
      process.shrinkingConePFTauDiscriminationByTrackIsolation*
      process.shrinkingConePFTauDiscriminationByECALIsolation*
      process.shrinkingConePFTauDiscriminationByIsolationUsingLeadingPion*
      process.shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion*
      process.shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion*
      process.shrinkingConePFTauDiscriminationAgainstElectron*
      process.shrinkingConePFTauDiscriminationAgainstMuon
)

# include particle flow based MET
from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')

# include particle flow based jets
from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process,cms.InputTag('ak5PFJets'),
		'AK5', 'PF',
		doJTA	      = True,
		doBTagging    = True,
		jetCorrLabel  = ('AK5','PF'),
		doType1MET    = False,
		doL1Cleaning  = False,  	       
		doL1Counters  = False,
		genJetCollection=cms.InputTag("ak5GenJets"),
		doJetID       = False
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


if(signal):
  process.p = cms.Path(
    process.producePatTuple +
    process.savePatTuple	       
  )
else:
  process.p = cms.Path(
    process.hltFilter +       
    process.theSkim +	       
    process.producePatTuple +
    process.savePatTuple	       
  )
	

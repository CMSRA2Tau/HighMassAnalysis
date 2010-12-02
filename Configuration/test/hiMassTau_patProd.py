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
process.GlobalTag.globaltag = 'GR10_P_V11::All'
#process.GlobalTag.globaltag = 'GR10_P_V7::All'
#process.GlobalTag.globaltag = 'START36_V10::All'
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
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('muTauSkimPat.root')
)

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32( -1 )
)

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
#        '/store/mc/Fall10/WtoMuNu_TuneD6T_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0013/F8D51F20-7DCF-DF11-9B87-E0CB4E55365A.root',
#        '/store/mc/Fall10/WtoMuNu_TuneD6T_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0013/F46C7DA7-7DCF-DF11-B814-001E0B8D7834.root',
#        '/store/mc/Fall10/WtoMuNu_TuneD6T_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0013/EEBEE73D-8BCF-DF11-8187-001E682F8738.root',
#        '/store/mc/Fall10/WtoMuNu_TuneD6T_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0013/EC7C93E7-FCCE-DF11-80C9-00261834B5B4.root',
#        '/store/mc/Fall10/WtoMuNu_TuneD6T_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0013/E25F4674-7FCF-DF11-AAD0-0015172C08D4.root',
#        '/store/mc/Fall10/WtoMuNu_TuneD6T_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0013/D6D24D87-8ACF-DF11-9398-001E68A99494.root',
#        '/store/mc/Fall10/WtoMuNu_TuneD6T_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0013/D6967C7E-7DCF-DF11-B711-00261834B574.root',
#        '/store/mc/Fall10/WtoMuNu_TuneD6T_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0013/CEBB3102-7DCF-DF11-A426-00D0680B883B.root',
#        '/store/mc/Fall10/WtoMuNu_TuneD6T_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0013/CEA0AF0E-8BCF-DF11-8D32-001E682F84DE.root',
#        '/store/mc/Fall10/WtoMuNu_TuneD6T_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0013/B62C68A6-8ACF-DF11-B28D-001E68A993A4.root',
#        '/store/mc/Fall10/WtoMuNu_TuneD6T_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0013/A897F8F9-8ACF-DF11-A369-001E68A993EE.root',
#        '/store/mc/Fall10/WtoMuNu_TuneD6T_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0013/989A0F06-A3CE-DF11-B740-E0CB4E553688.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/711/F47056D9-8DE7-DF11-9355-003048F024FA.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/709/EC605D13-92E7-DF11-A30E-001617C3B79A.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/510/AAAB863A-2FE7-DF11-AF30-0030487CD16E.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/509/66337493-1DE7-DF11-8DF6-003048F024E0.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/461/52270E98-B3E6-DF11-BE88-0030487CD76A.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/459/7821C7A6-99E6-DF11-AE34-003048F11942.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/446/5ABFEAD7-6AE6-DF11-AE27-0030487CD7CA.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/442/F4CBBFB8-40E6-DF11-9FAD-0030487C2B86.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/442/E282E48F-4AE6-DF11-BDA4-00304879FA4A.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/442/B0642007-47E6-DF11-ACC0-0030487CD16E.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/442/A080AD3E-44E6-DF11-8C79-0030487C7828.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/442/98D8AFD9-42E6-DF11-BAC4-0030487CD6B4.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/442/8E396AD5-49E6-DF11-BDEF-000423D94908.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/442/46420456-46E6-DF11-9957-0030487CD716.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/442/2ED860B8-40E6-DF11-A060-0030487CD718.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/442/0C9D8607-47E6-DF11-996C-0030487C7392.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/441/54413500-18E6-DF11-A3F2-0030487C7392.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/392/0A91BA76-A3E5-DF11-B395-0030487CD17C.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/379/64F1872D-56E5-DF11-94D0-0030487CD6B4.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/378/04F1D08D-49E5-DF11-AF22-0030487CAEAC.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/375/963AC2C9-36E5-DF11-9B50-001D09F2426D.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/372/FCF91419-33E5-DF11-92EC-0030487CD7CA.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/304/220D14EB-DDE4-DF11-867F-001D09F24691.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/294/F40E37B8-F6E4-DF11-A54A-0030487CD906.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/294/D2EBA7B9-F6E4-DF11-9BCC-0030487CD6B4.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/294/BE272121-F8E4-DF11-B8B5-0030487CD17C.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/294/9EA0B5B2-F6E4-DF11-A1B7-000423D94494.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/294/9026ED9D-F7E4-DF11-ACFB-003048F11DE2.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/294/7C70DD9F-F7E4-DF11-B9DF-0016177CA7A0.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/294/66157B01-F6E4-DF11-B518-0030487CD6D2.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/294/44A6F1B1-F6E4-DF11-9228-000423D94908.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/294/1E47D920-F8E4-DF11-A64B-0030487CD7CA.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/293/289E56AF-C2E4-DF11-AF1F-001617E30E28.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/291/FC6AEC33-C6E4-DF11-862E-003048F024FE.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/291/F4FEDF51-E2E4-DF11-942D-001617C3B6CE.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/291/F2061C09-D5E4-DF11-8C3F-001D09F2532F.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/291/F08ACF28-D0E4-DF11-BBF8-001D09F23A34.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/291/E61D8328-D0E4-DF11-86AC-001D09F28F0C.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/291/E0BF91D0-D7E4-DF11-BB96-001D09F242EA.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/291/E0BDCCE2-D2E4-DF11-9047-003048CFB40C.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/291/DEF6D6AD-D5E4-DF11-8246-001D09F23174.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/291/DE1F5440-D9E4-DF11-9A83-001617E30E2C.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/291/DCD967F1-D9E4-DF11-AFB4-0030487C6090.root',
        '/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/291/CE8BB813-D7E4-DF11-9913-001D09F250AF.root',
    )
    #skipBadFiles = cms.untracked.bool(True) 
)
#process.source.inputCommands = cms.untracked.vstring(
#	"keep *", 
#	"drop *_MEtoEDMConverter_*_*", 
#	"drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__HLT"
#)

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

process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                     applyfilter = cms.untracked.bool(True),
                                     debugOn = cms.untracked.bool(False),
                                     numtrack = cms.untracked.uint32(10),
                                     thresh = cms.untracked.double(0.2)
                                     )

# Skimming step
process.load("HighMassAnalysis.Skimming.genLevelSequence_cff")
# Skim sequence
#process.load("HighMassAnalysis.Skimming.elecTauSkimSequence_cff")
process.load("HighMassAnalysis.Skimming.muTauSkimSequence_cff")
process.load("HighMassAnalysis.Skimming.triggerReq_cfi")

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
#                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID          = False
                 )

from PhysicsTools.PatAlgos.tools.electronTools import *
#addElectronUserIsolation(process,["Tracker"])
#addElectronUserIsolation(process,["Ecal"])
#addElectronUserIsolation(process,["Hcal"])

from PhysicsTools.PatAlgos.tools.muonTools import *
addMuonUserIsolation(process)

process.theSkim = cms.Sequence( 
#process.PFTauProducerSequence * 
#process.elecTauSkimSequence )
process.muTauSkimSequence )


if(signal):
	process.p = cms.Path(process.scrapingVeto + 
		      process.producePatTuple +
                      process.savePatTuple
	)
else:
	process.p = cms.Path(
                      process.scrapingVeto + 
                      process.muTauHLTFilter +
                      process.theSkim +
		      process.producePatTuple +
                      process.savePatTuple
	)
	

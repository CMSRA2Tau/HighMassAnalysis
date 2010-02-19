import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('hiMassTau')

# import of standard configurations for RECOnstruction
process.load('Configuration/StandardSequences/Services_cff')

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = 'INFO'
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

# standard sequences
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'MC_31X_V3::All'
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
    fileName = cms.untracked.string('test.root')
    #fileName = cms.untracked.string('CONDOR_OUTFILE')
    #fileName = cms.untracked.string("outputFILENAME")
)

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32( 100 )
)
#process.load("HighMassAnalysis.Configuration.FILESTOREAD")

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
       'dcache:/pnfs/cms/WAX/11/store/user/lpctau/HighMassTau/eluiggi/Ztautau/ZtautauSummer09MC31X7TeVZprimeTauTau_ETauSkim/96f5ca8867b4fcbe08570bd10fd38500/zttETauSkim_1.root',
    )
    #skipBadFiles = cms.untracked.bool(True) 
)

# import utility function for switching pat::Tau input to different reco::Tau collection stored on AOD
from PhysicsTools.PatAlgos.tools.tauTools import * 

# comment-out to take reco::CaloTaus instead of reco::PFTaus as input for pat::Tau production
#switchToCaloTau(process)

# comment-out to take shrinking dR = 5.0/Et(PFTau) signal cone instead of fixed dR = 0.07 signal cone reco::PFTaus
#switchToPFTauShrinkingCone(process)
switchToPFTauFixedCone(process)

# Use this to add the CaloTau side-by-side to the PAT PFTau (it will be called 'selectedLayer1CaloTaus')
# This doesn't take care of the trigger matching.  No iso deposits available yet
process.selectedLayer1CaloTaus = process.allLayer1Taus.clone(
    tauSource = cms.InputTag("caloRecoTauProducer"),
    addGenMatch = cms.bool(False),
    addGenJetMatch = cms.bool(False),
    #addTrigMatch = cms.bool(False),
    userIsolation    = cms.PSet(), # there is no path for calo tau isolation available at the moment
    isoDeposits  = cms.PSet(), # there is no path for calo tau isolation available at the moment
    tauIDSources = cms.PSet(  # all these are already present in 2.2.X AODSIM
            leadingTrackFinding = cms.InputTag("caloRecoTauDiscriminationByLeadingTrackFinding"),
            leadingTrackPtCut   = cms.InputTag("caloRecoTauDiscriminationByLeadingTrackPtCut"),
            byIsolation         = cms.InputTag("caloRecoTauDiscriminationByIsolation"),
            againstElectron     = cms.InputTag("caloRecoTauDiscriminationAgainstElectron"),  
    )
)

# Use this to add the tau collections with fixed signal cone of dR=0.15
process.selectedLayer1FixedConeHighEffPFTaus = process.allLayer1Taus.clone(
    tauSource = cms.InputTag("fixedConeHighEffPFTauProducer"),
    addGenMatch = cms.bool(True),
    genParticleMatch = cms.InputTag("tauMatchFixedConeHighEff"),
    addGenJetMatch = cms.bool(True),
    genJetMatch = cms.InputTag("tauGenJetMatchFixedConeHighEff"),
    #addTrigMatch = cms.bool(True),
    #trigPrimMatch = cms.VInputTag(cms.InputTag("tauTrigMatchHLT1TauFixedConeHighEff")),
    userIsolation    = cms.PSet(), # will not use
    isoDeposits  = cms.PSet(), # will not use
    tauIDSources = cms.PSet(  # all these are already present in 2.2.X AODSIM
        leadingTrackFinding = cms.InputTag("fixedConeHighEffPFTauDiscriminationByLeadingTrackFinding"),
        leadingTrackPtCut = cms.InputTag("fixedConeHighEffPFTauDiscriminationByLeadingTrackPtCut"),
        trackIsolation = cms.InputTag("fixedConeHighEffPFTauDiscriminationByTrackIsolation"),
        ecalIsolation = cms.InputTag("fixedConeHighEffPFTauDiscriminationByECALIsolation"),
        byIsolation = cms.InputTag("fixedConeHighEffPFTauDiscriminationByIsolation"),
        againstElectron = cms.InputTag("fixedConeHighEffPFTauDiscriminationAgainstElectron"),
        againstMuon = cms.InputTag("fixedConeHighEffPFTauDiscriminationAgainstMuon")
    )
)

# Use this to add the tau collections with fixed signal cone of dR=0.07
process.selectedLayer1FixedConePFTaus = process.allLayer1Taus.clone(
    tauSource = cms.InputTag("fixedConePFTauProducer"),
    addGenMatch = cms.bool(True),
    genParticleMatch = cms.InputTag("tauMatchFixedCone"),
    addGenJetMatch = cms.bool(True),
    genJetMatch = cms.InputTag("tauGenJetMatchFixedCone"),
    #addTrigMatch = cms.bool(True),
    #trigPrimMatch = cms.VInputTag(cms.InputTag("tauTrigMatchHLT1TauFixedCone")),
    userIsolation    = cms.PSet(), # will not use
    isoDeposits  = cms.PSet(), # will not use
    tauIDSources = cms.PSet(  # all these are already present in 2.2.X AODSIM
        leadingTrackFinding = cms.InputTag("fixedConePFTauDiscriminationByLeadingTrackFinding"),
        leadingTrackPtCut = cms.InputTag("fixedConePFTauDiscriminationByLeadingTrackPtCut"),
        trackIsolation = cms.InputTag("fixedConePFTauDiscriminationByTrackIsolation"),
        ecalIsolation = cms.InputTag("fixedConePFTauDiscriminationByECALIsolation"),
        byIsolation = cms.InputTag("fixedConePFTauDiscriminationByIsolation"),
        againstElectron = cms.InputTag("fixedConePFTauDiscriminationAgainstElectron"),
        againstMuon = cms.InputTag("fixedConePFTauDiscriminationAgainstMuon")
    )
)

# Use this to add the tau collections with shrinking signal cone of dR=5/ET
process.selectedLayer1ShrinkingConeHighEffPFTaus = process.allLayer1Taus.clone(
    tauSource = cms.InputTag("shrinkingConePFTauProducer"),
    addGenMatch = cms.bool(True),
    genParticleMatch = cms.InputTag("tauMatchShrinkingCone"),
    addGenJetMatch = cms.bool(True),
    genJetMatch = cms.InputTag("tauGenJetMatchShrinkingCone"),
    #addTrigMatch = cms.bool(True),
    #trigPrimMatch = cms.VInputTag(cms.InputTag("tauTrigMatchHLT1TauShrinkingCone")),
    userIsolation    = cms.PSet(), # will not use
    isoDeposits  = cms.PSet(), # will not use
    tauIDSources = cms.PSet(  # all these are already present in 2.2.X AODSIM
        leadingTrackFinding = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding"),
        leadingTrackPtCut = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackPtCut"),
        trackIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolation"),
        ecalIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolation"),
        byIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByIsolation"),
        againstElectron = cms.InputTag("shrinkingConePFTauDiscriminationAgainstElectron"),
        againstMuon = cms.InputTag("shrinkingConePFTauDiscriminationAgainstMuon")
    )
)

# Use this to add the tau collections with shrinking signal cone of dR=3/ET
process.selectedLayer1ShrinkingConePFTaus = process.allLayer1Taus.clone(
    tauSource = cms.InputTag("shrinkingTightConePFTauProducer"),
    addGenMatch = cms.bool(True),
    genParticleMatch = cms.InputTag("tauMatchShrinkingTightCone"),
    addGenJetMatch = cms.bool(True),
    genJetMatch = cms.InputTag("tauGenJetMatchShrinkingTightCone"),
    #addTrigMatch = cms.bool(True),
    #trigPrimMatch = cms.VInputTag(cms.InputTag("tauTrigMatchHLT1TauShrinkingTightCone")),
    userIsolation    = cms.PSet(), # will not use
    isoDeposits  = cms.PSet(), # will not use
    tauIDSources = cms.PSet(  # all these are already present in 2.2.X AODSIM
        leadingTrackFinding = cms.InputTag("shrinkingTightConePFTauDiscriminationByLeadingTrackFinding"),
        leadingTrackPtCut = cms.InputTag("shrinkingTightConePFTauDiscriminationByLeadingTrackPtCut"),
        trackIsolation = cms.InputTag("shrinkingTightConePFTauDiscriminationByTrackIsolation"),
        ecalIsolation = cms.InputTag("shrinkingTightConePFTauDiscriminationByECALIsolation"),
        byIsolation = cms.InputTag("shrinkingTightConePFTauDiscriminationByIsolation"),
        againstElectron = cms.InputTag("shrinkingTightConePFTauDiscriminationAgainstElectron"),
        againstMuon = cms.InputTag("shrinkingTightConePFTauDiscriminationAgainstMuon")
    )
)

# change the default pat sequence for taus - add creation of additional tau collections
process.allLayer1Objects.replace( process.allLayer1Taus, process.allLayer1Taus
                                                       + process.selectedLayer1CaloTaus
                                                       + process.selectedLayer1FixedConeHighEffPFTaus
                                                       + process.selectedLayer1FixedConePFTaus
                                                       + process.selectedLayer1ShrinkingConeHighEffPFTaus
                                                       + process.selectedLayer1ShrinkingConePFTaus
)

# change the default pat sequence for taus - add creation of additional tau collections
process.allLayer1Objects.replace( process.tauMatch, process.tauMatch
                                                       + process.tauMatchFixedConeHighEff
						       + process.tauMatchFixedCone
						       + process.tauMatchShrinkingCone
						       + process.tauMatchShrinkingTightCone
)

# change the default pat sequence for taus - add creation of additional tau collections
process.allLayer1Objects.replace( process.tauGenJetMatch, process.tauGenJetMatch
                                                       + process.tauGenJetMatchFixedConeHighEff
						       + process.tauGenJetMatchFixedCone
						       + process.tauGenJetMatchShrinkingCone
						       + process.tauGenJetMatchShrinkingTightCone
)

# import utility function for managing pat::METs
#from HighMassAnalysis.Configuration.metTools import *
# comment-out the addPFMet() function to add pfMET
#  - first Boolean switch on production of genMET with mu's 
#  - second Boolean switch on type-1 corrections
#addPFMet(process,False,False)
# comment-out to replace caloMET by pfMET in all di-tau objects
#replaceMETforDiTaus(process,
#                    cms.InputTag('layer1METs'),
#                    cms.InputTag('layer1PFMETs'))
## comment-out to add genMET with mu's to layer1MET (caloMET)  
#process.layer1METs.genMETSource = cms.InputTag('genMETWithMu')
## note: above line works when the first Boolean in the addPFMET is True, 
## otherwise you need comment-out the following:
### addGenMetWithMu(process) # comment-out only when the first Boolean in the addPFMET is False

# Modifications for 3X
from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')

# add a pf jet collection
#from PhysicsTools.PatAlgos.tools.jetTools import *
#addJetCollection(process, cms.InputTag('iterativeCone5PFJets'),  'PF',
#        doJTA=True,
#        doBTagging=False,
#        jetCorrLabel=None, # You may want to apply jet energy corrections
#        doType1MET=False)  # You don't want CaloMET with PFJets, do you?

# Modifications for 3X
from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process, cms.InputTag('iterativeCone5PFJets'),  'PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = None,
                 doType1MET   = False,
                 doL1Cleaning = False,                 
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("iterativeCone5GenJets"),
#                 doJetID          = False
                 )

# modify default pat sequence to include the matching sequences for the additional tau collections
#process.beforeLayer1Objects.replace( process.patTrigMatch, process.patTrigMatch
#                                                         + process.tauTrigMatchHLT1TauFixedConeHighEff
#                                                         + process.tauTrigMatchHLT1TauFixedCone
#                                                         + process.tauTrigMatchHLT1TauShrinkingCone
#                                                         + process.tauTrigMatchHLT1TauShrinkingTightCone)
#process.beforeLayer1Objects.replace( process.patMCTruth, process.patMCTruth
#                                                         + process.tauMatchFixedConeHighEff
#                                                         + process.tauGenJetMatchFixedConeHighEff
#                                                         + process.tauMatchFixedCone
#                                                         + process.tauGenJetMatchFixedCone
#                                                         + process.tauMatchShrinkingCone
#                                                         + process.tauGenJetMatchShrinkingCone
#                                                         + process.tauMatchShrinkingTightCone
#                                                         + process.tauGenJetMatchShrinkingTightCone)

#process.cteq65PdfWeights = cms.EDProducer("EwkPdfWeightProducer",
#      PdfInfoTag = cms.untracked.InputTag("genEventPdfInfo"),
#      PdfSetName = cms.untracked.string("cteq65.LHgrid")
#)
#process.MRST2006nnloPdfWeights = cms.EDProducer("EwkPdfWeightProducer",
#      PdfInfoTag = cms.untracked.InputTag("genEventPdfInfo"),
#      PdfSetName = cms.untracked.string("MRST2006nnlo.LHgrid")
#)
#process.MRST2007lomodPdfWeights = cms.EDProducer("EwkPdfWeightProducer",
#      PdfInfoTag = cms.untracked.InputTag("genEventPdfInfo"),
#      PdfSetName = cms.untracked.string("MRST2007lomod.LHgrid")
#)

process.p = cms.Path( process.producePatTuple
#                    + process.cteq65PdfWeights
#                    + process.MRST2006nnloPdfWeights
#                    + process.MRST2007lomodPdfWeights
                    + process.savePatTuple )

# print-out all python configuration parameter information
#print process.dumpPython()


import FWCore.ParameterSet.Config as cms

# import customized settings
from HighMassAnalysis.Configuration.customizePAT_cff import *

# produce collections of dR = 0.07 and dR = 0.15 fixed
# and dR = 5.0/Et shrinking signal cone taus using latest tags
from RecoTauTag.Configuration.RecoPFTauTag_cff import *
pfRecoTauTagInfoProducer.ChargedHadrCand_AssociationCone   = cms.double(1.0)

ic5PFJetTracksAssociatorAtVertex.j2tParametersVX = cms.PSet(
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(1.0)
)

ak5PFJetTracksAssociatorAtVertex.j2tParametersVX = cms.PSet(
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(1.0)
)

#--- tracker isolation
shrinkingConePFTauProducer.TrackerIsolConeSizeFormula   = cms.string('1.00') ## **
shrinkingConePFTauProducer.TrackerIsolConeSize_max   = cms.double(1.1)
#--- ecal isolation
shrinkingConePFTauProducer.ECALIsolConeSizeFormula      = cms.string('1.00') ## **
shrinkingConePFTauProducer.ECALIsolConeSize_max      = cms.double(1.1)
#--- hcal isolation
shrinkingConePFTauProducer.HCALIsolConeSizeFormula      = cms.string('1.00') ## **
shrinkingConePFTauProducer.HCALIsolConeSize_max      = cms.double(1.1)
#--- tracker signal
shrinkingConePFTauProducer.TrackerSignalConeSize_min = cms.double(0.07)
shrinkingConePFTauProducer.TrackerSignalConeSize_max = cms.double(0.15)
#--- ecal signal
shrinkingConePFTauProducer.ECALSignalConeSizeFormula    = cms.string('0.15') ## **
shrinkingConePFTauProducer.ECALSignalConeSize_min    = cms.double(0.00)
shrinkingConePFTauProducer.ECALSignalConeSize_max    = cms.double(0.15)
#--- hcal signal
shrinkingConePFTauProducer.HCALSignalConeSizeFormula    = cms.string('0.15') ## **
shrinkingConePFTauProducer.HCALSignalConeSize_min    = cms.double(0.00)
shrinkingConePFTauProducer.HCALSignalConeSize_max    = cms.double(0.15)

# gen met
from RecoMET.Configuration.GenMETParticles_cff import *
from RecoMET.METProducers.genMetTrue_cfi import *
#from RecoMET.METProducers.genMet_cfi import *
from RecoMET.METProducers.genMetCalo_cfi import *

# define sequence for gen jet production
from PhysicsTools.JetMCAlgos.TauGenJets_cfi import *

produceAndDiscriminatePFTausCustomized = cms.Sequence(
      ic5PFJetTracksAssociatorAtVertex *
      ak5PFJetTracksAssociatorAtVertex *
      recoTauAK5PFJets08Region*
      recoTauPileUpVertices*
      pfRecoTauTagInfoProducer *
      ak5PFJetsRecoTauPiZeros *
      ak5PFJetsLegacyTaNCPiZeros *
      shrinkingConePFTauProducer*
      shrinkingConePFTauDecayModeProducer*
      shrinkingConePFTauDecayModeIndexProducer*
      shrinkingConePFTauDiscriminationByLeadingTrackFinding*
      shrinkingConePFTauDiscriminationByLeadingTrackPtCut*
      shrinkingConePFTauDiscriminationByLeadingPionPtCut*
      shrinkingConePFTauDiscriminationByIsolation*
      shrinkingConePFTauDiscriminationByTrackIsolation*
      shrinkingConePFTauDiscriminationByECALIsolation*
      shrinkingConePFTauDiscriminationByIsolationUsingLeadingPion*
      shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion*
      shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion*
      shrinkingConePFTauDiscriminationAgainstElectron*
      shrinkingConePFTauDiscriminationAgainstMuon*
      ak5PFJetsLegacyHPSPiZeros *
      combinatoricRecoTaus *
      produceAndDiscriminateHPSPFTaus
)

produceThePatCands = cms.Sequence(
    recoElectronIsolation *
    eIdSequence *
#    ic5PFJetTracksAssociatorAtVertex *
#    ak5PFJetTracksAssociatorAtVertex *
#    recoTauAK5PFJets08Region*
#    recoTauPileUpVertices*
#    pfRecoTauTagInfoProducer *
#    ak5PFJetsRecoTauPiZeros *
#    ak5PFJetsLegacyTaNCPiZeros *
#    produceAndDiscriminateShrinkingConePFTausCustomized *
#    produceShrinkingConeDiscriminationByTauNeuralClassifier *
#    ak5PFJetsLegacyHPSPiZeros *
#    combinatoricRecoTaus *
#    produceAndDiscriminateHPSPFTaus *
    patCustomizedCandidates *
    selectedPatCustomizedCandidates 
    )

produceGenMETInfo = cms.Sequence(
    genMETParticles *
    genMetCalo *
    genMetTrue
)

if(data):
  producePatTuple = cms.Sequence(
    produceThePatCands
  )
else:
  producePatTuple = cms.Sequence(
    produceThePatCands
    + produceGenMETInfo
  )


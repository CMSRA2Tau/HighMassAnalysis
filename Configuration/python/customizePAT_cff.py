import FWCore.ParameterSet.Config as cms
import copy

# Standard pat sequences
from PhysicsTools.PatAlgos.patSequences_cff import * 

# --------------------Modifications for taus--------------------

shrinkingConeTauIsoDepositPFCandidates = copy.deepcopy(tauIsoDepositPFCandidates)
shrinkingConeTauIsoDepositPFChargedHadrons = copy.deepcopy(tauIsoDepositPFChargedHadrons)
shrinkingConeTauIsoDepositPFNeutralHadrons = copy.deepcopy(tauIsoDepositPFNeutralHadrons)
shrinkingConeTauIsoDepositPFGammas = copy.deepcopy(tauIsoDepositPFGammas)
fixedConeHighEffTauIsoDepositPFCandidates = copy.deepcopy(tauIsoDepositPFCandidates)
fixedConeHighEffTauIsoDepositPFChargedHadrons = copy.deepcopy(tauIsoDepositPFChargedHadrons)
fixedConeHighEffTauIsoDepositPFNeutralHadrons = copy.deepcopy(tauIsoDepositPFNeutralHadrons)
fixedConeHighEffTauIsoDepositPFGammas = copy.deepcopy(tauIsoDepositPFGammas)
fixedConeHighEffTauIsoDepositPFCandidates.src = cms.InputTag("fixedConeHighEffPFTauProducer")
fixedConeHighEffTauIsoDepositPFChargedHadrons.src = cms.InputTag("fixedConeHighEffPFTauProducer")
fixedConeHighEffTauIsoDepositPFNeutralHadrons.src = cms.InputTag("fixedConeHighEffPFTauProducer")
fixedConeHighEffTauIsoDepositPFGammas.src = cms.InputTag("fixedConeHighEffPFTauProducer")
fixedConeHighEffTauIsoDepositPFCandidates.ExtractorPSet.tauSource = cms.InputTag("fixedConeHighEffPFTauProducer")
fixedConeHighEffTauIsoDepositPFChargedHadrons.ExtractorPSet.tauSource = cms.InputTag("fixedConeHighEffPFTauProducer")
fixedConeHighEffTauIsoDepositPFNeutralHadrons.ExtractorPSet.tauSource = cms.InputTag("fixedConeHighEffPFTauProducer")
fixedConeHighEffTauIsoDepositPFGammas.ExtractorPSet.tauSource = cms.InputTag("fixedConeHighEffPFTauProducer")
fixedConeTauIsoDepositPFCandidates = copy.deepcopy(tauIsoDepositPFCandidates)
fixedConeTauIsoDepositPFChargedHadrons = copy.deepcopy(tauIsoDepositPFChargedHadrons)
fixedConeTauIsoDepositPFNeutralHadrons = copy.deepcopy(tauIsoDepositPFNeutralHadrons)
fixedConeTauIsoDepositPFGammas = copy.deepcopy(tauIsoDepositPFGammas)
fixedConeTauIsoDepositPFCandidates.src = cms.InputTag("fixedConePFTauProducer")
fixedConeTauIsoDepositPFChargedHadrons.src = cms.InputTag("fixedConePFTauProducer")
fixedConeTauIsoDepositPFNeutralHadrons.src = cms.InputTag("fixedConePFTauProducer")
fixedConeTauIsoDepositPFGammas.src = cms.InputTag("fixedConePFTauProducer")
fixedConeTauIsoDepositPFCandidates.ExtractorPSet.tauSource = cms.InputTag("fixedConePFTauProducer")
fixedConeTauIsoDepositPFChargedHadrons.ExtractorPSet.tauSource = cms.InputTag("fixedConePFTauProducer")
fixedConeTauIsoDepositPFNeutralHadrons.ExtractorPSet.tauSource = cms.InputTag("fixedConePFTauProducer")
fixedConeTauIsoDepositPFGammas.ExtractorPSet.tauSource = cms.InputTag("fixedConePFTauProducer")
shrinkingTightConeTauIsoDepositPFCandidates = copy.deepcopy(tauIsoDepositPFCandidates)
shrinkingTightConeTauIsoDepositPFChargedHadrons = copy.deepcopy(tauIsoDepositPFChargedHadrons)
shrinkingTightConeTauIsoDepositPFNeutralHadrons = copy.deepcopy(tauIsoDepositPFNeutralHadrons)
shrinkingTightConeTauIsoDepositPFGammas = copy.deepcopy(tauIsoDepositPFGammas)
shrinkingTightConeTauIsoDepositPFCandidates.src = cms.InputTag("shrinkingTightConePFTauProducer")
shrinkingTightConeTauIsoDepositPFChargedHadrons.src = cms.InputTag("shrinkingTightConePFTauProducer")
shrinkingTightConeTauIsoDepositPFNeutralHadrons.src = cms.InputTag("shrinkingTightConePFTauProducer")
shrinkingTightConeTauIsoDepositPFGammas.src = cms.InputTag("shrinkingTightConePFTauProducer")
shrinkingTightConeTauIsoDepositPFCandidates.ExtractorPSet.tauSource = cms.InputTag("shrinkingTightConePFTauProducer")
shrinkingTightConeTauIsoDepositPFChargedHadrons.ExtractorPSet.tauSource = cms.InputTag("shrinkingTightConePFTauProducer")
shrinkingTightConeTauIsoDepositPFNeutralHadrons.ExtractorPSet.tauSource = cms.InputTag("shrinkingTightConePFTauProducer")
shrinkingTightConeTauIsoDepositPFGammas.ExtractorPSet.tauSource = cms.InputTag("shrinkingTightConePFTauProducer")

patCustomizedPFTauIsolation = cms.Sequence(  shrinkingConeTauIsoDepositPFCandidates
                                 * shrinkingConeTauIsoDepositPFChargedHadrons
                                 * shrinkingConeTauIsoDepositPFNeutralHadrons
                                 * shrinkingConeTauIsoDepositPFGammas
                                 * shrinkingTightConeTauIsoDepositPFCandidates
                                 * shrinkingTightConeTauIsoDepositPFChargedHadrons
                                 * shrinkingTightConeTauIsoDepositPFNeutralHadrons
                                 * shrinkingTightConeTauIsoDepositPFGammas
                                 * fixedConeTauIsoDepositPFCandidates
                                 * fixedConeTauIsoDepositPFChargedHadrons
                                 * fixedConeTauIsoDepositPFNeutralHadrons
                                 * fixedConeTauIsoDepositPFGammas                                   
                                 * fixedConeHighEffTauIsoDepositPFCandidates
                                 * fixedConeHighEffTauIsoDepositPFChargedHadrons
                                 * fixedConeHighEffTauIsoDepositPFNeutralHadrons
                                 * fixedConeHighEffTauIsoDepositPFGammas
                                )

shrinkingConeTauMatch = copy.deepcopy(tauMatch)
shrinkingConeTauGenJetMatch = copy.deepcopy(tauGenJetMatch)
fixedConeHighEffTauMatch = copy.deepcopy(tauMatch)
fixedConeHighEffTauMatch.src         = cms.InputTag("fixedConeHighEffPFTauProducer")
fixedConeHighEffTauGenJetMatch = copy.deepcopy(tauGenJetMatch)
fixedConeHighEffTauGenJetMatch.src         = cms.InputTag("fixedConeHighEffPFTauProducer")
fixedConeTauMatch = copy.deepcopy(tauMatch)
fixedConeTauMatch.src         = cms.InputTag("fixedConePFTauProducer")
fixedConeTauGenJetMatch = copy.deepcopy(tauGenJetMatch)
fixedConeTauGenJetMatch.src         = cms.InputTag("fixedConePFTauProducer")
shrinkingTightConeTauMatch = copy.deepcopy(tauMatch)
shrinkingTightConeTauMatch.src         = cms.InputTag("shrinkingTightConePFTauProducer")
shrinkingTightConeTauGenJetMatch = copy.deepcopy(tauGenJetMatch)
shrinkingTightConeTauGenJetMatch.src         = cms.InputTag("shrinkingTightConePFTauProducer")

shrinkingConePatTaus = copy.deepcopy(patTaus)
fixedConeHighEffPatTaus = copy.deepcopy(patTaus)
fixedConePatTaus = copy.deepcopy(patTaus)
shrinkingTightConePatTaus = copy.deepcopy(patTaus)

shrinkingConePatTaus.tauSource = cms.InputTag("shrinkingConePFTauProducer")
fixedConeHighEffPatTaus.tauSource = cms.InputTag("fixedConeHighEffPFTauProducer")
fixedConePatTaus.tauSource = cms.InputTag("fixedConePFTauProducer")
shrinkingTightConePatTaus.tauSource = cms.InputTag("shrinkingTightConePFTauProducer")

shrinkingConePatTaus.isoDeposits.pfAllParticles = cms.InputTag("shrinkingConeTauIsoDepositPFCandidates")
fixedConeHighEffPatTaus.isoDeposits.pfAllParticles = cms.InputTag("fixedConeHighEffTauIsoDepositPFCandidates")
fixedConePatTaus.isoDeposits.pfAllParticles = cms.InputTag("fixedConeTauIsoDepositPFCandidates")
shrinkingTightConePatTaus.isoDeposits.pfAllParticles = cms.InputTag("shrinkingTightConeTauIsoDepositPFCandidates")

shrinkingConePatTaus.isoDeposits.pfChargedHadron = cms.InputTag("shrinkingConeTauIsoDepositPFChargedHadrons")
fixedConeHighEffPatTaus.isoDeposits.pfChargedHadron = cms.InputTag("fixedConeHighEffTauIsoDepositPFChargedHadrons")
fixedConePatTaus.isoDeposits.pfChargedHadron = cms.InputTag("fixedConeTauIsoDepositPFChargedHadrons")
shrinkingTightConePatTaus.isoDeposits.pfChargedHadron = cms.InputTag("shrinkingTightConeTauIsoDepositPFChargedHadrons")

shrinkingConePatTaus.isoDeposits.pfNeutralHadron = cms.InputTag("shrinkingConeTauIsoDepositPFNeutralHadrons")
fixedConeHighEffPatTaus.isoDeposits.pfNeutralHadron = cms.InputTag("fixedConeHighEffTauIsoDepositPFNeutralHadrons")
fixedConePatTaus.isoDeposits.pfNeutralHadron = cms.InputTag("fixedConeTauIsoDepositPFNeutralHadrons")
shrinkingTightConePatTaus.isoDeposits.pfNeutralHadron = cms.InputTag("shrinkingTightConeTauIsoDepositPFNeutralHadrons")

shrinkingConePatTaus.isoDeposits.pfGamma = cms.InputTag("shrinkingConeTauIsoDepositPFGammas")
fixedConeHighEffPatTaus.isoDeposits.pfGamma = cms.InputTag("fixedConeHighEffTauIsoDepositPFGammas")
fixedConePatTaus.isoDeposits.pfGamma = cms.InputTag("fixedConeTauIsoDepositPFGammas")
shrinkingTightConePatTaus.isoDeposits.pfGamma = cms.InputTag("shrinkingTightConeTauIsoDepositPFGammas")

shrinkingConePatTaus.userIsolation.pfAllParticles.src = cms.InputTag("shrinkingConeTauIsoDepositPFCandidates")
fixedConeHighEffPatTaus.userIsolation.pfAllParticles.src = cms.InputTag("fixedConeHighEffTauIsoDepositPFCandidates")
fixedConePatTaus.userIsolation.pfAllParticles.src = cms.InputTag("fixedConeTauIsoDepositPFCandidates")
shrinkingTightConePatTaus.userIsolation.pfAllParticles.src = cms.InputTag("shrinkingTightConeTauIsoDepositPFCandidates")

shrinkingConePatTaus.userIsolation.pfChargedHadron.src = cms.InputTag("shrinkingConeTauIsoDepositPFChargedHadrons")
fixedConeHighEffPatTaus.userIsolation.pfChargedHadron.src = cms.InputTag("fixedConeHighEffTauIsoDepositPFChargedHadrons")
fixedConePatTaus.userIsolation.pfChargedHadron.src = cms.InputTag("fixedConeTauIsoDepositPFChargedHadrons")
shrinkingTightConePatTaus.userIsolation.pfChargedHadron.src = cms.InputTag("shrinkingTightConeTauIsoDepositPFChargedHadrons")

shrinkingConePatTaus.userIsolation.pfNeutralHadron.src = cms.InputTag("shrinkingConeTauIsoDepositPFNeutralHadrons")
fixedConeHighEffPatTaus.userIsolation.pfNeutralHadron.src = cms.InputTag("fixedConeHighEffTauIsoDepositPFNeutralHadrons")
fixedConePatTaus.userIsolation.pfNeutralHadron.src = cms.InputTag("fixedConeTauIsoDepositPFNeutralHadrons")
shrinkingTightConePatTaus.userIsolation.pfNeutralHadron.src = cms.InputTag("shrinkingTightConeTauIsoDepositPFNeutralHadrons")

shrinkingConePatTaus.userIsolation.pfGamma.src = cms.InputTag("shrinkingConeTauIsoDepositPFGammas")
fixedConeHighEffPatTaus.userIsolation.pfGamma.src = cms.InputTag("fixedConeHighEffTauIsoDepositPFGammas")
fixedConePatTaus.userIsolation.pfGamma.src = cms.InputTag("fixedConeTauIsoDepositPFGammas")
shrinkingTightConePatTaus.userIsolation.pfGamma.src = cms.InputTag("shrinkingTightConeTauIsoDepositPFGammas")

fixedConeHighEffPatTaus.tauIDSources = cms.PSet(
        leadingTrackFinding = cms.InputTag("fixedConeHighEffPFTauDiscriminationByLeadingTrackFinding"),
        leadingTrackPtCut = cms.InputTag("fixedConeHighEffPFTauDiscriminationByLeadingTrackPtCut"),
        leadingPionPtCut = cms.InputTag("fixedConeHighEffPFTauDiscriminationByLeadingPionPtCut"),
        trackIsolation = cms.InputTag("fixedConeHighEffPFTauDiscriminationByTrackIsolation"),
        trackIsolationUsingLeadingPion = cms.InputTag("fixedConeHighEffPFTauDiscriminationByTrackIsolationUsingLeadingPion"),
        ecalIsolation = cms.InputTag("fixedConeHighEffPFTauDiscriminationByECALIsolation"),
        ecalIsolationUsingLeadingPion = cms.InputTag("fixedConeHighEffPFTauDiscriminationByECALIsolationUsingLeadingPion"),
        byIsolation = cms.InputTag("fixedConeHighEffPFTauDiscriminationByIsolation"),
        byIsolationUsingLeadingPion = cms.InputTag("fixedConeHighEffPFTauDiscriminationByIsolationUsingLeadingPion"),
        againstElectron = cms.InputTag("fixedConeHighEffPFTauDiscriminationAgainstElectron"),
        againstMuon = cms.InputTag("fixedConeHighEffPFTauDiscriminationAgainstMuon"),
)
fixedConePatTaus.tauIDSources = cms.PSet(
        leadingTrackFinding = cms.InputTag("fixedConePFTauDiscriminationByLeadingTrackFinding"),
        leadingTrackPtCut = cms.InputTag("fixedConePFTauDiscriminationByLeadingTrackPtCut"),
        leadingPionPtCut = cms.InputTag("fixedConePFTauDiscriminationByLeadingPionPtCut"),
        trackIsolation = cms.InputTag("fixedConePFTauDiscriminationByTrackIsolation"),
        trackIsolationUsingLeadingPion = cms.InputTag("fixedConePFTauDiscriminationByTrackIsolationUsingLeadingPion"),
        ecalIsolation = cms.InputTag("fixedConePFTauDiscriminationByECALIsolation"),
        ecalIsolationUsingLeadingPion = cms.InputTag("fixedConePFTauDiscriminationByECALIsolationUsingLeadingPion"),
        byIsolation = cms.InputTag("fixedConePFTauDiscriminationByIsolation"),
        byIsolationUsingLeadingPion = cms.InputTag("fixedConePFTauDiscriminationByIsolationUsingLeadingPion"),
        againstElectron = cms.InputTag("fixedConePFTauDiscriminationAgainstElectron"),
        againstMuon = cms.InputTag("fixedConePFTauDiscriminationAgainstMuon"),
)
shrinkingTightConePatTaus.tauIDSources = cms.PSet(
        leadingTrackFinding = cms.InputTag("shrinkingTightConePFTauDiscriminationByLeadingTrackFinding"),
        leadingTrackPtCut = cms.InputTag("shrinkingTightConePFTauDiscriminationByLeadingTrackPtCut"),
        leadingPionPtCut = cms.InputTag("shrinkingTightConePFTauDiscriminationByLeadingPionPtCut"),
        trackIsolation = cms.InputTag("shrinkingTightConePFTauDiscriminationByTrackIsolation"),
        trackIsolationUsingLeadingPion = cms.InputTag("shrinkingTightConePFTauDiscriminationByTrackIsolationUsingLeadingPion"),
        ecalIsolation = cms.InputTag("shrinkingTightConePFTauDiscriminationByECALIsolation"),
        ecalIsolationUsingLeadingPion = cms.InputTag("shrinkingTightConePFTauDiscriminationByECALIsolationUsingLeadingPion"),
        byIsolation = cms.InputTag("shrinkingTightConePFTauDiscriminationByIsolation"),
        byIsolationUsingLeadingPion = cms.InputTag("shrinkingTightConePFTauDiscriminationByIsolationUsingLeadingPion"),
        againstElectron = cms.InputTag("shrinkingTightConePFTauDiscriminationAgainstElectron"),
        againstMuon = cms.InputTag("shrinkingTightConePFTauDiscriminationAgainstMuon"),
)
shrinkingConePatTaus.tauIDSources = cms.PSet(
        leadingTrackFinding = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding"),
        leadingTrackPtCut = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackPtCut"),
        leadingPionPtCut = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingPionPtCut"),
        trackIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolation"),
        trackIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion"),
        ecalIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolation"),
        ecalIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion"),
        byIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByIsolation"),
        byIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByIsolationUsingLeadingPion"),
        againstElectron = cms.InputTag("shrinkingConePFTauDiscriminationAgainstElectron"),
        againstMuon = cms.InputTag("shrinkingConePFTauDiscriminationAgainstMuon"),
#        byTaNC = cms.InputTag("shrinkingConePFTauDiscriminationByTaNC"),
#        byTaNCfrOnePercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrOnePercent"),
#        byTaNCfrHalfPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrHalfPercent"),
#        byTaNCfrQuarterPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrQuarterPercent"),
#        byTaNCfrTenthPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrTenthPercent")
)

fixedConeHighEffPatTaus.addDecayMode = cms.bool(False)
fixedConePatTaus.addDecayMode = cms.bool(False)
shrinkingTightConePatTaus.addDecayMode = cms.bool(False)
shrinkingConePatTaus.addDecayMode = cms.bool(False)

fixedConeHighEffPatTaus.decayModeSrc = cms.InputTag("fixedConeHighEffPFTauDecayModeProducer")
fixedConePatTaus.decayModeSrc = cms.InputTag("fixedConePFTauDecayModeProducer")
shrinkingTightConePatTaus.decayModeSrc = cms.InputTag("shrinkingTightConePFTauDecayModeProducer")
shrinkingConePatTaus.decayModeSrc = cms.InputTag("shrinkingConePFTauDecayModeProducer")

fixedConeHighEffPatTaus.addGenMatch = cms.bool(False)
fixedConePatTaus.addGenMatch = cms.bool(False)
shrinkingTightConePatTaus.addGenMatch = cms.bool(False)
shrinkingConePatTaus.addGenMatch = cms.bool(False)

fixedConeHighEffPatTaus.embedGenMatch = cms.bool(False)
fixedConePatTaus.embedGenMatch = cms.bool(False)
shrinkingTightConePatTaus.embedGenMatch = cms.bool(False)
shrinkingConePatTaus.embedGenMatch = cms.bool(False)

shrinkingConePatTaus.genParticleMatch = cms.InputTag("shrinkingConeTauMatch")
fixedConeHighEffPatTaus.genParticleMatch = cms.InputTag("fixedConeHighEffTauMatch")
fixedConePatTaus.genParticleMatch = cms.InputTag("fixedConeTauMatch")
shrinkingTightConePatTaus.genParticleMatch = cms.InputTag("shrinkingTightConeTauMatch")

fixedConeHighEffPatTaus.addGenJetMatch = cms.bool(False)
fixedConePatTaus.addGenJetMatch = cms.bool(False)
shrinkingTightConePatTaus.addGenJetMatch = cms.bool(False)
shrinkingConePatTaus.addGenJetMatch = cms.bool(False)

fixedConeHighEffPatTaus.embedGenJetMatch = cms.bool(False)
fixedConePatTaus.embedGenJetMatch = cms.bool(False)
shrinkingTightConePatTaus.embedGenJetMatch = cms.bool(False)
shrinkingConePatTaus.embedGenJetMatch = cms.bool(False)

shrinkingConePatTaus.genJetMatch = cms.InputTag("shrinkingConeTauGenJetMatch")
fixedConeHighEffPatTaus.genJetMatch = cms.InputTag("fixedConeHighEffTauGenJetMatch")
fixedConePatTaus.genJetMatch = cms.InputTag("fixedConeTauGenJetMatch")
shrinkingTightConePatTaus.genJetMatch = cms.InputTag("shrinkingTightConeTauGenJetMatch")

selectedLayer1ShrinkingConeHighEffPFTaus = copy.deepcopy(selectedPatTaus)
selectedLayer1ShrinkingConeHighEffPFTaus.src = cms.InputTag("shrinkingConePatTaus")
selectedLayer1ShrinkingConePFTaus = copy.deepcopy(selectedPatTaus)
selectedLayer1ShrinkingConePFTaus.src = cms.InputTag("shrinkingTightConePatTaus")
selectedLayer1FixedConePFTaus = copy.deepcopy(selectedPatTaus)
selectedLayer1FixedConePFTaus.src = cms.InputTag("fixedConePatTaus")
selectedLayer1FixedConeHighEffPFTaus = copy.deepcopy(selectedPatTaus)
selectedLayer1FixedConeHighEffPFTaus.src = cms.InputTag("fixedConeHighEffPatTaus")

selectedLayer1CaloTaus = copy.deepcopy(patTaus)
selectedLayer1CaloTaus.tauSource = cms.InputTag("caloRecoTauProducer")
selectedLayer1CaloTaus.isoDeposits  = cms.PSet() # there is no path for calo tau isolation available at the moment
selectedLayer1CaloTaus.userIsolation    = cms.PSet() # there is no path for calo tau isolation available at the moment
selectedLayer1CaloTaus.tauIDSources = cms.PSet(
        leadingTrackFinding = cms.InputTag("caloRecoTauDiscriminationByLeadingTrackFinding"),
        leadingTrackPtCut   = cms.InputTag("caloRecoTauDiscriminationByLeadingTrackPtCut"),
        byIsolation         = cms.InputTag("caloRecoTauDiscriminationByIsolation"),
        againstElectron     = cms.InputTag("caloRecoTauDiscriminationAgainstElectron"),
)
selectedLayer1CaloTaus.addDecayMode = cms.bool(False)
selectedLayer1CaloTaus.addGenMatch = cms.bool(False)
selectedLayer1CaloTaus.embedGenMatch = cms.bool(False)
selectedLayer1CaloTaus.addGenJetMatch = cms.bool(False)
selectedLayer1CaloTaus.addGenJetMatch = cms.bool(False)
#selectedLayer1CaloTaus.addTrigMatch = cms.bool(False)

makeCustomizedPatTaus = cms.Sequence(
    patPFCandidateIsoDepositSelection *
    patCustomizedPFTauIsolation *
#    shrinkingConeTauMatch *
#    fixedConeHighEffTauMatch *
#    fixedConeTauMatch *
#    shrinkingTightConeTauMatch *    
#    tauGenJets *
#    tauGenJetsSelectorAllHadrons *
#    shrinkingConeTauGenJetMatch *
#    fixedConeHighEffTauGenJetMatch *
#    fixedConeTauGenJetMatch *
#    shrinkingTightConeTauGenJetMatch *
    # object production
    fixedConeHighEffPatTaus *
    fixedConePatTaus *
    shrinkingTightConePatTaus *
    shrinkingConePatTaus *
    selectedLayer1CaloTaus
    )

# --------------------Modifications for muons--------------------

patMuons.userIsolation = cms.PSet(
  hcal = cms.PSet(
    src = cms.InputTag("muIsoDepositCalByAssociatorTowers","hcal"),
    deltaR = cms.double(0.6)
  ),
  tracker = cms.PSet(
    src = cms.InputTag("muIsoDepositTk"),
    deltaR = cms.double(0.6)
  ),
  user = cms.VPSet(
    cms.PSet(
      src = cms.InputTag("muIsoDepositCalByAssociatorTowers","ho"),
      deltaR = cms.double(0.6)
    ), 
    cms.PSet(
      src = cms.InputTag("muIsoDepositJets"),
      deltaR = cms.double(0.6)
    )
  ),
  ecal = cms.PSet(
    src = cms.InputTag("muIsoDepositCalByAssociatorTowers","ecal"),
    deltaR = cms.double(0.6)
  )
)
patMuons.isoDeposits = cms.PSet(
  tracker         = patMuons.userIsolation.tracker.src,
  ecal            = patMuons.userIsolation.ecal.src,
  hcal            = patMuons.userIsolation.hcal.src,
)
patMuons.addGenMatch = cms.bool(False)
patMuons.embedGenMatch = cms.bool(False)

# --------------------Modifications for electrons--------------------

from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositHcalFromTowers_cff import *
electronIsoDeposits = cms.Sequence( eleIsoDepositTk
                                   *eleIsoDepositEcalFromHits 
                                   *eleIsoDepositHcalFromTowers 
				   *eleIsoDepositHcalDepth1FromTowers
				   *eleIsoDepositHcalDepth2FromTowers
				   )
recoElectronIsolation = cms.Sequence( electronIsoDeposits )

patElectrons.addGenMatch = cms.bool(False)
patElectrons.embedGenMatch = cms.bool(False)
patElectrons.genParticleMatch = cms.InputTag("electronMatch")

cleanPatElectrons.checkOverlaps = cms.PSet()

from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
heepPatElectrons = cms.EDProducer("HEEPAttStatusToPAT",
                                          eleLabel = cms.InputTag("selectedPatElectrons"),
                                          barrelCuts = cms.PSet(heepBarrelCuts),
                                          endcapCuts = cms.PSet(heepEndcapCuts)
                                          )

# --------------------Modifications for jets--------------------

patJets.addGenPartonMatch   = cms.bool(False)                           ## switch on/off matching to quarks from hard scatterin
patJets.embedGenPartonMatch = cms.bool(False)                           ## switch on/off embedding of the GenParticle parton for this jet
patJets.addGenJetMatch      = cms.bool(False)                           ## switch on/off matching to GenJet's
patJets.embedGenJetMatch    = cms.bool(False)                           ## switch on/off embedding of matched genJet's
patJets.addPartonJetMatch   = cms.bool(False)                          ## switch on/off matching to PartonJet's (not implemented yet)
patJets.getJetMCFlavour     = cms.bool(False)


cleanPatJets.checkOverlaps = cms.PSet()

makeCustomizedPatJets = cms.Sequence(
    patJetCorrections *
    patJetCharge *
    patJets
)


# --------------------Modifications for photons--------------------

patPhotons.addGenMatch = cms.bool(False)
patPhotons.embedGenMatch = cms.bool(False)

# --------------------Modifications for MET--------------------

patMETs.addGenMET    = cms.bool(False)
patMETs.metSource  = cms.InputTag("met")

# --------------------Modified PAT sequences--------------------

patCustomizedCandidates = cms.Sequence(
    patElectrons +
    patMuons     +
    makeCustomizedPatTaus      +
    patPhotons   +
    makeCustomizedPatJets      +
    patMETs )

selectedPatCustomizedCandidates = cms.Sequence(
    selectedPatElectrons +
    heepPatElectrons	 +
    selectedPatMuons     +
    selectedLayer1ShrinkingConeHighEffPFTaus +
    selectedLayer1ShrinkingConePFTaus +
    selectedLayer1FixedConePFTaus +
    selectedLayer1FixedConeHighEffPFTaus +    
    selectedPatPhotons   +
    selectedPatJets )

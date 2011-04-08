import FWCore.ParameterSet.Config as cms
import copy

from HighMassAnalysis.Configuration.hiMassSetup_cfi import *
if(data):
  doGen = False
else:
  doGen = True

# Standard pat sequences
from PhysicsTools.PatAlgos.patSequences_cff import * 

# --------------------Modifications for taus--------------------

shrinkingConeTauIsoDepositPFCandidates = copy.deepcopy(tauIsoDepositPFCandidates)
shrinkingConeTauIsoDepositPFChargedHadrons = copy.deepcopy(tauIsoDepositPFChargedHadrons)
shrinkingConeTauIsoDepositPFNeutralHadrons = copy.deepcopy(tauIsoDepositPFNeutralHadrons)
shrinkingConeTauIsoDepositPFGammas = copy.deepcopy(tauIsoDepositPFGammas)

patCustomizedPFTauIsolation = cms.Sequence(  shrinkingConeTauIsoDepositPFCandidates
                                 * shrinkingConeTauIsoDepositPFChargedHadrons
                                 * shrinkingConeTauIsoDepositPFNeutralHadrons
                                 * shrinkingConeTauIsoDepositPFGammas
                                )

shrinkingConeTauMatch = copy.deepcopy(tauMatch)
shrinkingConeTauGenJetMatch = copy.deepcopy(tauGenJetMatch)

shrinkingConePatTaus = copy.deepcopy(patTaus)

shrinkingConePatTaus.tauSource = cms.InputTag("shrinkingConePFTauProducer")

shrinkingConePatTaus.isoDeposits = cms.PSet()

shrinkingConePatTaus.userIsolation = cms.PSet()

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

shrinkingConePatTaus.addDecayMode = cms.bool(doGen)

shrinkingConePatTaus.decayModeSrc = cms.InputTag("shrinkingConePFTauDecayModeProducer")

shrinkingConePatTaus.addGenMatch = cms.bool(doGen)

shrinkingConePatTaus.embedGenMatch = cms.bool(doGen)

shrinkingConePatTaus.genParticleMatch = cms.InputTag("shrinkingConeTauMatch")

shrinkingConePatTaus.addGenJetMatch = cms.bool(doGen)

shrinkingConePatTaus.embedGenJetMatch = cms.bool(doGen)

shrinkingConePatTaus.genJetMatch = cms.InputTag("shrinkingConeTauGenJetMatch")

selectedLayer1ShrinkingConeHighEffPFTaus = copy.deepcopy(selectedPatTaus)
selectedLayer1ShrinkingConeHighEffPFTaus.src = cms.InputTag("shrinkingConePatTaus")

makeCustomizedPatTauIso = cms.Sequence(
    patPFCandidateIsoDepositSelection
)

makeCustomizedPatTauMatch = cms.Sequence(
    shrinkingConeTauMatch *
    tauGenJets *
    tauGenJetsSelectorAllHadrons *
    shrinkingConeTauGenJetMatch
) 

makePatTaus = cms.Sequence(
    shrinkingConePatTaus
)

if(data):
  makeCustomizedPatTaus = cms.Sequence(
    makeCustomizedPatTauIso
    + makePatTaus
  )
else:
  makeCustomizedPatTaus = cms.Sequence(
    makeCustomizedPatTauIso
    + makeCustomizedPatTauMatch
    + makePatTaus
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
patMuons.addGenMatch = cms.bool(doGen)
patMuons.embedGenMatch = cms.bool(doGen)

# --------------------Modifications for electrons--------------------

from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositHcalFromTowers_cff import *
electronIsoDeposits = cms.Sequence( eleIsoDepositTk
                                   *eleIsoDepositEcalFromHits 
                                   *eleIsoDepositHcalFromTowers 
				   *eleIsoDepositHcalDepth1FromTowers
				   *eleIsoDepositHcalDepth2FromTowers
				   )
recoElectronIsolation = cms.Sequence( electronIsoDeposits )

patElectrons.addGenMatch = cms.bool(doGen)
patElectrons.embedGenMatch = cms.bool(doGen)
patElectrons.genParticleMatch = cms.InputTag("electronMatch")

cleanPatElectrons.checkOverlaps = cms.PSet()

from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
heepPatElectrons = cms.EDProducer("HEEPAttStatusToPAT",
                                          eleLabel = cms.InputTag("selectedPatElectrons"),
                                          barrelCuts = cms.PSet(heepBarrelCuts),
                                          endcapCuts = cms.PSet(heepEndcapCuts)
                                          )

# --------------------Modifications for jets--------------------

patJets.addJetCorrFactors    = cms.bool(False)
patJets.addGenPartonMatch   = cms.bool(doGen)                           ## switch on/off matching to quarks from hard scatterin
patJets.embedGenPartonMatch = cms.bool(doGen)                           ## switch on/off embedding of the GenParticle parton for this jet
patJets.addGenJetMatch      = cms.bool(doGen)                           ## switch on/off matching to GenJet's
patJets.embedGenJetMatch    = cms.bool(doGen)                           ## switch on/off embedding of matched genJet's
patJets.addPartonJetMatch   = cms.bool(doGen)                          ## switch on/off matching to PartonJet's (not implemented yet)
patJets.getJetMCFlavour     = cms.bool(doGen)

cleanPatJets.checkOverlaps = cms.PSet()

patJetsPF = copy.deepcopy(patJets)
patJetsPF.jetSource = cms.InputTag("ak5PFJets")
patJetsPF.addJetCorrFactors    = cms.bool(False)
patJetsPF.addJetID = cms.bool(False)
patJetsPF.addBTagInfo          = cms.bool(False)
patJetsPF.addDiscriminators    = cms.bool(False)
patJetsPF.addAssociatedTracks    = cms.bool(False)
patJetsPF.addJetCharge    = cms.bool(False)
patJetsPF.getJetMCFlavour    = cms.bool(False)
patJetsPF.addGenPartonMatch   = cms.bool(False)
patJetsPF.embedGenPartonMatch = cms.bool(False)
patJetsPF.addGenJetMatch      = cms.bool(False)
patJetsPF.embedGenJetMatch    = cms.bool(False)
patJetsPF.addPartonJetMatch   = cms.bool(False)
patJetsPF.getJetMCFlavour     = cms.bool(False)

makeCustomizedPatJets = cms.Sequence(
#    patJetCorrections *
    patJetCharge *
    patJets *
    patJetsPF
)
makeCustomizedPatJetsForMC = cms.Sequence(
#    patJetCorrections *
    patJetCharge *
    patJetPartonMatch *
    patJetGenJetMatch *
    patJetFlavourId *
    patJets *
    patJetsPF
)

# --------------------Modifications for photons--------------------

patPhotons.addGenMatch = cms.bool(doGen)
patPhotons.embedGenMatch = cms.bool(doGen)

# --------------------Modifications for MET--------------------

patMETs.addGenMET    = cms.bool(doGen)
patMETs.metSource  = cms.InputTag("met")
patMETsPF = copy.deepcopy(patMETs)
patMETsPF.metSource  = cms.InputTag("pfMet")
patMETsPF.addMuonCorrections = cms.bool(False)
patMETsPF.addGenMET    = cms.bool(doGen)

# --------------------Modified PAT sequences--------------------


selectedPatCustomizedCandidates = cms.Sequence(
    selectedPatElectrons +
    heepPatElectrons	 +
    selectedPatMuons     +
    selectedLayer1ShrinkingConeHighEffPFTaus +
    selectedPatPhotons
#    selectedPatJets
)

if(data):
  patCustomizedCandidates = cms.Sequence(
    patElectrons +
    patMuons     +
    makeCustomizedPatTaus      +
    patPhotons   +
    makeCustomizedPatJets      +
#    patMETs + 
    patMETsPF 
)
else:
  patCustomizedCandidates = cms.Sequence(
    makePatElectrons +
    makePatMuons     +
    makeCustomizedPatTaus      +
    makePatPhotons   +
    makeCustomizedPatJetsForMC      +
#    patMETs + 
    patMETsPF
)

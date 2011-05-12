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

# turn off matching for data
patTaus.addDecayMode = cms.bool(doGen)
patTaus.addGenMatch = cms.bool(doGen)
patTaus.embedGenMatch = cms.bool(doGen)
patTaus.addGenJetMatch = cms.bool(doGen)
patTaus.embedGenJetMatch = cms.bool(doGen)

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

hpsTauMatch = copy.deepcopy(tauMatch)
hpsTauMatch.src         = cms.InputTag("hpsPFTauProducer")
hpsTauGenJetMatch = copy.deepcopy(tauGenJetMatch)
hpsTauGenJetMatch.src         = cms.InputTag("hpsPFTauProducer")
hpsPatTaus = copy.deepcopy(patTaus)
hpsPatTaus.tauSource = cms.InputTag("hpsPFTauProducer")
hpsPatTaus.isoDeposits = cms.PSet()
hpsPatTaus.userIsolation = cms.PSet()

hpsPatTaus.tauIDSources = cms.PSet(
        byVLooseIsolation = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolation"),
        byLooseIsolation = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation"),
        byMediumIsolation = cms.InputTag("hpsPFTauDiscriminationByMediumIsolation"),
        byTightIsolation = cms.InputTag("hpsPFTauDiscriminationByTightIsolation"),
        byDecayModeFinding = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
        againstLooseElectron = cms.InputTag("hpsPFTauDiscriminationByLooseElectronRejection"),
        againstMediumElectron = cms.InputTag("hpsPFTauDiscriminationByMediumElectronRejection"),
        againstTightElectron = cms.InputTag("hpsPFTauDiscriminationByTightElectronRejection"),
        againstLooseMuon = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection"),
        againstTightMuon = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection"),
)

hpsPatTaus.addDecayMode = cms.bool(False)
hpsPatTaus.decayModeSrc = cms.InputTag("")
hpsPatTaus.addGenMatch = cms.bool(doGen)
hpsPatTaus.embedGenMatch = cms.bool(doGen)
hpsPatTaus.genParticleMatch = cms.InputTag("hpsTauMatch")
hpsPatTaus.addGenJetMatch = cms.bool(doGen)
hpsPatTaus.embedGenJetMatch = cms.bool(doGen)
hpsPatTaus.genJetMatch = cms.InputTag("hpsTauGenJetMatch")

selectedLayer1HPSPFTaus = copy.deepcopy(selectedPatTaus)
selectedLayer1HPSPFTaus.src = cms.InputTag("hpsPatTaus")

makeCustomizedPatTauIso = cms.Sequence(
    patPFCandidateIsoDepositSelection
)

makeCustomizedPatTauMatch = cms.Sequence(
    shrinkingConeTauMatch *
    hpsTauMatch *  
    tauGenJets *
    tauGenJetsSelectorAllHadrons *
    shrinkingConeTauGenJetMatch *
    hpsTauGenJetMatch
) 

makePatTaus = cms.Sequence(
    shrinkingConePatTaus *
    hpsPatTaus
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

# modified jet corrections for calojets
patJetCorrFactorsCalo = copy.deepcopy(patJetCorrFactors)
patJetCorrFactorsCalo.levels = ['L1Offset', 'L2Relative', 'L3Absolute']
# modified jet corrections for pfjets
patJetCorrFactorsPF = copy.deepcopy(patJetCorrFactors)
patJetCorrFactorsPF.src = cms.InputTag("ak5PFJets")
patJetCorrFactorsPF.payload = cms.string('AK5PF')
patJetCorrFactorsPF.levels = ['L1Offset', 'L2Relative', 'L3Absolute']

# jet charge association for pfjets
patJetChargePF = copy.deepcopy(patJetCharge)
patJetChargePF.src = cms.InputTag("ak5PFJetTracksAssociatorAtVertex")

patJets.addJetCorrFactors    = cms.bool(True)
patJets.jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsCalo") )
patJets.addGenPartonMatch   = cms.bool(doGen)                           ## switch on/off matching to quarks from hard scatterin
patJets.embedGenPartonMatch = cms.bool(doGen)                           ## switch on/off embedding of the GenParticle parton for this jet
patJets.addGenJetMatch      = cms.bool(doGen)                           ## switch on/off matching to GenJet's
patJets.embedGenJetMatch    = cms.bool(doGen)                           ## switch on/off embedding of matched genJet's
patJets.addPartonJetMatch   = cms.bool(doGen)                          ## switch on/off matching to PartonJet's (not implemented yet)
patJets.getJetMCFlavour     = cms.bool(doGen)

patJetsPF = copy.deepcopy(patJets)
patJetsPF.jetSource = cms.InputTag("ak5PFJets")
patJetsPF.addJetCorrFactors    = cms.bool(True)
patJetsPF.jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsPF") )
patJetsPF.addBTagInfo          = cms.bool(True)
patJetsPF.addDiscriminators    = cms.bool(True)
patJetsPF.discriminatorSources = cms.VInputTag(
        cms.InputTag("combinedSecondaryVertexBJetTagsPF"),
        cms.InputTag("combinedSecondaryVertexMVABJetTagsPF"),
        cms.InputTag("jetBProbabilityBJetTagsPF"),
        cms.InputTag("jetProbabilityBJetTagsPF"),
        cms.InputTag("simpleSecondaryVertexHighEffBJetTagsPF"),
        cms.InputTag("simpleSecondaryVertexHighPurBJetTagsPF"),
        cms.InputTag("trackCountingHighEffBJetTagsPF"),
        cms.InputTag("trackCountingHighPurBJetTagsPF"),
)
patJetsPF.addJetID = cms.bool(False)
patJetsPF.addAssociatedTracks    = cms.bool(True)
patJetsPF.trackAssociationSource = cms.InputTag("ak5PFJetTracksAssociatorAtVertex")
patJetsPF.addJetCharge    = cms.bool(True)
patJetsPF.jetChargeSource = cms.InputTag("patJetChargePF")
patJetsPF.getJetMCFlavour    = cms.bool(False)
patJetsPF.addGenPartonMatch   = cms.bool(False)
patJetsPF.embedGenPartonMatch = cms.bool(False)
patJetsPF.addGenJetMatch      = cms.bool(False)
patJetsPF.embedGenJetMatch    = cms.bool(False)
patJetsPF.addPartonJetMatch   = cms.bool(False)
patJetsPF.getJetMCFlavour     = cms.bool(False)

from RecoBTag.Configuration.RecoBTag_cff import *
impactParameterTagInfosPF = copy.deepcopy(impactParameterTagInfos)
impactParameterTagInfosPF.jetTracks = cms.InputTag("ak5PFJetTracksAssociatorAtVertex")
impactParameterMVABJetTagsPF = copy.deepcopy(impactParameterMVABJetTags)
impactParameterMVABJetTagsPF.tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosPF"))
jetBProbabilityBJetTagsPF = copy.deepcopy(jetBProbabilityBJetTags)
jetBProbabilityBJetTagsPF.tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosPF"))
jetProbabilityBJetTagsPF = copy.deepcopy(jetProbabilityBJetTags)
jetProbabilityBJetTagsPF.tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosPF"))
trackCountingHighEffBJetTagsPF = copy.deepcopy(trackCountingHighEffBJetTags)
trackCountingHighEffBJetTagsPF.tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosPF"))
trackCountingHighPurBJetTagsPF = copy.deepcopy(trackCountingHighPurBJetTags)
trackCountingHighPurBJetTagsPF.tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosPF"))
secondaryVertexTagInfosPF = copy.deepcopy(secondaryVertexTagInfos)
secondaryVertexTagInfosPF.trackIPTagInfos = cms.InputTag("impactParameterTagInfosPF")
simpleSecondaryVertexHighEffBJetTagsPF = copy.deepcopy(simpleSecondaryVertexHighEffBJetTags)
simpleSecondaryVertexHighEffBJetTagsPF.tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosPF"))
simpleSecondaryVertexHighPurBJetTagsPF = copy.deepcopy(simpleSecondaryVertexHighPurBJetTags)
simpleSecondaryVertexHighPurBJetTagsPF.tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosPF"))
combinedSecondaryVertexBJetTagsPF = copy.deepcopy(combinedSecondaryVertexBJetTags)
combinedSecondaryVertexBJetTagsPF.tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosPF"),
                                                           cms.InputTag("secondaryVertexTagInfosPF"))
combinedSecondaryVertexMVABJetTagsPF = copy.deepcopy(combinedSecondaryVertexMVABJetTags)
combinedSecondaryVertexMVABJetTagsPF.tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosPF"),
                                                              cms.InputTag("secondaryVertexTagInfosPF"))
ghostTrackVertexTagInfosPF = copy.deepcopy(ghostTrackVertexTagInfos)
ghostTrackVertexTagInfosPF.trackIPTagInfos = cms.InputTag("impactParameterTagInfosPF")
ghostTrackBJetTagsPF = copy.deepcopy(ghostTrackBJetTags)
ghostTrackBJetTagsPF.tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosPF"),
                                              cms.InputTag("ghostTrackVertexTagInfosPF"))
simpleSecondaryVertexBJetTagsPF = simpleSecondaryVertexHighEffBJetTagsPF.clone()

makeCustomizedPatJets = cms.Sequence(
    patJetCorrFactorsCalo *
    patJetCorrFactorsPF *
    patJetCharge *
    patJetChargePF *
    impactParameterTagInfosPF *
    ( trackCountingHighEffBJetTagsPF +
      trackCountingHighPurBJetTagsPF +
      jetProbabilityBJetTagsPF +
      jetBProbabilityBJetTagsPF +
      secondaryVertexTagInfosPF *
      ( simpleSecondaryVertexHighEffBJetTagsPF +
        simpleSecondaryVertexHighPurBJetTagsPF +
        combinedSecondaryVertexBJetTagsPF + 
        combinedSecondaryVertexMVABJetTagsPF
      ) +
      ghostTrackVertexTagInfosPF *
        ghostTrackBJetTagsPF
    ) *
    patJets *
    patJetsPF
)
makeCustomizedPatJetsForMC = cms.Sequence(
    patJetCorrFactorsCalo *
    patJetCorrFactorsPF *
    patJetCharge *
    patJetChargePF *
    impactParameterTagInfosPF *
    ( trackCountingHighEffBJetTagsPF +
      trackCountingHighPurBJetTagsPF +
      jetProbabilityBJetTagsPF +
      jetBProbabilityBJetTagsPF +
      secondaryVertexTagInfosPF *
      ( simpleSecondaryVertexHighEffBJetTagsPF +
        simpleSecondaryVertexHighPurBJetTagsPF +
        combinedSecondaryVertexBJetTagsPF + 
        combinedSecondaryVertexMVABJetTagsPF
      ) +
      ghostTrackVertexTagInfosPF *
        ghostTrackBJetTagsPF
    ) *
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

# Compute TypeI/TypeII MET
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from RecoJets.Configuration.RecoPFJets_cff import *
kt6PFJets.doRhoFastjet = True
kt6PFJets.Rho_EtaMax = cms.double(2.5)

from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
selectedPFJetsReco = ak5PFJets.clone()
selectedPFJetsReco.src = cms.InputTag('particleFlow')
selectedPFJetsReco.doAreaFastjet = True
selectedPFJetsReco.Rho_EtaMax = cms.double(5.0)

from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
ak5PFJetsC = ak5PFJets.clone()
ak5PFJetsC.src = cms.InputTag('particleFlow')
ak5PFJetsC.doAreaFastjet = True
ak5PFJetsC.Rho_EtaMax = cms.double(5.0)

from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5PFJet
metJESCorPFAK5 = metJESCorAK5PFJet.clone()
metJESCorPFAK5.inputUncorJetsLabel = "ak5PFJetsC"
metJESCorPFAK5.metType = "PFMET"
metJESCorPFAK5.inputUncorMetLabel = "pfMet"
metJESCorPFAK5.useTypeII = False
metJESCorPFAK5.jetPTthreshold = cms.double(10.0)
metJESCorPFAK5.corrector = cms.string('ak5PFL2L3')

metJES10CorPFAK5 = metJESCorAK5PFJet.clone()
metJES10CorPFAK5.inputUncorJetsLabel = "ak5PFJetsC"
metJES10CorPFAK5.metType = "PFMET"
metJES10CorPFAK5.inputUncorMetLabel = "pfMet"
metJES10CorPFAK5.useTypeII = False
metJES10CorPFAK5.jetPTthreshold = cms.double(10.0)
metJES10CorPFAK5.corrector = cms.string('ak5PFL1FastL2L3')

patMETs.addGenMET    = cms.bool(doGen)
patMETs.metSource  = cms.InputTag("met")
patMETsPFL1L2L3Cor = copy.deepcopy(patMETs)
patMETsPFL1L2L3Cor.metSource  = cms.InputTag("metJES10CorPFAK5")
patMETsPFL1L2L3Cor.addMuonCorrections = cms.bool(False)
patMETsPFL1L2L3Cor.addGenMET    = cms.bool(doGen)

patMETsPFL2L3Cor = copy.deepcopy(patMETs)
patMETsPFL2L3Cor.metSource  = cms.InputTag("metJESCorPFAK5")
patMETsPFL2L3Cor.addMuonCorrections = cms.bool(False)
patMETsPFL2L3Cor.addGenMET    = cms.bool(doGen)

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
    selectedLayer1HPSPFTaus +
    selectedPatPhotons
)

if(data):
  patCustomizedCandidates = cms.Sequence(
    patElectrons +
    patMuons     +
    makeCustomizedPatTaus      +
    patPhotons   +
    makeCustomizedPatJets      +
#    patMETs + 
    kt6PFJets + # needed for L1FastJet corrections
    ak5PFJetsC + # jet clustering with a given value of Rho_EtaMax for L1FastJet
    ak5PFJetsL2L3 + # compute JEC
    metJESCorPFAK5 + # compute typeI pfmet ak5PFL2L3
    metJES10CorPFAK5 + # compute typeI pfmet ak5PFL1FastL2L3
    patMETsPF +
    patMETsPFL1L2L3Cor +
    patMETsPFL2L3Cor
)
else:
  patCustomizedCandidates = cms.Sequence(
    makePatElectrons +
    makePatMuons     +
    makeCustomizedPatTaus      +
    makePatPhotons   +
    makeCustomizedPatJetsForMC      +
#    patMETs + 
    kt6PFJets + # needed for L1FastJet corrections
    ak5PFJetsC + # jet clustering with a given value of Rho_EtaMax for L1FastJet
    ak5PFJetsL2L3 + # compute JEC
    metJESCorPFAK5 + # compute typeI pfmet ak5PFL2L3
    metJES10CorPFAK5 + # compute typeI pfmet ak5PFL1FastL2L3
    patMETsPF +
    patMETsPFL1L2L3Cor +
    patMETsPFL2L3Cor
)

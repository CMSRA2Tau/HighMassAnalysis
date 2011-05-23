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

from RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi import *
eidLooseNoCuts = eidCutBasedExt.clone()
eidLooseNoCuts.electronIDType = 'classbased'
eidLooseNoCuts.electronQuality = 'loose'
eidLooseNoCuts.classbasedlooseEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutdetainl = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cutdphiin = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutdphiinl = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cuteseedopcor = cms.vdouble(-9e+9, -9e+9, -9e+9, -9e+9, -9e+9, -9e+9, -9e+9, -9e+9, -9e+9), ## Alfredo's modifications
        cutfmishits = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutdcotdist = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cuthoe = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cuthoel = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cutip_gsf = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutip_gsfl = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cutiso_sum = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutiso_sumoet = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutiso_sumoetl = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cutsee = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutseel = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5) ## Alfredo's modifications
)
eidLooseOnlyID = eidCutBasedExt.clone()
eidLooseOnlyID.electronIDType = 'classbased'
eidLooseOnlyID.electronQuality = 'loose'
eidLooseOnlyID.classbasedlooseEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(1.37e-02, 6.78e-03, 2.41e-02, 1.87e-02, 1.61e-02, 2.24e-02, 2.52e-02, 3.08e-02, 2.73e-02),
        cutdetainl = cms.vdouble(1.24e-02, 5.03e-03, 2.57e-02, 2.28e-02, 1.18e-02, 1.78e-02, 1.88e-02, 1.40e-01, 2.40e-02),
        cutdphiin = cms.vdouble(8.97e-02, 2.62e-01, 3.53e-01, 1.16e-01, 3.57e-01, 3.19e-01, 3.42e-01, 4.04e-01, 3.36e-01),
        cutdphiinl = cms.vdouble(7.47e-02, 2.50e-01, 3.56e-01, 9.56e-02, 3.47e-01, 3.26e-01, 3.33e-01, 6.47e-01, 2.89e-01),
        cuteseedopcor = cms.vdouble(6.30e-01, 8.20e-01, 4.01e-01, 7.18e-01, 4.00e-01, 4.58e-01, 1.50e-01, 6.64e-01, 3.73e-01),
        cutfmishits = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutdcotdist = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cuthoe = cms.vdouble(2.47e-01, 1.37e-01, 1.47e-01, 3.71e-01, 5.88e-02, 1.47e-01, 5.20e-01, 4.52e-01, 4.04e-01),
        cuthoel = cms.vdouble(2.36e-01, 1.26e-01, 1.47e-01, 3.75e-01, 3.92e-02, 1.45e-01, 3.65e-01, 3.83e-01, 3.84e-01),
        cutip_gsf = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutip_gsfl = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cutiso_sum = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutiso_sumoet = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutiso_sumoetl = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cutsee = cms.vdouble(1.76e-02, 1.25e-02, 1.81e-02, 4.15e-02, 3.64e-02, 4.18e-02, 1.46e-02, 6.78e-02, 1.33e-01),
        cutseel = cms.vdouble(1.64e-02, 1.18e-02, 1.50e-02, 5.23e-02, 3.26e-02, 4.56e-02, 1.85e-02, 5.89e-02, 5.44e-02)
)
eidLooseOnlyConversionRejection = eidCutBasedExt.clone()
eidLooseOnlyConversionRejection.electronIDType = 'classbased'
eidLooseOnlyConversionRejection.electronQuality = 'loose'
eidLooseOnlyConversionRejection.classbasedlooseEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutdetainl = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cutdphiin = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutdphiinl = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cuteseedopcor = cms.vdouble(-9e+9, -9e+9, -9e+9, -9e+9, -9e+9, -9e+9, -9e+9, -9e+9, -9e+9), ## Alfredo's modifications
        cutfmishits = cms.vdouble(4.50e+00, 1.50e+00, 1.50e+00, 2.50e+00, 2.50e+00, 1.50e+00, 4.50e+00, 3.50e+00, 3.50e+00),
        cutdcotdist = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cuthoe = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cuthoel = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cutip_gsf = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutip_gsfl = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cutiso_sum = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutiso_sumoet = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutiso_sumoetl = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cutsee = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutseel = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5) ## Alfredo's modifications
)
eidLooseOnlyIP = eidCutBasedExt.clone()
eidLooseOnlyIP.electronIDType = 'classbased'
eidLooseOnlyIP.electronQuality = 'loose'
eidLooseOnlyIP.classbasedlooseEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutdetainl = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cutdphiin = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutdphiinl = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cuteseedopcor = cms.vdouble(-9e+9, -9e+9, -9e+9, -9e+9, -9e+9, -9e+9, -9e+9, -9e+9, -9e+9), ## Alfredo's modifications
        cutfmishits = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutdcotdist = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cuthoe = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cuthoel = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cutip_gsf = cms.vdouble(5.51e-02, 7.65e-02, 1.43e-01, 8.74e-02, 5.94e-01, 3.70e-01, 9.13e-02, 1.15e+00, 2.31e-01),
        cutip_gsfl = cms.vdouble(1.86e-02, 7.59e-02, 1.38e-01, 4.73e-02, 6.20e-01, 3.04e-01, 1.09e-01, 7.75e-01, 4.79e-02),
        cutiso_sum = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutiso_sumoet = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutiso_sumoetl = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cutsee = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutseel = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5) ## Alfredo's modifications
)
eidLooseOnlyIso = eidCutBasedExt.clone()
eidLooseOnlyIso.electronIDType = 'classbased'
eidLooseOnlyIso.electronQuality = 'loose'
eidLooseOnlyIso.classbasedlooseEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutdetainl = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cutdphiin = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutdphiinl = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cuteseedopcor = cms.vdouble(-9e+9, -9e+9, -9e+9, -9e+9, -9e+9, -9e+9, -9e+9, -9e+9, -9e+9), ## Alfredo's modifications
        cutfmishits = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutdcotdist = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cuthoe = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cuthoel = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cutip_gsf = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutip_gsfl = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5), ## Alfredo's modifications
        cutiso_sum = cms.vdouble(3.30e+01, 1.70e+01, 1.79e+01, 1.88e+01, 8.55e+00, 1.25e+01, 1.76e+01, 1.85e+01, 2.98e+00),
        cutiso_sumoet = cms.vdouble(3.45e+01, 1.27e+01, 1.21e+01, 1.99e+01, 6.35e+00, 8.85e+00, 1.40e+01, 1.05e+01, 9.74e+00),
        cutiso_sumoetl = cms.vdouble(1.13e+01, 9.05e+00, 9.07e+00, 9.94e+00, 5.25e+00, 6.15e+00, 1.07e+01, 1.08e+01, 4.40e+00),
        cutsee = cms.vdouble(9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9, 9e+9), ## Alfredo's modifications
        cutseel = cms.vdouble(9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5, 9e+5) ## Alfredo's modifications
)

eIdSequence = cms.Sequence(eidRobustLoose +
                           eidRobustTight +
                           eidRobustHighEnergy +
                           eidLoose +
                           eidLooseNoCuts +
                           eidLooseOnlyID + 
                           eidLooseOnlyConversionRejection +
                           eidLooseOnlyIP +
                           eidLooseOnlyIso +
                           eidTight)

patElectrons.addElectronID = cms.bool(True)
patElectrons.electronIDSources = cms.PSet(
    eidRobustLoose      = cms.InputTag("eidRobustLoose"),
    eidRobustTight      = cms.InputTag("eidRobustTight"),
    eidLoose            = cms.InputTag("eidLoose"),
    eidLooseNoCuts      = cms.InputTag("eidLooseNoCuts"),
    eidLooseOnlyID      = cms.InputTag("eidLooseOnlyID"),
    eidLooseOnlyConversionRejection            = cms.InputTag("eidLooseOnlyConversionRejection"),
    eidLooseOnlyIP      = cms.InputTag("eidLooseOnlyIP"),
    eidLooseOnlyIso     = cms.InputTag("eidLooseOnlyIso"),
    eidTight            = cms.InputTag("eidTight"),
    eidRobustHighEnergy = cms.InputTag("eidRobustHighEnergy"),
)

from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
heepPatElectrons = cms.EDProducer("HEEPAttStatusToPAT",
                                          eleLabel = cms.InputTag("selectedPatElectrons"),
                                          barrelCuts = cms.PSet(heepBarrelCuts),
                                          endcapCuts = cms.PSet(heepEndcapCuts)
                                          )

# --------------------Modifications for jets--------------------

from JetMETCorrections.Configuration.DefaultJEC_cff import *
from RecoJets.Configuration.RecoPFJets_cff import *

ak5CaloL1Offset.useCondDB = False
ak5PFL1Offset.useCondDB = False
ak5JPTL1Offset.useCondDB = False

ak5PFL1Fastjet.useCondDB = False
ak5CaloL1Fastjet.useCondDB = False
ak5JPTL1Fastjet.useCondDB = False

# modified jet corrections for calojets
patJetCorrFactorsCalo = copy.deepcopy(patJetCorrFactors)
patJetCorrFactorsCalo.levels = ['L1FastJet', 'L2Relative', 'L3Absolute']
# modified jet corrections for pfjets
patJetCorrFactorsPF = copy.deepcopy(patJetCorrFactors)
patJetCorrFactorsPF.src = cms.InputTag("ak5PFJets")
patJetCorrFactorsPF.payload = cms.string('AK5PF')
patJetCorrFactorsPF.levels = ['L1FastJet', 'L2Relative', 'L3Absolute']

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
    kt6PFJets + # needed for L1FastJet corrections
    ak5PFJetsC + # jet clustering with a given value of Rho_EtaMax for L1FastJet
    ak5PFJetsL2L3 + # compute JEC
    makeCustomizedPatJets      +
#    patMETs + 
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
    kt6PFJets + # needed for L1FastJet corrections
    ak5PFJetsC + # jet clustering with a given value of Rho_EtaMax for L1FastJet
    ak5PFJetsL2L3 + # compute JEC
    makeCustomizedPatJetsForMC      +
#    patMETs + 
    metJESCorPFAK5 + # compute typeI pfmet ak5PFL2L3
    metJES10CorPFAK5 + # compute typeI pfmet ak5PFL1FastL2L3
    patMETsPF +
    patMETsPFL1L2L3Cor +
    patMETsPFL2L3Cor
)

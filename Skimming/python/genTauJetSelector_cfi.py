import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# selection of tau --> e nu_tau nu_e decays
#--------------------------------------------------------------------------------

# require generated tau to decay into electrons
genTauDecaysToElectronCands = cms.EDFilter("TauGenJetDecayModeSelector",
     src = cms.InputTag("tauGenJets"),
     select = cms.vstring('electron'),
     filter = cms.bool(True)
)

# require generator level electron produced in tau-decay to be within muon acceptance
selectedGenTauDecaysToElectronEta = cms.EDFilter("TauGenJetSelector",
     src = cms.InputTag("genTauDecaysToElectronCands"),
     cut = cms.string('abs(eta) < 2.5'),
     filter = cms.bool(True)
)

# require generator level electron produced in tau-decay to have transverse momentum above threshold
selectedGenTauDecaysToElectronPt = cms.EDFilter("TauGenJetSelector",
     src = cms.InputTag("selectedGenTauDecaysToElectronEta"),
     cut = cms.string('pt > 8.'),
     filter = cms.bool(True)
)

#require 2 gen level electrons with pt and eta cuts
selectedGenTauElectronPairs = cms.EDFilter("CandViewCountFilter",
    src =  cms.InputTag("selectedGenTauDecaysToElectronPt"),
    minNumber = cms.uint32(2)
)

#--------------------------------------------------------------------------------
# selection of tau --> mu nu_tau nu_mu decays
#--------------------------------------------------------------------------------

# require generated tau to decay into muon
genTauDecaysToMuonCands = cms.EDFilter("TauGenJetDecayModeSelector",
     src = cms.InputTag("tauGenJets"),
     select = cms.vstring('muon'),
     filter = cms.bool(True)
)

# require generator level muon produced in tau-decay to be within muon acceptance
selectedGenTauDecaysToMuonEta = cms.EDFilter("TauGenJetSelector",
     src = cms.InputTag("genTauDecaysToMuonCands"),
     cut = cms.string('abs(eta) < 2.5'),
     filter = cms.bool(True)
)

# require generator level muon produced in tau-decay to have transverse momentum above threshold
selectedGenTauDecaysToMuonPt = cms.EDFilter("TauGenJetSelector",
     src = cms.InputTag("selectedGenTauDecaysToMuonEta"),
     cut = cms.string('pt > 8.'),
     filter = cms.bool(True)
)

#require 2 gen level muons with pt and eta cuts
selectedGenTauMuonPairs = cms.EDFilter("CandViewCountFilter",
    src =  cms.InputTag("selectedGenTauDecaysToMuonPt"),
    minNumber = cms.uint32(2)
)

#--------------------------------------------------------------------------------
# selection of tau --> hadrons nu_tau decays
#--------------------------------------------------------------------------------

# require generated tau to decay hadronically
genTauDecaysToHadronsCands = cms.EDFilter("TauGenJetDecayModeSelector",
     src = cms.InputTag("tauGenJets"),
     select = cms.vstring('oneProng0Pi0', 'oneProng1Pi0', 'oneProng2Pi0', 'oneProngOther',
                          'threeProng0Pi0', 'threeProng1Pi0', 'threeProngOther', 'rare'),
     filter = cms.bool(True)
)

TwoGenTauDecaysToHadronsCandsFilter = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('genTauDecaysToHadronsCands'),
  minNumber = cms.uint32(2)
)

# require generator level hadrons produced in tau-decay to be within muon acceptance
selectedGenTauDecaysToHadronsEta = cms.EDFilter("TauGenJetSelector",
     src = cms.InputTag("genTauDecaysToHadronsCands"),
     cut = cms.string('abs(eta) < 2.5'),
     filter = cms.bool(True)
)

TwoGenTauDecaysToHadronsEtaFilter = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('selectedGenTauDecaysToHadronsEta'),
  minNumber = cms.uint32(2)
)

# require generator level hadrons produced in tau-decay to have transverse momentum above threshold
selectedGenTauDecaysToHadronsPt = cms.EDFilter("TauGenJetSelector",
     src = cms.InputTag("selectedGenTauDecaysToHadronsEta"),
     cut = cms.string('pt > 8.'),
     filter = cms.bool(True)
)

TwoGenTauDecaysToHadronsPtFilter = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('selectedGenTauDecaysToHadronsPt'),
  minNumber = cms.uint32(2)
)


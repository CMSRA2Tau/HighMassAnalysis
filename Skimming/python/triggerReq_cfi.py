import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt

emuHLTFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::HLT'),
    triggerConditions = (
      'HLT_Mu17_Ele8_CaloIdL*',
      'HLT_Mu17_Ele8*',
      'HLT_Mu8_Ele17*',
      'HLT_Mu8_Ele17_CaloIdL*',
      'HLT_IsoMu24_eta2p1*',
    ),
    l1tResults = '',
    throw = False
)

eeHLTFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::HLT'),
    triggerConditions = (
      'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL*',
    ),
    l1tResults = '',
    throw = False
)

mumuHLTFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::HLT'),
    triggerConditions = (
      'HLT_Mu17_Mu8*',
    ),
    l1tResults = '',
    throw = False
)

etauHLTFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::HLT'),
    triggerConditions = (
      'HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20*',
      'HLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20*',
      'HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20*',
      'HLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20*',
    ),
    l1tResults = '',
    throw = False
)

mutauHLTFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::HLT'),
    triggerConditions = (
      'HLT_IsoMu24*',
      'HLT_Mu24*',
    ),
    l1tResults = '',
    throw = False
)

tautauHLTFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::HLT'),
    triggerConditions = (
      'HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1_Prong1*',
      'HLT_DoubleMediumIsoPFTau*_Trk*_eta2p1_Prong1_*',
      'HLT_DoubleMediumIsoPFTau*_Trk*_eta2p1_*',
      'HLT_Double*PFTau*_Trk*_eta2p1_*',
    ),
    l1tResults = '',
    throw = False
)

SusyHLTFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::HLT'),
    triggerConditions = (
      'HLT_IsoMu24_eta2p1_v*',
      'HLT_IsoMu24_eta2p1*',
    ),
    l1tResults = '',
    throw = False
)


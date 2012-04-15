import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt

emuHLTFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::HLT'),
    triggerConditions = (
      'HLT_Mu17_Ele8_CaloIdL',
      'HLT_Mu17_Ele8_CaloIdL*',
      'HLT_Mu8_Ele17_CaloIdL',
      'HLT_Mu8_Ele17_CaloIdL*',
      'HLT_Mu10_Ele10_CaloIdL',
      'HLT_Mu10_Ele10_CaloIdL*',
      'HLT_Mu24',
      'HLT_Mu24*',
      'HLT_IsoMu17',
      'HLT_IsoMu17*',
    ),
    l1tResults = '',
    throw = False
)

etauHLTFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::HLT'),
    triggerConditions = (
      'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15*',
      'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20*',
      'HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20*'
    ),
    l1tResults = '',
    throw = False
)

mutauHLTFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::HLT'),
    triggerConditions = (
      'HLT_Mu24',
      'HLT_Mu24*',
      'HLT_IsoMu17',
      'HLT_IsoMu17*',
      'HLT_IsoMu12_LooseIsoPFTau10',
      'HLT_IsoMu12_LooseIsoPFTau10*',
      'HLT_Mu15_LooseIsoPFTau20',
      'HLT_Mu15_LooseIsoPFTau20*',
      'HLT_IsoMu15_LooseIsoPFTau15',
      'HLT_IsoMu15_LooseIsoPFTau15*',
      'HLT_IsoMu15_LooseIsoPFTau20',
      'HLT_IsoMu15_LooseIsoPFTau20*',
    ),
    l1tResults = '',
    throw = False
)

tautauHLTFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::HLT'),
    triggerConditions = (
      'HLT_DoubleIsoPFTau20_Trk5',
      'HLT_DoubleIsoPFTau20_Trk5*',
      'HLT_DoubleIsoPFTau*',
      'HLT_IsoPFTau*_IsoPFTau*',
    ),
    l1tResults = '',
    throw = False
)

SusyHLTFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::HLT'),
    triggerConditions = (
      'HLT_PFMHT150',
      'HLT_PFMHT150*',
      'HLT_PFMET150',
      'HLT_PFMET150*',
      'HLT_HT300_MHT75',
      'HLT_HT300_MHT75*',
    ),
    l1tResults = '',
    throw = False
)


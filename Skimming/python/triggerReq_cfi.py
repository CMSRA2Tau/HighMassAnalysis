import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt

emuHLTFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::HLT'),
    triggerConditions = (
      'HLT_Mu9',
      'HLT_Mu11',
      'HLT_IsoMu9*',
      'HLT_IsoMu11*',
      'HLT_IsoMu13*',
      'HLT_IsoMu15*',
      'HLT_Mu15*',
      'HLT_Mu5_Elec9*',
      'HLT_Mu11_Ele8*',
      'HLT_Mu8_Ele8*'
    ),
    l1tResults = '',
    throw = False
)

etauHLTFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::HLT'),
    triggerConditions = (
      'HLT_Ele10_LW_EleId_L1R*',
      'HLT_Ele15_LW_L1R*',
      'HLT_Ele15_SW_L1R*',
      'HLT_Ele12_SW_TightEleIdIsol_L1R*',
      'HLT_Ele17_SW_EleId_L1R*',
      'HLT_Ele17_SW_LooseEleId_L1R*',
      'HLT_Ele12_SW_TighterEleIdIsol_L1R*',
      'HLT_Ele17_SW_TightEleId_L1R*',
      'HLT_Ele17_SW_TighterEleId_L1R*',
      'HLT_Ele17_SW_TighterEleIdIsol_L1R*',
      'HLT_Ele22_SW_TighterEleId_L1R*',
      'HLT_IsoEle12_PFTau15*',
      'HLT_IsoEle12_PFTau15*',
      'HLT_IsoEle12_PFTau15*'
    ),
    l1tResults = '',
    throw = False
)

mutauHLTFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::HLT'),
    triggerConditions = (
      'HLT_Mu9',
      'HLT_Mu11',
      'HLT_IsoMu9',
      'HLT_IsoMu11*',
      'HLT_IsoMu13*',
      'HLT_IsoMu15*',
      'HLT_Mu15*',
      'HLT_IsoMu9_PFTau15*'
      'HLT_Mu11_PFTau15*'
    ),
    l1tResults = '',
    throw = False
)

tautauHLTFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::HLT'),
    triggerConditions = (
      'HLT_DoubleLooseIsoTau15',
      'HLT_DoubleIsoTau15_Trk5'
    ),
    l1tResults = '',
    throw = False
)

import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt

muTauHLTFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::HLT'),
    triggerConditions = (
      'HLT_Mu9',
      'HLT_Mu11'
    ),
    l1tResults = '',
    throw = False
)

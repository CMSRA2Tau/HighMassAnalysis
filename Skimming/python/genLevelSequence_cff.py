import FWCore.ParameterSet.Config as cms

from PhysicsTools.JetMCAlgos.TauGenJets_cfi import * 
from HiMassTauAnalyzer.Skimming.genTauJetSelector_cfi import *  

genLevelElecTauSequence = cms.Sequence(
    tauGenJets    
 *  genTauDecaysToElectronCands
 *  genTauDecaysToHadronsCands
 *  selectedGenTauDecaysToElectronEta
 *  selectedGenTauDecaysToElectronPt   
 *  selectedGenTauDecaysToHadronsEta
 *  selectedGenTauDecaysToHadronsPt
)

genLevelMuTauSequence = cms.Sequence(
    tauGenJets    
 *  genTauDecaysToMuonCands
 *  genTauDecaysToHadronsCands
 *  selectedGenTauDecaysToMuonEta
 *  selectedGenTauDecaysToMuonPt 
 *  selectedGenTauDecaysToHadronsEta
 *  selectedGenTauDecaysToHadronsPt
)

genLevelElecMuSequence = cms.Sequence(
    tauGenJets
 *  genTauDecaysToMuonCands
 *  genTauDecaysToElectronCands
 *  selectedGenTauDecaysToMuonEta
 *  selectedGenTauDecaysToMuonPt 
 *  selectedGenTauDecaysToElectronEta
 *  selectedGenTauDecaysToElectronPt   
)

genLevelMuMuSequence = cms.Sequence(
    tauGenJets
 *  genTauDecaysToMuonCands
 *  selectedGenTauDecaysToMuonEta
 *  selectedGenTauDecaysToMuonPt
 *  selectedGenTauMuonPairs
)

genLevelElecElecSequence = cms.Sequence(
    tauGenJets
 *  genTauDecaysToElectronCands
 *  selectedGenTauDecaysToElectronEta
 *  selectedGenTauDecaysToElectronPt
 *  selectedGenTauElectronPairs
)

genLevelTauTauSequence = cms.Sequence(
    tauGenJets
 *  genTauDecaysToHadronsCands
 *  TwoGenTauDecaysToHadronsCandsFilter
 *  selectedGenTauDecaysToHadronsEta
 *  TwoGenTauDecaysToHadronsEtaFilter
 *  selectedGenTauDecaysToHadronsPt
 *  TwoGenTauDecaysToHadronsPtFilter
)

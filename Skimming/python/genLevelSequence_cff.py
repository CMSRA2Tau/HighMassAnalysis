import FWCore.ParameterSet.Config as cms

from PhysicsTools.JetMCAlgos.TauGenJets_cfi import * 
from HighMassAnalysis.Skimming.genTauJetSelector_cfi import *  

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
 *  TwoGenTauDecaysToMuonCandsFilter
 *  selectedGenTauDecaysToMuonEta
 *  TwoGenTauDecaysToMuonCandsEtaFilter
 *  selectedGenTauDecaysToMuonPt
 *  TwoGenTauDecaysToMuonCandsPtFilter
)

genLevelElecElecSequence = cms.Sequence(
    tauGenJets
 *  genTauDecaysToElectronCands
 *  TwoGenTauDecaysToElectronCandsFilter
 *  selectedGenTauDecaysToElectronEta
 *  TwoGenTauDecaysToElectronCandsEtaFilter
 *  selectedGenTauDecaysToElectronPt
 *  TwoGenTauDecaysToElectronCandsPtFilter
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

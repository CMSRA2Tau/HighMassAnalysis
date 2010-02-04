import FWCore.ParameterSet.Config as cms

# import customized settings
from HighMassAnalysis.Configuration.customizePAT_cff import *

# produce collections of dR = 0.07 and dR = 0.15 fixed
# and dR = 5.0/Et shrinking signal cone taus using latest tags
from RecoTauTag.Configuration.RecoPFTauTag_cff import *
pfRecoTauTagInfoProducer.ChargedHadrCand_AssociationCone   = cms.double(1.0)

ic5PFJetTracksAssociatorAtVertex.j2tParametersVX = cms.PSet(
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(1.0)
)

#--- tracker isolation
fixedConeHighEffPFTauProducer.TrackerIsolConeSizeFormula   = cms.string('1.00') ## **
fixedConeHighEffPFTauProducer.TrackerIsolConeSize_max = cms.double(1.1)
#--- ecal isolation
fixedConeHighEffPFTauProducer.ECALIsolConeSizeFormula      = cms.string('1.00') ## **
fixedConeHighEffPFTauProducer.ECALIsolConeSize_max = cms.double(1.1)
#--- hcal isolation
fixedConeHighEffPFTauProducer.HCALIsolConeSizeFormula      = cms.string('1.00') ## **
fixedConeHighEffPFTauProducer.HCALIsolConeSize_max = cms.double(1.1)
#--- tracker signal
fixedConeHighEffPFTauProducer.TrackerSignalConeSize_min = cms.double(0.00)
fixedConeHighEffPFTauProducer.TrackerSignalConeSize_max = cms.double(0.15)
#--- ecal signal
fixedConeHighEffPFTauProducer.ECALSignalConeSizeFormula    = cms.string('0.15') ## **
fixedConeHighEffPFTauProducer.ECALSignalConeSize_min    = cms.double(0.00)
fixedConeHighEffPFTauProducer.ECALSignalConeSize_max    = cms.double(0.15)
#--- hcal signal
fixedConeHighEffPFTauProducer.HCALSignalConeSizeFormula    = cms.string('0.15') ## **
fixedConeHighEffPFTauProducer.HCALSignalConeSize_min    = cms.double(0.00)
fixedConeHighEffPFTauProducer.HCALSignalConeSize_max    = cms.double(0.15)


#--- tracker isolation
fixedConePFTauProducer.TrackerIsolConeSizeFormula   = cms.string('1.00') ## **
fixedConePFTauProducer.TrackerIsolConeSize_max = cms.double(1.1)
#--- ecal isolation
fixedConePFTauProducer.ECALIsolConeSizeFormula      = cms.string('1.00') ## **
fixedConePFTauProducer.ECALIsolConeSize_max = cms.double(1.1)
#--- hcal isolation
fixedConePFTauProducer.HCALIsolConeSizeFormula      = cms.string('1.00') ## **
fixedConePFTauProducer.HCALIsolConeSize_max = cms.double(1.1)
#--- tracker signal
fixedConePFTauProducer.TrackerSignalConeSize_min = cms.double(0.00)
fixedConePFTauProducer.TrackerSignalConeSize_max = cms.double(0.07)
#--- ecal signal
fixedConePFTauProducer.ECALSignalConeSizeFormula    = cms.string('0.15') ## **
fixedConePFTauProducer.ECALSignalConeSize_min    = cms.double(0.00)
fixedConePFTauProducer.ECALSignalConeSize_max    = cms.double(0.15)
#--- hcal signal
fixedConePFTauProducer.HCALSignalConeSizeFormula    = cms.string('0.15') ## **
fixedConePFTauProducer.HCALSignalConeSize_min    = cms.double(0.00)
fixedConePFTauProducer.HCALSignalConeSize_max    = cms.double(0.15)


#--- tracker isolation
shrinkingConePFTauProducer.TrackerIsolConeSizeFormula   = cms.string('1.00') ## **
shrinkingConePFTauProducer.TrackerIsolConeSize_max   = cms.double(1.1)
#--- ecal isolation
shrinkingConePFTauProducer.ECALIsolConeSizeFormula      = cms.string('1.00') ## **
shrinkingConePFTauProducer.ECALIsolConeSize_max      = cms.double(1.1)
#--- hcal isolation
shrinkingConePFTauProducer.HCALIsolConeSizeFormula      = cms.string('1.00') ## **
shrinkingConePFTauProducer.HCALIsolConeSize_max      = cms.double(1.1)
#--- tracker signal
shrinkingConePFTauProducer.TrackerSignalConeSize_min = cms.double(0.07)
shrinkingConePFTauProducer.TrackerSignalConeSize_max = cms.double(0.15)
#--- ecal signal
shrinkingConePFTauProducer.ECALSignalConeSizeFormula    = cms.string('0.15') ## **
shrinkingConePFTauProducer.ECALSignalConeSize_min    = cms.double(0.00)
shrinkingConePFTauProducer.ECALSignalConeSize_max    = cms.double(0.15)
#--- hcal signal
shrinkingConePFTauProducer.HCALSignalConeSizeFormula    = cms.string('0.15') ## **
shrinkingConePFTauProducer.HCALSignalConeSize_min    = cms.double(0.00)
shrinkingConePFTauProducer.HCALSignalConeSize_max    = cms.double(0.15)


# tau collection w/ shrinking cone of 3/ET
shrinkingTightConePFTauProducer = copy.deepcopy(pfRecoTauProducer)
shrinkingTightConePFTauProducer.LeadPFCand_minPt      = cms.double(5.0)
#--- tracker signal
shrinkingTightConePFTauProducer.TrackerSignalConeSizeFormula = cms.string('3/ET') ## **
shrinkingTightConePFTauProducer.TrackerSignalConeSize_min = cms.double(0.07)
shrinkingTightConePFTauProducer.TrackerSignalConeSize_max = cms.double(0.15)
#--- ecal signal
shrinkingTightConePFTauProducer.ECALSignalConeSizeFormula    = cms.string('0.15') ## **
shrinkingTightConePFTauProducer.ECALSignalConeSize_min    = cms.double(0.00)
shrinkingTightConePFTauProducer.ECALSignalConeSize_max    = cms.double(0.15)
#--- hcal signal
shrinkingTightConePFTauProducer.HCALSignalConeSizeFormula    = cms.string('0.15') ## **
shrinkingTightConePFTauProducer.HCALSignalConeSize_min    = cms.double(0.00)
shrinkingTightConePFTauProducer.HCALSignalConeSize_max    = cms.double(0.15)
#--- tracker isolation
shrinkingTightConePFTauProducer.TrackerIsolConeSizeFormula   = cms.string('1.00') ## **
shrinkingTightConePFTauProducer.TrackerIsolConeSize_min   = cms.double(0.0)
shrinkingTightConePFTauProducer.TrackerIsolConeSize_max   = cms.double(1.1)
#--- ecal isolation
shrinkingTightConePFTauProducer.ECALIsolConeSizeFormula      = cms.string('1.00') ## **
shrinkingTightConePFTauProducer.ECALIsolConeSize_min      = cms.double(0.0)
shrinkingTightConePFTauProducer.ECALIsolConeSize_max      = cms.double(1.1)
#--- hcal isolation
shrinkingTightConePFTauProducer.HCALIsolConeSizeFormula      = cms.string('1.00') ## **
shrinkingTightConePFTauProducer.HCALIsolConeSize_min      = cms.double(0.0)
shrinkingTightConePFTauProducer.HCALIsolConeSize_max      = cms.double(1.1)

from RecoTauTag.RecoTau.PFRecoTauDecayModeDeterminator_cfi                          import *
shrinkingTightConePFTauDecayModeProducer               = copy.deepcopy(pfTauDecayMode)
shrinkingTightConePFTauDecayModeProducer.PFTauProducer = 'shrinkingTightConePFTauProducer'
from RecoTauTag.RecoTau.PFRecoTauDecayModeIndexProducer_cfi                             import *
shrinkingTightConePFTauDecayModeIndexProducer                        = copy.deepcopy(pfTauDecayModeIndexProducer)
shrinkingTightConePFTauDecayModeIndexProducer.PFTauProducer          = cms.InputTag("shrinkingTightConePFTauProducer")
shrinkingTightConePFTauDecayModeIndexProducer.PFTauDecayModeProducer = cms.InputTag("shrinkingTightConePFTauDecayModeProducer")
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByIsolation_cfi                      import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByLeadingTrackFinding_cfi            import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByLeadingTrackPtCut_cfi              import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByTrackIsolation_cfi                 import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByECALIsolation_cfi                  import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationAgainstElectron_cfi                  import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationAgainstMuon_cfi                      import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByIsolationUsingLeadingPion_cfi      import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByLeadingPionPtCut_cfi               import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByTrackIsolationUsingLeadingPion_cfi import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByECALIsolationUsingLeadingPion_cfi  import *
shrinkingTightConePFTauDiscriminationByLeadingTrackFinding                          = copy.deepcopy(pfRecoTauDiscriminationByLeadingTrackFinding)
shrinkingTightConePFTauDiscriminationByLeadingTrackFinding.PFTauProducer            = 'shrinkingTightConePFTauProducer'
shrinkingTightConePFTauDiscriminationByLeadingTrackPtCut                            = copy.deepcopy(pfRecoTauDiscriminationByLeadingTrackPtCut)
shrinkingTightConePFTauDiscriminationByLeadingTrackPtCut.PFTauProducer              = 'shrinkingTightConePFTauProducer'
shrinkingTightConePFTauDiscriminationByLeadingPionPtCut                             = copy.deepcopy(pfRecoTauDiscriminationByLeadingPionPtCut)
shrinkingTightConePFTauDiscriminationByLeadingPionPtCut.PFTauProducer               = 'shrinkingTightConePFTauProducer'
shrinkingTightConePFTauDiscriminationByIsolation                                    = copy.deepcopy(pfRecoTauDiscriminationByIsolation)
shrinkingTightConePFTauDiscriminationByIsolation.PFTauProducer                      = 'shrinkingTightConePFTauProducer'
shrinkingTightConePFTauDiscriminationByTrackIsolation                               = copy.deepcopy(pfRecoTauDiscriminationByTrackIsolation)
shrinkingTightConePFTauDiscriminationByTrackIsolation.PFTauProducer                 = 'shrinkingTightConePFTauProducer'
shrinkingTightConePFTauDiscriminationByECALIsolation                                = copy.deepcopy(pfRecoTauDiscriminationByECALIsolation)
shrinkingTightConePFTauDiscriminationByECALIsolation.PFTauProducer                  = 'shrinkingTightConePFTauProducer'
shrinkingTightConePFTauDiscriminationByIsolationUsingLeadingPion                    = copy.deepcopy(pfRecoTauDiscriminationByIsolationUsingLeadingPion)
shrinkingTightConePFTauDiscriminationByIsolationUsingLeadingPion.PFTauProducer      = 'shrinkingTightConePFTauProducer'
shrinkingTightConePFTauDiscriminationByTrackIsolationUsingLeadingPion               = copy.deepcopy(pfRecoTauDiscriminationByTrackIsolationUsingLeadingPion)
shrinkingTightConePFTauDiscriminationByTrackIsolationUsingLeadingPion.PFTauProducer = 'shrinkingTightConePFTauProducer'
shrinkingTightConePFTauDiscriminationByECALIsolationUsingLeadingPion                = copy.deepcopy(pfRecoTauDiscriminationByECALIsolationUsingLeadingPion)
shrinkingTightConePFTauDiscriminationByECALIsolationUsingLeadingPion.PFTauProducer  = 'shrinkingTightConePFTauProducer'
shrinkingTightConePFTauDiscriminationAgainstElectron                                = copy.deepcopy(pfRecoTauDiscriminationAgainstElectron)
shrinkingTightConePFTauDiscriminationAgainstElectron.PFTauProducer                  = 'shrinkingTightConePFTauProducer'
shrinkingTightConePFTauDiscriminationAgainstMuon                                    = copy.deepcopy(pfRecoTauDiscriminationAgainstMuon)
shrinkingTightConePFTauDiscriminationAgainstMuon.PFTauProducer                      = 'shrinkingTightConePFTauProducer'

# gen met
from RecoMET.Configuration.GenMETParticles_cff import *
from RecoMET.METProducers.genMetTrue_cfi import *
#from RecoMET.METProducers.genMet_cfi import *
from RecoMET.METProducers.genMetCalo_cfi import *

# define sequence for gen jet production
from PhysicsTools.JetMCAlgos.TauGenJets_cfi import *

producePatTuple = cms.Sequence(
    tauGenJets *
    electronIdCutBased *
    recoElectronIsolation *
    ic5PFJetTracksAssociatorAtVertex *
    pfRecoTauTagInfoProducer *
    produceAndDiscriminateShrinkingConePFTaus *
    produceAndDiscriminateFixedConeHighEffPFTaus * 
    produceAndDiscriminateFixedConePFTaus *
    shrinkingTightConePFTauProducer*
    shrinkingTightConePFTauDecayModeProducer*
    shrinkingTightConePFTauDecayModeIndexProducer*
    shrinkingTightConePFTauDiscriminationByLeadingTrackFinding*
    shrinkingTightConePFTauDiscriminationByLeadingTrackPtCut*
    shrinkingTightConePFTauDiscriminationByLeadingPionPtCut*
    shrinkingTightConePFTauDiscriminationByIsolation*
    shrinkingTightConePFTauDiscriminationByTrackIsolation*
    shrinkingTightConePFTauDiscriminationByECALIsolation*
    shrinkingTightConePFTauDiscriminationByIsolationUsingLeadingPion*
    shrinkingTightConePFTauDiscriminationByTrackIsolationUsingLeadingPion*
    shrinkingTightConePFTauDiscriminationByECALIsolationUsingLeadingPion*
    shrinkingTightConePFTauDiscriminationAgainstElectron*
    shrinkingTightConePFTauDiscriminationAgainstMuon *
    patDefaultSequence *
    genMETParticles *
#    #genMet *
    genMetCalo *
    genMetTrue
    )


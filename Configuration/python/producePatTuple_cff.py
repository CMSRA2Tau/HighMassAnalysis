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

ak5PFJetTracksAssociatorAtVertex.j2tParametersVX = cms.PSet(
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(1.0)
)

# tau collection w/ shrinking cone of 3/ET
fixedConeHighEffPFTauProducer = copy.deepcopy(pfRecoTauProducer)
fixedConeHighEffPFTauProducer.LeadPFCand_minPt      = cms.double(5.0)
#--- tracker signal
fixedConeHighEffPFTauProducer.TrackerSignalConeSizeFormula = cms.string('0.15') ## **
fixedConeHighEffPFTauProducer.TrackerSignalConeSize_min = cms.double(0.07)
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
fixedConeHighEffPFTauProducer.TrackerIsolConeSizeFormula   = cms.string('1.00') ## **
fixedConeHighEffPFTauProducer.TrackerIsolConeSize_min   = cms.double(0.0)
fixedConeHighEffPFTauProducer.TrackerIsolConeSize_max   = cms.double(1.1)
#--- ecal isolation
fixedConeHighEffPFTauProducer.ECALIsolConeSizeFormula      = cms.string('1.00') ## **
fixedConeHighEffPFTauProducer.ECALIsolConeSize_min      = cms.double(0.0)
fixedConeHighEffPFTauProducer.ECALIsolConeSize_max      = cms.double(1.1)
#--- hcal isolation
fixedConeHighEffPFTauProducer.HCALIsolConeSizeFormula      = cms.string('1.00') ## **
fixedConeHighEffPFTauProducer.HCALIsolConeSize_min      = cms.double(0.0)
fixedConeHighEffPFTauProducer.HCALIsolConeSize_max      = cms.double(1.1)

from RecoTauTag.RecoTau.PFRecoTauDecayModeDeterminator_cfi                          import *
fixedConeHighEffPFTauDecayModeProducer               = copy.deepcopy(pfTauDecayMode)
fixedConeHighEffPFTauDecayModeProducer.PFTauProducer = 'fixedConeHighEffPFTauProducer'
from RecoTauTag.RecoTau.PFRecoTauDecayModeIndexProducer_cfi                             import *
fixedConeHighEffPFTauDecayModeIndexProducer                        = copy.deepcopy(pfTauDecayModeIndexProducer)
fixedConeHighEffPFTauDecayModeIndexProducer.PFTauProducer          = cms.InputTag("fixedConeHighEffPFTauProducer")
fixedConeHighEffPFTauDecayModeIndexProducer.PFTauDecayModeProducer = cms.InputTag("fixedConeHighEffPFTauDecayModeProducer")
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
from RecoTauTag.RecoTau.TauDiscriminatorTools import *
fixedConeHighEffPFTauDiscriminationByLeadingTrackFinding                          = copy.deepcopy(pfRecoTauDiscriminationByLeadingTrackFinding)
setTauSource(fixedConeHighEffPFTauDiscriminationByLeadingTrackFinding, 'fixedConeHighEffPFTauProducer')
fixedConeHighEffPFTauDiscriminationByLeadingTrackPtCut                            = copy.deepcopy(pfRecoTauDiscriminationByLeadingTrackPtCut)
setTauSource(fixedConeHighEffPFTauDiscriminationByLeadingTrackPtCut, 'fixedConeHighEffPFTauProducer')
fixedConeHighEffPFTauDiscriminationByLeadingPionPtCut                             = copy.deepcopy(pfRecoTauDiscriminationByLeadingPionPtCut)
setTauSource(fixedConeHighEffPFTauDiscriminationByLeadingPionPtCut, 'fixedConeHighEffPFTauProducer')
fixedConeHighEffPFTauDiscriminationByIsolation                                    = copy.deepcopy(pfRecoTauDiscriminationByIsolation)
setTauSource(fixedConeHighEffPFTauDiscriminationByIsolation, 'fixedConeHighEffPFTauProducer')
fixedConeHighEffPFTauDiscriminationByTrackIsolation                               = copy.deepcopy(pfRecoTauDiscriminationByTrackIsolation)
setTauSource(fixedConeHighEffPFTauDiscriminationByTrackIsolation, 'fixedConeHighEffPFTauProducer')
fixedConeHighEffPFTauDiscriminationByECALIsolation                                = copy.deepcopy(pfRecoTauDiscriminationByECALIsolation)
setTauSource(fixedConeHighEffPFTauDiscriminationByECALIsolation, 'fixedConeHighEffPFTauProducer')
fixedConeHighEffPFTauDiscriminationByIsolationUsingLeadingPion                    = copy.deepcopy(pfRecoTauDiscriminationByIsolationUsingLeadingPion)
setTauSource(fixedConeHighEffPFTauDiscriminationByIsolationUsingLeadingPion, 'fixedConeHighEffPFTauProducer')
fixedConeHighEffPFTauDiscriminationByTrackIsolationUsingLeadingPion               = copy.deepcopy(pfRecoTauDiscriminationByTrackIsolationUsingLeadingPion)
setTauSource(fixedConeHighEffPFTauDiscriminationByTrackIsolationUsingLeadingPion, 'fixedConeHighEffPFTauProducer')
fixedConeHighEffPFTauDiscriminationByECALIsolationUsingLeadingPion                = copy.deepcopy(pfRecoTauDiscriminationByECALIsolationUsingLeadingPion)
setTauSource(fixedConeHighEffPFTauDiscriminationByECALIsolationUsingLeadingPion, 'fixedConeHighEffPFTauProducer')
fixedConeHighEffPFTauDiscriminationAgainstElectron                                = copy.deepcopy(pfRecoTauDiscriminationAgainstElectron)
setTauSource(fixedConeHighEffPFTauDiscriminationAgainstElectron, 'fixedConeHighEffPFTauProducer')
fixedConeHighEffPFTauDiscriminationAgainstMuon                                    = copy.deepcopy(pfRecoTauDiscriminationAgainstMuon)
setTauSource(fixedConeHighEffPFTauDiscriminationAgainstMuon, 'fixedConeHighEffPFTauProducer')

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
from RecoTauTag.RecoTau.TauDiscriminatorTools import *
shrinkingTightConePFTauDiscriminationByLeadingTrackFinding                          = copy.deepcopy(pfRecoTauDiscriminationByLeadingTrackFinding)
setTauSource(shrinkingTightConePFTauDiscriminationByLeadingTrackFinding, 'shrinkingTightConePFTauProducer')
shrinkingTightConePFTauDiscriminationByLeadingTrackPtCut                            = copy.deepcopy(pfRecoTauDiscriminationByLeadingTrackPtCut)
setTauSource(shrinkingTightConePFTauDiscriminationByLeadingTrackPtCut, 'shrinkingTightConePFTauProducer')
shrinkingTightConePFTauDiscriminationByLeadingPionPtCut                             = copy.deepcopy(pfRecoTauDiscriminationByLeadingPionPtCut)
setTauSource(shrinkingTightConePFTauDiscriminationByLeadingPionPtCut, 'shrinkingTightConePFTauProducer')
shrinkingTightConePFTauDiscriminationByIsolation                                    = copy.deepcopy(pfRecoTauDiscriminationByIsolation)
setTauSource(shrinkingTightConePFTauDiscriminationByIsolation, 'shrinkingTightConePFTauProducer')
shrinkingTightConePFTauDiscriminationByTrackIsolation                               = copy.deepcopy(pfRecoTauDiscriminationByTrackIsolation)
setTauSource(shrinkingTightConePFTauDiscriminationByTrackIsolation, 'shrinkingTightConePFTauProducer')
shrinkingTightConePFTauDiscriminationByECALIsolation                                = copy.deepcopy(pfRecoTauDiscriminationByECALIsolation)
setTauSource(shrinkingTightConePFTauDiscriminationByECALIsolation, 'shrinkingTightConePFTauProducer')
shrinkingTightConePFTauDiscriminationByIsolationUsingLeadingPion                    = copy.deepcopy(pfRecoTauDiscriminationByIsolationUsingLeadingPion)
setTauSource(shrinkingTightConePFTauDiscriminationByIsolationUsingLeadingPion, 'shrinkingTightConePFTauProducer')
shrinkingTightConePFTauDiscriminationByTrackIsolationUsingLeadingPion               = copy.deepcopy(pfRecoTauDiscriminationByTrackIsolationUsingLeadingPion)
setTauSource(shrinkingTightConePFTauDiscriminationByTrackIsolationUsingLeadingPion, 'shrinkingTightConePFTauProducer')
shrinkingTightConePFTauDiscriminationByECALIsolationUsingLeadingPion                = copy.deepcopy(pfRecoTauDiscriminationByECALIsolationUsingLeadingPion)
setTauSource(shrinkingTightConePFTauDiscriminationByECALIsolationUsingLeadingPion, 'shrinkingTightConePFTauProducer')
shrinkingTightConePFTauDiscriminationAgainstElectron                                = copy.deepcopy(pfRecoTauDiscriminationAgainstElectron)
setTauSource(shrinkingTightConePFTauDiscriminationAgainstElectron, 'shrinkingTightConePFTauProducer')
shrinkingTightConePFTauDiscriminationAgainstMuon                                    = copy.deepcopy(pfRecoTauDiscriminationAgainstMuon)
setTauSource(shrinkingTightConePFTauDiscriminationAgainstMuon, 'shrinkingTightConePFTauProducer')

# gen met
from RecoMET.Configuration.GenMETParticles_cff import *
from RecoMET.METProducers.genMetTrue_cfi import *
#from RecoMET.METProducers.genMet_cfi import *
from RecoMET.METProducers.genMetCalo_cfi import *

# define sequence for gen jet production
from PhysicsTools.JetMCAlgos.TauGenJets_cfi import *

produceAndDiscriminateShrinkingConePFTausCustomized = cms.Sequence(
      shrinkingConePFTauProducer*
      shrinkingConePFTauDiscriminationByLeadingTrackFinding*
      shrinkingConePFTauDiscriminationByLeadingTrackPtCut*
      shrinkingConePFTauDiscriminationByLeadingPionPtCut*
      shrinkingConePFTauDiscriminationByIsolation*
      shrinkingConePFTauDiscriminationByTrackIsolation*
      shrinkingConePFTauDiscriminationByECALIsolation*
      shrinkingConePFTauDiscriminationByIsolationUsingLeadingPion*
      shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion*
      shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion*
      shrinkingConePFTauDiscriminationAgainstElectron*
      shrinkingConePFTauDiscriminationAgainstMuon
)

produceThePatCands = cms.Sequence(
    recoElectronIsolation *
    ic5PFJetTracksAssociatorAtVertex *
    ak5PFJetTracksAssociatorAtVertex *
    pfRecoTauTagInfoProducer *
    produceAndDiscriminateShrinkingConePFTaus *
    produceShrinkingConeDiscriminationByTauNeuralClassifier *
    fixedConeHighEffPFTauProducer*
    fixedConeHighEffPFTauDecayModeProducer*
    fixedConeHighEffPFTauDecayModeIndexProducer*
    fixedConeHighEffPFTauDiscriminationByLeadingTrackFinding*
    fixedConeHighEffPFTauDiscriminationByLeadingTrackPtCut*
    fixedConeHighEffPFTauDiscriminationByLeadingPionPtCut*
    fixedConeHighEffPFTauDiscriminationByIsolation*
    fixedConeHighEffPFTauDiscriminationByTrackIsolation*
    fixedConeHighEffPFTauDiscriminationByECALIsolation*
    fixedConeHighEffPFTauDiscriminationByIsolationUsingLeadingPion*
    fixedConeHighEffPFTauDiscriminationByTrackIsolationUsingLeadingPion*
    fixedConeHighEffPFTauDiscriminationByECALIsolationUsingLeadingPion*
    fixedConeHighEffPFTauDiscriminationAgainstElectron*
    fixedConeHighEffPFTauDiscriminationAgainstMuon *    
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
    patCustomizedCandidates *
    selectedPatCustomizedCandidates
)

produceGenMETInfo = cms.Sequence(
    genMETParticles *
    genMetCalo *
    genMetTrue
)

if(data):
  producePatTuple = cms.Sequence(
    produceThePatCands
  )
else:
  producePatTuple = cms.Sequence(
    produceThePatCands
    + produceGenMETInfo
  )



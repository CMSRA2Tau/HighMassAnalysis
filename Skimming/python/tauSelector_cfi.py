import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# select CaloTaus
#--------------------------------------------------------------------------------

selectedCaloTaus = cms.EDFilter("CaloTauSelector",
  src = cms.InputTag('caloRecoTauProducer'),
  discriminators = cms.VPSet(
    cms.PSet(
      discriminator = cms.InputTag("caloRecoTauDiscriminationByLeadingTrackPtCut"),
      selectionCut = cms.double(0.5)
    )
  ),
  filter = cms.bool(True)
)

#--------------------------------------------------------------------------------
# select PFTaus
#--------------------------------------------------------------------------------

selectedPFTaus = cms.EDFilter("PFTauSelector",
  src = cms.InputTag('shrinkingConePFTauProducer'),
  discriminators = cms.VPSet(
    cms.PSet(
      discriminator = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingPionPtCut"),
      selectionCut = cms.double(0.5)
    )
  ),
  filter = cms.bool(True)
)


import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('HiMassTau')

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( -1 )
)
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
'/store/user/eluiggi/zprimeTauTau500_7TeV_STARTUP31X_V4_GEN-SIM-RAW/ZprimeTauTau500PatSummer10Start36V10/a93249a823fb0dc819ca18ebc3447f9a/zprimeTauTauPat_8_1_RCT.root'
#'/store/user/florez/Zprime500_PAT_STARTUP_36X/florez/zprimeTauTau500_7TeV_STARTUP31X_V4_GEN-SIM-RAW/Zprime500_7TeV_PAT_eTau_STARTUP_36X/8b0dd7d13c5c85a732deb529e37114fd/eTauPatTuple_9_1_8vh.root',

    )
)

process.TFileService = cms.Service("TFileService", 
#    fileName = cms.string("outputFILENAME")
    fileName = cms.string("muTauAnalysis.root")
)

process.analyzeHiMassTau = cms.EDAnalyzer('HiMassTauAnalysis',

    #-----Generator level Inputs
    GenParticleSource = cms.untracked.InputTag('genParticles'),				# gen particle collection

    #-----Inputs to determine which channel to analyze
    AnalyzeTauForLeg1		= cms.bool(False),					# if true, taus will be used for leg1
    AnalyzeMuonForLeg1		= cms.bool(False),					# if true, muons will be used for leg1
    AnalyzeElectronForLeg1	= cms.bool(True),					# if true, electrons will be used for leg1
    AnalyzeTauForLeg2		= cms.bool(True),					# if true, taus will be used for leg2
    AnalyzeMuonForLeg2		= cms.bool(False),					# if true, muons will be used for leg2
    AnalyzeElectronForLeg2	= cms.bool(False),					# if true, electrons will be used for leg2

    #-----Reco Tau Inputs
    RecoTauSource = cms.InputTag('selectedLayer1FixedConePFTaus'),			# other choices include:
											# selectedLayer1FixedConeHighEffPFTaus
											# selectedLayer1ShrinkingConeHighEffPFTaus
											# selectedLayer1ShrinkingConePFTaus
    RecoTauEtaCut = cms.double(999.1),							# require tau |eta|<=X
    RecoTauPtMinCut = cms.double(0.),							# require tau pt>=X
    RecoTauPtMaxCut = cms.double(9999.),						# require tau pt<=X
    DoRecoTauDiscrByLeadTrack = cms.bool(False),					# if true, tau is required to pass a lead track pt cut
    UseRecoTauDiscrByLeadTrackFlag = cms.bool(True), 					# if true, default seed track discriminator is used
                                                     					# if false, seed track cut will be recalculated using the parameters below
    RecoTauDiscrByLeadTrack = cms.untracked.string('leadingTrackPtCut'),		# name of the lead track discriminator flag
    DoRecoTauDiscrByLeadTrackNhits = cms.bool(False),					# if true, tau leading track is required to have >= X hits
    RecoTauLeadTrackMinHits = cms.int32(10),						# tau leading track hits >= X
    DoRecoTauDiscrByIsolation = cms.bool(False), 					# if true, isolation will be applied
    UseRecoTauDiscrByIsolationFlag = cms.bool(False), 					# if true, the default isolation discriminator is used
                                                      					# if false, isolation is recalculated using the parameters below
    RecoTauDiscrByIsolation = cms.untracked.string('byIsolation'),			# name of the isolation discriminator flag
    UseRecoTauIsoSumPtInsteadOfNiso = cms.bool(True),					# if true, sum pt is used for tau isolation instead
											# of the number of isolation candidates
    UseRecoTauEllipseForEcalIso = cms.bool(False),					# if true, an ellipse in eta-phi space will be used to define
											# the signal to isolation annulus for ECAL isolation
    RecoTauEcalIsoRphiForEllipse = cms.double(0.15),					# a:  dphi^2 / a^2 + deta^2 / b^2  (ECAL ellipse)
    RecoTauEcalIsoRetaForEllipse = cms.double(0.07),					# b:  dphi^2 / a^2 + deta^2 / b^2  (ECAL ellipse)
    RecoTauNisoMax = cms.int32(0),							# number of isolation candidates <=X
    RecoTauIsoSumPtMaxCutValue = cms.double(6.0),					# sum pt of tracks & gammas < X
    RecoTauIsoSumPtMinCutValue = cms.double(0.0),					# sum pt of tracks & gammas >= X
    RecoTauDiscrByProngType = cms.string('NONE'),		 			# if '1or3', taus will be required to have 1 || 3 prongs
											# if '1', taus will be required to have 1 prongs
											# if '3', taus will be required to have 3 prongs
											# if anything else is used, no prong cuts are applied
    DoRecoTauDiscrBySignalTracksAndGammasMass = cms.bool(False),
    RecoTauSignal3ProngAndGammasMassMinCutValue = cms.double(0.9),
    RecoTauSignal3ProngAndGammasMassMaxCutValue = cms.double(999.0),
    RecoTauSignal1ProngAndGammasMassForPionMinCutValue = cms.double(0.08),
    RecoTauSignal1ProngAndGammasMassForPionMaxCutValue = cms.double(0.2),
    RecoTauSignal1ProngAndGammasMassForKaonVetoMinCutValue = cms.double(0.5),
    RecoTauSignal1ProngAndGammasMassForKaonVetoMaxCutValue = cms.double(999.0),
    RecoTauLeadTrackThreshold = cms.double(5.0), 					# seed track pt > X
    RecoTauSigGamThreshold = cms.double(5.0),						# gamma pt > X
    RecoTauIsoDeltaRCone = cms.double(0.5),	 					# tau outer isolation conesize
    RecoTauTrackIsoTrkThreshold = cms.double(1.0), 					# min pt requirement for isolation tracks
    RecoTauGammaIsoGamThreshold = cms.double(1.5), 					# min pt requirement for isolation gammas
    DoRecoTauDiscrAgainstElectron = cms.bool(False),					# if true, electron veto will be applied
    RecoTauDiscrAgainstElectron = cms.untracked.string('againstElectron'),		# name of electron veto discriminator flag
    DoRecoTauDiscrByCrackCut = cms.bool(False),						# if true, taus that fall on the cracks will not be considered
    DoRecoTauDiscrAgainstMuon = cms.bool(False),					# if true, muon veto will be applied
    RecoTauDiscrAgainstMuon = cms.untracked.string('againstMuon'),			# name of muon veto discriminator flag
    SetTANC = cms.bool(False),								# set true if wanting to fill TanC info in the Ntuple

    #-----Reco Muon Inputs
    RecoMuonSource = cms.InputTag('selectedPatMuons'),				# muon collection
    RecoMuonEtaCut = cms.double(999.1),							# require muon |eta|<=X
    RecoMuonPtMinCut = cms.double(0.),							# require muon pt>=X
    RecoMuonPtMaxCut = cms.double(9999.),						# require muon pt<=X
    DoRecoMuonDiscrByGlobal = cms.bool(False),						# if true, muon will be required to be a 'global muon
    DoRecoMuonDiscrByIsolation = cms.bool(False),					# if true, muon isolation will be applied
    RecoMuonIsoSumPtMaxCutValue = cms.double(3.0),					# sum pt of isolation tracks & ecal rechits < X
    RecoMuonIsoSumPtMinCutValue = cms.double(0.0),					# sum pt of isolation tracks & ecal rechits >= X
    RecoMuonIsoDeltaRCone = cms.double(0.6),						# outer conesize used for isolation
    RecoMuonTrackIsoTrkThreshold = cms.double(0.5),					# isolation tracks are required to have pt>X
    RecoMuonEcalIsoRecHitThreshold = cms.double(0.3),					# isolation rechits are required to have pt>X
    DoRecoMuonDiscrByIp = cms.bool(False),						# if true, muon will be required to have |d0|<X
    RecoMuonIpCut = cms.double(0.028),							# |d0|<X
    DoRecoMuonDiscrByPionVeto = cms.bool(False),                                        # if true, muon will be required to pass pion veto cut
    RecoMuonCaloCompCoefficient = cms.double(0.8),                                      # a -> pion veto: a*caloComp + b*segmComp
    RecoMuonSegmCompCoefficient = cms.double(1.2),                                      # b -> pion veto: a*caloComp + b*segmComp
    RecoMuonAntiPionCut = cms.double(1.0),                                              # pion veto > X

    #-----Reco Electron Inputs
    RecoElectronSource = cms.InputTag('selectedPatElectrons'),			# electron collection
    RecoElectronEtaCut = cms.double(999.1),						# require electron |eta|<=X
    RecoElectronPtMinCut = cms.double(0.),						# require electron pt>=X
    RecoElectronPtMaxCut = cms.double(9999.),						# require electron pt<=X
    DoRecoElectronDiscrByTrackIsolation = cms.bool(False), 				# if true, electrons will be required to pass track isolation
    RecoElectronTrackIsoSumPtCutValue = cms.double(1.0), 				# sum pt of tracks < X
    RecoElectronTrackIsoDeltaRCone = cms.double(0.6), 					# isolation conesize used to calculate sum pt of tracks
    RecoElectronTrackIsoTrkThreshold = cms.double(1.0), 				# min pt requirement for isolation tracks
    DoRecoElectronDiscrByEcalIsolation = cms.bool(False), 				# if true, electrons will be required to pass ecal isolation
    RecoElectronEcalIsoSumPtCutValue = cms.double(1.0), 				# sum pt of ecal rechits < X
    RecoElectronEcalIsoDeltaRCone = cms.double(0.6), 					# isolation conesize used to calculate sum pt of ecal rechits
    RecoElectronEcalIsoRecHitThreshold = cms.double(1.0),				# min pt requirement for isolation rechits
    DoRecoElectronDiscrByIp = cms.bool(False), 						# if true, electron will be required to have |d0|<X
    RecoElectronIpCut = cms.double(999.), 						# |d0(impact parameter)| < X
    DoRecoElectronDiscrByEoverP = cms.bool(False), 					# if true, require the electron to have X < E/p < Y
    RecoElectronEoverPMax = cms.double(1.3), 						# E/p < Y
    RecoElectronEoverPMin = cms.double(0.7), 						# E/p > X
    DoRecoElectronDiscrByHoverEm = cms.bool(False), 					# if true, require the electron to have H/Em < X
    RecoElectronHoverEmCut = cms.double(0.008), 					# H/Em < X
    DoRecoElectronDiscrBySigmaEtaEta = cms.bool(False),					# at the moment, does nothing
    RecoElectronSigmaEtaEtaCut = cms.double(0.0),					# at the moment, does nothing
    DoRecoElectronDiscrBySigmaIEtaIEta = cms.bool(False),				# at the moment, does nothing
    RecoElectronSigmaIEtaIEtaCut = cms.double(0.0),					# at the moment, does nothing
    DoRecoElectronDiscrBySCE5by5 = cms.bool(False),					# at the moment, does nothing
    RecoElectronSCE5by5Cut = cms.double(0.0),						# at the moment, does nothing
    DoRecoElectronDiscrByEcalDrivenSeed = cms.bool(False),
    DoRecoElectronDiscrByTrackerDrivenSeed = cms.bool(False),

    #-----Reco Jet Inputs
    RecoJetSource                       = cms.InputTag('selectedPatJets'),           # jet collection
    RecoJetEtaMinCut                    = cms.double(0.0),                              # require jet |eta|>=X
    RecoJetEtaMaxCut                    = cms.double(999.5),                            # require jet |eta|<=X
    RecoJetPtCut                        = cms.double(0.0),                             # require jet pt>=X
    UseCorrectedJet                     = cms.bool(False),                              # if true, jet corrections are used
    RemoveJetOverlapWithMuons           = cms.bool(False),                               # if true, jets w/ dR(muon,jet)<X will not be considered
    JetMuonMatchingDeltaR               = cms.double(0.5),                              # dR(muon,jet)<X used for removing jets from the "good jet" list
    RemoveJetOverlapWithElectrons       = cms.bool(False),                              # if true, jets w/ dR(electron,jet)<X will not be considered
    JetElectronMatchingDeltaR           = cms.double(0.5),                              # dR(electron,jet)<X used for removing jets from the "good jet" list
    RemoveJetOverlapWithTaus            = cms.bool(False),                               # if true, jets w/ dR(tau,jet)<X will not be considered
    JetTauMatchingDeltaR                = cms.double(0.25),                             # dR(tau,jet)<X used for removing jets from the "good jet" list

    #-----Vertex Inputs
    RecoVertexSource = cms.InputTag('offlinePrimaryVertices'), 				# vertex collection
    RecoVertexMaxZposition = cms.double(999.0),						# vertex |z| < X
    RecoVertexMinTracks = cms.int32(-1),						# vertex must have >= 2 "good" tracks used to reconstruct it
    RecoVertexTrackWeight = cms.double(0.5),						# weight used to define "good" tracks used to reconstruct vertex

    #-----Trigger Inputs
    RecoTriggerSource = cms.InputTag("TriggerResults","","REDIGI36X"),			# trigger collection
    TriggerRequirements = cms.vstring('HLT_Mu9'),					# trigger path name

    #-----Topology Inputs
    RecoMetSource = cms.InputTag('patMETsPF'),					# met collection
											# particle flow met for 2X	= layer1PFMETs
											# particle flow met for 3X	= layer1METsPF
											# standard calo met for 2X & 3X	= layer1METs
    DoDiscrByMet = cms.bool(False), 							# if true, met will be required to be > X
    RecoMetCut = cms.double(30.0), 							# met > X
    DoDiTauDiscrByDeltaR = cms.bool(False), 						# if true, ditau pairs must have dR(leg1,leg2) > X
    DiTauDeltaRCut = cms.double(0.7),	 						# dR(leg1,leg2) > X
    DiTauDiscrByOSLSType = cms.string('NONE'),		 				# if 'OS', product of leg1 charge and leg2 charge < 0
											# if 'LS', product of leg1 charge and leg2 charge > 0
											# if anything else is used, no OSLS requirement is applied
    UseTauSeedTrackForDiTauDiscrByOSLS = cms.bool(False), 				# if true, then use q(tau) = q(tau seed track)
                                                          				# if false, then use q(tau) = sum of q(tau jet tracks)
    DoDiTauDiscrByCosDphi = cms.bool(False), 						# if true, require X <= cosine dphi(leg1,leg2) <= Y
    DiTauCosDphiMaxCut = cms.double(-0.95), 						# cosine dphi(leg1, leg2) <= Y
    DiTauCosDphiMinCut = cms.double(-1.00), 						# cosine dphi(leg1, leg2) >= X
    DoDiscrByMassReco = cms.bool(False), 						# if true, apply a mass cut on leg1+leg2+met combinations
    UseVectorSumOfVisProductsAndMetMassReco = cms.bool(False), 				# if true, calculate invariant mass using
                                                               				# vector sum of leg1, leg2, and met
    UseCollinerApproxMassReco = cms.bool(False), 					# if true, use the collinear approximation for the invariant mass
                                                 					# if both are set to false, then use only the visible products
											# to calculate the mass (leg1 + leg2)
    MassMinCut = cms.double(150.0), 							# mass > X
    MassMaxCut = cms.double(1000.0), 							# mass < Y
    DoDiTauDiscrByCDFzeta2D = cms.bool(False),						# if true, apply 2D zeta cut ( a*pzeta + b*pzetaVis > X )
    PZetaCutCoefficient = cms.double(1.0),						# a -- ( a*pzeta + b*pzetaVis > X )
    PZetaVisCutCoefficient = cms.double(-0.875),					# b -- ( a*pzeta + b*pzetaVis > X )
    CDFzeta2DCutValue = cms.double(-7.00),						# X -- ( a*pzeta + b*pzetaVis > X )
    DoDiTauDiscrByDeltaPtDivSumPt = cms.bool(False),					# if true, apply cut on deltapt(leg2,leg1) / sumpt(leg1,leg2)
    DiTauDeltaPtDivSumPtMinCutValue = cms.double(0.1),					# ( pt_leg2 - pt_leg1 ) / ( pt_leg2 + pt_leg1 ) >= X
    DiTauDeltaPtDivSumPtMaxCutValue = cms.double(1.0),					# ( pt_leg2 - pt_leg1 ) / ( pt_leg2 + pt_leg1 ) <= X
    DoDiscrByLeg1MetDphi = cms.bool(False),						# if true, require X <= cosine dphi(leg1,met) <= Y
    Leg1MetDphiMinCut = cms.double(1.30),						# cosine dphi(leg1, met) >= X
    Leg1MetDphiMaxCut = cms.double(3.15),						# cosine dphi(leg1, met) <= Y
    DoDiscrByLeg2MetDphi = cms.bool(False),						# if true, require X <= cosine dphi(leg2,met) <= Y
    Leg2MetDphiMinCut = cms.double(1.30),						# cosine dphi(leg2, met) >= X
    Leg2MetDphiMaxCut = cms.double(3.15),						# cosine dphi(leg2, met) <= Y
    DoTauDiscrByIsZeeCut = cms.bool(False),

    #-----do matching to gen?
    MatchLeptonToGen = cms.bool(True),							# if true, match reco lepton to a gen lepton
    UseLeptonMotherId = cms.bool(True),							# if true, require the matched lepton to come from a certain 
											# 'mother' particle
    UseLeptonGrandMotherId = cms.bool(True),						# if true, require the matched lepton to come from a certain
											# 'grandmother' particle
    LeptonMotherId = cms.int32(15),							# pdgId of the 'mother' particle
    LeptonGrandMotherId = cms.int32(32),						# pdgId of the 'grandmother' particle
    MatchTauToGen = cms.bool(True),							# if true, match reco tau to a gen had tau
    UseTauMotherId = cms.bool(True),							# if true, require the matched tau to come from a certain
											# 'mother' particle ('mother' here is NOT 15!!!!! Matching
											# for the had tau leg already requires the vis had tau to come
											# from a tau lepton)
    UseTauGrandMotherId = cms.bool(False),						# if true, require the matched tau to come from a certain
											# 'grandmother' particle
    TauMotherId = cms.int32(32),							# pdgId of the 'mother' particle
    TauGrandMotherId = cms.int32(1),							# pdgId of the 'grandmother' particle
    TauToGenMatchingDeltaR = cms.double(0.25),						# matching dR:  dR(vis gen tau,reco tau)<X

    #-----ntuple Inputs
    DoProduceNtuple = cms.bool(True),							# if true, Ntuple will be filled
    NtupleTreeName = cms.untracked.string('HMTTree'),					# name of the Ntuple tree

    #-----Fill Histograms? Histograms are filled for events passing the specified cuts
    FillRecoVertexHists = cms.bool(True),						# if true, fill histograms for vertices
    FillGenTauHists = cms.bool(True),							# if true, fill histograms for gen had taus
    FillRecoTauHists = cms.bool(True),							# if true, fill histograms for reco taus
    FillRecoMuonHists = cms.bool(False),							# if true, fill histograms for reco muons
    FillRecoElectronHists = cms.bool(True),						# if true, fill histograms for reco electrons
    FillRecoJetHists = cms.bool(True),							# if true, fill histograms for reco jets
    FillTopologyHists = cms.bool(True),							# if true, fill topology histograms (e.g. met, mass, ...)

    #-----Event Sequence inputs
    RecoTriggersNmin = cms.int32(0),							# require event to pass >=X trigger paths defined above
    RecoVertexNmin = cms.int32(0),							# require event to have >=X vertices passing specified cuts
    RecoVertexNmax = cms.int32(1000),							# require event to have <=X vertices passing specified cuts
    RecoLeg1Nmin = cms.int32(0),							# require event to have >=X leg1 objects passing specified cuts
    RecoLeg1Nmax = cms.int32(1000),							# require event to have <=X leg1 objects passing specified cuts
    RecoLeg2Nmin = cms.int32(0),							# require event to have >=X leg2 objects passing specified cuts
    RecoLeg2Nmax = cms.int32(1000),							# require event to have <=X leg2 objects passing specified cuts
    RecoJetNmin = cms.int32(0),								# require event to have >=X "jets" passing specified cuts
    RecoJetNmax = cms.int32(1000),							# require event to have <=X "jets" passing specified cuts
    CombinationsNmin = cms.int32(0),							# require event to have >=X leg1+leg2+met combinations 
											# passing specified cuts
    CombinationsNmax = cms.int32(1000),							# require event to have <=X leg1+leg2+met combinations
											# passing specified cuts

    EventSelectionSequence = cms.vstring('RecoTriggersNmin','RecoVertexNmin','RecoVertexNmax','RecoLeg1Nmin','RecoLeg1Nmax','RecoLeg2Nmin','RecoLeg2Nmax',
                                         'RecoJetNmin','RecoJetNmax','CombinationsNmin','CombinationsNmax'),

    #-----Inputs for systematic uncertainties
    CalculatePdfSystematicUncertanties = cms.bool(False),				# if true, pdf systematic uncertanties will be calculated
    PdfWeightTags = cms.untracked.VInputTag("cteq65PdfWeights"),			# collection of weights for systematic uncertanties
    #PdfWeightTags = cms.untracked.VInputTag("cteq65PdfWeights", "MRST2006nnloPdfWeights", "MRST2007lomodPdfWeights"),
    CalculateFSRSystematics = cms.bool(False),
    CalculateISRGluonSystematics = cms.bool(False),
    CalculateISRGammaSystematics = cms.bool(False),
    SmearTheMuon = cms.bool(False),
    MuonPtScaleOffset = cms.double(1.0),
    MuonPtSigmaOffset = cms.double(1.0),
    MuonEtaScaleOffset = cms.double(1.0),
    MuonEtaSigmaOffset = cms.double(1.0),
    MuonPhiScaleOffset = cms.double(1.0),
    MuonPhiSigmaOffset = cms.double(1.0),
    SmearTheElectron = cms.bool(False),
    ElectronPtScaleOffset = cms.double(1.0),
    ElectronPtSigmaOffset = cms.double(1.0),
    ElectronEtaScaleOffset = cms.double(1.0),
    ElectronEtaSigmaOffset = cms.double(1.0),
    ElectronPhiScaleOffset = cms.double(1.0),
    ElectronPhiSigmaOffset = cms.double(1.0),
    SmearTheTau = cms.bool(False),
    TauPtScaleOffset = cms.double(1.0),
    TauPtSigmaOffset = cms.double(1.0),
    TauEtaScaleOffset = cms.double(1.0),
    TauEtaSigmaOffset = cms.double(1.0),
    TauPhiScaleOffset = cms.double(1.0),
    TauPhiSigmaOffset = cms.double(1.0),
    SmearTheJet = cms.bool(False),
    JetEnergyScaleOffset = cms.double(1.0),    
    SmearThePt = cms.bool(False),
    SmearTheEta = cms.bool(False),
    SmearThePhi = cms.bool(False),
    BosonPtBinEdges = cms.untracked.vdouble(
               0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.
            , 10., 11., 12., 13., 14., 15., 16., 17., 18., 19.
            , 20., 21., 22., 23., 24., 25., 26., 27., 28., 29.
            , 30., 31., 32., 33., 34., 35., 36., 37., 38., 39.
            , 40., 41., 42., 43., 44., 45., 46., 47., 48., 49.
            , 999999.
    ),
    PtWeights = cms.untracked.vdouble( 
              0.800665, 0.822121, 0.851249, 0.868285, 0.878733
            , 0.953853, 0.928108, 0.982021, 1.00659 , 1.00648
            , 1.03218 , 1.04924 , 1.03621 , 1.08743 , 1.01951
            , 1.10519 , 0.984263, 1.04853 , 1.06724 , 1.10183
            , 1.0503  , 1.13162 , 1.03837 , 1.12936 , 0.999173
            , 1.01453 , 1.11435 , 1.10545 , 1.07199 , 1.04542
            , 1.00828 , 1.0822  , 1.09667 , 1.16144 , 1.13906
            , 1.27974 , 1.14936 , 1.23235 , 1.06667 , 1.06363
            , 1.14225 , 1.22955 , 1.12674 , 1.03944 , 1.04639
            , 1.13667 , 1.20493 , 1.09349 , 1.2107  , 1.21073
    )
)

process.p = cms.Path(
    process.analyzeHiMassTau
    )

# print-out all python configuration parameter information
#print process.dumpPython()


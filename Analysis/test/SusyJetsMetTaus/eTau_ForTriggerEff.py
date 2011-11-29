import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('HiMassTau')

process.load('Configuration.StandardSequences.Services_cff')
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( -1 )
)
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
'/store/user/amarotta/amarotta/TTJets_TuneZ2_7TeV-madgraph-tauola/TTJets_mutauSkim_06222011/d2196243fc6117c733eff4717010f8b7/skimPat_9_1_RLi.root',
'/store/user/amarotta/amarotta/TTJets_TuneZ2_7TeV-madgraph-tauola/TTJets_mutauSkim_06222011/d2196243fc6117c733eff4717010f8b7/skimPat_99_1_IZ8.root',
) 
)

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("analysis.root")
)

process.analyzeHiMassTau = cms.EDAnalyzer('HiMassTauAnalysis')

    #-----Generator level Inputs
process.analyzeHiMassTau.GenParticleSource = cms.untracked.InputTag('genParticles')				# gen particle collection

    #-----Inputs to determine which channel to analyze
process.analyzeHiMassTau.AnalyzeTauForLeg1		= cms.bool(False)					# if true, taus will be used for leg1
process.analyzeHiMassTau.AnalyzeMuonForLeg1		= cms.bool(False)					# if true, muons will be used for leg1
process.analyzeHiMassTau.AnalyzeElectronForLeg1	= cms.bool(True)					# if true, electrons will be used for leg1
process.analyzeHiMassTau.AnalyzeTauForLeg2		= cms.bool(True)					# if true, taus will be used for leg2
process.analyzeHiMassTau.AnalyzeMuonForLeg2		= cms.bool(False)					# if true, muons will be used for leg2
process.analyzeHiMassTau.AnalyzeElectronForLeg2	= cms.bool(False)					# if true, electrons will be used for leg2
    #-----Reco Tau Inputs
process.analyzeHiMassTau.RecoTauSource = cms.InputTag('selectedLayer1HPSPFTaus')		# choices include:
											# selectedLayer1FixedConeHighEffPFTaus
											# selectedLayer1HPSPFTaus
											# selectedLayer1ShrinkingConePFTaus
											# selectedLayer1FixedConePFTaus
process.analyzeHiMassTau.RecoTauEtaCut = cms.double(2.1)							# require tau |eta|<=X
process.analyzeHiMassTau.RecoTauPtMinCut = cms.double(15.)							# require tau pt>=X
process.analyzeHiMassTau.RecoTauPtMaxCut = cms.double(9999.)						# require tau pt<=X
process.analyzeHiMassTau.DoRecoTauDiscrByLeadTrack = cms.bool(False) 					# if true, tau is required to pass a lead track pt cut
process.analyzeHiMassTau.UseRecoTauDiscrByLeadTrackFlag = cms.bool(True) 					# if true, default seed track discriminator is used
                                                     					# if false, seed track cut will be recalculated using the parameters below
process.analyzeHiMassTau.RecoTauDiscrByLeadTrack = cms.untracked.string('leadingPionPtCut')			# name of the lead track discriminator flag
process.analyzeHiMassTau.DoRecoTauDiscrByLeadTrackNhits = cms.bool(False)					# if true, tau leading track is required to have >= X hits
process.analyzeHiMassTau.RecoTauLeadTrackMinHits = cms.int32(12)						# tau leading track hits >= X
process.analyzeHiMassTau.DoRecoTauDiscrByH3x3OverP = cms.bool(False)					# if true, tau will be required
											# to pass H(3x3)/P(lead) cut
process.analyzeHiMassTau.RecoTauH3x3OverP = cms.double(0.03)						# H(3x3)/P(lead) > X
process.analyzeHiMassTau.DoRecoTauDiscrByIsolation = cms.bool(True) 					# if true, isolation will be applied
process.analyzeHiMassTau.UseRecoTauDiscrByIsolationFlag = cms.bool(True) 					# if true, the default isolation discriminator is used
                                                      					# if false, isolation is recalculated using the parameters below
process.analyzeHiMassTau.RecoTauDiscrByIsolation = cms.untracked.string('byVLooseIsolation')		# name of the isolation discriminator flag
process.analyzeHiMassTau.UseRecoTauIsoSumPtInsteadOfNiso = cms.bool(False)					# if true, sum pt is used for tau isolation instead
											# of the number of isolation candidates
process.analyzeHiMassTau.UseRecoTauEllipseForEcalIso = cms.bool(False)					# if true, an ellipse in eta-phi space will be used to define
											# the signal to isolation annulus for ECAL isolation
process.analyzeHiMassTau.RecoTauEcalIsoRphiForEllipse = cms.double(0.15)					# a:  dphi^2 / a^2 + deta^2 / b^2  (ECAL ellipse)
process.analyzeHiMassTau.RecoTauEcalIsoRetaForEllipse = cms.double(0.07)					# b:  dphi^2 / a^2 + deta^2 / b^2  (ECAL ellipse)
process.analyzeHiMassTau.RecoTauTrackNisoMax = cms.int32(0)							# number of isolation candidates <=X
process.analyzeHiMassTau.RecoTauEcalNisoMax = cms.int32(0)							# number of isolation candidates <=X
process.analyzeHiMassTau.RecoTauTrackIsoSumPtMaxCutValue = cms.double(1.0)					# sum pt of tracks & gammas < X
process.analyzeHiMassTau.RecoTauTrackIsoSumPtMinCutValue = cms.double(0.0)					# sum pt of tracks & gammas >= X
process.analyzeHiMassTau.RecoTauEcalIsoSumPtMaxCutValue = cms.double(1.0)					# sum pt of tracks & gammas < X
process.analyzeHiMassTau.RecoTauEcalIsoSumPtMinCutValue = cms.double(0.0)					# sum pt of tracks & gammas >= X
process.analyzeHiMassTau.RecoTauDiscrByProngType = cms.string('1or3hps')  		 			# if '1or3', taus will be required to have 1 || 3 prongs
											# if '1', taus will be required to have 1 prongs
											# if '3', taus will be required to have 3 prongs
											# if anything else is used, no prong cuts are applied
process.analyzeHiMassTau.RecoTauLeadTrackThreshold = cms.double(5.0) 					# seed track pt > X (used only
											# if default seed track discriminator
											# is NOT used
process.analyzeHiMassTau.RecoTauSigGamThreshold = cms.double(1.0)						# signal gammas pt > X
process.analyzeHiMassTau.RecoTauIsoDeltaRCone = cms.double(0.3)	 					# tau outer isolation conesize
process.analyzeHiMassTau.RecoTauTrackIsoTrkThreshold = cms.double(1.5) 					# min pt requirement for isolation tracks
process.analyzeHiMassTau.RecoTauGammaIsoGamThreshold = cms.double(2.0) 					# min pt requirement for isolation gammas
process.analyzeHiMassTau.DoRecoTauDiscrAgainstElectron = cms.bool(True)					# if true, default tau POG electron
											# discriminator veto will be applied
process.analyzeHiMassTau.RecoTauDiscrAgainstElectron = cms.untracked.string('againstTightElectron')		# name of electron veto discriminator flag
process.analyzeHiMassTau.DoRecoTauDiscrByCrackCut = cms.bool(False)						# if true, taus that fall on the cracks will not be considered
process.analyzeHiMassTau.DoRecoTauDiscrAgainstMuon = cms.bool(True) 					# if true, muon veto will be applied
process.analyzeHiMassTau.RecoTauDiscrAgainstMuon = cms.untracked.string('againstTightMuon')			# name of muon veto discriminator flag
process.analyzeHiMassTau.SelectTausThatAreMuons = cms.bool(False)
process.analyzeHiMassTau.SelectTausThatAreElectrons = cms.bool(False)

    #-----Reco Muon Inputs
process.analyzeHiMassTau.RecoMuonSource = cms.InputTag('selectedPatMuons')				        # muon collection
process.analyzeHiMassTau.RecoMuonEtaCut = cms.double(2.1)							# require muon |eta|<=X
process.analyzeHiMassTau.RecoMuonPtMinCut = cms.double(10.)							# require muon pt>=X
process.analyzeHiMassTau.RecoMuonPtMaxCut = cms.double(9999.)						# require muon pt<=X
process.analyzeHiMassTau.DoRecoMuonDiscrByGlobal = cms.bool(True)						# if true, muon will be required to be a 'global muon
process.analyzeHiMassTau.DoRecoMuonDiscrByIsolation = cms.bool(True)					# if true, muon isolation will be applied
process.analyzeHiMassTau.RecoMuonTrackIsoSumPtMaxCutValue = cms.double(1.0)					# sum pt of isolation tracks & ecal rechits < X
process.analyzeHiMassTau.RecoMuonTrackIsoSumPtMinCutValue = cms.double(0.0)					# sum pt of isolation tracks & ecal rechits >= X
process.analyzeHiMassTau.RecoMuonEcalIsoSumPtMaxCutValue = cms.double(1.0)					# sum pt of isolation tracks & ecal rechits < X
process.analyzeHiMassTau.RecoMuonEcalIsoSumPtMinCutValue = cms.double(0.0)					# sum pt of isolation tracks & ecal rechits >= X
process.analyzeHiMassTau.RecoMuonIsoDeltaRCone = cms.double(0.3)						# outer conesize used for isolation
process.analyzeHiMassTau.RecoMuonTrackIsoTrkThreshold = cms.double(1.5)					# isolation tracks are required to have pt>X
process.analyzeHiMassTau.RecoMuonEcalIsoRecHitThreshold = cms.double(2.0)					# isolation rechits are required to have pt>X
process.analyzeHiMassTau.DoRecoMuonDiscrByIp = cms.bool(True)						# if true, muon will be required to have |d0|<X
process.analyzeHiMassTau.RecoMuonIpCut = cms.double(0.2)							# |d0|<X
process.analyzeHiMassTau.DoRecoMuonDiscrByPionVeto = cms.bool(False)                                        # if true, muon will be required to pass pion veto cut
process.analyzeHiMassTau.RecoMuonCaloCompCoefficient = cms.double(0.8)                                      # a -> pion veto: a*caloComp + b*segmComp
process.analyzeHiMassTau.RecoMuonSegmCompCoefficient = cms.double(1.2)                                      # b -> pion veto: a*caloComp + b*segmComp
process.analyzeHiMassTau.RecoMuonAntiPionCut = cms.double(1.0)                                              # pion veto > X

    #-----Reco Electron Inputs
process.analyzeHiMassTau.RecoElectronSource = cms.InputTag('heepPatElectrons')			        # electron collection
process.analyzeHiMassTau.UseHeepInfo = cms.bool(True)							# if true, heep discriminators and
											# variables will be used instead of pat defaults
process.analyzeHiMassTau.RecoElectronEtaCut = cms.double(2.1)						# require electron |eta|<=X
process.analyzeHiMassTau.RecoElectronPtMinCut = cms.double(10.)						# require electron pt>=X
process.analyzeHiMassTau.RecoElectronPtMaxCut = cms.double(9999.)						# require electron pt<=X
process.analyzeHiMassTau.DoRecoElectronDiscrByTrackIsolation = cms.bool(True) 				# if true, electrons will be required to pass track isolation
process.analyzeHiMassTau.RecoElectronTrackIsoSumPtMaxCutValue = cms.double(1.0) 				# sum pt of tracks < X
process.analyzeHiMassTau.RecoElectronTrackIsoSumPtMinCutValue = cms.double(0.0) 				# sum pt of tracks >= X
process.analyzeHiMassTau.RecoElectronTrackIsoDeltaRCone = cms.double(0.4) 					# isolation conesize used to calculate sum pt of tracks
process.analyzeHiMassTau.RecoElectronTrackIsoTrkThreshold = cms.double(1.0) 				# min pt requirement for isolation tracks
process.analyzeHiMassTau.DoRecoElectronDiscrByEcalIsolation = cms.bool(True) 				# if true, electrons will be required to pass ecal isolation
process.analyzeHiMassTau.RecoElectronEcalIsoSumPtMaxCutValue = cms.double(1.0) 				# sum pt of ecal rechits < X
process.analyzeHiMassTau.RecoElectronEcalIsoSumPtMinCutValue = cms.double(0.0) 				# sum pt of ecal rechits >= X
process.analyzeHiMassTau.RecoElectronEcalIsoDeltaRCone = cms.double(0.4) 					# isolation conesize used to calculate sum pt of ecal rechits
process.analyzeHiMassTau.RecoElectronEcalIsoRecHitThreshold = cms.double(1.0)				# min pt requirement for isolation rechits
process.analyzeHiMassTau.DoRecoElectronDiscrByIp = cms.bool(True) 						# if true, electron will be required to have |d0|<X
#    process.analyzeHiMassTau.RecoElectronIpCut = cms.double(999.) 						# |d0(impact parameter)| < X
process.analyzeHiMassTau.DoRecoElectronDiscrByEoverP = cms.bool(False) 					# if true, require the electron to have X < E/p < Y
#    process.analyzeHiMassTau.RecoElectronEoverPMax = cms.double(1.3) 						# E/p < Y
#    process.analyzeHiMassTau.RecoElectronEoverPMin = cms.double(0.7) 						# E/p > X
process.analyzeHiMassTau.DoRecoElectronDiscrByHoverEm = cms.bool(False) 					# if true, require the electron to have H/Em < X
#    process.analyzeHiMassTau.RecoElectronHoverEmCut = cms.double(0.05)	 					# H/Em < X
process.analyzeHiMassTau.DoRecoElectronDiscrBySigmaIEtaIEta = cms.bool(False)				# if true, apply sigmaIEtaIEta cut
#    process.analyzeHiMassTau.RecoElectronSigmaIEtaIEta = cms.double(0.03)					# sigmaIEtaIEta <= X (will be used
											# ONLY if heep info is NOT used
process.analyzeHiMassTau.DoRecoElectronDiscrByDEtaIn = cms.bool(False)					# if true, apply 
#    process.analyzeHiMassTau.RecoElectronEEDEtaIn = cms.double(0.007)
#    process.analyzeHiMassTau.RecoElectronEBDEtaIn = cms.double(0.005)
process.analyzeHiMassTau.DoRecoElectronDiscrByDPhiIn = cms.bool(False)
#    process.analyzeHiMassTau.RecoElectronEEDPhiIn = cms.double(0.09)
#    process.analyzeHiMassTau.RecoElectronEBDPhiIn = cms.double(0.09)
process.analyzeHiMassTau.DoRecoElectronDiscrBySCE2by5Over5by5 = cms.bool(True)			
#    process.analyzeHiMassTau.RecoElectronEBscE1by5Over5by5 = cms.double(0.83)					
#    process.analyzeHiMassTau.RecoElectronEBscE2by5Over5by5 = cms.double(0.94)					
process.analyzeHiMassTau.DoRecoElectronDiscrByMissingHits = cms.bool(True)
#    process.analyzeHiMassTau.RecoElectronMissingHits = cms.int32(1)
#    process.analyzeHiMassTau.DoRecoElectronDiscrByEcalDrivenSeed = cms.bool(False)				# if true, require electron to be ecal driven
#    process.analyzeHiMassTau.DoRecoElectronDiscrByTrackerDrivenSeed = cms.bool(False)				# if true, require electron to be tracker driven

    #-----Reco Jet Inputs
process.analyzeHiMassTau.RecoJetSource                       = cms.InputTag('patJetsPF')              # jet collection
process.analyzeHiMassTau.RecoJetEtaMinCut                    = cms.double(0.0)                              # require jet |eta|>=X
process.analyzeHiMassTau.RecoJetEtaMaxCut                    = cms.double(999.5)                            # require jet |eta|<=X
process.analyzeHiMassTau.RecoJetPtCut                        = cms.double(30.0)                             # require jet pt>=X
process.analyzeHiMassTau.UseCorrectedJet                     = cms.bool(True)                              # if true, jet corrections are used
process.analyzeHiMassTau.RemoveJetOverlapWithMuons           = cms.bool(False)                              # if true, jets w/ dR(muon,jet)<X will not be considered
process.analyzeHiMassTau.JetMuonMatchingDeltaR               = cms.double(0.5)                              # dR(muon,jet)<X used for removing jets from the "good jet" list
process.analyzeHiMassTau.RemoveJetOverlapWithElectrons       = cms.bool(False)                              # if true, jets w/ dR(electron,jet)<X will not be considered
process.analyzeHiMassTau.JetElectronMatchingDeltaR           = cms.double(0.5)                              # dR(electron,jet)<X used for removing jets from the "good jet" list
process.analyzeHiMassTau.RemoveJetOverlapWithTaus            = cms.bool(False)                               # if true, jets w/ dR(tau,jet)<X will not be considered
process.analyzeHiMassTau.JetTauMatchingDeltaR                = cms.double(0.25)                             # dR(tau,jet)<X used for removing jets from the "good jet" list
process.analyzeHiMassTau.DoDiscrByFirstLeadingJet		= cms.bool(False)
process.analyzeHiMassTau.RecoFirstLeadingJetPt		= cms.double(100.0)
process.analyzeHiMassTau.RecoFirstLeadingJetEtaMinCut	= cms.double(0.0)
process.analyzeHiMassTau.RecoFirstLeadingJetEtaMaxCut	= cms.double(3.0)
process.analyzeHiMassTau.DoDiscrBySecondLeadingJet		= cms.bool(False)
process.analyzeHiMassTau.RecoSecondLeadingJetPt		= cms.double(100.0)
process.analyzeHiMassTau.RecoSecondLeadingJetEtaMinCut	= cms.double(0.0)
process.analyzeHiMassTau.RecoSecondLeadingJetEtaMaxCut	= cms.double(3.0)
process.analyzeHiMassTau.RemoveFirstLeadingJetOverlapWithMuons           = cms.bool(False)
process.analyzeHiMassTau.FirstLeadingJetMuonMatchingDeltaR               = cms.double(0.3)
process.analyzeHiMassTau.RemoveFirstLeadingJetOverlapWithElectrons       = cms.bool(True)
process.analyzeHiMassTau.FirstLeadingJetElectronMatchingDeltaR           = cms.double(0.3)
process.analyzeHiMassTau.RemoveFirstLeadingJetOverlapWithTaus            = cms.bool(True)
process.analyzeHiMassTau.FirstLeadingJetTauMatchingDeltaR                = cms.double(0.3)
process.analyzeHiMassTau.RemoveSecondLeadingJetOverlapWithMuons           = cms.bool(False)
process.analyzeHiMassTau.SecondLeadingJetMuonMatchingDeltaR               = cms.double(0.3)
process.analyzeHiMassTau.RemoveSecondLeadingJetOverlapWithElectrons       = cms.bool(True)
process.analyzeHiMassTau.SecondLeadingJetElectronMatchingDeltaR           = cms.double(0.3)
process.analyzeHiMassTau.RemoveSecondLeadingJetOverlapWithTaus            = cms.bool(True)
process.analyzeHiMassTau.SecondLeadingJetTauMatchingDeltaR                = cms.double(0.3)

    #-----Reco b-Jet Inputs
process.analyzeHiMassTau.RecoBJetEtaMinCut                    = cms.double(0.0)                              # require jet |eta|>=X
process.analyzeHiMassTau.RecoBJetEtaMaxCut                    = cms.double(999.4)                            # require jet |eta|<=X
process.analyzeHiMassTau.RecoBJetPtCut                        = cms.double(15.0)                             # require jet pt>=X
process.analyzeHiMassTau.RemoveBJetOverlapWithMuons           = cms.bool(False)                              # if true, jets w/ dR(muon,jet)<X will not be considered
process.analyzeHiMassTau.BJetMuonMatchingDeltaR               = cms.double(0.3)                              # dR(muon,jet)<X used for removing jets from the "good jet" list
process.analyzeHiMassTau.RemoveBJetOverlapWithElectrons       = cms.bool(True)                              # if true, jets w/ dR(electron,jet)<X will not be considered
process.analyzeHiMassTau.BJetElectronMatchingDeltaR           = cms.double(0.3)                              # dR(electron,jet)<X used for removing jets from the "good jet" list
process.analyzeHiMassTau.RemoveBJetOverlapWithTaus            = cms.bool(True)                               # if true, jets w/ dR(tau,jet)<X will not be considered
process.analyzeHiMassTau.BJetTauMatchingDeltaR                = cms.double(0.3)                             # dR(tau,jet)<X used for removing jets from the "good jet" list
process.analyzeHiMassTau.ApplyJetBTagging                    = cms.bool(True)                              # if true, apply track counting high$
process.analyzeHiMassTau.JetBTaggingTCHEcut                  = cms.double(1.7)                              # tagged as b-jet if TCHE > X

    #-----Vertex Inputs
process.analyzeHiMassTau.RecoVertexSource = cms.InputTag('offlinePrimaryVertices') 				# vertex collection
process.analyzeHiMassTau.RecoVertexMaxZposition = cms.double(20.0)						# vertex |z| < X
process.analyzeHiMassTau.RecoVertexMinTracks = cms.int32(2)						# vertex must have >= 2 "good" tracks used to reconstruct it
process.analyzeHiMassTau.RecoVertexTrackWeight = cms.double(0.5)						# weight used to define "good" tracks used to reconstruct vertex

    #-----Trigger Inputs
process.analyzeHiMassTau.RecoTriggerSource = cms.InputTag("TriggerResults","","HLT")			# trigger collection
process.analyzeHiMassTau.TriggerRequirements = cms.vstring('HLT_PFMHT150','HLT_PFMHT150_v1','HLT_PFMHT150_v2','HLT_PFMHT150_v3','HLT_PFMHT150_v4')  # trigger path name

    #-----Susy Topology Inputs
process.analyzeHiMassTau.DoSUSYDiscrByMHT = cms.bool(False)
process.analyzeHiMassTau.MhtCut = cms.double(200.0)
process.analyzeHiMassTau.DoSUSYDiscrByR1 = cms.bool(False)
process.analyzeHiMassTau.R1MinCut = cms.double(0.85)
process.analyzeHiMassTau.R1MaxCut = cms.double(999.0)
process.analyzeHiMassTau.DoSUSYDiscrByR2 = cms.bool(False)
process.analyzeHiMassTau.R2MinCut = cms.double(0.0)
process.analyzeHiMassTau.R2MaxCut = cms.double(3.6)
process.analyzeHiMassTau.DoSUSYDiscrByAlpha = cms.bool(False)
process.analyzeHiMassTau.AlphaMinCut = cms.double(0.5)
process.analyzeHiMassTau.AlphaMaxCut = cms.double(9999999999.9)
process.analyzeHiMassTau.DoSUSYDiscrByDphi1 = cms.bool(False)
process.analyzeHiMassTau.Dphi1MinCut = cms.double(0.9)
process.analyzeHiMassTau.Dphi1MaxCut = cms.double(999.9)
process.analyzeHiMassTau.DoSUSYDiscrByDphi2 = cms.bool(False)
process.analyzeHiMassTau.Dphi2MinCut = cms.double(-0.5)
process.analyzeHiMassTau.Dphi2MaxCut = cms.double(0.5)

    #-----Topology Inputs
process.analyzeHiMassTau.RecoMetSource = cms.InputTag('patMETsPFL1L2L3Cor')				         	# met collection
											# particle flow met for 2X	= layer1PFMETs
											# particle flow met for 3X	= layer1METsPF
											# standard calo met for 2X & 3X	= layer1METs
process.analyzeHiMassTau.DoDiscrByMet = cms.bool(False) 							# if true, met will be required to be > X
process.analyzeHiMassTau.CalculateMetUsingOnlyLeg1AndLeg2 = cms.bool(False)					# if true, recalculate met using leg1 and leg2 momenta
process.analyzeHiMassTau.RecoMetCut = cms.double(30.0) 							# met > X
process.analyzeHiMassTau.DoDiTauDiscrByDeltaR = cms.bool(True) 						# if true, ditau pairs must have dR(leg1,leg2) > X
process.analyzeHiMassTau.DiTauDeltaRCut = cms.double(0.3)	 						# dR(leg1,leg2) > X
process.analyzeHiMassTau.DiTauDiscrByOSLSType = cms.string('NONE')		 				# if 'OS', product of leg1 charge and leg2 charge < 0
											# if 'LS', product of leg1 charge and leg2 charge > 0
											# if anything else is used, no OSLS requirement is applied
process.analyzeHiMassTau.UseTauSeedTrackForDiTauDiscrByOSLS = cms.bool(True) 				# if true, then use q(tau) = q(tau seed track)
                                                          				# if false, then use q(tau) = sum of q(tau jet tracks)
process.analyzeHiMassTau.DoDiTauDiscrByCosDphi = cms.bool(False) 						# if true, require X <= cosine dphi(leg1,leg2) <= Y
process.analyzeHiMassTau.DiTauCosDphiMaxCut = cms.double(-0.95) 						# cosine dphi(leg1, leg2) <= Y
process.analyzeHiMassTau.DiTauCosDphiMinCut = cms.double(-1.00) 						# cosine dphi(leg1, leg2) >= X
process.analyzeHiMassTau.DoDiscrByMassReco = cms.bool(False) 						# if true, apply a mass cut on leg1+leg2+met combinations
process.analyzeHiMassTau.UseVectorSumOfVisProductsAndMetMassReco = cms.bool(False) 				# if true, calculate invariant mass using
                                                               				# vector sum of leg1, leg2, and met
process.analyzeHiMassTau.UseCollinerApproxMassReco = cms.bool(False) 					# if true, use the collinear approximation for the invariant mass
                                                 					# if both are set to false, then use only the visible products
											# to calculate the mass (leg1 + leg2)
process.analyzeHiMassTau.MassMinCut = cms.double(150.0) 							# mass > X
process.analyzeHiMassTau.MassMaxCut = cms.double(1000.0) 							# mass < Y
process.analyzeHiMassTau.DoDiTauDiscrByCDFzeta2D = cms.bool(False)						# if true, apply 2D zeta cut ( a*pzeta + b*pzetaVis > X )
process.analyzeHiMassTau.PZetaCutCoefficient = cms.double(1.0)						# a -- ( a*pzeta + b*pzetaVis > X )
process.analyzeHiMassTau.PZetaVisCutCoefficient = cms.double(-0.875)					# b -- ( a*pzeta + b*pzetaVis > X )
process.analyzeHiMassTau.CDFzeta2DCutValue = cms.double(-7.00)						# X -- ( a*pzeta + b*pzetaVis > X )
process.analyzeHiMassTau.DoDiTauDiscrByDeltaPtDivSumPt = cms.bool(False)					# if true, apply cut on deltapt(leg2,leg1) / sumpt(leg1,leg2)
process.analyzeHiMassTau.DiTauDeltaPtDivSumPtMinCutValue = cms.double(0.1)					# ( pt_leg2 - pt_leg1 ) / ( pt_leg2 + pt_leg1 ) >= X
process.analyzeHiMassTau.DiTauDeltaPtDivSumPtMaxCutValue = cms.double(1.0)					# ( pt_leg2 - pt_leg1 ) / ( pt_leg2 + pt_leg1 ) <= X
process.analyzeHiMassTau.DoDiTauDiscrByDeltaPt = cms.bool(False)						# if true, apply cut on deltapt(leg2,leg1)
process.analyzeHiMassTau.DiTauDeltaPtMinCutValue = cms.double(30.0)						# ( pt_leg2 - pt_leg1 ) >= X
process.analyzeHiMassTau.DiTauDeltaPtMaxCutValue = cms.double(9999.0)					# ( pt_leg2 - pt_leg1 ) <= X
process.analyzeHiMassTau.DoDiscrByLeg1MetDphi = cms.bool(False)						# if true, require X <= cosine dphi(leg1,met) <= Y
process.analyzeHiMassTau.Leg1MetDphiMinCut = cms.double(1.30)						# cosine dphi(leg1, met) >= X
process.analyzeHiMassTau.Leg1MetDphiMaxCut = cms.double(3.15)						# cosine dphi(leg1, met) <= Y
process.analyzeHiMassTau.DoDiscrByLeg2MetDphi = cms.bool(False)						# if true, require X <= cosine dphi(leg2,met) <= Y
process.analyzeHiMassTau.Leg2MetDphiMinCut = cms.double(1.30)						# cosine dphi(leg2, met) >= X
process.analyzeHiMassTau.Leg2MetDphiMaxCut = cms.double(3.15)						# cosine dphi(leg2, met) <= Y
process.analyzeHiMassTau.DoTauDiscrByIsZeeCut = cms.bool(False)

    #-----do matching to gen?
process.analyzeHiMassTau.MatchLeptonToGen = cms.bool(False)                                                 # if true, match reco lepton to a gen lepton
process.analyzeHiMassTau.UseLeptonMotherId = cms.bool(False)                                                # if true, require the matched lepton to come from a certain
                                                                                        # 'mother' particle
process.analyzeHiMassTau.UseLeptonGrandMotherId = cms.bool(False)                                           # if true, require the matched lepton to come from a certain
                                                                                        # 'grandmother' particle
process.analyzeHiMassTau.LeptonMotherId = cms.int32(24)                                                     # pdgId of the 'mother' particle
process.analyzeHiMassTau.LeptonGrandMotherId = cms.int32(32)                                                # pdgId of the 'grandmother' particle
process.analyzeHiMassTau.MatchTauToGen = cms.bool(False)                                                    # if true, match reco tau to a gen had tau
process.analyzeHiMassTau.UseTauMotherId = cms.bool(False)                                                   # if true, require the matched tau to come from a certain
                                                                                        # 'mother' particle ('mother' here is NOT 15!!!!! Matching
                                                                                        # for the had tau leg already requires the vis had tau to come
                                                                                        # from a tau lepton)
process.analyzeHiMassTau.UseTauGrandMotherId = cms.bool(False)                                              # if true, require the matched tau to come from a certain
                                                                                        # 'grandmother' particle
process.analyzeHiMassTau.TauMotherId = cms.int32(32)                                                        # pdgId of the 'mother' particle
process.analyzeHiMassTau.TauGrandMotherId = cms.int32(1)                                                    # pdgId of the 'grandmother' particle
process.analyzeHiMassTau.TauToGenMatchingDeltaR = cms.double(0.25)                                          # matching dR:  dR(vis gen tau,reco tau)<X

    #-----Event Sequence inputs
process.analyzeHiMassTau.RecoTriggersNmin = cms.int32(0)							# require event to pass >=X trigger paths defined above
process.analyzeHiMassTau.RecoVertexNmin = cms.int32(1)							# require event to have >=X vertices passing specified cuts
process.analyzeHiMassTau.RecoVertexNmax = cms.int32(1000)							# require event to have <=X vertices passing specified cuts
process.analyzeHiMassTau.RecoLeg1Nmin = cms.int32(1)							# require event to have >=X leg1 objects passing specified cuts
process.analyzeHiMassTau.RecoLeg1Nmax = cms.int32(1000)							# require event to have <=X leg1 objects passing specified cuts
process.analyzeHiMassTau.RecoLeg2Nmin = cms.int32(1)							# require event to have >=X leg2 objects passing specified cuts
process.analyzeHiMassTau.RecoLeg2Nmax = cms.int32(1000)							# require event to have <=X leg2 objects passing specified cuts
process.analyzeHiMassTau.RecoJetNmin = cms.int32(0)								# require event to have >=X "jets" passing specified cuts
process.analyzeHiMassTau.RecoJetNmax = cms.int32(1000)						        # require event to have <=X "jets" passing specified cuts
process.analyzeHiMassTau.RecoFirstLeadingJetNmin = cms.int32(0)
process.analyzeHiMassTau.RecoSecondLeadingJetNmin = cms.int32(0)
process.analyzeHiMassTau.RecoBJetNmin = cms.int32(1)							# require event to have >=X "jets" passing specified cuts
process.analyzeHiMassTau.RecoBJetNmax = cms.int32(1000)						        # require event to have <=X "jets" passing specified cuts
process.analyzeHiMassTau.SusyCombinationsNmin = cms.int32(0)
process.analyzeHiMassTau.CombinationsNmin = cms.int32(1)							# require event to have >=X leg1+leg2+met combinations 
											# passing specified cuts
process.analyzeHiMassTau.CombinationsNmax = cms.int32(10000)							# require event to have <=X leg1+leg2+met combinations
											# passing specified cuts

process.analyzeHiMassTau.EventSelectionSequence = cms.vstring('RecoTriggersNmin','RecoVertexNmin','RecoVertexNmax','RecoLeg1Nmin','RecoLeg1Nmax','RecoLeg2Nmin',
                                         'RecoLeg2Nmax','RecoJetNmin','RecoJetNmax','RecoFirstLeadingJetNmin','RecoSecondLeadingJetNmin',
                                         'RecoBJetNmin','RecoBJetNmax','SusyCombinationsNmin','CombinationsNmin','CombinationsNmax')

    #-----Inputs for systematic uncertainties
process.analyzeHiMassTau.CalculatePdfSystematicUncertanties = cms.bool(False)				# if true, pdf systematic uncertanties will be calculated
process.analyzeHiMassTau.PdfWeightTags = cms.untracked.VInputTag("cteq65PdfWeights")			# collection of weights for systematic uncertanties
process.analyzeHiMassTau.CalculateFSRSystematics = cms.bool(False)
process.analyzeHiMassTau.CalculateISRGluonSystematics = cms.bool(False)
process.analyzeHiMassTau.CalculateISRGammaSystematics = cms.bool(False)
process.analyzeHiMassTau.SmearTheMuon = cms.bool(False)
process.analyzeHiMassTau.MuonPtScaleOffset = cms.double(1.0)
process.analyzeHiMassTau.MuonPtSigmaOffset = cms.double(1.0)
process.analyzeHiMassTau.MuonEtaScaleOffset = cms.double(1.0)
process.analyzeHiMassTau.MuonEtaSigmaOffset = cms.double(1.0)
process.analyzeHiMassTau.MuonPhiScaleOffset = cms.double(1.0)
process.analyzeHiMassTau.MuonPhiSigmaOffset = cms.double(1.0)
process.analyzeHiMassTau.SmearTheElectron = cms.bool(False)
process.analyzeHiMassTau.ElectronPtScaleOffset = cms.double(1.0)
process.analyzeHiMassTau.ElectronPtSigmaOffset = cms.double(1.0)
process.analyzeHiMassTau.ElectronEtaScaleOffset = cms.double(1.0)
process.analyzeHiMassTau.ElectronEtaSigmaOffset = cms.double(1.0)
process.analyzeHiMassTau.ElectronPhiScaleOffset = cms.double(1.0)
process.analyzeHiMassTau.ElectronPhiSigmaOffset = cms.double(1.0)
process.analyzeHiMassTau.SmearTheTau = cms.bool(False)
process.analyzeHiMassTau.TauPtScaleOffset = cms.double(1.0)
process.analyzeHiMassTau.TauPtSigmaOffset = cms.double(1.0)
process.analyzeHiMassTau.TauEtaScaleOffset = cms.double(1.0)
process.analyzeHiMassTau.TauEtaSigmaOffset = cms.double(1.0)
process.analyzeHiMassTau.TauPhiScaleOffset = cms.double(1.0)
process.analyzeHiMassTau.TauPhiSigmaOffset = cms.double(1.0)
process.analyzeHiMassTau.SmearTheJet = cms.bool(False)
process.analyzeHiMassTau.JetEnergyScaleOffset = cms.double(1.0)    
process.analyzeHiMassTau.SmearThePt = cms.bool(False)
process.analyzeHiMassTau.SmearTheEta = cms.bool(False)
process.analyzeHiMassTau.SmearThePhi = cms.bool(False)
process.analyzeHiMassTau.CalculatePUSystematics = cms.bool(True)
process.analyzeHiMassTau.DataHistos = cms.string("AllData_NoCuts.root")
process.analyzeHiMassTau.MCHistos = cms.string("TTJets_NoCuts.root")
#    process.analyzeHiMassTau.BosonPtBinEdges = cms.untracked.vdouble(
#               0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.
#            , 10., 11., 12., 13., 14., 15., 16., 17., 18., 19.
#            , 20., 21., 22., 23., 24., 25., 26., 27., 28., 29.
#            , 30., 31., 32., 33., 34., 35., 36., 37., 38., 39.
#            , 40., 41., 42., 43., 44., 45., 46., 47., 48., 49.
#            , 999999.
#    )
#    process.analyzeHiMassTau.PtWeights = cms.untracked.vdouble( 
#              0.800665, 0.822121, 0.851249, 0.868285, 0.878733
#            , 0.953853, 0.928108, 0.982021, 1.00659 , 1.00648
#            , 1.03218 , 1.04924 , 1.03621 , 1.08743 , 1.01951
#            , 1.10519 , 0.984263, 1.04853 , 1.06724 , 1.10183
#            , 1.0503  , 1.13162 , 1.03837 , 1.12936 , 0.999173
#            , 1.01453 , 1.11435 , 1.10545 , 1.07199 , 1.04542
#            , 1.00828 , 1.0822  , 1.09667 , 1.16144 , 1.13906
#            , 1.27974 , 1.14936 , 1.23235 , 1.06667 , 1.06363
#            , 1.14225 , 1.22955 , 1.12674 , 1.03944 , 1.04639
#            , 1.13667 , 1.20493 , 1.09349 , 1.2107  , 1.21073
#    )
process.analyzeHiMassTau.ApplyMuonTriggerScaleFactors = cms.bool(False)
process.analyzeHiMassTau.ApplyElectronTriggerScaleFactors = cms.bool(False)
process.analyzeHiMassTau.ApplyTauTriggerScaleFactors = cms.bool(False)
#    process.analyzeHiMassTau.MuonTrigPtBinEdges = cms.untracked.vdouble(0., 20., 999999.)
#    process.analyzeHiMassTau.MuonTrigPtWeights = cms.untracked.vdouble(1., 1.)
#    process.analyzeHiMassTau.MuonTrigEtaBinEdges = cms.untracked.vdouble(-999999., -2.1, -1.2, -0.8, 0.8, 1.2, 2.1, 999999.)
#    process.analyzeHiMassTau.MuonTrigEtaWeights = cms.untracked.vdouble(1., 0.9661, 0.9498, 0.9837, 0.9498, 0.9661, 1.)
#    process.analyzeHiMassTau.ElectronTrigPtBinEdges = cms.untracked.vdouble(0., 20., 999999.)
#    process.analyzeHiMassTau.ElectronTrigPtWeights = cms.untracked.vdouble(1., 1.)
#    process.analyzeHiMassTau.ElectronTrigEtaBinEdges = cms.untracked.vdouble(-999999., -2.5, -1.479, 1.479, 2.5, 999999.)
#    process.analyzeHiMassTau.ElectronTrigEtaWeights = cms.untracked.vdouble(1., 0.9713, 0.9757, 0.9713, 1.)
#    process.analyzeHiMassTau.TauTrigPtBinEdges = cms.untracked.vdouble(0., 20., 999999.)
#    process.analyzeHiMassTau.TauTrigPtWeights = cms.untracked.vdouble(1., 1.)
#    process.analyzeHiMassTau.TauTrigEtaBinEdges = cms.untracked.vdouble(-999999., -2.1, -1.2, -0.8, 0.8, 1.2, 2.1, 999999.)
#    process.analyzeHiMassTau.TauTrigEtaWeights = cms.untracked.vdouble(1., 1., 1., 1., 1., 1., 1.)

process.p = cms.Path(
  process.analyzeHiMassTau
)

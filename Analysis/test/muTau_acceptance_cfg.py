import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('HiMassTau')

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("HiMassTauAnalyzer.Analysis.FILESTOREAD")

process.eTauCandidates = cms.EDProducer("DeltaRMinCandCombiner",
  decay = cms.string('selectedLayer1FixedConePFTaus@+ selectedLayer1Electrons@-'),
  checkCharge = cms.bool(False),
  cut = cms.string(''),
  name = cms.string('etauCandidates'),
  deltaRMin = cms.double(0.7),
  roles = cms.vstring('tau', 'electron')
)    

process.analyzeHiMassTau = cms.EDAnalyzer('HiMassTauAnalysis',

    #-----Generator level Inputs
    GenParticleSource = cms.untracked.InputTag('genParticles'),				# gen particle collection

    #-----Inputs to determine which channel to analyze
    AnalyzeTauForLeg1		= cms.bool(False),					# if true, taus will be used for leg1
    AnalyzeMuonForLeg1		= cms.bool(True),					# if true, muons will be used for leg1
    AnalyzeElectronForLeg1	= cms.bool(False),					# if true, electrons will be used for leg1
    AnalyzeTauForLeg2		= cms.bool(True),					# if true, taus will be used for leg2
    AnalyzeMuonForLeg2		= cms.bool(False),					# if true, muons will be used for leg2
    AnalyzeElectronForLeg2	= cms.bool(False),					# if true, electrons will be used for leg2

    #-----Reco Tau Inputs
    RecoTauSource = cms.InputTag('selectedLayer1FixedConePFTaus'),			# other choices include:
											# selectedLayer1FixedConeHighEffPFTaus
											# selectedLayer1ShrinkingConeHighEffPFTaus
											# selectedLayer1ShrinkingConePFTaus
    RecoTauEtaCut = cms.double(2.1),							# require tau |eta|<=X
    RecoTauPtMinCut = cms.double(15.),							# require tau pt>=X
    RecoTauPtMaxCut = cms.double(9999.),						# require tau pt<=X
    DoRecoTauDiscrByLeadTrack = cms.bool(True),						# if true, tau is required to pass a lead track pt cut
    UseRecoTauDiscrByLeadTrackFlag = cms.bool(True), 					# if true, default seed track discriminator is used
                                                     					# if false, seed track cut will be recalculated using the parameters below
    RecoTauDiscrByLeadTrack = cms.untracked.string('leadingTrackPtCut'),		# name of the lead track discriminator flag
    DoRecoTauDiscrByLeadTrackNhits = cms.bool(False),					# if true, tau leading track is required to have >= X hits
    RecoTauLeadTrackMinHits = cms.int32(12),						# tau leading track hits >= X
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
    RecoTauIsoSumPtMaxCutValue = cms.double(1.0),					# sum pt of tracks & gammas < X
    RecoTauIsoSumPtMinCutValue = cms.double(0.0),					# sum pt of tracks & gammas >= X
    RecoTauDiscrByProngType = cms.string('NONE'),		 			# if '1or3', taus will be required to have 1 || 3 prongs
											# if '1', taus will be required to have 1 prongs
											# if '3', taus will be required to have 3 prongs
											# if anything else is used, no prong cuts are applied
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

    #-----Reco Muon Inputs
    RecoMuonSource = cms.InputTag('selectedLayer1Muons'),				# muon collection
    RecoMuonEtaCut = cms.double(2.1),							# require muon |eta|<=X
    RecoMuonPtMinCut = cms.double(15.),							# require muon pt>=X
    RecoMuonPtMaxCut = cms.double(9999.),						# require muon pt<=X
    DoRecoMuonDiscrByGlobal = cms.bool(True),						# if true, muon will be required to be a 'global muon
    DoRecoMuonDiscrByIsolation = cms.bool(False),					# if true, muon isolation will be applied
    RecoMuonIsoSumPtMaxCutValue = cms.double(3.5),					# sum pt of isolation tracks & ecal rechits < X
    RecoMuonIsoSumPtMinCutValue = cms.double(0.0),					# sum pt of isolation tracks & ecal rechits >= X
    RecoMuonIsoDeltaRCone = cms.double(0.6),						# outer conesize used for isolation
    RecoMuonTrackIsoTrkThreshold = cms.double(1.0),					# isolation tracks are required to have pt>X
    RecoMuonEcalIsoRecHitThreshold = cms.double(1.0),					# isolation rechits are required to have pt>X
    DoRecoMuonDiscrByIp = cms.bool(False),						# if true, muon will be required to have |d0|<X
    RecoMuonIpCut = cms.double(0.028),							# |d0|<X

    #-----Reco Electron Inputs
    RecoElectronSource = cms.InputTag('selectedLayer1Electrons'),			# electron collection
    RecoElectronEtaCut = cms.double(2.1),						# require electron |eta|<=X
    RecoElectronPtMinCut = cms.double(15.),						# require electron pt>=X
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

    #-----Reco Jet Inputs
    RecoJetSource			= cms.InputTag('selectedLayer1Jets'),		# jet collection
    RecoJetEtaMinCut			= cms.double(0.0),				# require jet |eta|>=X
    RecoJetEtaMaxCut			= cms.double(999.5),				# require jet |eta|<=X
    RecoJetPtCut			= cms.double(15.0),				# require jet pt>=X
    UseCorrectedJet			= cms.bool(False),				# if true, jet corrections are used
    RemoveJetOverlapWithMuons		= cms.bool(True),				# if true, jets w/ dR(muon,jet)<X will not be considered
    JetMuonMatchingDeltaR		= cms.double(0.5),				# dR(muon,jet)<X used for removing jets from the "good jet" list
    RemoveJetOverlapWithElectrons	= cms.bool(False),				# if true, jets w/ dR(electron,jet)<X will not be considered
    JetElectronMatchingDeltaR		= cms.double(0.5),				# dR(electron,jet)<X used for removing jets from the "good jet" list
    RemoveJetOverlapWithTaus		= cms.bool(True),				# if true, jets w/ dR(tau,jet)<X will not be considered
    JetTauMatchingDeltaR		= cms.double(0.25),				# dR(tau,jet)<X used for removing jets from the "good jet" list

    #-----Vertex Inputs
    RecoVertexSource = cms.InputTag('offlinePrimaryVertices'), 				# vertex collection
    RecoVertexMaxZposition = cms.double(999.0),						# vertex |z| < X
    RecoVertexMinTracks = cms.int32(-1),						# vertex must have >= 2 "good" tracks used to reconstruct it
    RecoVertexTrackWeight = cms.double(0.5),						# weight used to define "good" tracks used to reconstruct vertex

    #-----Trigger Inputs
    RecoTriggerSource = cms.InputTag("TriggerResults","","HLT"),			# trigger collection
    TriggerRequirements = cms.vstring('HLT_Mu9'),					# trigger path name

    #-----Topology Inputs
    RecoDiTauSource = cms.InputTag('eTauCandidates'),					# ditau collection needed for the filling of the Ntuple
    RecoMetSource = cms.InputTag('layer1METsPF'),					# met collection
											# particle flow met for 2X	= layer1PFMETs
											# particle flow met for 3X	= layer1METsPF
											# standard calo met for 2X & 3X	= layer1METs
    DoDiscrByMet = cms.bool(False), 							# if true, met will be required to be > X
    RecoMetCut = cms.double(30.0), 							# met > X
    DoDiTauDiscrByDeltaR = cms.bool(True), 						# if true, ditau pairs must have dR(leg1,leg2) > X
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
    UseVectorSumOfVisProductsAndMetMassReco = cms.bool(True), 				# if true, calculate invariant mass using
                                                               				# vector sum of leg1, leg2, and met
    UseCollinerApproxMassReco = cms.bool(False), 					# if true, use the collinear approximation for the invariant mass
                                                 					# if both are set to false, then use only the visible products
											# to calculate the mass (leg1 + leg2)
    MassMinCut = cms.double(260.0), 							# mass > X
    MassMaxCut = cms.double(1000.0), 							# mass < Y
    DoDiTauDiscrByCDFzeta2D = cms.bool(False),						# if true, apply 2D zeta cut ( a*pzeta + b*pzetaVis > X )
    PZetaCutCoefficient = cms.double(1.0),						# a -- ( a*pzeta + b*pzetaVis > X )
    PZetaVisCutCoefficient = cms.double(-0.875),					# b -- ( a*pzeta + b*pzetaVis > X )
    CDFzeta2DCutValue = cms.double(-7.00),						# X -- ( a*pzeta + b*pzetaVis > X )
    DoDiTauDiscrByDeltaPtDivSumPt = cms.bool(False),					# if true, apply cut on deltapt(leg2,leg1) / sumpt(leg1,leg2)
    DiTauDeltaPtDivSumPtMinCutValue = cms.double(0.0),					# ( pt_leg2 - pt_leg1 ) / ( pt_leg2 + pt_leg1 ) >= X
    DiTauDeltaPtDivSumPtMaxCutValue = cms.double(1.0),					# ( pt_leg2 - pt_leg1 ) / ( pt_leg2 + pt_leg1 ) <= X
    DoDiscrByLeg1MetDphi = cms.bool(False),						# if true, require X <= cosine dphi(leg1,met) <= Y
    Leg1MetDphiMinCut = cms.double(1.30),						# cosine dphi(leg1, met) >= X
    Leg1MetDphiMaxCut = cms.double(3.15),						# cosine dphi(leg1, met) <= Y
    DoDiscrByLeg2MetDphi = cms.bool(False),						# if true, require X <= cosine dphi(leg2,met) <= Y
    Leg2MetDphiMinCut = cms.double(1.30),						# cosine dphi(leg2, met) >= X
    Leg2MetDphiMaxCut = cms.double(3.15),						# cosine dphi(leg2, met) <= Y

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
    TauToGenMatchingDeltaR = cms.double(0.2),						# matching dR:  dR(vis gen tau,reco tau)<X

    #-----ntuple Inputs
    DoProduceNtuple = cms.bool(False),							# if true, Ntuple will be filled
    NtupleTreeName = cms.untracked.string('HMTTree'),					# name of the Ntuple tree

    #-----Fill Histograms? Histograms are filled for events passing the specified cuts
    FillRecoVertexHists = cms.bool(True),						# if true, fill histograms for vertices
    FillGenTauHists = cms.bool(False),							# if true, fill histograms for gen had taus
    FillRecoTauHists = cms.bool(True),							# if true, fill histograms for reco taus
    FillRecoMuonHists = cms.bool(True),							# if true, fill histograms for reco muons
    FillRecoElectronHists = cms.bool(False),						# if true, fill histograms for reco electrons
    FillRecoJetHists = cms.bool(True),							# if true, fill histograms for reco jets
    FillTopologyHists = cms.bool(True),							# if true, fill topology histograms (e.g. met, mass, ...)

    #-----Event Sequence inputs
    RecoTriggersNmin = cms.int32(0),							# require event to pass >=X trigger paths defined above
    RecoVertexNmin = cms.int32(0),							# require event to have >=X vertices passing specified cuts
    RecoVertexNmax = cms.int32(1000),							# require event to have <=X vertices passing specified cuts
    RecoLeg1Nmin = cms.int32(1),							# require event to have >=X leg1 objects passing specified cuts
    RecoLeg1Nmax = cms.int32(1000),							# require event to have <=X leg1 objects passing specified cuts
    RecoLeg2Nmin = cms.int32(1),							# require event to have >=X leg2 objects passing specified cuts
    RecoLeg2Nmax = cms.int32(1000),							# require event to have <=X leg2 objects passing specified cuts
    RecoJetNmin = cms.int32(0),								# require event to have >=X "jets" passing specified cuts
    RecoJetNmax = cms.int32(1000),							# require event to have <=X "jets" passing specified cuts
    CombinationsNmin = cms.int32(1),							# require event to have >=X leg1+leg2+met combinations 
											# passing specified cuts
    CombinationsNmax = cms.int32(1000),							# require event to have <=X leg1+leg2+met combinations
											# passing specified cuts

    EventSelectionSequence = cms.vstring('RecoTriggersNmin','RecoVertexNmin','RecoVertexNmax','RecoLeg1Nmin','RecoLeg1Nmax','RecoLeg2Nmin','RecoLeg2Nmax',
                                         'RecoJetNmin','RecoJetNmax','CombinationsNmin','CombinationsNmax'),

    #-----Inputs for systematic uncertainties
    CalculatePdfSystematicUncertanties = cms.bool(False),				# if true, pdf systematic uncertanties will be calculated
    PdfWeightTags = cms.untracked.VInputTag("cteq65PdfWeights"),			# collection of weights for systematic uncertanties
    #PdfWeightTags = cms.untracked.VInputTag("cteq65PdfWeights", "MRST2006nnloPdfWeights", "MRST2007lomodPdfWeights"),
    SmearTheMuon = cms.bool(False),
    RelativeMuonPtOffset = cms.double(0.0),
    RelativeMuonPtSigma = cms.double(0.0309692),
    AbsoluteMuonEtaOffset = cms.double(0.0),
    AbsoluteMuonEtaSigma = cms.double(0.000706872),
    AbsoluteMuonPhiOffset = cms.double(0.0),
    AbsoluteMuonPhiSigma = cms.double(0.000406238),
    SmearTheElectron = cms.bool(False),
    RelativeElectronPtOffset = cms.double(0.0),
    RelativeElectronPtSigma = cms.double(0.0),
    AbsoluteElectronEtaOffset = cms.double(0.0),
    AbsoluteElectronEtaSigma = cms.double(0.0),
    AbsoluteElectronPhiOffset = cms.double(0.0),
    AbsoluteElectronPhiSigma = cms.double(0.0),
)

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("outputFILENAME")
)

process.p = cms.Path(
    process.eTauCandidates*
    process.analyzeHiMassTau
    )

# print-out all python configuration parameter information
#print process.dumpPython()


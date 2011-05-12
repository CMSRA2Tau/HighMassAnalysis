import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/user/lpctau/HighMassTau/florez/TTJets_TuneZ2_7TeV-madgraph-tauola/Ttbar413MuTauSkimPat/9798e5ab1d4a976e6f05eba49b220175/skimPat_9_1_QV8.root'
    )
)

process.elecMatchMap = cms.EDProducer("MCTruthDeltaRMatcherNew",
    src = cms.InputTag('selectedPatElectrons'),
    matched = cms.InputTag('genParticles'),
    distMin = cms.double(0.15),
    matchPDGId = cms.vint32(11)
)

process.muonMatchMap = cms.EDProducer("MCTruthDeltaRMatcherNew",
    src = cms.InputTag('selectedPatMuons'),
    matched = cms.InputTag('genParticles'),
    distMin = cms.double(0.15),
    matchPDGId = cms.vint32(13)
)

process.smearElectron = cms.EDProducer('SmearedElectronProducer',
     ElectronSource = cms.InputTag('selectedPatElectrons'),
     GenMatchMapTag = cms.untracked.InputTag('elecMatchMap'),
     SmearThePt = cms.bool(True),
     PtScaleOffset = cms.double(1.05),
     PtSigmaOffset = cms.double(1.05),
     SmearTheEta = cms.bool(False),
     EtaScaleOffset = cms.double(1.0),
     EtaSigmaOffset = cms.double(1.0),
     SmearThePhi = cms.bool(False),
     PhiScaleOffset = cms.double(1.0),
     PhiSigmaOffset = cms.double(1.0)
)

process.smearMuon = cms.EDProducer('SmearedMuonProducer',
     MuonSource = cms.InputTag('selectedPatMuons'),
     GenMatchMapTag = cms.untracked.InputTag('muonMatchMap'),
     SmearThePt = cms.bool(True),
     PtScaleOffset = cms.double(1.05),
     PtSigmaOffset = cms.double(1.05),
     SmearTheEta = cms.bool(False),
     EtaScaleOffset = cms.double(1.0),
     EtaSigmaOffset = cms.double(1.0),
     SmearThePhi = cms.bool(False),
     PhiScaleOffset = cms.double(1.0),
     PhiSigmaOffset = cms.double(1.0)
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

  
process.p = cms.Path(
     (process.elecMatchMap * process.smearElectron) + 
     (process.muonMatchMap*process.smearMuon)
)

process.e = cms.EndPath(process.out)

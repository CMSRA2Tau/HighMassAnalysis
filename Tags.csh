cvs co -r V01-02-00 RecoTauTag/RecoTau 
cvs co -r V01-02-00 RecoTauTag/TauTagTools
cvs co -r V01-02-00 RecoTauTag/Configuration
cvs co -r -d SHarper/HEEPAnalyzer UserCode/SHarper/HEEPAnalyzer
cvs co -r V03-01-21  CondFormats/JetMETObjects
cvs co -r V02-02-03  JetMETCorrections/Algorithms
cvs co -r V05-00-17  JetMETCorrections/Modules
cvs co -r V03-02-07  JetMETCorrections/Configuration
cvs co -r V08-06-16 PhysicsTools/PatAlgos
cvs co -r V04-04-04 JetMETCorrections/Type1MET
cvs co -d HighMassAnalysis/Configuration -r for423_05232011 UserCode/AlfredoGurrola/HighMassAnalysis/Configuration
cvs co -d HighMassAnalysis/Skimming -r for413_05112011  UserCode/AlfredoGurrola/HighMassAnalysis/Skimming
scram build -c
scram b -j8

cvs co -r V01-02-00 RecoTauTag/RecoTau 
cvs co -r V01-02-00 RecoTauTag/TauTagTools
cvs co -r V01-02-00 RecoTauTag/Configuration
cvs co -d SHarper/HEEPAnalyzer UserCode/SHarper/HEEPAnalyzer
cvs co -r V03-01-21  CondFormats/JetMETObjects
cvs co -r V02-02-03  JetMETCorrections/Algorithms
cvs co -r V05-00-17  JetMETCorrections/Modules
cvs co -r V03-02-07  JetMETCorrections/Configuration
cvs co -r V08-06-16 PhysicsTools/PatAlgos
cvs co -r V04-04-04 JetMETCorrections/Type1MET
cvs co -d HighMassAnalysis/Configuration -r for423_06222011 UserCode/AlfredoGurrola/HighMassAnalysis/Configuration
cvs co -d HighMassAnalysis/Skimming -r for423_06222011  UserCode/AlfredoGurrola/HighMassAnalysis/Skimming
cvs co -d HighMassAnalysis/Analysis -r forSusy_423_06222011_b UserCode/AlfredoGurrola/HighMassAnalysis/Analysis
scram build -c
scram b -j8

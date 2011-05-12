#include "HighMassAnalysis/Utilities/interface/SmearedElectronProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

SmearedElectronProducer::SmearedElectronProducer(const edm::ParameterSet& iConfig){
  produces<std::vector<pat::Electron> >();
   
  _ElectronSource = iConfig.getParameter<edm::InputTag>("ElectronSource");
  _GenMatchMap = iConfig.getUntrackedParameter<edm::InputTag>("GenMatchMapTag");
  
  _SmearThePt = iConfig.getParameter<bool>("SmearThePt");
  _PtScaleOffset = iConfig.getParameter<double>("PtScaleOffset");
  _PtSigmaOffset = iConfig.getParameter<double>("PtSigmaOffset");
  
  _SmearTheEta = iConfig.getParameter<bool>("SmearTheEta");
  _EtaScaleOffset = iConfig.getParameter<double>("EtaScaleOffset");
  _EtaSigmaOffset = iConfig.getParameter<double>("EtaSigmaOffset");
  
  _SmearThePhi = iConfig.getParameter<bool>("SmearThePhi");
  _PhiScaleOffset = iConfig.getParameter<double>("PhiScaleOffset");
  _PhiSigmaOffset = iConfig.getParameter<double>("PhiSigmaOffset");
}


SmearedElectronProducer::~SmearedElectronProducer(){
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void SmearedElectronProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  													     
  if(iEvent.isRealData()) return;
 
  Handle<pat::ElectronCollection> theElectronCollection; 
  iEvent.getByLabel(_ElectronSource, theElectronCollection);

  //Handle<reco::GenParticleCollection> theGenParticleCollection;
  Handle<reco::GenParticleMatch> genMatchMap;
  if(!iEvent.getByLabel(_GenMatchMap, genMatchMap)){
    edm::LogError("") << ">>> Electron-GenParticle match map does not exist !!!"; 
    return; 
  }
  
  std::auto_ptr<pat::ElectronCollection> smearedElectrons (new pat::ElectronCollection);
												     
  for (unsigned int vIt = 0; vIt < theElectronCollection->size(); vIt++) {
    edm::Ref<std::vector<pat::Electron> > elec(theElectronCollection, vIt);
		    	 
    double smearedPt;
    double smearedEta;
    double smearedPhi;
  
    double ptgen = elec->pt();
    double etagen = elec->eta(); 
    double phigen = elec->phi();  									    	 
    int chrgen = elec->charge();  									    	 
    reco::GenParticleRef gen = (*genMatchMap)[elec];							    	 
    if( !gen.isNull()) {										    	 
      ptgen = gen->pt();										      
      etagen = gen->eta();										      
      phigen = gen->phi();										      
      chrgen = gen->charge();										      
      LogTrace("") << ">>> Electron-GenParticle match found; ptelec= " << elec->pt() << ", ptgen= " << ptgen;	      
    } 
    else{
      LogTrace("") << ">>> Electron-GENPARTICLE MATCH NOT FOUND!!!";
    }
    
    if(_SmearThePt){
      smearedPt = (ptgen * _PtScaleOffset) + ((elec->pt() -  ptgen) * _PtSigmaOffset);
    }
    else
      smearedPt = elec->pt();
    
    if(_SmearTheEta){
      smearedEta = (etagen - _EtaScaleOffset) + ((elec->eta() - etagen) * _EtaSigmaOffset);
    }
    else 
      smearedEta = elec->eta();
    
    if(_SmearThePhi){
      smearedPhi = (phigen - _PhiScaleOffset) + ((elec->phi() - phigen) * _PhiSigmaOffset);
    }
    else
      smearedPhi = elec->phi();
      
    pat::Electron *theSmearedElectron =  elec->clone();
    theSmearedElectron->setP4(reco::Particle::PolarLorentzVector(smearedPt, smearedEta, smearedPhi, elec->mass()));
    smearedElectrons->push_back(*theSmearedElectron);

  }
  iEvent.put(smearedElectrons);													     
}

// ------------ method called once each job just before starting event loop  ------------
void SmearedElectronProducer::beginJob(){
}

// ------------ method called once each job just after ending the event loop  ------------
void SmearedElectronProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void SmearedElectronProducer::beginRun(edm::Run&, edm::EventSetup const&){
}

// ------------ method called when ending the processing of a run  ------------
void SmearedElectronProducer::endRun(edm::Run&, edm::EventSetup const&){
}

// ------------ method called when starting to processes a luminosity block  ------------
void SmearedElectronProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&){
}

// ------------ method called when ending the processing of a luminosity block  ------------
void SmearedElectronProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&){
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SmearedElectronProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SmearedElectronProducer);

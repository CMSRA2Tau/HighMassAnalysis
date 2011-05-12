#include "HighMassAnalysis/Utilities/interface/SmearedMuonProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

SmearedMuonProducer::SmearedMuonProducer(const edm::ParameterSet& iConfig){
  produces<std::vector<pat::Muon> >();
   
  _MuonSource = iConfig.getParameter<edm::InputTag>("MuonSource");
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


SmearedMuonProducer::~SmearedMuonProducer(){
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void SmearedMuonProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  													     
  if(iEvent.isRealData()) return;
 
  Handle<pat::MuonCollection> theMuonCollection; 
  iEvent.getByLabel(_MuonSource, theMuonCollection);

  //Handle<reco::GenParticleCollection> theGenParticleCollection;
  Handle<reco::GenParticleMatch> genMatchMap;
  if(!iEvent.getByLabel(_GenMatchMap, genMatchMap)){
    edm::LogError("") << ">>> Muon-GenParticle match map does not exist !!!"; 
    return; 
  }
  
  std::auto_ptr<pat::MuonCollection> smearedMuons (new pat::MuonCollection);
												     
  for (unsigned int vIt = 0; vIt < theMuonCollection->size(); vIt++) {
    edm::Ref<std::vector<pat::Muon> > mu(theMuonCollection, vIt);
		    	 
    double smearedPt;
    double smearedEta;
    double smearedPhi;
  
    double ptgen = mu->pt();
    double etagen = mu->eta(); 
    double phigen = mu->phi();  									    	 
    int chrgen = mu->charge();  									    	 
    reco::GenParticleRef gen = (*genMatchMap)[mu];							    	 
    if( !gen.isNull()) {										    	 
      ptgen = gen->pt();										      
      etagen = gen->eta();										      
      phigen = gen->phi();										      
      chrgen = gen->charge();										      
      LogTrace("") << ">>> Muon-GenParticle match found; ptmu= " << mu->pt() << ", ptgen= " << ptgen;	      
    } 
    else{
      LogTrace("") << ">>> MUON-GENPARTICLE MATCH NOT FOUND!!!";
    }
    
    if(_SmearThePt){
      smearedPt = (ptgen * _PtScaleOffset) + ((mu->pt() -  ptgen) * _PtSigmaOffset);
    }
    else
      smearedPt = mu->pt();
    
    if(_SmearTheEta){
      smearedEta = (etagen - _EtaScaleOffset) + ((mu->eta() - etagen) * _EtaSigmaOffset);
    }
    else 
      smearedEta = mu->eta();
    
    if(_SmearThePhi){
      smearedPhi = (phigen - _PhiScaleOffset) + ((mu->phi() - phigen) * _PhiSigmaOffset);
    }
    else
      smearedPhi = mu->phi();
      
    pat::Muon *theSmearedMuon =  mu->clone();
    theSmearedMuon->setP4(reco::Particle::PolarLorentzVector(smearedPt, smearedEta, smearedPhi, mu->mass()));
    smearedMuons->push_back(*theSmearedMuon);

  }
  iEvent.put(smearedMuons);													     
}

// ------------ method called once each job just before starting event loop  ------------
void SmearedMuonProducer::beginJob(){
}

// ------------ method called once each job just after ending the event loop  ------------
void SmearedMuonProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void SmearedMuonProducer::beginRun(edm::Run&, edm::EventSetup const&){
}

// ------------ method called when ending the processing of a run  ------------
void SmearedMuonProducer::endRun(edm::Run&, edm::EventSetup const&){
}

// ------------ method called when starting to processes a luminosity block  ------------
void SmearedMuonProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&){
}

// ------------ method called when ending the processing of a luminosity block  ------------
void SmearedMuonProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&){
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SmearedMuonProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SmearedMuonProducer);

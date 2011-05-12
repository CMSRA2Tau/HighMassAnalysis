// -*- C++ -*-
//
// Package:    SmearedElectronProducer
// Class:      SmearedElectronProducer
// 
/**\class SmearedElectronProducer SmearedElectronProducer.cc HighMassAnalysis/SmearedElectronProducer/src/SmearedElectronProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Eduardo Luiggi
//         Created:  Fri May  6 16:01:23 CDT 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


//
// class declaration
//

class SmearedElectronProducer : public edm::EDProducer {
   public:
      explicit SmearedElectronProducer(const edm::ParameterSet&);
      ~SmearedElectronProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      
      edm::InputTag _ElectronSource;
      edm::InputTag _GenMatchMap;
      
      bool _SmearThePt;
      double _PtScaleOffset;
      double _PtSigmaOffset;
      
      bool _SmearTheEta;
      double _EtaScaleOffset;
      double _EtaSigmaOffset;
      
      bool _SmearThePhi;
      double _PhiScaleOffset;
      double _PhiSigmaOffset;

      // ----------member data ---------------------------
};


// -*- C++ -*-
//
// Package:    BetaFilter
// Class:      BetaFilter
//
/**\class BetaFilter BetaFilter.cc psi2s1s/BetaFilter/src/BetaFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
//
//

// user include files
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include <cmath>
#include <cstdlib>
#include <string>

//
// class declaration
//

class BetaFilter : public edm::global::EDFilter<> {
public:
  explicit BetaFilter(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  bool filter(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  // ----------member data ---------------------------

  const edm::EDGetToken token_;
  const double minBeta;
  const double maxBeta;
  const double maxEta;
  const int particleID;
};

BetaFilter::BetaFilter(const edm::ParameterSet& iConfig)
    : token_(consumes<edm::HepMCProduct>(
          edm::InputTag(iConfig.getUntrackedParameter("moduleLabel", std::string("generator")), "unsmeared"))),
      minBeta(iConfig.getUntrackedParameter("MinBeta", 0.)),
      maxBeta(iConfig.getUntrackedParameter("MaxBeta", 10.)),
      maxEta(iConfig.getUntrackedParameter("MaxEta", 10.)),
      particleID(iConfig.getUntrackedParameter("ParticleID", 0)) {}

// ------------ method called on each new Event  ------------
bool BetaFilter::filter(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const {
  bool accepted = false;

  edm::Handle<edm::HepMCProduct> evt;
  iEvent.getByToken(token_, evt);
  const HepMC::GenEvent* myGenEvent = evt->GetEvent();

  for (HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end();
       ++p) {
    //if ((*p)->status() != 1)
    //  continue;
    //if ((*p)->momentum().perp() > minPt && std::fabs((*p)->momentum().eta()) < maxEta &&
    //    (*p)->momentum().perp() < maxPt && std::fabs((*p)->momentum().eta()) > minEta) {
    //  if (std::abs((*p)->pdg_id()) == particleID)
    //    nLeptons++;
    //}
    //if (nLeptons >= 4) {
    //  accepted = true;
    //  break;
    //}
    if (std::abs((*p)->pdg_id()) == particleID){
       double custom_p = sqrt((*p)->momentum().px()*(*p)->momentum().px()+(*p)->momentum().py()*(*p)->momentum().py()+(*p)->momentum().pz()*(*p)->momentum().pz());
       if (custom_p/(*p)->momentum().e() > minBeta and custom_p/(*p)->momentum().e() < maxBeta and fabs((*p)->momentum().eta()) < maxEta) 
	       accepted=true;
    }
  }
  return accepted;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void BetaFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(BetaFilter);

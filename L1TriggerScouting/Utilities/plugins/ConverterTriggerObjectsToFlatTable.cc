#include "FWCore/Framework/interface/MakerMacros.h"

#include <fstream>
#include <memory>
#include <string>
#include <cmath>
#include <cstdint>

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

using namespace edm;
using namespace reco;
using namespace std;

class ConverterTriggerObjectsToFlatTable : public edm::stream::EDProducer<> {
public:
  // constructor and destructor
  explicit ConverterTriggerObjectsToFlatTable(const edm::ParameterSet&);
  ~ConverterTriggerObjectsToFlatTable() override;

  void produce(edm::Event&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  void beginRun(edm::Run const&, edm::EventSetup const&) override;
  void endRun(edm::Run const&, edm::EventSetup const&) override;

  // the tokens to access the data
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjects_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;

  std::string name_, doc_;

};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ConverterTriggerObjectsToFlatTable::ConverterTriggerObjectsToFlatTable(const edm::ParameterSet& iConfig)
    : triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone> >(iConfig.getParameter<edm::InputTag>("objects"))),
      triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<nanoaod::FlatTable>();

}
// -----------------------------------------------------------------------------


ConverterTriggerObjectsToFlatTable::~ConverterTriggerObjectsToFlatTable() {

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}


// ----------------------- method called for each orbit  -----------------------
void ConverterTriggerObjectsToFlatTable::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::vector<float> pt;
  std::vector<float> eta;
  std::vector<float> phi;

  edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
     obj.unpackPathNames(names);
     std::vector pathNamesAll = obj.pathNames(false);
     for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
	 string hltName = pathNamesAll[h];
         if(!(hltName.find("HLT_IsoMu24_v") == string::npos) or !(hltName.find("HLT_Mu50_v") == string::npos)){
            bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
	    if (isLF){
               pt.push_back(obj.pt());
               eta.push_back(obj.eta());
               phi.push_back(obj.phi());
            }
	 }
     }
  }
  
  auto out = std::make_unique<nanoaod::FlatTable>(pt.size(), name_, false, false);
  out->setDoc(doc_);

  out->addColumn<float>("pt", pt, "pt");
  out->addColumn<float>("eta", eta, "eta");
  out->addColumn<float>("phi", phi, "phi");

  iEvent.put(std::move(out));
}


// ------------ method called when starting to processes a run  ------------
void ConverterTriggerObjectsToFlatTable::beginRun(edm::Run const& run, edm::EventSetup const& iSetup) {

}

// ------------ method called when ending to processes a run  ------------
void ConverterTriggerObjectsToFlatTable::endRun(edm::Run const&, edm::EventSetup const&) {}

void ConverterTriggerObjectsToFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ConverterTriggerObjectsToFlatTable);

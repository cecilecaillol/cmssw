#include "FWCore/Framework/interface/MakerMacros.h"

#include <fstream>
#include <memory>
#include <string>
#include <cmath>
#include <cstdint>

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingBMTFStub.h"
#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

#include "L1TriggerScouting/Utilities/interface/convertToL1TFormat.h"

using namespace l1ScoutingRun3;

class ConverterEtSumToFlatTable : public edm::stream::EDProducer<> {
public:
  // constructor and destructor
  explicit ConverterEtSumToFlatTable(const edm::ParameterSet&);
  ~ConverterEtSumToFlatTable() override;

  void produce(edm::Event&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  void beginRun(edm::Run const&, edm::EventSetup const&) override;
  void endRun(edm::Run const&, edm::EventSetup const&) override;

  // the tokens to access the data
  edm::EDGetTokenT<BXVector<l1t::EtSum>> src_;

  std::string name_, doc_;

};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ConverterEtSumToFlatTable::ConverterEtSumToFlatTable(const edm::ParameterSet& iConfig)
    : src_(consumes<BXVector<l1t::EtSum>>(iConfig.getParameter<edm::InputTag>("src"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")){
  produces<nanoaod::FlatTable>();

}
// -----------------------------------------------------------------------------


ConverterEtSumToFlatTable::~ConverterEtSumToFlatTable() {

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}


// ----------------------- method called for each orbit  -----------------------
void ConverterEtSumToFlatTable::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<BXVector<l1t::EtSum>> src;
  iEvent.getByToken(src_, src);

  std::vector<float> pt;
  std::vector<float> phi;
  std::vector<int> bx;

  int startBX = src->getFirstBX();

  for (int ibx = startBX; ibx <= src->getLastBX(); ++ibx) {
     for (auto itr = src->begin(ibx); itr != src->end(ibx); ++itr) {
       l1t::EtSum::EtSumType type = itr->getType();
       if (type == l1t::EtSum::EtSumType::kMissingEt) {
	 pt.push_back(itr->pt());
         phi.push_back(itr->phi());
	 bx.push_back(ibx);
       }
     }
  }

  auto out = std::make_unique<nanoaod::FlatTable>(pt.size(), name_, false, false);
  out->setDoc(doc_);

  out->addColumn<float>("pt", pt, "pt (physical units)");
  out->addColumn<float>("phi", phi, "phi (physical units)");
  out->addColumn<int>("bx", bx, "bx");

  iEvent.put(std::move(out));
}


// ------------ method called when starting to processes a run  ------------
void ConverterEtSumToFlatTable::beginRun(edm::Run const& run, edm::EventSetup const& iSetup) {

}

// ------------ method called when ending to processes a run  ------------
void ConverterEtSumToFlatTable::endRun(edm::Run const&, edm::EventSetup const&) {}

void ConverterEtSumToFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ConverterEtSumToFlatTable);

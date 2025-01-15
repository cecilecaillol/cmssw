#include "FWCore/Framework/interface/MakerMacros.h"

#include <fstream>
#include <memory>
#include <string>
#include <cmath>
#include <cstdint>

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include "DataFormats/L1Scouting/interface/L1ScoutingBMTFStub.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "L1Trigger/L1TMuonBarrel/interface/L1TMuonBarrelKalmanStubProcessor.h"

//For masks

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CondFormats/L1TObjects/interface/L1TMuonBarrelParams.h"
#include "CondFormats/DataRecord/interface/L1TMuonBarrelParamsRcd.h"
#include "CondFormats/L1TObjects/interface/L1MuDTTFMasks.h"
#include "CondFormats/DataRecord/interface/L1MuDTTFMasksRcd.h"


using namespace l1ScoutingRun3;

class ConverterStubsToFlatTable : public edm::global::EDProducer<> {
public:
  // constructor and destructor
  explicit ConverterStubsToFlatTable(const edm::ParameterSet&);
  ~ConverterStubsToFlatTable() override{};

  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // the tokens to access the data
  edm::EDGetTokenT<L1MuKBMTCombinedStubCollection> src_;

  std::string name_, doc_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ConverterStubsToFlatTable::ConverterStubsToFlatTable(const edm::ParameterSet& iConfig)
    : src_(consumes<L1MuKBMTCombinedStubCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<nanoaod::FlatTable>();
}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
void ConverterStubsToFlatTable::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  edm::Handle<L1MuKBMTCombinedStubCollection> src;
  iEvent.getByToken(src_, src);
  //auto out = std::make_unique<OrbitFlatTable>(src->bxOffsets(), name_);
  //out->setDoc(doc_);

  std::vector<int16_t> phi;
  std::vector<int16_t> phiB;
  std::vector<int16_t> qual;
  std::vector<int16_t> eta;
  std::vector<int16_t> qeta;
  std::vector<int16_t> station;
  std::vector<int16_t> wheel;
  std::vector<int16_t> sector;
  std::vector<int16_t> tag;
  std::vector<int16_t> bx;

  unsigned int i = 0;
  for (const L1MuKBMTCombinedStub& stub : *src) {
    phi.push_back(stub.phi());
    phiB.push_back(stub.phiB());
    qual.push_back(stub.quality());
    eta.push_back(stub.eta1());
    qeta.push_back(stub.qeta1());
    station.push_back(stub.stNum());
    wheel.push_back(stub.whNum());
    sector.push_back(stub.scNum());
    tag.push_back(stub.tag());
    bx.push_back(stub.bxNum());
    ++i;
  }

  auto out = std::make_unique<nanoaod::FlatTable>(tag.size(), name_, false, false);
  out->setDoc(doc_);

  out->addColumn<int16_t>("phi", phi, "phi (raw L1T units)");
  out->addColumn<int16_t>("phiB", phiB, "phiB (raw L1T units)");
  out->addColumn<int16_t>("qual", qual, "qual (raw L1T units)");
  out->addColumn<int16_t>("eta", eta, "eta (raw L1T units)");
  out->addColumn<int16_t>("qeta", qeta, "qeta (raw L1T units)");
  out->addColumn<int16_t>("station", station, "station (raw L1T units)");
  out->addColumn<int16_t>("wheel", wheel, "wheel (raw L1T units)");
  out->addColumn<int16_t>("sector", sector, "sector (raw L1T units)");
  out->addColumn<int16_t>("tag", tag, "tag (raw L1T units)");
  out->addColumn<int16_t>("bx", bx, "bx (raw L1T units)");

  iEvent.put(std::move(out));
}

void ConverterStubsToFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("src");
  desc.add<std::string>("name");
  desc.add<std::string>("doc");

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ConverterStubsToFlatTable);

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
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "L1Trigger/L1TMuonBarrel/interface/L1TMuonBarrelKalmanAlgo.h"
#include "L1Trigger/L1TMuonBarrel/interface/L1TMuonBarrelKalmanTrackFinder.h"
#include "L1Trigger/L1TMuonBarrel/interface/L1TMuonBarrelKalmanStubProcessor.h"

#include "DataFormats/L1TMuon/interface/L1MuKBMTrack.h"
#include "DataFormats/L1TMuon/interface/L1MuKBMTCombinedStub.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "CondFormats/L1TObjects/interface/L1TMuonGlobalParams.h"
#include "CondFormats/DataRecord/interface/L1TMuonGlobalParamsRcd.h"
#include "L1Trigger/L1TMuon/interface/L1TMuonGlobalParamsHelper.h"
#include "L1Trigger/L1TMuon/interface/MicroGMTLUTFactories.h"

#include "L1TriggerScouting/Utilities/interface/convertToL1TFormat.h"

using namespace l1ScoutingRun3;

class ConverterMuonsToFlatTable : public edm::stream::EDProducer<> {
public:
  // constructor and destructor
  explicit ConverterMuonsToFlatTable(const edm::ParameterSet&);
  ~ConverterMuonsToFlatTable() override;

  void produce(edm::Event&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  void beginRun(edm::Run const&, edm::EventSetup const&) override;
  void endRun(edm::Run const&, edm::EventSetup const&) override;

  // the tokens to access the data
  edm::EDGetTokenT<BXVector<l1t::Muon>> src_;

  std::string name_, doc_;

};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ConverterMuonsToFlatTable::ConverterMuonsToFlatTable(const edm::ParameterSet& iConfig)
    : src_(consumes<BXVector<l1t::Muon>>(iConfig.getParameter<edm::InputTag>("src"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")){
  produces<nanoaod::FlatTable>();

}
// -----------------------------------------------------------------------------


ConverterMuonsToFlatTable::~ConverterMuonsToFlatTable() {

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}


// ----------------------- method called for each orbit  -----------------------
void ConverterMuonsToFlatTable::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<BXVector<l1t::Muon>> src;
  iEvent.getByToken(src_, src);

  std::vector<float> pt;
  std::vector<float> eta;
  std::vector<float> phi;
  std::vector<int16_t> charge;
  std::vector<int16_t> quality;
  std::vector<int16_t> dxy;
  std::vector<int16_t> index;
  std::vector<float> ptUnconstrained;
  std::vector<float> etaAtVtx;
  std::vector<float> phiAtVtx;

  unsigned int i = 0;
  for (const l1t::Muon& mu : *src) {

    pt.push_back(mu.pt());
    eta.push_back(mu.eta());
    phi.push_back(mu.phi());
    charge.push_back(mu.charge());
    quality.push_back(mu.hwQual());
    dxy.push_back(mu.hwDXY());
    index.push_back(1); // wrong for now
    ptUnconstrained.push_back(mu.ptUnconstrained());
    etaAtVtx.push_back(mu.etaAtVtx());
    phiAtVtx.push_back(mu.phiAtVtx());

    /*pt[i] = ugmt::fPt(muon.hwPt());
    eta[i] = ugmt::fEta(muon.hwEta());
    phi[i] = ugmt::fPhi(muon.hwPhi());
    charge[i] = muon.hwCharge();
    quality[i] = muon.hwQual();
    dxy[i] = muon.hwDXY();
    index[i] = muon.tfMuonIndex();
    ptUnconstrained[i] = ugmt::fPtUnconstrained(muon.hwPtUnconstrained());
    etaAtVtx[i] = ugmt::fEtaAtVtx(muon.hwEtaAtVtx());
    phiAtVtx[i] = ugmt::fPhiAtVtx(muon.hwPhiAtVtx());*/

    ++i;
  }

  auto out = std::make_unique<nanoaod::FlatTable>(pt.size(), name_, false, false);
  out->setDoc(doc_);

  out->addColumn<float>("pt", pt, "pt (physical units)");
  out->addColumn<float>("eta", eta, "eta at second muon station (physical units)");
  out->addColumn<float>("phi", phi, "phi at second muon station (physical units)");
  out->addColumn<int16_t>("hwCharge", charge, "hwCharge (hw units)");
  out->addColumn<int16_t>("hwQual", quality, "hwQual (hw units)");
  out->addColumn<int16_t>("hwDXY", dxy, "untruncated transverse impact parameter (hw units)");
  out->addColumn<int16_t>("processor", index, "processor ([0-11])");
  out->addColumn<float>("ptUnconstrained", ptUnconstrained, "pt without vertex constraint (physical units)");
  out->addColumn<float>("etaAtVtx", etaAtVtx, "eta re-extrapolated at vertex (physical units)");
  out->addColumn<float>("phiAtVtx", phiAtVtx, "phi re-extrapolated at vertex (physical units)");

  iEvent.put(std::move(out));
}


// ------------ method called when starting to processes a run  ------------
void ConverterMuonsToFlatTable::beginRun(edm::Run const& run, edm::EventSetup const& iSetup) {

}

// ------------ method called when ending to processes a run  ------------
void ConverterMuonsToFlatTable::endRun(edm::Run const&, edm::EventSetup const&) {}

void ConverterMuonsToFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ConverterMuonsToFlatTable);

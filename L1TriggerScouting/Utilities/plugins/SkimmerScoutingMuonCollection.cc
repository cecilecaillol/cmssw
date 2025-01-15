#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "L1Trigger/L1TMuonBarrel/interface/L1TMuonBarrelKalmanAlgo.h"
#include "L1Trigger/L1TMuonBarrel/interface/L1TMuonBarrelKalmanTrackFinder.h"
#include "L1Trigger/L1TMuonBarrel/interface/L1TMuonBarrelKalmanStubProcessor.h"

#include "DataFormats/L1TMuon/interface/L1MuKBMTCombinedStub.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingMuon.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingBMTFStub.h"

#include "L1TriggerScouting/Utilities/interface/conversion.h"
#include "L1TriggerScouting/Utilities/interface/convertToL1TFormat.h"

using namespace l1ScoutingRun3;

//
// class declaration
//

class SkimmerScoutingMuonCollection : public edm::stream::EDProducer<> {
public:
  explicit SkimmerScoutingMuonCollection(const edm::ParameterSet&);
  ~SkimmerScoutingMuonCollection() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  edm::EDGetTokenT<MuonOrbitCollection> gmtSrc_;
  double ptmin_;
  double etamax_;
};
SkimmerScoutingMuonCollection::SkimmerScoutingMuonCollection(const edm::ParameterSet& iConfig)
    : gmtSrc_(consumes<MuonOrbitCollection>(iConfig.getParameter<edm::InputTag>("gmtSrc"))),
      ptmin_(iConfig.getParameter<double>("ptmin")),
      etamax_(iConfig.getParameter<double>("etamax")) {
  produces<MuonOrbitCollection>("L1MuonSkimmed").setBranchAlias("L1MuonSkimmedOrbitCollection");
}

SkimmerScoutingMuonCollection::~SkimmerScoutingMuonCollection() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void SkimmerScoutingMuonCollection::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  Handle<MuonOrbitCollection> muonsCollection;
  iEvent.getByToken(gmtSrc_, muonsCollection);

  std::vector<int> seenBxs;
  /*for (int i=0; i<muonsCollection->size(); ++i) {
    seenBxs.push_back((*muonsCollection)[i].bxNum());
  }*/
  for (int i=0; i<3565; ++i) {
    seenBxs.push_back(i);
  }

  std::unique_ptr<MuonOrbitCollection> skimmedMuonCollection(new MuonOrbitCollection);
  std::vector<std::vector<l1ScoutingRun3::Muon>> skimmedMuonBuffer(3565);
  unsigned nSkimmedMuon = 0;

  std::sort(seenBxs.begin(), seenBxs.end());
  seenBxs.erase(std::unique(seenBxs.begin(), seenBxs.end()), seenBxs.end());

  for (const auto& bx : seenBxs) {
      const auto& muons = muonsCollection->bxIterator(bx);
      for (const auto& muon : muons) {
          //const l1t::Muon gmt_m = getL1TMuon(muon);
	  if (ugmt::fPt(muon.hwPt())>ptmin_ and fabs(ugmt::fEta(muon.hwEta()))<etamax_){
	     const l1ScoutingRun3::Muon gmt_m = muon;
             skimmedMuonBuffer[bx].push_back(gmt_m);
             nSkimmedMuon++;
	  }
      }
  }

  skimmedMuonCollection->fillAndClear(skimmedMuonBuffer, nSkimmedMuon);
  iEvent.put(std::move(skimmedMuonCollection), "L1MuonSkimmed");
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void SkimmerScoutingMuonCollection::beginStream(edm::StreamID) {}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void SkimmerScoutingMuonCollection::endStream() {}

void SkimmerScoutingMuonCollection::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SkimmerScoutingMuonCollection);

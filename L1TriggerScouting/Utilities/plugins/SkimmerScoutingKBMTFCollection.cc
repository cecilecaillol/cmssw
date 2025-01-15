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

class SkimmerScoutingKBMTFCollection : public edm::stream::EDProducer<> {
public:
  explicit SkimmerScoutingKBMTFCollection(const edm::ParameterSet&);
  ~SkimmerScoutingKBMTFCollection() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  edm::EDGetTokenT<L1MuKBMTrackOrbitCollection> src_;
  L1TMuonBarrelKalmanAlgo* algo_;
  double ptmin_;
  double etamax_;
};
SkimmerScoutingKBMTFCollection::SkimmerScoutingKBMTFCollection(const edm::ParameterSet& iConfig)
    : src_(consumes<L1MuKBMTrackOrbitCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      algo_(new L1TMuonBarrelKalmanAlgo(iConfig.getParameter<edm::ParameterSet>("algoSettings"))),
      ptmin_(iConfig.getParameter<double>("ptmin")),
      etamax_(iConfig.getParameter<double>("etamax")) {
  produces<L1MuKBMTrackOrbitCollection>("L1MuKBMTrackSkimmed").setBranchAlias("L1MuKBMTrackSkimmedOrbitCollection");
}

SkimmerScoutingKBMTFCollection::~SkimmerScoutingKBMTFCollection() {
  if (algo_ != nullptr)
  delete algo_;

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void SkimmerScoutingKBMTFCollection::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  Handle<L1MuKBMTrackOrbitCollection> kbmtfCollection;
  iEvent.getByToken(src_, kbmtfCollection);

  std::vector<int> seenBxs;
  /*for (int i=0; i<muonsCollection->size(); ++i) {
    seenBxs.push_back((*muonsCollection)[i].bxNum());
  }*/
  for (int i=0; i<3565; ++i) {
    seenBxs.push_back(i);
  }

  std::unique_ptr<L1MuKBMTrackOrbitCollection> kbmTrackCollection(new L1MuKBMTrackOrbitCollection);
  std::vector<std::vector<L1MuKBMTrack>> kbmTrackBuffer(3565);
  unsigned nKbmTrack = 0;


  std::sort(seenBxs.begin(), seenBxs.end());
  seenBxs.erase(std::unique(seenBxs.begin(), seenBxs.end()), seenBxs.end());

  for (const auto& bx : seenBxs) {
      const auto& tracks = kbmtfCollection->bxIterator(bx);
      for (const auto& track : tracks) {
	   l1t::RegionalMuonCand bmtf_m = algo_->convertToBMTF(track);
	  if (ugmt::fPt(bmtf_m.hwPt())>ptmin_ and fabs(ugmt::fEta(bmtf_m.hwEta()))<etamax_){
	     const L1MuKBMTrack trk = track;
             kbmTrackBuffer[bx].push_back(trk);
             nKbmTrack++;
	  }
      }
  }
  kbmTrackCollection->fillAndClear(kbmTrackBuffer, nKbmTrack);
  iEvent.put(std::move(kbmTrackCollection), "L1MuKBMTrackSkimmed");

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void SkimmerScoutingKBMTFCollection::beginStream(edm::StreamID) {}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void SkimmerScoutingKBMTFCollection::endStream() {}

void SkimmerScoutingKBMTFCollection::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SkimmerScoutingKBMTFCollection);

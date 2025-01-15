// system include files

#include <memory>


// user include files

#include "FWCore/Framework/interface/global/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingBMTFStub.h"
#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"

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

//

// class declaration

//


class TrackFilter : public edm::global::EDFilter<> {

public:
  explicit TrackFilter(const edm::ParameterSet&);
  ~TrackFilter() override;

private:
  bool filter(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
    edm::EDGetTokenT<L1MuKBMTrackOrbitCollection> muSrc;

  // ----------member data ---------------------------

};


TrackFilter::TrackFilter(const edm::ParameterSet& iConfig)

    : 
    muSrc(consumes<L1MuKBMTrackOrbitCollection>(iConfig.getParameter<edm::InputTag>("muons")))
      {}


TrackFilter::~TrackFilter() {}


bool TrackFilter::filter(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {

    edm::Handle<L1MuKBMTrackOrbitCollection> src;
    iEvent.getByToken(muSrc, src);

    bool is_highptmu=false;

    unsigned int nmu = src->size();
    if (nmu>0) is_highptmu=true;

    /*for (unsigned int i = 0; i < nmu; ++i) {
        const pat::Muon & mu = (*muons)[i];
        if (fabs(mu.eta())<2.4 and mu.isPFMuon() and mu.isMediumMuon()){
	  if (mu.pt()>9) nmu9++;
          if (mu.pt()>18) nmu18++;
          if (mu.pt()>24) nmu24++;
	}
    }*/

    return is_highptmu;
}

#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(TrackFilter);

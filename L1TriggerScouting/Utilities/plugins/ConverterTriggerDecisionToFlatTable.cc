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

class ConverterTriggerDecisionToFlatTable : public edm::stream::EDProducer<> {
public:
  // constructor and destructor
  explicit ConverterTriggerDecisionToFlatTable(const edm::ParameterSet&);
  ~ConverterTriggerDecisionToFlatTable() override;

  void produce(edm::Event&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  void beginRun(edm::Run const&, edm::EventSetup const&) override;
  void endRun(edm::Run const&, edm::EventSetup const&) override;

  // the tokens to access the data
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;

  std::string name_, doc_;

};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ConverterTriggerDecisionToFlatTable::ConverterTriggerDecisionToFlatTable(const edm::ParameterSet& iConfig)
    : triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<nanoaod::FlatTable>();

}
// -----------------------------------------------------------------------------


ConverterTriggerDecisionToFlatTable::~ConverterTriggerDecisionToFlatTable() {

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}


// ----------------------- method called for each orbit  -----------------------
void ConverterTriggerDecisionToFlatTable::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::vector<int> HLTIsoMu24;
  int is_HLTIsoMu24=0;

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
     if (!(names.triggerName(i).find("HLT_IsoMu24_v") == string::npos) and triggerBits->accept(i)) is_HLTIsoMu24=1;
  }
  HLTIsoMu24.push_back(is_HLTIsoMu24);

  auto out = std::make_unique<nanoaod::FlatTable>(HLTIsoMu24.size(), name_, false, false);
  out->setDoc(doc_);

  out->addColumn<int>("HLTIsoMu24", HLTIsoMu24, "HLTIsoMu24");

  iEvent.put(std::move(out));
}


// ------------ method called when starting to processes a run  ------------
void ConverterTriggerDecisionToFlatTable::beginRun(edm::Run const& run, edm::EventSetup const& iSetup) {

}

// ------------ method called when ending to processes a run  ------------
void ConverterTriggerDecisionToFlatTable::endRun(edm::Run const&, edm::EventSetup const&) {}

void ConverterTriggerDecisionToFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ConverterTriggerDecisionToFlatTable);

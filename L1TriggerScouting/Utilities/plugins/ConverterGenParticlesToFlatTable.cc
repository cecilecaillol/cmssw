#include "FWCore/Framework/interface/MakerMacros.h"

#include <fstream>
#include <memory>
#include <string>
#include <cmath>
#include <cstdint>

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

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

class ConverterGenParticlesToFlatTable : public edm::stream::EDProducer<> {
public:
  // constructor and destructor
  explicit ConverterGenParticlesToFlatTable(const edm::ParameterSet&);
  ~ConverterGenParticlesToFlatTable() override;

  void produce(edm::Event&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  void beginRun(edm::Run const&, edm::EventSetup const&) override;
  void endRun(edm::Run const&, edm::EventSetup const&) override;

  // the tokens to access the data
  edm::EDGetTokenT<GenParticleCollection> src_;

  std::string name_, doc_;

};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ConverterGenParticlesToFlatTable::ConverterGenParticlesToFlatTable(const edm::ParameterSet& iConfig)
    : src_(consumes<GenParticleCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<nanoaod::FlatTable>();

}
// -----------------------------------------------------------------------------


ConverterGenParticlesToFlatTable::~ConverterGenParticlesToFlatTable() {

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}


// ----------------------- method called for each orbit  -----------------------
void ConverterGenParticlesToFlatTable::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::vector<float> pt;
  std::vector<float> eta;
  std::vector<float> phi;
  std::vector<int16_t> charge;
  std::vector<int16_t> pdgid;
  std::vector<float> beta;
  std::vector<float> mass;

  edm::Handle<GenParticleCollection> pruned;
  iEvent.getByToken(src_,pruned);
  for (size_t i = 0; i < pruned->size(); ++i) {
      const GenParticle &p = (*pruned)[i];
      const Candidate * mom = p.mother();
      //if (p.status()==1) cout<<p.pdgId()<<"    "<<p.pt()<<"    "<<p.eta()<<" "<<p.phi()<<" "<<mom->pdgId()<<endl;
      //if (i>1) cout<<p.pdgId()<<"    "<<p.pt()<<"    "<<p.eta()<<" "<<p.phi()<<" "<<mom->pdgId()<<" "<<p.status()<<endl;
      if (p.status()==1){
         pt.push_back(p.pt());
         eta.push_back(p.eta());
         phi.push_back(p.phi());
         beta.push_back(p.p()/p.energy());
	 mass.push_back(p.mass());
         charge.push_back(p.charge());
         pdgid.push_back(p.pdgId());
      }
  }
  //cout<<endl;

  auto out = std::make_unique<nanoaod::FlatTable>(pt.size(), name_, false, false);
  out->setDoc(doc_);

  out->addColumn<float>("pt", pt, "pt");
  out->addColumn<float>("eta", eta, "eta");
  out->addColumn<float>("phi", phi, "phi");
  out->addColumn<float>("beta", beta, "beta");
  out->addColumn<float>("mass", mass, "mass");
  out->addColumn<int16_t>("pdgid", pdgid, "pdgid");
  out->addColumn<int16_t>("charge", charge, "charge");

  iEvent.put(std::move(out));
}


// ------------ method called when starting to processes a run  ------------
void ConverterGenParticlesToFlatTable::beginRun(edm::Run const& run, edm::EventSetup const& iSetup) {

}

// ------------ method called when ending to processes a run  ------------
void ConverterGenParticlesToFlatTable::endRun(edm::Run const&, edm::EventSetup const&) {}

void ConverterGenParticlesToFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ConverterGenParticlesToFlatTable);

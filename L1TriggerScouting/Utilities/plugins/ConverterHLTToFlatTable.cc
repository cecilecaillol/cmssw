#include "FWCore/Framework/interface/MakerMacros.h"

#include <fstream>
#include <memory>
#include <string>
#include <cmath>
#include <cstdint>

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

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
//using namespace reco;
using namespace std;

class ConverterHLTToFlatTable : public edm::stream::EDProducer<> {
public:
  // constructor and destructor
  explicit ConverterHLTToFlatTable(const edm::ParameterSet&);
  ~ConverterHLTToFlatTable() override;

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

ConverterHLTToFlatTable::ConverterHLTToFlatTable(const edm::ParameterSet& iConfig)
    : triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<nanoaod::FlatTable>();

}
// -----------------------------------------------------------------------------


ConverterHLTToFlatTable::~ConverterHLTToFlatTable() {

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}


// ----------------------- method called for each orbit  -----------------------
void ConverterHLTToFlatTable::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::vector<int16_t> v_HLT_IsoMu24;
  std::vector<int16_t> v_HLT_Mu50;
  std::vector<int16_t> v_HLT_PFMET120_PFMHT120_IDTight;
  std::vector<int16_t> v_HLT_PFHT500_PFMET100_PFMHT100_IDTight;
  std::vector<int16_t> v_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
  std::vector<int16_t> v_HLT_MET105_IsoTrk50;
  int HLT_IsoMu24=0;
  int HLT_Mu50=0;
  int HLT_PFMET120_PFMHT120_IDTight=0;
  int HLT_PFHT500_PFMET100_PFMHT100_IDTight=0;
  int HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60=0;
  int HLT_MET105_IsoTrk50=0;

  /*edm::Handle<edm::TriggerResults> handle;
  iEvent.getByToken(m_token, handle);
  const edm::TriggerResults& triggers = *handle;
  const edm::TriggerNames& names = triggerNames(triggers);*/

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);

  const edm::TriggerNames names = iEvent.triggerNames(*triggerBits);

  for (size_t i_hlt = 0; i_hlt != triggerBits->size(); ++i_hlt){
      string hltName = names.triggerName(i_hlt);
      if(!(hltName.find("HLT_IsoMu24_v") == string::npos)){ 
         if( triggerBits->wasrun(i_hlt) && !triggerBits->error(i_hlt) && triggerBits->accept(i_hlt )) HLT_IsoMu24 = 1;
      }
      if(!(hltName.find("HLT_Mu50_v") == string::npos)){
         if( triggerBits->wasrun(i_hlt) && !triggerBits->error(i_hlt) && triggerBits->accept(i_hlt )) HLT_Mu50 = 1;
      }
      if(!(hltName.find("HLT_PFMET120_PFMHT120_IDTight_v") == string::npos)){
         if( triggerBits->wasrun(i_hlt) && !triggerBits->error(i_hlt) && triggerBits->accept(i_hlt )) HLT_PFMET120_PFMHT120_IDTight = 1;
      }
      if(!(hltName.find("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v") == string::npos)){
         if( triggerBits->wasrun(i_hlt) && !triggerBits->error(i_hlt) && triggerBits->accept(i_hlt )) HLT_PFHT500_PFMET100_PFMHT100_IDTight = 1;
      }
      if(!(hltName.find("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v") == string::npos)){
         if( triggerBits->wasrun(i_hlt) && !triggerBits->error(i_hlt) && triggerBits->accept(i_hlt )) HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = 1;
      }
      if(!(hltName.find("HLT_MET105_IsoTrk50_v") == string::npos)){
         if( triggerBits->wasrun(i_hlt) && !triggerBits->error(i_hlt) && triggerBits->accept(i_hlt )) HLT_MET105_IsoTrk50 = 1;
      }
  }
  v_HLT_IsoMu24.push_back(HLT_IsoMu24);
  v_HLT_Mu50.push_back(HLT_Mu50);
  v_HLT_PFMET120_PFMHT120_IDTight.push_back(HLT_PFMET120_PFMHT120_IDTight);
  v_HLT_PFHT500_PFMET100_PFMHT100_IDTight.push_back(HLT_PFHT500_PFMET100_PFMHT100_IDTight);
  v_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60.push_back(HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
  v_HLT_MET105_IsoTrk50.push_back(HLT_MET105_IsoTrk50);

  //std::cout << "\n == TRIGGER PATHS= " << std::endl;
  /*for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
        std::cout << "Trigger " << names.triggerName(i) <<
                ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
                ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
                << std::endl;
  }*/

  /*edm::Handle<GenParticleCollection> pruned;
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
         charge.push_back(p.charge());
         pdgid.push_back(p.pdgId());
      }
  }*/
  //cout<<endl;

  auto out = std::make_unique<nanoaod::FlatTable>(v_HLT_IsoMu24.size(), name_, false, false);
  out->setDoc(doc_);

  out->addColumn<int>("IsoMu24", v_HLT_IsoMu24, "IsoMu24");
  out->addColumn<int>("Mu50", v_HLT_Mu50, "Mu50");
  out->addColumn<int>("PFMET120_PFMHT120_IDTight", v_HLT_PFMET120_PFMHT120_IDTight, "PFMET120_PFMHT120_IDTight");
  out->addColumn<int>("PFHT500_PFMET100_PFMHT100_IDTight", v_HLT_PFHT500_PFMET100_PFMHT100_IDTight, "PFHT500_PFMET100_PFMHT100_IDTight");
  out->addColumn<int>("PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", v_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, "PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60");
  out->addColumn<int>("MET105_IsoTrk50", v_HLT_MET105_IsoTrk50, "MET105_IsoTrk50");

  iEvent.put(std::move(out));
}


// ------------ method called when starting to processes a run  ------------
void ConverterHLTToFlatTable::beginRun(edm::Run const& run, edm::EventSetup const& iSetup) {

}

// ------------ method called when ending to processes a run  ------------
void ConverterHLTToFlatTable::endRun(edm::Run const&, edm::EventSetup const&) {}

void ConverterHLTToFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ConverterHLTToFlatTable);

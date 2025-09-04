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

class ConverterKbmtfTracksToFlatTable : public edm::stream::EDProducer<> {
public:
  // constructor and destructor
  explicit ConverterKbmtfTracksToFlatTable(const edm::ParameterSet&);
  ~ConverterKbmtfTracksToFlatTable() override;

  void produce(edm::Event&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  void beginRun(edm::Run const&, edm::EventSetup const&) override;
  void endRun(edm::Run const&, edm::EventSetup const&) override;

  unsigned calcGlobalPhi(l1t::RegionalMuonCand&);

  // the tokens to access the data
  edm::EDGetTokenT<BXVector<L1MuKBMTrack>> src_;

  std::string name_, doc_;

  L1TMuonBarrelKalmanAlgo* algo_;

  bool addStubs_;

  std::shared_ptr<l1t::MicroGMTExtrapolationLUT> m_BEtaExtrapolation_;
  std::shared_ptr<l1t::MicroGMTExtrapolationLUT> m_BPhiExtrapolation_;

  std::unique_ptr<L1TMuonGlobalParamsHelper> microGMTParamsHelper_;
  edm::ESGetToken<L1TMuonGlobalParams, L1TMuonGlobalParamsRcd> m_microGMTParamsToken_;

  int fwRev_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ConverterKbmtfTracksToFlatTable::ConverterKbmtfTracksToFlatTable(const edm::ParameterSet& iConfig)
    : src_(consumes<BXVector<L1MuKBMTrack>>(iConfig.getParameter<edm::InputTag>("src"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")),
      algo_(new L1TMuonBarrelKalmanAlgo(iConfig.getParameter<edm::ParameterSet>("algoSettings"))),
      addStubs_(iConfig.getParameter<bool>("addStubs")) {
  produces<nanoaod::FlatTable>();

  m_microGMTParamsToken_ = esConsumes<L1TMuonGlobalParams, L1TMuonGlobalParamsRcd, edm::Transition::BeginRun>();
  microGMTParamsHelper_ = std::make_unique<L1TMuonGlobalParamsHelper>();
  fwRev_ = 0x8010000;
}
// -----------------------------------------------------------------------------


ConverterKbmtfTracksToFlatTable::~ConverterKbmtfTracksToFlatTable() {
  if (algo_ != nullptr)
    delete algo_;

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}


// ----------------------- method called for each orbit  -----------------------
void ConverterKbmtfTracksToFlatTable::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<BXVector<L1MuKBMTrack>> src;
  iEvent.getByToken(src_, src);
  //auto out = std::make_unique<nanoaod::FlatTable>(1, name_, false, false);
  //out->setDoc(doc_);

  int outputShiftPhi = 3;
  int outputShiftEta = 3;
  if (fwRev_ >= 0x4010000) {
    outputShiftPhi = 2;
    outputShiftEta = 0;
  }

  /*std::vector<float> pt(out->size());
  std::vector<float> eta(out->size());
  std::vector<float> phi(out->size());
  std::vector<int16_t> charge(out->size());
  std::vector<int16_t> quality(out->size());
  std::vector<int16_t> dxy(out->size());
  std::vector<int16_t> index(out->size());
  std::vector<float> ptUnconstrained(out->size());
  std::vector<float> etaAtVtx(out->size());
  std::vector<float> phiAtVtx(out->size());*/
  std::vector<float> pt;
  std::vector<float> eta;
  std::vector<float> phi;
  std::vector<int16_t> charge;
  std::vector<int16_t> quality;
  std::vector<int16_t> dxy;
  std::vector<int16_t> curvature;
  std::vector<int16_t> index;
  std::vector<float> ptUnconstrained;
  std::vector<float> etaAtVtx;
  std::vector<float> phiAtVtx;

  /*std::vector<int16_t> nStub(out->size());
  std::vector<std::vector<int16_t>> sStation(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sSector(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sWheel(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sHwQual(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sHwPhi(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sHwPhiB(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sHwEta1(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sHwQEta1(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sHwEta2(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sHwQEta2(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sTag(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sBx(4, std::vector<int16_t>(out->size(), 0));*/

  std::vector<int16_t> nStub;

  std::vector<int16_t> s1Station;
  std::vector<int16_t> s1Sector;
  std::vector<int16_t> s1Wheel;
  std::vector<int16_t> s1HwQual;
  std::vector<int16_t> s1HwPhi;
  std::vector<int16_t> s1HwPhiB;
  std::vector<int16_t> s1HwEta1;
  std::vector<int16_t> s1HwQEta1;
  std::vector<int16_t> s1HwEta2;
  std::vector<int16_t> s1HwQEta2;
  std::vector<int16_t> s1Tag;
  std::vector<int16_t> s1Bx;

  std::vector<int16_t> s2Station;
  std::vector<int16_t> s2Sector;
  std::vector<int16_t> s2Wheel;
  std::vector<int16_t> s2HwQual;
  std::vector<int16_t> s2HwPhi;
  std::vector<int16_t> s2HwPhiB;
  std::vector<int16_t> s2HwEta1;
  std::vector<int16_t> s2HwQEta1;
  std::vector<int16_t> s2HwEta2;
  std::vector<int16_t> s2HwQEta2;
  std::vector<int16_t> s2Tag;
  std::vector<int16_t> s2Bx;

  std::vector<int16_t> s3Station;
  std::vector<int16_t> s3Sector;
  std::vector<int16_t> s3Wheel;
  std::vector<int16_t> s3HwQual;
  std::vector<int16_t> s3HwPhi;
  std::vector<int16_t> s3HwPhiB;
  std::vector<int16_t> s3HwEta1;
  std::vector<int16_t> s3HwQEta1;
  std::vector<int16_t> s3HwEta2;
  std::vector<int16_t> s3HwQEta2;
  std::vector<int16_t> s3Tag;
  std::vector<int16_t> s3Bx;

  std::vector<int16_t> s4Station;
  std::vector<int16_t> s4Sector;
  std::vector<int16_t> s4Wheel;
  std::vector<int16_t> s4HwQual;
  std::vector<int16_t> s4HwPhi;
  std::vector<int16_t> s4HwPhiB;
  std::vector<int16_t> s4HwEta1;
  std::vector<int16_t> s4HwQEta1;
  std::vector<int16_t> s4HwEta2;
  std::vector<int16_t> s4HwQEta2;
  std::vector<int16_t> s4Tag;
  std::vector<int16_t> s4Bx;

  unsigned int i = 0;
  for (const L1MuKBMTrack& track : *src) {

    l1t::RegionalMuonCand bmtf_m = algo_->convertToBMTF(track);
    /*pt[i] = ugmt::fPt(bmtf_m.hwPt());
    eta[i] = ugmt::fEta(bmtf_m.hwEta());
    phi[i] = ugmt::fPhi(calcGlobalPhi(bmtf_m));
    charge[i] = bmtf_m.hwSign()==1? -1 : 1;
    quality[i] = bmtf_m.hwQual();
    dxy[i] = track.dxy();
    index[i] = bmtf_m.processor(); // wrong for now
    ptUnconstrained[i] = ugmt::fPtUnconstrained(bmtf_m.hwPtUnconstrained());*/
    pt.push_back(ugmt::fPt(bmtf_m.hwPt()));
    eta.push_back(ugmt::fEta(bmtf_m.hwEta()));
    phi.push_back(ugmt::fPhi(calcGlobalPhi(bmtf_m)));
    charge.push_back(bmtf_m.hwSign()==1? -1 : 1);
    quality.push_back(bmtf_m.hwQual());
    //dxy.push_back(track.dxy()); // Do we need this quantity?
    dxy.push_back(bmtf_m.hwDXY());
    curvature.push_back(bmtf_m.hwK());
    index.push_back(bmtf_m.processor()); // wrong for now
    ptUnconstrained.push_back(ugmt::fPtUnconstrained(bmtf_m.hwPtUnconstrained()));

    int ptRedInWidth = m_BPhiExtrapolation_->getPtRedInWidth();
    int ptMask = (1 << ptRedInWidth) - 1;
    int etaRedInWidth = m_BPhiExtrapolation_->getEtaRedInWidth();
    int redEtaShift = 8 - etaRedInWidth;

    // only use LSBs of pt:
    int ptRed = bmtf_m.hwPt() & ptMask;
    // here we drop the LSBs and mask the MSB
    int etaAbsRed = (std::abs(bmtf_m.hwEta()) >> redEtaShift) & ((1 << etaRedInWidth) - 1);
    int deltaPhi = 0;
    int deltaEta = 0;

    if (bmtf_m.hwPt() < (1 << ptRedInWidth)) {  // extrapolation only for "low" pT muons
      int sign = 1;
      if (bmtf_m.hwSign() == 1) {
        sign = -1;
      }
      deltaPhi = (m_BPhiExtrapolation_->lookup(etaAbsRed, ptRed) << outputShiftPhi) * sign;
      deltaEta = (m_BEtaExtrapolation_->lookup(etaAbsRed, ptRed) << outputShiftEta);
      if (bmtf_m.hwEta() > 0) {
        deltaEta *= -1;
      }
    }

    /*etaAtVtx[i] = ugmt::fEta(bmtf_m.hwEta() + deltaEta);
    phiAtVtx[i] = ugmt::fPhi(calcGlobalPhi(bmtf_m) + deltaPhi);*/

    etaAtVtx.push_back(ugmt::fEta(bmtf_m.hwEta() + deltaEta));
    phiAtVtx.push_back(ugmt::fPhi(calcGlobalPhi(bmtf_m) + deltaPhi));

    if (addStubs_) {
      /*nStub[i] = track.stubs().size();
      unsigned j = 0;
      for (const auto& stub : track.stubs()) {
        sStation[j][i] = (*stub).stNum();
        sSector[j][i] = (*stub).scNum();
        sWheel[j][i] = (*stub).whNum();
        sHwQual[j][i] = (*stub).quality();
        sHwPhi[j][i] = (*stub).phi();
        sHwPhiB[j][i] = (*stub).phiB();
        sHwEta1[j][i] = (*stub).eta1();
        sHwQEta1[j][i] = (*stub).qeta1();
        sHwEta2[j][i] = (*stub).eta2();
        sHwQEta2[j][i] = (*stub).qeta2();
        sTag[j][i] = (*stub).tag();
	sBx[j][i] = (*stub).bxNum();
        ++j;
      }*/
      unsigned j = 0;
      nStub.push_back(track.stubs().size());
      for (const auto& stub : track.stubs()) {
	++j;
	if (j==1){
           s1Station.push_back((*stub).stNum());
           s1Sector.push_back((*stub).scNum());
           s1Wheel.push_back((*stub).whNum());
           s1HwQual.push_back((*stub).quality());
           s1HwPhi.push_back((*stub).phi());
           s1HwPhiB.push_back((*stub).phiB());
           s1HwEta1.push_back((*stub).eta1());
           s1HwQEta1.push_back((*stub).qeta1());
           s1HwEta2.push_back((*stub).eta2());
           s1HwQEta2.push_back((*stub).qeta2());
           s1Tag.push_back((*stub).tag());
           s1Bx.push_back((*stub).bxNum());
	}
	if (j==2){
           s2Station.push_back((*stub).stNum());
           s2Sector.push_back((*stub).scNum());
           s2Wheel.push_back((*stub).whNum());
           s2HwQual.push_back((*stub).quality());
           s2HwPhi.push_back((*stub).phi());
           s2HwPhiB.push_back((*stub).phiB());
           s2HwEta1.push_back((*stub).eta1());
           s2HwQEta1.push_back((*stub).qeta1());
           s2HwEta2.push_back((*stub).eta2());
           s2HwQEta2.push_back((*stub).qeta2());
           s2Tag.push_back((*stub).tag());
           s2Bx.push_back((*stub).bxNum());
        }
	if (j==3){
           s3Station.push_back((*stub).stNum());
           s3Sector.push_back((*stub).scNum());
           s3Wheel.push_back((*stub).whNum());
           s3HwQual.push_back((*stub).quality());
           s3HwPhi.push_back((*stub).phi());
           s3HwPhiB.push_back((*stub).phiB());
           s3HwEta1.push_back((*stub).eta1());
           s3HwQEta1.push_back((*stub).qeta1());
           s3HwEta2.push_back((*stub).eta2());
           s3HwQEta2.push_back((*stub).qeta2());
           s3Tag.push_back((*stub).tag());
           s3Bx.push_back((*stub).bxNum());
        }
	if (j==4){
           s4Station.push_back((*stub).stNum());
           s4Sector.push_back((*stub).scNum());
           s4Wheel.push_back((*stub).whNum());
           s4HwQual.push_back((*stub).quality());
           s4HwPhi.push_back((*stub).phi());
           s4HwPhiB.push_back((*stub).phiB());
           s4HwEta1.push_back((*stub).eta1());
           s4HwQEta1.push_back((*stub).qeta1());
           s4HwEta2.push_back((*stub).eta2());
           s4HwQEta2.push_back((*stub).qeta2());
           s4Tag.push_back((*stub).tag());
           s4Bx.push_back((*stub).bxNum());
        }
      }
      if (track.stubs().size()<4){
	s4Station.push_back(-1);
        s4Sector.push_back(-1);
        s4Wheel.push_back(-1);
        s4HwQual.push_back(-1);
        s4HwPhi.push_back(-1);
        s4HwPhiB.push_back(-1);
        s4HwEta1.push_back(-1);
        s4HwQEta1.push_back(-1);
        s4HwEta2.push_back(-1);
        s4HwQEta2.push_back(-1);
        s4Tag.push_back(-1);
        s4Bx.push_back(-99);
      }
      if (track.stubs().size()<3){
        s3Station.push_back(-1);
        s3Sector.push_back(-1);
        s3Wheel.push_back(-1);
        s3HwQual.push_back(-1);
        s3HwPhi.push_back(-1);
        s3HwPhiB.push_back(-1);
        s3HwEta1.push_back(-1);
        s3HwQEta1.push_back(-1);
        s3HwEta2.push_back(-1);
        s3HwQEta2.push_back(-1);
        s3Tag.push_back(-1);
        s3Bx.push_back(-99);
      }
    }

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
  out->addColumn<int16_t>("hwK", curvature, "curvature");
  out->addColumn<int16_t>("processor", index, "processor ([0-11])");
  out->addColumn<float>("ptUnconstrained", ptUnconstrained, "pt without vertex constraint (physical units)");
  out->addColumn<float>("etaAtVtx", etaAtVtx, "eta re-extrapolated at vertex (physical units)");
  out->addColumn<float>("phiAtVtx", phiAtVtx, "phi re-extrapolated at vertex (physical units)");

  if (addStubs_) {
    out->addColumn<int16_t>("nStub", nStub, "number of stubs used to reconstruct KBMTF track");
    /*for (int i=0; i<4; ++i) {
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"Station", sStation[i], "stub station");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"Sector", sSector[i], "stub sector");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"Wheel", sWheel[i], "stub wheel");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"HwQual", sHwQual[i], "stub quality (hw units)");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"HwPhi", sHwPhi[i], "stub local phi position (hw units)");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"HwPhiB", sHwPhiB[i], "stub phi bending (hw units)");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"HwEta1", sHwEta1[i], "eta of first stub in chamber (hw units)");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"HwQEta1", sHwQEta1[i], "eta quality of first stub in chamber (hw units)");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"HwEta2", sHwEta2[i], "eta of second stub in chamber (hw units)");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"HwQEta2", sHwQEta2[i], "eta quality of second stub in chamber (hw units)");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"Tag", sTag[i], "tag=0 is for second stub in chamber");
	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"Bx", sBx[i], "bx");
    }*/
    out->addColumn<int16_t>("s1Station", s1Station, "stub station");
    out->addColumn<int16_t>("s1Sector", s1Sector, "stub sector");
    out->addColumn<int16_t>("s1Wheel", s1Wheel, "stub wheel");
    /*out->addColumn<int16_t>("s1HwQual", s1HwQual, "stub quality (hw units)");
    out->addColumn<int16_t>("s1HwPhi", s1HwPhi, "stub local phi position (hw units)");
    out->addColumn<int16_t>("s1HwPhiB", s1HwPhiB, "stub phi bending (hw units)");
    out->addColumn<int16_t>("s1HwEta1", s1HwEta1, "eta of first stub in chamber (hw units)");
    out->addColumn<int16_t>("s1HwQEta1", s1HwQEta1, "eta quality of first stub in chamber (hw units)");
    out->addColumn<int16_t>("s1HwEta2", s1HwEta2, "eta of second stub in chamber (hw units)");
    out->addColumn<int16_t>("s1HwQEta2", s1HwQEta2, "eta quality of second stub in chamber (hw units)");
    out->addColumn<int16_t>("s1Tag", s1Tag, "tag=0 is for second stub in chamber");*/
    out->addColumn<int16_t>("s1Bx", s1Bx, "bx");

    out->addColumn<int16_t>("s2Station", s2Station, "stub station");
    out->addColumn<int16_t>("s2Sector", s2Sector, "stub sector");
    out->addColumn<int16_t>("s2Wheel", s2Wheel, "stub wheel");
    /*out->addColumn<int16_t>("s2HwQual", s2HwQual, "stub quality (hw units)");
    out->addColumn<int16_t>("s2HwPhi", s2HwPhi, "stub local phi position (hw units)");
    out->addColumn<int16_t>("s2HwPhiB", s2HwPhiB, "stub phi bending (hw units)");
    out->addColumn<int16_t>("s2HwEta1", s2HwEta1, "eta of first stub in chamber (hw units)");
    out->addColumn<int16_t>("s2HwQEta1", s2HwQEta1, "eta quality of first stub in chamber (hw units)");
    out->addColumn<int16_t>("s2HwEta2", s2HwEta2, "eta of second stub in chamber (hw units)");
    out->addColumn<int16_t>("s2HwQEta2", s2HwQEta2, "eta quality of second stub in chamber (hw units)");
    out->addColumn<int16_t>("s2Tag", s2Tag, "tag=0 is for second stub in chamber");*/
    out->addColumn<int16_t>("s2Bx", s2Bx, "bx");

    out->addColumn<int16_t>("s3Station", s3Station, "stub station");
    out->addColumn<int16_t>("s3Sector", s3Sector, "stub sector");
    out->addColumn<int16_t>("s3Wheel", s3Wheel, "stub wheel");
    /*out->addColumn<int16_t>("s3HwQual", s3HwQual, "stub quality (hw units)");
    out->addColumn<int16_t>("s3HwPhi", s3HwPhi, "stub local phi position (hw units)");
    out->addColumn<int16_t>("s3HwPhiB", s3HwPhiB, "stub phi bending (hw units)");
    out->addColumn<int16_t>("s3HwEta1", s3HwEta1, "eta of first stub in chamber (hw units)");
    out->addColumn<int16_t>("s3HwQEta1", s3HwQEta1, "eta quality of first stub in chamber (hw units)");
    out->addColumn<int16_t>("s3HwEta2", s3HwEta2, "eta of second stub in chamber (hw units)");
    out->addColumn<int16_t>("s3HwQEta2", s3HwQEta2, "eta quality of second stub in chamber (hw units)");
    out->addColumn<int16_t>("s3Tag", s3Tag, "tag=0 is for second stub in chamber");*/
    out->addColumn<int16_t>("s3Bx", s3Bx, "bx");

    out->addColumn<int16_t>("s4Station", s4Station, "stub station");
    out->addColumn<int16_t>("s4Sector", s4Sector, "stub sector");
    out->addColumn<int16_t>("s4Wheel", s4Wheel, "stub wheel");
    /*out->addColumn<int16_t>("s4HwQual", s4HwQual, "stub quality (hw units)");
    out->addColumn<int16_t>("s4HwPhi", s4HwPhi, "stub local phi position (hw units)");
    out->addColumn<int16_t>("s4HwPhiB", s4HwPhiB, "stub phi bending (hw units)");
    out->addColumn<int16_t>("s4HwEta1", s4HwEta1, "eta of first stub in chamber (hw units)");
    out->addColumn<int16_t>("s4HwQEta1", s4HwQEta1, "eta quality of first stub in chamber (hw units)");
    out->addColumn<int16_t>("s4HwEta2", s4HwEta2, "eta of second stub in chamber (hw units)");
    out->addColumn<int16_t>("s4HwQEta2", s4HwQEta2, "eta quality of second stub in chamber (hw units)");
    out->addColumn<int16_t>("s4Tag", s4Tag, "tag=0 is for second stub in chamber");*/
    out->addColumn<int16_t>("s4Bx", s4Bx, "bx");
  }

  iEvent.put(std::move(out));
}


// ------------ method called when starting to processes a run  ------------
void ConverterKbmtfTracksToFlatTable::beginRun(edm::Run const& run, edm::EventSetup const& iSetup) {

  /*
  edm::ESHandle<L1TMuonGlobalParams> microGMTParamsHandle = iSetup.getHandle(m_microGMTParamsToken);

  std::unique_ptr<L1TMuonGlobalParams_PUBLIC> microGMTParams(
      new L1TMuonGlobalParams_PUBLIC(cast_to_L1TMuonGlobalParams_PUBLIC(*microGMTParamsHandle.product())));
  if (microGMTParams->pnodes_.empty()) {
    edm::ESHandle<L1TMuonGlobalParams> o2oProtoHandle = iSetup.getHandle(m_o2oProtoToken);
    microGMTParamsHelper = std::make_unique<L1TMuonGlobalParamsHelper>(*o2oProtoHandle.product());
  } else
    microGMTParamsHelper =
        std::make_unique<L1TMuonGlobalParamsHelper>(cast_to_L1TMuonGlobalParams(*microGMTParams.get()));
  */

  edm::ESHandle<L1TMuonGlobalParams> microGMTParamsHandle_ = iSetup.getHandle(m_microGMTParamsToken_);
  microGMTParamsHelper_ = std::make_unique<L1TMuonGlobalParamsHelper>(*microGMTParamsHandle_.product());
  if (!microGMTParamsHelper_) {
    edm::LogError("L1TMicroGMTLUTDumper") << "Could not retrieve parameters from Event Setup" << std::endl;
  }
  
  // int fwRev_ = 0x8010000; // microGMTParamsHelper_->fwVersion();

  m_BEtaExtrapolation_ = l1t::MicroGMTExtrapolationLUTFactory::create(microGMTParamsHelper_->bEtaExtrapolationLUT(), l1t::MicroGMTConfiguration::ETA_OUT, fwRev_);
  m_BPhiExtrapolation_ = l1t::MicroGMTExtrapolationLUTFactory::create(microGMTParamsHelper_->bPhiExtrapolationLUT(), l1t::MicroGMTConfiguration::PHI_OUT, fwRev_);
}

// ------------ method called when ending to processes a run  ------------
void ConverterKbmtfTracksToFlatTable::endRun(edm::Run const&, edm::EventSetup const&) {}

void ConverterKbmtfTracksToFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

unsigned ConverterKbmtfTracksToFlatTable::calcGlobalPhi(l1t::RegionalMuonCand& l1_reg_m) {

  unsigned globalPhi = l1_reg_m.processor()*48 + l1_reg_m.hwPhi();
  globalPhi += 576 - 24;      // first processor starts at -15degrees in cms phi
  globalPhi = globalPhi%576;  // wrap around

  return globalPhi;
}

DEFINE_FWK_MODULE(ConverterKbmtfTracksToFlatTable);

// system include files
#include <memory>
#include <algorithm>
#include <fstream>
#include <cstdio>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TCorrelator/interface/TkPrimaryVertex.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "L1Trigger/Phase2L1ParticleFlow/src/newfirmware/dataformats/layer1_emulator.h"
#include "L1Trigger/Phase2L1ParticleFlow/src/newfirmware/regionizer/regionizer_base_ref.h"
#include "L1Trigger/Phase2L1ParticleFlow/src/newfirmware/pf/pfalgo2hgc_ref.h"
#include "L1Trigger/Phase2L1ParticleFlow/src/newfirmware/pf/pfalgo3_ref.h"
#include "L1Trigger/Phase2L1ParticleFlow/src/newfirmware/pf/pfalgo_dummy_ref.h"
#include "L1Trigger/Phase2L1ParticleFlow/src/newfirmware/puppi/linpuppi_ref.h"
#include "L1Trigger/Phase2L1ParticleFlow/src/newfirmware/egamma/pftkegalgo_ref.h"

#include "DataFormats/L1TCorrelator/interface/TkMuon.h"
#include "DataFormats/L1TCorrelator/interface/TkMuonFwd.h"

#include "DataFormats/L1TCorrelator/interface/TkElectron.h"
#include "DataFormats/L1TCorrelator/interface/TkElectronFwd.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1TCorrelator/interface/TkEm.h"
#include "DataFormats/L1TCorrelator/interface/TkEmFwd.h"

//--------------------------------------------------------------------------------------------------
class L1TCorrelatorLayer1Producer : public edm::stream::EDProducer<> {
public:
  explicit L1TCorrelatorLayer1Producer(const edm::ParameterSet&);
  ~L1TCorrelatorLayer1Producer() override;

private:
  edm::ParameterSet config_;
  int debug_;

  bool useStandaloneMuons_;
  bool useTrackerMuons_;

  bool hasTracks_;
  edm::EDGetTokenT<l1t::PFTrackCollection> tkCands_;
  float trkPt_;
  edm::EDGetTokenT<std::vector<l1t::TkPrimaryVertex>> extTkVtx_;

  edm::EDGetTokenT<l1t::MuonBxCollection> muCands_;    // standalone muons
  edm::EDGetTokenT<l1t::TkMuonCollection> tkMuCands_;  // tk muons

  std::vector<edm::EDGetTokenT<l1t::PFClusterCollection>> emCands_;
  std::vector<edm::EDGetTokenT<l1t::PFClusterCollection>> hadCands_;

  float emPtCut_, hadPtCut_;

  bool sortOutputs_;

  l1ct::Event event_; 
  std::unique_ptr<l1ct::RegionizerEmulator> regionizer_;
  std::unique_ptr<l1ct::PFAlgoEmulatorBase> l1pfalgo_;
  std::unique_ptr<l1ct::LinPuppiEmulator> l1pualgo_;
  std::unique_ptr<l1ct::PFTkEGAlgoEmulator> l1tkegalgo_;

  bool writeEgSta_;
  // Region dump
  const std::string regionDumpName_;
  std::fstream fRegionDump_;

  // region of interest debugging
  float debugEta_, debugPhi_, debugR_;

  // these are used to link items back
  std::unordered_map<const l1t::PFCluster *, l1t::PFClusterRef> clusterRefMap_;
  std::unordered_map<const l1t::PFTrack *, l1t::PFTrackRef> trackRefMap_;
  std::unordered_map<const l1t::Muon *, l1t::PFCandidate::MuonRef> muonRefMap_;

  // main methods
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void addUInt(unsigned int value, std::string iLabel, edm::Event& iEvent);
 
  void initSectorsAndRegions(const edm::ParameterSet & iConfig);
  void initEvent(const edm::Event &e);
  // add object, tracking references
  void addTrack(const l1t::PFTrack &t, l1t::PFTrackRef ref);
  void addMuon(const l1t::Muon &t, l1t::PFCandidate::MuonRef ref);
  void addHadCalo(const l1t::PFCluster &t, l1t::PFClusterRef ref);
  void addEmCalo(const l1t::PFCluster &t, l1t::PFClusterRef ref);
  // add objects in already-decoded format
  void addDecodedTrack(l1ct::DetectorSector<l1ct::TkObjEmu> & sec, const l1t::PFTrack &t);
  void addDecodedMuon(l1ct::DetectorSector<l1ct::MuObjEmu> & sec, const l1t::Muon &t);
  void addDecodedHadCalo(l1ct::DetectorSector<l1ct::HadCaloObjEmu> & sec, const l1t::PFCluster &t);
  void addDecodedEmCalo(l1ct::DetectorSector<l1ct::EmCaloObjEmu> & sec, const l1t::PFCluster &t);
  // fetching outputs
  std::unique_ptr<l1t::PFCandidateCollection> fetchHadCalo() const;
  std::unique_ptr<l1t::PFCandidateCollection> fetchEmCalo() const;
  std::unique_ptr<l1t::PFCandidateCollection> fetchTracks() const;
  std::unique_ptr<l1t::PFCandidateCollection> fetchPF() const;
  std::unique_ptr<l1t::PFCandidateCollection> fetchPuppi() const;

  void putEgObjects(edm::Event &iEvent,
                    const bool writeEgSta,
                    const std::string &egLablel,
                    const std::string &tkEmLabel,
                    const std::string &tkEleLabel) const;
                    
  template<typename T>
  void setRefs_(l1t::PFCandidate &pf, const T & p) const ; 
};

//
// constructors and destructor
//
L1TCorrelatorLayer1Producer::L1TCorrelatorLayer1Producer(const edm::ParameterSet& iConfig)
    : config_(iConfig),
      debug_(iConfig.getUntrackedParameter<int>("debug", 0)),
      useStandaloneMuons_(true),//iConfig.getParameter<bool>("useStandaloneMuons")),
      useTrackerMuons_(false), //iConfig.getParameter<bool>("useTrackerMuons")),
      hasTracks_(!iConfig.getParameter<edm::InputTag>("tracks").label().empty()),
      tkCands_(hasTracks_ ? consumes<l1t::PFTrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))
                          : edm::EDGetTokenT<l1t::PFTrackCollection>()),
      trkPt_(iConfig.getParameter<double>("trkPtCut")),
      muCands_(consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      //tkMuCands_(consumes<l1t::TkMuonCollection>(iConfig.getParameter<edm::InputTag>("tkMuons"))),
      emPtCut_(iConfig.getParameter<double>("emPtCut")),
      hadPtCut_(iConfig.getParameter<double>("hadPtCut")),
      sortOutputs_(iConfig.getParameter<bool>("sortOutputs")),
      regionizer_(nullptr),
      l1pfalgo_(nullptr),
      l1pualgo_(nullptr),
      l1tkegalgo_(nullptr),
      writeEgSta_(iConfig.getParameter<bool>("writeEGSta")),
      regionDumpName_(iConfig.getUntrackedParameter<std::string>("dumpFileName", "")),
      debugEta_(iConfig.getUntrackedParameter<double>("debugEta", 0)),
      debugPhi_(iConfig.getUntrackedParameter<double>("debugPhi", 0)),
      debugR_(iConfig.getUntrackedParameter<double>("debugR", -1)) {
  produces<l1t::PFCandidateCollection>("PF");
  produces<l1t::PFCandidateCollection>("Puppi");

  produces<l1t::PFCandidateCollection>("EmCalo");
  produces<l1t::PFCandidateCollection>("Calo");
  produces<l1t::PFCandidateCollection>("TK");
#if 0 // LATER
  produces<l1t::PFCandidateCollection>("TKVtx");
#endif

  if (writeEgSta_)
    produces<BXVector<l1t::EGamma>>("L1Eg");
  produces<l1t::TkElectronCollection>("L1TkEle");
  produces<l1t::TkEmCollection>("L1TkEm");


  for (const auto& tag : iConfig.getParameter<std::vector<edm::InputTag>>("emClusters")) {
    emCands_.push_back(consumes<l1t::PFClusterCollection>(tag));
  }
  for (const auto& tag : iConfig.getParameter<std::vector<edm::InputTag>>("hadClusters")) {
    hadCands_.push_back(consumes<l1t::PFClusterCollection>(tag));
  }


  const std::string& regalgo = iConfig.getParameter<std::string>("regionizerAlgo");
  if (regalgo == "Ideal") {
    regionizer_ = std::make_unique<l1ct::RegionizerEmulator>();
  } else
    throw cms::Exception("Configuration", "Unsupported regionizerAlgo");

  const std::string& algo = iConfig.getParameter<std::string>("pfAlgo");
  if (algo == "PFAlgo3") {
    l1pfalgo_ = std::make_unique<l1ct::PFAlgo3Emulator>(iConfig.getParameter<edm::ParameterSet>("pfAlgoParameters"));
  } else if (algo == "PFAlgo2HGC") {
    l1pfalgo_ = std::make_unique<l1ct::PFAlgo2HGCEmulator>(iConfig.getParameter<edm::ParameterSet>("pfAlgoParameters"));
  } else if (algo == "PFAlgoDummy") {
    l1pfalgo_ = std::make_unique<l1ct::PFAlgoDummyEmulator>(iConfig.getParameter<edm::ParameterSet>("pfAlgoParameters"));
  } else
    throw cms::Exception("Configuration", "Unsupported pfAlgo");

  const std::string& pualgo = iConfig.getParameter<std::string>("puAlgo");
  if (pualgo == "LinearizedPuppi") {
    l1pualgo_ = std::make_unique<l1ct::LinPuppiEmulator>(iConfig.getParameter<edm::ParameterSet>("puAlgoParameters"));
  } else
    throw cms::Exception("Configuration", "Unsupported puAlgo");


  l1tkegalgo_ = std::make_unique<l1ct::PFTkEGAlgoEmulator>(l1ct::pftkegalgo_config(iConfig.getParameter<edm::ParameterSet>("tkEgAlgoParameters")));

  extTkVtx_ = consumes<std::vector<l1t::TkPrimaryVertex>>(iConfig.getParameter<edm::InputTag>("vtxCollection"));

#if 0
  for (int tot = 0; tot <= 1; ++tot) {
    for (int i = 0; i < l1tpf_impl::Region::n_input_types; ++i) {
      produces<unsigned int>(std::string(tot ? "totNL1" : "maxNL1") + l1tpf_impl::Region::inputTypeName(i));
    }
    for (int i = 0; i < l1tpf_impl::Region::n_output_types; ++i) {
      produces<unsigned int>(std::string(tot ? "totNL1PF" : "maxNL1PF") + l1tpf_impl::Region::outputTypeName(i));
      produces<unsigned int>(std::string(tot ? "totNL1Puppi" : "maxNL1Puppi") + l1tpf_impl::Region::outputTypeName(i));
    }
  }
  for (int i = 0; i < l1tpf_impl::Region::n_input_types; ++i) {
    produces<std::vector<unsigned>>(std::string("vecNL1") + l1tpf_impl::Region::inputTypeName(i));
  }
  for (int i = 0; i < l1tpf_impl::Region::n_output_types; ++i) {
    produces<std::vector<unsigned>>(std::string("vecNL1PF") + l1tpf_impl::Region::outputTypeName(i));
    produces<std::vector<unsigned>>(std::string("vecNL1Puppi") + l1tpf_impl::Region::outputTypeName(i));
  }
#endif

  initSectorsAndRegions(iConfig);
}

L1TCorrelatorLayer1Producer::~L1TCorrelatorLayer1Producer() {
}

void L1TCorrelatorLayer1Producer::beginStream(edm::StreamID id) {
  if (!regionDumpName_.empty()) {
    if (id == 0) {
      fRegionDump_.open(regionDumpName_.c_str(), std::ios::out | std::ios::binary);
    } else {
      edm::LogWarning("L1TCorrelatorLayer1Producer")
          << "Job running with multiple streams, but dump file will have only events on stream zero.";
    }
  }
}

// ------------ method called to produce the data  ------------
void L1TCorrelatorLayer1Producer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // clear the regions also at the beginning, in case one event didn't complete but the job continues on
  initEvent(iEvent);

  /// ------ READ TRACKS ----
  if (hasTracks_) {
    edm::Handle<l1t::PFTrackCollection> htracks;
    iEvent.getByToken(tkCands_, htracks);
    const auto& tracks = *htracks;
    for (unsigned int itk = 0, ntk = tracks.size(); itk < ntk; ++itk) {
      const auto& tk = tracks[itk];
      // adding objects to PF
      if (debugR_ > 0 && deltaR(tk.eta(), tk.phi(), debugEta_, debugPhi_) > debugR_)
        continue;
      if (tk.pt() > trkPt_ && tk.quality() > 0) {
        addTrack(tk, l1t::PFTrackRef(htracks, itk));
      }
    }
  }

  /// ------ READ MUONS ----
  /// ------- first check that not more than one version of muons (standaloneMu or trackerMu) is set to be used in l1pflow
  if (useStandaloneMuons_ && useTrackerMuons_) {
    throw cms::Exception(
        "Configuration",
        "setting useStandaloneMuons=True && useTrackerMuons=True is not to be done, as it would duplicate all muons\n");
  }

  if (useStandaloneMuons_) {
    edm::Handle<l1t::MuonBxCollection> muons;
    iEvent.getByToken(muCands_, muons);
    for (auto it = muons->begin(0), ed = muons->end(0); it != ed; ++it) {
      const l1t::Muon& mu = *it;
      if (debugR_ > 0 && deltaR(mu.eta(), mu.phi(), debugEta_, debugPhi_) > debugR_)
        continue;
      addMuon(mu, l1t::PFCandidate::MuonRef(muons, muons->key(it)));
    }
  }

  if (useTrackerMuons_) {
      throw cms::Exception("Configuration", "Unsupported for now");
  }

  // ------ READ CALOS -----
  edm::Handle<l1t::PFClusterCollection> caloHandle;
  for (const auto& tag : emCands_) {
    iEvent.getByToken(tag, caloHandle);
    const auto& calos = *caloHandle;
    for (unsigned int ic = 0, nc = calos.size(); ic < nc; ++ic) {
      const auto& calo = calos[ic];
      if (debugR_ > 0 && deltaR(calo.eta(), calo.phi(), debugEta_, debugPhi_) > debugR_)
        continue;
      if (calo.pt() > emPtCut_)
        addEmCalo(calo, l1t::PFClusterRef(caloHandle, ic));
    }
  }
  for (const auto& tag : hadCands_) {
    iEvent.getByToken(tag, caloHandle);
    const auto& calos = *caloHandle;
    for (unsigned int ic = 0, nc = calos.size(); ic < nc; ++ic) {
      const auto& calo = calos[ic];
      if (debugR_ > 0 && deltaR(calo.eta(), calo.phi(), debugEta_, debugPhi_) > debugR_)
        continue;
      if (calo.pt() > hadPtCut_)
        addHadCalo(calo, l1t::PFClusterRef(caloHandle, ic));
    }
  }

  regionizer_->run(event_.decoded, event_.pfinputs);

  // First, get a copy of the discretized and corrected inputs, and write them out
  iEvent.put(fetchEmCalo(), "EmCalo");
  iEvent.put(fetchHadCalo(), "Calo");
  iEvent.put(fetchTracks(), "TK");

  // Then do the vertexing, and save it out
  float z0 = 0; double ptsum = 0;
  edm::Handle<std::vector<l1t::TkPrimaryVertex>> vtxHandle;
  iEvent.getByToken(extTkVtx_, vtxHandle);
  for (const l1t::TkPrimaryVertex& vtx : *vtxHandle) {
      if (ptsum == 0 || vtx.sum() > ptsum) {
          z0 = vtx.zvertex();
          ptsum = vtx.sum();
          l1ct::PVObjEmu hwpv;
          hwpv.hwZ0 = l1ct::Scales::makeZ0(z0);
          if (event_.pvs.empty()) { event_.pvs.push_back(hwpv); }
          else { event_.pvs[0] = hwpv; }
      }
  }

  // Then also save the tracks with a vertex cut
#if 0
  iEvent.put(l1regions_.fetchTracks(/*ptmin=*/0.0, /*fromPV=*/true), "TKVtx");
#endif

  // Then run PF in each region
  event_.out.resize(event_.pfinputs.size());
  for (unsigned int ir = 0, nr = event_.pfinputs.size(); ir < nr; ++ir) {
    l1pfalgo_->run(event_.pfinputs[ir], event_.out[ir]);
    l1pfalgo_->mergeNeutrals(event_.out[ir]);
    l1tkegalgo_->run(event_.pfinputs[ir], event_.out[ir]);
    //l1tkegalgo_->runTkIso(l1region, z0);
    //l1tkegalgo_->runPFIso(l1region, z0);
  }

  // Then run puppi (regionally)
  for (unsigned int ir = 0, nr = event_.pfinputs.size(); ir < nr; ++ir) {
    l1pualgo_->run(event_.pfinputs[ir], event_.pvs, event_.out[ir]);
    //l1pualgo_->runNeutralsPU(l1region, z0, -1., puGlobals);
    //l1region.outputCrop(sortOutputs_);
  }

  // save PF into the event
  iEvent.put(fetchPF(), "PF");

  // and save puppi
  iEvent.put(fetchPuppi(), "Puppi");

  // save the EG objects
  putEgObjects(iEvent, writeEgSta_, "L1Eg", "L1TkEm", "L1TkEle");
#if 0

  // Then go do the multiplicities

  for (int i = 0; i < l1tpf_impl::Region::n_input_types; ++i) {
    auto totAndMax = l1regions_.totAndMaxInput(i);
    addUInt(totAndMax.first, std::string("totNL1") + l1tpf_impl::Region::inputTypeName(i), iEvent);
    addUInt(totAndMax.second, std::string("maxNL1") + l1tpf_impl::Region::inputTypeName(i), iEvent);
    iEvent.put(l1regions_.vecInput(i), std::string("vecNL1") + l1tpf_impl::Region::inputTypeName(i));
  }
  for (int i = 0; i < l1tpf_impl::Region::n_output_types; ++i) {
    auto totAndMaxPF = l1regions_.totAndMaxOutput(i, false);
    auto totAndMaxPuppi = l1regions_.totAndMaxOutput(i, true);
    addUInt(totAndMaxPF.first, std::string("totNL1PF") + l1tpf_impl::Region::outputTypeName(i), iEvent);
    addUInt(totAndMaxPF.second, std::string("maxNL1PF") + l1tpf_impl::Region::outputTypeName(i), iEvent);
    addUInt(totAndMaxPuppi.first, std::string("totNL1Puppi") + l1tpf_impl::Region::outputTypeName(i), iEvent);
    addUInt(totAndMaxPuppi.second, std::string("maxNL1Puppi") + l1tpf_impl::Region::outputTypeName(i), iEvent);
    iEvent.put(l1regions_.vecOutput(i, false), std::string("vecNL1PF") + l1tpf_impl::Region::outputTypeName(i));
    iEvent.put(l1regions_.vecOutput(i, true), std::string("vecNL1Puppi") + l1tpf_impl::Region::outputTypeName(i));
  }

#endif
  if (fRegionDump_.is_open()) {
    event_.write(fRegionDump_);
  }


  // finally clear the regions
  event_.clear();
}

void L1TCorrelatorLayer1Producer::addUInt(unsigned int value, std::string iLabel, edm::Event& iEvent) {
  iEvent.put(std::make_unique<unsigned>(value), iLabel);
}

void L1TCorrelatorLayer1Producer::initSectorsAndRegions(const edm::ParameterSet & iConfig) {
    // the track finder geometry is fixed
    unsigned int TF_phiSlices = 9;
    float TF_phiWidth = 2 * M_PI / TF_phiSlices;
    event_.decoded.track.clear();
    for (unsigned int ieta = 0, neta = 2; ieta < neta; ++ieta) {
          for (unsigned int iphi = 0; iphi < TF_phiSlices; ++iphi) {
              float phiCenter = reco::reduceRange(iphi * TF_phiWidth); 
              event_.decoded.track.emplace_back(
                          (ieta ? 0. : -2.5), 
                          (ieta ? 2.5 : 0.0),
                          phiCenter,
                          TF_phiWidth);
          }
    }
 
    event_.decoded.emcalo.clear();
    event_.decoded.hadcalo.clear();
    for (const edm::ParameterSet &preg : iConfig.getParameter<std::vector<edm::ParameterSet>>("caloSectors")) {
      std::vector<double> etaBoundaries = preg.getParameter<std::vector<double>>("etaBoundaries");
      if (!std::is_sorted(etaBoundaries.begin(), etaBoundaries.end())) throw cms::Exception("Configuration", "caloSectors.etaBoundaries not sorted\n");
      unsigned int phiSlices = preg.getParameter<uint32_t>("phiSlices");
      float phiWidth = 2 * M_PI / phiSlices;
      if (phiWidth > 2*l1ct::Scales::maxAbsPhi()) throw cms::Exception("Configuration", "caloSectors phi range too large for phi_t data type");
      for (unsigned int ieta = 0, neta = etaBoundaries.size() - 1; ieta < neta; ++ieta) {
          float etaWidth = etaBoundaries[ieta + 1]-etaBoundaries[ieta];
          if (etaWidth > 2*l1ct::Scales::maxAbsEta()) throw cms::Exception("Configuration", "caloSectors eta range too large for eta_t data type");
          for (unsigned int iphi = 0; iphi < phiSlices; ++iphi) {
              float phiCenter = reco::reduceRange(iphi * phiWidth);  //align with L1 TrackFinder phi sector indexing for now
              event_.decoded.hadcalo.emplace_back(etaBoundaries[ieta],
                          etaBoundaries[ieta + 1],
                          phiCenter,
                          phiWidth);
              event_.decoded.emcalo.emplace_back(etaBoundaries[ieta],
                          etaBoundaries[ieta + 1],
                          phiCenter,
                          phiWidth);
        }
      }
    }
    
    event_.decoded.muon.region = l1ct::PFRegionEmu(
            -l1ct::Scales::maxAbsGlbEta(), l1ct::Scales::maxAbsGlbEta(),
            0.f, 2*l1ct::Scales::maxAbsGlbPhi(), 0.f, 0.f);

    event_.pfinputs.clear();
    for (const edm::ParameterSet &preg : iConfig.getParameter<std::vector<edm::ParameterSet>>("regions")) {
      std::vector<double> etaBoundaries = preg.getParameter<std::vector<double>>("etaBoundaries");
      if (!std::is_sorted(etaBoundaries.begin(), etaBoundaries.end())) throw cms::Exception("Configuration", "regions.etaBoundaries not sorted\n");
      unsigned int phiSlices = preg.getParameter<uint32_t>("phiSlices");
      float etaExtra = preg.getParameter<double>("etaExtra");
      float phiExtra = preg.getParameter<double>("phiExtra");
      float phiWidth = 2 * M_PI / phiSlices;
      for (unsigned int ieta = 0, neta = etaBoundaries.size() - 1; ieta < neta; ++ieta) {
          for (unsigned int iphi = 0; iphi < phiSlices; ++iphi) {
              float phiCenter = reco::reduceRange(iphi * phiWidth);  //align with L1 TrackFinder phi sector indexing
              event_.pfinputs.emplace_back(etaBoundaries[ieta],
                          etaBoundaries[ieta + 1],
                          phiCenter,
                          phiWidth,
                          etaExtra,
                          phiExtra);
        }
      }
    }

}

void L1TCorrelatorLayer1Producer::initEvent(const edm::Event &iEvent) {
  event_.clear();
  event_.run = iEvent.id().run(); 
  event_.lumi = iEvent.id().luminosityBlock();
  event_.event = iEvent.id().event();
  clusterRefMap_.clear();
  trackRefMap_.clear();
  muonRefMap_.clear();
}

void L1TCorrelatorLayer1Producer::addTrack(const l1t::PFTrack &t, l1t::PFTrackRef ref) {
    auto & sectors = event_.decoded.track;
    assert(sectors.size() == 18);
    int isec = t.track()->phiSector() + (t.eta() >= 0 ? 9 : 0);
    addDecodedTrack(sectors[isec], t);
    trackRefMap_[&t] = ref;
}
void L1TCorrelatorLayer1Producer::addMuon(const l1t::Muon &mu, l1t::PFCandidate::MuonRef ref) {
    addDecodedMuon(event_.decoded.muon, mu);
    muonRefMap_[&mu] = ref;
}
void L1TCorrelatorLayer1Producer::addHadCalo(const l1t::PFCluster &c, l1t::PFClusterRef ref) {
    for (auto & sec : event_.decoded.hadcalo) {
        if (sec.region.contains(c.eta(), c.phi())) {
            addDecodedHadCalo(sec, c);
        }
    }
    clusterRefMap_[&c] = ref;
}
void L1TCorrelatorLayer1Producer::addEmCalo(const l1t::PFCluster &c, l1t::PFClusterRef ref) {
    for (auto & sec : event_.decoded.emcalo) {
        if (sec.region.contains(c.eta(), c.phi())) {
            addDecodedEmCalo(sec, c);
        }
    }
  clusterRefMap_[&c] = ref;
}

void L1TCorrelatorLayer1Producer::addDecodedTrack(l1ct::DetectorSector<l1ct::TkObjEmu> & sec, const l1t::PFTrack &t) {
    l1ct::TkObjEmu tk;
    tk.hwPt = l1ct::Scales::makePtFromFloat(t.pt());
    tk.hwEta = l1ct::Scales::makeEta(sec.region.localEta(t.caloEta()));
    tk.hwPhi = l1ct::Scales::makePhi(sec.region.localPhi(t.caloPhi()));
    tk.hwCharge = t.charge() > 0;
    tk.hwQuality = t.quality();
    tk.hwDEta = l1ct::Scales::makeEta(t.eta() - t.caloEta());
    tk.hwDPhi = l1ct::Scales::makePhi(std::abs(reco::deltaPhi(t.phi(), t.caloPhi())));
    tk.hwZ0 =   l1ct::Scales::makeZ0(t.vertex().Z());
    tk.hwDxy  = 0;
    tk.hwChi2 = round(t.chi2() * 10);
    tk.hwStubs = t.nStubs();
    tk.src = & t;
    sec.obj.push_back(tk);
}


void L1TCorrelatorLayer1Producer::addDecodedMuon(l1ct::DetectorSector<l1ct::MuObjEmu> & sec, const l1t::Muon &t) {
    l1ct::MuObjEmu mu;
    mu.hwPt = l1ct::Scales::makePtFromFloat(t.pt());
    mu.hwEta = l1ct::Scales::makeEta(t.eta());
    mu.hwPhi = l1ct::Scales::makePhi(t.phi());
    mu.hwCharge = t.charge() > 0;
    mu.hwQuality = t.hwQual();
    mu.hwDEta = 0;
    mu.hwDPhi = 0;
    mu.hwZ0 =   l1ct::Scales::makeZ0(t.vertex().Z());
    mu.hwDxy  = 0;
    mu.src = & t;
    sec.obj.push_back(mu);
}

void L1TCorrelatorLayer1Producer::addDecodedHadCalo(l1ct::DetectorSector<l1ct::HadCaloObjEmu> & sec, const l1t::PFCluster &c) {
    l1ct::HadCaloObjEmu calo;
    calo.hwPt = l1ct::Scales::makePtFromFloat(c.pt());
    calo.hwEta = l1ct::Scales::makeEta(sec.region.localEta(c.eta()));
    calo.hwPhi = l1ct::Scales::makePhi(sec.region.localPhi(c.phi()));
    calo.hwEmPt = l1ct::Scales::makePtFromFloat(c.emEt());
    calo.hwIsEM = c.isEM();
    calo.src = & c;
    sec.obj.push_back(calo);
}

void L1TCorrelatorLayer1Producer::addDecodedEmCalo(l1ct::DetectorSector<l1ct::EmCaloObjEmu> & sec, const l1t::PFCluster &c) {
    l1ct::EmCaloObjEmu calo;
    calo.hwPt = l1ct::Scales::makePtFromFloat(c.pt());
    calo.hwEta = l1ct::Scales::makeEta(sec.region.localEta(c.eta()));
    calo.hwPhi = l1ct::Scales::makePhi(sec.region.localPhi(c.phi()));
    calo.hwPtErr = l1ct::Scales::makePtFromFloat(c.ptError());
    calo.hwFlags = c.hwQual();
    calo.src = & c;
    sec.obj.push_back(calo);
}


template<typename T>
void L1TCorrelatorLayer1Producer::setRefs_(l1t::PFCandidate &pf, const T & p) const {
    if (p.srcCluster) {
        auto match = clusterRefMap_.find(p.srcCluster);
        if (match == clusterRefMap_.end()) {
            throw cms::Exception("CorruptData") << "Invalid cluster pointer in PF candidate id " << p.intId() << " pt "
                << p.floatPt() << " eta " << p.floatEta() << " phi " << p.floatPhi();
        }
        pf.setPFCluster(match->second);
    }
    if (p.srcTrack) {
        auto match = trackRefMap_.find(p.srcTrack);
        if (match == trackRefMap_.end()) {
            throw cms::Exception("CorruptData") << "Invalid track pointer in PF candidate id " << p.intId() << " pt "
                << p.floatPt() << " eta " << p.floatEta() << " phi " << p.floatPhi();
        }
        pf.setPFTrack(match->second);
    }
    if (p.srcMu) {
        auto match = muonRefMap_.find(p.srcMu);
        if (match == muonRefMap_.end()) {
            throw cms::Exception("CorruptData") << "Invalid muon pointer in PF candidate id " << p.intId() << " pt "
                << p.floatPt() << " eta " << p.floatEta() << " phi " << p.floatPhi();
        }
        pf.setMuon(match->second);
    }
}

template<>
void L1TCorrelatorLayer1Producer::setRefs_<l1ct::PFNeutralObjEmu>(l1t::PFCandidate &pf, const l1ct::PFNeutralObjEmu & p) const {
    if (p.srcCluster) {
        auto match = clusterRefMap_.find(p.srcCluster);
        if (match == clusterRefMap_.end()) {
            throw cms::Exception("CorruptData") << "Invalid cluster pointer in PF candidate id " << p.intId() << " pt "
                << p.floatPt() << " eta " << p.floatEta() << " phi " << p.floatPhi();
        }
        pf.setPFCluster(match->second);
    }
}

template<>
void L1TCorrelatorLayer1Producer::setRefs_<l1ct::HadCaloObjEmu>(l1t::PFCandidate &pf, const l1ct::HadCaloObjEmu & p) const {
    if (p.src) {
        auto match = clusterRefMap_.find(p.src);
        if (match == clusterRefMap_.end()) {
            throw cms::Exception("CorruptData") << "Invalid cluster pointer in hadcalo candidate  pt "
                << p.floatPt() << " eta " << p.floatEta() << " phi " << p.floatPhi();
        }
        pf.setPFCluster(match->second);
    }
}

template<>
void L1TCorrelatorLayer1Producer::setRefs_<l1ct::EmCaloObjEmu>(l1t::PFCandidate &pf, const l1ct::EmCaloObjEmu & p) const {
    if (p.src) {
        auto match = clusterRefMap_.find(p.src);
        if (match == clusterRefMap_.end()) {
            throw cms::Exception("CorruptData") << "Invalid cluster pointer in emcalo candidate  pt "
                << p.floatPt() << " eta " << p.floatEta() << " phi " << p.floatPhi();
        }
        pf.setPFCluster(match->second);
    }
}

template<>
void L1TCorrelatorLayer1Producer::setRefs_<l1ct::TkObjEmu>(l1t::PFCandidate &pf, const l1ct::TkObjEmu & p) const {
    if (p.src) {
        auto match = trackRefMap_.find(p.src);
        if (match == trackRefMap_.end()) {
            throw cms::Exception("CorruptData") << "Invalid track pointer in track candidate  pt "
                << p.floatPt() << " eta " << p.floatEta() << " phi " << p.floatPhi();
        }
        pf.setPFTrack(match->second);
    }
}



std::unique_ptr<l1t::PFCandidateCollection> L1TCorrelatorLayer1Producer::fetchHadCalo() const {
  auto ret = std::make_unique<l1t::PFCandidateCollection>();
  for (const auto r : event_.pfinputs) {
      const auto & reg = r.region;
      for (const auto &p : r.hadcalo) {
          if (p.hwPt == 0 || !reg.isFiducial(p)) continue;
          reco::Particle::PolarLorentzVector p4(
                  p.floatPt(), reg.floatGlbEtaOf(p), reg.floatGlbPhiOf(p), 0.13f);
          l1t::PFCandidate::ParticleType type = p.hwIsEM ? l1t::PFCandidate::Photon : l1t::PFCandidate::NeutralHadron;
          ret->emplace_back(type, 0, p4, 1, p.intPt(), p.intEta(), p.intPhi());
          setRefs_(ret->back(), p);
      }
  }
  return ret;
}
std::unique_ptr<l1t::PFCandidateCollection> L1TCorrelatorLayer1Producer::fetchEmCalo() const {
  auto ret = std::make_unique<l1t::PFCandidateCollection>();
  for (const auto r : event_.pfinputs) {
      const auto & reg = r.region;
      for (const auto &p : r.emcalo) {
          if (p.hwPt == 0 || !reg.isFiducial(p)) continue;
          reco::Particle::PolarLorentzVector p4(
                  p.floatPt(), reg.floatGlbEtaOf(p), reg.floatGlbPhiOf(p), 0.13f);
          ret->emplace_back(l1t::PFCandidate::Photon, 0, p4, 1, p.intPt(), p.intEta(), p.intPhi());
          setRefs_(ret->back(), p);
      }
  }
  return ret;
}
std::unique_ptr<l1t::PFCandidateCollection> L1TCorrelatorLayer1Producer::fetchTracks() const {
  auto ret = std::make_unique<l1t::PFCandidateCollection>();
  for (const auto r : event_.pfinputs) {
      const auto & reg = r.region;
      for (const auto &p : r.track) {
          if (p.hwPt == 0 || !reg.isFiducial(p)) continue;
          reco::Particle::PolarLorentzVector p4(
                  p.floatPt(), reg.floatGlbEta(p.hwVtxEta()), reg.floatGlbPhi(p.hwVtxPhi()), 0.13f);
          ret->emplace_back(l1t::PFCandidate::ChargedHadron, p.intCharge(), p4, 1, p.intPt(), p.intEta(), p.intPhi());
          setRefs_(ret->back(), p);
      }
  }
  return ret;
}



std::unique_ptr<l1t::PFCandidateCollection> L1TCorrelatorLayer1Producer::fetchPF() const {
  auto ret = std::make_unique<l1t::PFCandidateCollection>();
  for (unsigned int ir = 0, nr = event_.pfinputs.size(); ir < nr; ++ir) {
    const auto & reg = event_.pfinputs[ir].region;
      for (const auto &p : event_.out[ir].pfcharged) {
          if (p.hwPt == 0 || !reg.isFiducial(p)) continue;
          reco::Particle::PolarLorentzVector p4(
                  p.floatPt(), reg.floatGlbEta(p.hwVtxEta()), reg.floatGlbPhi(p.hwVtxPhi()), 0.13f);
          l1t::PFCandidate::ParticleType type = l1t::PFCandidate::ChargedHadron;
          if (p.hwId.isMuon()) type = l1t::PFCandidate::Muon;
          else if (p.hwId.isElectron()) type = l1t::PFCandidate::Electron;
          ret->emplace_back(type, p.intCharge(), p4, 1, p.intPt(), p.intEta(), p.intPhi());
          ret->back().setVertex(reco::Particle::Point(0, 0, p.floatZ0()));
          setRefs_(ret->back(), p);
      }
      for (const auto &p : event_.out[ir].pfneutral) {
          if (p.hwPt == 0 || !reg.isFiducial(p)) continue;
          reco::Particle::PolarLorentzVector p4(
                  p.floatPt(), reg.floatGlbEtaOf(p), reg.floatGlbPhiOf(p), 0.13f);
          l1t::PFCandidate::ParticleType type = p.hwId.isPhoton() ? l1t::PFCandidate::Photon : l1t::PFCandidate::NeutralHadron;
          ret->emplace_back(type, 0, p4, 1, p.intPt(), p.intEta(), p.intPhi());
          setRefs_(ret->back(), p);
      }
  }
  return ret;
}

std::unique_ptr<l1t::PFCandidateCollection> L1TCorrelatorLayer1Producer::fetchPuppi() const {
    auto ret = std::make_unique<l1t::PFCandidateCollection>();
    for (unsigned int ir = 0, nr = event_.pfinputs.size(); ir < nr; ++ir) {
        for (const auto &p : event_.out[ir].puppi) {
            if (p.hwPt == 0) continue;
            // note: Puppi candidates are already in global coordinates & fiducial-only!
            l1t::PFCandidate::ParticleType type; float mass = 0.13f;
            if (p.hwId.charged()) {
                if (p.hwId.isMuon()) { type = l1t::PFCandidate::Muon; mass = 0.105; }
                else if (p.hwId.isElectron()) { type = l1t::PFCandidate::Electron; mass = 0.005; }
                else type = l1t::PFCandidate::ChargedHadron;
            } else {
                type = p.hwId.isPhoton() ? l1t::PFCandidate::Photon : l1t::PFCandidate::NeutralHadron;
                mass = p.hwId.isPhoton() ? 0.0 : 0.5;
            }
            reco::Particle::PolarLorentzVector p4(p.floatPt(), p.floatEta(), p.floatPhi(), mass);
            ret->emplace_back(type, p.intCharge(), p4, p.floatPuppiW(), p.intPt(), p.intEta(), p.intPhi());
            if (p.hwId.charged()) {
                ret->back().setVertex(reco::Particle::Point(0, 0, p.floatZ0()));
            }
        }
    }
    return ret;
}

void L1TCorrelatorLayer1Producer::putEgObjects(edm::Event &iEvent,
                                               const bool writeEgSta,
                                               const std::string &egLablel,
                                               const std::string &tkEmLabel,
                                               const std::string &tkEleLabel) const {
  auto egs = std::make_unique<BXVector<l1t::EGamma>>();
  auto tkems = std::make_unique<l1t::TkEmCollection>();
  auto tkeles = std::make_unique<l1t::TkElectronCollection>();

  // edm::RefProd<BXVector<l1t::EGamma>> ref_egs;
  // if (writeEgSta)
  //   ref_egs = iEvent.getRefBeforePut<BXVector<l1t::EGamma>>(egLablel);
  // 
  // edm::Ref<BXVector<l1t::EGamma>>::key_type idx = 0;
  // 
  // for (const Region &r : regions_) {
  //   for (const auto &egphoton : r.egphotons) {
  //     if (egphoton.floatPt() < ptMin)
  //       continue;
  // 
  //     if (!r.fiducialLocal(egphoton.floatEta(), egphoton.floatPhi()))
  //       continue;
  // 
  //     edm::Ref<BXVector<l1t::EGamma>> reg;
  //     auto mom = reco::Candidate::PolarLorentzVector(
  //         egphoton.floatPt(), r.globalEta(egphoton.floatEta()), r.globalPhi(egphoton.floatPhi()), 0.);
  //     if (writeEgSta) {
  //       l1t::EGamma eg(mom);
  //       eg.setHwQual(egphoton.hwQual);
  //       egs->push_back(0, eg);
  //       reg = edm::Ref<BXVector<l1t::EGamma>>(ref_egs, idx++);
  //     } else {
  //       auto egptr = egphoton.cluster.src->constituentsAndFractions()[0].first;
  //       reg = edm::Ref<BXVector<l1t::EGamma>>(egptr.id(), dynamic_cast<const l1t::EGamma *>(egptr.get()), egptr.key());
  //     }
  // 
  //     l1t::TkEm tkem(reco::Candidate::LorentzVector(mom), reg, egphoton.floatIso(), egphoton.floatIsoPV());
  //     tkem.setHwQual(egphoton.hwQual);
  //     tkem.setPFIsol(egphoton.floatPFIso());
  //     tkem.setPFIsolPV(egphoton.floatPFIsoPV());
  //     tkems->push_back(tkem);
  // 
  //     if (egphoton.ele_idx == -1)
  //       continue;
  // 
  //     const auto &egele = r.egeles[egphoton.ele_idx];
  // 
  //     if (!r.fiducialLocal(egele.floatEta(), egele.floatPhi()))
  //       continue;
  // 
  //     auto mom_ele = reco::Candidate::PolarLorentzVector(
  //         egele.floatPt(), r.globalEta(egele.floatEta()), r.globalPhi(egele.floatPhi()), 0.);
  // 
  //     l1t::TkElectron tkele(
  //         reco::Candidate::LorentzVector(mom_ele), reg, edm::refToPtr(egele.track.src->track()), egele.floatIso());
  //     tkele.setHwQual(egele.hwQual);
  //     tkele.setPFIsol(egele.floatPFIso());
  //     tkeles->push_back(tkele);
  //   }
  // }
  if (writeEgSta)
    iEvent.put(std::move(egs), egLablel);
  iEvent.put(std::move(tkems), tkEmLabel);
  iEvent.put(std::move(tkeles), tkEleLabel);
}



//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TCorrelatorLayer1Producer);

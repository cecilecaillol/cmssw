#include "L1Trigger/Phase2L1ParticleFlow/interface/egamma/l2egsorter_ref.h"

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <memory>
#include <iostream>

using namespace l1ct;

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"

l1ct::L2EgSorterEmulator::L2EgSorterEmulator(const edm::ParameterSet &pset)
    : L2EgSorterEmulator(pset.getParameter<uint32_t>("nBOARDS"),
                         pset.getParameter<uint32_t>("nEGPerBoard"),
                         pset.getParameter<uint32_t>("nEGOut")) {}

#endif

void L2EgSorterEmulator::toFirmware(const std::vector<EGIsoObjEmu> &out_photons,
                                    const std::vector<EGIsoEleObjEmu> &out_eles,
                                    EGIsoObj out_egphs[/*nObjOut*/],
                                    EGIsoEleObj out_egeles[/*nObjOut*/]) const {
  for (unsigned int io = 0; io < nEGOut; io++) {
    EGIsoObj pho;
    EGIsoEleObj ele;
    if (io < out_photons.size())
      pho = out_photons[io];
    else
      pho.clear();
    if (io < out_eles.size())
      ele = out_eles[io];
    else
      ele.clear();

    out_egphs[io] = pho;
    out_egeles[io] = ele;
  }
}

void L2EgSorterEmulator::run(const std::vector<l1ct::OutputBoard> &in,
                             std::vector<EGIsoObjEmu> &out_photons,
                             std::vector<EGIsoEleObjEmu> &out_eles) const {
  // we copy to be able to resize them
  std::vector<std::vector<EGIsoObjEmu>> photons_in;
  std::vector<std::vector<EGIsoEleObjEmu>> eles_in;
  photons_in.reserve(in.size());
  eles_in.reserve(in.size());
  for (const auto &board : in) {
    std::vector<EGIsoObjEmu> photons = board.egphoton;
    std::vector<EGIsoEleObjEmu> eles = board.egelectron;
    resize_input(photons);
    resize_input(eles);

    photons_in.push_back(photons);
    eles_in.push_back(eles);
  }
  merge(photons_in, out_photons);
  merge(eles_in, out_eles);
}

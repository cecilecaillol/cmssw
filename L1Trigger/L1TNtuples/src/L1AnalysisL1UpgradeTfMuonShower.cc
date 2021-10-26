#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeTfMuonShower.h"
#include <cmath>
L1Analysis::L1AnalysisL1UpgradeTfMuonShower::L1AnalysisL1UpgradeTfMuonShower() {}

L1Analysis::L1AnalysisL1UpgradeTfMuonShower::~L1AnalysisL1UpgradeTfMuonShower() {}

void L1Analysis::L1AnalysisL1UpgradeTfMuonShower::SetTfMuonShower(const l1t::RegionalMuonShowerBxCollection& muonShower,
                                                      unsigned maxL1UpgradeTfMuonShower) {
  for (int ibx = muonShower.getFirstBX(); ibx <= muonShower.getLastBX(); ++ibx) {
    for (auto it = muonShower.begin(ibx);
         it != muonShower.end(ibx) && l1upgradetfmuonshower_.nTfMuonShowers < maxL1UpgradeTfMuonShower;
         ++it) {
      if (it->isValid()) {
        l1upgradetfmuonshower_.tfMuonShowerBx.push_back(ibx);
        int type;
        if (it->isOneNominalInTime())
          type = L1AnalysisL1UpgradeTfMuonShowerDataFormat::kOneNominal;
        else if (it->isOneTightInTime())
          type = L1AnalysisL1UpgradeTfMuonShowerDataFormat::kOneTight;
        else if (it->isTwoLooseInTime())
          type = L1AnalysisL1UpgradeTfMuonShowerDataFormat::kTwoLoose;
        else
          type = L1AnalysisL1UpgradeTfMuonShowerDataFormat::kInvalid;
        l1upgradetfmuonshower_.tfMuonShowerType.push_back(type);
        l1upgradetfmuonshower_.nTfMuonShowers++;
      }
    }
  }
}

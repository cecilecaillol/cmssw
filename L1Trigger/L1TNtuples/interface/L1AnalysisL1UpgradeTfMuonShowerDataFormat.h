#ifndef __L1Analysis_L1AnalysisL1UpgradeTfMuonShowerDataFormat_H__
#define __L1Analysis_L1AnalysisL1UpgradeTfMuonShowerDataFormat_H__

#include <vector>
#include <map>

namespace L1Analysis {

  struct L1AnalysisL1UpgradeTfMuonShowerDataFormat {

    enum TfMuonShowerType {
      kInvalid,
      kOneNominal,
      kOneTight,
      kTwoLoose
    };

    L1AnalysisL1UpgradeTfMuonShowerDataFormat() { Reset(); };
    ~L1AnalysisL1UpgradeTfMuonShowerDataFormat(){};

    void Reset() {
      nTfMuonShowers = 0;
      tfMuonShowerBx.clear();
      tfMuonShowerType.clear();
    }

    unsigned short int nTfMuonShowers;
    std::vector<short int> tfMuonShowerBx;
    std::vector<short int> tfMuonShowerType;
  };
}  // namespace L1Analysis
#endif

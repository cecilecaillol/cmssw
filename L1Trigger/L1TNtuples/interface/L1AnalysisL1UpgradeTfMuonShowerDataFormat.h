#ifndef __L1Analysis_L1AnalysisL1UpgradeTfMuonShowerDataFormat_H__
#define __L1Analysis_L1AnalysisL1UpgradeTfMuonShowerDataFormat_H__

#include <vector>
#include <map>

namespace L1Analysis {
  struct L1AnalysisL1UpgradeTfMuonShowerDataFormat {
    L1AnalysisL1UpgradeTfMuonShowerDataFormat() { Reset(); };
    ~L1AnalysisL1UpgradeTfMuonShowerDataFormat(){};

    void Reset() {
      nTfMuonShowers = 0;
      tfMuonShowerOneNominal.clear();
      tfMuonShowerOneTight.clear();
      tfMuonShowerTwoLoose.clear();
    }

    unsigned short int nTfMuonShowers;
    std::vector<short int> tfMuonShowerOneNominal;
    std::vector<short int> tfMuonShowerOneTight;
    std::vector<short int> tfMuonShowerTwoLoose;
  };
}  // namespace L1Analysis
#endif

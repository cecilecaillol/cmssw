#ifndef UCTGeometry_hh
#define UCTGeometry_hh

// UCT Geometric mapping between EB+EE/HB+HE & HF TPGs to CPT7s
// Access methods check if the mapping is correct
// Access methods also allow conversion of TPG to CTP7 mapping
// This class only has functions and defintions, no local memory

/*
    
  Crate Eta mapping:

  eta:      -5 ... -4 ... -3 ... -2 ... -1 ... 0 ... 1 ... 2 ...  3 ...  4 ...  5
  absIEta:  41   ...  30(29)28     ...        1 1     ...       28(29)30   ... 41
  region:   12 11 10  9  8  7  6  5  4 3 2 1  0 0 1 2  3  4   5  6 7 8 9 10 11 12

  For EB/EE towers we have 8-bit ET + 1-bit finegrain input
  For HB/HE towers we have 8-bit ET + 5-bit feature (TBD)
  For HF only 8-bits ET per input tower

  Each EB/HB+EE/HE regions are 4x4 == ~0.35 x 0.35 iEta x iPhi wide

  Calo TPG eta index is 1 through 28 for EB/HB+EE/HE
  Note that tower 29 is unused in trigger as it shadows portion of HE
  Corresponding negative values are for negative eta

  The uHTR HF eta index is 30-41 for 12 divisions in eta between 3-5
  Note that eta=30-39 are 10-degree phi towers.  However, Layer2 wants the same
  5-degree phi towers as in 1-28.  Therefore, we make artificial splitting for
  output.  For eta=40,42 we take 20-degree phi tower and split four ways !

  In the HF region size is proposed as 2x4 HF towers in etaxphi == 0.33 x 0.35,
  where the phi division of x4 assumes above /2 and /4 splits.

*/

#include <iostream>

#include <utility>

namespace l1tcalo {
  constexpr uint32_t NCrates{3};
  constexpr uint32_t NCardsInCrate{6};
  constexpr uint32_t NRegionsInCard{7};
  constexpr uint32_t NEtaInRegion{4};
  constexpr uint32_t NPhiInRegion{4};
  constexpr uint32_t NPhiInCard{
      NPhiInRegion};  // For convenience of Layer-2 we always have the same number of phi divisions

  constexpr uint32_t HFEtaOffset{NRegionsInCard * NEtaInRegion + 1};
  constexpr uint32_t NHFRegionsInCard{6};
  constexpr uint32_t NHFEtaInRegion{2};

  constexpr uint32_t NSides{2};  // Positive and Negative Eta sides
  constexpr uint32_t NEta{NRegionsInCard * NEtaInRegion};
  constexpr uint32_t NPhi{NCrates * NCardsInCrate * NPhiInRegion};
  constexpr uint32_t NHFEta{NHFRegionsInCard * NHFEtaInRegion};
  constexpr uint32_t NHFPhi{NCrates * NCardsInCrate * NPhiInRegion};

  constexpr uint32_t MaxCrateNumber{(NCrates - 1)};
  constexpr uint32_t MaxCardNumber{(NCardsInCrate - 1)};
  constexpr uint32_t MaxRegionNumber{(NRegionsInCard - 1)};
  constexpr uint32_t MaxEtaInRegion{(NEtaInRegion - 1)};
  constexpr uint32_t MaxPhiInRegion{(NPhiInRegion - 1)};

  constexpr int MaxCaloEta{41};
  constexpr int MaxCaloPhi{72};
  constexpr uint32_t CaloHFRegionStart{7};
  constexpr uint32_t CaloVHFRegionStart{12};

  constexpr uint32_t MaxUCTRegionsPhi{MaxCaloPhi / NPhiInRegion};
  constexpr uint32_t MaxUCTRegionsEta{NRegionsInCard + NHFRegionsInCard}; // Labelled -MaxUCTRegionsEta to +MaxUCTRegionsEta skipping 0

  // Binning for Layer1 calibration LUTs
  const uint32_t nEtBins = 256;
  const uint32_t nCalEtaBins = 28;
  const uint32_t nCalSideBins = 2;
  const uint32_t nHfEtaBins = 12;

}  // namespace l1tcalo

typedef std::pair<int, uint32_t> UCTRegionIndex;
typedef std::pair<int, int> UCTTowerIndex;

class UCTGeometry {
public:
  UCTGeometry();
  ~UCTGeometry() = default;

  // Calorimeter indices are defined to be ints and do not count from zero
  // Eta index sign indicates negative or positive eta
  // The zero eta value is illegal
  // Phi indices go 1 through 72 and other values are illegal

  int getCaloEtaIndex(bool negativeSide, uint32_t region, uint32_t iEta);
  int getCaloPhiIndex(uint32_t crate, uint32_t card, uint32_t region, uint32_t iPhi);

  uint32_t getLinkNumber(bool negativeSide, uint32_t region, uint32_t iEta, uint32_t iPhi);
  uint32_t getChannelNumber(bool negativeSide, uint32_t iEta, uint32_t iPhi);

  uint32_t getNCrates() { return l1tcalo::NCrates; }
  uint32_t getNCards() { return l1tcalo::NCardsInCrate; }
  uint32_t getNRegions() { return (l1tcalo::NRegionsInCard + l1tcalo::NHFRegionsInCard); }
  uint32_t getNEta(uint32_t region);
  uint32_t getNPhi(uint32_t region);

  uint32_t getCrate(int caloEta, int caloPhi);
  uint32_t getCard(int caloEta, int caloPhi);
  uint32_t getRegion(int caloEta, int caloPhi);
  uint32_t getiEta(int caloEta);
  uint32_t getiPhi(int caloPhi);
  bool getNegativeSide(int caloEta) {return (caloEta < 0);}

  bool checkCrate(uint32_t crate) { return !(crate < l1tcalo::NCrates); }
  bool checkCard(uint32_t card) { return !(card < l1tcalo::NCardsInCrate); }
  bool checkRegion(uint32_t region) { return !(region < (l1tcalo::NRegionsInCard + l1tcalo::NHFRegionsInCard)); }
  bool checkEtaIndex(uint32_t region, uint32_t iEta) {
    if (region < l1tcalo::NRegionsInCard)
      return !(iEta < l1tcalo::NEtaInRegion);
    else
      return !(iEta < l1tcalo::NHFEtaInRegion);
  }
  bool checkPhiIndex(uint32_t region, uint32_t iPhi) { return !(iPhi < l1tcalo::NPhiInRegion); }

  // We label regions by phi and eta indices
  //  UCTRegionPhiIndices are 0-17
  //  UCTRegionEtaIndices are 1-7 for EB/HB+EE/HE and 8-13 for HF,
  //   with negative values for negative eta, and zero being illegal
  //  GCTRegionPhiIndices are the same
  //  GCTRegionEtaIndices are 4-10 for -7 thru -1 and 11-17 for 1 thru 7
  //  GCTRegionEtaIndices 0-3 are for -HF and 18-22 are for +HF, which are now unused!
  // We label by the pair (UCTRegionPhiIndex, UCTRegionEtaIndex)

  uint32_t getUCTRegionPhiIndex(int caloPhi) {
    if (caloPhi < 71)
      return ((caloPhi + 1) / 4);
    else if (caloPhi < 73)
      return 0;
    else
      return 0xDEADBEEF;
  }
  uint32_t getGCTRegionPhiIndex(int caloPhi) { return getUCTRegionPhiIndex(caloPhi); }
  uint32_t getGCTRegionPhiIndex(UCTRegionIndex r) { return r.second; }

  int getUCTRegionEtaIndex(int caloEta) {
    // Region index is same for all phi; so get for phi = 1
    uint32_t rgn = getRegion(caloEta, 1);
    if (caloEta < 0)
      return -(rgn + 1);
    return (rgn + 1);
  }
  uint32_t getGCTRegionEtaIndex(int caloEta) {
    uint32_t gctEta = 0xDEADBEEF;
    // Region index is same for all phi; so get for phi = 1
    uint32_t region = getRegion(caloEta, 1);
    if (caloEta > 0 && region < l1tcalo::NRegionsInCard)
      gctEta = region + 11;
    else if (caloEta < 0 && region < l1tcalo::NRegionsInCard)
      gctEta = 10 - region;
    return gctEta;
  }
  uint32_t getGCTRegionEtaIndex(UCTRegionIndex r) {
    uint32_t gctEta = 0xDEADBEEF;
    uint32_t region = std::abs(r.first) - 1;
    if (r.first > 0 && region < l1tcalo::NRegionsInCard)
      gctEta = region + 11;
    else if (r.first < 0 && region < l1tcalo::NRegionsInCard)
      gctEta = 10 - region;
    return gctEta;
  }

  uint32_t getUCTRegionPhiIndex(uint32_t crate, uint32_t card);

  int getUCTRegionEtaIndex(bool negativeSide, uint32_t region) {
    if (negativeSide)
      return -(region + 1);
    else
      return (region + 1);
  }

  UCTRegionIndex getUCTRegionIndex(UCTTowerIndex caloTower) {
    return getUCTRegionIndex(caloTower.first, caloTower.second);
  }

  UCTRegionIndex getUCTRegionIndex(int caloEta, int caloPhi);

  UCTRegionIndex getUCTRegionIndex(bool negativeSide, uint32_t crate, uint32_t card, uint32_t region);

  UCTTowerIndex getUCTTowerIndex(UCTRegionIndex r, uint32_t iEta = 0, uint32_t iPhi = 0);

  double getUCTTowerEta(int caloEta);
  double getUCTTowerPhi(int caloPhi);

  UCTRegionIndex getUCTRegionIndexFromL1CaloRegion(uint32_t caloRegionEta, uint32_t caloRegionPhi);

  UCTTowerIndex getUCTTowerIndexFromL1CaloRegion(UCTRegionIndex r, uint32_t rawData);

private:
};

#endif

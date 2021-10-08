#ifndef DataFormats_L1Trigger_MuonShower_h
#define DataFormats_L1Trigger_MuonShower_h

/*
  This class is derived from the L1Candidate primarily to interface easily
  with the Global Muon Trigger. In the trigger system the MuonShower object
  carries only up to 4 bits of information, 2 for in-time showers,
  2 for out-of-time showers.
*/

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/L1TObjComparison.h"

namespace l1t {

  class MuonShower;
  typedef BXVector<MuonShower> MuonShowerBxCollection;
  typedef edm::Ref<MuonShowerBxCollection> MuonShowerRef;
  typedef edm::RefVector<MuonShowerBxCollection> MuonShowerRefVector;
  typedef std::vector<MuonShowerRef> MuonShowerVectorRef;

  typedef ObjectRefBxCollection<MuonShower> MuonShowerRefBxCollection;
  typedef ObjectRefPair<MuonShower> MuonShowerRefPair;
  typedef ObjectRefPairBxCollection<MuonShower> MuonShowerRefPairBxCollection;

  class MuonShower : public L1Candidate {
  public:
    MuonShower(bool oneNominalInTime = false,
               bool oneNominalOutOfTime = false,
               bool twoLooseInTime = false,
               bool twoLooseOutOfTime = false,
               bool oneTightInTime = false,
               bool oneTightOutOfTime = false);

    ~MuonShower() override;

    /*
      In CMSSW we consider 3 valid cases:
      - 1 nominal shower (baseline trigger for physics at Run-3)
      - 2 loose showers (to extend the physics reach)
      - 1 tight shower (backup trigger)

      In the uGT and UTM library, the hadronic shower trigger data is split
      over 4 bits: 2 for in-time trigger data, 2 for out-of-time trigger data
      - mus0, mus1 for in-time
      - musOutOfTime0, musOutOfTime1 for out-of-time

      The mapping for Run-3 is as follows:
      - 1 nominal shower -> 0b01
      - 2 loose showers -> 0b10
      - 1 tight shower -> 0b11
      This is done separately for the in-time and out-of-time trigger data
    */

    void setMus0(const bool bit) {mus0_ = bit;}
    void setMus1(const bool bit) {mus1_ = bit;}
    void setMusOutOfTime0(const bool bit) {musOutOfTime0_ = bit;}
    void setMusOutOfTime1(const bool bit) {musOutOfTime1_ = bit;}

    bool mus0() const {return mus0_;}
    bool mus1() const {return mus1_;}
    bool musOutOfTime0() const {return musOutOfTime0_;}
    bool musOutOfTime1() const {return musOutOfTime1_;}

    // at least one bit must be valid
    bool isValid() const;

    // useful members for trigger performance studies
    bool isOneNominalInTime() const {return mus0_;}
    bool isOneNominalOutOfTime() const {return musOutOfTime0_;}
    bool isTwoLooseInTime() const {return mus1_;}
    bool isTwoLooseOutOfTime() const {return musOutOfTime1_;}
    bool isOneTightInTime() const {return mus0_ and mus1_;}
    bool isOneTightOutOfTime() const {return musOutOfTime0_ and musOutOfTime1_;}

    virtual bool operator==(const l1t::MuonShower& rhs) const;
    virtual inline bool operator!=(const l1t::MuonShower& rhs) const { return !(operator==(rhs)); };

  private:
    // Run-3 definitions as provided in DN-20-033
    // in time and out-of-time qualities. only 2 bits each.
    bool mus0_;
    bool mus1_;
    bool musOutOfTime0_;
    bool musOutOfTime1_;
  };

}  // namespace l1t

#endif

#ifndef L1_TRACK_TRIGGER_TRACK_WORD_H
#define L1_TRACK_TRIGGER_TRACK_WORD_H

////////
//
// class to store the 96-bit track word produced by the L1 Track Trigger.  Intended to be inherited by L1 TTTrack.
// packing scheme given below.
//
// author:      Mike Hildreth
// modified by: Alexx Perloff
// created:     April 9, 2019
// modified:    March 9, 2021
//
///////

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include <ap_int.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <string>
#include <vector>

class TTTrack_TrackWord {
public:
  // ----------constants, enums and typedefs ---------
  enum TrackBitWidths {
    // The sizes of the track word components
    kMVAOtherSize = 6,    // Space for two specialized MVA selections
    kMVAQualitySize = 3,  // Width of track quality MVA
    kHitPatternSize = 7,  // Width of the hit pattern for stubs
    kBendChi2Size = 3,    // Width of the Bend-Chi2
    kD0Size = 13,         // Width of D0
    kChi2RZSize = 4,      // Width of Chi2 for r-z
    kZ0Size = 12,         // Width of z-position (40cm / 0.1)
    kTanlSize = 16,       // Width of tan(lambda)
    kChi2RPhiSize = 4,    // Width of Chi2 for r-phi
    kPhiSize = 12,        // Width of phi
    kRinvSize = 15,       // Width of Rinv
    kValidSize = 1,       // Valid bit

    kTrackWordSize = kValidSize + kRinvSize + kPhiSize + kChi2RPhiSize + kTanlSize + kZ0Size + kChi2RZSize 
                   + kD0Size + kBendChi2Size + kHitPatternSize + kMVAQualitySize + kMVAOtherSize,  // Width of the track word in bits
  };

  enum TrackBitLocations {
    // The location of the least significant bit (LSB) and most significant bit (MSB) in the track word for different fields
    kMVAOtherLSB = 0,
    kMVAOtherMSB = kMVAOtherLSB + TrackBitWidths::kMVAOtherSize - 1,
    kMVAQualityLSB = kMVAOtherMSB + 1,
    kMVAQualityMSB = kMVAQualityLSB + TrackBitWidths::kMVAQualitySize - 1,
    kHitPatternLSB = kMVAQualityMSB + 1,
    kHitPatternMSB = kHitPatternLSB + TrackBitWidths::kHitPatternSize - 1,
    kBendChi2LSB = kHitPatternMSB + 1,
    kBendChi2MSB = kBendChi2LSB + TrackBitWidths::kBendChi2Size - 1,
    kD0LSB = kBendChi2MSB + 1,
    kD0MSB = kD0LSB + TrackBitWidths::kD0Size - 1,
    kChi2RZLSB = kD0MSB + 1,
    kChi2RZMSB = kChi2RZLSB + TrackBitWidths::kChi2RZSize - 1,
    kZ0LSB = kChi2RZMSB + 1,
    kZ0MSB = kZ0LSB + TrackBitWidths::kZ0Size - 1,
    kTanlLSB = kZ0MSB + 1,
    kTanlMSB = kTanlLSB + TrackBitWidths::kTanlSize - 1,
    kChi2RPhiLSB = kTanlMSB + 1,
    kChi2RPhiMSB = kChi2RPhiLSB + TrackBitWidths::kChi2RPhiSize - 1,
    kPhiLSB = kChi2RPhiMSB + 1,
    kPhiMSB = kPhiLSB + TrackBitWidths::kPhiSize - 1,
    kRinvLSB = kPhiMSB + 1,
    kRinvMSB = kRinvLSB + TrackBitWidths::kRinvSize - 1,
    kValidLSB = kRinvMSB + 1,
    kValidMSB = kValidLSB + TrackBitWidths::kValidSize - 1,
  };

  // Binning constants
  static constexpr double minRinv = -0.006;
  static constexpr double minPhi0 = -0.7853981696;  // relative to the center of the sector
  static constexpr double minTanl = -8.;
  static constexpr double minZ0 = -20.46912512;
  static constexpr double minD0 = -16.;

  static constexpr double stepRinv = (2. * std::abs(minRinv)) / (1 << TrackBitWidths::kRinvSize);
  static constexpr double stepPhi0 = (2. * std::abs(minPhi0)) / (1 << TrackBitWidths::kPhiSize);
  static constexpr double stepTanL = (1. / (1 << 12));
  static constexpr double stepZ0 = (2. * std::abs(minZ0)) / (1 << TrackBitWidths::kZ0Size);
  static constexpr double stepD0 = (1. / (1 << 8));

  static constexpr std::array<double, 1 << TrackBitWidths::kChi2RPhiSize> chi2RPhiBins = {
      {0., 0.25, 0.5, 1., 2., 3., 5., 7., 10., 20., 40., 100., 200., 500., 1000., 3000.}};
  static constexpr std::array<double, 1 << TrackBitWidths::kChi2RZSize> chi2RZBins = {
      {0., 0.25, 0.5, 1., 2., 3., 5., 7., 10., 20., 40., 100., 200., 500., 1000., 3000.}};
  static constexpr std::array<double, 1 << TrackBitWidths::kBendChi2Size> bendChi2Bins = {
      {0., 0.5, 1.25, 2., 3., 5., 10., 50.}};

  // Track flags
  typedef ap_uint<TrackBitWidths::kValidSize> valid_t;  // valid bit

  // Track parameters types
  typedef ap_uint<TrackBitWidths::kRinvSize> rinv_t;  // Track Rinv
  typedef ap_uint<TrackBitWidths::kPhiSize> phi_t;    // Track phi
  typedef ap_uint<TrackBitWidths::kTanlSize> tanl_t;  // Track tan(l)
  typedef ap_uint<TrackBitWidths::kZ0Size> z0_t;      // Track z
  typedef ap_uint<TrackBitWidths::kD0Size> d0_t;      // D0

  // Track quality types
  typedef ap_uint<TrackBitWidths::kChi2RPhiSize> chi2rphi_t;      // Chi2 r-phi
  typedef ap_uint<TrackBitWidths::kChi2RZSize> chi2rz_t;          // Chi2 r-z
  typedef ap_uint<TrackBitWidths::kBendChi2Size> bendChi2_t;      // Bend-Chi2
  typedef ap_uint<TrackBitWidths::kHitPatternSize> hit_t;         // Hit mask bits
  typedef ap_uint<TrackBitWidths::kMVAQualitySize> qualityMVA_t;  // Track quality MVA
  typedef ap_uint<TrackBitWidths::kMVAOtherSize> otherMVA_t;      // Specialized MVA selection

  // Track word types
  typedef ap_uint<TrackBitWidths::kTrackWordSize> tkword_t;  // Entire track word;

public:
  // ----------Constructors --------------------------
  TTTrack_TrackWord() {}
  TTTrack_TrackWord(unsigned int valid,
                    const GlobalVector& momentum,
                    const GlobalPoint& POCA,
                    double rInv,
                    double chi2RPhi,  // would be xy chisq if chi2Z is non-zero
                    double chi2RZ,
                    double bendChi2,
                    unsigned int hitPattern,
                    unsigned int mvaQuality,
                    unsigned int mvaOther);
  TTTrack_TrackWord(unsigned int valid,
                    unsigned int rInv,
                    unsigned int phi0,
                    unsigned int tanl,
                    unsigned int z0,
                    unsigned int d0,
                    unsigned int chi2RPhi,  // would be total chisq if chi2Z is zero
                    unsigned int chi2RZ,
                    unsigned int bendChi2,
                    unsigned int hitPattern,
                    unsigned int mvaQuality,
                    unsigned int mvaOther);

  // ----------copy constructor ----------------------
  TTTrack_TrackWord(const TTTrack_TrackWord& word) { trackWord_ = word.trackWord_; }

  // ----------operators -----------------------------
  TTTrack_TrackWord& operator=(const TTTrack_TrackWord& word) {
    trackWord_ = word.trackWord_;
    return *this;
  }

  // ----------member functions (getters) ------------
  // These functions return arbitarary precision unsigned int words (lists of bits) for each quantity
  // Signed quantities have the sign enconded in the left-most bit.
  ap_uint<TrackBitWidths::kValidSize> getValidWord() const {
    return trackWord_(TrackBitLocations::kValidMSB, TrackBitLocations::kValidLSB);
  }
  ap_uint<TrackBitWidths::kRinvSize> getRinvWord() const {
    return trackWord_(TrackBitLocations::kRinvMSB, TrackBitLocations::kRinvLSB);
  }
  ap_uint<TrackBitWidths::kPhiSize> getPhiWord() const {
    return trackWord_(TrackBitLocations::kPhiMSB, TrackBitLocations::kPhiLSB);
  }
  ap_uint<TrackBitWidths::kTanlSize> getTanlWord() const {
    return trackWord_(TrackBitLocations::kTanlMSB, TrackBitLocations::kTanlLSB);
  }
  ap_uint<TrackBitWidths::kZ0Size> getZ0Word() const {
    return trackWord_(TrackBitLocations::kZ0MSB, TrackBitLocations::kZ0LSB);
  }
  ap_uint<TrackBitWidths::kD0Size> getD0Word() const {
    return trackWord_(TrackBitLocations::kD0MSB, TrackBitLocations::kD0LSB);
  }
  ap_uint<TrackBitWidths::kChi2RPhiSize> getChi2RPhiWord() const {
    return trackWord_(TrackBitLocations::kChi2RPhiMSB, TrackBitLocations::kChi2RPhiLSB);
  }
  ap_uint<TrackBitWidths::kChi2RZSize> getChi2RZWord() const {
    return trackWord_(TrackBitLocations::kChi2RZMSB, TrackBitLocations::kChi2RZLSB);
  }
  ap_uint<TrackBitWidths::kBendChi2Size> getBendChi2Word() const {
    return trackWord_(TrackBitLocations::kBendChi2MSB, TrackBitLocations::kBendChi2LSB);
  }
  ap_uint<TrackBitWidths::kHitPatternSize> getHitPatternWord() const {
    return trackWord_(TrackBitLocations::kHitPatternMSB, TrackBitLocations::kHitPatternLSB);
  }
  ap_uint<TrackBitWidths::kMVAQualitySize> getMVAQualityWord() const {
    return trackWord_(TrackBitLocations::kMVAQualityMSB, TrackBitLocations::kMVAQualityLSB);
  }
  ap_uint<TrackBitWidths::kMVAOtherSize> getMVAOtherWord() const {
    return trackWord_(TrackBitLocations::kMVAOtherMSB, TrackBitLocations::kMVAOtherLSB);
  }
  ap_uint<TrackBitWidths::kTrackWordSize> getTrackWord() const { return trackWord_; }

  // These functions return the packed bits in integer format for each quantity
  // Signed quantities have the sign enconded in the left-most bit.
  unsigned int getValidBits() const { return getValidWord().to_uint(); }
  unsigned int getRinvBits() const { return getRinvWord().to_uint(); }
  unsigned int getPhiBits() const { return getPhiWord().to_uint(); }
  unsigned int getTanlBits() const { return getTanlWord().to_uint(); }
  unsigned int getZ0Bits() const { return getZ0Word().to_uint(); }
  unsigned int getD0Bits() const { return getD0Word().to_uint(); }
  unsigned int getChi2RPhiBits() const { return getChi2RPhiWord().to_uint(); }
  unsigned int getChi2RZBits() const { return getChi2RZWord().to_uint(); }
  unsigned int getBendChi2Bits() const { return getBendChi2Word().to_uint(); }
  unsigned int getHitPatternBits() const { return getHitPatternWord().to_uint(); }
  unsigned int getMVAQualityBits() const { return getMVAQualityWord().to_uint(); }
  unsigned int getMVAOtherBits() const { return getMVAOtherWord().to_uint(); }

  // These functions return the unpacked and converted values
  // These functions return real numbers converted from the digitized quantities by unpacking the 96-bit track word
  bool getValid() const { return getValidWord().to_bool(); }
  double getRinv() const { return unpackSignedValue(getRinvBits(), TrackBitWidths::kRinvSize, stepRinv); }
  double getPhi() const { return unpackSignedValue(getPhiBits(), TrackBitWidths::kPhiSize, stepRinv); }
  double getTanl() const { return unpackSignedValue(getTanlBits(), TrackBitWidths::kTanlSize, stepRinv); }
  double getZ0() const { return unpackSignedValue(getZ0Bits(), TrackBitWidths::kZ0Size, stepRinv); }
  double getD0() const { return unpackSignedValue(getD0Bits(), TrackBitWidths::kD0Size, stepRinv); }
  double getChi2RPhi() const { return chi2RPhiBins[getChi2RPhiBits()]; }
  double getChi2RZ() const { return chi2RZBins[getChi2RZBits()]; }
  double getBendChi2() const { return bendChi2Bins[getBendChi2Bits()]; }
  unsigned int getHitPattern() const { return getHitPatternBits(); }
  unsigned int getMVAQuality() const { return getMVAQualityBits(); }
  unsigned int getMVAOther() const { return getMVAOtherBits(); }

  // ----------member functions (setters) ------------
  void setTrackWord(unsigned int valid,
                    const GlobalVector& momentum,
                    const GlobalPoint& POCA,
                    double rInv,
                    double chi2RPhi,  // would be total chisq if chi2Z is zero
                    double chi2RZ,
                    double bendChi2,
                    unsigned int hitPattern,
                    unsigned int mvaQuality,
                    unsigned int mvaOther);

  void setTrackWord(unsigned int valid,
                    unsigned int rInv,
                    unsigned int phi0,
                    unsigned int tanl,
                    unsigned int z0,
                    unsigned int d0,
                    unsigned int chi2RPhi,  // would be total chisq if chi2Z is zero
                    unsigned int chi2RZ,
                    unsigned int bendChi2,
                    unsigned int hitPattern,
                    unsigned int mvaQuality,
                    unsigned int mvaOther);

  void setTrackWord(ap_uint<TrackBitWidths::kValidSize> valid,
                    ap_uint<TrackBitWidths::kRinvSize> rInv,
                    ap_uint<TrackBitWidths::kPhiSize> phi0,
                    ap_uint<TrackBitWidths::kTanlSize> tanl,
                    ap_uint<TrackBitWidths::kZ0Size> z0,
                    ap_uint<TrackBitWidths::kD0Size> d0,
                    ap_uint<TrackBitWidths::kChi2RPhiSize> chi2RPhi,  // would be total chisq if chi2Z is zero
                    ap_uint<TrackBitWidths::kChi2RZSize> chi2RZ,
                    ap_uint<TrackBitWidths::kBendChi2Size> bendChi2,
                    ap_uint<TrackBitWidths::kHitPatternSize> hitPattern,
                    ap_uint<TrackBitWidths::kMVAQualitySize> mvaQuality,
                    ap_uint<TrackBitWidths::kMVAOtherSize> mvaOther);

private:
  // ----------private member functions --------------
  unsigned int digitizeSignedValue(double value, unsigned int nBits, double lsb) const {
    unsigned int digitized_value = std::floor(std::abs(value) / lsb);
    unsigned int digitized_maximum = (1 << (nBits - 1)) - 1;  // The remove 1 bit from nBits to account for the sign
    if (digitized_value > digitized_maximum)
      digitized_value = digitized_maximum;
    if (value < 0)
      digitized_value = (1 << nBits) - digitized_value;  // two's complement encoding
    return digitized_value;
  }

  template <typename T>
  constexpr unsigned int getBin(double value, const T& bins) const {
    auto up = std::upper_bound(bins.begin(), bins.end(), value);
    return (up - bins.begin() - 1);
  }

  double unpackSignedValue(unsigned int bits, unsigned int nBits, double lsb) const {
    int isign = 1;
    unsigned int digitized_maximum = (1 << nBits) - 1;
    if (bits & (1 << (nBits - 1))) {  // check the sign
      isign = -1;
      bits = (1 << (nBits + 1)) - bits;  // if negative, flip everything for two's complement encoding
    }
    return (double(bits & digitized_maximum) + 0.5) * lsb * isign;
  }

  // ----------member data ---------------------------
  tkword_t trackWord_;
};

#endif

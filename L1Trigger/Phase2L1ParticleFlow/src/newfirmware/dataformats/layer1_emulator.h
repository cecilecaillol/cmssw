#ifndef FIRMWARE_dataformats_layer1_emulator_h
#define FIRMWARE_dataformats_layer1_emulator_h

#include <fstream>
#include <vector>
#include "layer1_objs.h"
#include "pf.h"
#include "puppi.h"

namespace l1t { class PFTrack; class PFCluster; class PFCandidate; class Muon; }

namespace l1ct {
    
    struct HadCaloObjEmu : public HadCaloObj {
        const l1t::PFCluster *src;
        bool read(std::fstream & from) ;
        bool write(std::fstream & to) const ;
        void clear() {
            HadCaloObj::clear();
            src = nullptr;
        }
    };

    struct EmCaloObjEmu : public EmCaloObj { 
        const l1t::PFCluster *src;
        bool read(std::fstream & from) ;
        bool write(std::fstream & to) const ;
        void clear() {
            EmCaloObj::clear();
            src = nullptr;
        }
    };

    struct TkObjEmu : public TkObj {
        uint16_t hwChi2, hwStubs;
        const l1t::PFTrack *src;
        bool read(std::fstream & from) ;
        bool write(std::fstream & to) const ;
        void clear() {
            TkObj::clear();
            src = nullptr;
            hwChi2 = 0;
            hwStubs = 0;
        }
    };

    struct MuObjEmu : public MuObj {
        const l1t::Muon *src;
        bool read(std::fstream & from) ;
        bool write(std::fstream & to) const ;
        void clear() {
            MuObj::clear();
            src = nullptr;
        }
   };

    struct PFChargedObjEmu : public PFChargedObj {
        const l1t::PFCluster *srcCluster;
        const l1t::PFTrack *srcTrack;
        const l1t::Muon *srcMu;
        bool read(std::fstream & from) ;
        bool write(std::fstream & to) const ;
        void clear() {
            PFChargedObj::clear();
            srcCluster = nullptr;
            srcTrack = nullptr;
            srcMu = nullptr;
        }
    };

    struct PFNeutralObjEmu : public PFNeutralObj {
        const l1t::PFCluster *srcCluster;
        bool read(std::fstream & from) ;
        bool write(std::fstream & to) const ;
        void clear() {
            PFNeutralObj::clear();
            srcCluster = nullptr;
        }
    };

    struct PFRegionEmu : public PFRegion {
        float etaExtra, phiExtra;

        PFRegionEmu() : PFRegion() {}
        PFRegionEmu(float etamin,
           float etamax,
           float phicenter,
           float phiwidth,
           float etaextra,
           float phiextra) ;

        // global coordinates
        bool contains(float eta, float phi) const ;
        float localEta(float globalEta) const ;
        float localPhi(float globalPhi) const ;

        bool read(std::fstream & from) ;
        bool write(std::fstream & to) const ;
    };


    struct PuppiObjEmu : public PuppiObj {
        const l1t::PFCluster *srcCluster;
        const l1t::PFTrack *srcTrack;
        const l1t::Muon *srcMu;
        bool read(std::fstream & from) ;
        bool write(std::fstream & to) const ;
        void clear() {
            PuppiObj::clear();
            srcCluster = nullptr;
            srcTrack = nullptr;
            srcMu = nullptr;
        }
        inline void fill(const PFRegionEmu & region, const PFChargedObjEmu &src) {
            PuppiObj::fill(region, src);
            srcCluster = src.srcCluster;
            srcTrack = src.srcTrack;
            srcMu = src.srcMu;
        }
        inline void fill(const PFRegionEmu & region, const PFNeutralObjEmu &src, pt_t puppiPt,
                puppiWgt_t puppiWgt) {
            PuppiObj::fill(region, src, puppiPt, puppiWgt);
            srcCluster = src.srcCluster;
            srcTrack = nullptr;
            srcMu = nullptr;
        }
        inline void fill(const PFRegionEmu & region, const HadCaloObjEmu &src, pt_t puppiPt,
                puppiWgt_t puppiWgt) {
            PuppiObj::fill(region, src, puppiPt, puppiWgt);
            srcCluster = src.src;
            srcTrack = nullptr;
            srcMu = nullptr;
        }

    };

   
    struct PVObjEmu {
        z0_t hwZ0;
        bool read(std::fstream & from) ;
        bool write(std::fstream & to) const ;
    };

    template<typename T>
    struct DetectorSector {
        PFRegionEmu region;
        std::vector<T> obj;
        DetectorSector() {}
        DetectorSector(float etamin, float etamax, float phicenter, float phiwidth, float etaextra=0, float phiextra=0) : 
               region(etamin, etamax, phicenter, phiwidth, etaextra, phiextra) {}
        // convenience forwarding of some methods
        typedef typename std::vector<T>::const_iterator const_iterator;
        typedef typename std::vector<T>::iterator iterator;
        inline const T & operator[](unsigned int i) const { return obj[i]; }
        inline T & operator[](unsigned int i) { return obj[i]; }
        inline const_iterator begin() const { return obj.begin(); }
        inline iterator begin() { return obj.begin(); }
        inline const_iterator end() const { return obj.end(); }
        inline iterator end() { return obj.end(); }
        inline unsigned int size() const { return obj.size(); }
        inline void resize(unsigned int size) { obj.resize(size); }
        inline void clear() { obj.clear(); }

    };

    struct RegionizerDecodedInputs {
        std::vector<DetectorSector<HadCaloObjEmu>> hadcalo;
        std::vector<DetectorSector<EmCaloObjEmu>> emcalo;
        std::vector<DetectorSector<TkObjEmu>> track;
        std::vector<MuObjEmu> muon; // muons are global

        bool read(std::fstream & from) ;
        bool write(std::fstream & to) const ;
        void clear() ;
    }; 

    struct PFInputRegion {
        PFRegionEmu region;
        std::vector<HadCaloObjEmu> hadcalo;
        std::vector<EmCaloObjEmu> emcalo;
        std::vector<TkObjEmu> track;
        std::vector<MuObjEmu> muon;

        PFInputRegion() {}
        PFInputRegion(float etamin, float etamax, float phicenter, float phiwidth, float etaextra, float phiextra) : 
               region(etamin, etamax, phicenter, phiwidth, etaextra, phiextra) {}
        bool read(std::fstream & from) ;
        bool write(std::fstream & to) const ;
        void clear() ;
    }; 

    struct OutputRegion {
        std::vector<PFChargedObjEmu> pfcharged;
        std::vector<PFNeutralObjEmu> pfphoton;
        std::vector<PFNeutralObjEmu> pfneutral;
        std::vector<PFChargedObjEmu> pfmuon;
        std::vector<PuppiObjEmu> puppi;

        bool read(std::fstream & from) ;
        bool write(std::fstream & to) const ;
        void clear() ;
    };


    struct Event {
        static const int VERSION = 2;
        uint32_t run, lumi; uint64_t event;
        RegionizerDecodedInputs decoded;
        std::vector<PFInputRegion> pfinputs;
        std::vector<PVObjEmu> pvs;
        std::vector<OutputRegion> out;

        Event() : run(0), lumi(0), event(0) {}

        bool read(std::fstream & from) ;
        bool write(std::fstream & to) const ;
        void clear() ;
        void init(uint32_t run, uint32_t lumi, uint64_t event) ;
    };
  
    template<typename T1, typename T2>
    void toFirmware(const std::vector<T1> &in, unsigned int NMAX, T2 out[/*NMAX*/]) {
        unsigned int n = std::min<unsigned>(in.size(), NMAX);
        for (unsigned int i = 0; i < n; ++i) out[i] = in[i];
        for (unsigned int i = n; i < NMAX; ++i) out[i].clear();
    }

} // namespace

#endif

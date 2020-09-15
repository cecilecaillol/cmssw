#ifndef __L1Analysis_L1AnalysisPhaseIIStep1DataFormat_H__
#define __L1Analysis_L1AnalysisPhaseIIStep1DataFormat_H__

//-------------------------------------------------------------------------------
// Created 20/04/2010 - E. Conte, A.C. Le Bihan
// 
// 
// Original code : UserCode/L1TriggerDPG/L1ExtraTreeProducer - Jim Brooke
//-------------------------------------------------------------------------------


#include <vector>

namespace L1Analysis
{
  struct L1AnalysisPhaseIIStep1DataFormat
  {
    L1AnalysisPhaseIIStep1DataFormat(){Reset();};
    ~L1AnalysisPhaseIIStep1DataFormat(){};
    
    void Reset()
    {

      z0Puppi=0;
      z0VertexTDR=0;
      z0Vertices.clear();
      z0L1TkPV.clear();
      sumL1TkPV.clear();
      nL1TkPVs=0;
      nVertices=0;

      nCaloTaus = 0;
      caloTauEt.clear();
      caloTauEta.clear();
      caloTauPhi.clear(); 
      caloTauIEt.clear();
      caloTauIEta.clear();
      caloTauIPhi.clear(); 
      caloTauIso.clear();
      caloTauBx.clear();
      caloTauTowerIPhi.clear();
      caloTauTowerIEta.clear();
      caloTauRawEt.clear();
      caloTauIsoEt.clear();
      caloTauNTT.clear();
      caloTauHasEM.clear();
      caloTauIsMerged.clear();
      caloTauHwQual.clear();

      nPfPhase1L1Jets = 0;
      pfPhase1L1JetEt.clear();
      pfPhase1L1JetEta.clear();
      pfPhase1L1JetPhi.clear();

      pfPhase1L1HT.clear();
      pfPhase1L1MHTEt.clear();
      pfPhase1L1MHTPhi.clear();
      nPfPhase1L1MHT=0;
 
      nEG = 0;
      EGEt.clear();
      EGEta.clear();
      EGPhi.clear();
      EGBx.clear();
      EGIso.clear();
      EGzVtx.clear();
      EGHwQual.clear();      
      EGHGC.clear();
      EGPassesLooseTrackID.clear();
      EGPassesPhotonID.clear();

      nTkElectrons = 0;
      tkElectronEt.clear();
      tkElectronEta.clear();
      tkElectronPhi.clear();
      tkElectronChg.clear();
      tkElectronBx.clear();
      tkElectronTrkIso.clear();
      tkElectronzVtx.clear();
      tkElectronHwQual.clear();
      tkElectronEGRefPt.clear(); 
      tkElectronEGRefEta.clear();
      tkElectronEGRefPhi.clear();
      tkElectronHGC.clear();
      tkElectronPassesLooseTrackID.clear();
      tkElectronPassesPhotonID.clear();

      nTkPhotons = 0;
      tkPhotonEt.clear();
      tkPhotonEta.clear();
      tkPhotonPhi.clear();
      tkPhotonBx.clear();
      tkPhotonTrkIso.clear();
      tkPhotonTrkIsoPV.clear();
      tkPhotonzVtx.clear();
      tkPhotonHwQual.clear();
      tkPhotonEGRefPt.clear();
      tkPhotonEGRefEta.clear();
      tkPhotonEGRefPhi.clear();
      tkPhotonHGC.clear();
      tkPhotonPassesLooseTrackID.clear();
      tkPhotonPassesPhotonID.clear();

      nTkMuons = 0;
      tkMuonPt.clear();
      tkMuonEta.clear();
      tkMuonPhi.clear();
      tkMuonChg.clear();
      tkMuonTrkIso.clear();
      tkMuonBx.clear();
      tkMuonQual.clear();
      tkMuonzVtx.clear();
      tkMuonMuRefPt.clear();
      tkMuonTrkRefPt.clear();
      tkMuonMuRefPhi.clear(); 
      tkMuonMuRefEta.clear();
      tkMuonDRMuTrack.clear();
      tkMuonNMatchedTracks.clear();
      tkMuonMuRefChg.clear();
      tkMuonRegion.clear();

      nTkMuonStubs = 0;
      tkMuonStubsPt.clear();
      tkMuonStubsEta.clear();
      tkMuonStubsPhi.clear();
      tkMuonStubsChg.clear();
      tkMuonStubsTrkIso.clear();
      tkMuonStubsBx.clear();
      tkMuonStubsQual.clear();
      tkMuonStubszVtx.clear();
      tkMuonStubsBarrelStubs.clear();
      tkMuonStubsRegion.clear();

      puppiMETEt=0;
      puppiMETPhi=0;

      nNNTaus = 0;
      nnTauEt.clear();
      nnTauEta.clear();
      nnTauPhi.clear();
      nnTauChg.clear();
      nnTauChargedIso.clear();
      nnTauFullIso.clear();    
      nnTauID.clear();         
      nnTauPassLooseNN.clear();   
      nnTauPassLoosePF.clear(); 
      nnTauPassTightPF.clear();   
      nnTauPassTightNN.clear();

      nNNTauPFs = 0;
      nnTauPFEt.clear();
      nnTauPFEta.clear();
      nnTauPFPhi.clear();
      nnTauPFChg.clear();
      nnTauPFChargedIso.clear();
      nnTauPFFullIso.clear();
      nnTauPFID.clear();
      nnTauPFPassLooseNN.clear();
      nnTauPFPassLoosePF.clear();
      nnTauPFPassTightPF.clear();
      nnTauPFPassTightNN.clear();


      nTkBsCands=0;
      tkBsCandPt.clear();
      tkBsCandEta.clear();
      tkBsCandPhi.clear();
      tkBsCandMass.clear();
      tkBsCandPhi1Pt.clear();
      tkBsCandPhi2Pt.clear();
      tkBsCandPhi1Eta.clear();
      tkBsCandPhi2Eta.clear();
      tkBsCandPhi1Phi.clear();
      tkBsCandPhi2Phi.clear();
      tkBsCandPhi1Mass.clear();
      tkBsCandPhi2Mass.clear();
      tkBsCandDRPhiPair.clear(); 
      tkBsCandDxyPhiPair.clear();
      tkBsCandDzPhiPair.clear(); 
      tkBsCandKind.clear();

    }
 
    double z0Puppi;
    double z0VertexTDR;
    unsigned short int nVertices;
    std::vector<double>z0Vertices;  
    unsigned short int nL1TkPVs;
    std::vector<double>  z0L1TkPV;
    std::vector<double>  sumL1TkPV;
 
    unsigned short int nCaloTaus;
    std::vector<double> caloTauEt;
    std::vector<double> caloTauEta;
    std::vector<double> caloTauPhi;
    std::vector<short int> caloTauIEt;
    std::vector<short int> caloTauIEta;
    std::vector<short int> caloTauIPhi;
    std::vector<short int> caloTauIso;
    std::vector<short int> caloTauBx;
    std::vector<short int> caloTauTowerIPhi;
    std::vector<short int> caloTauTowerIEta;
    std::vector<short int> caloTauRawEt;    
    std::vector<short int> caloTauIsoEt;
    std::vector<short int> caloTauNTT;
    std::vector<short int> caloTauHasEM;
    std::vector<short int> caloTauIsMerged;
    std::vector<short int> caloTauHwQual;


    unsigned short int nPfPhase1L1Jets;
    std::vector<double> pfPhase1L1JetEt;
    std::vector<double> pfPhase1L1JetEta;
    std::vector<double> pfPhase1L1JetPhi;

    std::vector<double> pfPhase1L1HT;
    std::vector<double> pfPhase1L1MHTEt;
    std::vector<double> pfPhase1L1MHTPhi;
     unsigned int nPfPhase1L1MHT;


    unsigned int nEG;
    std::vector<double> EGEt;
    std::vector<double> EGEta;
    std::vector<double> EGPhi;
    std::vector<int>    EGBx;
    std::vector<double> EGIso;
    std::vector<double> EGzVtx;
    std::vector<int>    EGHwQual;
    std::vector<unsigned int> EGHGC;
    std::vector<unsigned int>   EGPassesLooseTrackID;
    std::vector<unsigned int>   EGPassesPhotonID;

    unsigned int nTkElectrons;
    std::vector<double> tkElectronEt;
    std::vector<double> tkElectronEta;
    std::vector<double> tkElectronPhi;
    std::vector<int>    tkElectronChg;
    std::vector<int>    tkElectronBx;
    std::vector<double> tkElectronTrkIso;
    std::vector<double> tkElectronzVtx;
    std::vector<double> tkElectronHwQual;
    std::vector<double>   tkElectronEGRefPt;
    std::vector<double>   tkElectronEGRefEta;
    std::vector<double>   tkElectronEGRefPhi;
    std::vector<unsigned int> tkElectronHGC;
    std::vector<unsigned int> tkElectronPassesLooseTrackID;
    std::vector<unsigned int> tkElectronPassesPhotonID;

    unsigned int nTkPhotons;
    std::vector<double> tkPhotonEt;
    std::vector<double> tkPhotonEta;
    std::vector<double> tkPhotonPhi;
    std::vector<int>    tkPhotonBx;
    std::vector<double> tkPhotonTrkIso;
    std::vector<double> tkPhotonTrkIsoPV;
    std::vector<double> tkPhotonzVtx;
    std::vector<double> tkPhotonHwQual;
    std::vector<double>   tkPhotonEGRefPt;
    std::vector<double>   tkPhotonEGRefEta;
    std::vector<double>   tkPhotonEGRefPhi;
    std::vector<unsigned int> tkPhotonHGC;
    std::vector<unsigned int> tkPhotonPassesLooseTrackID;
    std::vector<unsigned int> tkPhotonPassesPhotonID;


    unsigned int nTkMuons;
    std::vector<double>   tkMuonPt;
    std::vector<double>   tkMuonEta;
    std::vector<double>   tkMuonPhi;
    std::vector<int>      tkMuonChg;
    std::vector<unsigned int> tkMuonIso;
    std::vector<double> tkMuonTrkIso;
    std::vector<unsigned int> tkMuonFwd;
    std::vector<unsigned int> tkMuonMip;
    std::vector<unsigned int> tkMuonRPC;
    std::vector<int>      tkMuonBx;
    std::vector<unsigned int>      tkMuonQual;
    std::vector<double>   tkMuonzVtx;
    std::vector<double> tkMuonMuRefPt;
    std::vector<double> tkMuonTrkRefPt;
    std::vector<double>  tkMuonMuRefPhi;
    std::vector<double>  tkMuonMuRefEta;
    std::vector<double>  tkMuonDRMuTrack;
    std::vector<double>  tkMuonNMatchedTracks;
    std::vector<int>  tkMuonMuRefChg;
    std::vector<unsigned int>   tkMuonRegion;

    unsigned int nTkMuonStubs;
    std::vector<double>   tkMuonStubsPt;
    std::vector<double>   tkMuonStubsEta;
    std::vector<double>   tkMuonStubsPhi;
    std::vector<int>      tkMuonStubsChg;
    std::vector<int>      tkMuonStubsBx;
    std::vector<double>   tkMuonStubsTrkIso;
    std::vector<unsigned int>      tkMuonStubsQual;
    std::vector<double>   tkMuonStubszVtx;
    std::vector<double>   tkMuonStubsBarrelStubs;
    std::vector<unsigned int>   tkMuonStubsRegion;

    double puppiMETEt;
    double puppiMETPhi;

    unsigned int nNNTaus;
    std::vector<double>   nnTauEt;
    std::vector<double>   nnTauEta;
    std::vector<double>   nnTauPhi;
    std::vector<int>   nnTauChg;
    std::vector<double>   nnTauChargedIso;
    std::vector<double>   nnTauFullIso;
    std::vector<unsigned int>   nnTauID;
    std::vector<unsigned int>   nnTauPassLooseNN;
    std::vector<unsigned int>   nnTauPassLoosePF;
    std::vector<unsigned int>   nnTauPassTightPF;
    std::vector<unsigned int>   nnTauPassTightNN;


    unsigned int nNNTauPFs;
    std::vector<double>   nnTauPFEt;
    std::vector<double>   nnTauPFEta;
    std::vector<double>   nnTauPFPhi;
    std::vector<int>   nnTauPFChg;
    std::vector<double>   nnTauPFChargedIso;
    std::vector<double>   nnTauPFFullIso;
    std::vector<unsigned int>   nnTauPFID;
    std::vector<unsigned int>   nnTauPFPassLooseNN;
    std::vector<unsigned int>   nnTauPFPassLoosePF;
    std::vector<unsigned int>   nnTauPFPassTightPF;
    std::vector<unsigned int>   nnTauPFPassTightNN;


   unsigned int nTkBsCands;
   std::vector<double> tkBsCandPt;
   std::vector<double> tkBsCandEta;
   std::vector<double> tkBsCandPhi;
   std::vector<double> tkBsCandMass;
   std::vector<double> tkBsCandPhi1Pt;
   std::vector<double> tkBsCandPhi2Pt;
   std::vector<double> tkBsCandPhi1Eta;
   std::vector<double> tkBsCandPhi2Eta;
   std::vector<double> tkBsCandPhi1Phi;
   std::vector<double> tkBsCandPhi2Phi;
   std::vector<double> tkBsCandPhi1Mass;
   std::vector<double> tkBsCandPhi2Mass;
   std::vector<double> tkBsCandDRPhiPair;    
   std::vector<double> tkBsCandDxyPhiPair;   
   std::vector<double> tkBsCandDzPhiPair;    
   std::vector<int>    tkBsCandKind;

  }; 
}
#endif



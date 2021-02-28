#include "L1Trigger/L1TTrackMatch/interface/L1TkEtMissEmuAlgo.h"

namespace L1TkEtMissEmuAlgo{
  //Function to count stubs in hitpattern
  iNstub CountNStub(unsigned int Hitpattern){  
    iNstub Nstub = 0;
    for (int i = (N_hitpatternbits-1); i >= 0; i--) {
      int k = Hitpattern >> i;
        if (k & 1)
          Nstub ++;
    }
    return Nstub;
  }

  //Function to take float phi to local integer phi -> would be part of tttrack
  iPhi FloatPhi_to_iPhi(float phi,unsigned int sector){   
    float temp_phi = 0.0;
    if (sector < 4 ){ 
      temp_phi = phi - (sector * (2*M_PI)/9); 
    } else if (sector > 5 ){
      temp_phi = phi + ((9-sector) * (2*M_PI)/9); 
    } else if (sector == 4 ){
      if (phi > 0){
        temp_phi = phi - (sector *  (2*M_PI)/9);
      } else{ 
        temp_phi = phi + ((9-sector) * (2*M_PI)/9);
      }
    } else if (sector == 5 ){
      if (phi < 0){
        temp_phi = phi + ((9-sector) * (2*M_PI)/9); 
      } else{
         temp_phi = phi - (sector *  (2*M_PI)/9);
      }
    }
    return digitize_Signed <iPhi> (temp_phi,-maxPhi,maxPhi,N_phiBins);
  }


  void FillCosLUT(iglobPhi cosLUT[LUT_bins],float max_LUT_phi){  //Fill cosine LUT with integer values
    float phi = 0;
    for (int LUT_idx = 0; LUT_idx < LUT_bins; LUT_idx++){
      cosLUT[LUT_idx] = digitize_Signed <iglobPhi> (cos(phi),0,1,N_globPhiBins);
      phi += (2*maxPhi)/(N_globPhiBins -1 );
    }
  }

    // For generating phi slices for LUTs
  std::vector<iglobPhi> generate_phi_slice_LUTs(unsigned int Num_entries){
    std::vector<iglobPhi> table;
    float slice_centre = 0.0;
    for (unsigned int q=0;q<Num_entries;q++){
      table.push_back(slice_centre / (2*maxPhi / (N_globPhiBins-1)));
      slice_centre += 2*M_PI/(Num_entries-1);
    }
    return table;
  }

  //Converts local int phi to global int phi
  iglobPhi local_to_global(iPhi local_phi,iglobPhi sector_shift,std::vector<iglobPhi> phi_quadrants){ 
      int PhiShift = N_globPhiBins/2;
      int PhiMin = phi_quadrants.front();
      int PhiMax = phi_quadrants.back();
      int phi_multiplier = N_phiBits - N_globphiBits;

      int tempPhi = 0;
      iglobPhi globalPhi = 0;

      tempPhi = (local_phi / pow(2,phi_multiplier)) + sector_shift - PhiShift;
      if (tempPhi < PhiMin) {tempPhi = tempPhi + PhiMax;}
      else if(tempPhi > PhiMax) {tempPhi = tempPhi - PhiMax;}
      else tempPhi = tempPhi;

      globalPhi = iglobPhi (tempPhi);

      return globalPhi;

  }

    // Generate Eta LUT for track to vertex association
  std::vector<iEta> generate_EtaRegions(const float EtaRegions[N_etaregions+1]) {
    std::vector<iEta> iEtaRegions;
    for (int q=0;q<=N_etaregions;q++){
      iEtaRegions.push_back(digitize_Signed <iEta> (EtaRegions[q], -maxEta, maxEta, N_etaBins));
    }
    return iEtaRegions;
  }

  std::vector<iZ0> generate_DeltaZBins(const float DeltaZ[N_etaregions]) {
    std::vector<iZ0> iDeltaZ;
    for (int q=0;q<N_etaregions;q++){
      iDeltaZ.push_back(digitize_Signed <iZ0> (DeltaZ[q], 0, maxZ0, N_z0Bins/2));  //Only half range otherwise deltaZ values midway through integer representation
    } 
    return iDeltaZ;

  }





}
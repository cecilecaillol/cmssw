#include <L1Trigger/L1TMuonBayes/interface/Omtf/GoldenPattern.h>
#include <L1Trigger/L1TMuonBayes/interface/Omtf/GoldenPatternWithStat.h>
#include <L1Trigger/L1TMuonBayes/interface/Omtf/OMTFConfiguration.h>
#include <L1Trigger/L1TMuonBayes/interface/Omtf/OMTFSorter.h>
#include <cassert>
#include <iostream>
#include <strstream>
#include <algorithm>
#include <bitset>


#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/RPCTrigger/interface/RPCConst.h"
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
template <class GoldenPatternType>
AlgoMuons::value_type OMTFSorter<GoldenPatternType>::sortRefHitResults(unsigned int procIndx, unsigned int iRefHit, const std::vector< std::shared_ptr<GoldenPatternType> >& gPatterns,
					  int charge){
  GoldenPatternType* bestGP = 0; //the GoldenPattern with the best result for this iRefHit
//  std::cout <<" ====== sortRefHitResults: " << std::endl;

  for(auto& itGP: gPatterns) {
    if(!itGP->getResults()[procIndx][iRefHit].isValid())
      continue;

    if(charge!=0 && itGP->key().theCharge != charge)
      continue; //charge==0 means ignore charge

    ///Accept only candidates with >2 hits
    if(itGP->getResults()[procIndx][iRefHit].getFiredLayerCnt() < 3) //TODO - move 3 to the configuration??
      continue;

    if(bestGP == 0) {
      bestGP = itGP.get();
    }
    else if(itGP->getResults()[procIndx][iRefHit].getFiredLayerCnt() > bestGP->getResults()[procIndx][iRefHit].getFiredLayerCnt() ){
      bestGP = itGP.get();
      //std::cout <<" sorter, byQual, now best is: "<<bestKey << " RefLayer "<<itKey.second.getRefLayer()<<" FiredLayerCn "<<itKey.second.getFiredLayerCnt()<<std::endl;
    }
    else if(itGP->getResults()[procIndx][iRefHit].getFiredLayerCnt() == bestGP->getResults()[procIndx][iRefHit].getFiredLayerCnt() ) {
      if(itGP->getResults()[procIndx][iRefHit].getPdfSum() > bestGP->getResults()[procIndx][iRefHit].getPdfSum()) {
        //if the PdfWeigtSum is equal, we take the GP with the lower number, i.e. lower pt = check if this is ok for physics FIXME (KB)
        bestGP = itGP.get();
        //std::cout <<" sorter, byDisc, now best is: "<<bestKey << " "<<itKey.second.getRefLayer()<<" FiredLayerCn "<<itKey.second.getFiredLayerCnt()<< std::endl;
      }
    }
  }
  if(bestGP) {
    AlgoMuons::value_type candidate(new AlgoMuon(bestGP->getResults()[procIndx][iRefHit], bestGP, iRefHit) );
     //std::cout<<__FUNCTION__<<" line "<<__LINE__ <<" return: " << candidate << std::endl;
     return candidate;
  }
  else {
    AlgoMuons::value_type candidate(new AlgoMuon() );
    candidate->setRefHitNumber(iRefHit);
    return candidate;
  }
}

template class OMTFSorter<GoldenPattern>;
template class OMTFSorter<GoldenPatternWithStat>;
template class OMTFSorter<GoldenPatternWithThresh>;


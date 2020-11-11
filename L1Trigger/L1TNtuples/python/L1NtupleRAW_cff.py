import FWCore.ParameterSet.Config as cms

from L1Trigger.L1TNtuples.l1EventTree_cfi import *
from L1Trigger.L1TNtuples.l1ExtraTree_cfi import *
from L1Trigger.L1TNtuples.l1CaloTowerTree_cfi import *
from L1Trigger.L1TNtuples.l1UpgradeTfMuonTree_cfi import *
from L1Trigger.L1TNtuples.l1UpgradeTree_cfi import *
from L1Trigger.L1TNtuples.l1uGTTree_cfi import *
from L1Trigger.L1TNtuples.l1HOTree_cfi import *

L1NtupleRAW = cms.Sequence(
  l1EventTree
  #+l1ExtraTree
  +l1CaloTowerTree
  +l1UpgradeTfMuonTree
  +l1UpgradeTree
  +l1uGTTree
  +l1HOTree
)

#  do not have l1t::CaloTowerBxCollection in Stage1 
from Configuration.Eras.Modifier_stage1L1Trigger_cff import stage1L1Trigger
_stage1_L1NTupleRAW = L1NtupleRAW.copyAndExclude([l1CaloTowerTree,l1UpgradeTfMuonTree])
stage1L1Trigger.toReplaceWith(L1NtupleRAW,_stage1_L1NTupleRAW)

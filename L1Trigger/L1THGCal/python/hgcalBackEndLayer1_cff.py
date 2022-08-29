import FWCore.ParameterSet.Config as cms

from L1Trigger.L1THGCal.hgcalTriggerGeometryESProducer_cfi import *
from L1Trigger.L1THGCal.hgcalBackEndLayer1Producer_cfi import *


<<<<<<< HEAD
hgcalBackEndLayer1 = cms.Task(hgcalBackEndLayer1Producer)

hgcalBackEndLayer1HFNose = cms.Task(hgcalBackEndLayer1ProducerHFNose)

=======
L1THGCalBackEndLayer1 = cms.Task(l1tHGCalBackEndLayer1Producer)
L1THGCalBackEndStage1 = cms.Task(l1tHGCalBackEndStage1Producer)

L1THGCalBackEndLayer1HFNose = cms.Task(l1tHGCalBackEndLayer1ProducerHFNose)
>>>>>>> e68eee3787a... Rename L1T modules/sequences/tasks to follow conventions

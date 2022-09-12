import FWCore.ParameterSet.Config as cms

from L1Trigger.L1THGCal.l1tHGCalTriggerGeometryESProducer_cfi import *
from L1Trigger.L1THGCal.l1tHGCalBackEndLayer2Producer_cfi import *


<<<<<<< HEAD
hgcalBackEndLayer2 = cms.Task(hgcalBackEndLayer2Producer)

hgcalBackEndLayer2HFNose = cms.Task(hgcalBackEndLayer2ProducerHFNose)

=======
L1THGCalBackEndLayer2 = cms.Task(l1tHGCalBackEndLayer2Producer)
L1THGCalBackEndStage2 = cms.Task(l1tHGCalBackEndStage2Producer)

L1THGCalBackEndLayer2HFNose = cms.Task(l1tHGCalBackEndLayer2ProducerHFNose)
>>>>>>> e68eee3787a... Rename L1T modules/sequences/tasks to follow conventions

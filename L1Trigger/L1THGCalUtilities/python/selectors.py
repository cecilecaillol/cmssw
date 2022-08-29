import FWCore.ParameterSet.Config as cms

<<<<<<< HEAD
def create_genmatch(process, inputs,
        distance=0.3
        ):
    producer = process.hgc3DClusterGenMatchSelector.clone(
            dR = cms.double(distance),
            src = cms.InputTag('{}:HGCalBackendLayer2Processor3DClustering'.format(inputs))
            )
    return producer
=======
class CreateGenMatch(object):
    def __init__(self,
            distance=0.3
            ):
        self.dR = distance

    def __call__(self, process, inputs):
        producer = process.l1tHGCal3DClusterGenMatchSelector.clone(
                dR = cms.double(self.dR),
                src = cms.InputTag(inputs)
                )
        return producer
>>>>>>> e68eee3787a... Rename L1T modules/sequences/tasks to follow conventions

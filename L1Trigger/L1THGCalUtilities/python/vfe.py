import FWCore.ParameterSet.Config as cms

from L1Trigger.L1THGCal.l1tHGCalVFEProducer_cfi import vfe_proc

def create_vfe(process,
        linearization_si=vfe_proc.linearizationCfg_si,
        linearization_sc=vfe_proc.linearizationCfg_sc,
        compression_ldm=vfe_proc.compressionCfg_ldm,
        compression_hdm=vfe_proc.compressionCfg_hdm,
        ):
    producer = process.hgcalVFEProducer.clone(
        ProcessorParameters = vfe_proc.clone(
            linearizationCfg_si = linearization_si,
            linearizationCfg_sc = linearization_sc,
            compressionCfg_ldm = compression_ldm,
            compressionCfg_hdm = compression_hdm,
        )
<<<<<<< HEAD
    )
    return producer
=======

    def __call__(self, process):
        producer = process.l1tHGCalVFEProducer.clone(
            ProcessorParameters = self.processor
        )
        return producer
>>>>>>> e68eee3787a... Rename L1T modules/sequences/tasks to follow conventions

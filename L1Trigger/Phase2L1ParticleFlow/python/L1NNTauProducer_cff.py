import FWCore.ParameterSet.Config as cms

from L1Trigger.Phase2L1ParticleFlow.l1ctLayer1_cff import l1tLayer1Barrel,l1tLayer1HGCal,l1tLayer1

#from L1Trigger.Phase2L1ParticleFlow.L1NNTauProducer_cfi import *

#L1NNTauProducerPuppi = L1NNTauProducer.clone(
#                                NNFileName      = cms.string("L1Trigger/Phase2L1ParticleFlow/data/tau_3layer_puppi.pb")
#                                )


l1tNNTauProducerPuppi = cms.EDProducer("L1NNTauProducer",
                                      seedpt          = cms.double(10),
                                      conesize        = cms.double(0.4),
                                      tausize         = cms.double(0.1),
                                      maxtaus         = cms.int32(5),
                                      nparticles      = cms.int32(10),
<<<<<<< HEAD
                                      HW              = cms.bool(True),
                                      debug           = cms.bool(False),
                                      L1PFObjects     = cms.InputTag("l1ctLayer1:Puppi"), #1pfCandidates:Puppi"),#l1pfCandidates
=======
                                      L1PFObjects     = cms.InputTag("l1tLayer1:Puppi"), #1pfCandidates:Puppi"),#l1pfCandidates
>>>>>>> e68eee3787a... Rename L1T modules/sequences/tasks to follow conventions
                                      NNFileName      = cms.string("L1Trigger/Phase2L1ParticleFlow/data/tau_3layer_puppi.pb")
)

l1tNNTauProducerPF = cms.EDProducer("L1NNTauProducer",
                                      seedpt          = cms.double(10),
                                      conesize        = cms.double(0.4),
                                      tausize         = cms.double(0.1),
                                      maxtaus         = cms.int32(5),
                                      nparticles      = cms.int32(10),
<<<<<<< HEAD
                                      HW              = cms.bool(True),
                                      debug           = cms.bool(False),
                                      L1PFObjects     = cms.InputTag("l1ctLayer1:PF"),#l1pfCandidates
=======
                                      L1PFObjects     = cms.InputTag("l1tLayer1:PF"),#l1pfCandidates
>>>>>>> e68eee3787a... Rename L1T modules/sequences/tasks to follow conventions
                                      NNFileName      = cms.string("L1Trigger/Phase2L1ParticleFlow/data/tau_3layer.pb")
)


l1ctLayer1Barrel2Vtx = l1ctLayer1Barrel.clone()
l1ctLayer1Barrel2Vtx.nVtx = 2
l1ctLayer1Barrel2Vtx.puAlgoParameters.nVtx = 2
l1ctLayer1HGCal2Vtx  = l1ctLayer1HGCal.clone()
l1ctLayer1HGCal2Vtx.nVtx = 2
l1ctLayer1HGCal2Vtx.puAlgoParameters.nVtx = 2
l1ctLayer12Vtx       = l1ctLayer1.clone()
l1ctLayer12Vtx.pfProducers = cms.VInputTag(
    cms.InputTag("l1ctLayer1Barrel2Vtx"),
    cms.InputTag("l1ctLayer1HGCal2Vtx"),
    cms.InputTag("l1ctLayer1HGCalNoTK"),
    cms.InputTag("l1ctLayer1HF")
)
L1NNTauProducerPuppi2Vtx = L1NNTauProducerPuppi.clone()
L1NNTauProducerPuppi2Vtx.L1PFObjects =  cms.InputTag("l1ctLayer12Vtx:Puppi")
tau2VtxTaskHW = cms.Task(
    l1ctLayer1Barrel2Vtx,
    l1ctLayer1HGCal2Vtx,
    l1ctLayer12Vtx,
    L1NNTauProducerPuppi2Vtx
)

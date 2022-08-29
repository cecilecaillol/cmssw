import FWCore.ParameterSet.Config as cms

l1tSeedConePFJetProducer = cms.EDProducer("L1SeedConePFJetProducer",
                           L1PFObjects = cms.InputTag("l1tLayer1","Puppi"),
                           nJets       = cms.uint32(10),
                           coneSize    = cms.double(0.4),
                           HW          = cms.bool(False),
                           debug       = cms.bool(False)
                         )

<<<<<<< HEAD
L1SeedConePFJetEmulatorProducer = cms.EDProducer("L1SeedConePFJetProducer",
                                   L1PFObjects = cms.InputTag("l1ctLayer1","Puppi"),
                                   nJets       = cms.uint32(10),
                                   coneSize    = cms.double(0.4),
                                   HW          = cms.bool(True),
                                   debug       = cms.bool(False) 
                                  )
=======
l1tSeedConePFJetEmulatorProducer = l1tSeedConePFJetProducer.clone(HW = True)

>>>>>>> e68eee3787a... Rename L1T modules/sequences/tasks to follow conventions

import FWCore.ParameterSet.Config as cms

<<<<<<< HEAD
L1TrackerEtMiss = cms.EDProducer('L1TrackerEtMissProducer',
    L1TrackInputTag = cms.InputTag("L1TrackSelectionProducer", L1TrackSelectionProducer.outputCollectionName.value()),
    L1TrackAssociatedInputTag = cms.InputTag("L1TrackSelectionProducer", L1TrackSelectionProducer.outputCollectionName.value() + "Associated"),
    L1MetCollectionName = cms.string("L1TrackerEtMiss"),
    maxPt = cms.double( -10. ),	    # in GeV. When maxPt > 0, tracks with PT above maxPt are considered as
=======
l1tTrackerEtMiss = cms.EDProducer('L1TrackerEtMissProducer',
    L1TrackInputTag = cms.InputTag("l1tTrackSelectionProducer", "Level1TTTracksSelected"),
    L1TrackAssociatedInputTag = cms.InputTag("l1tTrackSelectionProducer", "Level1TTTracksSelectedAssociated"),
    L1MetCollectionName = cms.string("l1tTrackerEtMiss"),
    maxPt = cms.double(-10.) ,	    # in GeV. When maxPt > 0, tracks with PT above maxPt are considered as
>>>>>>> e68eee3787a... Rename L1T modules/sequences/tasks to follow conventions
                                    # mismeasured and are treated according to highPtTracks below.
                                    # When maxPt < 0, no special treatment is done for high PT tracks.
    highPtTracks = cms.int32( 1 ),  # when = 0 : truncation. Tracks with PT above maxPt are ignored
                                    # when = 1 : saturation. Tracks with PT above maxPt are set to PT=maxPt.
                                    # When maxPt < 0, no special treatment is done for high PT tracks.
    debug     = cms.bool(False)
)

<<<<<<< HEAD
L1TrackerEtMissExtended = L1TrackerEtMiss.clone( #NOT OPTIMIZED, STUDIED, OR USED
    L1TrackInputTag = cms.InputTag("L1TrackSelectionProducerExtended", L1TrackSelectionProducerExtended.outputCollectionName.value()),
    L1TrackAssociatedInputTag = cms.InputTag("L1TrackSelectionProducerExtended", L1TrackSelectionProducerExtended.outputCollectionName.value() + "Associated"),
    L1MetCollectionName = cms.string("L1TrackerExtendedEtMiss"),
=======
l1tTrackerEtMissExtended = l1tTrackerEtMiss.clone( #NOT OPTIMIZED, STUDIED, OR USED
    L1TrackInputTag = ("l1tTrackSelectionProducerExtended", "Level1TTTracksExtendedSelected"),
    L1TrackAssociatedInputTag = ("l1tTrackSelectionProducerExtended", "Level1TTTracksExtendedSelectedAssociated"),
    L1MetCollectionName = "l1tTrackerExtendedEtMiss",
>>>>>>> e68eee3787a... Rename L1T modules/sequences/tasks to follow conventions
)

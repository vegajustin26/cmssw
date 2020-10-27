import FWCore.ParameterSet.Config as cms

TrackFindingTrackletProducer_params = cms.PSet (

  InputTag         = cms.InputTag( "TTTracksFromTrackletEmulation", "Level1TTTracks"), #
  LabelKFin        = cms.string  ( "TrackFindingTrackletProducerKFin"  ),              #
  LabelKF          = cms.string  ( "TrackFindingTrackletProducerKF"    ),              #
  LabelTT          = cms.string  ( "TrackFindingTrackletProducerKFout" ),              #
  LabelAS          = cms.string  ( "TrackFindingTrackletProducerAS"    ),              #
  BranchAccepted   = cms.string  ( "StubAccepted" ),                                   #
  BranchLost       = cms.string  ( "StubLost"     ),                                   #
  CheckHistory     = cms.bool    ( True  ),                                            # checks if input sample production is configured as current process
  EnableTruncation = cms.bool    ( True  )                                             # enable emulation of truncation, lost stubs are filled in BranchLost

)
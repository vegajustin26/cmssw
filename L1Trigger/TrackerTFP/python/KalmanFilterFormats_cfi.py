import FWCore.ParameterSet.Config as cms

TrackTriggerKalmanFilterFormats_params = cms.PSet (

    BaseShiftx0           = cms.int32(  -3 ),
    BaseShiftx1           = cms.int32(  -9 ),
    BaseShiftx2           = cms.int32(  -5 ),
    BaseShiftx3           = cms.int32(  -6 ),
    BaseShiftv0           = cms.int32(  -6 ),
    BaseShiftv1           = cms.int32(   7 ),
    BaseShiftr0           = cms.int32( -12 ),
    BaseShiftr1           = cms.int32(  -4 ),
    BaseShiftS00          = cms.int32(  -2 ),
    BaseShiftS01          = cms.int32(  -7 ),
    BaseShiftS12          = cms.int32(   3 ),
    BaseShiftS13          = cms.int32(   1 ),
    BaseShiftK00          = cms.int32( -16 ),
    BaseShiftK10          = cms.int32( -23 ),
    BaseShiftK21          = cms.int32( -22 ),
    BaseShiftK31          = cms.int32( -23 ),
    BaseShiftR00          = cms.int32(  -5 ),
    BaseShiftR11          = cms.int32(   6 ),
    BaseShiftInvR00Approx = cms.int32( -26 ),
    BaseShiftInvR11Approx = cms.int32( -37 ),
    BaseShiftInvR00Cor    = cms.int32( -15 ),
    BaseShiftInvR11Cor    = cms.int32( -15 ),
    BaseShiftInvR00       = cms.int32( -26 ),
    BaseShiftInvR11       = cms.int32( -37 ),
    BaseShiftC00          = cms.int32(   5 ),
    BaseShiftC01          = cms.int32(  -3 ),
    BaseShiftC11          = cms.int32(  -7 ),
    BaseShiftC22          = cms.int32(   3 ),
    BaseShiftC23          = cms.int32(   0 ),
    BaseShiftC33          = cms.int32(   1 )

)
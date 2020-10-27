import FWCore.ParameterSet.Config as cms

from L1Trigger.TrackerDTC.ProducerES_cff import TrackTriggerSetup
from L1Trigger.TrackerTFP.Producer_cfi import TrackerTFPProducer_params
from L1Trigger.TrackerTFP.ProducerES_cff import TrackTriggerDataFormats
from L1Trigger.TrackerTFP.ProducerLayerEncoding_cff import TrackTriggerLayerEncoding
from L1Trigger.TrackerTFP.KalmanFilterFormats_cff import TrackTriggerKalmanFilterFormats
from L1Trigger.TrackFindingTracklet.Producer_cfi import TrackFindingTrackletProducer_params

TrackFindingTrackletProducerKFin = cms.EDProducer( 'trackFindingTracklet::ProducerKFin', TrackFindingTrackletProducer_params )
TrackFindingTrackletProducerKF = cms.EDProducer( 'trackerTFP::ProducerKF', TrackFindingTrackletProducer_params )
TrackFindingTrackletProducerKFout = cms.EDProducer( 'trackFindingTracklet::ProducerKFout', TrackFindingTrackletProducer_params )
TrackFindingTrackletProducerAS = cms.EDProducer( 'trackerTFP::ProducerAS', TrackFindingTrackletProducer_params )
#ifndef L1Trigger_TrackFindingTracklet_interface_HybridFit_h
#define L1Trigger_TrackFindingTracklet_interface_HybridFit_h

#ifdef USEHYBRID
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "L1Trigger/TrackFindingTMTT/interface/L1track3D.h"
#include "L1Trigger/TrackFindingTMTT/interface/Stub.h"
#include "L1Trigger/TrackFindingTMTT/interface/KFParamsComb.h"
#include "L1Trigger/TrackFindingTracklet/interface/HybridFit.h"
#include "L1Trigger/TrackFindingTMTT/interface/Settings.h"
#include "L1Trigger/TrackFindingTMTT/interface/L1fittedTrack.h"
#include "L1Trigger/TrackFindingTMTT/interface/KFTrackletTrack.h"
#endif

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <vector>

namespace Trklet {

  class Settings;

  class HybridFit {
  public:
    HybridFit(unsigned int iSector, const Settings* const settings);

    void Fit(Tracklet* tracklet, std::vector<std::pair<Stub*, L1TStub*>>& trackstublist);

  private:
    unsigned int iSector_;

    const Settings* const settings_;
  };
};  // namespace Trklet
#endif

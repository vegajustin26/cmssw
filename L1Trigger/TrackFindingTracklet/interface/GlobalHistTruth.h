#ifndef L1Trigger_TrackFindingTracklet_interface_GlobalHistTruth_h
#define L1Trigger_TrackFindingTracklet_interface_GlobalHistTruth_h

#include "HistBase.h"

using namespace std;

class GlobalHistTruth {
public:
  static SLHCEvent*& event() {
    static SLHCEvent* theEvent = 0;
    return theEvent;
  }

  static HistBase*& histograms() {
    static HistBase dummy;
    static HistBase* theHistBase = &dummy;
    return theHistBase;
  }
};

#endif

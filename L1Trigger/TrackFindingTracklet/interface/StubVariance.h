#ifndef L1Trigger_TrackFindingTracklet_interface_StubVariance_h
#define L1Trigger_TrackFindingTracklet_interface_StubVariance_h

class GlobalHistTruth;
class SLHCEvent;

namespace Trklet {
  
  class StubVariance {
  public:
    StubVariance(SLHCEvent& ev, GlobalHistTruth* globals);
    
    void process(SLHCEvent& ev, GlobalHistTruth* globals);
    
  private:
    
  };
};
#endif

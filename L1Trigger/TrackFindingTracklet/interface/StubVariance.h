#ifndef L1Trigger_TrackFindingTracklet_interface_StubVariance_h
#define L1Trigger_TrackFindingTracklet_interface_StubVariance_h

namespace Trklet {

  class Globals;
  class SLHCEvent;
  
  class StubVariance {
  public:
    StubVariance(SLHCEvent& ev, Globals* globals);
    
    void process(SLHCEvent& ev, Globals* globals);
    
  private:
    
  };
};
#endif

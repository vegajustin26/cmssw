#ifndef L1Trigger_TrackFindingTracklet_interface_StubVariance_h
#define L1Trigger_TrackFindingTracklet_interface_StubVariance_h

namespace trklet {

  class Globals;
  class SLHCEvent;

  class StubVariance {
  public:
    StubVariance(SLHCEvent& ev, Globals* globals);

    ~StubVariance() {}

    void process(SLHCEvent& ev, Globals* globals);

  private:
  };
};  // namespace trklet
#endif

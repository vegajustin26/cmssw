#ifndef L1Trigger_TrackFindingTracklet_interface_StubVariance_h
#define L1Trigger_TrackFindingTracklet_interface_StubVariance_h

namespace trklet {

  class Globals;
  class Settings;
  class SLHCEvent;

  class StubVariance {
  public:
    StubVariance(SLHCEvent& ev, Globals* globals, const Settings* settings);

    ~StubVariance() = default;

    void process(SLHCEvent& ev, Globals* globals, const Settings* settings);

  private:
  };
};  // namespace trklet
#endif

#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletEventProcessor_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletEventProcessor_h

#include "L1Trigger/TrackFindingTracklet/interface/CPUTimer.h"
#include "L1Trigger/TrackFindingTracklet/interface/Cabling.h"

#include <map>
#include <vector>
#include <string>

namespace Trklet {

  class Settings;
  class SLHCEvent;
  class Globals;
  class Sector;
  class HistImp;
  class Track;

  class TrackletEventProcessor {
  public:
    TrackletEventProcessor(){};

    virtual ~TrackletEventProcessor();

    void init(const Settings* theSettings);

    void event(SLHCEvent& ev);

    void printSummary();

    std::vector<Track*>& tracks() { return tracks_; }

  private:
    const Settings* settings_{0};

    Globals* globals_{};

    Sector** sectors_{};

    HistImp* histimp_{};

    int eventnum_ = {0};

    Cabling cabling_;

    CPUTimer cleanTimer_;
    CPUTimer addStubTimer_;
    CPUTimer VMRouterTimer_;
    CPUTimer TETimer_;
    CPUTimer TEDTimer_;
    CPUTimer TRETimer_;
    CPUTimer TCTimer_;
    CPUTimer TCDTimer_;
    CPUTimer PRTimer_;
    CPUTimer METimer_;
    CPUTimer MCTimer_;
    CPUTimer MPTimer_;
    CPUTimer FTTimer_;
    CPUTimer PDTimer_;

    std::vector<Track*> tracks_;

    std::map<std::string, std::vector<int> > dtclayerdisk_;
  };

};  // namespace Trklet
#endif

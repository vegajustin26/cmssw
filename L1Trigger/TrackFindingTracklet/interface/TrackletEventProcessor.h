#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletEventProcessor_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletEventProcessor_h

#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"
#include "L1Trigger/TrackFindingTracklet/interface/CPUTimer.h"
#include "L1Trigger/TrackFindingTracklet/interface/Cabling.h"
#include "L1Trigger/TrackFindingTracklet/interface/slhcevent.h"
#include "L1Trigger/TrackFindingTracklet/interface/Sector.h"


namespace Trklet {

  class TrackletEventProcessor {
    
  public:
    
    TrackletEventProcessor() { };

    ~TrackletEventProcessor();

    void init(const Settings* theSettings);

    void event(SLHCEvent& ev);

    void printSummary();

    std::vector<Track*>& tracks() { return tracks_; }

  private:
    
    const Settings* settings_{0};
    
    Globals* globals_{};
    
    Sector** sectors_{};
    
    int eventnum_={0};
    
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

    std::map<string,vector<int> > dtclayerdisk_;
    
  };
  
};    
#endif
      

#ifndef L1Trigger_TrackFindingTracklet_interface_MatchEngine_h
#define L1Trigger_TrackFindingTracklet_interface_MatchEngine_h

#include "L1Trigger/TrackFindingTracklet/interface/ProcessBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletCalculator.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubsMEMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMProjectionsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/CandidateMatchMemory.h"

#include <vector>

class MemoryBase;

namespace Trklet {
  
  class Settings;
  class Globals;

  class MatchEngine:public ProcessBase{
    
  public:
    
    MatchEngine(string name, const Settings* settings, Globals* global, unsigned int iSector);
    
    void addOutput(MemoryBase* memory,string output);

    void addInput(MemoryBase* memory,string input);
    
    void execute();

    double bend(double r, double rinv);
  
  private:
    
    VMStubsMEMemory* vmstubs_;
    VMProjectionsMemory* vmprojs_;
    
    CandidateMatchMemory* candmatches_;
    
    int layer_;
    int disk_;
    
    //used in the layers
    vector<bool> table_;
    
    //used in the disks
    vector<bool> tablePS_;
    vector<bool> table2S_;
  };
  
};
#endif

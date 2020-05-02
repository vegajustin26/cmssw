#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletCalculator_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletCalculator_h

#include "L1Trigger/TrackFindingTracklet/interface/TrackletCalculatorBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletProjectionsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/StubPairsMemory.h"

#include <vector>

namespace Trklet {

  class Settings;
  class Globals;
  
  class TrackletCalculator:public TrackletCalculatorBase{
    
  public:
    
    TrackletCalculator(std::string name, const Settings* const settings, Globals* globals, unsigned int iSector);
    
    void addOutputProjection(TrackletProjectionsMemory* &outputProj, MemoryBase* memory);
    
    void addOutput(MemoryBase* memory,std::string output);
    
    void addInput(MemoryBase* memory,std::string input);
    
    void execute();
    
  private:
    
    int iTC_;
    
    unsigned int maxtracklet_; //maximum numbor of tracklets that be stored

    std::vector<AllStubsMemory*> innerallstubs_;
    std::vector<AllStubsMemory*> outerallstubs_;
    std::vector<StubPairsMemory*> stubpairs_;
    
  };
  
};



#endif

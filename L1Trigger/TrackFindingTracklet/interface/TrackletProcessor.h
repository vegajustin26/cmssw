#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletProcessor_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletProcessor_h

#include "L1Trigger/TrackFindingTracklet/interface/TrackletCalculatorBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubsTEMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/StubPairsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletProjectionsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"

#include <vector>

class MemoryBase;

namespace Trklet {

  class Settings;
  class Globals;

  class TrackletProcessor:public TrackletCalculatorBase{
    
  public:
    
    TrackletProcessor(string name, const Settings* const settings,Globals* globals, unsigned int iSector);
    
    void addOutputProjection(TrackletProjectionsMemory* &outputProj, MemoryBase* memory);
  
    void addOutput(MemoryBase* memory,string output);
    
    void addInput(MemoryBase* memory,string input);

    void execute();

    void setVMPhiBin();
    
    double rinv(double phi1, double phi2,double r1, double r2);

    double bend(double r, double rinv);

    void writeTETable();
    
  private:
    
    int iTC_;
    
    unsigned int maxtracklet_; //maximum numbor of tracklets that be stored
    
    std::vector<VMStubsTEMemory*> innervmstubs_;
    std::vector<VMStubsTEMemory*> outervmstubs_;
    
    std::vector<AllStubsMemory*> innerallstubs_;
    std::vector<AllStubsMemory*> outerallstubs_;
    
    bool extra_;
    
    std::map<unsigned int, std::vector<bool> > phitable_;
    std::map<unsigned int, std::vector<bool> > pttableinner_;
    std::map<unsigned int, std::vector<bool> > pttableouter_;
    
    int innerphibits_;
    int outerphibits_;
    
  };

};
#endif

#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletEngineDisplaced_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletEngineDisplaced_h

#include "L1Trigger/TrackFindingTracklet/interface/ProcessBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubsTEMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/StubPairsMemory.h"

#include <vector>

class MemoryBase;

namespace Trklet {

  class Settings;
  class Globals;

  class TrackletEngineDisplaced:public ProcessBase{
    
  public:
    
    TrackletEngineDisplaced(string name, const Settings* settings, Globals* global, unsigned int iSector);

    ~TrackletEngineDisplaced();
    
    void addOutput(MemoryBase* memory,string output);
    
    void addInput(MemoryBase* memory,string input);
    
    void execute();
    
    void readTables();
  
  private:
    
    double phimax_;
    double phimin_;
    
    int layer1_;
    int layer2_;
    int disk1_;
    int disk2_;
    
    vector<VMStubsTEMemory*> firstvmstubs_;
    VMStubsTEMemory* secondvmstubs_;
    
    std::vector<StubPairsMemory*> stubpairs_;
    
    std::vector<set<string> > table_;
    
    int firstphibits_;
    int secondphibits_;
    
    int iSeed_;
    
  };

};
#endif

#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletEngineDisplaced_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletEngineDisplaced_h

#include "L1Trigger/TrackFindingTracklet/interface/ProcessBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubsTEMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/StubPairsMemory.h"

#include <vector>

namespace Trklet {

  class Settings;
  class Globals;
  class MemoryBase;

  class TrackletEngineDisplaced : public ProcessBase {
    
  public:
    
    TrackletEngineDisplaced(std::string name, const Settings* settings, Globals* global, unsigned int iSector);

    ~TrackletEngineDisplaced();
    
    void addOutput(MemoryBase* memory,std::string output);
    void addInput(MemoryBase* memory,std::string input);
    
    void execute();
    
    void readTables();
  
  private:
    
    double phimax_;
    double phimin_;
    
    int layer1_;
    int layer2_;
    int disk1_;
    int disk2_;
    
    std::vector<VMStubsTEMemory*> firstvmstubs_;
    VMStubsTEMemory* secondvmstubs_;
    
    std::vector<StubPairsMemory*> stubpairs_;
    
    std::vector<std::set<std::string> > table_;
    
    int firstphibits_;
    int secondphibits_;
    
    int iSeed_;
    
  };
};
#endif

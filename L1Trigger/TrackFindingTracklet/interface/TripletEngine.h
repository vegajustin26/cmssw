#ifndef L1Trigger_TrackFindingTracklet_interface_TripletEngine_h
#define L1Trigger_TrackFindingTracklet_interface_TripletEngine_h

#include "L1Trigger/TrackFindingTracklet/interface/ProcessBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/StubTripletsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubsTEMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/StubPairsMemory.h"

#include <vector>

namespace Trklet {

  class Settings;
  class Globals;

  class TripletEngine:public ProcessBase{
    
  public:
    
    TripletEngine(std::string name, const Settings* settings, Globals* global, unsigned int iSector);

    ~TripletEngine();

    void addOutput(MemoryBase* memory,std::string output);

    void addInput(MemoryBase* memory,std::string input);
    
    void execute();

    void readTables();
    void writeTables();

  private:

    int count_;
    
    double phimax_;
    double phimin_;
    
    int layer1_;
    int layer2_;
    int layer3_;
    int disk1_;
    int disk2_;
    int disk3_;
    int dct1_;
    int dct2_;
    int dct3_;
    int phi1_;
    int phi2_;
    int phi3_;
    int z1_;
    int z2_;
    int z3_;
    int r1_;
    int r2_;
    int r3_;
  
    std::vector<VMStubsTEMemory*> thirdvmstubs_;
    std::vector<StubPairsMemory*> stubpairs_;
    
    StubTripletsMemory* stubtriplets_;

    std::map<std::string, std::vector<std::vector<std::string> > > tmpSPTable_;
    std::map<std::string, std::vector<std::map<std::string, unsigned> > > spTable_;
    std::vector<bool> table_;

    int secondphibits_;
    int thirdphibits_;
    
    int iSeed_;
  };

};
#endif

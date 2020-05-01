#ifndef L1Trigger_TrackFindingTracklet_interface_PurgeDuplicate_h
#define L1Trigger_TrackFindingTracklet_interface_PurgeDuplicate_h

//#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackFitMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/CleanTrackMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/ProcessBase.h"

#include <vector>

class GlobalHistTruth;
class MemoryBase;
class Tracklet;
class L1TStub;

namespace Trklet {

  class Settings;
  class Stub;
  class Track;

  class PurgeDuplicate:public ProcessBase{
    
  public:
    
    PurgeDuplicate(string name, const Settings* settings, GlobalHistTruth* global, unsigned int iSector);

    void addOutput(MemoryBase* memory,string output);

    void addInput(MemoryBase* memory,string input);
    
    void execute(std::vector<Track*>& outputtracks_);

  private:
    
    double getPhiRes(Tracklet* curTracklet, std::pair<Stub*, L1TStub*> curStub);
    
    std::vector<Track*> inputtracks_;
    std::vector<std::vector<std::pair<Stub*,L1TStub*>>> inputstublists_;
    std::vector<std::vector<std::pair<int,int>>> inputstubidslists_;
    std::vector<std::vector<std::pair<int,int>>> mergedstubidslists_;
    std::vector<TrackFitMemory*> inputtrackfits_;
    std::vector<Tracklet*> inputtracklets_;
    std::vector<CleanTrackMemory*> outputtracklets_;
    
  };

};
#endif

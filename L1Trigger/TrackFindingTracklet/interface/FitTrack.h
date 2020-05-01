#ifndef L1Trigger_TrackFindingTracklet_interface_FitTrack_H
#define L1Trigger_TrackFindingTracklet_interface_FitTrack_H

#include "L1Trigger/TrackFindingTracklet/interface/ProcessBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletParametersMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/FullMatchMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackFitMemory.h"

#include <vector>

class L1TStub;
class Stub;
class GlobalHistTruth;


namespace Trklet {

  class Settings;
  
  class FitTrack:public ProcessBase{
    
  public:
    
    FitTrack(string name, const Settings* settings, GlobalHistTruth* global, unsigned int iSector);
    
    void addOutput(MemoryBase* memory,string output);

    void addInput(MemoryBase* memory,string input);
    
    // used if USEHYBRID is not defined
    void trackFitChisq(Tracklet* tracklet, std::vector<std::pair<Stub*,L1TStub*>> &, std::vector<std::pair<int,int>> &);
    
    // used if USEHYBRID is defined
    void trackFitKF(Tracklet* tracklet, std::vector<std::pair<Stub*,L1TStub*>> &trackstublist, std::vector<std::pair<int,int>> &stubidslist);

    // used for propagating tracklet without fitting
    void trackFitFake(Tracklet* tracklet, std::vector<std::pair<Stub*,L1TStub*>> &, std::vector<std::pair<int,int>> &);
    
    std::vector<Tracklet*> orderedMatches(vector<FullMatchMemory*>& fullmatch);

    void execute();

  private:
    
    std::vector<TrackletParametersMemory*> seedtracklet_;
    std::vector<FullMatchMemory*> fullmatch1_;
    std::vector<FullMatchMemory*> fullmatch2_;
    std::vector<FullMatchMemory*> fullmatch3_;
    std::vector<FullMatchMemory*> fullmatch4_;
    
    TrackFitMemory* trackfit_;
    
  };

};
#endif

#ifndef L1Trigger_TrackFindingTracklet_interface_FullMatchMemory_h
#define L1Trigger_TrackFindingTracklet_interface_FullMatchMemory_h

#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"

#include <vector>

class Tracklet;
class L1TStub;
class Stub;

namespace Trklet {

  class Settings;

  
  class FullMatchMemory : public MemoryBase {
  public:
    FullMatchMemory(string name, const Settings* const settings, unsigned int iSector, double phimin, double phimax);

    void addMatch(Tracklet* tracklet, std::pair<Stub*, L1TStub*> stub);

    void addMatch(std::pair<Tracklet*, std::pair<Stub*, L1TStub*> > match);
    
    unsigned int nMatches() const { return matches_.size(); }
    
    Tracklet* getFPGATracklet(unsigned int i) const { return matches_[i].first; }
    
    std::pair<Tracklet*, std::pair<Stub*, L1TStub*> > getMatch(unsigned int i) const { return matches_[i]; }
    
    void clean() { matches_.clear(); }
    
    void writeMC(bool first);
    
    int layer() const { return layer_; }
    int disk() const { return disk_; }
    
  private:
    double phimin_;
    double phimax_;
    std::vector<std::pair<Tracklet*, std::pair<Stub*, L1TStub*> > > matches_;
    
    int layer_;
    int disk_;
  };

};
#endif

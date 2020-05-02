#ifndef L1Trigger_TrackFindingTracklet_interface_AllProjectionsMemory_h
#define L1Trigger_TrackFindingTracklet_interface_AllProjectionsMemory_h

#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"

#include <vector>


namespace Trklet {

  class Settings;
  
  class AllProjectionsMemory : public MemoryBase {
  public:
    AllProjectionsMemory(std::string name, const Settings* const settings, unsigned int iSector, double phimin, double phimax);
    
    void addTracklet(Tracklet* tracklet) { tracklets_.push_back(tracklet); }
    
    unsigned int nTracklets() const { return tracklets_.size(); }
    
    Tracklet* getFPGATracklet(unsigned int i) const { return tracklets_[i]; }
    
    void clean() { tracklets_.clear(); }
    
    void writeAP(bool first);
    
    //int layer() const { return layer_; }
    //int disk() const { return disk_; }
    
    
  private:
    double phimin_;
    double phimax_;
    std::vector<Tracklet*> tracklets_;
    
    int layer_;
    int disk_;
    
  };
  
};

#endif

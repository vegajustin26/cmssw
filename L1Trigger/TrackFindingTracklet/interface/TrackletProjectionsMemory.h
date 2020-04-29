#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletProjectionsMemory_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletProjectionsMemory_h

#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"

#include <vector>

class Tracklet;

namespace Trklet {

  class Settings;

  
  class TrackletProjectionsMemory : public MemoryBase {
  public:
    TrackletProjectionsMemory(string name, const Settings* const settings, unsigned int iSector, double phimin, double phimax);

    void addProj(Tracklet* tracklet);

    unsigned int nTracklets() const { return tracklets_.size(); }

    Tracklet* getFPGATracklet(unsigned int i) const { return tracklets_[i]; }

    void clean();

    void writeTPROJ(bool first);

    int layer() const { return layer_; }
    int disk() const { return disk_; }
    
  private:
    double phimin_;
    double phimax_;
    std::vector<Tracklet*> tracklets_;
    
    int layer_;
    int disk_;
  };

};
#endif

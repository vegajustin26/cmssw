// This class holds the tracklet parameters for the selected stub pairs
// This class owns the tracklets. Further modules only holds pointers
#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletParametersMemory_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletParametersMemory_h

#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"


#include <vector>


namespace Trklet {

  class Settings;
  class Globals;

  
  class TrackletParametersMemory : public MemoryBase {
  public:
    TrackletParametersMemory(string name, const Settings* const settings, unsigned int iSector, double phimin, double phimax);
    
    void addTracklet(Tracklet *tracklet) { tracklets_.push_back(tracklet); }
    
    unsigned int nTracklets() const { return tracklets_.size(); }
    
    Tracklet *getFPGATracklet(unsigned int i) const { return tracklets_[i]; }
    
    void clean() {
      for (unsigned int i = 0; i < tracklets_.size(); i++) {
	delete tracklets_[i];
      }
      tracklets_.clear();
    }
    
    void writeMatches(Globals* globals, int &matchesL1, int &matchesL3, int &matchesL5);

    void writeTPAR(bool first);
    
  private:
    double phimin_;
    double phimax_;
    std::vector<Tracklet *> tracklets_;
  };

};
#endif

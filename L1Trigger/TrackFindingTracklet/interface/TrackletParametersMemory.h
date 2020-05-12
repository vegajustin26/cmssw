// This class holds the tracklet parameters for the selected stub pairs
// This class owns the tracklets. Further modules only holds pointers
#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletParametersMemory_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletParametersMemory_h

#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"

#include <string>
#include <vector>

namespace trklet {

  class Settings;
  class Globals;
  class Tracklet;

  class TrackletParametersMemory : public MemoryBase {
  public:
    TrackletParametersMemory(std::string name, const Settings *const settings, unsigned int iSector);

    void addTracklet(Tracklet *tracklet) { tracklets_.push_back(tracklet); }

    unsigned int nTracklets() const { return tracklets_.size(); }

    Tracklet *getFPGATracklet(unsigned int i) const { return tracklets_[i]; }

    void clean();

    void writeMatches(Globals *globals, int &matchesL1, int &matchesL3, int &matchesL5);

    void writeTPAR(bool first);

  private:
    std::vector<Tracklet *> tracklets_;
  };

};  // namespace trklet
#endif

#ifndef L1Trigger_TrackFindingTracklet_interface_TETableOuterDisk_h
#define L1Trigger_TrackFindingTracklet_interface_TETableOuterDisk_h

#include "L1Trigger/TrackFindingTracklet/interface/TETableBase.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <vector>

namespace Trklet {

  class Settings;

  class TETableOuterDisk : public TETableBase {
  public:
    TETableOuterDisk(const Settings* settings);
    TETableOuterDisk(const Settings* settings, int disk, int zbits, int rbits);

    ~TETableOuterDisk() {}

    void init(const Settings* settings, int disk, int zbits, int rbits);

    // negative return means that seed can not be formed
    int getLookupValue(int irbin, int izbin);

    int lookup(int zbin, int rbin);

  private:
    double zmean_;

    double rmin_;
    double rmax_;

    double zmin_;
    double zmax_;

    double dr_;
    double dz_;

    int zbits_;
    int rbits_;

    int zbins_;
    int rbins_;

    int disk_;
  };
};  // namespace Trklet
#endif

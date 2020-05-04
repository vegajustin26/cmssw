#ifndef L1Trigger_TrackFindingTracklet_interface_TETableOuter_h
#define L1Trigger_TrackFindingTracklet_interface_TETableOuter_h

#include "L1Trigger/TrackFindingTracklet/interface/TETableBase.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <vector>

namespace Trklet {

  class Settings;

  class TETableOuter : public TETableBase {
  public:
    TETableOuter(const Settings* settings);
    TETableOuter(const Settings* settings, int layer, int zbits, int rbits);

    ~TETableOuter() {}

    void init(const Settings* settings, int layer, int zbits, int rbits);

    // negative return means that seed can not be formed
    int getLookupValue(const Settings* settings, int izbin, int irbin);

    int lookup(int zbin, int rbin);

  private:
    double rmean_;

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

    int layer_;
  };
};  // namespace Trklet
#endif

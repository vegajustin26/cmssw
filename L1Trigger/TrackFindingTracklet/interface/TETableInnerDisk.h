#ifndef L1Trigger_TrackFindingTracklet_interface_TETableInnerDisk_h
#define L1Trigger_TrackFindingTracklet_interface_TETableInnerDisk_h

#include "L1Trigger/TrackFindingTracklet/interface/TETableBase.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <vector>

namespace Trklet {

  class Settings;

  class TETableInnerDisk : public TETableBase {
  public:
    TETableInnerDisk(const Settings* settings);
    TETableInnerDisk(const Settings* settings, int disk1, int disk2, int layer3, int zbits, int rbits);

    ~TETableInnerDisk() {}

    void init(const Settings* settings, int disk1, int disk2, int zbits, int rbits);
    void init(const Settings* settings, int disk1, int disk2, int layer3, int zbits, int rbits);

    // negative return means that seed can not be formed
    int getLookupValue(const Settings* settings, int irbin, int izbin);

    void findr(double r, double z, double& rmind2, double& rmaxd2);
    double rintercept(double zcut, double r, double z);

    void findz(double z, double r, double& zminl3, double& zmaxl3);
    double zintercept(double zcut, double z, double r);

    int lookup(int zbin, int rbin);

  private:
    int disk1_;
    int disk2_;
    int layer3_;
    int zbits_;
    int rbits_;

    int rbins_;
    double rmind1_;
    double rmaxd1_;
    double dr_;

    int zbins_;
    double zmind1_;
    double zmaxd1_;
    double dz_;

    double zmeand2_;
    double rmeanl3_;
  };
};  // namespace Trklet
#endif

#ifndef L1Trigger_TrackFindingTracklet_interface_TETableInner_h
#define L1Trigger_TrackFindingTracklet_interface_TETableInner_h

#include "L1Trigger/TrackFindingTracklet/interface/TETableBase.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <vector>

namespace Trklet {
  
  class Settings;
  
  class TETableInner : public TETableBase {
  public:

    TETableInner(const Settings* settings);
    TETableInner(const Settings* settings, int layer1, int layer2, int layer3, int zbits, int rbits, bool thirdLayerIsDisk = false);
    
    ~TETableInner() {}

    void init(const Settings* settings, int layer1, int layer2, int zbits, int rbits);
    void init(const Settings* settings, int layer1, int layer2, int layer3, int zbits, int rbits, bool thirdLayerIsDisk = false);
    
    // negative return means that seed can not be formed
    int getLookupValue(const Settings* settings, int izbin, int irbin, bool extra);
    
    void findzL2(double z, double r, double& zminl2, double& zmaxl2);
    double zinterceptL2(double zcut, double z, double r);

    void findzL3(double z, double r, double& zminl3, double& zmaxl3);
    double zinterceptL3(double zcut, double z, double r);

    void findr(double r, double z, double& rmind2, double& rmaxd2);
    double rintercept(double zcut, double r, double z);

    int lookup(int zbin, int rbin);
    
  private:
    bool thirdLayerIsDisk_;
    
    int layer1_;
    int layer2_;
    int layer3_;
    int zbits_;
    int rbits_;
    
    int rbins_;
    double rminl1_;
    double rmaxl1_;
    double dr_;
    
    int zbins_;
    double zminl1_;
    double zminl2_;
    double dz_;
    
    double rmeanl2_;
    double rmeanl3_;
    double zmeand3_;
    
    double rmaxdisk_;
    double rmindisk_;
  };
  
};
#endif

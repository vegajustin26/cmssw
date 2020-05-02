#ifndef L1Trigger_TrackFindingTracklet_interface_TETableInnerOverlap_h
#define L1Trigger_TrackFindingTracklet_interface_TETableInnerOverlap_h

#include "L1Trigger/TrackFindingTracklet/interface/TETableBase.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <vector>

namespace Trklet {
  
  class Settings;
  
  class TETableInnerOverlap : public TETableBase{
    
  public:
    
    TETableInnerOverlap(const Settings* settings);
    TETableInnerOverlap(const Settings* settings, int layer1, int disk2, int zbits, int rbits);

    ~TETableInnerOverlap() { }

    void init(const Settings* settings, int layer1, int disk2, int zbits, int rbits);
    
    // negative return means that seed can not be formed
    int getLookupValue(const Settings* settings, int izbin, int irbin);

    void findr(double r, double z, double& rmind2, double& rmaxd2);
    double rintercept(double zcut, double r, double z);

    int lookup(int zbin, int rbin);
    
  private:
    
    int layer1_;
    int disk2_;
    int zbits_;
    int rbits_;
    
    int rbins_;
    double rminl1_;
    double rmaxl1_;
    double dr_;
    
    int zbins_;
    double zminl1_;
    double zmaxl1_;
    double dz_;
    
    double zmeand2_;
    
    double rmaxdisk_;
    double rmindisk_;
  };
};
#endif

#ifndef L1Trigger_TrackFindingTracklet_interface_ProjectionTemp_h
#define L1Trigger_TrackFindingTracklet_interface_ProjectionTemp_h

#include <cassert>
#include "Tracklet.h"

using namespace std;
using namespace Trklet;

class ProjectionTemp {
public:
  ProjectionTemp(Tracklet* proj,
                 unsigned int slot,
                 unsigned int projrinv,
                 int projfinerz,
                 unsigned int projfinephi,
                 unsigned int iphi,
                 bool isPSseed) {
    proj_ = proj;
    slot_ = slot;
    projrinv_ = projrinv;
    projfinerz_ = projfinerz;
    projfinephi_ = projfinephi;
    iphi_ = iphi;
    isPSseed_ = isPSseed;
  }

  ProjectionTemp() {
    proj_ = 0;
    slot_ = 0;
    projrinv_ = 0;
    projfinerz_ = 0;
    projfinephi_ = 0;
    iphi_ = 0;
    isPSseed_ = false;
  }

  ~ProjectionTemp() {}

  Tracklet* proj() const { return proj_; }
  unsigned int slot() const { return slot_; }
  unsigned int projrinv() const { return projrinv_; }
  int projfinerz() const { return projfinerz_; }
  unsigned int projfinephi() const { return projfinephi_; }
  unsigned int iphi() const { return iphi_; }
  bool isPSseed() const { return isPSseed_; }

private:
  Tracklet* proj_;
  unsigned int slot_;
  unsigned int projrinv_;
  unsigned int projfinerz_;
  unsigned int projfinephi_;
  unsigned int iphi_;
  bool isPSseed_;
};

#endif

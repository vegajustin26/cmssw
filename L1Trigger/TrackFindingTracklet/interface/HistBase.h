#ifndef L1Trigger_TrackFindingTracklet_interface_HistBase_h
#define L1Trigger_TrackFindingTracklet_interface_HistBase_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <bitset>
#include <cassert>
#include <cmath>

using namespace std;

class HistBase {
public:
  HistBase() {}

  virtual ~HistBase() {}

  virtual void FillLayerResidual(int, int, double, double, double, double, bool) {
    return;  //default implementation does nothing
  }

  virtual void FillDiskResidual(int, int, double, double, double, double, bool) {
    return;  //default implementation does nothing
  }

  //arguments are
  // int seedIndex
  // int iSector
  // double irinv, rinv
  // double iphi0, phi0
  // double ieta, eta
  // double iz0, z0
  // int tp
  virtual void fillTrackletParams(const Settings*, int, int, double, double, double, double, double, double, double, double, int) {
    return;  //default implementation does nothing
  }

  //int seedIndex
  //double etaTP
  //bool eff
  virtual void fillSeedEff(int, double, bool) {
    return;  //default implementation does nothing
  }

private:
};

#endif

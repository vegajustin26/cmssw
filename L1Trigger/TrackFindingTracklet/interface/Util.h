#ifndef L1Trigger_TrackFindingTracklet_interface_Util_h
#define L1Trigger_TrackFindingTracklet_interface_Util_h

namespace Trklet{

  //method return phi in the -pi to +pi range
  inline double phiRange(double phi) {
    //catch if phi is very out of range, not a number etc
    assert(std::abs(phi) < 100.0);
    while (phi < -M_PI)
      phi += 2 * M_PI;
    while (phi > M_PI)
      phi -= 2 * M_PI;
    return phi;
  }
  
};

#endif

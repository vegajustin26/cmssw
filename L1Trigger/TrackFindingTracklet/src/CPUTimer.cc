#include "L1Trigger/TrackFindingTracklet/interface/CPUTimer.h"

using namespace trklet;

void CPUTimer::start() { gettimeofday(&tstart_, 0); }

void CPUTimer::stop() {
  timeval tstop;
  gettimeofday(&tstop, 0);
  float tsec = tstop.tv_sec - tstart_.tv_sec;
  float tusec = tstop.tv_usec - tstart_.tv_usec;
  if (tusec < 0) {
    tusec += 1000000.0;
    tsec -= 1.0;
  }
  double tmp = tsec + tusec / 1000000.0;
  ttot_ += tmp;
  ttotsq_ += tmp * tmp;
  ntimes_++;
}

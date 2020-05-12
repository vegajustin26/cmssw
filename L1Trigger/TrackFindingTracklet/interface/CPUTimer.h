#ifndef L1Trigger_TrackFindingTracklet_interface_CPUTimer_h
#define L1Trigger_TrackFindingTracklet_interface_CPUTimer_h

#include <cmath>
#include <sys/time.h>

namespace Trklet {

  class CPUTimer {
  public:
    CPUTimer(){}

    ~CPUTimer(){}

    void start();
    void stop();
    unsigned int ntimes() const { return ntimes_; }
    double avgtime() const { return ttot_ / ntimes_; }
    double rms() const { return sqrt((ttot_ * ttot_ - ttotsq_)) / ntimes_; }
    double tottime() const { return ttot_; }

  private:
    unsigned int ntimes_{0};
    double ttot_{0.0};
    double ttotsq_{0.0};

    timeval tstart_;
  };
};  // namespace Trklet
#endif

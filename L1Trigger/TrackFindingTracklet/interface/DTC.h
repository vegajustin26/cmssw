#ifndef L1Trigger_TrackFindingTracklet_interface_DTC_h
#define L1Trigger_TrackFindingTracklet_interface_DTC_h

#include "L1Trigger/TrackFindingTracklet/interface/DTCLink.h"

class L1TStub;

namespace Trklet {

  class Stub;
  
  class DTC {
  public:
    DTC(std::string name = "");

    void init(std::string name);
    
    void addSec(int sector);
    
    void addphi(double phi, int layerdisk);
    
    void addLink(double phimin, double phimax);
    
    int addStub(std::pair<Stub*, L1TStub*> stub);

    unsigned int nLinks() const { return links_.size(); }

    const DTCLink& link(unsigned int i) const { return links_[i]; }
    
    void clean();

    double min(unsigned int i) const { return phimin_[i]; }
    double max(unsigned int i) const { return phimax_[i]; }

  private:
    std::string name_;
    std::vector<DTCLink> links_;
    std::vector<int> sectors_;
    
    double phimin_[11];
    double phimax_[11];
  };
};
#endif

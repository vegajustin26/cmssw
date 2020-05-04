#include "L1Trigger/TrackFindingTracklet/interface/Track.h"

using namespace std;
using namespace Trklet;

Track::Track(int irinv,
             int iphi0,
             int id0,
             int it,
             int iz0,
             int ichisqrphi,
             int ichisqrz,
             double chisqrphi,
             double chisqrz,
             int hitpattern,
             std::map<int, int> stubID,
             std::vector<L1TStub*> l1stub,
             int seed) {
  irinv_ = irinv;
  iphi0_ = iphi0;
  id0_ = id0;
  iz0_ = iz0;
  it_ = it;
  ichisqrphi_ = ichisqrphi;
  ichisqrz_ = ichisqrz;

  chisqrphi_ = chisqrphi;
  chisqrz_ = chisqrz;

  hitpattern_ = hitpattern;

  nstubs_ = l1stub.size();
  if (nstubs_ > 6)
    nstubs_ = 6;  //maximum used in fit

  stubID_ = stubID;
  l1stub_ = l1stub;

  seed_ = seed;
  duplicate_ = false;
  sector_ = -1;
}

double Track::phi0(const Settings* settings) const {
  double dphi = 2 * M_PI / settings->NSector();
  double dphiHG = 0.5 * settings->dphisectorHG() - M_PI / settings->NSector();
  double phimin = sector_ * dphi - dphiHG;
  double phimax = phimin + dphi + 2 * dphiHG;
  phimin -= M_PI / settings->NSector();
  phimax -= M_PI / settings->NSector();
  phimin = Trklet::phiRange(phimin);
  phimax = Trklet::phiRange(phimax);
  if (phimin > phimax)
    phimin -= 2 * M_PI;
  double phioffset = phimin;

  return iphi0_ * settings->kphi0pars() + phioffset;
}

#include "L1Trigger/TrackFindingTracklet/interface/DTC.h"
#include "L1Trigger/TrackFindingTracklet/interface/L1TStub.h"
#include "L1Trigger/TrackFindingTracklet/interface/Stub.h"

using namespace std;
using namespace Trklet;


DTC::DTC(string name) {
  name_ = name;
  for (unsigned int i = 0; i < 11; i++) {
    phimin_[i] = 10.0;
    phimax_[i] = -10.0;
  }
}

void DTC::init(string name) {
  name_ = name;
}

void DTC::addSec(int sector) {
  sectors_.push_back(sector);
}

void DTC::addphi(double phi, int layerdisk) {
  assert(layerdisk >= 0);
  assert(layerdisk < 11);
  if (phi < phimin_[layerdisk])
    phimin_[layerdisk] = phi;
  if (phi > phimax_[layerdisk])
    phimax_[layerdisk] = phi;
}

void DTC::addLink(double phimin, double phimax) {
  DTCLink link(phimin, phimax);
  links_.push_back(link);
}

int DTC::addStub(std::pair<Stub*, L1TStub*> stub) {
  double phi = Trklet::phiRange(stub.second->phi());
  bool overlaplayer = ((stub.second->layer() + 1) % 2 == 0);
  int added = 0;
  for (unsigned int i = 0; i < links_.size(); i++) {
    if (links_[i].inRange(phi, overlaplayer)) {
      added++;
      links_[i].addStub(stub);
    }
  }
  return added;
}

void DTC::clean() {
  for (unsigned int i = 0; i < links_.size(); i++) {
    links_[i].clean();
  }
}

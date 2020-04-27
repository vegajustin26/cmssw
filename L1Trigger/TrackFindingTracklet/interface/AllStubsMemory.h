// This class holds all the stubs in a DTC region for a give layer
#ifndef L1Trigger_TrackFindingTracklet_interface_AllStubsMemory_h
#define L1Trigger_TrackFindingTracklet_interface_AllStubsMemory_h


#include "L1Trigger/TrackFindingTracklet/interface/L1TStub.h"
#include "L1Trigger/TrackFindingTracklet/interface/Stub.h"
#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"

#include <vector>


class AllStubsMemory : public MemoryBase {
public:
  AllStubsMemory(string name, const Settings* const settings, unsigned int iSector, double phimin, double phimax);

  void addStub(std::pair<Stub*, L1TStub*> stub) { stubs_.push_back(stub); }

  unsigned int nStubs() const { return stubs_.size(); }

  Stub* getFPGAStub(unsigned int i) const { return stubs_[i].first; }
  L1TStub* getL1TStub(unsigned int i) const { return stubs_[i].second; }
  std::pair<Stub*, L1TStub*> getStub(unsigned int i) const { return stubs_[i]; }

  void clean() { stubs_.clear(); }

  void writeStubs(bool first);

  int layer() const { return layer_; }
  int disk() const { return disk_; }

private:
  double phimin_;
  double phimax_;
  std::vector<std::pair<Stub*, L1TStub*> > stubs_;

  int layer_;
  int disk_;
};

#endif

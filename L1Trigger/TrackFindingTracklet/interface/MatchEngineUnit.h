#ifndef L1Trigger_TrackFindingTracklet_interface_MatchEngineUnit_h
#define L1Trigger_TrackFindingTracklet_interface_MatchEngineUnit_h

#include <cassert>
#include <vector>

#include "VMStubsTEMemory.h"
#include "CircularBuffer.h"

using namespace std;

class MatchEngineUnit {
public:
  MatchEngineUnit(bool barrel, vector<bool> table, vector<bool> tablePS, vector<bool> table2S) : candmatches_(5) {
    idle_ = true;
    barrel_ = barrel;
    table_ = table;
    tablePS_ = tablePS;
    table2S_ = table2S;
    slot_ = 1;  //This makes it idle until initialized
  }

  ~MatchEngineUnit() {}

  void init(VMStubsMEMemory* vmstubsmemory,
            unsigned int slot,
            int projrinv,
            int projfinerz,
            int projfinephi,
            bool isPSseed,
            Tracklet* proj) {
    vmstubsmemory_ = vmstubsmemory;
    idle_ = false;
    slot_ = slot;
    istub_ = 0;
    projrinv_ = projrinv;
    projfinerz_ = projfinerz;
    projfinephi_ = projfinephi;
    isPSseed_ = isPSseed;
    proj_ = proj;
  }

  bool empty() const { return candmatches_.empty(); }

  std::pair<Tracklet*, std::pair<Stub*, L1TStub*> > read() { return candmatches_.read(); }

  std::pair<Tracklet*, std::pair<Stub*, L1TStub*> > peek() const { return candmatches_.peek(); }

  Tracklet* currentProj() const { return proj_; }

  void step() {
    if (idle())
      return;
    if (candmatches_.almostfull())
      return;

    //std::pair<Stub*,L1TStub*> stub=vmstubsmemory_->getStubBin(slot_,istub_);
    VMStubME vmstub = vmstubsmemory_->getVMStubMEBin(slot_, istub_);

    istub_++;
    if (istub_ >= vmstubsmemory_->nStubsBin(slot_))
      idle_ = true;

    bool isPSmodule = vmstub.isPSmodule();
    int stubfinerz = vmstub.finerz().value();
    int stubfinephi = vmstub.finephi().value();

    int deltaphi = stubfinephi - projfinephi_;

    bool dphicut = (abs(deltaphi) < 3) || (abs(deltaphi) > 5);  //FIXME... ugly

    if (!barrel_)
      dphicut = true;

    int nbits = isPSmodule ? 3 : 4;

    unsigned int index = (projrinv_ << nbits) + vmstub.bend().value();

    //Check if stub z position consistent
    int idrz = stubfinerz - projfinerz_;
    bool pass;

    if (barrel_) {
      if (isPSseed_) {
        pass = idrz >= -2 && idrz <= 2;
      } else {
        pass = idrz >= -5 && idrz <= 5;
      }
    } else {
      if (isPSmodule) {
        pass = idrz >= -1 && idrz <= 1;
      } else {
        pass = idrz >= -5 && idrz <= 5;
      }
    }

    //Check if stub bend and proj rinv consistent
    if (pass && dphicut) {
      if (barrel_ ? table_[index] : (isPSmodule ? tablePS_[index] : table2S_[index])) {
        std::pair<Tracklet*, std::pair<Stub*, L1TStub*> > tmp(proj_, vmstub.stub());
        candmatches_.store(tmp);
      }
    }
  }

  bool idle() const { return idle_; }

private:
  VMStubsMEMemory* vmstubsmemory_;

  //unsigned int memory slot
  unsigned int slot_;
  unsigned int istub_;

  bool barrel_;
  int projrinv_;
  int projfinerz_;
  int projfinephi_;
  bool isPSseed_;
  Tracklet* proj_;

  bool idle_;

  //used in the layers
  vector<bool> table_;

  //used in the disks
  vector<bool> tablePS_;
  vector<bool> table2S_;

  //save the candidate matches
  CircularBuffer<std::pair<Tracklet*, std::pair<Stub*, L1TStub*> > > candmatches_;
};

#endif

#include "L1Trigger/TrackFindingTracklet/interface/MatchEngineUnit.h"

using namespace std;
using namespace trklet;

MatchEngineUnit::MatchEngineUnit(bool barrel, unsigned int layerdisk, vector<bool> table)
    : candmatches_(5) {
  idle_ = true;
  barrel_ = barrel;
  table_ = table;
  slot_ = 1;  //This makes it idle until initialized
  layerdisk_=layerdisk;
}

void MatchEngineUnit::init(VMStubsMEMemory* vmstubsmemory,
                           unsigned int slot,
                           int projrinv,
                           int projfinerz,
                           int projfinephi,
			   int shift,
                           bool usesecond,
                           bool isPSseed,
                           Tracklet* proj) {
  vmstubsmemory_ = vmstubsmemory;
  idle_ = false;
  slot_ = slot;
  istub_ = 0;
  projrinv_ = projrinv;
  projfinerz_ = projfinerz;
  projfinephi_ = projfinephi;
  shift_ = shift;
  usesecond_ = usesecond;
  isPSseed_ = isPSseed;
  proj_ = proj;
}

void MatchEngineUnit::step() {
  if (idle() || candmatches_.almostfull())
    return;

  const VMStubME& vmstub = vmstubsmemory_->getVMStubMEBin(slot_, istub_);

  bool isPSmodule = vmstub.isPSmodule();
  int stubfinerz = vmstub.finerz().value();
  int stubfinephi = vmstub.finephi().value();

  int deltaphi = stubfinephi - projfinephi_ + (1<<NFINERZBITS)*shift_;

  bool dphicut = (abs(deltaphi) < 3);

  int nbits = isPSmodule ? 3 : 4;

  int diskps = (!barrel_)&&isPSmodule;
  
  unsigned int index = (diskps<<(4+5)) + (projrinv_ << nbits) + vmstub.bend().value();

  //Check if stub z position consistent
  int idrz = stubfinerz - projfinerz_;
  bool pass;

  if (barrel_) {
    if (isPSseed_) {
      pass = idrz >= -1 && idrz <= 1;
    } else {
      pass = idrz >= -5 && idrz <= 5;
    }
  } else {
    if (isPSmodule) {
      pass = idrz >= -1 && idrz <= 1;
    } else {
      pass = idrz >= -3 && idrz <= 3;
    }
  }

  //Check if stub bend and proj rinv consistent
  if ((pass && dphicut) && table_[index] ) {
    std::pair<Tracklet*, const Stub*> tmp(proj_, vmstub.stub());
    candmatches_.store(tmp);
  }

  istub_++;
  if (istub_ >= vmstubsmemory_->nStubsBin(slot_)) {
    if (usesecond_) {
      usesecond_ = false;
      istub_ = 0;
      slot_++;
      projfinerz_ -= 8;
    } else {
      idle_ = true;
    }
  }
}

void MatchEngineUnit::reset() {
  candmatches_.reset();
  idle_ = true;
  istub_ = 0;
}

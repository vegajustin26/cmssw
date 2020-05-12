#include "L1Trigger/TrackFindingTracklet/interface/FullMatchMemory.h"

#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"
#include "L1Trigger/TrackFindingTracklet/interface/Stub.h"
#include "L1Trigger/TrackFindingTracklet/interface/L1TStub.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;
using namespace trklet;

FullMatchMemory::FullMatchMemory(string name, const Settings* const settings, unsigned int iSector)
    : MemoryBase(name, settings, iSector) {
  if (settings_->extended()) {
    initLayerDisk(10, layer_, disk_);
  } else {
    initLayerDisk(8, layer_, disk_);
  }
}

void FullMatchMemory::addMatch(Tracklet* tracklet, std::pair<const Stub*, const L1TStub*> stub) {
  if (!settings_->doKF()) {  //When using KF we allow multiple matches
    for (unsigned int i = 0; i < matches_.size(); i++) {
      if (matches_[i].first == tracklet) {  //Better match, replace
        matches_[i].second = stub;
        return;
      }
    }
  }
  std::pair<Tracklet*, std::pair<const Stub*, const L1TStub*> > tmp(tracklet, stub);
  //Check that we have the right TCID order
  if (matches_.size() > 0) {
    if ((!settings_->doKF() && matches_[matches_.size() - 1].first->TCID() >= tracklet->TCID()) ||
        (settings_->doKF() && matches_[matches_.size() - 1].first->TCID() > tracklet->TCID())) {
      edm::LogPrint("Tracklet") << "Wrong TCID ordering in " << getName() << " : "
                                << matches_[matches_.size() - 1].first->TCID() << " " << tracklet->TCID() << " "
                                << matches_[matches_.size() - 1].first->trackletIndex() << " "
                                << tracklet->trackletIndex();
    }
  }
  matches_.push_back(tmp);
}

void FullMatchMemory::writeMC(bool first) {
  std::string fname = "../data/MemPrints/Matches/FullMatches_";
  fname += getName();
  fname += "_";
  if (iSector_ + 1 < 10)
    fname += "0";
  fname += std::to_string(iSector_ + 1);
  fname += ".dat";
  if (first) {
    bx_ = 0;
    event_ = 1;
    out_.open(fname.c_str());
  } else
    out_.open(fname.c_str(), std::ofstream::app);

  out_ << "BX = " << (bitset<3>)bx_ << " Event : " << event_ << endl;

  for (unsigned int j = 0; j < matches_.size(); j++) {
    string match = (layer_ > 0) ? matches_[j].first->fullmatchstr(layer_) : matches_[j].first->fullmatchdiskstr(disk_);
    out_ << "0x";
    if (j < 16)
      out_ << "0";
    out_ << hex << j << dec;
    out_ << " " << match << " " << trklet::hexFormat(match) << endl;
  }
  out_.close();

  bx_++;
  event_++;
  if (bx_ > 7)
    bx_ = 0;
}

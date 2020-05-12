#include "L1Trigger/TrackFindingTracklet/interface/TrackletProjectionsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;
using namespace Trklet;

TrackletProjectionsMemory::TrackletProjectionsMemory(string name, const Settings* const settings, unsigned int iSector)
    : MemoryBase(name, settings, iSector) {
  if (settings_->extended()) {
    initLayerDisk(14, layer_, disk_);
  } else {
    initLayerDisk(12, layer_, disk_);
  }
}

void TrackletProjectionsMemory::addProj(Tracklet* tracklet) {
  if (layer_ != 0 && disk_ == 0)
    assert(tracklet->validProj(layer_));
  if (layer_ == 0 && disk_ != 0)
    assert(tracklet->validProjDisk(disk_));
  if (layer_ != 0 && disk_ != 0)
    assert(tracklet->validProj(layer_) || tracklet->validProjDisk(disk_));

  for (unsigned int i = 0; i < tracklets_.size(); i++) {
    if (tracklets_[i] == tracklet) {
      edm::LogPrint("Tracklet") << "Adding same tracklet " << tracklet << " twice in " << getName();
    }
    assert(tracklets_[i] != tracklet);
  }

  tracklets_.push_back(tracklet);
}

void TrackletProjectionsMemory::clean() { tracklets_.clear(); }

void TrackletProjectionsMemory::writeTPROJ(bool first) {
  std::string fname = "../data/MemPrints/TrackletProjections/TrackletProjections_";
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

  for (unsigned int j = 0; j < tracklets_.size(); j++) {
    string proj = (layer_ > 0 && tracklets_[j]->validProj(layer_)) ? tracklets_[j]->trackletprojstrlayer(layer_)
                                                                   : tracklets_[j]->trackletprojstrdisk(disk_);
    out_ << "0x";
    if (j < 16)
      out_ << "0";
    out_ << hex << j << dec;
    out_ << " " << proj << "  " << Trklet::hexFormat(proj) << endl;
  }
  out_.close();

  bx_++;
  event_++;
  if (bx_ > 7)
    bx_ = 0;
}

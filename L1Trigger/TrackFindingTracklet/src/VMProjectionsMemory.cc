#include "L1Trigger/TrackFindingTracklet/interface/VMProjectionsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;
using namespace Trklet;

VMProjectionsMemory::VMProjectionsMemory(string name, const Settings* const settings, unsigned int iSector)
    : MemoryBase(name, settings, iSector) {

  initLayerDisk(7,layer_,disk_);

}

void VMProjectionsMemory::addTracklet(Tracklet* tracklet, unsigned int allprojindex) {
  std::pair<Tracklet*, unsigned int> tmp(tracklet, allprojindex);
  //Check that order of TCID is correct
  if (tracklets_.size() > 0) {
    assert(tracklets_[tracklets_.size() - 1].first->TCID() <= tracklet->TCID());
  }
  tracklets_.push_back(tmp);
}

void VMProjectionsMemory::writeVMPROJ(bool first) {
  std::string fname = "../data/MemPrints/VMProjections/VMProjections_";
  fname += getName();
  //get rid of duplicates
  int len = fname.size();
  if (fname[len - 2] == 'n' && fname[len - 1] > '1' && fname[len - 1] <= '9')
    return;

  fname += "_";
  ostringstream oss;
  oss << iSector_ + 1;
  if (iSector_ + 1 < 10)
    fname += "0";
  fname += oss.str();
  fname += ".dat";
  if (first) {
    bx_ = 0;
    event_ = 1;
    out_.open(fname.c_str());
  } else
    out_.open(fname.c_str(), std::ofstream::app);

  out_ << "BX = " << (bitset<3>)bx_ << " Event : " << event_ << endl;

  for (unsigned int j = 0; j < tracklets_.size(); j++) {
    string vmproj = (layer_ > 0) ? tracklets_[j].first->vmstrlayer(layer_, tracklets_[j].second)
                                 : tracklets_[j].first->vmstrdisk(disk_, tracklets_[j].second);
    out_ << "0x";
    if (j < 16)
      out_ << "0";
    out_ << hex << j << dec;
    out_ << " " << vmproj << " " << Trklet::hexFormat(vmproj) << endl;
  }
  out_.close();

  bx_++;
  event_++;
  if (bx_ > 7)
    bx_ = 0;
}

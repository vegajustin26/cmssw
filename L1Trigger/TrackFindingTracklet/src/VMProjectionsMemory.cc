#include "L1Trigger/TrackFindingTracklet/interface/VMProjectionsMemory.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;
using namespace Trklet;

VMProjectionsMemory::VMProjectionsMemory(string name, const Settings* const settings, unsigned int iSector, double phimin, double phimax) :
  MemoryBase(name, settings, iSector) {
    phimin_ = phimin;
    phimax_ = phimax;
    string subname = name.substr(7, 2);
    layer_ = 0;
    disk_ = 0;
    if (subname == "L1")
      layer_ = 1;
    if (subname == "L2")
      layer_ = 2;
    if (subname == "L3")
      layer_ = 3;
    if (subname == "L4")
      layer_ = 4;
    if (subname == "L5")
      layer_ = 5;
    if (subname == "L6")
      layer_ = 6;
    if (subname == "D1")
      disk_ = 1;
    if (subname == "D2")
      disk_ = 2;
    if (subname == "D3")
      disk_ = 3;
    if (subname == "D4")
      disk_ = 4;
    if (subname == "D5")
      disk_ = 5;
    if (layer_ == 0 && disk_ == 0) {
      edm::LogPrint("Tracklet") << name << " subname = " << subname << " " << layer_ << " " << disk_;
    }
    assert((layer_ != 0) || (disk_ != 0));
}

void VMProjectionsMemory::writeVMPROJ(bool first) {
  
  std::string fname = "../data/MemPrints/VMProjections/VMProjections_";
  fname += getName();
  //get rid of duplicates
  int len = fname.size();
  if (fname[len - 2] == 'n' && fname[len - 1] > '1' && fname[len - 1] <= '9')
    return;
  //
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

#include "L1Trigger/TrackFindingTracklet/interface/VMStubsMEMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"

using namespace std;
using namespace trklet;

VMStubsMEMemory::VMStubsMEMemory(string name, const Settings* const settings, unsigned int iSector)
    : MemoryBase(name, settings, iSector) {

  unsigned int layerdisk=initLayerDisk(6);
  if (layerdisk<6) {
    binnedstubs_.resize(settings_->NLONGVMBINS());
  } else {
    //For disks we have NLONGVMBITS on each disk
    binnedstubs_.resize(2*settings_->NLONGVMBINS());
  }
  
}

void VMStubsMEMemory::writeStubs(bool first) {
  std::string fname = "../data/MemPrints/VMStubsME/VMStubs_";
  fname += getName();
  //get rid of duplicates
  int len = fname.size();
  if (fname[len - 2] == 'n' && fname[len - 1] > '1' && fname[len - 1] <= '9')
    return;

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

  for (unsigned int i = 0; i < binnedstubs_.size(); i++) {
    for (unsigned int j = 0; j < binnedstubs_[i].size(); j++) {
      string stub = binnedstubs_[i][j].stubindex().str();
      stub += "|" + binnedstubs_[i][j].bend().str();

      FPGAWord finepos = binnedstubs_[i][j].finerz();
      stub += "|" + finepos.str();
      out_ << hex << i << " " << j << dec << " " << stub << " " << trklet::hexFormat(stub) << endl;
    }
  }
  out_.close();

  bx_++;
  event_++;
  if (bx_ > 7)
    bx_ = 0;
}

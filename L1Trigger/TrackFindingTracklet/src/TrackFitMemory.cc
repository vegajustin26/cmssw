#include "L1Trigger/TrackFindingTracklet/interface/TrackFitMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/slhcevent.h"
#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"

using namespace std;
using namespace trklet;

TrackFitMemory::TrackFitMemory(
    string name, const Settings* const settings, unsigned int iSector, double phimin, double phimax)
    : MemoryBase(name, settings, iSector) {
  phimin_ = phimin;
  phimax_ = phimax;
}

void TrackFitMemory::writeTF(bool first) {
  std::string fname = "../data/MemPrints/FitTrack/TrackFit_";
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

  for (unsigned int j = 0; j < tracks_.size(); j++) {
    out_ << "0x";
    if (j < 16)
      out_ << "0";
    out_ << hex << j << dec << " ";
    out_ << tracks_[j]->trackfitstr() << " " << trklet::hexFormat(tracks_[j]->trackfitstr());
    out_ << "\n";
  }
  out_.close();

  bx_++;
  event_++;
  if (bx_ > 7)
    bx_ = 0;
}

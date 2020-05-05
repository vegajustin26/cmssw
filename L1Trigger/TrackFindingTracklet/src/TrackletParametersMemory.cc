#include "L1Trigger/TrackFindingTracklet/interface/TrackletParametersMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"
#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"


using namespace std;
using namespace Trklet;

TrackletParametersMemory::TrackletParametersMemory(string name, const Settings* const settings, unsigned int iSector)
    : MemoryBase(name, settings, iSector) {}

void TrackletParametersMemory::clean() {
  for (unsigned int i = 0; i < tracklets_.size(); i++) {
    delete tracklets_[i];
  }
  tracklets_.clear();
}



void TrackletParametersMemory::writeMatches(Globals* globals, int& matchesL1, int& matchesL3, int& matchesL5) {
  ofstream& out = globals->ofstream("nmatches.txt");
  for (unsigned int i = 0; i < tracklets_.size(); i++) {
    if ((tracklets_[i]->nMatches() + tracklets_[i]->nMatchesDisk()) > 0) {
      if (tracklets_[i]->layer() == 1)
        matchesL1++;
      if (tracklets_[i]->layer() == 3)
        matchesL3++;
      if (tracklets_[i]->layer() == 5)
        matchesL5++;
    }
    out << tracklets_[i]->layer() << " " << tracklets_[i]->disk() << " " << tracklets_[i]->nMatches() << " "
        << tracklets_[i]->nMatchesDisk() << endl;
  }
}

void TrackletParametersMemory::writeTPAR(bool first) {
  std::string fname = "../data/MemPrints/TrackletParameters/TrackletParameters_";
  fname += getName();
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
    string tpar = tracklets_[j]->trackletparstr();
    out_ << "0x";
    if (j < 16)
      out_ << "0";
    out_ << hex << j << dec;
    out_ << " " << tpar << " " << Trklet::hexFormat(tpar) << endl;
  }
  out_.close();

  bx_++;
  event_++;
  if (bx_ > 7)
    bx_ = 0;
}

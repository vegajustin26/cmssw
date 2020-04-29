// This class holds the tracklet parameters for the selected stub pairs
// This class owns the tracklets. Further modules only holds pointers
#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletParametersMemory_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletParametersMemory_h

#include "Tracklet.h"
#include "MemoryBase.h"
#include "GlobalHistTruth.h"

using namespace std;

class TrackletParametersMemory : public MemoryBase {
public:
  TrackletParametersMemory(string name, const Settings* const settings, unsigned int iSector, double phimin, double phimax)
    : MemoryBase(name, settings, iSector) {
    phimin_ = phimin;
    phimax_ = phimax;
  }

  void addTracklet(Tracklet *tracklet) {
    tracklets_.push_back(tracklet);
  }

  unsigned int nTracklets() const { return tracklets_.size(); }

  Tracklet *getFPGATracklet(unsigned int i) const { return tracklets_[i]; }

  void writeMatches(GlobalHistTruth* globals, int &matchesL1, int &matchesL3, int &matchesL5) {
    ofstream& out=globals->ofstream("nmatches.txt");
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

  void clean() {
    for (unsigned int i = 0; i < tracklets_.size(); i++) {
      delete tracklets_[i];
    }
    tracklets_.clear();
  }

  void writeTPAR(bool first) {
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

private:
  double phimin_;
  double phimax_;
  std::vector<Tracklet *> tracklets_;
};

#endif

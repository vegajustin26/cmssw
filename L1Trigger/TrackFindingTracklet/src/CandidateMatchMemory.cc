#include "L1Trigger/TrackFindingTracklet/interface/CandidateMatchMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;
using namespace Trklet;

CandidateMatchMemory::CandidateMatchMemory(string name, const Settings* const settings, unsigned int iSector)
    : MemoryBase(name, settings, iSector) {
  string subname = name.substr(3, 2);
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

void CandidateMatchMemory::writeCM(bool first) {
  std::string fname = "../data/MemPrints/Matches/CandidateMatches_";
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

  for (unsigned int j = 0; j < matches_.size(); j++) {
    string stubid = matches_[j].second.first->stubindex().str();                         // stub ID
    int projindex = (layer_ > 0) ? matches_[j].first.second : matches_[j].first.second;  // Allproj index
    FPGAWord tmp;
    if (projindex >= (1 << 7)) {
      projindex = (1 << 7) - 1;
    }
    tmp.set(projindex, 7, true, __LINE__, __FILE__);
    out_ << "0x";
    if (j < 16)
      out_ << "0";
    out_ << hex << j << dec;
    out_ << " " << tmp.str() << "|" << stubid << " " << Trklet::hexFormat(tmp.str() + stubid) << endl;
  }
  out_.close();

  bx_++;
  event_++;
  if (bx_ > 7)
    bx_ = 0;
}

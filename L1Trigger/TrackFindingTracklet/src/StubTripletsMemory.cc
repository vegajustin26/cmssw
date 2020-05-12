#include "L1Trigger/TrackFindingTracklet/interface/StubTripletsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/L1TStub.h"
#include "L1Trigger/TrackFindingTracklet/interface/Stub.h"

using namespace std;
using namespace Trklet;

StubTripletsMemory::StubTripletsMemory(string name, const Settings* const settings, unsigned int iSector)
    : MemoryBase(name, settings, iSector) {}

void StubTripletsMemory::writeST(bool first) {
  std::string fname = "../data/MemPrints/StubPairs/StubTriplets_";
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

  for (unsigned int j = 0; j < stubs1_.size(); j++) {
    string stub1index = stubs1_[j].first->stubindex().str();
    string stub2index = stubs2_[j].first->stubindex().str();
    string stub3index = stubs3_[j].first->stubindex().str();
    if (j < 16)
      out_ << "0";
    out_ << hex << j << dec;
    out_ << " " << stub1index << "|" << stub2index << "|" << stub3index << endl;
  }
  out_.close();

  bx_++;
  event_++;
  if (bx_ > 7)
    bx_ = 0;
}

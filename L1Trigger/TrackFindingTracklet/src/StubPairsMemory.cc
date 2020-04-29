#include "L1Trigger/TrackFindingTracklet/interface/StubPairsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/L1TStub.h"
#include "L1Trigger/TrackFindingTracklet/interface/Stub.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubTE.h"

using namespace std;
using namespace Trklet;


StubPairsMemory::StubPairsMemory(string name, const Settings* const settings, unsigned int iSector, double phimin, double phimax) :
  MemoryBase(name, settings, iSector) {
  phimin_ = phimin;
  phimax_ = phimax;
  }

void StubPairsMemory::writeSP(bool first) {
  std::string fname = "../data/MemPrints/StubPairs/StubPairs_";
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
  
  for (unsigned int j = 0; j < stubs1_.size(); j++) {
    string stub1index = stubs1_[j].stub().first->stubindex().str();
    string stub2index = stubs2_[j].stub().first->stubindex().str();
    out_ << "0x";
    if (j < 16)
      out_ << "0";
    out_ << hex << j << dec;
    out_ << " " << stub1index << "|" << stub2index << " " << Trklet::hexFormat(stub1index + stub2index) << endl;
  }
  out_.close();
  
  bx_++;
  event_++;
  if (bx_ > 7)
    bx_ = 0;
}

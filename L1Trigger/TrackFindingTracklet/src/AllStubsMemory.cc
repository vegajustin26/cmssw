#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"
//#include "L1Trigger/TrackFindingTracklet/interface/Stub.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"

using namespace std;
using namespace Trklet;


AllStubsMemory::AllStubsMemory(string name, const Settings* const settings, unsigned int iSector) :
  MemoryBase(name, settings, iSector) {
  
  //set the layer or disk that the memory is in
  initLayerDisk(3, layer_, disk_);
  
  assert(name.substr(5, 3) == "PHI");
}

void AllStubsMemory::writeStubs(bool first) {
  openFile(first, "../data/MemPrints/Stubs/AllStubs_");
  
  for (unsigned int j = 0; j < stubs_.size(); j++) {
    string stub = stubs_[j].first->str();
    out_ << "0x";
    if (j < 16)
      out_ << "0";
    out_ << hex << j << dec;
    out_ << " " << stub << " " << hexFormat(stub) << endl;
  }
  out_.close();
}


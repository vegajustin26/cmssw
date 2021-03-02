#include "L1Trigger/TrackFindingTracklet/interface/InputLinkMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/L1TStub.h"
#include "L1Trigger/TrackFindingTracklet/interface/Stub.h"

#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"

#include <iomanip>
#include <cmath>
#include <sstream>
#include <cctype>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace trklet;
using namespace std;

InputLinkMemory::InputLinkMemory(string name, Settings const& settings, unsigned int iSector, double, double)
    : MemoryBase(name, settings, iSector) {
  string subname = name.substr(5, 7);
  phiregion_ = subname[3] - 'A';
  assert(phiregion_ >= 0 && phiregion_ < 8);

  layerdisk_ = initLayerDisk(3);
}

void InputLinkMemory::addStub(const L1TStub& al1stub, const Stub& stub) {


  //Various consistency checks
  unsigned int stublayerdisk = stub.layerdisk();
  assert(stublayerdisk == layerdisk_);
  
  FPGAWord iphi = stub.phicorr();
  unsigned int nallbits = settings_.nbitsallstubs(layerdisk_);
  int phibin = iphi.bits(iphi.nbits() - nallbits, nallbits);
  int iphivmRaw = iphi.bits(iphi.nbits() - 5, 5);

  assert ( phibin==phiregion_);

  if (settings_.debugTracklet()) {
    edm::LogVerbatim("Tracklet") << "Will add stub in " << getName() << " "
                                 << "iphiwmRaw = " << iphivmRaw << " phi=" << al1stub.phi() << " z=" << al1stub.z()
                                 << " r=" << al1stub.r();
  }

  //Make new objects owned by the inputlink memory and save in list of stubs
  if (stubs_.size() < settings_.maxStep("Link")) {
    Stub* stubptr = new Stub(stub);
    stubptr->setl1tstub(new L1TStub(al1stub));

    stubs_.emplace_back(stubptr);
  }
  
}

void InputLinkMemory::writeStubs(bool first) {
  const string dirIS = settings_.memPath() + "InputStubs/";
  openFile(first, dirIS, "InputStubs_");

  for (unsigned int j = 0; j < stubs_.size(); j++) {
    string stub = stubs_[j]->str();
    out_ << std::setfill('0') << std::setw(2);
    out_ << hex << j << dec;
    out_ << " " << stub << " " << trklet::hexFormat(stub) << endl;
  }
  out_.close();
}

void InputLinkMemory::clean() {
  for (auto& stub : stubs_) {
    delete stub->l1tstub();
    delete stub;
  }
  stubs_.clear();
}

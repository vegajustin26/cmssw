#include "L1Trigger/TrackFindingTracklet/interface/DTCLinkMemory.h"
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

DTCLinkMemory::DTCLinkMemory(string name, Settings const& settings, double, double) : MemoryBase(name, settings) {}

void DTCLinkMemory::addStub(const L1TStub& al1stub, const Stub& stub) {
  //Make new objects owned by the dtclink memory and save in list of stubs
  if (stubs_.size() < settings_.maxStep("IR")) {
    Stub* stubptr = new Stub(stub);
    stubptr->setl1tstub(new L1TStub(al1stub));

    stubs_.emplace_back(stubptr);
  }
}

void DTCLinkMemory::writeStubs(bool first, unsigned int iSector) {
  iSector_ = iSector;

  static map<string, vector<int> > dtclayers{{"PS10G_1", {0, 6, 8, 10}},
                                             {"PS10G_2", {0, 7, 9}},
                                             {"PS10G_3", {1, 7}},
                                             {"PS10G_4", {6, 8, 10}},
                                             {"PS_1", {2, 7}},
                                             {"PS_2", {2, 9}},
                                             {"2S_1", {3, 4}},
                                             {"2S_2", {4}},
                                             {"2S_3", {5}},
                                             {"2S_4", {5, 8}},
                                             {"2S_5", {6, 9}},
                                             {"2S_6", {7, 10}}};

  const string dirIS = settings_.memPath() + "InputStubs/";
  openFile(first, dirIS, "Link_");

  for (unsigned int j = 0; j < stubs_.size(); j++) {
    string dtcname = stubs_[j]->l1tstub()->DTClink();
    int layerdisk = stubs_[j]->l1tstub()->layerdisk();

    //If the string starts with 'neg' skip the first three character
    int start = dtcname.substr(0, 3) == "neg" ? 3 : 0;

    //For the dtcbase name remove the leading 'neg' if in the name and the trailing '_A' or '_B'
    string dtcbase = dtcname.substr(start, dtcname.size() - 2 - start);

    vector<int> layers = dtclayers[dtcbase];

    int lcode = -1;
    for (unsigned int index = 0; index < layers.size(); index++) {
      if (layerdisk == layers[index]) {
        lcode = index;
      }
    }
    assert(lcode != -1);

    FPGAWord ldcode(lcode, 2, true);

    string stub = stubs_[j]->str() + "|" + ldcode.str() + "|1";
    out_ << std::setfill('0') << std::setw(2);
    out_ << hex << j << dec;
    out_ << " " << stub << " " << trklet::hexFormat(stub) << endl;
  }
  out_.close();
}

void DTCLinkMemory::clean() {
  for (auto& stub : stubs_) {
    delete stub->l1tstub();
    delete stub;
  }
  stubs_.clear();
}

#include "L1Trigger/TrackFindingTracklet/interface/TETableBase.h"

using namespace std;
using namespace trklet;

TETableBase::TETableBase(const Settings* settings) : settings_(settings) {}

int TETableBase::lookup(int, int) {
  assert(0);  //Should never get here - this should be a pure virtual fcn
};

void TETableBase::writeVMTable(std::string name, bool positive) {
  ofstream out;
  out.open(name.c_str());
  out << "{" << endl;
  for (unsigned int i = 0; i < table_.size(); i++) {
    if (i != 0) {
      out << "," << endl;
    }

    assert(nbits_ > 0);

    int itable = table_[i];
    if (positive) {
      if (table_[i] < 0) {
        itable = (1 << nbits_) - 1;
      }
    }

    out << itable;
  }
  out << endl << "};" << endl;
  out.close();
}

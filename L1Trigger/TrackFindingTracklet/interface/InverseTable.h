#ifndef L1Trigger_TrackFindingTracklet_interface_InverseTable_h
#define L1Trigger_TrackFindingTracklet_interface_InverseTable_h

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <vector>

class InverseTable {
public:
  InverseTable() {}

  ~InverseTable() {}

  void initR(int nbits, int offset, int invbits, bool pos) {
    nbits_ = nbits;
    entries_ = (1 << nbits);

    for (int i = 0; i < entries_; i++) {
      int idrrel = i;
      if (!pos) {
        if (i > ((1 << (nbits - 1)) - 1)) {
          idrrel = i - (1 << nbits);
        }
      }
      int idr = offset + idrrel;
      table_.push_back(round_int((1 << invbits) / (1.0 * idr)));
    }
  }
  void initT(int nbits, int offset, int invbits, bool pos) {
    nbits_ = nbits;
    entries_ = (1 << nbits);

    for (int i = 0; i < entries_; i++) {
      int itrel = i;
      if (!pos)
        itrel = i - entries_;
      int it = itrel << offset;
      int invt = round_int((1 << invbits) / (1.0 * it));
      table_.push_back(invt);
    }
  }

  void write(std::string fname) {
    ofstream out(fname.c_str());

    for (int i = 0; i < entries_; i++) {
      unsigned int tt = table_[i];
      out << std::hex << tt << endl;
    }
    out.close();
  }

  int lookup(int drrel) const {
    assert(drrel >= 0);
    assert(drrel < (1 << nbits_));
    return table_[drrel];
  }

private:
  int nbits_;
  int entries_;
  std::vector<int> table_;
};

#endif

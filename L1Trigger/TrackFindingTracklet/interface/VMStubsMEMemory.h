//This class holds the reduced VM stubs
#ifndef L1Trigger_TrackFindingTracklet_interface_VMStubsMEMemory_h
#define L1Trigger_TrackFindingTracklet_interface_VMStubsMEMemory_h

#include "L1TStub.h"
#include "Stub.h"
#include "VMStubME.h"
#include "MemoryBase.h"

using namespace std;

class VMStubsMEMemory : public MemoryBase {
public:
  VMStubsMEMemory(string name, const Settings* const settings, unsigned int iSector, double phimin, double phimax) :
  MemoryBase(name, settings, iSector) {
    phimin_ = phimin;
    phimax_ = phimax;
  }

  void addStub(VMStubME stub, unsigned int bin) {
    stubs_.push_back(stub);
    binnedstubs_[bin].push_back(stub);
  }

  unsigned int nStubs() const { return stubs_.size(); }

  VMStubME getVMStubME(unsigned int i) const { return stubs_[i]; }
  Stub* getFPGAStub(unsigned int i) const { return stubs_[i].stub().first; }
  L1TStub* getL1TStub(unsigned int i) const { return stubs_[i].stub().second; }
  std::pair<Stub*, L1TStub*> getStub(unsigned int i) const { return stubs_[i].stub(); }

  unsigned int nStubsBin(unsigned int bin) const {
    assert(bin < settings_->MEBinsDisks() * 2);
    return binnedstubs_[bin].size();
  }

  VMStubME getVMStubMEBin(unsigned int bin, unsigned int i) const {
    assert(bin < settings_->MEBinsDisks() * 2);
    assert(i < binnedstubs_[bin].size());
    return binnedstubs_[bin][i];
  }

  std::pair<Stub*, L1TStub*> getStubBin(unsigned int bin, unsigned int i) const {
    assert(bin < settings_->MEBinsDisks() * 2);
    assert(i < binnedstubs_[bin].size());
    return binnedstubs_[bin][i].stub();
  }

  void clean() {
    stubs_.clear();
    for (unsigned int i = 0; i < settings_->MEBinsDisks() * 2; i++) {
      binnedstubs_[i].clear();
    }
  }

  void writeStubs(bool first) {
    std::string fname = "../data/MemPrints/VMStubsME/VMStubs_";
    fname += getName();
    //get rid of duplicates
    int len = fname.size();
    if (fname[len - 2] == 'n' && fname[len - 1] > '1' && fname[len - 1] <= '9')
      return;
    //
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

    for (unsigned int i = 0; i < settings_->NLONGVMBINS(); i++) {
      for (unsigned int j = 0; j < binnedstubs_[i].size(); j++) {
        string stub = binnedstubs_[i][j].stubindex().str();
        stub += "|" + binnedstubs_[i][j].bend().str();

        FPGAWord finepos = binnedstubs_[i][j].finerz();
        stub += "|" + finepos.str();
        out_ << hex << i << " " << j << dec << " " << stub << " " << Trklet::hexFormat(stub) << endl;
      }
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
  std::vector<VMStubME> stubs_;
  //std::vector<VMStubME> binnedstubs_[settings_->MEBinsDisks() * 2];
  // LS HACK: the above doesn't work :(
  // settings_->MEBinsDisks() * 2 = 16
  std::vector<VMStubME> binnedstubs_[16];
};

#endif

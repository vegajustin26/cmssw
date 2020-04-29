#ifndef L1Trigger_TrackFindingTracklet_interface_VMStubsMEMemory_h
#define L1Trigger_TrackFindingTracklet_interface_VMStubsMEMemory_h

#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubME.h"

#include <vector>

class L1TStub;
class Stub;

namespace Trklet {

  class Settings;

  
  class VMStubsMEMemory : public MemoryBase {
  public:
    VMStubsMEMemory(string name, const Settings* const settings, unsigned int iSector, double phimin, double phimax);

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
    
    void writeStubs(bool first);

  private:
    double phimin_;
    double phimax_;
    std::vector<VMStubME> stubs_;
    //std::vector<VMStubME> binnedstubs_[settings_->MEBinsDisks() * 2];
    // LS HACK: the above doesn't work :(
    // settings_->MEBinsDisks() * 2 = 16
    std::vector<VMStubME> binnedstubs_[16];
  };

};
#endif

#ifndef L1Trigger_TrackFindingTracklet_interface_VMStubsMEMemory_h
#define L1Trigger_TrackFindingTracklet_interface_VMStubsMEMemory_h

#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubME.h"

#include <string>
#include <vector>

namespace trklet {

  class Settings;
  class Stub;
  class L1TStub;

  class VMStubsMEMemory : public MemoryBase {
  public:
    VMStubsMEMemory(std::string name, const Settings* const settings, unsigned int iSector);

    virtual ~VMStubsMEMemory() {}

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
      assert(bin < binnedstubs_.size());
      return binnedstubs_[bin].size();
    }

    VMStubME getVMStubMEBin(unsigned int bin, unsigned int i) const {
      assert(bin < binnedstubs_.size());
      assert(i < binnedstubs_[bin].size());
      return binnedstubs_[bin][i];
    }

    std::pair<Stub*, L1TStub*> getStubBin(unsigned int bin, unsigned int i) const {
      assert(bin < binnedstubs_.size());
      assert(i < binnedstubs_[bin].size());
      return binnedstubs_[bin][i].stub();
    }

    void clean() {
      stubs_.clear();
      for (unsigned int i = 0; i < binnedstubs_.size(); i++) {
        binnedstubs_[i].clear();
      }
    }

    void writeStubs(bool first);

  private:
    std::vector<VMStubME> stubs_; //TODO - not used and should be removed
    std::vector<std::vector<VMStubME> >binnedstubs_;
  };

};  // namespace trklet
#endif

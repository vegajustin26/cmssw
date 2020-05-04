#ifndef L1Trigger_TrackFindingTracklet_interface_StubPairsMemory_h
#define L1Trigger_TrackFindingTracklet_interface_StubPairsMemory_h

#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubTE.h"

#include <vector>

namespace Trklet {

  class Settings;
  class Stub;
  class L1TStub;

  class StubPairsMemory : public MemoryBase {
  public:
    StubPairsMemory(std::string name, const Settings* const settings, unsigned int iSector);

    void addStubPair(const VMStubTE& stub1,
                     const VMStubTE& stub2,
                     const unsigned index = 0,
                     const std::string& tedName = "") {
      stubs_.push_back(std::pair<VMStubTE, VMStubTE>(stub1, stub2));
      indices_.push_back(index);
      tedNames_.push_back(tedName);
    }

    unsigned int nStubPairs() const { return stubs_.size(); }

    VMStubTE getVMStub1(unsigned int i) const { return stubs_[i].first; }
    Stub* getFPGAStub1(unsigned int i) const { return stubs_[i].first.stub().first; }
    L1TStub* getL1TStub1(unsigned int i) const { return stubs_[i].first.stub().second; }
    std::pair<Stub*, L1TStub*> getStub1(unsigned int i) const { return stubs_[i].first.stub(); }

    VMStubTE getVMStub2(unsigned int i) const { return stubs_[i].second; }
    Stub* getFPGAStub2(unsigned int i) const { return stubs_[i].second.stub().first; }
    L1TStub* getL1TStub2(unsigned int i) const { return stubs_[i].second.stub().second; }
    std::pair<Stub*, L1TStub*> getStub2(unsigned int i) const { return stubs_[i].second.stub(); }

    unsigned getIndex(const unsigned i) const { return indices_.at(i); }
    const std::string& getTEDName(const unsigned i) const { return tedNames_.at(i); }

    void clean() {
      stubs_.clear();
      indices_.clear();
      tedNames_.clear();
    }

    void writeSP(bool first);

  private:
    std::vector<std::pair<VMStubTE, VMStubTE> > stubs_;

    std::vector<unsigned> indices_;
    std::vector<std::string> tedNames_;
  };

};  // namespace Trklet
#endif

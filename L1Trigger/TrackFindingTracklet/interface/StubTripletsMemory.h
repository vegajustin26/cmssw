#ifndef L1Trigger_TrackFindingTracklet_interface_StubTripletsMemory_h
#define L1Trigger_TrackFindingTracklet_interface_StubTripletsMemory_h

#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"

#include <vector>

namespace trklet {

  class Settings;
  class Stub;
  class L1TStub;

  class StubTripletsMemory : public MemoryBase {
  public:
    StubTripletsMemory(std::string name, const Settings* const settings, unsigned int iSector);

    ~StubTripletsMemory() = default;

    void addStubs(const Stub* stub1,
                  const Stub* stub2,
                  const Stub* stub3) {
      stubs1_.push_back(stub1);
      stubs2_.push_back(stub2);
      stubs3_.push_back(stub3);
    }

    unsigned int nStubTriplets() const { return stubs1_.size(); }

    const Stub* getFPGAStub1(unsigned int i) const { return stubs1_[i]; }
    //const L1TStub* getL1TStub1(unsigned int i) const { return stubs1_[i]->l1tstub(); }

    const Stub* getFPGAStub2(unsigned int i) const { return stubs2_[i]; }
    //const L1TStub* getL1TStub2(unsigned int i) const { return stubs2_[i]->l1tstub(); }

    const Stub* getFPGAStub3(unsigned int i) const { return stubs3_[i]; }
    //const L1TStub* getL1TStub3(unsigned int i) const { return stubs3_[i]->l1tstub(); }

    void clean() override {
      stubs1_.clear();
      stubs2_.clear();
      stubs3_.clear();
    }

    void writeST(bool first);

  private:
    std::vector<const Stub* > stubs1_;
    std::vector<const Stub* > stubs2_;
    std::vector<const Stub* > stubs3_;
  };

};  // namespace trklet
#endif

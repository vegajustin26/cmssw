#ifndef L1Trigger_TrackFindingTracklet_interface_AllStubsMemory_h
#define L1Trigger_TrackFindingTracklet_interface_AllStubsMemory_h

#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"

#include <utility>
#include <string>
#include <vector>

namespace trklet {

  class Settings;
  class Stub;
  class L1TStub;

  class AllStubsMemory : public MemoryBase {
  public:
    AllStubsMemory(std::string name, const Settings* const settings, unsigned int iSector);

    ~AllStubsMemory() = default;

    void addStub(std::pair<Stub*, L1TStub*> stub) { stubs_.push_back(stub); }

    unsigned int nStubs() const { return stubs_.size(); }

    std::pair<Stub*, L1TStub*> getStub(unsigned int i) const { return stubs_[i]; }

    void clean() { stubs_.clear(); }

    void writeStubs(bool first);

  private:
    std::vector<std::pair<Stub*, L1TStub*> > stubs_;
  };

};  // namespace trklet
#endif

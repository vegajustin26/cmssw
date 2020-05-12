// This class holds a list of stubs for an input link.
// This modules 'owns' the pointers to the stubs. All subsequent modules that handles stubs uses a pointer to the original stored here.
#ifndef L1Trigger_TrackFindingTracklet_interface_InputLinkMemory_h
#define L1Trigger_TrackFindingTracklet_interface_InputLinkMemory_h

#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"

#include <vector>

namespace trklet {

  class Settings;
  class Globals;
  class Stub;
  class L1TStub;
  class VMRouterPhiCorrTable;

  class InputLinkMemory : public MemoryBase {
  public:
    InputLinkMemory(std::string name, const Settings* const settings, unsigned int iSector, double, double);

    bool addStub(const Settings* settings, Globals* globals, L1TStub& al1stub, Stub& stub, std::string dtc);

    unsigned int nStubs() const { return stubs_.size(); }

    std::pair<Stub*, L1TStub*> getStub(unsigned int i) const { return stubs_[i]; }

    void writeStubs(bool first);

    void clean();

  private:
    std::vector<std::pair<Stub*, L1TStub*> > stubs_;
    int phiregion_;
    unsigned int layerdisk_;
  };

};  // namespace trklet
#endif

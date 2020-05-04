#ifndef L1Trigger_TrackFindingTracklet_interface_VMRouter_h
#define L1Trigger_TrackFindingTracklet_interface_VMRouter_h

#include "L1Trigger/TrackFindingTracklet/interface/ProcessBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/InputLinkMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubsMEMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubsTEMemory.h"

#include <vector>

namespace Trklet {

  class Settings;
  class Globals;
  class MemoryBase;

  class VMRouter : public ProcessBase {
  public:
    VMRouter(std::string name, const Settings* settings, Globals* global, unsigned int iSector);

    void addOutput(MemoryBase* memory, std::string output);

    void addInput(MemoryBase* memory, std::string input);

    void execute();

    void initFineBinTable();

    FPGAWord lookup(unsigned int iseed, unsigned int inner, FPGAWord z, FPGAWord r, bool negdisk, bool isPSmodule);

  private:
    //0-5 are the layers and 6-10 are the disks
    unsigned int layerdisk_;

    //overlapbits_ is the top bits of phicorr used to add or subtract one to see if stub should be added to
    //two VMs. nextrabits_ is the number of bits beyond the bits for the phivm that is used by overlapbits_
    unsigned int overlapbits_;
    unsigned int nextrabits_;

    int nbitsfinebintable_;
    std::vector<int> finebintable_;

    //The input stub memories
    std::vector<InputLinkMemory*> stubinputs_;

    //The all stub memories
    std::vector<AllStubsMemory*> allstubs_;

    //The VM stubs memories used by the MEs
    std::vector<VMStubsMEMemory*> vmstubsMEPHI_;

    //The VM stubs memories used by the TEs
    // vmstubsTEPHI_[i].first.first is the seed number [0,11] for
    // vmstubsTEPHI_[i].first.second is the stub position in the seed 0 is inner 1 is outer 2 is
    //                               third stub in extended tracking
    // vmstubsTEPHI_[i].second[iVM][n] is the VMStubsTEMemory for iVM and the nth copy
    std::vector<std::pair<std::pair<unsigned int, unsigned int>, std::vector<std::vector<VMStubsTEMemory*> > > >
        vmstubsTEPHI_;
  };
};  // namespace Trklet
#endif

#ifndef L1Trigger_TrackFindingTracklet_interface_StubPairsMemory_h
#define L1Trigger_TrackFindingTracklet_interface_StubPairsMemory_h

#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubTE.h"

#include <vector>

class L1TStub;

namespace Trklet {

  class Settings;
  class Stub;

  class StubPairsMemory : public MemoryBase {
  public:
    StubPairsMemory(string name, const Settings* const settings, unsigned int iSector, double phimin, double phimax);

    void addStubPair(const VMStubTE& stub1,
		     const VMStubTE& stub2,
		     const unsigned index = 0,
		     const std::string& tedName = "") {
      stubs1_.push_back(stub1);
      stubs2_.push_back(stub2);
      indices_.push_back(index);
      tedNames_.push_back(tedName);
    }

    unsigned int nStubPairs() const { return stubs1_.size(); }
    
    VMStubTE getVMStub1(unsigned int i) const { return stubs1_[i]; }
    Stub* getFPGAStub1(unsigned int i) const { return stubs1_[i].stub().first; }
    L1TStub* getL1TStub1(unsigned int i) const { return stubs1_[i].stub().second; }
    std::pair<Stub*, L1TStub*> getStub1(unsigned int i) const { return stubs1_[i].stub(); }
    
    VMStubTE getVMStub2(unsigned int i) const { return stubs2_[i]; }
    Stub* getFPGAStub2(unsigned int i) const { return stubs2_[i].stub().first; }
    L1TStub* getL1TStub2(unsigned int i) const { return stubs2_[i].stub().second; }
    std::pair<Stub*, L1TStub*> getStub2(unsigned int i) const { return stubs2_[i].stub(); }
    
    unsigned getIndex(const unsigned i) const { return indices_.at(i); }
    const std::string& getTEDName(const unsigned i) const { return tedNames_.at(i); }
    
    void clean() {
      stubs1_.clear();
      stubs2_.clear();
      indices_.clear();
      tedNames_.clear();
    }
    
    void writeSP(bool first);

    
  private:
    double phimin_;
    double phimax_;
    //FIXME should not be two vectors
    std::vector<VMStubTE> stubs1_;
    std::vector<VMStubTE> stubs2_;
    
    std::vector<unsigned> indices_;
    std::vector<std::string> tedNames_;
  };
  
};
#endif

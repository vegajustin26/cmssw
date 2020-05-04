#ifndef L1Trigger_TrackFindingTracklet_interface_CandidateMatchMemory_h
#define L1Trigger_TrackFindingTracklet_interface_CandidateMatchMemory_h

#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"

#include <vector>

namespace Trklet {

  class Settings;
  class Stub;
  class L1TStub;
  
  class CandidateMatchMemory : public MemoryBase {
  public:
    CandidateMatchMemory(std::string name, const Settings* const settings, unsigned int iSector);
    
    void addMatch(std::pair<Tracklet*, int> tracklet, std::pair<Stub*, L1TStub*> stub) {
      std::pair<std::pair<Tracklet*, int>, std::pair<Stub*, L1TStub*> > tmp(tracklet, stub);
      
      //Check for consistency
      for (unsigned int i = 0; i < matches_.size(); i++) {
	if (tracklet.first->TCID() < matches_[i].first.first->TCID()) {
	  edm::LogPrint("Tracklet") << "In " << getName() << " adding tracklet " << tracklet.first
				    << " with lower TCID : " << tracklet.first->TCID() << " than earlier TCID "
				    << matches_[i].first.first->TCID();
	  assert(0);
	}
      }
      matches_.push_back(tmp);
    }
    
    unsigned int nMatches() const { return matches_.size(); }
    
    Tracklet* getFPGATracklet(unsigned int i) const { return matches_[i].first.first; }
    std::pair<Stub*, L1TStub*> getStub(unsigned int i) const { return matches_[i].second; }
    std::pair<std::pair<Tracklet*, int>, std::pair<Stub*, L1TStub*> > getMatch(unsigned int i) const {
      return matches_[i];
    }
    
    void clean() { matches_.clear(); }
    
    void writeCM(bool first);
    
    int layer() const { return layer_; }
    int disk() const { return disk_; }
    
  private:

    std::vector<std::pair<std::pair<Tracklet*, int>, std::pair<Stub*, L1TStub*> > > matches_;
    
    int layer_;
    int disk_;
  };

};
#endif

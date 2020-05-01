#ifndef L1Trigger_TrackFindingTracklet_interface_MatchCalculator_h
#define L1Trigger_TrackFindingTracklet_interface_MatchCalculator_h

#include "L1Trigger/TrackFindingTracklet/interface/ProcessBase.h"

#include "L1Trigger/TrackFindingTracklet/interface/CandidateMatchMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/FullMatchMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllProjectionsMemory.h"


class L1TStub;

#include <vector>

namespace Trklet {

  class Settings; 
  class Globals;
  class Stub;
  class Tracklet;

  class MatchCalculator:public ProcessBase{
    
  public:
    
    MatchCalculator(string name, const Settings* settings, Globals* global, unsigned int iSector);
    
    void addOutput(MemoryBase* memory,string output);

    void addInput(MemoryBase* memory,string input);

    void execute();

    std::vector<std::pair<std::pair<Tracklet*,int>,std::pair<Stub*,L1TStub*> > > mergeMatches(std::vector<CandidateMatchMemory*>& candmatch);

  private:
    
    unsigned int layerdisk_;
    
    int fact_;
    int icorrshift_;
    int icorzshift_;
    int phi0shift_;
    double phimin_;
    double phimax_;
    double phioffset_;
    
    unsigned int phimatchcut_[12];
    unsigned int zmatchcut_[12];
    unsigned int rphicutPS_[12];
    unsigned int rphicut2S_[12];
    unsigned int rcutPS_[12];
    unsigned int rcut2S_[12];
    
    int ialphafactinner_[10];
    int ialphafactouter_[10];
    
    AllStubsMemory* allstubs_;
    AllProjectionsMemory* allprojs_;
    
    std::vector<CandidateMatchMemory*> matches_;
    
    std::vector<FullMatchMemory*> fullMatches_;
  };

};
#endif

//This class holds functional blocks of a sector
#ifndef L1Trigger_TrackFindingTracklet_interface_Sector_h
#define L1Trigger_TrackFindingTracklet_interface_Sector_h

#include "L1Trigger/TrackFindingTracklet/interface/InputLinkMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubsTEMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubsMEMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/StubPairsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/StubTripletsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletParametersMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletProjectionsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllProjectionsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMProjectionsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/CandidateMatchMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/FullMatchMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackFitMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/CleanTrackMemory.h"

#include "L1Trigger/TrackFindingTracklet/interface/VMRouter.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletEngine.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletEngineDisplaced.h"
#include "L1Trigger/TrackFindingTracklet/interface/TripletEngine.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletCalculator.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletProcessor.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletCalculatorDisplaced.h"
#include "L1Trigger/TrackFindingTracklet/interface/ProjectionRouter.h"
#include "L1Trigger/TrackFindingTracklet/interface/MatchEngine.h"
#include "L1Trigger/TrackFindingTracklet/interface/MatchCalculator.h"
#include "L1Trigger/TrackFindingTracklet/interface/MatchProcessor.h"
#include "L1Trigger/TrackFindingTracklet/interface/FitTrack.h"
#include "L1Trigger/TrackFindingTracklet/interface/PurgeDuplicate.h"
#include "L1Trigger/TrackFindingTracklet/interface/Util.h"
//#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"

namespace Trklet {

  class Settings;

  class Sector{
  public:
    Sector(unsigned int i, const Settings* settings, Globals* globals);
    ~Sector();
    
    bool addStub(L1TStub stub, string dtc);

    // Creates all required memory modules based on wiring map (args: module type, module instance)
    void addMem(string memType,string memName);

    // Creates all required processing modules based on wiring map (args: module type, module instance)
    void addProc(string procType,string procName);

    //--- Create all required proc -> mem module connections, based on wiring map
    //--- (args: memory instance & input/output proc modules it connects to in format procName.pinName)
    void addWire(string mem,string procinfull,string procoutfull);

    ProcessBase* getProc(string procName);
    MemoryBase* getMem(string memName);

    void writeInputStubs(bool first);
    void writeVMSTE(bool first);
    void writeVMSME(bool first);
    void writeAS(bool first);
    void writeSP(bool first); 
    void writeST(bool first);
    void writeTPAR(bool first);
    void writeTPROJ(bool first);
    void writeAP(bool first);
    void writeVMPROJ(bool first);
    void writeCM(bool first);
    void writeMC(bool first);
    void writeTF(bool first);
    void writeCT(bool first);
    
    void clean();

    // execute the different tracklet processing modules
    void executeVMR();
    void executeTE();
    void executeTED();
    void executeTRE();
    void executeTP();
    void executeTC();
    void executeTCD();
    void executePR();
    void executeME();
    void executeMC();
    void executeMP();
    void executeFT();
    void executePD(std::vector<Track*>& tracks);
    
    bool foundTrack(ofstream& outres, L1SimTrack simtrk);

    std::vector<Tracklet*> getAllTracklets();
    std::vector<std::pair<Stub*,L1TStub*> > getStubs() const;
    
    std::set<int> seedMatch(int itp);
    
    double phimin() const {return phimin_;}
    double phimax() const {return phimax_;}

  private:
    
    int isector_;
    const Settings* const settings_;
    Globals* globals_;
    double phimin_;
    double phimax_;
    
    std::map<string, MemoryBase*> Memories_;
    std::vector<MemoryBase*> MemoriesV_;
    std::vector<InputLinkMemory*> IL_;
    std::vector<AllStubsMemory*> AS_;
    std::vector<VMStubsTEMemory*> VMSTE_;
    std::vector<VMStubsMEMemory*> VMSME_;
    std::vector<StubPairsMemory*> SP_;
    std::vector<StubTripletsMemory*> ST_;
    std::vector<TrackletParametersMemory*> TPAR_;
    std::vector<TrackletProjectionsMemory*> TPROJ_;
    std::vector<AllProjectionsMemory*> AP_;
    std::vector<VMProjectionsMemory*> VMPROJ_;
    std::vector<CandidateMatchMemory*> CM_;
    std::vector<FullMatchMemory*> FM_;
    std::vector<TrackFitMemory*> TF_;
    std::vector<CleanTrackMemory*> CT_;
    
    std::map<string, ProcessBase*> Processes_;
    std::vector<VMRouter*> VMR_;
    std::vector<TrackletEngine*> TE_;
    std::vector<TrackletEngineDisplaced*> TED_;
    std::vector<TripletEngine*> TRE_;
    std::vector<TrackletProcessor*> TP_;
    std::vector<TrackletCalculator*> TC_;
    std::vector<TrackletCalculatorDisplaced*> TCD_;
    std::vector<ProjectionRouter*> PR_;
    std::vector<MatchEngine*> ME_;
    std::vector<MatchCalculator*> MC_;
    std::vector<MatchProcessor*> MP_;
    std::vector<FitTrack*> FT_;
    std::vector<PurgeDuplicate*> PD_;
    
  };
};
#endif

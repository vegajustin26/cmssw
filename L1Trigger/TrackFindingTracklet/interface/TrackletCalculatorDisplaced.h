#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletCalculatorDisplaced_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletCalculatorDisplaced_h

#include "L1Trigger/TrackFindingTracklet/interface/ProcessBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletProjectionsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/StubTripletsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletParametersMemory.h"

#include <vector>

class GlobalHistTruth;
class MemoryBase;
class Stub;
class L1TStub;

namespace Trklet {

  class Settings;
  
  class TrackletCalculatorDisplaced:public ProcessBase{
    
  public:
    
    TrackletCalculatorDisplaced(string name, const Settings* settings, GlobalHistTruth* global, unsigned int iSector);
    
    void addOutputProjection(TrackletProjectionsMemory* &outputProj, MemoryBase* memory);

    void addOutput(MemoryBase* memory,string output);
    
    void addInput(MemoryBase* memory,string input);

    void execute(); 

    void addDiskProj(Tracklet* tracklet, int disk);
    
    bool addLayerProj(Tracklet* tracklet, int layer); 

    void addProjection(int layer,int iphi,TrackletProjectionsMemory* trackletprojs, Tracklet* tracklet);

    void addProjectionDisk(int disk,int iphi,TrackletProjectionsMemory* trackletprojs, Tracklet* tracklet);

    bool LLLSeeding(Stub* innerFPGAStub, L1TStub* innerStub, Stub* middleFPGAStub, L1TStub* middleStub, Stub* outerFPGAStub, L1TStub* outerStub);

    bool DDLSeeding(Stub* innerFPGAStub, L1TStub* innerStub, Stub* middleFPGAStub, L1TStub* middleStub, Stub* outerFPGAStub, L1TStub* outerStub);

    bool LLDSeeding(Stub* innerFPGAStub, L1TStub* innerStub, Stub* middleFPGAStub, L1TStub* middleStub, Stub* outerFPGAStub, L1TStub* outerStub);

    int round_int( double r ) {
      return (r > 0.0) ? (r + 0.5) : (r - 0.5); 
    }
    
    void exactproj(double rproj,double rinv, double phi0, double d0,
		   double t, double z0, double r0,
		   double &phiproj, double &zproj,
		   double &phider, double &zder); 

    
    void exactprojdisk(double zproj, double rinv, double, double,  //phi0 and d0 are not used.
		       double t, double z0,
		       double x0, double y0,
		       double &phiproj, double &rproj,
		       double &phider, double &rder);
    
    void exacttracklet(double r1, double z1, double phi1,
		       double r2, double z2, double phi2,
		       double r3, double z3, double phi3,
		       int take3,
		       double& rinv, double& phi0, double &d0,
		       double& t, double& z0,
		       double phiproj[5], double zproj[5], 
		       double phiprojdisk[5], double rprojdisk[5],
		       double phider[5], double zder[5],
		       double phiderdisk[5], double rderdisk[5]);
    
  private:
    
    int TCIndex_;
    int layer_;
    int disk_;
    double phimin_;
    double phimax_;
    double rproj_[4];
    int lproj_[4];
    double zproj_[3];
    int dproj_[3];
    
    vector<double> toR_;
    vector<double> toZ_;
    
    unsigned int maxtracklet_; //maximum numbor of tracklets that be stored
    
    vector<AllStubsMemory*> innerallstubs_;
    vector<AllStubsMemory*> middleallstubs_;
    vector<AllStubsMemory*> outerallstubs_;
    vector<StubTripletsMemory*> stubtriplets_;
    
    TrackletParametersMemory* trackletpars_;
    
    TrackletProjectionsMemory* trackletproj_L1PHI1_;
    TrackletProjectionsMemory* trackletproj_L1PHI2_;
    TrackletProjectionsMemory* trackletproj_L1PHI3_;
    TrackletProjectionsMemory* trackletproj_L1PHI4_;
    TrackletProjectionsMemory* trackletproj_L1PHI5_;
    TrackletProjectionsMemory* trackletproj_L1PHI6_;
    TrackletProjectionsMemory* trackletproj_L1PHI7_;
    TrackletProjectionsMemory* trackletproj_L1PHI8_;
    
    TrackletProjectionsMemory* trackletproj_L2PHI1_;
    TrackletProjectionsMemory* trackletproj_L2PHI2_;
    TrackletProjectionsMemory* trackletproj_L2PHI3_;
    TrackletProjectionsMemory* trackletproj_L2PHI4_;
    
    TrackletProjectionsMemory* trackletproj_L3PHI1_;
    TrackletProjectionsMemory* trackletproj_L3PHI2_;
    TrackletProjectionsMemory* trackletproj_L3PHI3_;
    TrackletProjectionsMemory* trackletproj_L3PHI4_;
    
    TrackletProjectionsMemory* trackletproj_L4PHI1_;
    TrackletProjectionsMemory* trackletproj_L4PHI2_;
    TrackletProjectionsMemory* trackletproj_L4PHI3_;
    TrackletProjectionsMemory* trackletproj_L4PHI4_;
    
    TrackletProjectionsMemory* trackletproj_L5PHI1_;
    TrackletProjectionsMemory* trackletproj_L5PHI2_;
    TrackletProjectionsMemory* trackletproj_L5PHI3_;
    TrackletProjectionsMemory* trackletproj_L5PHI4_;
    
    TrackletProjectionsMemory* trackletproj_L6PHI1_;
    TrackletProjectionsMemory* trackletproj_L6PHI2_;
    TrackletProjectionsMemory* trackletproj_L6PHI3_;
    TrackletProjectionsMemory* trackletproj_L6PHI4_;
    
    TrackletProjectionsMemory* trackletproj_D1PHI1_;
    TrackletProjectionsMemory* trackletproj_D1PHI2_;
    TrackletProjectionsMemory* trackletproj_D1PHI3_;
    TrackletProjectionsMemory* trackletproj_D1PHI4_;
    
    TrackletProjectionsMemory* trackletproj_D2PHI1_;
    TrackletProjectionsMemory* trackletproj_D2PHI2_;
    TrackletProjectionsMemory* trackletproj_D2PHI3_;
    TrackletProjectionsMemory* trackletproj_D2PHI4_;
    
    TrackletProjectionsMemory* trackletproj_D3PHI1_;
    TrackletProjectionsMemory* trackletproj_D3PHI2_;
    TrackletProjectionsMemory* trackletproj_D3PHI3_;
    TrackletProjectionsMemory* trackletproj_D3PHI4_;
    
    TrackletProjectionsMemory* trackletproj_D4PHI1_;
    TrackletProjectionsMemory* trackletproj_D4PHI2_;
    TrackletProjectionsMemory* trackletproj_D4PHI3_;
    TrackletProjectionsMemory* trackletproj_D4PHI4_;
    
    TrackletProjectionsMemory* trackletproj_D5PHI1_;
    TrackletProjectionsMemory* trackletproj_D5PHI2_;
    TrackletProjectionsMemory* trackletproj_D5PHI3_;
    TrackletProjectionsMemory* trackletproj_D5PHI4_;
    
    TrackletProjectionsMemory* trackletproj_L1Plus_; 
    TrackletProjectionsMemory* trackletproj_L1Minus_;
    
    TrackletProjectionsMemory* trackletproj_L2Plus_; 
    TrackletProjectionsMemory* trackletproj_L2Minus_;
    
    TrackletProjectionsMemory* trackletproj_L3Plus_; 
    TrackletProjectionsMemory* trackletproj_L3Minus_;
    
    TrackletProjectionsMemory* trackletproj_L4Plus_; 
    TrackletProjectionsMemory* trackletproj_L4Minus_;
    
    TrackletProjectionsMemory* trackletproj_L5Plus_; 
    TrackletProjectionsMemory* trackletproj_L5Minus_;
    
    TrackletProjectionsMemory* trackletproj_L6Plus_; 
    TrackletProjectionsMemory* trackletproj_L6Minus_;
    
    
    TrackletProjectionsMemory* trackletproj_D1Plus_; 
    TrackletProjectionsMemory* trackletproj_D1Minus_;
    
    TrackletProjectionsMemory* trackletproj_D2Plus_; 
    TrackletProjectionsMemory* trackletproj_D2Minus_;
    
    TrackletProjectionsMemory* trackletproj_D3Plus_; 
    TrackletProjectionsMemory* trackletproj_D3Minus_;
    
    TrackletProjectionsMemory* trackletproj_D4Plus_; 
    TrackletProjectionsMemory* trackletproj_D4Minus_;
    
    TrackletProjectionsMemory* trackletproj_D5Plus_; 
    TrackletProjectionsMemory* trackletproj_D5Minus_;
  };

};
#endif

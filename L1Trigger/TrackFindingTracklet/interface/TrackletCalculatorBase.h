#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletCalculatorBase_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletCalculatorBase_h

#include "L1Trigger/TrackFindingTracklet/interface/ProcessBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletParametersMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletProjectionsMemory.h"

#include <vector>

class GlobalHistTruth;
class Tracklet;
class L1TStub;

namespace Trklet {

  class Settings;
  class Stub;
  
  class TrackletCalculatorBase : public ProcessBase{

  public:
    
    TrackletCalculatorBase(string name, const Settings* const settings, GlobalHistTruth* global, unsigned int iSector);
  
    void exacttracklet(double r1, double z1, double phi1,
		       double r2, double z2, double phi2, double,
		       double& rinv, double& phi0,
		       double& t, double& z0,
		       double phiproj[4], double zproj[4], 
		       double phider[4], double zder[4],
		       double phiprojdisk[5], double rprojdisk[5], 
		       double phiderdisk[5], double rderdisk[5]);
    
    void exacttrackletdisk(double r1, double z1, double phi1,
			   double r2, double z2, double phi2, double,
			   double& rinv, double& phi0,
			   double& t, double& z0,
			   double phiprojLayer[3], double zprojLayer[3], 
			   double phiderLayer[3], double zderLayer[3],
			   double phiproj[3], double rproj[3], 
			   double phider[3], double rder[3]);
    
    void exacttrackletOverlap(double r1, double z1, double phi1,
			      double r2, double z2, double phi2, double,
			      double& rinv, double& phi0,
			      double& t, double& z0,
			      double phiprojLayer[3], double zprojLayer[3], 
			      double phiderLayer[3], double zderLayer[3],
			      double phiproj[3], double rproj[3], 
			      double phider[3], double rder[3]);
    
    
    void exactproj(double rproj,double rinv,double phi0,double t, double z0,
		   double &phiproj, double &zproj, double &phider, double &zder);
    
    void exactprojdisk(double zproj,double rinv,double phi0, double t, double z0,
		       double &phiproj, double &rproj, double &phider, double &rder);
    
    void addDiskProj(Tracklet* tracklet, int disk);
    bool addLayerProj(Tracklet* tracklet, int layer);
    
    void addProjection(int layer,int iphi,TrackletProjectionsMemory* trackletprojs, Tracklet* tracklet);
    void addProjectionDisk(int disk,int iphi,TrackletProjectionsMemory* trackletprojs, Tracklet* tracklet);
    
    bool goodTrackPars(bool goodrinv, bool goodz0); 
    
    bool inSector(int iphi0, int irinv, double phi0approx, double rinvapprox);
    
    bool barrelSeeding(Stub* innerFPGAStub, L1TStub* innerStub, Stub* outerFPGAStub, L1TStub* outerStub);
    bool diskSeeding(Stub* innerFPGAStub,L1TStub* innerStub,Stub* outerFPGAStub,L1TStub* outerStub); 
    bool overlapSeeding(Stub* innerFPGAStub, L1TStub* innerStub, Stub* outerFPGAStub, L1TStub* outerStub);
    
  protected:
    
    int iSeed_;
    int TCIndex_;
    
    double phimin_;
    double phimax_;
    double phioffset_;
    
    int layer_;
    int disk_;
    
    int lproj_[4];
    int dproj_[3];
    
    double rproj_[4];
    double zproj_[3];
    double zprojoverlap_[4];
    
    TrackletParametersMemory* trackletpars_;
    
    //First index is layer/disk second is phi region
    std::vector<std::vector<TrackletProjectionsMemory*> > trackletprojlayers_;
    std::vector<std::vector<TrackletProjectionsMemory*> > trackletprojdisks_;  
  };
  
};
#endif

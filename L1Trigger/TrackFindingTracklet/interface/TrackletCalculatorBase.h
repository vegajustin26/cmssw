#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletCalculatorBase_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletCalculatorBase_h

#include "L1Trigger/TrackFindingTracklet/interface/ProcessBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletParametersMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletProjectionsMemory.h"

#include <vector>

namespace trklet {

  class Settings;
  class Globals;
  class Stub;
  class L1TStub;
  class Tracklet;

  class TrackletCalculatorBase : public ProcessBase {
  public:
    TrackletCalculatorBase(std::string name, Settings const& settings, Globals* global, unsigned int iSector);

    ~TrackletCalculatorBase() override = default;

    void exacttracklet(double r1,
                       double z1,
                       double phi1,
                       double r2,
                       double z2,
                       double phi2,
                       double,
                       double& rinv,
                       double& phi0,
                       double& t,
                       double& z0,
                       double phiproj[N_PROJLAYER],
                       double zproj[N_PROJLAYER],
                       double phider[N_PROJLAYER],
                       double zder[N_PROJLAYER],
                       double phiprojdisk[N_PROJDISK],
                       double rprojdisk[N_PROJDISK],
                       double phiderdisk[N_PROJDISK],
                       double rderdisk[N_PROJDISK]);

    void exacttrackletdisk(double r1,
                           double z1,
                           double phi1,
                           double r2,
                           double z2,
                           double phi2,
                           double,
                           double& rinv,
                           double& phi0,
                           double& t,
                           double& z0,
                           double phiprojLayer[N_DPROJMAX],
                           double zprojLayer[N_DPROJMAX],
                           double phiderLayer[N_DPROJMAX],
                           double zderLayer[N_DPROJMAX],
                           double phiproj[N_DPROJMAX],
                           double rproj[N_DPROJMAX],
                           double phider[N_DPROJMAX],
                           double rder[N_DPROJMAX]);

    void exacttrackletOverlap(double r1,
                              double z1,
                              double phi1,
                              double r2,
                              double z2,
                              double phi2,
                              double,
                              double& rinv,
                              double& phi0,
                              double& t,
                              double& z0,
                              double phiprojLayer[N_DPROJMAX],
                              double zprojLayer[N_DPROJMAX],
                              double phiderLayer[N_DPROJMAX],
                              double zderLayer[N_DPROJMAX],
                              double phiproj[N_DPROJMAX],
                              double rproj[N_DPROJMAX],
                              double phider[N_DPROJMAX],
                              double rder[N_DPROJMAX]);

    void exactproj(double rproj,
                   double rinv,
                   double phi0,
                   double t,
                   double z0,
                   double& phiproj,
                   double& zproj,
                   double& phider,
                   double& zder);

    void exactprojdisk(double zproj,
                       double rinv,
                       double phi0,
                       double t,
                       double z0,
                       double& phiproj,
                       double& rproj,
                       double& phider,
                       double& rder);

    void addDiskProj(Tracklet* tracklet, int disk);
    bool addLayerProj(Tracklet* tracklet, int layer);

    void addProjection(int layer, int iphi, TrackletProjectionsMemory* trackletprojs, Tracklet* tracklet);
    void addProjectionDisk(int disk, int iphi, TrackletProjectionsMemory* trackletprojs, Tracklet* tracklet);

    bool goodTrackPars(bool goodrinv, bool goodz0);

    bool inSector(int iphi0, int irinv, double phi0approx, double rinvapprox);

    bool barrelSeeding(const Stub* innerFPGAStub,
                       const L1TStub* innerStub,
                       const Stub* outerFPGAStub,
                       const L1TStub* outerStub);
    bool diskSeeding(const Stub* innerFPGAStub,
                     const L1TStub* innerStub,
                     const Stub* outerFPGAStub,
                     const L1TStub* outerStub);
    bool overlapSeeding(const Stub* innerFPGAStub,
                        const L1TStub* innerStub,
                        const Stub* outerFPGAStub,
                        const L1TStub* outerStub);

  protected:
    unsigned int iSeed_;
    unsigned int layerdisk1_;
    unsigned int layerdisk2_;

    int TCIndex_;

    double phioffset_;

    int layer_;
    int disk_;

    TrackletParametersMemory* trackletpars_;

    //First index is layer/disk second is phi region
    std::vector<std::vector<TrackletProjectionsMemory*> > trackletprojlayers_;
    std::vector<std::vector<TrackletProjectionsMemory*> > trackletprojdisks_;
  };

};  // namespace trklet
#endif

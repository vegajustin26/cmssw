// This class holds a list of stubs that are in a given layer and DCT region
#ifndef L1Trigger_TrackFindingTracklet_interface_Cabling_h
#define L1Trigger_TrackFindingTracklet_interface_Cabling_h

#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/DTC.h"
#include "L1Trigger/TrackFindingTracklet/interface/DTCLink.h"

#include <vector>

namespace trklet {

  class Settings;

  class Cabling {
  public:
    Cabling();

    virtual ~Cabling() {}
    
    void init(std::string dtcconfig, std::string moduleconfig);

    std::string dtc(int layer, int ladder, int module);

    void addphi(std::string dtc, double phi, int layer, int module);

    void writephirange();

    std::vector<std::string> DTCs() const;

  private:
    std::vector<DTCLink> links_;
    std::map<std::string, DTC> dtcranges;
    std::map<std::string, DTC> dtcs;
    std::map<int, std::map<int, std::map<int, std::string> > > modules;
  };

};  // namespace trklet
#endif

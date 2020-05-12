#ifndef L1Trigger_TrackFindingTracklet_interface_Stub_h
#define L1Trigger_TrackFindingTracklet_interface_Stub_h

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>

#include "L1Trigger/TrackFindingTracklet/interface/FPGAWord.h"
#include "L1Trigger/TrackFindingTracklet/interface/L1TStub.h"
#include "L1Trigger/TrackFindingTracklet/interface/Util.h"

namespace Trklet {

  class Settings;

  class Stub {
  public:
    Stub(const Settings* const settings);

    Stub(const L1TStub& stub, const Settings* const settings, double phiminsec, double phimaxsec);

    ~Stub() {}

    FPGAWord iphivmFineBins(int VMbits, int finebits) const;

    std::string str() const {
      if (layer_.value() != -1) {
        return r_.str() + "|" + z_.str() + "|" + phi_.str() + "|" + bend_.str();
      } else {
        if (isPSmodule()) {
          return r_.str() + "|" + z_.str() + "|" + phi_.str() + "|" + bend_.str();
        } else {
          return "000" + r_.str() + "|" + z_.str() + "|" + phi_.str() + "|" + alphanew_.str() + "|" + bend_.str();
        }
      }
    }

    std::string strbare() const { return bend_.str() + r_.str() + z_.str() + phi_.str(); }

    std::string phiregionaddressstr() const;  //TODO - should migrate away from suing this method

    FPGAWord phiregion() const;  //TODO - should migrate away from using this method

    void setAllStubIndex(int nstub);  //TODO - should migrate away from using this method

    void setPhiCorr(int phiCorr);

    FPGAWord bend() const { return bend_; }

    FPGAWord r() const { return r_; }
    FPGAWord z() const { return z_; }
    FPGAWord phi() const { return phi_; }
    FPGAWord phicorr() const { return phicorr_; }
    FPGAWord alphanew() const { return alphanew_; }

    FPGAWord stubindex() const { return stubindex_; }
    FPGAWord layer() const { return layer_; }
    FPGAWord disk() const { return disk_; }
    unsigned int layerdisk() const;

    bool isBarrel() const { return layer_.value() != -1; }
    bool isDisk() const { return disk_.value() != 0; }

    bool isPSmodule() const { return isPSmodule_; }

    double rapprox() const;
    double zapprox() const;
    double phiapprox(double phimin, double) const;

    void setfiner(int finer) { finer_.set(finer, 4, true, __LINE__, __FILE__); }

    FPGAWord finer() const { return finer_; }

    void setfinez(int finez) { finez_.set(finez, 4, true, __LINE__, __FILE__); }

    FPGAWord finez() const { return finez_; }

  private:
    bool isPSmodule_;  //TODO should not be used can be removed
    FPGAWord layer_;
    FPGAWord disk_;
    FPGAWord r_;
    FPGAWord z_;
    FPGAWord phi_;
    FPGAWord alphanew_;

    FPGAWord bend_;

    FPGAWord phicorr_;  //Corrected for bend to nominal radius

    FPGAWord stubindex_;

    FPGAWord finer_;  //TODO should not be member data
    FPGAWord finez_;  //TODO should not be member data

    const Settings* const settings_;
  };

};  // namespace Trklet
#endif

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

  class Stub{
    
  public:
    
    Stub(const Settings* const settings);

    Stub(const L1TStub& stub,const Settings* const settings, double phiminsec, double phimaxsec);

    ~Stub() {  }

    //Returns a number from 0 to 31 //FIXME should not be used
    unsigned int iphivmRaw() const;

    FPGAWord iphivmFineBins(int VMbits, int finebits) const;

    std::string str() const {
      if (layer_.value()!=-1) {
	return r_.str()+"|"+z_.str()+"|"+phi_.str()+"|"+bend_.str();
      }
      else {
	if (isPSmodule()) {
	  return r_.str()+"|"+z_.str()+"|"+phi_.str()+"|"+bend_.str();
	}
	else {
	  return "000"+r_.str()+"|"+z_.str()+"|"+phi_.str()+"|"+alphanew_.str()+"|"+bend_.str();
	}
      }
    }
    
    std::string strbare() const {
      return bend_.str()+r_.str()+z_.str()+phi_.str();
    }
    
    std::string phiregionaddressstr();
    
    FPGAWord phiregion() const;

    void setAllStubIndex(int nstub);
    
    void setPhiCorr(int phiCorr);


    FPGAWord bend() const {return bend_; }
    
    FPGAWord r() const { return r_; }
    FPGAWord z() const { return z_; }
    FPGAWord phi() const { return phi_; }
    FPGAWord phicorr() const { return phicorr_; }
    FPGAWord alphanew() const { return alphanew_; }
    
    //FIXME should remove these...
    int ir() const { return r_.value(); }
    int iz() const { return z_.value(); }
    int iphi() const { return phi_.value(); }
    
    FPGAWord stubindex() const {return stubindex_;}
    FPGAWord layer() const {return layer_;}
    FPGAWord disk() const {return disk_;}
    
    bool isBarrel() const {return layer_.value()!=-1;}
    bool isDisk() const {return disk_.value()!=0;}
    
    bool isPSmodule() const {return isPSmodule_;}
    
    double rapprox() const;
    double zapprox() const;
    double phiapprox(double phimin, double) const;
    
    void setfiner(int finer) {
      finer_.set(finer,4,true,__LINE__,__FILE__);
    }
    
    FPGAWord finer() const {
      return finer_;
    }
    
    void setfinez(int finez) {
      finez_.set(finez,4,true,__LINE__,__FILE__);
    }
    
    FPGAWord finez() const {
      return finez_;
    }
    
    
    
  private:
    
    bool isPSmodule_;  //FIXME can be removed
    FPGAWord layer_;  
    FPGAWord disk_;  
    FPGAWord r_;
    FPGAWord z_;
    FPGAWord phi_;
    FPGAWord alphanew_;
    
    FPGAWord bend_;
    
    FPGAWord phicorr_;  //Corrected for bend to nominal radius
    
    FPGAWord stubindex_;
    
    FPGAWord finer_;   //FIXME should not be member data
    FPGAWord finez_;   //FIXME should not be member data
    
    const Settings* const settings_;
    
  };
  
};
#endif




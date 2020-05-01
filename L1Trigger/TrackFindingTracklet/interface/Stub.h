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
      std::ostringstream oss;
      if (layer_.value()!=-1) {
	oss << r_.str()<<"|"<< z_.str()<<"|"<< phi_.str()<<"|"<<bend_.str();
      }
      else {
	if (isPSmodule()) {
	  oss <<r_.str()<<"|"<< z_.str()<<"|"<< phi_.str()<<"|"<<bend_.str();
	}
	else {
	  oss << "000"<<r_.str()<<"|"<< z_.str()<<"|"<< phi_.str()<<"|"<<alphanew_.str()<<"|"<<bend_.str();
	}
      }
      return oss.str(); 
    }
    
    std::string strbare() const {
      std::ostringstream oss;
      oss << bend_.str()<<r_.str()<< z_.str()<< phi_.str();
      return oss.str();
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
    
    int round_int( double r ) {
      return (r > 0.0) ? (r + 0.5) : (r - 0.5); 
    }
    
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
    
    
    //Should be optimized by layer - now first implementation to make sure it works OK
static int bendencode(double bend, bool isPS) {
  
  int ibend=2.0*bend;
  
  assert(std::abs(ibend-2.0*bend)<0.1);
  
  if (isPS) {
    
    if (ibend==0||ibend==1) return 0;
    if (ibend==2||ibend==3) return 1;
    if (ibend==4||ibend==5) return 2;
    if (ibend>=6) return 3;
    if (ibend==-1||ibend==-2) return 4;
    if (ibend==-3||ibend==-4) return 5;
    if (ibend==-5||ibend==-6) return 6;
    if (ibend<=-7) return 7;
    
    assert(0);
    
  }
  
  
  if (ibend==0||ibend==1) return 0;
  if (ibend==2||ibend==3) return 1;
  if (ibend==4||ibend==5) return 2;
  if (ibend==6||ibend==7) return 3;
  if (ibend==8||ibend==9) return 4;
  if (ibend==10||ibend==11) return 5;
  if (ibend==12||ibend==13) return 6;
  if (ibend>=14) return 7;
  if (ibend==-1||ibend==-2) return 8;
  if (ibend==-3||ibend==-4) return 9;
  if (ibend==-5||ibend==-6) return 10;
  if (ibend==-7||ibend==-8) return 11;
  if (ibend==-9||ibend==-10) return 12;
  if (ibend==-11||ibend==-12) return 13;
  if (ibend==-13||ibend==-14) return 14;
  if (ibend<=-15) return 15;
  
  edm::LogVerbatim("Tracklet") << "bend ibend : "<<bend<<" "<<ibend;
  
  assert(0);
  
  
}
      
    //Should be optimized by layer - now first implementation to make sure it works OK
static double benddecode(int ibend, bool isPS) {
  
  if (isPS) {
    
    if (ibend==0) return 0.25;
    if (ibend==1) return 1.25;
    if (ibend==2) return 2.25;
    if (ibend==3) return 3.25;
    if (ibend==4) return -0.75;
    if (ibend==5) return -1.75;
    if (ibend==6) return -2.75;
    if (ibend==7) return -3.75;
    
    assert(0);
  }
  
  if (ibend==0) return 0.25;
  if (ibend==1) return 1.25;
  if (ibend==2) return 2.25;
  if (ibend==3) return 3.25;
  if (ibend==4) return 4.25;
  if (ibend==5) return 5.25;
  if (ibend==6) return 6.25;
  if (ibend==7) return 7.25;
  if (ibend==8) return -0.75;
  if (ibend==9) return -1.75;
  if (ibend==10) return -2.75;
  if (ibend==11) return -3.75;
  if (ibend==12) return -4.75;
  if (ibend==13) return -5.75;
  if (ibend==14) return -6.75;
  if (ibend==15) return -7.75;
  
  assert(0);
  
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




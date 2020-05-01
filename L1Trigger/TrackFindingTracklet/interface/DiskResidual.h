#ifndef L1Trigger_TrackFindingTracklet_interface_DiskResidual_h
#define L1Trigger_TrackFindingTracklet_interface_DiskResidual_h

#include "L1Trigger/TrackFindingTracklet/interface/FPGAWord.h"

class L1TStub;
class Stub;

namespace Trklet {

  class Settings;

  class DiskResidual {
  public:

    DiskResidual() { valid_ = false; }

    void init(const Settings* settings,
	      int disk,
	      int iphiresid,
	      int irresid,
	      int istubid,
	      double phiresid,
	      double rresid,
	      double phiresidapprox,
	      double rresidapprox,
	      double zstub,
	      double alpha,
	      FPGAWord ialpha,
	      std::pair<Stub*, L1TStub*> stubptrs);

    bool valid() const { return valid_; }
    
    FPGAWord fpgaphiresid() const {
      assert(valid_);
      return fpgaphiresid_;
    };
    
    FPGAWord fpgarresid() const {
      assert(valid_);
      return fpgarresid_;
    };
    
    FPGAWord fpgastubid() const {
      assert(valid_);
      return fpgastubid_;
    };
    
    double phiresid() const {
      assert(valid_);
      return phiresid_;
    };
    
    double rresid() const {
      assert(valid_);
      return rresid_;
    };
    
    double phiresidapprox() const {
      assert(valid_);
      return phiresidapprox_;
    };
    
    double rresidapprox() const {
      assert(valid_);
      return rresidapprox_;
    };
    
    double zstub() const {
      assert(valid_);
      return zstub_;
    };
    
    double alpha() const {
      assert(valid_);
      return alpha_;
    };
    
    FPGAWord ialpha() const {
      assert(valid_);
      return ialpha_;
    };
    
    std::pair<Stub*, L1TStub*> stubptrs() const {
      assert(valid_);
      return stubptrs_;
    };
    
  protected:
    bool valid_;
    
    int disk_;
    
    FPGAWord fpgaphiresid_;
    FPGAWord fpgarresid_;
    FPGAWord fpgastubid_;
    
    double phiresid_;
    double rresid_;
    
    double phiresidapprox_;
    double rresidapprox_;
    
    double zstub_;
    double alpha_;
    FPGAWord ialpha_;
    std::pair<Stub*, L1TStub*> stubptrs_;
  };

};
#endif

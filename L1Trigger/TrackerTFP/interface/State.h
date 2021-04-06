#ifndef L1Trigger_TrackerTFP_State_h
#define L1Trigger_TrackerTFP_State_h

#include "L1Trigger/TrackerDTC/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/DataFormats.h"

#include <vector>
#include <numeric>


namespace trackerTFP {

  // 
  class State {
  public:
    //
    State(State* state);
    // proto state constructor
    State(const DataFormats* dataFormats, TrackKFin* track, int trackId);
    // combinatoric state constructor
    State(State* state, StubKFin* stub);
    // updated state constructor
    State(State* state, const std::vector<double>& doubles);
    ~State(){}

    //
    void finish();
    //
    int numSkippedLayers() const { return numSkippedLayers_; }
    //
    int numConsistentLayers() const { return numConsistentLayers_; }
    //
    TrackKFin* track() const { return track_; }
    //
    State* parent() const { return parent_; }
    //
    StubKFin*  stub() const { return stub_; }
    //
    double r() const { return stub_->r(); }
    //
    double phi() const { return stub_->phi(); }
    //
    double z() const { return stub_->z(); }
    //
    int sectorPhi() const { return track_->sectorPhi(); }
    //
    int sectorEta() const { return track_->sectorEta(); }
    //
    const TTBV& hitPattern() const { return hitPattern_; }
    //
    int trackId() const { return trackId_; }
    //
    TTBV maybePattern() const { return track_->maybePattern(); }
    //
    int sector() const { return sectorPhi() * setup_->numSectorsEta() + sectorEta(); }
    //
    const std::vector<int>& layerMap() const { return layerMap_; }
    //
    bool barrel() const { return setup_->barrel(stub_->ttStubRef()); }
    //
    bool psModule() const { return setup_->psModule(stub_->ttStubRef()); }
    //
    int layer() const { return stub_->layer(); }
    //
    void x0(double d) { x0_ = d; }
    //
    void x1(double d) { x1_ = d; }
    //
    void x2(double d) { x2_ = d; }
    //
    void x3(double d) { x3_ = d; }
    //
    double x0() const { return x0_; }
    //
    double x1() const { return x1_; }
    //
    double x2() const { return x2_; }
    //
    double x3() const { return x3_; }
    //
    double C00() const { return C00_; }
    //
    double C01() const { return C01_; }
    //
    double C11() const { return C11_; }
    //
    double C22() const { return C22_; }
    //
    double C23() const { return C23_; }
    //
    double C33() const { return C33_; }
    //
    double H12() const { return r() + (dataFormats_->hybrid() ? setup_->hybridChosenRofPhi() : setup_->chosenRofPhi()) - setup_->chosenRofZ(); }
    //
    double H00() const { return r(); }
    //
    double m0() const { return stub_->phi(); }
    //
    double m1() const { return stub_->z(); }
    //
    double dPhi() const { return stub_->dPhi(); }
    //
    double dZ() const { return stub_->dZ(); }
    //
    double v0() const { return pow(stub_->dPhi(), 2); }
    //double v0() const { return setup_->v0(stub_->ttStubRef(), track_->qOverPt()); }
    //
    double v1() const { return pow(stub_->dZ(), 2); }
    //double v1() const { return setup_->v1(stub_->ttStubRef(), track_->cotGlobal()); }
    //
    FrameTrack frame() const {
      TrackKF track(*track_, x1_, x0_, x3_, x2_);
      //std::cout << "KF " << track.inv2R() << " " << track.phiT() << std::endl;
      return track.frame();
      }
    //
    std::vector<StubKF> stubs() const;
    //
    std::vector<StubKFin*> stubsIn() const;

  private:
    //
    const DataFormats* dataFormats_;
    //
    const trackerDTC::Setup* setup_;
    // found mht track
    TrackKFin* track_;
    //
    int trackId_;
    // previous state, nullptr for first states
    State* parent_;
    // stub to add
    StubKFin* stub_;
    // shows which stub on each layer has been added so far
    std::vector<int> layerMap_;
    // shows which layer has been added so far
    TTBV hitPattern_;
    //
    double x0_;
    //
    double x1_;
    //
    double x2_;
    //
    double x3_;
    //
    double C00_;
    //
    double C01_;
    //
    double C11_;
    //
    double C22_;
    //
    double C23_;
    //
    double C33_;
    //
    int numSkippedLayers_;
    //
    int numConsistentLayers_;
  };

}

#endif
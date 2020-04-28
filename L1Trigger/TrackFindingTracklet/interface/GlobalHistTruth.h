#ifndef L1Trigger_TrackFindingTracklet_interface_GlobalHistTruth_h
#define L1Trigger_TrackFindingTracklet_interface_GlobalHistTruth_h

#include "HistBase.h"

#include "L1Trigger/TrackFindingTracklet/interface/IMATH_TrackletCalculator.h"
#include "L1Trigger/TrackFindingTracklet/interface/IMATH_TrackletCalculatorDisk.h"
#include "L1Trigger/TrackFindingTracklet/interface/IMATH_TrackletCalculatorOverlap.h"

using namespace std;

class GlobalHistTruth {
public:

  GlobalHistTruth(Trklet::Settings* settings){
  // tracklet calculators 
  ITC_L1L2_ = new IMATH_TrackletCalculator(settings,1,2);
  ITC_L2L3_ = new IMATH_TrackletCalculator(settings,2,3);
  ITC_L3L4_ = new IMATH_TrackletCalculator(settings,3,4);
  ITC_L5L6_ = new IMATH_TrackletCalculator(settings,5,6);
  
  ITC_F1F2_ = new IMATH_TrackletCalculatorDisk(settings,1,2);
  ITC_F3F4_ = new IMATH_TrackletCalculatorDisk(settings,3,4);
  ITC_B1B2_ = new IMATH_TrackletCalculatorDisk(settings,-1,-2);
  ITC_B3B4_ = new IMATH_TrackletCalculatorDisk(settings,-3,-4);
  
  ITC_L1F1_ = new IMATH_TrackletCalculatorOverlap(settings,1,1);
  ITC_L2F1_ = new IMATH_TrackletCalculatorOverlap(settings,2,1);
  ITC_L1B1_ = new IMATH_TrackletCalculatorOverlap(settings,1,-1);
  ITC_L2B1_ = new IMATH_TrackletCalculatorOverlap(settings,2,-1);




  }
  
  static SLHCEvent*& event() {
    static SLHCEvent* theEvent = 0;
    return theEvent;
  }

  static HistBase*& histograms() {
    static HistBase dummy;
    static HistBase* theHistBase = &dummy;
    return theHistBase;
  }

  IMATH_TrackletCalculator* ITC_L1L2() { return ITC_L1L2_; }

  IMATH_TrackletCalculator* ITC_L2L3() { return ITC_L2L3_; }

  IMATH_TrackletCalculator* ITC_L3L4() { return ITC_L3L4_; }

  IMATH_TrackletCalculator* ITC_L5L6() { return ITC_L5L6_; }

  IMATH_TrackletCalculatorDisk* ITC_F1F2() { return ITC_F1F2_; }

  IMATH_TrackletCalculatorDisk* ITC_F3F4() { return ITC_F3F4_; }

  IMATH_TrackletCalculatorDisk* ITC_B1B2() { return ITC_B1B2_; }

  IMATH_TrackletCalculatorDisk* ITC_B3B4() { return ITC_B3B4_; }

  IMATH_TrackletCalculatorOverlap* ITC_L1F1() { return ITC_L1F1_; }

  IMATH_TrackletCalculatorOverlap* ITC_L1B1() { return ITC_L1B1_; }

  IMATH_TrackletCalculatorOverlap* ITC_L2F1() { return ITC_L2F1_; }

  IMATH_TrackletCalculatorOverlap* ITC_L2B1() { return ITC_L2B1_; }

  std::ofstream& ofstream(std::string fname) {
    if (ofstreams_.find(fname)!=ofstreams_.end()) {
      return *(ofstreams_[fname]);
    }
    std::ofstream* outptr=new std::ofstream(fname.c_str());
    ofstreams_[fname]=outptr;
    return *outptr;
  }
  
private:

  std::map<std::string, std::ofstream*> ofstreams_;

  // tracklet calculators 
  IMATH_TrackletCalculator* ITC_L1L2_{0};
  IMATH_TrackletCalculator* ITC_L2L3_{0};
  IMATH_TrackletCalculator* ITC_L3L4_{0};
  IMATH_TrackletCalculator* ITC_L5L6_{0};
  
  IMATH_TrackletCalculatorDisk* ITC_F1F2_{0};
  IMATH_TrackletCalculatorDisk* ITC_F3F4_{0};
  IMATH_TrackletCalculatorDisk* ITC_B1B2_{0};
  IMATH_TrackletCalculatorDisk* ITC_B3B4_{0};
  
  IMATH_TrackletCalculatorOverlap* ITC_L1F1_{0};
  IMATH_TrackletCalculatorOverlap* ITC_L2F1_{0};
  IMATH_TrackletCalculatorOverlap* ITC_L1B1_{0};
  IMATH_TrackletCalculatorOverlap* ITC_L2B1_{0};


  
};

#endif

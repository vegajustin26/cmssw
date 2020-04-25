#ifndef L1Trigger_TrackFindingTracklet_interface_GlobalHistTruth_h
#define L1Trigger_TrackFindingTracklet_interface_GlobalHistTruth_h

#include "HistBase.h"

using namespace std;

class GlobalHistTruth {
public:
  static SLHCEvent*& event() {
    static SLHCEvent* theEvent = 0;
    return theEvent;
  }

  static HistBase*& histograms() {
    static HistBase dummy;
    static HistBase* theHistBase = &dummy;
    return theHistBase;
  }

  static  IMATH_TrackletCalculator*& ITC_L1L2() {
    static IMATH_TrackletCalculator* theITC = 0;
    return theITC;
  }

  static  IMATH_TrackletCalculator*& ITC_L2L3() {
    static IMATH_TrackletCalculator* theITC = 0;
    return theITC;
  }

  static  IMATH_TrackletCalculator*& ITC_L3L4() {
    static IMATH_TrackletCalculator* theITC = 0;
    return theITC;
  }

  static  IMATH_TrackletCalculator*& ITC_L5L6() {
    static IMATH_TrackletCalculator* theITC = 0;
    return theITC;
  }

  static  IMATH_TrackletCalculatorDisk*& ITC_F1F2() {
    static IMATH_TrackletCalculatorDisk* theITC = 0;
    return theITC;
  }

  static  IMATH_TrackletCalculatorDisk*& ITC_F3F4() {
    static IMATH_TrackletCalculatorDisk* theITC = 0;
    return theITC;
  }

  static  IMATH_TrackletCalculatorDisk*& ITC_B1B2() {
    static IMATH_TrackletCalculatorDisk* theITC = 0;
    return theITC;
  }

  static  IMATH_TrackletCalculatorDisk*& ITC_B3B4() {
    static IMATH_TrackletCalculatorDisk* theITC = 0;
    return theITC;
  }

  static  IMATH_TrackletCalculatorOverlap*& ITC_L1F1() {
    static IMATH_TrackletCalculatorOverlap* theITC = 0;
    return theITC;
  }

  static  IMATH_TrackletCalculatorOverlap*& ITC_L1B1() {
    static IMATH_TrackletCalculatorOverlap* theITC = 0;
    return theITC;
  }

  static  IMATH_TrackletCalculatorOverlap*& ITC_L2F1() {
    static IMATH_TrackletCalculatorOverlap* theITC = 0;
    return theITC;
  }

  static  IMATH_TrackletCalculatorOverlap*& ITC_L2B1() {
    static IMATH_TrackletCalculatorOverlap* theITC = 0;
    return theITC;
  }




  
};

#endif

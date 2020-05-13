#ifndef L1Trigger_TrackFindingTracklet_interface_Globals_h
#define L1Trigger_TrackFindingTracklet_interface_Globals_h

#include "HistBase.h"

#include "L1Trigger/TrackFindingTracklet/interface/IMATH_TrackletCalculator.h"
#include "L1Trigger/TrackFindingTracklet/interface/IMATH_TrackletCalculatorDisk.h"
#include "L1Trigger/TrackFindingTracklet/interface/IMATH_TrackletCalculatorOverlap.h"

namespace tmtt {
  class Settings;
  class KFParamsComb;
}  // namespace tmtt

namespace trklet {

  class TETableBase;
  class TrackDerTable;
  class ProjectionRouterBendTable;
  class VMRouterPhiCorrTable;
  class SLHCEvent;

  class Globals {
  public:
    
    Globals(const Settings* settings) {
      imathGlobals_ = new imathGlobals();

      // tracklet calculators
      ITC_L1L2_ = new IMATH_TrackletCalculator(settings, imathGlobals_, 1, 2);
      ITC_L2L3_ = new IMATH_TrackletCalculator(settings, imathGlobals_, 2, 3);
      ITC_L3L4_ = new IMATH_TrackletCalculator(settings, imathGlobals_, 3, 4);
      ITC_L5L6_ = new IMATH_TrackletCalculator(settings, imathGlobals_, 5, 6);

      ITC_F1F2_ = new IMATH_TrackletCalculatorDisk(settings, imathGlobals_, 1, 2);
      ITC_F3F4_ = new IMATH_TrackletCalculatorDisk(settings, imathGlobals_, 3, 4);
      ITC_B1B2_ = new IMATH_TrackletCalculatorDisk(settings, imathGlobals_, -1, -2);
      ITC_B3B4_ = new IMATH_TrackletCalculatorDisk(settings, imathGlobals_, -3, -4);

      ITC_L1F1_ = new IMATH_TrackletCalculatorOverlap(settings, imathGlobals_, 1, 1);
      ITC_L2F1_ = new IMATH_TrackletCalculatorOverlap(settings, imathGlobals_, 2, 1);
      ITC_L1B1_ = new IMATH_TrackletCalculatorOverlap(settings, imathGlobals_, 1, -1);
      ITC_L2B1_ = new IMATH_TrackletCalculatorOverlap(settings, imathGlobals_, 2, -1);
    }

    ~Globals() {
      delete imathGlobals_;
    }
    
    SLHCEvent*& event() { return theEvent_; }

    HistBase*& histograms() { return theHistBase_; }

    TrackDerTable*& trackDerTable() { return trackDerTable_; }

    tmtt::Settings*& tmttSettings() { return tmttSettings_; }

    tmtt::KFParamsComb*& tmttKFParamsComb() { return tmttKFParamsComb_; }

    VMRouterPhiCorrTable*& phiCorr(unsigned int layer) { return thePhiCorr_[layer]; }

    TETableBase*& teTable(unsigned int inner, unsigned int iSeed) { return theTETable_[inner][iSeed]; }

    ProjectionRouterBendTable*& projectionRouterBendTable() { return projectionRouterBendTable_; }

    std::map<std::string, std::vector<int> >& ILindex() { return ILindex_; }

    std::map<std::string, int>& layerdiskmap() { return layerdiskmap_; }

    double& Vfull(int i, int j, int ptbin, int index) { return Vfull_[i][j][ptbin][index]; }

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
      if (ofstreams_.find(fname) != ofstreams_.end()) {
        return *(ofstreams_[fname]);
      }
      std::ofstream* outptr = new std::ofstream(fname.c_str());
      ofstreams_[fname] = outptr;
      return *outptr;
    }

  private:
    std::map<std::string, std::ofstream*> ofstreams_;

    // tracklet calculators
    IMATH_TrackletCalculator* ITC_L1L2_{nullptr};
    IMATH_TrackletCalculator* ITC_L2L3_{nullptr};
    IMATH_TrackletCalculator* ITC_L3L4_{nullptr};
    IMATH_TrackletCalculator* ITC_L5L6_{nullptr};

    IMATH_TrackletCalculatorDisk* ITC_F1F2_{nullptr};
    IMATH_TrackletCalculatorDisk* ITC_F3F4_{nullptr};
    IMATH_TrackletCalculatorDisk* ITC_B1B2_{nullptr};
    IMATH_TrackletCalculatorDisk* ITC_B3B4_{nullptr};

    IMATH_TrackletCalculatorOverlap* ITC_L1F1_{nullptr};
    IMATH_TrackletCalculatorOverlap* ITC_L2F1_{nullptr};
    IMATH_TrackletCalculatorOverlap* ITC_L1B1_{nullptr};
    IMATH_TrackletCalculatorOverlap* ITC_L2B1_{nullptr};

    imathGlobals* imathGlobals_{nullptr};

    SLHCEvent* theEvent_{nullptr};

    HistBase* theHistBase_{nullptr};

    TrackDerTable* trackDerTable_{nullptr};

    ProjectionRouterBendTable* projectionRouterBendTable_{nullptr};

    tmtt::Settings* tmttSettings_{nullptr};

    tmtt::KFParamsComb* tmttKFParamsComb_{nullptr};

    std::array<VMRouterPhiCorrTable*, 6> thePhiCorr_{{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}};

    std::array<std::array<TETableBase*, 12>, 3> theTETable_{{{{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}},
                                                             {{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}},
                                                             {{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}}}};

    std::map<std::string, std::vector<int> > ILindex_;

    std::map<std::string, int> layerdiskmap_;

    double Vfull_[11][11][4][1000];
  };
};  // namespace trklet

#endif

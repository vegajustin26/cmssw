#include "L1Trigger/TrackFindingTracklet/interface/TrackletCalculator.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletProjectionsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/StubPairsMemory.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;
using namespace Trklet;

TrackletCalculator::TrackletCalculator(string name,
                                       const Settings* const settings,
                                       Globals* globals,
                                       unsigned int iSector)
    : TrackletCalculatorBase(name, settings, globals, iSector) {
  phioffset_ = phimin_;

  for (unsigned int ilayer = 0; ilayer < 6; ilayer++) {
    vector<TrackletProjectionsMemory*> tmp(settings->nallstubs(ilayer), 0);
    trackletprojlayers_.push_back(tmp);
  }

  for (unsigned int idisk = 0; idisk < 5; idisk++) {
    vector<TrackletProjectionsMemory*> tmp(settings->nallstubs(idisk + 6), 0);
    trackletprojdisks_.push_back(tmp);
  }

  layer_ = 0;
  disk_ = 0;

  if (name_[3] == 'L')
    layer_ = name_[4] - '0';
  if (name_[3] == 'D')
    disk_ = name_[4] - '0';

  // set TC index
  iTC_=name_[7]-'A';

  if (name_.substr(3, 4) == "L1L2")
    iSeed_ = 0;
  else if (name_.substr(3, 4) == "L3L4")
    iSeed_ = 2;
  else if (name_.substr(3, 4) == "L5L6")
    iSeed_ = 3;
  else if (name_.substr(3, 4) == "D1D2")
    iSeed_ = 4;
  else if (name_.substr(3, 4) == "D3D4")
    iSeed_ = 5;
  else if (name_.substr(3, 4) == "D1L1")
    iSeed_ = 6;
  else if (name_.substr(3, 4) == "D1L2")
    iSeed_ = 7;
  else if (name_.substr(3, 4) == "L1D1")
    iSeed_ = 6;
  else if (name_.substr(3, 4) == "L2D1")
    iSeed_ = 7;
  else if (name_.substr(3, 4) == "L2L3")
    iSeed_ = 1;

  assert(iSeed_ != -1);

  TCIndex_ = (iSeed_ << 4) + iTC_;
  assert(TCIndex_ >= 0 && TCIndex_ < 128);

  assert((layer_ != 0) || (disk_ != 0));

  if (iSeed_ == 0 || iSeed_ == 1 || iSeed_ == 2 || iSeed_ == 3) {
    if (layer_ == 1) {
      lproj_[0] = 3;
      lproj_[1] = 4;
      lproj_[2] = 5;
      lproj_[3] = 6;
    }

    if (layer_ == 2) {
      lproj_[0] = 1;
      lproj_[1] = 4;
      lproj_[2] = 5;
      lproj_[3] = 6;
    }

    if (layer_ == 3) {
      lproj_[0] = 1;
      lproj_[1] = 2;
      lproj_[2] = 5;
      lproj_[3] = 6;
    }

    if (layer_ == 5) {
      lproj_[0] = 1;
      lproj_[1] = 2;
      lproj_[2] = 3;
      lproj_[3] = 4;
    }
  }

  if (iSeed_ == 4 || iSeed_ == 5) {
    if (disk_ == 1) {
      dproj_[0] = 3;
      dproj_[1] = 4;
      dproj_[2] = 5;
    }

    if (disk_ == 3) {
      dproj_[0] = 1;
      dproj_[1] = 2;
      dproj_[2] = 5;
    }
  }

  if (settings_->usephicritapprox()) {
    double phicritFactor =
        0.5 * settings_->rcrit() * globals_->ITC_L1L2()->rinv_final.get_K() / globals_->ITC_L1L2()->phi0_final.get_K();
    if (std::abs(phicritFactor - 2.) > 0.25)
      edm::LogPrint("Tracklet")
          << "TrackletCalculator::TrackletCalculator phicrit approximation may be invalid! Please check.";
  }
}

void TrackletCalculator::addOutputProjection(TrackletProjectionsMemory*& outputProj, MemoryBase* memory) {
  outputProj = dynamic_cast<TrackletProjectionsMemory*>(memory);
  assert(outputProj != 0);
}

void TrackletCalculator::addOutput(MemoryBase* memory, string output) {
  if (settings_->writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding output to " << memory->getName() << " to output "
                                 << output;
  }
  if (output == "trackpar") {
    TrackletParametersMemory* tmp = dynamic_cast<TrackletParametersMemory*>(memory);
    assert(tmp != 0);
    trackletpars_ = tmp;
    return;
  }

  if (output.substr(0, 7) == "projout") {
    //output is on the form 'projoutL2PHIC' or 'projoutD3PHIB'
    TrackletProjectionsMemory* tmp = dynamic_cast<TrackletProjectionsMemory*>(memory);
    assert(tmp != 0);

    unsigned int layerdisk = output[8] - '1';   //layer or disk counting from 0
    unsigned int phiregion = output[12] - 'A';  //phiregion counting from 0

    if (output[7] == 'L') {
      assert(layerdisk < 6);
      assert(phiregion < trackletprojlayers_[layerdisk].size());
      //check that phiregion not already initialized
      assert(trackletprojlayers_[layerdisk][phiregion] == 0);
      trackletprojlayers_[layerdisk][phiregion] = tmp;
      return;
    }

    if (output[7] == 'D') {
      assert(layerdisk < 5);
      assert(phiregion < trackletprojdisks_[layerdisk].size());
      //check that phiregion not already initialized
      assert(trackletprojdisks_[layerdisk][phiregion] == 0);
      trackletprojdisks_[layerdisk][phiregion] = tmp;
      return;
    }
  }

  edm::LogProblem("Tracklet") << "Could not find output : " << output;
  assert(0);
}

void TrackletCalculator::addInput(MemoryBase* memory, string input) {
  if (settings_->writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding input from " << memory->getName() << " to input "
                                 << input;
  }
  if (input == "innerallstubin") {
    AllStubsMemory* tmp = dynamic_cast<AllStubsMemory*>(memory);
    assert(tmp != 0);
    innerallstubs_.push_back(tmp);
    return;
  }
  if (input == "outerallstubin") {
    AllStubsMemory* tmp = dynamic_cast<AllStubsMemory*>(memory);
    assert(tmp != 0);
    outerallstubs_.push_back(tmp);
    return;
  }
  if (input.substr(0, 8) == "stubpair") {
    StubPairsMemory* tmp = dynamic_cast<StubPairsMemory*>(memory);
    assert(tmp != 0);
    stubpairs_.push_back(tmp);
    return;
  }
  assert(0);
}

void TrackletCalculator::execute() {
  unsigned int countall = 0;
  unsigned int countsel = 0;

  for (unsigned int l = 0; l < stubpairs_.size(); l++) {
    if (trackletpars_->nTracklets() >= settings_->ntrackletmax()) {
      edm::LogVerbatim("Tracklet") << "Will break on too many tracklets in " << getName();
      break;
    }
    for (unsigned int i = 0; i < stubpairs_[l]->nStubPairs(); i++) {
      countall++;
      L1TStub* innerStub = stubpairs_[l]->getL1TStub1(i);
      Stub* innerFPGAStub = stubpairs_[l]->getFPGAStub1(i);

      L1TStub* outerStub = stubpairs_[l]->getL1TStub2(i);
      Stub* outerFPGAStub = stubpairs_[l]->getFPGAStub2(i);

      if (settings_->debugTracklet()) {
        edm::LogVerbatim("Tracklet") << "TrackletCalculator execute " << getName() << "[" << iSector_ << "]";
      }

      if (innerFPGAStub->isBarrel() && (getName() != "TC_D1L2A" && getName() != "TC_D1L2B")) {
        if (outerFPGAStub->isDisk()) {
          //overlap seeding
          bool accept = overlapSeeding(outerFPGAStub, outerStub, innerFPGAStub, innerStub);
          if (accept)
            countsel++;
        } else {
          //barrel+barrel seeding
          bool accept = barrelSeeding(innerFPGAStub, innerStub, outerFPGAStub, outerStub);
          if (accept)
            countsel++;
        }
      } else {
        if (outerFPGAStub->isDisk()) {
          //disk+disk seeding
          bool accept = diskSeeding(innerFPGAStub, innerStub, outerFPGAStub, outerStub);
          if (accept)
            countsel++;
        } else if (innerFPGAStub->isDisk()) {
          //layer+disk seeding
          bool accept = overlapSeeding(innerFPGAStub, innerStub, outerFPGAStub, outerStub);
          if (accept)
            countsel++;
        } else {
          assert(0);
        }
      }

      if (trackletpars_->nTracklets() >= settings_->ntrackletmax()) {
        edm::LogVerbatim("Tracklet") << "Will break on number of tracklets in " << getName();
        break;
      }

      if (countall >= settings_->maxStep("TC")) {
        if (settings_->debugTracklet())
          edm::LogVerbatim("Tracklet") << "Will break on MAXTC 1";
        break;
      }
      if (settings_->debugTracklet()) {
        edm::LogVerbatim("Tracklet") << "TrackletCalculator execute done";
      }
    }
    if (countall >= settings_->maxStep("TC")) {
      if (settings_->debugTracklet())
        edm::LogVerbatim("Tracklet") << "Will break on MAXTC 2";
      break;
    }
  }

  if (settings_->writeMonitorData("TC")) {
    globals_->ofstream("trackletcalculator.txt") << getName() << " " << countall << " " << countsel << endl;
  }
}

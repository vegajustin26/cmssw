
#include "L1Trigger/TrackFindingTracklet/interface/LayerProjection.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;
using namespace trklet;

void LayerProjection::init(Settings const& settings,
                           int projlayer,
                           int iphiproj,
                           int izproj,
                           int iphider,
                           int izder,
                           double phiproj,
                           double zproj,
                           double phiprojder,
                           double zprojder,
                           double phiprojapprox,
                           double zprojapprox,
                           double phiprojderapprox,
                           double zprojderapprox,
                           bool isPSseed) {
  assert(projlayer > 0);
  assert(projlayer <= N_LAYER);

  valid_ = true;

  projlayer_ = projlayer;
  unsigned int layerdisk = projlayer - 1;

  fpgaphiproj_.set(iphiproj, settings.nphibitsstub(layerdisk), true, __LINE__, __FILE__);
  fpgazproj_.set(izproj, settings.nzbitsstub(layerdisk), false, __LINE__, __FILE__);

  if (layerdisk < N_PSLAYER) {
    fpgaphiprojder_.set(iphider, settings.nbitsphiprojderL123(), false, __LINE__, __FILE__);
    fpgazprojder_.set(izder, settings.nbitszprojderL123(), false, __LINE__, __FILE__);
  } else {
    fpgaphiprojder_.set(iphider, settings.nbitsphiprojderL456(), false, __LINE__, __FILE__);
    fpgazprojder_.set(izder, settings.nbitszprojderL456(), false, __LINE__, __FILE__);
  }

  ////Separate the vm projections into zbins
  ////This determines the central bin:
  ////int zbin=4+(zproj.value()>>(zproj.nbits()-3));
  ////But we need some range (particularly for L5L6 seed projecting to L1-L3):
  int offset = isPSseed ? 1 : 4;

  int ztemp = fpgazproj_.value() >> (fpgazproj_.nbits() - settings.MEBinsBits() - NFINERZBITS);
  unsigned int zbin1 = (1 << (settings.MEBinsBits() - 1)) + ((ztemp - offset) >> NFINERZBITS);
  unsigned int zbin2 = (1 << (settings.MEBinsBits() - 1)) + ((ztemp + offset) >> NFINERZBITS);

  if (zbin1 >= settings.MEBins()) {
    zbin1 = 0;  //note that zbin1 is unsigned
  }
  if (zbin2 >= settings.MEBins()) {
    zbin2 = settings.MEBins() - 1;
  }

  assert(zbin1 <= zbin2);
  assert(zbin2 - zbin1 <= 1);

  fpgazbin1projvm_.set(zbin1, settings.MEBinsBits(), true, __LINE__, __FILE__);  // first z bin

  int nextbin = zbin1 != zbin2;
  fpgazbin2projvm_.set(nextbin, 1, true, __LINE__, __FILE__);  // need to check adjacent z bin?

  //fine vm z bits. Use 4 bits for fine position. starting at zbin 1
  int finez = ((1 << (settings.MEBinsBits() + NFINERZBITS - 1)) + ztemp) - (zbin1 << NFINERZBITS);

  fpgafinezvm_.set(finez, NFINERZBITS + 1, true, __LINE__, __FILE__);  // fine z postions starting at zbin1

  //fine phi bits
  int projfinephi =
      (fpgaphiproj_.value() >>
       (fpgaphiproj_.nbits() - (settings.nbitsallstubs(projlayer_) + settings.nbitsvmme(projlayer_) + NFINEPHIBITS))) &
      ((1 << NFINEPHIBITS) - 1);
  fpgafinephivm_.set(projfinephi, NFINEPHIBITS, true, __LINE__, __FILE__);  // fine phi postions

  phiproj_ = phiproj;
  zproj_ = zproj;
  phiprojder_ = phiprojder;
  zprojder_ = zprojder;

  phiprojapprox_ = phiprojapprox;
  zprojapprox_ = zprojapprox;
  phiprojderapprox_ = phiprojderapprox;
  zprojderapprox_ = zprojderapprox;
}

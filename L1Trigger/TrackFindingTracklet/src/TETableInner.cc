#include "L1Trigger/TrackFindingTracklet/interface/TETableInner.h"

using namespace std;
using namespace Trklet;

TETableInner::TETableInner(const Settings* settings) : TETableBase(settings) { nbits_ = 20; }

TETableInner::TETableInner(
    const Settings* settings, int layer1, int layer2, int layer3, int zbits, int rbits, bool thirdLayerIsDisk)
    : TETableBase(settings) {
  nbits_ = 20;
  init(settings, layer1, layer2, layer3, zbits, rbits, thirdLayerIsDisk);
}

void TETableInner::init(const Settings* settings, int layer1, int layer2, int zbits, int rbits) {
  init(settings, layer1, layer2, -1, zbits, rbits);
}

void TETableInner::init(
    const Settings* settings, int layer1, int layer2, int layer3, int zbits, int rbits, bool thirdLayerIsDisk) {
  thirdLayerIsDisk_ = thirdLayerIsDisk;

  layer1_ = layer1;
  layer2_ = layer2;
  layer3_ = layer3;
  zbits_ = zbits;
  rbits_ = rbits;

  bool extended = layer1 == 2 && layer2 == 3 && layer3_ == 1 && thirdLayerIsDisk_;
  bool extra = layer1 == 2 && layer2 == 3 && !extended;

  rbins_ = (1 << rbits);
  rminl1_ = settings->rmean(layer1 - 1) - settings->drmax();
  rmaxl1_ = settings->rmean(layer1 - 1) + settings->drmax();
  dr_ = 2 * settings->drmax() / rbins_;

  zbins_ = (1 << zbits);
  zminl1_ = -settings->zlength();
  zminl2_ = settings->zlength();
  dz_ = 2 * settings->zlength() / zbins_;

  if (layer1 == 1) {
    rmindisk_ = settings->rmindiskvm();
    rmaxdisk_ = settings->rmaxdiskl1overlapvm();
  }

  if (layer1 == 2) {
    rmindisk_ = settings->rmindiskl2overlapvm();
    rmaxdisk_ = (extended ? settings->rmaxdisk() : settings->rmaxdiskvm());
  }

  rmeanl2_ = settings->rmean(layer2 - 1);
  if (layer3_ > 0) {
    rmeanl3_ = settings->rmean(layer3 - 1);
    zmeand3_ = settings->zmean(layer3 - 1);
  } else {
    rmeanl3_ = 0.;
    zmeand3_ = 0.;
  }

  for (int izbin = 0; izbin < zbins_; izbin++) {
    for (int irbin = 0; irbin < rbins_; irbin++) {
      int value = getLookupValue(settings, izbin, irbin, extra);
      table_.push_back(value);
    }
  }

  if (settings->writeTable()) {
    writeVMTable("VMTableInnerL" + std::to_string(layer1_) + "L" + std::to_string(layer2_) + ".tab");
  }
}

int TETableInner::getLookupValue(const Settings* settings, int izbin, int irbin, bool extra) {
  double z1 = zminl1_ + izbin * dz_;
  double z2 = zminl1_ + (izbin + 1) * dz_;

  double r1 = rminl1_ + irbin * dr_;
  double r2 = rminl1_ + (irbin + 1) * dr_;

  if (extra and
      std::abs(0.5 * (z1 + z2)) < 52.0) {  //This seeding combinations should not be for central region of detector
    return -1;
  }

  double zmaxl2 = -2 * settings->zlength();
  double zminl2 = 2 * settings->zlength();

  findzL2(z1, r1, zminl2, zmaxl2);
  findzL2(z1, r2, zminl2, zmaxl2);
  findzL2(z2, r1, zminl2, zmaxl2);
  findzL2(z2, r2, zminl2, zmaxl2);

  assert(zminl2 < zmaxl2);

  if (zminl2 > settings->zlength())
    return -1;
  if (zmaxl2 < -settings->zlength())
    return -1;

  double zmaxl3 = -2 * settings->zlength();
  double zminl3 = 2 * settings->zlength();

  findzL3(z1, r1, zminl3, zmaxl3);
  findzL3(z1, r2, zminl3, zmaxl3);
  findzL3(z2, r1, zminl3, zmaxl3);
  findzL3(z2, r2, zminl3, zmaxl3);

  assert(zminl3 < zmaxl3);

  if (zminl3 > settings->zlength() && layer3_ > 0 && !thirdLayerIsDisk_)
    return -1;
  if (zmaxl3 < -settings->zlength() && layer3_ > 0 && !thirdLayerIsDisk_)
    return -1;

  int NBINS = settings->NLONGVMBINS() * settings->NLONGVMBINS();

  // first pack zbinmin and deltaz for second layer

  int zbinmin = NBINS * (zminl2 + settings->zlength()) / (2 * settings->zlength());
  int zbinmax = NBINS * (zmaxl2 + settings->zlength()) / (2 * settings->zlength());

  if (zbinmin < 0)
    zbinmin = 0;
  if (zbinmax >= NBINS)
    zbinmax = NBINS - 1;

  assert(zbinmin <= zbinmax);
  assert(zbinmax - zbinmin <= (int)settings->NLONGVMBINS());

  int valueL2 = zbinmin / 8;
  valueL2 *= 2;
  if (zbinmax / 8 - zbinmin / 8 > 0)
    valueL2 += 1;
  valueL2 *= 8;
  valueL2 += (zbinmin & 7);
  assert(valueL2 / 8 < 15);
  int deltaz = zbinmax - zbinmin;
  if (deltaz > 7) {
    deltaz = 7;
  }
  assert(deltaz < 8);
  valueL2 += (deltaz << 7);

  if (!settings->extended()) {
    return valueL2;
  }

  // then pack zbinmin and deltaz for third layer

  zbinmin = NBINS * (zminl3 + settings->zlength()) / (2 * settings->zlength());
  zbinmax = NBINS * (zmaxl3 + settings->zlength()) / (2 * settings->zlength());

  if (zbinmin < 0)
    zbinmin = 0;
  if (zbinmax >= NBINS)
    zbinmax = NBINS - 1;

  assert(zbinmin <= zbinmax);
  assert(zbinmax - zbinmin <= (int)settings->NLONGVMBINS());

  int valueL3 = zbinmin / 8;
  valueL3 *= 2;
  if (zbinmax / 8 - zbinmin / 8 > 0)
    valueL3 += 1;
  valueL3 *= 8;
  valueL3 += (zbinmin & 7);
  assert(valueL3 / 8 < 15);
  deltaz = zbinmax - zbinmin;
  if (deltaz > 7) {
    deltaz = 7;
  }
  assert(deltaz < 8);
  valueL3 += (deltaz << 7);

  int valueD3 = 0;
  if (layer3_ > 0 && thirdLayerIsDisk_) {
    if (std::abs(z1) <= settings->z0cut())
      return -1;
    if (std::abs(z2) <= settings->z0cut())
      return -1;

    double rmaxd3 = -2 * settings->rmaxdisk();
    double rmind3 = 2 * settings->rmaxdisk();

    findr(r1, z1, rmind3, rmaxd3);
    findr(r1, z2, rmind3, rmaxd3);
    findr(r2, z1, rmind3, rmaxd3);
    findr(r2, z2, rmind3, rmaxd3);

    assert(rmind3 < rmaxd3);

    if (rmind3 > rmaxdisk_)
      return -1;
    if (rmind3 < rmindisk_)
      rmind3 = rmindisk_;
    if (rmaxd3 > rmaxdisk_)
      rmaxd3 = rmaxdisk_;
    if (rmaxd3 < rmindisk_)
      return -1;

    int NBINS = settings->NLONGVMBINS() * settings->NLONGVMBINS() / 2;  //divide by two for + and - z

    int rbinmin = NBINS * (rmind3 - settings->rmindiskvm()) / (settings->rmaxdisk() - settings->rmindiskvm());
    int rbinmax = NBINS * (rmaxd3 - settings->rmindiskvm()) / (settings->rmaxdisk() - settings->rmindiskvm());

    if (rmind3 < settings->rmaxdiskvm())
      rbinmin = 0;
    if (rmaxd3 < settings->rmaxdiskvm())
      rbinmax = 0;

    if (rbinmin < 0)
      rbinmin = 0;
    if (rbinmax >= NBINS)
      rbinmax = NBINS - 1;

    assert(rbinmin <= rbinmax);

    valueD3 = rbinmin / 8;
    if (z1 < 0)
      valueD3 += 4;
    valueD3 *= 2;
    if (rbinmax / 8 - rbinmin / 8 > 0)
      valueD3 += 1;
    valueD3 *= 8;
    valueD3 += (rbinmin & 7);
    assert(valueD3 / 8 < 15);
    int deltar = rbinmax - rbinmin;
    if (deltar > 7) {
      deltar = 7;
    }
    assert(deltar < 8);
    valueD3 += (deltar << 7);
    assert(valueD3 < (1 << 10));
  }

  // mask out the values for third layer if this is not a table for a triplet seed
  if (layer3_ <= 0) {
    valueL3 = 0;
    valueD3 = 0;
  }

  // finally pack values for second and third layer together

  int value = (valueL3 << 10) + valueL2;
  if (thirdLayerIsDisk_)
    value = (valueD3 << 10) + valueL2;

  return value;
}

void TETableInner::findzL2(double z, double r, double& zminl2, double& zmaxl2) {
  double zl2 = zinterceptL2(settings_->z0cut(), z, r);

  if (zl2 < zminl2)
    zminl2 = zl2;
  if (zl2 > zmaxl2)
    zmaxl2 = zl2;

  zl2 = zinterceptL2(-settings_->z0cut(), z, r);

  if (zl2 < zminl2)
    zminl2 = zl2;
  if (zl2 > zmaxl2)
    zmaxl2 = zl2;
}

double TETableInner::zinterceptL2(double zcut, double z, double r) { return zcut + (z - zcut) * rmeanl2_ / r; }

void TETableInner::findzL3(double z, double r, double& zminl3, double& zmaxl3) {
  double zl3 = zinterceptL3(settings_->z0cut(), z, r);

  if (zl3 < zminl3)
    zminl3 = zl3;
  if (zl3 > zmaxl3)
    zmaxl3 = zl3;

  zl3 = zinterceptL3(-settings_->z0cut(), z, r);

  if (zl3 < zminl3)
    zminl3 = zl3;
  if (zl3 > zmaxl3)
    zmaxl3 = zl3;
}

double TETableInner::zinterceptL3(double zcut, double z, double r) { return zcut + (z - zcut) * rmeanl3_ / r; }

void TETableInner::findr(double r, double z, double& rmind2, double& rmaxd2) {
  double rd2 = rintercept(settings_->z0cut(), r, z);

  if (rd2 < rmind2)
    rmind2 = rd2;
  if (rd2 > rmaxd2)
    rmaxd2 = rd2;

  rd2 = rintercept(-settings_->z0cut(), r, z);

  if (rd2 < rmind2)
    rmind2 = rd2;
  if (rd2 > rmaxd2)
    rmaxd2 = rd2;
}

double TETableInner::rintercept(double zcut, double r, double z) {
  double zmean = (z > 0.0) ? zmeand3_ : -zmeand3_;
  return (zmean - zcut) * r / (z - zcut);
}

int TETableInner::lookup(int zbin, int rbin) {
  int index = zbin * rbins_ + rbin;
  return table_[index];
}

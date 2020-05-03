#include "L1Trigger/TrackFindingTracklet/interface/TETableInnerDisk.h"

using namespace std;
using namespace Trklet;

TETableInnerDisk::TETableInnerDisk(const Settings* settings) :
  TETableBase(settings) {
  nbits_ = 19;
}

TETableInnerDisk::TETableInnerDisk(const Settings* settings,int disk1, int disk2, int layer3, int zbits, int rbits) :
  TETableBase(settings) {
  nbits_ = 19;
  init(settings, disk1, disk2, layer3, zbits, rbits);
}

void TETableInnerDisk::init(const Settings* settings, int disk1, int disk2, int zbits, int rbits) {
  init(settings, disk1, disk2, -1, zbits, rbits);
}

void TETableInnerDisk::init(const Settings* settings, int disk1, int disk2, int layer3, int zbits, int rbits) {
  disk1_ = disk1;
  disk2_ = disk2;
  layer3_ = layer3;
  rbits_ = rbits;
  zbits_ = zbits;
  
  rbins_ = (1 << rbits);
  rmind1_ = 0.0;
  rmaxd1_ = settings->rmaxdisk();
  dr_ = settings->rmaxdisk() / rbins_;
  
  zbins_ = (1 << zbits);
  zmind1_ = settings->zmean(disk1 - 1) - settings->dzmax();
  zmaxd1_ = settings->zmean(disk1 - 1) + settings->dzmax();
  dz_ = 2 * settings->dzmax() / zbins_;
  
  zmeand2_ = settings->zmean(disk2 - 1);
  rmeanl3_ = 0.;
  if (layer3 > 0)
    rmeanl3_ = settings->rmean(layer3 - 1);
  
  for (int irbin = 0; irbin < rbins_; irbin++) {
    for (int izbin = 0; izbin < zbins_; izbin++) {
      int value = getLookupValue(settings, irbin, izbin);
      table_.push_back(value);
    }
  }
  if (settings->writeTable()) {
    writeVMTable("VMTableInnerD" + std::to_string(disk1_) + "D" + std::to_string(disk2_) + ".tab");
  }
}

int TETableInnerDisk::getLookupValue(const Settings* settings, int irbin, int izbin) {
  double r1 = rmind1_ + irbin * dr_;
  double r2 = rmind1_ + (irbin + 1) * dr_;
  
  double z1 = zmind1_ + izbin * dz_;
  double z2 = zmind1_ + (izbin + 1) * dz_;
  
  double rmaxd2 = -2 * settings->rmaxdisk();
  double rmind2 = 2 * settings->rmaxdisk();
  
  findr(r1, z1, rmind2, rmaxd2);
  findr(r1, z2, rmind2, rmaxd2);
  findr(r2, z1, rmind2, rmaxd2);
  findr(r2, z2, rmind2, rmaxd2);
  
  assert(rmind2 < rmaxd2);
  
  if (rmind2 > settings_->rmaxdiskvm())
    return -1;
  if (rmaxd2 < settings_->rmindiskvm())
    return -1;
  
  double zmaxl3 = -2 * settings->zlength();
  double zminl3 = 2 * settings->zlength();
  
  findz(z1, r1, zminl3, zmaxl3);
  findz(z1, r2, zminl3, zmaxl3);
  findz(z2, r1, zminl3, zmaxl3);
  findz(z2, r2, zminl3, zmaxl3);
  
  assert(zminl3 < zmaxl3);
  
  if (zminl3 > settings->zlength() && layer3_ > 0)
    return -1;
  if (zmaxl3 < -settings->zlength() && layer3_ > 0)
    return -1;
  
  int NBINS = settings_->NLONGVMBINS() * settings_->NLONGVMBINS() / 2;  //divide by two for + and - z
  
  // first pack rbinmin and deltar for second disk
  
  // This is a 9 bit word:
  // xxx|yy|z|rrr
  // xxx is the delta r window
  // yy is the r bin
  // z is flag to look in next bin
  // rrr fine r bin
  // NOTE : this encoding is not efficient z is one if xxx+rrr is greater than 8
  //        and xxx is only 1,2, or 3
  //        should also reject xxx=0 as this means projection is outside range
  
  int rbinmin = NBINS * (rmind2 - settings_->rmindiskvm()) / (settings_->rmaxdiskvm() - settings_->rmindiskvm());
  int rbinmax = NBINS * (rmaxd2 - settings_->rmindiskvm()) / (settings_->rmaxdiskvm() - settings_->rmindiskvm());
  
  if (rbinmin < 0)
    rbinmin = 0;
  if (rbinmax >= NBINS)
    rbinmax = NBINS - 1;
  
  assert(rbinmin <= rbinmax);
  assert(rbinmax - rbinmin <= (int)settings_->NLONGVMBINS());
  
  int valueD2 = rbinmin / 8;
  valueD2 *= 2;
  if (rbinmax / 8 - rbinmin / 8 > 0)
    valueD2 += 1;
  valueD2 *= 8;
  valueD2 += (rbinmin & 7);
  assert(valueD2 / 8 < 7);
  int deltar = rbinmax - rbinmin;
  assert(deltar < 8);
  valueD2 += (deltar << 6);
  assert(valueD2 < (1 << 9));
  
  // then pack zbinmin and deltaz for third layer
  
  NBINS = settings_->NLONGVMBINS() * settings_->NLONGVMBINS();
  
  int zbinmin = NBINS * (zminl3 + settings->zlength()) / (2 * settings->zlength());
  int zbinmax = NBINS * (zmaxl3 + settings->zlength()) / (2 * settings->zlength());
  
  if (zbinmin < 0)
    zbinmin = 0;
  if (zbinmax >= NBINS)
    zbinmax = NBINS - 1;
  
  assert(zbinmin <= zbinmax);
  assert(zbinmax - zbinmin <= (int)settings_->NLONGVMBINS());
  
  int valueL3 = zbinmin / 8;
  valueL3 *= 2;
  if (zbinmax / 8 - zbinmin / 8 > 0)
    valueL3 += 1;
  valueL3 *= 8;
  valueL3 += (zbinmin & 7);
  assert(valueL3 / 8 < 15);
  int deltaz = zbinmax - zbinmin;
  if (deltaz > 7) {
    deltaz = 7;
  }
  assert(deltaz < 8);
  valueL3 += (deltaz << 7);
  
  // mask out the values for third layer if this is not a table for a triplet seed
  if (layer3_ <= 0)
    valueL3 = 0;
  
  // finally pack values for second and third layer together
  int value = (valueL3 << 9) + valueD2;
  
  return value;
}

void TETableInnerDisk::findr(double r, double z, double& rmind2, double& rmaxd2) {
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

double TETableInnerDisk::rintercept(double zcut, double r, double z) {
  return (zmeand2_ - zcut) * r / (z - zcut);
}

void TETableInnerDisk::findz(double z, double r, double& zminl3, double& zmaxl3) {
  double zl3 = zintercept(settings_->z0cut(), z, r);
  
  if (zl3 < zminl3)
    zminl3 = zl3;
  if (zl3 > zmaxl3)
    zmaxl3 = zl3;
  
  zl3 = zintercept(-settings_->z0cut(), z, r);
  
  if (zl3 < zminl3)
    zminl3 = zl3;
  if (zl3 > zmaxl3)
    zmaxl3 = zl3;
}

double TETableInnerDisk::zintercept(double zcut, double z, double r) {
  return zcut + (z - zcut) * rmeanl3_ / r;
}

int TETableInnerDisk::lookup(int zbin, int rbin) {
  int index = rbin * zbins_ + zbin;
  assert(index < (int)table_.size());
  return table_[index];
}

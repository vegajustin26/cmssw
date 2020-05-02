#include "L1Trigger/TrackFindingTracklet/interface/TETableOuterDisk.h"

using namespace std;
using namespace Trklet;

TETableOuterDisk::TETableOuterDisk(const Settings* settings) :
  TETableBase(settings) {
  nbits_ = 5;
}

TETableOuterDisk::TETableOuterDisk(const Settings* settings, int disk, int zbits, int rbits) :
  TETableBase(settings) {
  nbits_ = 5;
  init(settings, disk, zbits, rbits);
}

void TETableOuterDisk::init(const Settings* settings, int disk, int zbits, int rbits) {
    disk_ = disk;
    rbits_ = rbits;
    zbits_ = zbits;

    rbins_ = (1 << rbits);
    rmin_ = 0;
    rmax_ = settings->rmaxdisk();
    dr_ = settings->rmaxdisk() / rbins_;

    zbins_ = (1 << zbits);
    zmin_ = settings->zmean(disk - 1) - settings->dzmax();
    zmax_ = settings->zmean(disk - 1) + settings->dzmax();
    dz_ = 2 * settings->dzmax() / zbins_;

    zmean_ = settings->zmean(disk - 1);

    for (int irbin = 0; irbin < rbins_; irbin++) {
      for (int izbin = 0; izbin < zbins_; izbin++) {
        int value = getLookupValue(irbin, izbin);
        table_.push_back(value);
      }
    }
    if (settings->writeTable()) {
      writeVMTable("VMTableOuterD" + std::to_string(disk_) + ".tab");
    }
}

int TETableOuterDisk::getLookupValue(int irbin, int izbin) {
    double r = rmin_ + (irbin + 0.5) * dr_;
    double z = zmin_ + (izbin + 0.5) * dz_;

    double rproj = r * zmean_ / z;

    int NBINS = settings_->NLONGVMBINS() * settings_->NLONGVMBINS() / 2;  //divide by two for + vs - z disks

    int rbin = NBINS * (rproj - settings_->rmindiskvm()) / (settings_->rmaxdiskvm() - settings_->rmindiskvm());

    if (rbin < 0)
      rbin = 0;
    if (rbin >= NBINS)
      rbin = NBINS - 1;

    return rbin;
}

int TETableOuterDisk::lookup(int zbin, int rbin) {
  int index = rbin * zbins_ + zbin;
  return table_[index];
}

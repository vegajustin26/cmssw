#include "L1Trigger/TrackFindingTracklet/interface/TETableOuter.h"

using namespace std;
using namespace Trklet;

TETableOuter::TETableOuter(const Settings* settings) :
  TETableBase(settings) {
  nbits_ = 6;
}

TETableOuter::TETableOuter(const Settings* settings, int layer, int zbits, int rbits) :
  TETableBase(settings) {
  nbits_ = 6;
  init(settings, layer, zbits, rbits);
}

void TETableOuter::init(const Settings* settings,int layer, int zbits, int rbits) {
    layer_ = layer;
    zbits_ = zbits;
    rbits_ = rbits;

    rbins_ = (1 << rbits);
    rmin_ = settings->rmean(layer - 1) - settings->drmax();
    rmax_ = settings->rmean(layer - 1) + settings->drmax();
    dr_ = 2 * settings->drmax() / rbins_;

    zbins_ = (1 << zbits);
    zmin_ = -settings->zlength();
    zmax_ = settings->zlength();
    dz_ = 2 * settings->zlength() / zbins_;

    rmean_ = settings->rmean(layer - 1);

    for (int izbin = 0; izbin < zbins_; izbin++) {
      for (int irbin = 0; irbin < rbins_; irbin++) {
        int value = getLookupValue(settings, izbin, irbin);
        table_.push_back(value);
      }
    }

    if (settings->writeTable()) {
      writeVMTable("VMTableOuterL" + std::to_string(layer_) + ".tab");
    }
}

int TETableOuter::getLookupValue(const Settings* settings, int izbin, int irbin) {
    double z = zmin_ + (izbin + 0.5) * dz_;
    double r = rmin_ + (irbin + 0.5) * dr_;

    double zproj = z * rmean_ / r;

    int NBINS = settings_->NLONGVMBINS() * settings_->NLONGVMBINS();

    int zbin = NBINS * (zproj + settings->zlength()) / (2 * settings->zlength());

    if (zbin < 0)
      zbin = 0;
    if (zbin >= NBINS)
      zbin = NBINS - 1;

    return zbin;
}

int TETableOuter::lookup(int zbin, int rbin) {
  int index = zbin * rbins_ + rbin;
  return table_[index];
}

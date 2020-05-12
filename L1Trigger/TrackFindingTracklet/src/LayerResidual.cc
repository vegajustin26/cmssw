#include "L1Trigger/TrackFindingTracklet/interface/LayerResidual.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"

using namespace std;
using namespace trklet;

void LayerResidual::init(const Settings* settings,
                         int layer,
                         int iphiresid,
                         int izresid,
                         int istubid,
                         double phiresid,
                         double zresid,
                         double phiresidapprox,
                         double zresidapprox,
                         double rstub,
                         std::pair<Stub*, L1TStub*> stubptrs) {
  assert(layer >= 1);
  assert(layer <= 6);

  if (valid_ && (std::abs(iphiresid) > std::abs(fpgaphiresid_.value())))
    return;

  valid_ = true;

  layer_ = layer;

  fpgaphiresid_.set(iphiresid, settings->phiresidbits(), false, __LINE__, __FILE__);
  fpgazresid_.set(izresid, settings->zresidbits(), false, __LINE__, __FILE__);
  int nbitsid = 10;
  fpgastubid_.set(istubid, nbitsid, true, __LINE__, __FILE__);
  assert(!fpgaphiresid_.atExtreme());

  phiresid_ = phiresid;
  zresid_ = zresid;

  phiresidapprox_ = phiresidapprox;
  zresidapprox_ = zresidapprox;

  rstub_ = rstub;
  stubptrs_ = stubptrs;
}

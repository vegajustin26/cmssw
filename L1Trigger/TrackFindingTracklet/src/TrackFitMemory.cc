#include "L1Trigger/TrackFindingTracklet/interface/TrackFitMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/slhcevent.h"
#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"

using namespace std;
using namespace Trklet;

TrackFitMemory::TrackFitMemory(
    string name, const Settings* const settings, unsigned int iSector, double phimin, double phimax)
    : MemoryBase(name, settings, iSector) {
  phimin_ = phimin;
  phimax_ = phimax;
}

bool TrackFitMemory::foundTrack(ofstream& outres, L1SimTrack simtrk) {
  bool match = false;
  double phioffset = phimin_;
  for (unsigned int i = 0; i < tracks_.size(); i++) {
    if (tracks_[i]->getTrack()->duplicate())
      continue;
    if (tracks_[i]->foundTrack(simtrk, phioffset))
      match = true;
    if (tracks_[i]->foundTrack(simtrk, phioffset)) {
      Tracklet* tracklet = tracks_[i];

      int charge = simtrk.trackid() / abs(simtrk.trackid());
      if (abs(simtrk.trackid()) < 100)
        charge = -charge;
      double simphi = simtrk.phi();
      if (simphi < 0.0)
        simphi += 2 * M_PI;
      int irinv = tracklet->irinvfit().value();
      if (irinv == 0)
        irinv = 1;
      int layerordisk = -1;
      if (tracklet->isBarrel()) {
        layerordisk = tracklet->layer();
      } else {
        layerordisk = tracklet->disk();
      }
      if (settings_->writeMonitorData("ResEff")) {
        outres << layerordisk << " " << tracklet->nMatches() << " " << simtrk.pt() * charge << " " << simphi << " "
               << simtrk.eta() << " " << simtrk.vz() << "   " << (0.3 * 3.8 / 100.0) / tracklet->rinvfit() << " "
               << tracklet->phi0fit() + phioffset << " " << asinh(tracklet->tfit()) << " " << tracklet->z0fit() << "   "
               << (0.3 * 3.8 / 100.0) / tracklet->rinvfitexact() << " " << tracklet->phi0fitexact() + phioffset << " "
               << asinh(tracklet->tfitexact()) << " " << tracklet->z0fitexact() << "   "
               << (0.3 * 3.8 / 100.0) / (irinv * settings_->krinvpars()) << " "
               << tracklet->iphi0fit().value() * settings_->kphi0pars() + phioffset << " "
               << asinh(tracklet->itfit().value() * settings_->ktpars()) << " "
               << tracklet->iz0fit().value() * settings_->kz() << "   "
               << (0.3 * 3.8 / 100.0) / (1e-20 + tracklet->fpgarinv().value() * settings_->krinvpars()) << " "
               << tracklet->fpgaphi0().value() * settings_->kphi0pars() + phioffset << " "
               << asinh(tracklet->fpgat().value() * settings_->ktpars()) << " "
               << tracklet->fpgaz0().value() * settings_->kz() << "               "
               << (0.3 * 3.8 / 100.0) / (1e-20 + tracklet->rinvapprox()) << " " << tracklet->phi0approx() + phioffset
               << " " << asinh(tracklet->tapprox()) << " " << tracklet->z0approx() << endl;
      }
    }
  }
  return match;
}

void TrackFitMemory::writeTF(bool first) {
  std::string fname = "../data/MemPrints/FitTrack/TrackFit_";
  fname += getName();
  fname += "_";
  ostringstream oss;
  oss << iSector_ + 1;
  if (iSector_ + 1 < 10)
    fname += "0";
  fname += oss.str();
  fname += ".dat";
  if (first) {
    bx_ = 0;
    event_ = 1;
    out_.open(fname.c_str());
  } else
    out_.open(fname.c_str(), std::ofstream::app);

  out_ << "BX = " << (bitset<3>)bx_ << " Event : " << event_ << endl;

  for (unsigned int j = 0; j < tracks_.size(); j++) {
    out_ << "0x";
    if (j < 16)
      out_ << "0";
    out_ << hex << j << dec << " ";
    out_ << tracks_[j]->trackfitstr() << " " << Trklet::hexFormat(tracks_[j]->trackfitstr());
    out_ << "\n";
  }
  out_.close();

  bx_++;
  event_++;
  if (bx_ > 7)
    bx_ = 0;
}

#include "L1Trigger/TrackFindingTracklet/interface/L1TStub.h"

using namespace std;
using namespace Trklet;

L1TStub::L1TStub() {}

L1TStub::L1TStub(int eventid,
                 vector<int> tps,
                 int iphi,
                 int iz,
                 int layer,
                 int ladder,
                 int module,
                 int strip,
                 double x,
                 double y,
                 double z,
                 double sigmax,
                 double sigmaz,
                 double pt,
                 double bend,
                 int isPSmodule,
                 int isFlipped) {
  eventid_ = eventid;
  tps_ = tps;
  iphi_ = iphi;
  iz_ = iz;
  layer_ = layer;
  ladder_ = ladder;
  module_ = module;
  strip_ = strip;
  x_ = x;
  y_ = y;
  z_ = z;
  sigmax_ = sigmax;
  sigmaz_ = sigmaz;
  pt_ = pt;
  bend_ = bend;
  isPSmodule_ = isPSmodule;
  isFlipped_ = isFlipped;

  allstubindex_ = 999;
}

void L1TStub::AddInnerDigi(int ladder, int module, int irphi, int iz) {
  pair<int, int> tmplm(ladder, module);
  innerdigisladdermodule_.push_back(tmplm);

  pair<int, int> tmp(irphi, iz);
  innerdigis_.push_back(tmp);
}

void L1TStub::AddOuterDigi(int ladder, int module, int irphi, int iz) {
  pair<int, int> tmplm(ladder, module);
  outerdigisladdermodule_.push_back(tmplm);

  pair<int, int> tmp(irphi, iz);
  outerdigis_.push_back(tmp);
}

void L1TStub::write(ofstream& out) {
  out << "Stub: " << layer_ + 1 << "\t" << ladder_ << "\t" << module_ << "\t" << strip_ << "\t" << eventid_ << "\t"
      << pt_ << "\t" << x_ << "\t" << y_ << "\t" << z_ << "\t" << bend_ << "\t" << isPSmodule_ << "\t" << isFlipped_
      << "\t" << tps_.size() << " \t";
  for (unsigned itps = 0; itps < tps_.size(); itps++) {
    out << tps_[itps] << " \t";
  }
  out << endl;
}

void L1TStub::write(ostream& out) {
  out << "Stub: " << layer_ + 1 << "\t" << ladder_ << "\t" << module_ << "\t" << strip_ << "\t" << eventid_ << "\t"
      << pt_ << "\t" << x_ << "\t" << y_ << "\t" << z_ << "\t" << bend_ << "\t" << isPSmodule_ << "\t" << isFlipped_
      << "\t" << tps_.size() << " \t";
  for (unsigned itps = 0; itps < tps_.size(); itps++) {
    out << tps_[itps] << " \t";
  }
  out << endl;
}

int L1TStub::ptsign() {
  int ptsgn = -1.0;
  if (diphi() < iphiouter())
    ptsgn = -ptsgn;
  return ptsgn;
}

double L1TStub::diphi() {
  if (innerdigis_.size() == 0) {
    return 0.0;
  }
  double phi_tmp = 0.0;
  for (unsigned int i = 0; i < innerdigis_.size(); i++) {
    phi_tmp += innerdigis_[i].first;
  }
  return phi_tmp / innerdigis_.size();
}

double L1TStub::iphiouter() {
  if (outerdigis_.size() == 0) {
    return 0.0;
  }
  double phi_tmp = 0.0;
  for (unsigned int i = 0; i < outerdigis_.size(); i++) {
    phi_tmp += outerdigis_[i].first;
  }
  return phi_tmp / outerdigis_.size();
}

double L1TStub::diz() {
  if (innerdigis_.size() == 0) {
    return 0.0;
  }
  double z_tmp = 0.0;
  for (unsigned int i = 0; i < innerdigis_.size(); i++) {
    z_tmp += innerdigis_[i].second;
  }
  return z_tmp / innerdigis_.size();
}

bool L1TStub::operator==(const L1TStub& other) const {
  if (other.iphi() == iphi_ && other.iz() == iz_ && other.layer() == layer_ && other.ladder() == ladder_ &&
      other.module() == module_)
    return true;

  else
    return false;
}

void L1TStub::lorentzcor(double shift) {
  double r = this->r();
  double phi = this->phi() - shift / r;
  this->x_ = r * cos(phi);
  this->y_ = r * sin(phi);
}

double L1TStub::alpha() const {
  if (isPSmodule())
    return 0.0;
  int flip = 1;
  if (isFlipped())
    flip = -1;
  if (z_ > 0.0) {
    return ((int)strip_ - 509.5) * 0.009 * flip / r2();
  }
  return -((int)strip_ - 509.5) * 0.009 * flip / r2();
}

double L1TStub::alphanew() const {
  if (isPSmodule())
    return 0.0;
  int flip = 1;
  if (isFlipped())
    flip = -1;
  if (z_ > 0.0) {
    return ((int)strip_ - 509.5) * flip / 510.0;
  }
  return -((int)strip_ - 509.5) * flip / 510.0;
}

double L1TStub::alphatruncated() const {
  if (isPSmodule())
    return 0.0;
  int flip = 1;
  if (isFlipped())
    flip = -1;
  int striptruncated = strip_ / 1;
  striptruncated *= 1;
  if (z_ > 0.0) {
    return (striptruncated - 509.5) * 0.009 * flip / r2();
  }
  return -(striptruncated - 509.5) * 0.009 * flip / r2();
}

void L1TStub::setXY(double x, double y) {
  x_ = x;
  y_ = y;
}

bool L1TStub::tpmatch(int tp) const {
  for (unsigned int i = 0; i < tps_.size(); i++) {
    if (tp == tps_[i])
      return true;
  }

  return false;
}

#include "L1Trigger/TrackFindingTracklet/interface/TrackletLUT.h"
#include "L1Trigger/TrackFindingTracklet/interface/Util.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"

#include <filesystem>

using namespace std;
using namespace trklet;

TrackletLUT::TrackletLUT(const Settings& settings) : settings_(settings) {}

void TrackletLUT::initBendMatch(unsigned int layerdisk) {
  
  unsigned int nrinv = NRINVBITS;
  double rinvhalf = 0.5 * ((1 << nrinv) - 1);

  bool barrel = layerdisk < N_LAYER;
  bool isPSmodule = layerdisk < N_PSLAYER;
  double stripPitch = settings_.stripPitch(isPSmodule);
  
  if (barrel) {

    unsigned int nbits = isPSmodule ? N_BENDBITS_PS : N_BENDBITS_2S;

    for (unsigned int irinv = 0; irinv < (1u << nrinv); irinv++) {
      double rinv = (irinv - rinvhalf) * (1 << (settings_.nbitsrinv() - nrinv)) * settings_.krinvpars();

      double projbend = bendstrip(settings_.rmean(layerdisk), rinv, stripPitch);
      for (unsigned int ibend = 0; ibend < (1u << nbits); ibend++) {
        double stubbend = settings_.benddecode(ibend, layerdisk, isPSmodule);
        bool pass = std::abs(stubbend - projbend) < settings_.bendcutme(ibend, layerdisk, isPSmodule);
        table_.push_back(pass);
      }
    }
  } else {
    for (unsigned int iprojbend = 0; iprojbend < (1u << nrinv); iprojbend++) {
      double projbend = 0.5 * (iprojbend - rinvhalf);
      for (unsigned int ibend = 0; ibend < (1 << N_BENDBITS_2S); ibend++) {
        double stubbend = settings_.benddecode(ibend, layerdisk, false);
        bool pass = std::abs(stubbend - projbend) < settings_.bendcutme(ibend, layerdisk, false);
        table_.push_back(pass);
      }
    }
    for (unsigned int iprojbend = 0; iprojbend < (1u << nrinv); iprojbend++) {
      double projbend = 0.5 * (iprojbend - rinvhalf);
      for (unsigned int ibend = 0; ibend < (1 << N_BENDBITS_PS); ibend++) {
        double stubbend = settings_.benddecode(ibend, layerdisk, true);
        bool pass = std::abs(stubbend - projbend) < settings_.bendcutme(ibend, layerdisk, true);
        table_.push_back(pass);
      }
    }
  }
}


void TrackletLUT::initVMRTable(unsigned int layerdisk, VMRTableType type) {

  unsigned int zbits = settings_.vmrlutzbits(layerdisk);
  unsigned int rbits = settings_.vmrlutrbits(layerdisk);

  unsigned int rbins = (1 << rbits);
  unsigned int zbins = (1 << zbits);

  double zmin, zmax, rmin, rmax;
  
  if (layerdisk < N_LAYER) {
    zmin = -settings_.zlength();
    zmax = settings_.zlength();
    rmin = settings_.rmean(layerdisk) - settings_.drmax();
    rmax = settings_.rmean(layerdisk) + settings_.drmax();
  } else {
    rmin = 0;
    rmax = settings_.rmaxdisk();
    zmin = settings_.zmean(layerdisk - N_LAYER) - settings_.dzmax();
    zmax = settings_.zmean(layerdisk - N_LAYER) + settings_.dzmax();
  }

  double dr = (rmax - rmin) / rbins;
  double dz = (zmax - zmin) / zbins;

  int NBINS = settings_.NLONGVMBINS() * settings_.NLONGVMBINS();

  for (unsigned int izbin = 0; izbin < zbins; izbin++) {
    for (unsigned int irbin = 0; irbin < rbins; irbin++) {
      double r = rmin + (irbin + 0.5) * dr;
      double z = zmin + (izbin + 0.5) * dz;

      if (settings_.combined()) {
        int iznew = izbin - (1 << (zbits - 1));
        if (iznew < 0)
          iznew += (1 << zbits);
        assert(iznew >= 0);
        assert(iznew < (1 << zbits));
        z = zmin + (iznew + 0.5) * dz;
        if (layerdisk < N_LAYER) {
          int irnew = irbin - (1 << (rbits - 1));
          if (irnew < 0)
            irnew += (1 << rbits);
          assert(irnew >= 0);
          assert(irnew < (1 << rbits));
          r = rmin + (irnew + 0.5) * dr;
        }
      }

      if (layerdisk >= N_LAYER  && irbin < 10)  //special case for the tabulated radii in 2S disks
        r = (layerdisk < N_LAYER+2) ? settings_.rDSSinner(irbin) : settings_.rDSSouter(irbin);

      int bin;
      if (layerdisk < N_LAYER) {
        double zproj = z * settings_.rmean(layerdisk) / r;
        bin = NBINS * (zproj + settings_.zlength()) / (2 * settings_.zlength());
      } else {
        double rproj = r * settings_.zmean(layerdisk - N_LAYER) / z;
        bin = NBINS * (rproj - settings_.rmindiskvm()) / (settings_.rmaxdisk() - settings_.rmindiskvm());
      }
      if (bin < 0)
        bin = 0;
      if (bin >= NBINS)
        bin = NBINS - 1;

      if (type == VMRTableType::me ) {
	table_.push_back(bin);
      }

      if (type == VMRTableType::disk ) {
	if (layerdisk >= N_LAYER) {
	  double rproj = r * settings_.zmean(layerdisk - N_LAYER) / z;
	  bin = 0.5 * NBINS * (rproj - settings_.rmindiskvm()) / (settings_.rmaxdiskvm() - settings_.rmindiskvm());
	  //bin value of zero indicates that stub is out of range
	  if (bin < 0)
	    bin = 0;
	  if (bin >= NBINS / 2)
	    bin = 0;
	  table_.push_back(bin);
	}
      }

      if (type == VMRTableType::inner ) {
	if (layerdisk == 0 || layerdisk == 2 || layerdisk == 4 || layerdisk == 6 || layerdisk == 8) {
	  table_.push_back(getVMRLookup(layerdisk + 1, z, r, dz, dr));
	}
	if (layerdisk == 1) {
	  table_.push_back(getVMRLookup(layerdisk + 1, z, r, dz, dr, 1));
	}
      }

      if (type == VMRTableType::inneroverlap ) {
	if (layerdisk == 0 || layerdisk == 1) {
	  table_.push_back(getVMRLookup(6, z, r, dz, dr, layerdisk + 6));
	}
      }
      
      
      if (type == VMRTableType::innerthird ) {
	if (layerdisk == 1) {  //projection from L2 to D1 for L2L3D1 seeding
	  table_.push_back(getVMRLookup(6, z, r, dz, dr, 10));
	}
	
	if (layerdisk == 4) {  //projection from L5 to L4 for L5L6L4 seeding
	  table_.push_back(getVMRLookup(3, z, r, dz, dr));
	}
	
	if (layerdisk == 2) {  //projection from L3 to L5 for L3L4L2 seeding
	  table_.push_back(getVMRLookup(1, z, r, dz, dr));
	}
	
	if (layerdisk == 6) {  //projection from D1 to L2 for D1D2L2 seeding
	  table_.push_back(getVMRLookup(1, z, r, dz, dr));
	}
      }

      if (type == VMRTableType::innerthird ) {
	if (layerdisk == 0 || layerdisk == 1) {
	  table_.push_back(getVMRLookup(6, z, r, dz, dr, layerdisk + 6));
	}
      }
    }
  }
}

int TrackletLUT::getVMRLookup(unsigned int layerdisk, double z, double r, double dz, double dr, int iseed) const {
  double z0cut = settings_.z0cut();

  if (layerdisk < N_LAYER) {
    if (iseed == 1 && std::abs(z) < 52.0)
      return -1;

    double rmean = settings_.rmean(layerdisk);

    double rratio1 = rmean / (r + 0.5 * dr);
    double rratio2 = rmean / (r - 0.5 * dr);

    double z1 = (z - 0.5 * dz) * rratio1 + z0cut * (rratio1 - 1.0);
    double z2 = (z + 0.5 * dz) * rratio1 + z0cut * (rratio1 - 1.0);
    double z3 = (z - 0.5 * dz) * rratio2 + z0cut * (rratio2 - 1.0);
    double z4 = (z + 0.5 * dz) * rratio2 + z0cut * (rratio2 - 1.0);
    double z5 = (z - 0.5 * dz) * rratio1 - z0cut * (rratio1 - 1.0);
    double z6 = (z + 0.5 * dz) * rratio1 - z0cut * (rratio1 - 1.0);
    double z7 = (z - 0.5 * dz) * rratio2 - z0cut * (rratio2 - 1.0);
    double z8 = (z + 0.5 * dz) * rratio2 - z0cut * (rratio2 - 1.0);

    double zmin = std::min({z1, z2, z3, z4, z5, z6, z7, z8});
    double zmax = std::max({z1, z2, z3, z4, z5, z6, z7, z8});

    int NBINS = settings_.NLONGVMBINS() * settings_.NLONGVMBINS();

    int zbin1 = NBINS * (zmin + settings_.zlength()) / (2 * settings_.zlength());
    int zbin2 = NBINS * (zmax + settings_.zlength()) / (2 * settings_.zlength());

    if (zbin1 >= NBINS)
      return -1;
    if (zbin2 < 0)
      return -1;

    if (zbin2 >= NBINS)
      zbin2 = NBINS - 1;
    if (zbin1 < 0)
      zbin1 = 0;

    // This is a 10 bit word:
    // xxx|yyy|z|rrr
    // xxx is the delta z window
    // yyy is the z bin
    // z is flag to look in next bin
    // rrr first fine z bin
    // NOTE : this encoding is not efficient z is one if xxx+rrr is greater than 8
    //        and xxx is only 1,2, or 3
    //        should also reject xxx=0 as this means projection is outside range

    int value = zbin1 / 8;
    value *= 2;
    if (zbin2 / 8 - zbin1 / 8 > 0)
      value += 1;
    value *= 8;
    value += (zbin1 & 7);
    assert(value / 8 < 15);
    int deltaz = zbin2 - zbin1;
    if (deltaz > 7) {
      deltaz = 7;
    }
    assert(deltaz < 8);
    value += (deltaz << 7);

    return value;

  } else {
    if (std::abs(z) < 2.0 * z0cut)
      return -1;

    double zmean = settings_.zmean(layerdisk - N_LAYER);
    if (z < 0.0)
      zmean = -zmean;

    double r1 = (r + 0.5 * dr) * (zmean + z0cut) / (z + 0.5 * dz + z0cut);
    double r2 = (r - 0.5 * dr) * (zmean - z0cut) / (z + 0.5 * dz - z0cut);
    double r3 = (r + 0.5 * dr) * (zmean + z0cut) / (z - 0.5 * dz + z0cut);
    double r4 = (r - 0.5 * dr) * (zmean - z0cut) / (z - 0.5 * dz - z0cut);
    double r5 = (r + 0.5 * dr) * (zmean - z0cut) / (z + 0.5 * dz - z0cut);
    double r6 = (r - 0.5 * dr) * (zmean + z0cut) / (z + 0.5 * dz + z0cut);
    double r7 = (r + 0.5 * dr) * (zmean - z0cut) / (z - 0.5 * dz - z0cut);
    double r8 = (r - 0.5 * dr) * (zmean + z0cut) / (z - 0.5 * dz + z0cut);

    double rmin = std::min({r1, r2, r3, r4, r5, r6, r7, r8});
    double rmax = std::max({r1, r2, r3, r4, r5, r6, r7, r8});

    int NBINS = settings_.NLONGVMBINS() * settings_.NLONGVMBINS() / 2;

    double rmindisk = settings_.rmindiskvm();
    double rmaxdisk = settings_.rmaxdiskvm();

    if (iseed == 6)
      rmaxdisk = settings_.rmaxdiskl1overlapvm();
    if (iseed == 7)
      rmindisk = settings_.rmindiskl2overlapvm();
    if (iseed == 10)
      rmaxdisk = settings_.rmaxdisk();

    if (rmin > rmaxdisk)
      return -1;
    if (rmax > rmaxdisk)
      rmax = rmaxdisk;

    if (rmax < rmindisk)
      return -1;
    if (rmin < rmindisk)
      rmin = rmindisk;

    int rbin1 = NBINS * (rmin - settings_.rmindiskvm()) / (settings_.rmaxdiskvm() - settings_.rmindiskvm());
    int rbin2 = NBINS * (rmax - settings_.rmindiskvm()) / (settings_.rmaxdiskvm() - settings_.rmindiskvm());

    if (iseed == 10) {
      constexpr double rminspec = 40.0;
      rbin1 = NBINS * (rmin - rminspec) / (settings_.rmaxdisk() - rminspec);
      rbin2 = NBINS * (rmax - rminspec) / (settings_.rmaxdisk() - rminspec);
    }

    if (rbin2 >= NBINS)
      rbin2 = NBINS - 1;
    if (rbin1 < 0)
      rbin1 = 0;

    // This is a 9 bit word:
    // xxx|yy|z|rrr
    // xxx is the delta r window
    // yy is the r bin yy is three bits for overlaps
    // z is flag to look in next bin
    // rrr fine r bin
    // NOTE : this encoding is not efficient z is one if xxx+rrr is greater than 8
    //        and xxx is only 1,2, or 3
    //        should also reject xxx=0 as this means projection is outside range

    bool overlap = iseed == 6 || iseed == 7 || iseed == 10;

    int value = rbin1 / 8;
    if (overlap) {
      if (z < 0.0)
        value += 4;
    }
    value *= 2;
    if (rbin2 / 8 - rbin1 / 8 > 0)
      value += 1;
    value *= 8;
    value += (rbin1 & 7);
    assert(value / 8 < 15);
    int deltar = rbin2 - rbin1;
    if (deltar > 7)
      deltar = 7;
    if (overlap) {
      value += (deltar << 7);
    } else {
      value += (deltar << 6);
    }

    return value;
  }
}

  


void TrackletLUT::initPhiCorrTable(unsigned int layerdisk, unsigned int rbits) {

  bool psmodule = layerdisk < N_PSLAYER;

  unsigned int bendbits = psmodule ? N_BENDBITS_PS : N_BENDBITS_2S;
  
  unsigned int rbins = (1 << rbits);

  double rmean = settings_.rmean(layerdisk);
  double drmax = settings_.drmax();
  
  double dr = 2.0 * drmax / rbins;

  unsigned int bendbins = (1 << bendbits);

  for (unsigned int ibend = 0; ibend < bendbins; ibend++) {
    for (unsigned int irbin = 0; irbin < rbins; irbin++) {
      int value = getphiCorrValue(layerdisk, ibend, irbin, rmean, dr, drmax);
      table_.push_back(value);
    }
  }

  name_ = "VMPhiCorrL" + std::to_string(layerdisk+1) + ".tab";
  nbits_ = 14;
  positive_ = false;
  
}

int TrackletLUT::getphiCorrValue(unsigned int layerdisk, unsigned int ibend, unsigned int irbin,
				 double rmean, double dr, double drmax) const {

  bool psmodule = layerdisk < N_PSLAYER;
  
  double bend = -settings_.benddecode(ibend, layerdisk, psmodule);

  //for the rbin - calculate the distance to the nominal layer radius
  double Delta = (irbin + 0.5) * dr - drmax;

  //calculate the phi correction - this is a somewhat approximate formula
  double dphi = (Delta / 0.18) * bend * settings_.stripPitch(psmodule) / rmean;

  double kphi =  psmodule ? settings_.kphi() : settings_.kphi1();
				  
  int idphi = dphi/kphi;

  return idphi;
}


void TrackletLUT::writeTable() const{
  // Write LUT table.

  ofstream out = openfile(settings_.tablePath(), name_, __FILE__, __LINE__);

  out << "{" << endl;
  for (unsigned int i = 0; i < table_.size(); i++) {
    if (i != 0) {
      out << "," << endl;
    }

    int itable = table_[i];
    if (positive_) {
      if (table_[i] < 0) {
        itable = (1 << nbits_) - 1;
      }
    }

    out << itable;
  }
  out << endl << "};" << endl;
  out.close();
}

int TrackletLUT::lookup(unsigned int index) const {
  assert(index < table_.size());
  return table_[index];
}

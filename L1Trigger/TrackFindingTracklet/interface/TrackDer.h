#ifndef L1Trigger_TrackFindingTracklet_interface_TrackDer_h
#define L1Trigger_TrackFindingTracklet_interface_TrackDer_h

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <vector>

namespace trklet {

  class TrackDer {
  public:
    TrackDer();

    ~TrackDer() = default;

    void setIndex(int layermask, int diskmask, int alphamask, int irinv);

    int getLayerMask() const { return layermask_; }
    int getDiskMask() const { return diskmask_; }
    int getAlphaMask() const { return alphamask_; }
    int getirinv() const { return irinv_; }

    void setirinvdphi(int i, int irinvdphi) { irinvdphi_[i] = irinvdphi; }
    void setirinvdzordr(int i, int irinvdzordr) { irinvdzordr_[i] = irinvdzordr; }
    void setiphi0dphi(int i, int iphi0dphi) { iphi0dphi_[i] = iphi0dphi; }
    void setiphi0dzordr(int i, int iphi0dzordr) { iphi0dzordr_[i] = iphi0dzordr; }
    void setitdphi(int i, int itdphi) { itdphi_[i] = itdphi; }
    void setitdzordr(int i, int itdzordr) { itdzordr_[i] = itdzordr; }
    void setiz0dphi(int i, int iz0dphi) { iz0dphi_[i] = iz0dphi; }
    void setiz0dzordr(int i, int iz0dzordr) { iz0dzordr_[i] = iz0dzordr; }

    void setitdzcorr(int i, int j, int itdzcorr) { itdzcorr_[i][j] = itdzcorr; }
    void setiz0dzcorr(int i, int j, int iz0dzcorr) { iz0dzcorr_[i][j] = iz0dzcorr; }

    void setrinvdphi(int i, double rinvdphi) { rinvdphi_[i] = rinvdphi; }
    void setrinvdzordr(int i, double rinvdzordr) { rinvdzordr_[i] = rinvdzordr; }
    void setphi0dphi(int i, double phi0dphi) { phi0dphi_[i] = phi0dphi; }
    void setphi0dzordr(int i, double phi0dzordr) { phi0dzordr_[i] = phi0dzordr; }
    void settdphi(int i, double tdphi) { tdphi_[i] = tdphi; }
    void settdzordr(int i, double tdzordr) { tdzordr_[i] = tdzordr; }
    void setz0dphi(int i, double z0dphi) { z0dphi_[i] = z0dphi; }
    void setz0dzordr(int i, double z0dzordr) { z0dzordr_[i] = z0dzordr; }

    void settdzcorr(int i, int j, double tdzcorr) { tdzcorr_[i][j] = tdzcorr; }
    void setz0dzcorr(int i, int j, double z0dzcorr) { z0dzcorr_[i][j] = z0dzcorr; }

    double getrinvdphi(int i) const { return rinvdphi_[i]; }
    double getrinvdzordr(int i) const { return rinvdzordr_[i]; }
    double getphi0dphi(int i) const { return phi0dphi_[i]; }
    double getphi0dzordr(int i) const { return phi0dzordr_[i]; }
    double gettdphi(int i) const { return tdphi_[i]; }
    double gettdzordr(int i) const { return tdzordr_[i]; }
    double getz0dphi(int i) const { return z0dphi_[i]; }
    double getz0dzordr(int i) const { return z0dzordr_[i]; }

    double gettdzcorr(int i, int j) const { return tdzcorr_[i][j]; }
    double getz0dzcorr(int i, int j) const { return z0dzcorr_[i][j]; }

    double getirinvdphi(int i) const { return irinvdphi_[i]; }
    double getirinvdzordr(int i) const { return irinvdzordr_[i]; }
    double getiphi0dphi(int i) const { return iphi0dphi_[i]; }
    double getiphi0dzordr(int i) const { return iphi0dzordr_[i]; }
    double getitdphi(int i) const { return itdphi_[i]; }
    double getitdzordr(int i) const { return itdzordr_[i]; }
    double getiz0dphi(int i) const { return iz0dphi_[i]; }
    double getiz0dzordr(int i) const { return iz0dzordr_[i]; }

    int getitdzcorr(int i, int j) const { return itdzcorr_[i][j]; }
    int getiz0dzcorr(int i, int j) const { return iz0dzcorr_[i][j]; }

    void sett(double t) { t_ = t; }
    double gett() const { return t_; }

    void fill(int t, double MinvDt[4][12], int iMinvDt[4][12]) const;

  private:
    int irinvdphi_[6];
    int irinvdzordr_[6];
    int iphi0dphi_[6];
    int iphi0dzordr_[6];
    int itdphi_[6];
    int itdzordr_[6];
    int iz0dphi_[6];
    int iz0dzordr_[6];

    int itdzcorr_[3][3];
    int iz0dzcorr_[3][3];

    double rinvdphi_[6];
    double rinvdzordr_[6];
    double phi0dphi_[6];
    double phi0dzordr_[6];
    double tdphi_[6];
    double tdzordr_[6];
    double z0dphi_[6];
    double z0dzordr_[6];

    double tdzcorr_[3][3];
    double z0dzcorr_[3][3];

    double t_;

    int layermask_;
    int diskmask_;
    int alphamask_;
    int irinv_;
  };

};  // namespace trklet
#endif

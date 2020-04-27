#ifndef L1Trigger_TrackFindingTracklet_interface_TETableInnerOverlap_h
#define L1Trigger_TrackFindingTracklet_interface_TETableInnerOverlap_h

#include "TETableBase.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <vector>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;

class TETableInnerOverlap:public TETableBase{

public:

  TETableInnerOverlap(const Settings* settings) :TETableBase(settings) {
    nbits_ = 10;
  }

  TETableInnerOverlap(const Settings* settings,
		      int layer1,
		      int disk2,
		      int zbits,
		      int rbits
		      ) : TETableBase(settings) {
    nbits_ = 10;
    init(settings,layer1,disk2,zbits,rbits);
  }

  ~TETableInnerOverlap() {

  }


  void init(const Settings* settings,
	    int layer1,
	    int disk2,
	    int zbits,
	    int rbits
	    ) {

    layer1_=layer1;
    disk2_=disk2;
    zbits_=zbits;
    rbits_=rbits;

    rbins_=(1<<rbits);
    rminl1_=settings->rmean(layer1-1)-settings->drmax();
    rmaxl1_=settings->rmean(layer1-1)-settings->drmax();
    dr_=2*settings->drmax()/rbins_;

    zbins_=(1<<zbits);
    zminl1_=-settings->zlength();
    zmaxl1_=settings->zlength();
    dz_=2*settings->zlength()/zbins_;

    assert(layer1==1||layer1==2);
    
    if (layer1==1){
      rmindisk_=settings_->rmindiskvm();
      rmaxdisk_=settings_->rmaxdiskl1overlapvm();
    }

    if (layer1==2){
      rmindisk_=settings_->rmindiskl2overlapvm();
      rmaxdisk_=settings_->rmaxdiskvm();
    }

      
    zmeand2_=settings->zmean(disk2-1);

    for (int izbin=0;izbin<zbins_;izbin++) {
      for (int irbin=0;irbin<rbins_;irbin++) {
	int value=getLookupValue(settings,izbin,irbin);
	table_.push_back(value);
      }
    }

    if (settings->writeTable()) {
      writeVMTable("VMTableInnerL"+std::to_string(layer1_)+"D"+std::to_string(disk2_)+".tab");
    }
    
  }

  // negative return means that seed can not be formed
  int getLookupValue(const Settings* settings, int izbin, int irbin){

    bool print=false;
    
    double r1=rminl1_+irbin*dr_;
    double r2=rminl1_+(irbin+1)*dr_;

    double z1=zminl1_+izbin*dz_;
    double z2=zminl1_+(izbin+1)*dz_;

    if (std::abs(z1)<=settings_->z0cut()) return -1;
    if (std::abs(z2)<=settings_->z0cut()) return -1;

    double rmaxd2=-2*settings->rmaxdisk();
    double rmind2=2*settings->rmaxdisk();

    findr(r1,z1,rmind2,rmaxd2);
    findr(r1,z2,rmind2,rmaxd2);
    findr(r2,z1,rmind2,rmaxd2);
    findr(r2,z2,rmind2,rmaxd2);

    if (print) edm::LogVerbatim("Tracklet") << "PRINT layer1 rmind2 rmaxd2 z2 r1 "<<layer1_<<" "
					    <<rmind2<<" "<<rmaxd2<<" "<<z2<<" "<<r1;
    
    assert(rmind2<rmaxd2);

    if (rmind2>rmaxdisk_) return -1;
    if (rmind2<rmindisk_) rmind2=rmindisk_;
    if (rmaxd2>rmaxdisk_) rmaxd2=rmaxdisk_;
    if (rmaxd2<rmindisk_) return -1;

    int NBINS=settings_->NLONGVMBINS()*settings_->NLONGVMBINS()/2; //divide by two for + and - z
    
    int rbinmin=NBINS*(rmind2-settings_->rmindiskvm())/(settings_->rmaxdiskvm()-settings_->rmindiskvm());
    int rbinmax=NBINS*(rmaxd2-settings_->rmindiskvm())/(settings_->rmaxdiskvm()-settings_->rmindiskvm());
    
    if (rbinmin<0) rbinmin=0;
    if (rbinmax>=NBINS) rbinmax=NBINS-1;

    if (print) edm::LogVerbatim("Tracklet") << "PRINT layer1 rmind2 rmaxd2 dr z2 r1 "<<layer1_<<" "
					    <<rmind2<<" "<<rmaxd2<<" "<<" "<<(settings_->rmaxdiskvm()-settings_->rmindiskvm())/NBINS<<z2<<" "<<r1;

    if (print) edm::LogVerbatim("Tracklet") <<"PRINT rbminmin rbinmax "<<rbinmin<<" "<<rbinmax;
    
    assert(rbinmin<=rbinmax);
    //assert(rbinmax-rbinmin<=(int)settings_->NLONGVMBINS());

    int value=rbinmin/8;
    if (z1<0) value+=4;
    value*=2;
    if (rbinmax/8-rbinmin/8>0) value+=1;
    value*=8;
    value+=(rbinmin&7);
    assert(value/8<15);
    int deltar=rbinmax-rbinmin;
    if (deltar>7) {
      deltar=7;
    }
    assert(deltar<8);
    value+=(deltar<<7);
    assert(value<(1<<10));
    return value;
    
  }


  void findr(double r, double z, double& rmind2, double& rmaxd2){

    double rd2=rintercept(settings_->z0cut(),r,z);

    if (rd2<rmind2) rmind2=rd2;
    if (rd2>rmaxd2) rmaxd2=rd2;
    
    rd2=rintercept(-settings_->z0cut(),r,z);

    if (rd2<rmind2) rmind2=rd2;
    if (rd2>rmaxd2) rmaxd2=rd2;

  }

  double rintercept(double zcut, double r, double z) {

    double zmean=(z>0.0)?zmeand2_:-zmeand2_;
    
    return (zmean-zcut)*r/(z-zcut);
    
  }

  int lookup(int zbin, int rbin) {

    int index=zbin*rbins_+rbin;
    assert(index<(int)table_.size());
    return table_[index];
    
  }
    
  /*
  
  void writephi(std::string fname) {

    ofstream out(fname.c_str());

    for (int i=0;i<phitableentries_;i++){
      FPGAWord entry;
      entry.set(i,phitablebits_);
      //out << entry.str()<<" "<<tablephi_[i]<<endl;
      out <<tablephi_[i]<<endl;
    }
    out.close();
  
  }

  */


private:

  int layer1_;
  int disk2_;
  int zbits_;
  int rbits_;
  
  int rbins_;
  double rminl1_;
  double rmaxl1_;
  double dr_;

  int zbins_;
  double zminl1_;
  double zmaxl1_;
  double dz_;
  
  double zmeand2_;

  double rmaxdisk_;
  double rmindisk_;

  
};



#endif




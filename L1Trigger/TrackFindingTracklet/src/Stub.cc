#include "L1Trigger/TrackFindingTracklet/interface/Stub.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;
using namespace Trklet;


Stub::Stub(const Trklet::Settings* const settings):
  settings_(settings){
}


Stub::Stub(const L1TStub& stub,const Trklet::Settings* const settings, double phiminsec, double phimaxsec):
  settings_(settings){
  
  double r=stub.r();
  double z=stub.z();
  double sbend = stub.bend();
  
  isPSmodule_ = false;
  if (stub.isPSmodule()) isPSmodule_=true;
  
  int ibend=bendencode(sbend,isPSmodule_);
  
  int bendbits=3;
  if (!isPSmodule_) bendbits=4;
  
  bend_.set(ibend,bendbits,true,__LINE__,__FILE__);
  
  int layer = stub.layer()+1; 
  
  // hold the real values from L1Stub	
  double stubphi=stub.phi();
  
  if (layer<999) {
    
    disk_.set(0,4,false,__LINE__,__FILE__);
    
    assert(layer>=1&&layer<=6);
    double rmin = settings_->rmean(layer-1) - settings_->drmax();
    double rmax = settings_->rmean(layer-1) + settings_->drmax();
    
    if (r<rmin||r>rmax) {
      edm::LogProblem("Tracklet") << "Error r, rmin, rmeas,  rmax :"<<r<<" "<<rmin<<" "<<0.5*(rmin+rmax)<<" "<<rmax;
    }
    
    int irbits=settings_->nrbitsstub(layer-1);
    
    int ir=round_int((1<<irbits)*((r-settings_->rmean(layer-1))/(rmax-rmin)));

    
    double zmin=-settings_->zlength();
    double zmax=settings_->zlength();
    
    if (z<zmin||z>zmax) {
      edm::LogProblem("Tracklet") << "Error z, zmin, zmax :"<<z<<" "<<zmin<<" "<<zmax;
    }
    
    int izbits=settings_->nzbitsstub(layer-1);
    
    int iz=round_int((1<<izbits)*z/(zmax-zmin));
    
    if (z<zmin||z>zmax) {
      edm::LogProblem("Tracklet") << "Error z, zmin, zmax :"<<z<<" "<<zmin<<" "<<zmax;
    }
    
    assert(phimaxsec-phiminsec>0.0);
    
    if (stubphi<phiminsec-(phimaxsec-phiminsec)/6.0) {
      stubphi+=2*M_PI;
    }
    assert((phimaxsec-phiminsec)>0.0);
    
    int iphibits=settings_->nphibitsstub(layer-1);
    
    double deltaphi=Trklet::phiRange(stubphi-phiminsec);
    
    int iphi=(1<<iphibits)*deltaphi/(phimaxsec-phiminsec);
    
    layer_.set(layer-1,3,true,__LINE__,__FILE__);
    r_.set(ir,irbits,false,__LINE__,__FILE__);
    z_.set(iz,izbits,false,__LINE__,__FILE__);
    phi_.set(iphi,iphibits,true,__LINE__,__FILE__);
    
    
    phicorr_.set(iphi,iphibits,true,__LINE__,__FILE__);
    
  } else {
    
    // Here we handle the hits on disks.
    
    int disk=stub.module();
    assert(disk>=1&&disk<=5);
    int sign=1;
    if (z<0.0) sign=-1;
    
    double zmin = sign*(settings_->zmean(disk-1) - sign*settings_->dzmax());
    double zmax = sign*(settings_->zmean(disk-1) + sign*settings_->dzmax());
    
    if ((z>zmax)||(z<zmin)) {
      edm::LogProblem("Tracklet") << "Error disk z, zmax, zmin: "<<z<<" "<<zmax<<" "<<zmin;
    }
    
    int iz=(1<<settings->nzbitsstub(disk+5))*((z-sign*settings_->zmean(disk-1))/std::abs(zmax-zmin));
    
    assert(phimaxsec-phiminsec>0.0);
    if (stubphi<phiminsec-(phimaxsec-phiminsec)/6.0) {
      stubphi+=2*M_PI;
    }
    
    assert(phimaxsec-phiminsec>0.0);
    if (stubphi<phiminsec-(phimaxsec-phiminsec)/6.0) {
      stubphi+=2*M_PI;
    }
    
    int iphibits=settings_->nphibitsstub(disk+5);
    
    double deltaphi=Trklet::phiRange(stubphi-phiminsec);
    
    int iphi=(1<<iphibits)*deltaphi/(phimaxsec-phiminsec);
    
    double rmin=0;
    double rmax=settings_->rmaxdisk();
    
    if (r<rmin||r>rmax) {
      edm::LogProblem("Tracklet") << "Error disk r, rmin, rmax :"<<r<<" "<<rmin<<" "<<rmax;
    }
    
    int ir=(1<<settings_->nrbitsstub(disk+5))*(r-rmin)/(rmax-rmin);
    
    int irSS = -1;
    if (!isPSmodule_) {
      for (int i=0; i<10; ++i){
	if (disk<=2) {
	  if (std::abs(r-settings_->rDSSinner(i)) < 0.2){
	    irSS = i;
	    break;
	  }
	} else {
	  if (std::abs(r-settings_->rDSSouter(i)) < 0.2){
	    irSS = i;
	    break;
	  }
	}
      }
      if (irSS<0) {
	edm::LogProblem("Tracklet") << "WARNING! didn't find rDSS value! r = " << r << " Check that correct geometry is used!";
	assert(0);
      }
    }
    if(irSS < 0){
      //PS modules
      r_.set(ir,settings_->nrbitsstub(disk+5),true,__LINE__,__FILE__);
    }
    else {
      //SS modules
      r_.set(irSS,4,true,__LINE__,__FILE__);  // in case of SS modules, store index, not r itself
    }
    
    z_.set(iz,settings->nzbitsstub(disk+5),false,__LINE__,__FILE__);
    phi_.set(iphi,iphibits,true,__LINE__,__FILE__);
    phicorr_.set(iphi,iphibits,true,__LINE__,__FILE__);
    
    disk_.set(sign*disk,4,false,__LINE__,__FILE__);    
    
    
    double alphanew=stub.alphanew();
    assert(std::abs(alphanew)<1.0);
    int ialphanew=alphanew*(1<<(settings->nbitsalpha()-1));
    assert(ialphanew<(1<<(settings->nbitsalpha()-1)));
    assert(ialphanew>=-(1<<(settings->nbitsalpha()-1)));
    alphanew_.set(ialphanew,settings->nbitsalpha(),false,__LINE__,__FILE__);
    
  }  
}

unsigned int Stub::iphivmRaw() const {
  unsigned int iphivm=(phicorr_.value()>>(phicorr_.nbits()-5));
  assert(iphivm<32);
  return iphivm;
}


FPGAWord Stub::iphivmFineBins(int VMbits, int finebits) const {
  
  unsigned int finephi=(phicorr_.value()>>(phicorr_.nbits()-VMbits-finebits))&((1<<finebits)-1);
  return FPGAWord(finephi,finebits,true,__LINE__,__FILE__);
  
}

std::string Stub::phiregionaddressstr() {
  assert(phiregion().value()>-1);
  return phiregion().str()+stubindex_.str();	
}

FPGAWord Stub::phiregion() const {
  // 3 bits
  if (layer_.value()>=0) {
    unsigned int nallstubs=settings_->nallstubs(layer_.value());
    int iphiregion=iphivmRaw()/(32/nallstubs);
    FPGAWord phi;
    phi.set(iphiregion,3);
    return phi;
  }
  if (abs(disk_.value())>=1) {
    unsigned int nallstubs=settings_->nallstubs(abs(disk_.value())+5);
    int iphiregion=iphivmRaw()/(32/nallstubs);
    FPGAWord phi;
    phi.set(iphiregion,3);
    return phi;
  }
  assert(0);     
}


void Stub::setAllStubIndex(int nstub){
  if (nstub>=(1<<7)){
    if (settings_->debugTracklet()) edm::LogPrint("Tracklet") << "Warning too large stubindex!";
    nstub=(1<<7)-1;
  }
  
  stubindex_.set(nstub,7);
}

void Stub::setPhiCorr(int phiCorr){
  
  int iphicorr=phi_.value()-phiCorr;
  
  if (iphicorr<0) iphicorr=0;
  if (iphicorr>=(1<<phi_.nbits())) iphicorr=(1<<phi_.nbits())-1;
  
  phicorr_.set(iphicorr,phi_.nbits(),true,__LINE__,__FILE__);
  
}


double Stub::rapprox() const {
  if (disk_.value()==0){
    int lr=1<<(8-settings_->nrbitsstub(layer_.value()));
    return r_.value()*settings_->kr()*lr+settings_->rmean(layer_.value());
  }
  return r_.value()*settings_->kr();
}

double Stub::zapprox() const {
  if (disk_.value()==0){
    int lz=1;
    if (layer_.value()>=3) {
      lz=16;
    }
    return z_.value()*settings_->kz()*lz;
  }
  int sign=1;
  if (disk_.value()<0) sign=-1;
  if (sign<0) {
    return (z_.value()+1)*settings_->kz()+sign*settings_->zmean(abs(disk_.value())-1);  //FIXME Not sure why this is needed to get agreement with integer calculations
  } else {
    return z_.value()*settings_->kz()+sign*settings_->zmean(abs(disk_.value())-1);
  }
}

double Stub::phiapprox(double phimin, double) const {
  int lphi=1;
  if (layer_.value()>=3) {
    lphi=8;
  }
  return Trklet::phiRange(phimin+phi_.value()*settings_->kphi()/lphi);
}

  

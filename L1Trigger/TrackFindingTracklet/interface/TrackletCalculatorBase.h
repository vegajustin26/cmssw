//This class implementes the tracklet engine
#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletCalculatorBase_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletCalculatorBase_h

#include "ProcessBase.h"
#include "Util.h"
#include "GlobalHistTruth.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;

class TrackletCalculatorBase:public ProcessBase{

public:

 TrackletCalculatorBase(string name, const Settings* const settings, GlobalHistTruth* global, unsigned int iSector):
   ProcessBase(name,settings,global,iSector) {
  }
  

  void exacttracklet(double r1, double z1, double phi1,
		     double r2, double z2, double phi2, double,
		     double& rinv, double& phi0,
		     double& t, double& z0,
		     double phiproj[4], double zproj[4], 
		     double phider[4], double zder[4],
		     double phiprojdisk[5], double rprojdisk[5], 
		     double phiderdisk[5], double rderdisk[5]) {

    double deltaphi=Trklet::phiRange(phi1-phi2);
    
    double dist=sqrt(r2*r2+r1*r1-2*r1*r2*cos(deltaphi));
    
    rinv=2*sin(deltaphi)/dist;

    double phi1tmp=phi1-phimin_;    
     
    phi0=Trklet::phiRange(phi1tmp+asin(0.5*r1*rinv));
    
    double rhopsi1=2*asin(0.5*r1*rinv)/rinv;
	    
    double rhopsi2=2*asin(0.5*r2*rinv)/rinv;
    
    t=(z1-z2)/(rhopsi1-rhopsi2);
    
    z0=z1-t*rhopsi1;

    for (int i=0;i<4;i++) {
      exactproj(rproj_[i],rinv,phi0,t,z0,
		phiproj[i],zproj[i],phider[i],zder[i]);
    }

    for (int i=0;i<5;i++) {
      //int sign=1;
      //if (t<0) sign=-1;
      exactprojdisk(settings_->zmean(i),rinv,phi0,t,z0,
		phiprojdisk[i],rprojdisk[i],phiderdisk[i],rderdisk[i]);
    }



  }


  void exacttrackletdisk(double r1, double z1, double phi1,
			 double r2, double z2, double phi2, double,
			 double& rinv, double& phi0,
			 double& t, double& z0,
			 double phiprojLayer[3], double zprojLayer[3], 
			 double phiderLayer[3], double zderLayer[3],
			 double phiproj[3], double rproj[3], 
			 double phider[3], double rder[3]) {

    double deltaphi=Trklet::phiRange(phi1-phi2);

    double dist=sqrt(r2*r2+r1*r1-2*r1*r2*cos(deltaphi));
    
    rinv=2*sin(deltaphi)/dist;

    double phi1tmp=phi1-phimin_;    

    phi0=Trklet::phiRange(phi1tmp+asin(0.5*r1*rinv));
    
    double rhopsi1=2*asin(0.5*r1*rinv)/rinv;
	    
    double rhopsi2=2*asin(0.5*r2*rinv)/rinv;
    
    t=(z1-z2)/(rhopsi1-rhopsi2);
    
    z0=z1-t*rhopsi1;


    for (int i=0;i<3;i++) {
      exactprojdisk(zproj_[i],rinv,phi0,t,z0,
		    phiproj[i],rproj[i],
		    phider[i],rder[i]);
    }


    for (int i=0;i<3;i++) {
      exactproj(settings_->rmean(i),rinv,phi0,t,z0,
		    phiprojLayer[i],zprojLayer[i],
		    phiderLayer[i],zderLayer[i]);
    }


  }

  void exacttrackletOverlap(double r1, double z1, double phi1,
			    double r2, double z2, double phi2, double,
			    double& rinv, double& phi0,
			    double& t, double& z0,
			    double phiprojLayer[3], double zprojLayer[3], 
			    double phiderLayer[3], double zderLayer[3],
			    double phiproj[3], double rproj[3], 
			    double phider[3], double rder[3]) {

    double deltaphi=Trklet::phiRange(phi1-phi2);

    double dist=sqrt(r2*r2+r1*r1-2*r1*r2*cos(deltaphi));
    
    rinv=2*sin(deltaphi)/dist;

    if (r1>r2) rinv=-rinv;

    double phi1tmp=phi1-phimin_;    

    phi0=Trklet::phiRange(phi1tmp+asin(0.5*r1*rinv));
    
    double rhopsi1=2*asin(0.5*r1*rinv)/rinv;
	    
    double rhopsi2=2*asin(0.5*r2*rinv)/rinv;
    
    t=(z1-z2)/(rhopsi1-rhopsi2);
    
    z0=z1-t*rhopsi1;

    for (int i=0;i<4;i++) {
      exactprojdisk(zprojoverlap_[i],rinv,phi0,t,z0,
		    phiproj[i],rproj[i],
		    phider[i],rder[i]);
    }


    for (int i=0;i<1;i++) {
      exactproj(settings_->rmean(i),rinv,phi0,t,z0,
		    phiprojLayer[i],zprojLayer[i],
		    phiderLayer[i],zderLayer[i]);
    }


  }


  void exactproj(double rproj,double rinv,double phi0,
		  double t, double z0,
		  double &phiproj, double &zproj,
		  double &phider, double &zder) {

    phiproj=phi0-asin(0.5*rproj*rinv);
    zproj=z0+(2*t/rinv)*asin(0.5*rproj*rinv);

    phider=-0.5*rinv/sqrt(1-pow(0.5*rproj*rinv,2));
    zder=t/sqrt(1-pow(0.5*rproj*rinv,2));

  }


  void exactprojdisk(double zproj,double rinv,double phi0,
		     double t, double z0,
		     double &phiproj, double &rproj,
		     double &phider, double &rder) {

    if (t<0) zproj=-zproj;
    
    double tmp=rinv*(zproj-z0)/(2.0*t);
    rproj=(2.0/rinv)*sin(tmp);
    phiproj=phi0-tmp;

    phider=-rinv/(2*t);
    rder=cos(tmp)/t;

  } 


  void addDiskProj(Tracklet* tracklet, int disk){

    FPGAWord fpgar=tracklet->fpgarprojdisk(disk);

    if (fpgar.value()*settings_->krprojshiftdisk()<12.0) return;
    if (fpgar.value()*settings_->krprojshiftdisk()>112.0) return;


    
    FPGAWord fpgaphi=tracklet->fpgaphiprojdisk(disk);
    
    int iphivmRaw=fpgaphi.value()>>(fpgaphi.nbits()-5);

    int iphi=iphivmRaw/(32/settings_->nallstubs(abs(disk)+5));

    addProjectionDisk(disk,iphi,trackletprojdisks_[abs(disk)-1][iphi],tracklet);

  }


  bool addLayerProj(Tracklet* tracklet, int layer){

    assert(layer>0);

    FPGAWord fpgaz=tracklet->fpgazproj(layer);
    FPGAWord fpgaphi=tracklet->fpgaphiproj(layer);


    if(fpgaphi.atExtreme()) edm::LogProblem("Tracklet")<<"at extreme! "<<fpgaphi.value();

    assert(!fpgaphi.atExtreme());
    
    if (fpgaz.atExtreme()) return false;

    if (std::abs(fpgaz.value()*settings_->kz())>settings_->zlength()) return false;

    int iphivmRaw=fpgaphi.value()>>(fpgaphi.nbits()-5);

    int iphi=iphivmRaw/(32/settings_->nallstubs(layer-1));

    addProjection(layer,iphi,trackletprojlayers_[layer-1][iphi],tracklet);
    
    return true;

  }

  void addProjection(int layer,int iphi,TrackletProjectionsMemory* trackletprojs, Tracklet* tracklet){
    if (trackletprojs==0) {
      if (settings_->warnNoMem()) {
	edm::LogVerbatim("Tracklet") << "No projection memory exists in "<<getName()<<" for layer = "<<layer<<" iphi = "<<iphi+1;
      }
      return;
    }
    assert(trackletprojs!=0);
    trackletprojs->addProj(tracklet);
  }

  void addProjectionDisk(int disk,int iphi,TrackletProjectionsMemory* trackletprojs, Tracklet* tracklet){
    if (iSeed_==1&&abs(disk)==4) return; //L3L4 projections to D3 are not used. Should be in configuration
    if (trackletprojs==0) {
      if (layer_==3&&abs(disk)==3) return; //L3L4 projections to D3 are not used.
      if (settings_->warnNoMem()) {       
	edm::LogVerbatim("Tracklet") << "No projection memory exists in "<<getName()<<" for disk = "<<abs(disk)<<" iphi = "<<iphi+1;
      }
      return;
    }
    assert(trackletprojs!=0);
    trackletprojs->addProj(tracklet);
  }

  //Check if rinv and z0 is good
  bool goodTrackPars(bool goodrinv, bool goodz0) {
  
    bool success = true;
    if(!goodrinv){
      if (settings_->debugTracklet()) {
	edm::LogVerbatim("Tracklet") << getName()<<" TrackletCalculatorBase irinv too large";
      }
      success = false;
    }
    if (!goodz0){
      if (settings_->debugTracklet()) {
	edm::LogVerbatim("Tracklet") << getName()<<" TrackletCalculatorBase z0 cut to large";
      }
      success = false;
    }
    return success;
  }


  bool inSector(int iphi0, int irinv, double phi0approx, double rinvapprox) {
      
    double phicritapprox=phi0approx-asin(0.5*settings_->rcrit()*rinvapprox);

    int ifactor=0.5*settings_->rcrit()*settings_->krinvpars()/settings_->kphi0pars()*(1<<8);
    int iphicrit=iphi0-(irinv>>8)*ifactor;
  
    int iphicritmincut=settings_->phicritminmc()/globals_->ITC_L1L2()->phi0_final.get_K();
    int iphicritmaxcut=settings_->phicritmaxmc()/globals_->ITC_L1L2()->phi0_final.get_K(); 

    bool keepapprox=(phicritapprox>settings_->phicritminmc())&&(phicritapprox<settings_->phicritmaxmc()),
         keep=(iphicrit>iphicritmincut)&&(iphicrit<iphicritmaxcut);
    if (settings_->debugTracklet())
      if (keepapprox && !keep)
        edm::LogVerbatim("Tracklet") << getName() << " Tracklet kept with exact phicrit cut but not approximate, phicritapprox: " << phicritapprox;
    if (settings_->usephicritapprox()) {
      return keepapprox;
    } else {
      return keep;
    }

    return true;
    
  }
    

  
  
  bool barrelSeeding(Stub* innerFPGAStub, L1TStub* innerStub, Stub* outerFPGAStub, L1TStub* outerStub){
	  
    if (settings_->debugTracklet()) {
      edm::LogVerbatim("Tracklet") << "TrackletCalculator "<<getName()<<" "<<layer_<<" trying stub pair in layer (inner outer): "
				   <<innerFPGAStub->layer().value()<<" "<<outerFPGAStub->layer().value();
    }
	    
    assert(outerFPGAStub->isBarrel());
    
    assert(layer_==innerFPGAStub->layer().value()+1);
    
    assert(layer_==1||layer_==2||layer_==3||layer_==5);

    	  
    double r1=innerStub->r();
    double z1=innerStub->z();
    double phi1=innerStub->phi();
    
    double r2=outerStub->r();
    double z2=outerStub->z();
    double phi2=outerStub->phi();
    
    
    double rinv,phi0,t,z0;
    
    double phiproj[4],zproj[4],phider[4],zder[4];
    double phiprojdisk[5],rprojdisk[5],phiderdisk[5],rderdisk[5];
    
    exacttracklet(r1,z1,phi1,r2,z2,phi2,outerStub->sigmaz(),
		  rinv,phi0,t,z0,
		  phiproj,zproj,phider,zder,
		  phiprojdisk,rprojdisk,phiderdisk,rderdisk);

    if (settings_->useapprox()) {
      phi1=innerFPGAStub->phiapprox(phimin_,phimax_);
      z1=innerFPGAStub->zapprox();
      r1=innerFPGAStub->rapprox();

      phi2=outerFPGAStub->phiapprox(phimin_,phimax_);
      z2=outerFPGAStub->zapprox();
      r2=outerFPGAStub->rapprox();
    }
    
    double rinvapprox,phi0approx,tapprox,z0approx;
    double phiprojapprox[4],zprojapprox[4];
    double phiprojdiskapprox[5],rprojdiskapprox[5];
    
    IMATH_TrackletCalculator *ITC;
    if(layer_==1)      ITC = globals_->ITC_L1L2();
    else if(layer_==2) ITC = globals_->ITC_L2L3();
    else if(layer_==3) ITC = globals_->ITC_L3L4();
    else               ITC = globals_->ITC_L5L6();
    
    ITC->r1.set_fval(r1-settings_->rmean(layer_-1));
    ITC->r2.set_fval(r2-settings_->rmean(layer_));
    ITC->z1.set_fval(z1);
    ITC->z2.set_fval(z2);
    double sphi1 = phi1 - phioffset_;
    if(sphi1<0) sphi1 += 8*atan(1.);
    if(sphi1>8*atan(1.)) sphi1 -= 8*atan(1.);
    double sphi2 = phi2 - phioffset_;
    if(sphi2<0) sphi2 += 8*atan(1.);
    if(sphi2>8*atan(1.)) sphi2 -= 8*atan(1.);
    ITC->phi1.set_fval(sphi1);
    ITC->phi2.set_fval(sphi2);

    ITC->rproj0.set_fval(rproj_[0]);
    ITC->rproj1.set_fval(rproj_[1]);
    ITC->rproj2.set_fval(rproj_[2]);
    ITC->rproj3.set_fval(rproj_[3]);

    ITC->zproj0.set_fval(t>0? settings_->zmean(0) : -settings_->zmean(0));
    ITC->zproj1.set_fval(t>0? settings_->zmean(1) : -settings_->zmean(1));
    ITC->zproj2.set_fval(t>0? settings_->zmean(2) : -settings_->zmean(2));
    ITC->zproj3.set_fval(t>0? settings_->zmean(3) : -settings_->zmean(3));
    ITC->zproj4.set_fval(t>0? settings_->zmean(4) : -settings_->zmean(4));

    ITC->rinv_final.calculate();
    ITC->phi0_final.calculate();
    ITC->t_final.calculate();
    ITC->z0_final.calculate();

    ITC->phiL_0_final.calculate();
    ITC->phiL_1_final.calculate();
    ITC->phiL_2_final.calculate();
    ITC->phiL_3_final.calculate();

    ITC->zL_0_final.calculate();
    ITC->zL_1_final.calculate();
    ITC->zL_2_final.calculate();
    ITC->zL_3_final.calculate();

    ITC->phiD_0_final.calculate();
    ITC->phiD_1_final.calculate();
    ITC->phiD_2_final.calculate();
    ITC->phiD_3_final.calculate();
    ITC->phiD_4_final.calculate();

    ITC->rD_0_final.calculate();
    ITC->rD_1_final.calculate();
    ITC->rD_2_final.calculate();
    ITC->rD_3_final.calculate();
    ITC->rD_4_final.calculate();

    ITC->der_phiL_final.calculate();
    ITC->der_zL_final.calculate();
    ITC->der_phiD_final.calculate();
    ITC->der_rD_final.calculate();

    //store the approcximate results
    rinvapprox = ITC->rinv_final.get_fval();
    phi0approx = ITC->phi0_final.get_fval();
    tapprox    = ITC->t_final.get_fval();
    z0approx   = ITC->z0_final.get_fval();
        
    phiprojapprox[0] = ITC->phiL_0_final.get_fval();
    phiprojapprox[1] = ITC->phiL_1_final.get_fval();
    phiprojapprox[2] = ITC->phiL_2_final.get_fval();
    phiprojapprox[3] = ITC->phiL_3_final.get_fval();

    zprojapprox[0]   = ITC->zL_0_final.get_fval();
    zprojapprox[1]   = ITC->zL_1_final.get_fval();
    zprojapprox[2]   = ITC->zL_2_final.get_fval();
    zprojapprox[3]   = ITC->zL_3_final.get_fval();

    phiprojdiskapprox[0] = ITC->phiD_0_final.get_fval();
    phiprojdiskapprox[1] = ITC->phiD_1_final.get_fval();
    phiprojdiskapprox[2] = ITC->phiD_2_final.get_fval();
    phiprojdiskapprox[3] = ITC->phiD_3_final.get_fval();
    phiprojdiskapprox[4] = ITC->phiD_4_final.get_fval();

    rprojdiskapprox[0] = ITC->rD_0_final.get_fval();
    rprojdiskapprox[1] = ITC->rD_1_final.get_fval();
    rprojdiskapprox[2] = ITC->rD_2_final.get_fval();
    rprojdiskapprox[3] = ITC->rD_3_final.get_fval();
    rprojdiskapprox[4] = ITC->rD_4_final.get_fval();

    //now binary
    
    int irinv,iphi0,it,iz0;
    LayerProjection layerprojs[4];
    DiskProjection diskprojs[5];
    int iphiproj[4],izproj[4];
    int iphiprojdisk[5],irprojdisk[5];
    
    int ir1=innerFPGAStub->ir();
    int iphi1=innerFPGAStub->iphi();
    int iz1=innerFPGAStub->iz();
      
    int ir2=outerFPGAStub->ir();
    int iphi2=outerFPGAStub->iphi();
    int iz2=outerFPGAStub->iz();
    
    iphi1<<=(settings_->nphibitsstub(5)-settings_->nphibitsstub(layer_-1));
    iphi2<<=(settings_->nphibitsstub(5)-settings_->nphibitsstub(layer_));
    ir1<<=(8-settings_->nrbitsstub(layer_-1));
    ir2<<=(8-settings_->nrbitsstub(layer_));

    iz1<<=(settings_->nzbitsstub(0)-settings_->nzbitsstub(layer_-1));
    iz2<<=(settings_->nzbitsstub(0)-settings_->nzbitsstub(layer_));  
      
    ITC->r1.set_ival(ir1);
    ITC->r2.set_ival(ir2);
    ITC->z1.set_ival(iz1);
    ITC->z2.set_ival(iz2);
    ITC->phi1.set_ival(iphi1);
    ITC->phi2.set_ival(iphi2);
    
    ITC->rinv_final.calculate();
    ITC->phi0_final.calculate();
    ITC->t_final.calculate();
    ITC->z0_final.calculate();

    ITC->phiL_0_final.calculate();
    ITC->phiL_1_final.calculate();
    ITC->phiL_2_final.calculate();
    ITC->phiL_3_final.calculate();

    ITC->zL_0_final.calculate();
    ITC->zL_1_final.calculate();
    ITC->zL_2_final.calculate();
    ITC->zL_3_final.calculate();

    ITC->phiD_0_final.calculate();
    ITC->phiD_1_final.calculate();
    ITC->phiD_2_final.calculate();
    ITC->phiD_3_final.calculate();
    ITC->phiD_4_final.calculate();

    ITC->rD_0_final.calculate();
    ITC->rD_1_final.calculate();
    ITC->rD_2_final.calculate();
    ITC->rD_3_final.calculate();
    ITC->rD_4_final.calculate();

    ITC->der_phiL_final.calculate();
    ITC->der_zL_final.calculate();
    ITC->der_phiD_final.calculate();
    ITC->der_rD_final.calculate();

    //store the binary results
    irinv = ITC->rinv_final.get_ival();
    iphi0 = ITC->phi0_final.get_ival();
    it    = ITC->t_final.get_ival();
    iz0   = ITC->z0_final.get_ival();

    iphiproj[0] = ITC->phiL_0_final.get_ival();
    iphiproj[1] = ITC->phiL_1_final.get_ival();
    iphiproj[2] = ITC->phiL_2_final.get_ival();
    iphiproj[3] = ITC->phiL_3_final.get_ival();
    
    izproj[0]   = ITC->zL_0_final.get_ival();
    izproj[1]   = ITC->zL_1_final.get_ival();
    izproj[2]   = ITC->zL_2_final.get_ival();
    izproj[3]   = ITC->zL_3_final.get_ival();


    if (!goodTrackPars(ITC->rinv_final.local_passes(),
		       ITC->z0_final.local_passes())) return false;


    if (!inSector(iphi0,irinv,phi0approx,rinvapprox)) return false;

    for(int i=0; i<4; ++i){

      //reject projection if z is out of range
      if (izproj[i]<-(1<<(settings_->nzbitsstub(0)-1))) continue;
      if (izproj[i]>=(1<<(settings_->nzbitsstub(0)-1))) continue;

      //reject projection if phi is out of range
      if (iphiproj[i]>=(1<<settings_->nphibitsstub(5))-1) continue;
      if (iphiproj[i]<=0) continue;
      
      //Adjust bits for r and z projection depending on layer
      if (lproj_[i]<=3) {  //FIXME clean up logic
	iphiproj[i]>>=(settings_->nphibitsstub(5)-settings_->nphibitsstub(lproj_[i]-1));
      }
      else {
	izproj[i]>>=(settings_->nzbitsstub(0)-settings_->nzbitsstub(5));
      }

      layerprojs[i].init(settings_,lproj_[i],rproj_[i],
			 iphiproj[i],izproj[i],
			 ITC->der_phiL_final.get_ival(),ITC->der_zL_final.get_ival(),
			 phiproj[i],zproj[i],
			 phider[i],zder[i],
			 phiprojapprox[i],zprojapprox[i],
			 ITC->der_phiL_final.get_fval(),ITC->der_zL_final.get_fval());
      
    }

    
      
    iphiprojdisk[0] = ITC->phiD_0_final.get_ival();
    iphiprojdisk[1] = ITC->phiD_1_final.get_ival();
    iphiprojdisk[2] = ITC->phiD_2_final.get_ival();
    iphiprojdisk[3] = ITC->phiD_3_final.get_ival();
    iphiprojdisk[4] = ITC->phiD_4_final.get_ival();

    irprojdisk[0]   = ITC->rD_0_final.get_ival();
    irprojdisk[1]   = ITC->rD_1_final.get_ival();
    irprojdisk[2]   = ITC->rD_2_final.get_ival();
    irprojdisk[3]   = ITC->rD_3_final.get_ival();
    irprojdisk[4]   = ITC->rD_4_final.get_ival();

   if(std::abs(it * ITC->t_final.get_K())>1.0) {
      for(int i=0; i<5; ++i){

	if (iphiprojdisk[i]<=0) continue;
	if (iphiprojdisk[i]>=(1<<settings_->nphibitsstub(0))-1) continue;
	
	if(irprojdisk[i]< 20. / ITC->rD_0_final.get_K() ||
	   irprojdisk[i] > 120. / ITC->rD_0_final.get_K() ) continue;

	diskprojs[i].init(settings_,i+1,rproj_[i],
			   iphiprojdisk[i],irprojdisk[i],
			   ITC->der_phiD_final.get_ival(),ITC->der_rD_final.get_ival(),
			   phiprojdisk[i],rprojdisk[i],
			   phiderdisk[i],rderdisk[i],
			   phiprojdiskapprox[i],rprojdiskapprox[i],
			   ITC->der_phiD_final.get_fval(),ITC->der_rD_final.get_fval());
	
      }
    }
 
    
    if (settings_->writeMonitorData("TPars")) {
      globals_->ofstream("trackletpars.txt")  <<"Trackpars "<<layer_
					     <<"   "<<rinv<<" "<<rinvapprox<<" "<<ITC->rinv_final.get_fval()
					     <<"   "<<phi0<<" "<<phi0approx<<" "<<ITC->phi0_final.get_fval()
					     <<"   "<<t<<" "<<tapprox<<" "<<ITC->t_final.get_fval()
					     <<"   "<<z0<<" "<<z0approx<<" "<<ITC->z0_final.get_fval()
					     <<endl;
    }	        
        
    Tracklet* tracklet=new Tracklet(settings_,innerStub,NULL,outerStub,
				    innerFPGAStub,NULL,outerFPGAStub,
				    rinv,phi0,0.0,z0,t,
				    rinvapprox,phi0approx,0.0,
				    z0approx,tapprox,
				    irinv,iphi0,0,iz0,it,
				    layerprojs,
				    diskprojs,
				    false);
    
    if (settings_->debugTracklet()) {
      edm::LogVerbatim("Tracklet") << "TrackletCalculator "<<getName()<<" Found tracklet in layer = "<<layer_<<" "<<iSector_<<" phi0 = "<<phi0;
    }
        

    tracklet->setTrackletIndex(trackletpars_->nTracklets());
    tracklet->setTCIndex(TCIndex_);

    if (settings_->writeMonitorData("Seeds")) {
      ofstream fout("seeds.txt", ofstream::app);
      fout << __FILE__ << ":" << __LINE__ << " " << name_ << "_" << iSector_ << " " << tracklet->getISeed() << endl;
      fout.close();
    }
    trackletpars_->addTracklet(tracklet);

    HistBase* hists=globals_->histograms();
    int tp=tracklet->tpseed();
    hists->fillTrackletParams(settings_,globals_,iSeed_,iSector_,
			      rinvapprox,irinv*ITC->rinv_final.get_K(),
			      phi0approx,iphi0*ITC->phi0_final.get_K(),
			      asinh(tapprox),asinh(it*ITC->t_final.get_K()),
			      z0approx,iz0*ITC->z0_final.get_K(),
			      tp);

    
    bool addL3=false;
    bool addL4=false;
    bool addL5=false;
    bool addL6=false;
    for(unsigned int j=0;j<4;j++){
      bool added=false;
      if (tracklet->validProj(lproj_[j])) {
	added=addLayerProj(tracklet,lproj_[j]);
	if (added&&lproj_[j]==3) addL3=true;
	if (added&&lproj_[j]==4) addL4=true;
	if (added&&lproj_[j]==5) addL5=true;
	if (added&&lproj_[j]==6) addL6=true;
      }
    }
    
    
    for(unsigned int j=0;j<4;j++){ //no projections to 5th disk!!
      int disk=j+1;
      if (disk==4&&addL3) continue;
      if (disk==3&&addL4) continue;
      if (disk==2&&addL5) continue;
      if (disk==1&&addL6) continue;
      if (it<0) disk=-disk;
      if (tracklet->validProjDisk(abs(disk))) {
	addDiskProj(tracklet,disk);
      }
    }
    
    return true;

  }
    
  bool diskSeeding(Stub* innerFPGAStub,L1TStub* innerStub,Stub* outerFPGAStub,L1TStub* outerStub){

	    
    if (settings_->debugTracklet()) {
      edm::LogVerbatim("Tracklet") <<  "TrackletCalculator::execute calculate disk seeds";
    }
	      
    int sign=1;
    if (innerFPGAStub->disk().value()<0) sign=-1;
    
    int disk=innerFPGAStub->disk().value();
    assert(abs(disk)==1||abs(disk)==3);
    
    
    assert(innerStub->isPSmodule());
    assert(outerStub->isPSmodule());
	    
    double r1=innerStub->r();
    double z1=innerStub->z();
    double phi1=innerStub->phi();
    
    double r2=outerStub->r();
    double z2=outerStub->z();
    double phi2=outerStub->phi();
	    
    
    if (r2<r1+2.0) {
      //assert(0);
      return false; //Protection... Should be handled cleaner
      //to avoid problem with floating point 
      //calculation
    }
    
    double rinv,phi0,t,z0;
    
    double phiproj[3],zproj[3],phider[3],zder[3];
    double phiprojdisk[3],rprojdisk[3],phiderdisk[3],rderdisk[3];
    
    exacttrackletdisk(r1,z1,phi1,r2,z2,phi2,outerStub->sigmaz(),
		      rinv,phi0,t,z0,
		      phiproj,zproj,phider,zder,
		      phiprojdisk,rprojdisk,phiderdisk,rderdisk);


    //Truncates floating point positions to integer
    //representation precision
    if (settings_->useapprox()) {
      phi1=innerFPGAStub->phiapprox(phimin_,phimax_);
      z1=innerFPGAStub->zapprox();
      r1=innerFPGAStub->rapprox();
      
      phi2=outerFPGAStub->phiapprox(phimin_,phimax_);
      z2=outerFPGAStub->zapprox();
      r2=outerFPGAStub->rapprox();
    }
    
    double rinvapprox,phi0approx,tapprox,z0approx;
    double phiprojapprox[3],zprojapprox[3];
    double phiprojdiskapprox[3],rprojdiskapprox[3];
	    
    IMATH_TrackletCalculatorDisk *ITC;
    if(disk==1)       ITC = globals_->ITC_F1F2();
    else if(disk==3)  ITC = globals_->ITC_F3F4();
    else if(disk==-1) ITC = globals_->ITC_B1B2();
    else               ITC = globals_->ITC_B3B4();
    
    ITC->r1.set_fval(r1);
    ITC->r2.set_fval(r2);
    int signt = t>0? 1 : -1;
    ITC->z1.set_fval(z1-signt*settings_->zmean(abs(disk_)-1));
    ITC->z2.set_fval(z2-signt*settings_->zmean(abs(disk_)));
    double sphi1 = phi1 - phioffset_;
    if(sphi1<0) sphi1 += 8*atan(1.);
    if(sphi1>8*atan(1.)) sphi1 -= 8*atan(1.);
    double sphi2 = phi2 - phioffset_;
    if(sphi2<0) sphi2 += 8*atan(1.);
    if(sphi2>8*atan(1.)) sphi2 -= 8*atan(1.);
    ITC->phi1.set_fval(sphi1);
    ITC->phi2.set_fval(sphi2);
	    
    ITC->rproj0.set_fval(settings_->rmean(0));
    ITC->rproj1.set_fval(settings_->rmean(1));
    ITC->rproj2.set_fval(settings_->rmean(2));

    ITC->zproj0.set_fval(t>0? zproj_[0] : -zproj_[0]);
    ITC->zproj1.set_fval(t>0? zproj_[1] : -zproj_[1]);
    ITC->zproj2.set_fval(t>0? zproj_[2] : -zproj_[2]);

    ITC->rinv_final.calculate();
    ITC->phi0_final.calculate();
    ITC->t_final.calculate();
    ITC->z0_final.calculate();

    ITC->phiL_0_final.calculate();
    ITC->phiL_1_final.calculate();
    ITC->phiL_2_final.calculate();

    ITC->zL_0_final.calculate();
    ITC->zL_1_final.calculate();
    ITC->zL_2_final.calculate();

    ITC->phiD_0_final.calculate();
    ITC->phiD_1_final.calculate();
    ITC->phiD_2_final.calculate();

    ITC->rD_0_final.calculate();
    ITC->rD_1_final.calculate();
    ITC->rD_2_final.calculate();

    ITC->der_phiL_final.calculate();
    ITC->der_zL_final.calculate();
    ITC->der_phiD_final.calculate();
    ITC->der_rD_final.calculate();

    //store the approximate results
    rinvapprox = ITC->rinv_final.get_fval();
    phi0approx = ITC->phi0_final.get_fval();
    tapprox    = ITC->t_final.get_fval();
    z0approx   = ITC->z0_final.get_fval();
    
    phiprojapprox[0] = ITC->phiL_0_final.get_fval();
    phiprojapprox[1] = ITC->phiL_1_final.get_fval();
    phiprojapprox[2] = ITC->phiL_2_final.get_fval();

    zprojapprox[0]   = ITC->zL_0_final.get_fval();
    zprojapprox[1]   = ITC->zL_1_final.get_fval();
    zprojapprox[2]   = ITC->zL_2_final.get_fval();

    phiprojdiskapprox[0] = ITC->phiD_0_final.get_fval();
    phiprojdiskapprox[1] = ITC->phiD_1_final.get_fval();
    phiprojdiskapprox[2] = ITC->phiD_2_final.get_fval();

    rprojdiskapprox[0] = ITC->rD_0_final.get_fval();
    rprojdiskapprox[1] = ITC->rD_1_final.get_fval();
    rprojdiskapprox[2] = ITC->rD_2_final.get_fval();

    //now binary
    
    int irinv,iphi0,it,iz0;
    int iphiproj[3],izproj[3];
    
    int iphiprojdisk[3],irprojdisk[3];

    int ir1=innerFPGAStub->ir();
    int iphi1=innerFPGAStub->iphi();
    int iz1=innerFPGAStub->iz();
    
    int ir2=outerFPGAStub->ir();
    int iphi2=outerFPGAStub->iphi();
    int iz2=outerFPGAStub->iz();
    
    //To get same precission as for layers.
    iphi1<<=(settings_->nphibitsstub(5)-settings_->nphibitsstub(0));
    iphi2<<=(settings_->nphibitsstub(5)-settings_->nphibitsstub(0));
    
    ITC->r1.set_ival(ir1);
    ITC->r2.set_ival(ir2);
    ITC->z1.set_ival(iz1);
    ITC->z2.set_ival(iz2);
    ITC->phi1.set_ival(iphi1);
    ITC->phi2.set_ival(iphi2);

    ITC->rinv_final.calculate();
    ITC->phi0_final.calculate();
    ITC->t_final.calculate();
    ITC->z0_final.calculate();

    ITC->phiL_0_final.calculate();
    ITC->phiL_1_final.calculate();
    ITC->phiL_2_final.calculate();

    ITC->zL_0_final.calculate();
    ITC->zL_1_final.calculate();
    ITC->zL_2_final.calculate();

    ITC->phiD_0_final.calculate();
    ITC->phiD_1_final.calculate();
    ITC->phiD_2_final.calculate();

    ITC->rD_0_final.calculate();
    ITC->rD_1_final.calculate();
    ITC->rD_2_final.calculate();

    ITC->der_phiL_final.calculate();
    ITC->der_zL_final.calculate();
    ITC->der_phiD_final.calculate();
    ITC->der_rD_final.calculate();

    //store the binary results
    irinv = ITC->rinv_final.get_ival();
    iphi0 = ITC->phi0_final.get_ival();
    it    = ITC->t_final.get_ival();
    iz0   = ITC->z0_final.get_ival();

    iphiproj[0] = ITC->phiL_0_final.get_ival();
    iphiproj[1] = ITC->phiL_1_final.get_ival();
    iphiproj[2] = ITC->phiL_2_final.get_ival();
    
    izproj[0]   = ITC->zL_0_final.get_ival();
    izproj[1]   = ITC->zL_1_final.get_ival();
    izproj[2]   = ITC->zL_2_final.get_ival();

    if (!goodTrackPars(ITC->rinv_final.local_passes(),
		       ITC->z0_final.local_passes())) return false;

    if (!inSector(iphi0,irinv,phi0approx,rinvapprox)) return false;

    LayerProjection layerprojs[4];
    DiskProjection diskprojs[3];

    
    for(int i=0; i<3; ++i){

      //Check is outside z range
      if (izproj[i]<-(1<<(settings_->nzbitsstub(0)-1))) continue;
      if (izproj[i]>=(1<<(settings_->nzbitsstub(0)-1))) continue;

      //Check if outside phi range
      if (iphiproj[i]>=(1<<settings_->nphibitsstub(5))-1) continue;
      if (iphiproj[i]<=0) continue;

      //shift bits - allways in PS modules for disk seeding
      iphiproj[i]>>=(settings_->nphibitsstub(5)-settings_->nphibitsstub(0));
      
      layerprojs[i].init(settings_,i+1,settings_->rmean(i),
			 iphiproj[i],izproj[i],
			 ITC->der_phiL_final.get_ival(),ITC->der_zL_final.get_ival(),
			 phiproj[i],zproj[i],
			 phider[i],zder[i],
			 phiprojapprox[i],zprojapprox[i],
			 ITC->der_phiL_final.get_fval(),ITC->der_zL_final.get_fval());
    }

 
    iphiprojdisk[0] = ITC->phiD_0_final.get_ival();
    iphiprojdisk[1] = ITC->phiD_1_final.get_ival();
    iphiprojdisk[2] = ITC->phiD_2_final.get_ival();

    irprojdisk[0]   = ITC->rD_0_final.get_ival();
    irprojdisk[1]   = ITC->rD_1_final.get_ival();
    irprojdisk[2]   = ITC->rD_2_final.get_ival();

    

    for(int i=0; i<3; ++i){

      //check that phi projection in range
      if (iphiprojdisk[i]<=0) continue;
      if (iphiprojdisk[i]>=(1<<settings_->nphibitsstub(0))-1) continue;

      //check that r projection in range
      if(irprojdisk[i]<=0 ||
	 irprojdisk[i] > 120. / ITC->rD_0_final.get_K() ) continue;

      diskprojs[i].init(settings_,i+1,rproj_[i],
			iphiprojdisk[i],irprojdisk[i],
			ITC->der_phiD_final.get_ival(),ITC->der_rD_final.get_ival(),
			phiprojdisk[i],rprojdisk[i],
			phiderdisk[i],rderdisk[i],
			phiprojdiskapprox[i],rprojdiskapprox[i],
			ITC->der_phiD_final.get_fval(),ITC->der_rD_final.get_fval());

      
    }

    
    if (settings_->writeMonitorData("TPars")) {
      globals_->ofstream("trackletparsdisk.txt")  <<"Trackpars         "<<disk_
				      <<"   "<<rinv<<" "<<rinvapprox<<" "<<ITC->rinv_final.get_fval()
				      <<"   "<<phi0<<" "<<phi0approx<<" "<<ITC->phi0_final.get_fval()
				      <<"   "<<t<<" "<<tapprox<<" "<<ITC->t_final.get_fval()
				      <<"   "<<z0<<" "<<z0approx<<" "<<ITC->z0_final.get_fval()
				      <<endl;
    }
	    
    Tracklet* tracklet=new Tracklet(settings_,innerStub,NULL,outerStub,
				    innerFPGAStub,NULL,outerFPGAStub,
				    rinv,phi0,0.0,z0,t,
				    rinvapprox,phi0approx,0.0,
				    z0approx,tapprox,
				    irinv,iphi0,0,iz0,it,
				    layerprojs,
				    diskprojs,
				    true);
    
    if (settings_->debugTracklet()) {
      edm::LogVerbatim("Tracklet") << "Found tracklet in disk = "<<disk_<<" "<<tracklet<<" "<<iSector_;
    }
        
    tracklet->setTrackletIndex(trackletpars_->nTracklets());
    tracklet->setTCIndex(TCIndex_);

    if (settings_->writeMonitorData("Seeds")) {
      ofstream fout("seeds.txt", ofstream::app);
      fout << __FILE__ << ":" << __LINE__ << " " << name_ << "_" << iSector_ << " " << tracklet->getISeed() << endl;
      fout.close();
    }
    trackletpars_->addTracklet(tracklet);
    
    if (tracklet->validProj(1)) {
      addLayerProj(tracklet,1);
    }
    
    if (tracklet->validProj(2)) {
      addLayerProj(tracklet,2);
    }
    
    for(unsigned int j=0;j<3;j++){
      if (tracklet->validProjDisk(sign*dproj_[j])) {
	addDiskProj(tracklet,sign*dproj_[j]);
      }
    }

    return true;
    
  }
  

  bool overlapSeeding(Stub* innerFPGAStub, L1TStub* innerStub, Stub* outerFPGAStub, L1TStub* outerStub){
    
    //Deal with overlap stubs here
    assert(outerFPGAStub->isBarrel());
    
    assert(innerFPGAStub->isDisk());
    
    int disk=innerFPGAStub->disk().value();

    if (settings_->debugTracklet()) {
      edm::LogVerbatim("Tracklet") << "trying to make overlap tracklet disk_ = "<<disk_<<" "<<getName();
    }
    
    double r1=innerStub->r();
    double z1=innerStub->z();
    double phi1=innerStub->phi();
    
    double r2=outerStub->r();
    double z2=outerStub->z();
    double phi2=outerStub->phi();
    
    //Protection... Should be handled cleaner
    //to avoid problem with floating point 
    //calculation and with overflows
    //in the integer calculation
    if (r1<r2+1.5) {
      //edm::LogVerbatim("Tracklet") << "in overlap tracklet: radii wrong";
      return false;
    }
    

    double rinv,phi0,t,z0;
	    
    double phiproj[3],zproj[3],phider[3],zder[3];
    double phiprojdisk[4],rprojdisk[4],phiderdisk[4],rderdisk[4];
    
    exacttrackletOverlap(r1,z1,phi1,r2,z2,phi2,outerStub->sigmaz(),
			 rinv,phi0,t,z0,
			 phiproj,zproj,phider,zder,
			 phiprojdisk,rprojdisk,phiderdisk,rderdisk);
    
    
    //Truncates floating point positions to integer
    //representation precision
    if (settings_->useapprox()) {
      phi1=innerFPGAStub->phiapprox(phimin_,phimax_);
      z1=innerFPGAStub->zapprox();
      r1=innerFPGAStub->rapprox();
	      
      phi2=outerFPGAStub->phiapprox(phimin_,phimax_);
      z2=outerFPGAStub->zapprox();
      r2=outerFPGAStub->rapprox();
    }

    double rinvapprox,phi0approx,tapprox,z0approx;
    double phiprojapprox[3],zprojapprox[3];
    double phiprojdiskapprox[4],rprojdiskapprox[4];

    IMATH_TrackletCalculatorOverlap *ITC;
    int ll = outerFPGAStub->layer().value()+1;
    if     (ll==1 && disk==1)  ITC = globals_->ITC_L1F1();
    else if(ll==2 && disk==1)  ITC = globals_->ITC_L2F1();
    else if(ll==1 && disk==-1) ITC = globals_->ITC_L1B1();
    else if(ll==2 && disk==-1) ITC = globals_->ITC_L2B1();
    else assert(0);
    
    ITC->r1.set_fval(r2-settings_->rmean(ll-1));
    ITC->r2.set_fval(r1);
    int signt = t>0? 1 : -1;
    ITC->z1.set_fval(z2);
    ITC->z2.set_fval(z1-signt*settings_->zmean(abs(disk_)));
    double sphi1 = phi1 - phioffset_;
    if(sphi1<0) sphi1 += 8*atan(1.);
    if(sphi1>8*atan(1.)) sphi1 -= 8*atan(1.);
    double sphi2 = phi2 - phioffset_;
    if(sphi2<0) sphi2 += 8*atan(1.);
    if(sphi2>8*atan(1.)) sphi2 -= 8*atan(1.);
    ITC->phi1.set_fval(sphi2);
    ITC->phi2.set_fval(sphi1);

    ITC->rproj0.set_fval(settings_->rmean(0));
    ITC->rproj1.set_fval(settings_->rmean(1));
    ITC->rproj2.set_fval(settings_->rmean(2));
    
    ITC->zproj0.set_fval(t>0? zprojoverlap_[0] : -zprojoverlap_[0]);
    ITC->zproj1.set_fval(t>0? zprojoverlap_[1] : -zprojoverlap_[1]);
    ITC->zproj2.set_fval(t>0? zprojoverlap_[2] : -zprojoverlap_[2]);
    ITC->zproj3.set_fval(t>0? zprojoverlap_[3] : -zprojoverlap_[3]);
    
    ITC->rinv_final.calculate();
    ITC->phi0_final.calculate();
    ITC->t_final.calculate();
    ITC->z0_final.calculate();

    ITC->phiL_0_final.calculate();
    ITC->phiL_1_final.calculate();
    ITC->phiL_2_final.calculate();

    ITC->zL_0_final.calculate();
    ITC->zL_1_final.calculate();
    ITC->zL_2_final.calculate();

    ITC->phiD_0_final.calculate();
    ITC->phiD_1_final.calculate();
    ITC->phiD_2_final.calculate();
    ITC->phiD_3_final.calculate();

    ITC->rD_0_final.calculate();
    ITC->rD_1_final.calculate();
    ITC->rD_2_final.calculate();
    ITC->rD_3_final.calculate();

    ITC->der_phiL_final.calculate();
    ITC->der_zL_final.calculate();
    ITC->der_phiD_final.calculate();
    ITC->der_rD_final.calculate();

    //store the approximate results
    rinvapprox = ITC->rinv_final.get_fval();
    phi0approx = ITC->phi0_final.get_fval();
    tapprox    = ITC->t_final.get_fval();
    z0approx   = ITC->z0_final.get_fval();

    phiprojapprox[0] = ITC->phiL_0_final.get_fval();
    phiprojapprox[1] = ITC->phiL_1_final.get_fval();
    phiprojapprox[2] = ITC->phiL_2_final.get_fval();

    zprojapprox[0]   = ITC->zL_0_final.get_fval();
    zprojapprox[1]   = ITC->zL_1_final.get_fval();
    zprojapprox[2]   = ITC->zL_2_final.get_fval();

    phiprojdiskapprox[0] = ITC->phiD_0_final.get_fval();
    phiprojdiskapprox[1] = ITC->phiD_1_final.get_fval();
    phiprojdiskapprox[2] = ITC->phiD_2_final.get_fval();
    phiprojdiskapprox[3] = ITC->phiD_3_final.get_fval();

    rprojdiskapprox[0] = ITC->rD_0_final.get_fval();
    rprojdiskapprox[1] = ITC->rD_1_final.get_fval();
    rprojdiskapprox[2] = ITC->rD_2_final.get_fval();
    rprojdiskapprox[3] = ITC->rD_3_final.get_fval();


    //now binary

    int irinv,iphi0,it,iz0;
    int iphiproj[3],izproj[3];
    
    int iphiprojdisk[4],irprojdisk[4];
    
    int ir2=innerFPGAStub->ir();
    int iphi2=innerFPGAStub->iphi();
    int iz2=innerFPGAStub->iz();
      
    int ir1=outerFPGAStub->ir();
    int iphi1=outerFPGAStub->iphi();
    int iz1=outerFPGAStub->iz();
      
    //To get global precission
    ir1<<=(8-settings_->nrbitsstub(ll-1));
    iphi1<<=(settings_->nphibitsstub(5)-settings_->nphibitsstub(0));
    iphi2<<=(settings_->nphibitsstub(5)-settings_->nphibitsstub(0));

    ITC->r1.set_ival(ir1);
    ITC->r2.set_ival(ir2);
    ITC->z1.set_ival(iz1);
    ITC->z2.set_ival(iz2);
    ITC->phi1.set_ival(iphi1);
    ITC->phi2.set_ival(iphi2);
      
    ITC->rinv_final.calculate();
    ITC->phi0_final.calculate();
    ITC->t_final.calculate();
    ITC->z0_final.calculate();

    ITC->phiL_0_final.calculate();
    ITC->phiL_1_final.calculate();
    ITC->phiL_2_final.calculate();

    ITC->zL_0_final.calculate();
    ITC->zL_1_final.calculate();
    ITC->zL_2_final.calculate();

    ITC->phiD_0_final.calculate();
    ITC->phiD_1_final.calculate();
    ITC->phiD_2_final.calculate();
    ITC->phiD_3_final.calculate();

    ITC->rD_0_final.calculate();
    ITC->rD_1_final.calculate();
    ITC->rD_2_final.calculate();
    ITC->rD_3_final.calculate();

    ITC->der_phiL_final.calculate();
    ITC->der_zL_final.calculate();
    ITC->der_phiD_final.calculate();
    ITC->der_rD_final.calculate();

    //store the binary results
    irinv = ITC->rinv_final.get_ival();
    iphi0 = ITC->phi0_final.get_ival();
    it    = ITC->t_final.get_ival();
    iz0   = ITC->z0_final.get_ival();

    iphiproj[0] = ITC->phiL_0_final.get_ival();
    iphiproj[1] = ITC->phiL_1_final.get_ival();
    iphiproj[2] = ITC->phiL_2_final.get_ival();
    
    izproj[0]   = ITC->zL_0_final.get_ival();
    izproj[1]   = ITC->zL_1_final.get_ival();
    izproj[2]   = ITC->zL_2_final.get_ival();

    iphiprojdisk[0] = ITC->phiD_0_final.get_ival();
    iphiprojdisk[1] = ITC->phiD_1_final.get_ival();
    iphiprojdisk[2] = ITC->phiD_2_final.get_ival();
    iphiprojdisk[3] = ITC->phiD_3_final.get_ival();

    irprojdisk[0]   = ITC->rD_0_final.get_ival();
    irprojdisk[1]   = ITC->rD_1_final.get_ival();
    irprojdisk[2]   = ITC->rD_2_final.get_ival();
    irprojdisk[3]   = ITC->rD_3_final.get_ival();

    if (!goodTrackPars(ITC->rinv_final.local_passes(),
		       ITC->z0_final.local_passes())) return false;

    
    if (!inSector(iphi0,irinv,phi0approx,rinvapprox)) return false;

    LayerProjection layerprojs[4];
    DiskProjection diskprojs[5];


    for(int i=0; i<3; ++i){

      //check that zproj is in range
      if (izproj[i]<-(1<<(settings_->nzbitsstub(0)-1))) continue;
      if (izproj[i]>=(1<<(settings_->nzbitsstub(0)-1))) continue;

      //check that phiproj is in range
      if (iphiproj[i]>=(1<<settings_->nphibitsstub(5))-1) continue;
      if (iphiproj[i]<=0) continue;

      //adjust bits for PS modules (no 2S modules in overlap seeds)
      iphiproj[i]>>=(settings_->nphibitsstub(5)-settings_->nphibitsstub(0));
      
      layerprojs[i].init(settings_,i+1,settings_->rmean(i),
			 iphiproj[i],izproj[i],
			 ITC->der_phiL_final.get_ival(),ITC->der_zL_final.get_ival(),
			 phiproj[i],zproj[i],
			 phider[i],zder[i],
			 phiprojapprox[i],zprojapprox[i],
			 ITC->der_phiL_final.get_fval(),ITC->der_zL_final.get_fval());
      
    }


    for(int i=0; i<4; ++i){

      //check that phi projection in range
      if (iphiprojdisk[i]<=0) continue;
      if (iphiprojdisk[i]>=(1<<settings_->nphibitsstub(0))-1) continue;

      //check that r projection in range
      if(irprojdisk[i]<=0 ||
	 irprojdisk[i] > 120. / ITC->rD_0_final.get_K() ) continue;

      diskprojs[i].init(settings_,i+1,rproj_[i],
			iphiprojdisk[i],irprojdisk[i],
			ITC->der_phiD_final.get_ival(),ITC->der_rD_final.get_ival(),
			phiprojdisk[i],rprojdisk[i],
			phiderdisk[i],rderdisk[i],
			phiprojdiskapprox[i],rprojdiskapprox[i],
			ITC->der_phiD_final.get_fval(),ITC->der_rD_final.get_fval());

      
    }

       
    if (settings_->writeMonitorData("TPars")) {
      globals_->ofstream("trackletparsoverlap.txt")  <<"Trackpars "<<disk_
	  <<"   "<<rinv<<" "<<irinv<<" "<<ITC->rinv_final.get_fval()
	  <<"   "<<phi0<<" "<<iphi0<<" "<<ITC->phi0_final.get_fval()
	  <<"   "<<t<<" "<<it<<" "<<ITC->t_final.get_fval()
	  <<"   "<<z0<<" "<<iz0<<" "<<ITC->z0_final.get_fval()
	  <<endl;
    }
	      
    Tracklet* tracklet=new Tracklet(settings_,innerStub,NULL,outerStub,
				    innerFPGAStub,NULL,outerFPGAStub,
				    rinv,phi0,0.0,z0,t,
				    rinvapprox,phi0approx,0.0,
				    z0approx,tapprox,
				    irinv,iphi0,0,iz0,it,
				    layerprojs,
				    diskprojs,
				    false,true);
    
    if (settings_->debugTracklet()) {
      edm::LogVerbatim("Tracklet") << "Found tracklet in overlap = "<<layer_<<" "<<disk_<<" "<<tracklet<<" "<<iSector_;
    }
    
        
    tracklet->setTrackletIndex(trackletpars_->nTracklets());
    tracklet->setTCIndex(TCIndex_);
    
    if (settings_->writeMonitorData("Seeds")) {
      ofstream fout("seeds.txt", ofstream::app);
      fout << __FILE__ << ":" << __LINE__ << " " << name_ << "_" << iSector_ << " " << tracklet->getISeed() << endl;
      fout.close();
    }
    trackletpars_->addTracklet(tracklet);
    
    int layer=outerFPGAStub->layer().value()+1;
    
    if (layer==2) {
      if (tracklet->validProj(1)) {
	addLayerProj(tracklet,1);
      }
    }
    
    
    for(unsigned int disk=2;disk<6;disk++){
      if (layer==2 && disk==5 ) continue;
      if (tracklet->validProjDisk(disk)) {
	addDiskProj(tracklet,disk);
      }
    }

    return true;
    
  }
  

  
    
protected:

  int iSeed_;
  int TCIndex_;
    
  double phimin_;
  double phimax_;
  double phioffset_;
  
  int layer_;
  int disk_;

  int lproj_[4];
  int dproj_[3];



  double rproj_[4];
  double zproj_[3];
  double zprojoverlap_[4];


  
  TrackletParametersMemory* trackletpars_;

  //First index is layer/disk second is phi region
  vector<vector<TrackletProjectionsMemory*> > trackletprojlayers_;
  vector<vector<TrackletProjectionsMemory*> > trackletprojdisks_;

 public:
    
};


#endif

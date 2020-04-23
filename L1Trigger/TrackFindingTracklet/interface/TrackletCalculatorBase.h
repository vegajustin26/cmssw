//This class implementes the tracklet engine
#ifndef TRACKLETCALCULATORBASE_H
#define TRACKLETCALCULATORBASE_H

#include "ProcessBase.h"
#include "Util.h"
#include "GlobalHistTruth.h"


using namespace std;

class TrackletCalculatorBase:public ProcessBase{

public:

 TrackletCalculatorBase(string name, const Settings* const settings, unsigned int iSector):
  ProcessBase(name,settings,iSector) {
  }
  

  void exacttracklet(double r1, double z1, double phi1,
		     double r2, double z2, double phi2, double,
		     double& rinv, double& phi0,
		     double& t, double& z0,
		     double phiproj[4], double zproj[4], 
		     double phider[4], double zder[4],
		     double phiprojdisk[5], double rprojdisk[5], 
		     double phiderdisk[5], double rderdisk[5]) {

    double deltaphi=Util::phiRange(phi1-phi2);
    
    double dist=sqrt(r2*r2+r1*r1-2*r1*r2*cos(deltaphi));
    
    rinv=2*sin(deltaphi)/dist;

    double phi1tmp=phi1-phimin_;    
     
    phi0=Util::phiRange(phi1tmp+asin(0.5*r1*rinv));
    
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
      exactprojdisk(zmean[i],rinv,phi0,t,z0,
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

    double deltaphi=Util::phiRange(phi1-phi2);

    double dist=sqrt(r2*r2+r1*r1-2*r1*r2*cos(deltaphi));
    
    rinv=2*sin(deltaphi)/dist;

    double phi1tmp=phi1-phimin_;    

    phi0=Util::phiRange(phi1tmp+asin(0.5*r1*rinv));
    
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
      exactproj(rmean[i],rinv,phi0,t,z0,
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

    double deltaphi=Util::phiRange(phi1-phi2);

    double dist=sqrt(r2*r2+r1*r1-2*r1*r2*cos(deltaphi));
    
    rinv=2*sin(deltaphi)/dist;

    if (r1>r2) rinv=-rinv;

    double phi1tmp=phi1-phimin_;    

    
    //cout << "phi1 phi2 phi1tmp : "<<phi1<<" "<<phi2<<" "<<phi1tmp<<endl;

    phi0=Util::phiRange(phi1tmp+asin(0.5*r1*rinv));
    
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
      exactproj(rmean[i],rinv,phi0,t,z0,
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

    if (fpgar.value()*krprojshiftdisk<12.0) return;
    if (fpgar.value()*krprojshiftdisk>112.0) return;


    
    FPGAWord fpgaphi=tracklet->fpgaphiprojdisk(disk);
    
    int iphivmRaw=fpgaphi.value()>>(fpgaphi.nbits()-5);

    int iphi=iphivmRaw/(32/nallstubsdisks[abs(disk)-1]);

    addProjectionDisk(disk,iphi,trackletprojdisks_[abs(disk)-1][iphi],tracklet);

  }


  bool addLayerProj(Tracklet* tracklet, int layer){

    assert(layer>0);

    FPGAWord fpgaz=tracklet->fpgazproj(layer);
    FPGAWord fpgaphi=tracklet->fpgaphiproj(layer);


    if(fpgaphi.atExtreme()) cout<<"at extreme! "<<fpgaphi.value()<<"\n";

    assert(!fpgaphi.atExtreme());
    
    if (fpgaz.atExtreme()) return false;

    if (std::abs(fpgaz.value()*kz)>zlength) return false;

    int iphivmRaw=fpgaphi.value()>>(fpgaphi.nbits()-5);

    int iphi=iphivmRaw/(32/nallstubslayers[layer-1]);

    addProjection(layer,iphi,trackletprojlayers_[layer-1][iphi],tracklet);
    
    return true;

  }

  void addProjection(int layer,int iphi,TrackletProjectionsMemory* trackletprojs, Tracklet* tracklet){
    if (trackletprojs==0) {
      if (warnNoMem) {
	cout << "No projection memory exists in "<<getName()<<" for layer = "<<layer<<" iphi = "<<iphi+1<<endl;
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
      if (warnNoMem) {       
	cout << "No projection memory exists in "<<getName()<<" for disk = "<<abs(disk)<<" iphi = "<<iphi+1<<endl;
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
      if (debug1) {
	cout << getName()<<" TrackletCalculatorBase irinv too large"<<endl;
      }
      success = false;
    }
    if (!goodz0){
      if (debug1) {
	cout << getName()<<" TrackletCalculatorBase z0 cut to large"<<endl;
      }
      success = false;
    }
    return success;
  }


  bool inSector(int iphi0, int irinv, double phi0approx, double rinvapprox) {
      
    double phicritapprox=phi0approx-asin(0.5*rcrit*rinvapprox);

    int ifactor=0.5*rcrit*krinvpars/kphi0pars*(1<<8);
    int iphicrit=iphi0-(irinv>>8)*ifactor;
  
    int iphicritmincut=phicritminmc/ITC_L1L2.phi0_final.get_K();
    int iphicritmaxcut=phicritmaxmc/ITC_L1L2.phi0_final.get_K(); 

    bool keepapprox=(phicritapprox>phicritminmc)&&(phicritapprox<phicritmaxmc),
         keep=(iphicrit>iphicritmincut)&&(iphicrit<iphicritmaxcut);
    if (debug1)
      if (keepapprox && !keep)
        cout << getName() << " Tracklet kept with exact phicrit cut but not approximate, phicritapprox: " << phicritapprox << endl;
    if (usephicritapprox) {
      return keepapprox;
    } else {
      return keep;
    }

    return true;
    
  }
    

  
  
  bool barrelSeeding(Stub* innerFPGAStub, L1TStub* innerStub, Stub* outerFPGAStub, L1TStub* outerStub){
	  
    if (debug1) {
      cout << "TrackletCalculator "<<getName()<<" "<<layer_<<" trying stub pair in layer (inner outer): "
	   <<innerFPGAStub->layer().value()<<" "<<outerFPGAStub->layer().value()<<endl;
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

    if (useapprox) {
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
    if(layer_==1)      ITC = &ITC_L1L2;
    else if(layer_==2) ITC = &ITC_L2L3;
    else if(layer_==3) ITC = &ITC_L3L4;
    else               ITC = &ITC_L5L6;
    
    ITC->r1.set_fval(r1-rmean[layer_-1]);
    ITC->r2.set_fval(r2-rmean[layer_]);
    ITC->z1.set_fval(z1);
    ITC->z2.set_fval(z2);
    double sphi1 = phi1 - phioffset_;
    if(sphi1<0) sphi1 += 8*atan(1.);
    if(sphi1>8*atan(1.)) sphi1 -= 8*atan(1.);
    //cout << "sphi1: "<<phi2<<endl;
    double sphi2 = phi2 - phioffset_;
    if(sphi2<0) sphi2 += 8*atan(1.);
    if(sphi2>8*atan(1.)) sphi2 -= 8*atan(1.);
    ITC->phi1.set_fval(sphi1);
    ITC->phi2.set_fval(sphi2);

    ITC->rproj0.set_fval(rproj_[0]);
    ITC->rproj1.set_fval(rproj_[1]);
    ITC->rproj2.set_fval(rproj_[2]);
    ITC->rproj3.set_fval(rproj_[3]);

    ITC->zproj0.set_fval(t>0? zmean[0] : -zmean[0]);
    ITC->zproj1.set_fval(t>0? zmean[1] : -zmean[1]);
    ITC->zproj2.set_fval(t>0? zmean[2] : -zmean[2]);
    ITC->zproj3.set_fval(t>0? zmean[3] : -zmean[3]);
    ITC->zproj4.set_fval(t>0? zmean[4] : -zmean[4]);

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
    
    /*
    cout << "texact:"<<t<<" tapprox:"<<tapprox<<endl;
    cout << "z0exact:"<<z0<<" z0approx:"<<z0approx<<" "
	 <<" zeroth order:"<<z1-tapprox*r1
	 <<" first order:"
	 <<z1-tapprox*r1*(1+r1*r1*rinvapprox*rinvapprox/24.0)
	 <<" wrong first order:"
	 <<z1-tapprox*r1*(1+r1*r1*rinvapprox*rinvapprox/6.0)
	 <<endl;
    */
    
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
    
    if (layer_<4) iphi1<<=(nbitsphistubL456-nbitsphistubL123);
    if (layer_<3) iphi2<<=(nbitsphistubL456-nbitsphistubL123);
    if (layer_<4) {
      ir1<<=(8-nbitsrL123);
    } else {
      ir1<<=(8-nbitsrL456);
    }
    if (layer_<3) {
      ir2<<=(8-nbitsrL123);
    } else {
      ir2<<=(8-nbitsrL456);
    }
    if (layer_>3) iz1<<=(nbitszL123-nbitszL456);
    if (layer_>2) iz2<<=(nbitszL123-nbitszL456);  
      
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
      if (izproj[i]<-(1<<(nbitszprojL123-1))) continue;
      if (izproj[i]>=(1<<(nbitszprojL123-1))) continue;

      //reject projection if phi is out of range
      if (iphiproj[i]>=(1<<nbitsphistubL456)-1) continue;
      if (iphiproj[i]<=0) continue;
      
      //Adjust bits for r and z projection depending on layer
      if (lproj_[i]<=3) {
	iphiproj[i]>>=(nbitsphistubL456-nbitsphistubL123);
      }
      else {
	izproj[i]>>=(nbitszprojL123-nbitszprojL456);
      }

      layerprojs[i].init(lproj_[i],rproj_[i],
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
	if (iphiprojdisk[i]>=(1<<nbitsphistubL123)-1) continue;
	
	if(irprojdisk[i]< 20. / ITC->rD_0_final.get_K() ||
	   irprojdisk[i] > 120. / ITC->rD_0_final.get_K() ) continue;

	diskprojs[i].init(i+1,rproj_[i],
			   iphiprojdisk[i],irprojdisk[i],
			   ITC->der_phiD_final.get_ival(),ITC->der_rD_final.get_ival(),
			   phiprojdisk[i],rprojdisk[i],
			   phiderdisk[i],rderdisk[i],
			   phiprojdiskapprox[i],rprojdiskapprox[i],
			   ITC->der_phiD_final.get_fval(),ITC->der_rD_final.get_fval());
	
      }
    }
 
    
    if (writeTrackletPars) {
      static ofstream out("trackletpars.txt");
      out <<"Trackpars "<<layer_
	  <<"   "<<rinv<<" "<<rinvapprox<<" "<<ITC->rinv_final.get_fval()
	  <<"   "<<phi0<<" "<<phi0approx<<" "<<ITC->phi0_final.get_fval()
	  <<"   "<<t<<" "<<tapprox<<" "<<ITC->t_final.get_fval()
	  <<"   "<<z0<<" "<<z0approx<<" "<<ITC->z0_final.get_fval()
	  <<endl;
    }	        
        
    Tracklet* tracklet=new Tracklet(innerStub,NULL,outerStub,
				    innerFPGAStub,NULL,outerFPGAStub,
				    rinv,phi0,0.0,z0,t,
				    rinvapprox,phi0approx,0.0,
				    z0approx,tapprox,
				    irinv,iphi0,0,iz0,it,
				    layerprojs,
				    diskprojs,
				    false);
    
    if (debug1) {
      cout << "TrackletCalculator "<<getName()<<" Found tracklet in layer = "<<layer_<<" "
	   <<iSector_<<" phi0 = "<<phi0<<endl;
    }
        

    tracklet->setTrackletIndex(trackletpars_->nTracklets());
    tracklet->setTCIndex(TCIndex_);

    if (writeSeeds) {
      ofstream fout("seeds.txt", ofstream::app);
      fout << __FILE__ << ":" << __LINE__ << " " << name_ << "_" << iSector_ << " " << tracklet->getISeed() << endl;
      fout.close();
    }
    trackletpars_->addTracklet(tracklet);

    HistBase* hists=GlobalHistTruth::histograms();
    int tp=tracklet->tpseed();
    hists->fillTrackletParams(iSeed_,iSector_,
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
      //	    cout<<" LL to L "<<lproj[j]<<"\n";
      bool added=false;
      if (tracklet->validProj(lproj_[j])) {
	added=addLayerProj(tracklet,lproj_[j]);
	//cout << "Add tracklet proj for layer "<<lproj_[j]<<": "<<phiproj[j]<<" "<<iphiproj[j]<<" added = "
	//     <<added<<endl;
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
      //	    cout<<" LL to disk "<<disk<<"\n";
      if (tracklet->validProjDisk(abs(disk))) {
	//cout << "Add tracklet "<<tracklet<<" for disk "<<disk<<endl;
	addDiskProj(tracklet,disk);
      }
    }
    
    return true;

  }
    
  bool diskSeeding(Stub* innerFPGAStub,L1TStub* innerStub,Stub* outerFPGAStub,L1TStub* outerStub){

	    
    if (debug1) {
      cout <<  "TrackletCalculator::execute calculate disk seeds" << endl;
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
    if (useapprox) {
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
    if(disk==1)       ITC = &ITC_F1F2;
    else if(disk==3)  ITC = &ITC_F3F4;
    else if(disk==-1) ITC = &ITC_B1B2;
    else               ITC = &ITC_B3B4;
    
    ITC->r1.set_fval(r1);
    ITC->r2.set_fval(r2);
    int signt = t>0? 1 : -1;
    ITC->z1.set_fval(z1-signt*zmean[abs(disk_)-1]);
    ITC->z2.set_fval(z2-signt*zmean[abs(disk_)]);
    double sphi1 = phi1 - phioffset_;
    if(sphi1<0) sphi1 += 8*atan(1.);
    if(sphi1>8*atan(1.)) sphi1 -= 8*atan(1.);
    double sphi2 = phi2 - phioffset_;
    if(sphi2<0) sphi2 += 8*atan(1.);
    if(sphi2>8*atan(1.)) sphi2 -= 8*atan(1.);
    ITC->phi1.set_fval(sphi1);
    ITC->phi2.set_fval(sphi2);
	    
    ITC->rproj0.set_fval(rmean[0]);
    ITC->rproj1.set_fval(rmean[1]);
    ITC->rproj2.set_fval(rmean[2]);

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

    /*
    cout << "texact:"<<t<<" tapprox:"<<tapprox<<endl;
    cout << "z0exact:"<<z0<<" z0approx:"<<z0approx<<" "
	 <<" zeroth order:"<<z1-tapprox*r1
	 <<" first order:"
	 <<z1-tapprox*r1*(1+r1*r1*rinvapprox*rinvapprox/24.0)
	 <<" wrong first order:"
	 <<z1-tapprox*r1*(1+r1*r1*rinvapprox*rinvapprox/6.0)
	 <<endl;
    */
    
    
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
    iphi1<<=(nbitsphistubL456-nbitsphistubL123);
    iphi2<<=(nbitsphistubL456-nbitsphistubL123);
    
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
      if (izproj[i]<-(1<<(nbitszprojL123-1))) continue;
      if (izproj[i]>=(1<<(nbitszprojL123-1))) continue;

      //Check if outside phi range
      if (iphiproj[i]>=(1<<nbitsphistubL456)-1) continue;
      if (iphiproj[i]<=0) continue;

      //shift bits - allways in PS modules for disk seeding
      iphiproj[i]>>=(nbitsphistubL456-nbitsphistubL123);
      
      layerprojs[i].init(i+1,rmean[i],
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
      if (iphiprojdisk[i]>=(1<<nbitsphistubL123)-1) continue;

      //check that r projection in range
      if(irprojdisk[i]<=0 ||
	 irprojdisk[i] > 120. / ITC->rD_0_final.get_K() ) continue;

      diskprojs[i].init(i+1,rproj_[i],
			iphiprojdisk[i],irprojdisk[i],
			ITC->der_phiD_final.get_ival(),ITC->der_rD_final.get_ival(),
			phiprojdisk[i],rprojdisk[i],
			phiderdisk[i],rderdisk[i],
			phiprojdiskapprox[i],rprojdiskapprox[i],
			ITC->der_phiD_final.get_fval(),ITC->der_rD_final.get_fval());

      
    }

    
    if (writeTrackletParsDisk) {
      static ofstream out("trackletparsdisk.txt");
      out <<"Trackpars         "<<disk_
	  <<"   "<<rinv<<" "<<rinvapprox<<" "<<ITC->rinv_final.get_fval()
	  <<"   "<<phi0<<" "<<phi0approx<<" "<<ITC->phi0_final.get_fval()
	  <<"   "<<t<<" "<<tapprox<<" "<<ITC->t_final.get_fval()
	  <<"   "<<z0<<" "<<z0approx<<" "<<ITC->z0_final.get_fval()
	  <<endl;
    }
	    
    Tracklet* tracklet=new Tracklet(innerStub,NULL,outerStub,
				    innerFPGAStub,NULL,outerFPGAStub,
				    rinv,phi0,0.0,z0,t,
				    rinvapprox,phi0approx,0.0,
				    z0approx,tapprox,
				    irinv,iphi0,0,iz0,it,
				    layerprojs,
				    diskprojs,
				    true);
    
    if (debug1) {
      cout << "Found tracklet in disk = "<<disk_<<" "<<tracklet
	   <<" "<<iSector_<<endl;
    }
        
    tracklet->setTrackletIndex(trackletpars_->nTracklets());
    tracklet->setTCIndex(TCIndex_);

    if (writeSeeds) {
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

    if (debug1) {
      cout << "trying to make overlap tracklet disk_ = "<<disk_<<" "<<getName()<<endl;
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
      //cout << "in overlap tracklet: radii wrong"<<endl;
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
    if (useapprox) {
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
    if     (ll==1 && disk==1)  ITC = &ITC_L1F1;
    else if(ll==2 && disk==1)  ITC = &ITC_L2F1;
    else if(ll==1 && disk==-1) ITC = &ITC_L1B1;
    else if(ll==2 && disk==-1) ITC = &ITC_L2B1;
    else assert(0);
    
    ITC->r1.set_fval(r2-rmean[ll-1]);
    ITC->r2.set_fval(r1);
    int signt = t>0? 1 : -1;
    ITC->z1.set_fval(z2);
    ITC->z2.set_fval(z1-signt*zmean[abs(disk_)-1]);
    double sphi1 = phi1 - phioffset_;
    if(sphi1<0) sphi1 += 8*atan(1.);
    if(sphi1>8*atan(1.)) sphi1 -= 8*atan(1.);
    double sphi2 = phi2 - phioffset_;
    if(sphi2<0) sphi2 += 8*atan(1.);
    if(sphi2>8*atan(1.)) sphi2 -= 8*atan(1.);
    ITC->phi1.set_fval(sphi2);
    ITC->phi2.set_fval(sphi1);

    ITC->rproj0.set_fval(rmean[0]);
    ITC->rproj1.set_fval(rmean[1]);
    ITC->rproj2.set_fval(rmean[2]);
    
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
    ir1<<=(8-nbitsrL123);
    iphi1<<=(nbitsphistubL456-nbitsphistubL123);
    iphi2<<=(nbitsphistubL456-nbitsphistubL123);

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
      if (izproj[i]<-(1<<(nbitszprojL123-1))) continue;
      if (izproj[i]>=(1<<(nbitszprojL123-1))) continue;

      //check that phiproj is in range
      if (iphiproj[i]>=(1<<nbitsphistubL456)-1) continue;
      if (iphiproj[i]<=0) continue;

      //adjust bits for PS modules (no 2S modules in overlap seeds)
      iphiproj[i]>>=(nbitsphistubL456-nbitsphistubL123);
      
      layerprojs[i].init(i+1,rmean[i],
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
      if (iphiprojdisk[i]>=(1<<nbitsphistubL123)-1) continue;

      //check that r projection in range
      if(irprojdisk[i]<=0 ||
	 irprojdisk[i] > 120. / ITC->rD_0_final.get_K() ) continue;

      diskprojs[i].init(i+1,rproj_[i],
			iphiprojdisk[i],irprojdisk[i],
			ITC->der_phiD_final.get_ival(),ITC->der_rD_final.get_ival(),
			phiprojdisk[i],rprojdisk[i],
			phiderdisk[i],rderdisk[i],
			phiprojdiskapprox[i],rprojdiskapprox[i],
			ITC->der_phiD_final.get_fval(),ITC->der_rD_final.get_fval());

      
    }

       
    if (writeTrackletParsOverlap) {
      static ofstream out("trackletparsoverlap.txt");
      out <<"Trackpars "<<disk_
	  <<"   "<<rinv<<" "<<irinv<<" "<<ITC->rinv_final.get_fval()
	  <<"   "<<phi0<<" "<<iphi0<<" "<<ITC->phi0_final.get_fval()
	  <<"   "<<t<<" "<<it<<" "<<ITC->t_final.get_fval()
	  <<"   "<<z0<<" "<<iz0<<" "<<ITC->z0_final.get_fval()
	  <<endl;
    }
	      
    Tracklet* tracklet=new Tracklet(innerStub,NULL,outerStub,
				    innerFPGAStub,NULL,outerFPGAStub,
				    rinv,phi0,0.0,z0,t,
				    rinvapprox,phi0approx,0.0,
				    z0approx,tapprox,
				    irinv,iphi0,0,iz0,it,
				    layerprojs,
				    diskprojs,
				    false,true);
    
    if (debug1) {
      cout << "Found tracklet in overlap = "<<layer_<<" "<<disk_
	   <<" "<<tracklet<<" "<<iSector_<<endl;
    }
    
        
    tracklet->setTrackletIndex(trackletpars_->nTracklets());
    tracklet->setTCIndex(TCIndex_);
    
    if (writeSeeds) {
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
  static IMATH_TrackletCalculator ITC_L1L2;
  static IMATH_TrackletCalculator ITC_L2L3;
  static IMATH_TrackletCalculator ITC_L3L4;
  static IMATH_TrackletCalculator ITC_L5L6;
  
  static IMATH_TrackletCalculatorDisk ITC_F1F2;
  static IMATH_TrackletCalculatorDisk ITC_F3F4;
  static IMATH_TrackletCalculatorDisk ITC_B1B2;
  static IMATH_TrackletCalculatorDisk ITC_B3B4;

  static IMATH_TrackletCalculatorOverlap ITC_L1F1;
  static IMATH_TrackletCalculatorOverlap ITC_L2F1;
  static IMATH_TrackletCalculatorOverlap ITC_L1B1;
  static IMATH_TrackletCalculatorOverlap ITC_L2B1;

    
};

IMATH_TrackletCalculator TrackletCalculatorBase::ITC_L1L2{1,2};
IMATH_TrackletCalculator TrackletCalculatorBase::ITC_L2L3{2,3};
IMATH_TrackletCalculator TrackletCalculatorBase::ITC_L3L4{3,4};
IMATH_TrackletCalculator TrackletCalculatorBase::ITC_L5L6{5,6};

IMATH_TrackletCalculatorDisk TrackletCalculatorBase::ITC_F1F2{1,2};
IMATH_TrackletCalculatorDisk TrackletCalculatorBase::ITC_F3F4{3,4};
IMATH_TrackletCalculatorDisk TrackletCalculatorBase::ITC_B1B2{-1,-2};
IMATH_TrackletCalculatorDisk TrackletCalculatorBase::ITC_B3B4{-3,-4};

IMATH_TrackletCalculatorOverlap TrackletCalculatorBase::ITC_L1F1{1,1};
IMATH_TrackletCalculatorOverlap TrackletCalculatorBase::ITC_L2F1{2,1};
IMATH_TrackletCalculatorOverlap TrackletCalculatorBase::ITC_L1B1{1,-1};
IMATH_TrackletCalculatorOverlap TrackletCalculatorBase::ITC_L2B1{2,-1};


#endif

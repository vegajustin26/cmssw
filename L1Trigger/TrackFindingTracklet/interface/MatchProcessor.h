//This class implementes the match processor
#ifndef L1Trigger_TrackFindingTracklet_interface_MatchProcessor_h
#define L1Trigger_TrackFindingTracklet_interface_MatchProcessor_h

#include "ProcessBase.h"
#include "CircularBuffer.h"
#include "MatchEngineUnit.h"
#include "ProjectionTemp.h"
#include "Util.h"

using namespace std;

class MatchProcessor:public ProcessBase{

public:

 MatchProcessor(string name, const Settings* settings, unsigned int iSector):
  ProcessBase(name,settings,iSector), fullmatches_(12), inputProjBuffer_(3){
    
    double dphi=2*M_PI/NSector;
    double dphiHG=0.5*dphisectorHG-M_PI/NSector;
    phimin_=iSector_*dphi-dphiHG;
    phimax_=phimin_+dphi+2*dphiHG;

    phimin_-=M_PI/NSector;
    phimax_-=M_PI/NSector;
    if (phimin_>M_PI) {
      phimin_-=2*M_PI;
      phimax_-=2*M_PI;
    }

    phioffset_=phimin_;

    initLayerDisk(3,layer_,disk_);

    //FIXME should sort out constants here
    icorrshift_=7+idrinvbits+phi0bitshift-rinvbitshift-phiderbitshift;
    if (layer_<=3) {
      icorzshift_=-1-PS_zderL_shift;
    } else {
      icorzshift_=-1-SS_zderL_shift;      
    }
    phi0shift_=3;
    fact_=1;
    if (layer_>=4) {
      fact_=(1<<(settings_->nzbitsstub(0)-settings_->nzbitsstub(5)));
      icorrshift_-=(10-settings_->nrbitsstub(layer_-1));
      icorzshift_+=(settings_->nzbitsstub(0)-settings_->nzbitsstub(5)+settings_->nrbitsstub(layer_-1)-settings_->nrbitsstub(0));
      phi0shift_=0;
    }


    nrbits_=5;
    nphiderbits_=6;

    
    //to adjust globaly the phi and rz matching cuts
    phifact_=1.0;
    rzfact_=1.0;


    for(unsigned int iSeed=0;iSeed<12;iSeed++) {
      if (layer_>0) {
	phimatchcut_[iSeed]=settings_->rphimatchcut(iSeed,layer_-1)/(kphi1*settings_->rmean(layer_-1));
	zmatchcut_[iSeed]=settings_->zmatchcut(iSeed,layer_-1)/kz;
      }
      if (disk_!=0) {
	rphicutPS_[iSeed]=settings_->rphicutPS(iSeed,abs(disk_)-1)/(kphiproj123*kr);
	rphicut2S_[iSeed]=settings_->rphicut2S(iSeed,abs(disk_)-1)/(kphiproj123*kr);
	rcut2S_[iSeed]=settings_->rcut2S(iSeed,abs(disk_)-1)/krprojshiftdisk;
	rcutPS_[iSeed]=settings_->rcutPS(iSeed,abs(disk_)-1)/krprojshiftdisk;
      }
    }

    
    
    if (iSector_==0&&layer_>0&&settings_->writeTable()) {

      ofstream outphicut;
      outphicut.open(getName()+"_phicut.tab");
      outphicut << "{"<<endl;
      for(unsigned int seedindex=0;seedindex<12;seedindex++){
	if (seedindex!=0) outphicut<<","<<endl;
	outphicut << phimatchcut_[seedindex];
      }
      outphicut <<endl<<"};"<<endl;
      outphicut.close();

      ofstream outzcut;
      outzcut.open(getName()+"_zcut.tab");
      outzcut << "{"<<endl;
      for(unsigned int seedindex=0;seedindex<12;seedindex++){
	if (seedindex!=0) outzcut<<","<<endl;
	outzcut << zmatchcut_[seedindex];
      }
      outzcut <<endl<<"};"<<endl;
      outzcut.close();
    }



    

    if (layer_>0) {

      unsigned int nbits=3;
      if (layer_>=4) nbits=4;
      
      for(unsigned int irinv=0;irinv<32;irinv++){
	double rinv=(irinv-15.5)*(1<<(nbitsrinv-5))*krinvpars;
	double projbend=bend(settings_->rmean(layer_-1),rinv);
	for(unsigned int ibend=0;ibend<(unsigned int)(1<<nbits);ibend++){
	  double stubbend=Stub::benddecode(ibend,layer_<=3);
	  bool pass=std::abs(stubbend-projbend)<mecut;
	  table_.push_back(pass);
	}
      }

      if (settings_->writeTable()){
	ofstream out;
	char layer='0'+layer_;
	string fname="METable_L";
	fname+=layer;
	fname+=".tab";
	out.open(fname.c_str());
	out << "{" <<endl;
	for(unsigned int i=0;i<table_.size();i++){
	  if (i!=0) {
	    out <<","<<endl;
	  }
	  out << table_[i] ;
	}
	out << "};"<<endl;
	out.close();
      }
      
    }

    if (disk_>0) {

      for(unsigned int iprojbend=0;iprojbend<32;iprojbend++){
	double projbend=0.5*(iprojbend-15.0);
	for(unsigned int ibend=0;ibend<8;ibend++){
	  double stubbend=Stub::benddecode(ibend,true);
	  bool pass=std::abs(stubbend-projbend)<mecutdisk;
	  tablePS_.push_back(pass);
	}
	for(unsigned int ibend=0;ibend<16;ibend++){
	  double stubbend=Stub::benddecode(ibend,false);
	  bool pass=std::abs(stubbend-projbend)<mecutdisk;
	  table2S_.push_back(pass);
	}
      }
      
    }


    for (unsigned int i=0; i<10; i++) {
      ialphafactinner_[i]= 
	(1<<settings_->alphashift())*krprojshiftdisk*half2SmoduleWidth/(1<<(settings_->nbitsalpha()-1))/(settings_->rDSSinner(i)*settings_->rDSSinner(i))/kphiproj123;
      ialphafactouter_[i]= 
	(1<<settings_->alphashift())*krprojshiftdisk*half2SmoduleWidth/(1<<(settings_->nbitsalpha()-1))/(settings_->rDSSouter(i)*settings_->rDSSouter(i))/kphiproj123;
    }

    barrel_=layer_>0;

    nvm_=barrel_?settings_->nvmme(layer_-1)*settings_->nallstubs(layer_-1):settings_->nvmme(disk_+5)*settings_->nallstubs(disk_+5);
    nvmbins_=barrel_?settings_->nvmme(layer_-1):settings_->nvmme(disk_+5);

    if (nvm_==32) nvmbits_=5;
    if (nvm_==16) nvmbits_=4;
    assert(nvmbits_!=-1);


    nMatchEngines_=4;
    for(unsigned int iME=0;iME<nMatchEngines_;iME++){
      MatchEngineUnit tmpME(barrel_,table_,tablePS_,table2S_);
      matchengines_.push_back(tmpME);
    }
 
    
  }

  void addOutput(MemoryBase* memory,string output){
    if (settings_->writetrace()) {
      cout << "In "<<name_<<" adding output to "<<memory->getName() << " to output "<<output<<endl;
    }
    if (output.find("matchout")!=std::string::npos){
      FullMatchMemory* tmp=dynamic_cast<FullMatchMemory*>(memory);
      assert(tmp!=0);
      unsigned int iSeed=getISeed(tmp->getName());
      assert(iSeed<fullmatches_.size());
      assert(fullmatches_[iSeed]==0);
      fullmatches_[iSeed]=tmp;
      return;
    }
    cout << "Could not find output = "<<output<<endl;
    assert(0);
  }

  void addInput(MemoryBase* memory,string input){
    if (settings_->writetrace()) {
      cout << "In "<<name_<<" adding input from "<<memory->getName() << " to input "<<input<<endl;
    }
    if (input=="allstubin"){
      AllStubsMemory* tmp=dynamic_cast<AllStubsMemory*>(memory);
      assert(tmp!=0);
      allstubs_=tmp;
      return;
    }
    if (input=="vmstubin"){
      VMStubsMEMemory* tmp=dynamic_cast<VMStubsMEMemory*>(memory);
      assert(tmp!=0);
      vmstubs_.push_back(tmp); //to allow more than one stub in?  vmstubs_=tmp;
      return;
    }
    if (input=="projin"){
      TrackletProjectionsMemory* tmp=dynamic_cast<TrackletProjectionsMemory*>(memory);
      assert(tmp!=0);
      inputprojs_.push_back(tmp);
      return;
    }
    cout << "MatchProcessor input = "<<input<<endl;
    assert(0);
  }

  void execute() {

    /*
      The code is organized in three 'steps' corresponding to the PR, ME, and MC functions. The output from
      the PR step is buffered in a 'circular' buffer, and similarly the ME output is put in a circular buffer. 

      The implementation is done in steps, emulating what can be done in firmware. One each step we do:

       1) A projection is read and if there is space it is insert into the inputProjBuffer_

       2) Process next match in the ME - if there is an idle ME the next projection is inserted

       3) Readout match from ME and send to match calculator

     */
    

    Tracklet* oldTracklet=0;
    
    unsigned int countall=0;
    unsigned int countsel=0;
    
    
    unsigned int iprojmem=0;
    unsigned int iproj=0;

    inputProjBuffer_.reset();
    
    for (unsigned int istep=0;istep<settings_->maxStep("MP");istep++) {

      //Step 1
      //First step here checks if we have more input projections to put into
      //the input puffer for projections
      if (iprojmem<inputprojs_.size()) {
	TrackletProjectionsMemory* projMem=inputprojs_[iprojmem];
	if (projMem->nTracklets()==0) {
	  iprojmem++;
	} else if (iproj<projMem->nTracklets()) {
	  if (!inputProjBuffer_.almostfull()) {

	    Tracklet* proj=projMem->getFPGATracklet(iproj);

	    FPGAWord fpgaphi=barrel_?proj->fpgaphiproj(layer_):proj->fpgaphiprojdisk(disk_);
	    FPGAWord fpgarz=barrel_?proj->fpgazproj(layer_):proj->fpgarprojdisk(disk_);

	    int iphi=(fpgaphi.value()>>(fpgaphi.nbits()-nvmbits_))&(nvmbins_-1);

	    int projrinv=-1;
	    if (barrel_) {
	      projrinv=16+(proj->fpgarinv().value()>>(proj->fpgarinv().nbits()-5));
	    }else{
	      //The next lines looks up the predicted bend based on:
	      // 1 - r projections
	      // 2 - phi derivative
	      // 3 - the sign - i.e. if track is forward or backward
	      int rindex=(proj->fpgarprojdisk(disk_).value()>>(proj->fpgarprojdisk(disk_).nbits()-nrbits_))&((1<<nrbits_)-1);
	      
	      int phiderindex=(proj->fpgaphiprojderdisk(disk_).value()>>(proj->fpgaphiprojderdisk(disk_).nbits()-nphiderbits_))&((1<<nphiderbits_)-1);
	      
	      int signindex=(proj->fpgarprojderdisk(disk_).value()<0);
	      
	      int bendindex=(signindex<<(nphiderbits_+nrbits_))+
		(rindex<<(nphiderbits_))+
		phiderindex;
	      
	      projrinv=bendTable(abs(disk_)-1,bendindex);
	      
	      proj->setBendIndex(projrinv,disk_);
	  
	    }
	    assert(projrinv>=0);

	    unsigned int slot = barrel_?proj->zbin1projvm(layer_):proj->rbin1projvm(disk_);
	    bool second=(barrel_?proj->zbin2projvm(layer_):proj->rbin2projvm(disk_))==1;
	    	    
	    unsigned int projfinephi=fpgaphi.value()>>(fpgaphi.nbits()-(nvmbits_+3))&7;
	    int projfinerz = barrel_?proj->finezvm(layer_):proj->finervm(disk_);
	  
	    bool isPSseed=proj->PSseed()==1;

	    VMStubsMEMemory* stubmem=vmstubs_[iphi];
	    if (stubmem->nStubsBin(slot)!=0) { 
	      ProjectionTemp tmpProj(proj,slot,projrinv,projfinerz,projfinephi,iphi,isPSseed);
	      inputProjBuffer_.store(tmpProj);
	    }
	    if (second&&(stubmem->nStubsBin(slot+1)!=0)) {
	      ProjectionTemp tmpProj(proj,slot+1,projrinv,projfinerz-8,projfinephi,iphi,isPSseed);
	      inputProjBuffer_.store(tmpProj);
	    }
	    iproj++;
	    if (iproj==projMem->nTracklets()){
	      iproj=0;
	      iprojmem++;
	    }
	  }
	}
      }

      //Step 2
      //Check if we have ME that can process projection

      bool addedProjection=false;
      for (unsigned int iME=0;iME<nMatchEngines_;iME++) {
	matchengines_[iME].step();
	//if match engine empty and we have queued projections add to match engine
	if ((!addedProjection)&&matchengines_[iME].idle()&&(!inputProjBuffer_.empty())) {

	  ProjectionTemp tmpProj=inputProjBuffer_.read();
	  VMStubsMEMemory* stubmem=vmstubs_[tmpProj.iphi()];
	  
	  matchengines_[iME].init(stubmem,
				  tmpProj.slot(),
				  tmpProj.projrinv(),
				  tmpProj.projfinerz(),
				  tmpProj.projfinephi(),
				  tmpProj.isPSseed(),
				  tmpProj.proj());
	  addedProjection=true;
	}
      }
    

      //Step 3
      //Check if we have candidate match to process

      unsigned int iMEbest=nMatchEngines_;
      int bestTCID=-1;
      bool bestInPipeline=false;
      for (unsigned int iME=0;iME<nMatchEngines_;iME++) {
	bool empty=matchengines_[iME].empty();
	if (empty&&matchengines_[iME].idle()) continue;
	int currentTCID=empty?matchengines_[iME].currentProj()->TCID():matchengines_[iME].peek().first->TCID();
	if ((iMEbest==nMatchEngines_)||(currentTCID<bestTCID)) {
	  iMEbest=iME;
	  bestTCID=currentTCID;
	  bestInPipeline=empty;
	}
      }

      if (iMEbest!=nMatchEngines_&&(!bestInPipeline)) {

	std::pair<Tracklet*,std::pair<Stub*,L1TStub*> > candmatch=matchengines_[iMEbest].read();
      
	L1TStub* stub=candmatch.second.second;
	Stub* fpgastub=candmatch.second.first;
	Tracklet* tracklet=candmatch.first;
	
	if (oldTracklet!=0) {
	//allow equal here since we can have more than one cadidate match per tracklet projection
	  assert(oldTracklet->TCID()<=tracklet->TCID());
	}
	oldTracklet=tracklet;
	
	bool match=matchCalculator(tracklet,fpgastub,stub);

	countall++;
	if (match) countsel++;;
	
      }

    }
    
    if (settings_->writeMonitorData("MC")) {
      static ofstream out("matchcalculator.txt");
      out << getName()<<" "<<countall<<" "<<countsel<<endl;
    }

    
  }


  

  bool matchCalculator(Tracklet* tracklet,Stub* fpgastub, L1TStub* stub){
      
    if (layer_!=0) {
	  
	
      int ir=fpgastub->r().value();
      int iphi=tracklet->fpgaphiproj(layer_).value();
      int icorr=(ir*tracklet->fpgaphiprojder(layer_).value())>>icorrshift_;	
      iphi+=icorr;
      
      int iz=tracklet->fpgazproj(layer_).value();
      int izcor=(ir*tracklet->fpgazprojder(layer_).value()+(1<<(icorzshift_-1)))>>icorzshift_;
      iz+=izcor;	
      
      int ideltaz=fpgastub->z().value()-iz;
      int ideltaphi=(fpgastub->phi().value()<<phi0shift_)-(iphi<<(phi0bitshift-1+phi0shift_)); 
      
      
      //Floating point calculations
      
      double phi=stub->phi();
      double r=stub->r();
      double z=stub->z();
      
      
      if (settings_->useapprox()) {
	double dphi=Util::phiRange(phi-fpgastub->phiapprox(phimin_,phimax_));
	assert(std::abs(dphi)<0.001);
	phi=fpgastub->phiapprox(phimin_,phimax_);
	z=fpgastub->zapprox();
	r=fpgastub->rapprox();
      }
      
      if (phi<0) phi+=2*M_PI;
      phi-=phioffset_;
      
      double dr=r-tracklet->rproj(layer_);
      assert(std::abs(dr)<drmax);
      
      double dphi=Util::phiRange(phi-(tracklet->phiproj(layer_)+dr*tracklet->phiprojder(layer_)));
      
      double dz=z-(tracklet->zproj(layer_)+dr*tracklet->zprojder(layer_));
	
      double dphiapprox=Util::phiRange(phi-(tracklet->phiprojapprox(layer_)+
					    dr*tracklet->phiprojderapprox(layer_)));
      
      double dzapprox=z-(tracklet->zprojapprox(layer_)+
			 dr*tracklet->zprojderapprox(layer_));
      
      int seedindex=tracklet->getISeed();
      
      assert(phimatchcut_[seedindex]>0);
      assert(zmatchcut_[seedindex]>0);
      
      bool truthmatch=tracklet->stubtruthmatch(stub);
      
      HistBase* hists=GlobalHistTruth::histograms();
      hists->FillLayerResidual(layer_, seedindex,
			       dphiapprox*settings_->rmean(layer_-1),
			       ideltaphi*kphi1*settings_->rmean(layer_-1),
			       ideltaz*fact_*kz, dz,truthmatch);
      
      
      
      
      if (settings_->writeMonitorData("Residuals")) {
	static ofstream out("layerresiduals.txt");
	
	double pt=0.003*3.8/std::abs(tracklet->rinv());
	  
	out << layer_<<" "<<seedindex<<" "<<pt<<" "<<ideltaphi*kphi1*settings_->rmean(layer_-1)
	    <<" "<<dphiapprox*settings_->rmean(layer_-1)
	    <<" "<<phimatchcut_[seedindex]*kphi1*settings_->rmean(layer_-1)
	    <<"   "<<ideltaz*fact_*kz<<" "<<dz<<" "<<zmatchcut_[seedindex]*kz<<endl;	  
      }
      
      
      bool imatch=(std::abs(ideltaphi)<=phifact_*phimatchcut_[seedindex])&&(std::abs(ideltaz*fact_)<=rzfact_*zmatchcut_[seedindex]);
      
      if (settings_->debugTracklet()) {
	cout << getName()<<" imatch = "<<imatch<<" ideltaphi cut "<<ideltaphi
	     <<" "<<phimatchcut_[seedindex]
	     <<" ideltaz*fact cut "<<ideltaz*fact_<<" "<<zmatchcut_[seedindex]<<endl;
      }

      if (std::abs(dphi)>0.2 || std::abs(dphiapprox)>0.2 ) {
	cout << "WARNING dphi and/or dphiapprox too large : "
	<<dphi<<" "<<dphiapprox<<endl;
      }
      
      assert(std::abs(dphi)<0.2);
      assert(std::abs(dphiapprox)<0.2);
      
      if (imatch) {
	
	std::pair<Stub*,L1TStub*> tmp(fpgastub,stub);
	
	tracklet->addMatch(layer_,ideltaphi,ideltaz,
			   dphi,dz,dphiapprox,dzapprox,
			   (fpgastub->phiregion().value()<<7)+fpgastub->stubindex().value(),
			   stub->r(),tmp);
	

	if (settings_->debugTracklet()) {
	  cout << "Accepted full match in layer " <<getName()
	       << " "<<tracklet
	       << " "<<iSector_<<endl;	   
	}

	int iSeed = tracklet->getISeed();
	assert(fullmatches_[iSeed]!=0);
	fullmatches_[iSeed]->addMatch(tracklet,tmp);

	return true;
      }
      else {
	return false;
      }
    } else {  //disk matches
	
      
      //check that stubs and projections in same half of detector
      assert(stub->z()*tracklet->t()>0.0);
      
      int sign=(tracklet->t()>0.0)?1:-1;
      int disk=sign*disk_;
      assert(disk!=0);
      
      //Perform integer calculations here
      
      int iz=fpgastub->z().value();
      int iphi=tracklet->fpgaphiprojdisk(disk).value();
      
      int shifttmp=t2bits+tbitshift+phi0bitshift+2-rinvbitshiftdisk-phiderdiskbitshift-PS_phiderD_shift;
      assert(shifttmp>=0);
      int iphicorr=(iz*tracklet->fpgaphiprojderdisk(disk).value())>>shifttmp;
      
      iphi+=iphicorr;
      
      int ir=tracklet->fpgarprojdisk(disk).value();
      
      
      int shifttmp2=rprojdiskbitshift+t3shift-rderdiskbitshift;
      
      assert(shifttmp2>=0);
      int ircorr=(iz*tracklet->fpgarprojderdisk(disk).value())>>shifttmp2;
      
      ir+=ircorr;
      
      int ideltaphi=fpgastub->phi().value()*kphi/kphiproj123-iphi; 
      
      
      int irstub = fpgastub->r().value();
      int ialphafact=0;
      if(!stub->isPSmodule()){
	assert(irstub<10);
	if (disk_<=2) {
	  ialphafact = ialphafactinner_[irstub];
	  irstub = settings_->rDSSinner(irstub)/kr;
	} else {
	  ialphafact = ialphafactouter_[irstub];
	  irstub = settings_->rDSSouter(irstub)/kr;
	}
      }
      
      int ideltar=(irstub*krdisk)/krprojshiftdisk-ir;
      
      if (!stub->isPSmodule()) {	  
	int ialphanew=fpgastub->alphanew().value();
	int iphialphacor=((ideltar*ialphanew*ialphafact)>>settings_->alphashift());
	ideltaphi+=iphialphacor;
      }
      
      
      
      
      //Perform floating point calculations here
      
      double phi=stub->phi();
      double z=stub->z();
      double r=stub->r();
      
      
      if (settings_->useapprox()) {
	double dphi=Util::phiRange(phi-fpgastub->phiapprox(phimin_,phimax_));
	assert(std::abs(dphi)<0.001);
	phi=fpgastub->phiapprox(phimin_,phimax_);
	z=fpgastub->zapprox();
	r=fpgastub->rapprox();
      }
      
      if (phi<0) phi+=2*M_PI;
      phi-=phioffset_;
      
      double dz=z-sign*zmean[disk_-1];
      
      if(std::abs(dz) > dzmax){
	cout << __FILE__ << ":" << __LINE__ << " " << name_ << "_" << iSector_ << " " << tracklet->getISeed() << endl;
	cout << "stub "<<stub->z() <<" disk "<<disk<<" "<<dz<<endl;
	assert(std::abs(dz)<dzmax);
      }	
      
      
      double phiproj=tracklet->phiprojdisk(disk)+dz*tracklet->phiprojderdisk(disk);
      
      double rproj=tracklet->rprojdisk(disk)+dz*tracklet->rprojderdisk(disk);
      
      double deltar=r-rproj;
      
	
      double dr=stub->r()-rproj;
      
      double dphi=Util::phiRange(phi-phiproj);
      
      double dphiapprox=Util::phiRange(phi-(tracklet->phiprojapproxdisk(disk)+
					    dz*tracklet->phiprojderapproxdisk(disk)));
      
      double drapprox=stub->r()-(tracklet->rprojapproxdisk(disk)+
				 dz*tracklet->rprojderapproxdisk(disk));
      
      double drphi=dphi*stub->r();
      double drphiapprox=dphiapprox*stub->r();
      
      
      
      if (!stub->isPSmodule()) {
	double alphanew=stub->alphanew();
	dphi+=dr*alphanew*half2SmoduleWidth/stub->r2();;
	dphiapprox+=drapprox*alphanew*half2SmoduleWidth/stub->r2();
	
	drphi+=dr*alphanew*half2SmoduleWidth/stub->r();
	drphiapprox+=dr*alphanew*half2SmoduleWidth/stub->r();
      }
      

      int seedindex=tracklet->getISeed();
      
      int idrphicut=rphicutPS_[seedindex];
      int idrcut=rcutPS_[seedindex]; 
      if (!stub->isPSmodule()) {
	idrphicut=rphicut2S_[seedindex];
	idrcut=rcut2S_[seedindex]; 
      }
      
      double drphicut=idrphicut*kphiproj123*kr;
      double drcut=idrcut*krprojshiftdisk;
      
      
      if (settings_->writeMonitorData("Residuals")) {
	static ofstream out("diskresiduals.txt");
	
	double pt=0.003*3.8/std::abs(tracklet->rinv());
	  
	out << disk_<<" "<<stub->isPSmodule()<<" "<<tracklet->layer()<<" "
	    <<abs(tracklet->disk())<<" "<<pt<<" "
	    <<ideltaphi*kphiproj123*stub->r()<<" "<<drphiapprox<<" "
	    <<drphicut<<" "
	    <<ideltar*krprojshiftdisk<<" "<<deltar<<" "
	    <<drcut<<" "
	    <<endl;	  
      }
      
      
      bool match=(std::abs(drphi)<drphicut)&&(std::abs(deltar)<drcut);
	
      bool imatch=(std::abs(ideltaphi*irstub)<idrphicut)&&(std::abs(ideltar)<idrcut);
      
      
      if (settings_->debugTracklet()) {
	cout << "imatch match disk: "<<imatch<<" "<<match<<" "
	     <<std::abs(ideltaphi)<<" "<<drphicut/(kphiproj123*stub->r())<<" "
	     <<std::abs(ideltar)<<" "<<drcut/krprojshiftdisk<<" r = "<<stub->r()<<endl;
      }
      
      
      if (imatch) {
	
	std::pair<Stub*,L1TStub*> tmp(fpgastub,stub);
	
	if (settings_->debugTracklet()) {
	  cout << "MatchCalculator found match in disk "<<getName()<<endl;
	}
	
	if(std::abs(dphi)>=0.25){
	  cout<<"dphi "<<dphi<<"\n";
	  cout<<"ISeed "<<tracklet->getISeed()<<"\n";
          }
	assert(std::abs(dphi)<0.25);
	assert(std::abs(dphiapprox)<0.25);
	
	tracklet->addMatchDisk(disk,ideltaphi,ideltar,
			       drphi/stub->r(),dr,drphiapprox/stub->r(),drapprox,
			       stub->alpha(),
			       (fpgastub->phiregion().value()<<7)+fpgastub->stubindex().value(),
			       stub->z(),tmp);
	if (settings_->debugTracklet()) {
	  cout << "Accepted full match in disk " <<getName()
	       << " "<<tracklet
	       << " "<<iSector_<<endl;	   
	}

	int iSeed = tracklet->getISeed();
	assert(fullmatches_[iSeed]!=0);
	fullmatches_[iSeed]->addMatch(tracklet,tmp);

	return true;
      } else {
	return false;
      }
    }
  }
  

  
  double bend(double r, double rinv) {

    double dr=0.18;
    
    double delta=r*dr*0.5*rinv;

    double bend=-delta/0.009;
    if (r<55.0) bend=-delta/0.01;
    
    return bend;
    
  }

  int bendTable(int diskindex,int bendindex) {

    static vector<int> bendtable[5];

    static bool first=true;

    if (first) {
      first=false;
    
      for (unsigned int idisk=0;idisk<5;idisk++) {

	unsigned int nsignbins=2;
	unsigned int nrbins=1<<(nrbits_);
	unsigned int nphiderbins=1<<(nphiderbits_);
      
	for(unsigned int isignbin=0;isignbin<nsignbins;isignbin++) {
	  for(unsigned int irbin=0;irbin<nrbins;irbin++) {
	    int ir=irbin;
	    if (ir>(1<<(nrbits_-1))) ir-=(1<<nrbits_);
	    ir=ir<<(settings_->nrbitsstub(6)-nrbits_);
	    for(unsigned int iphiderbin=0;iphiderbin<nphiderbins;iphiderbin++) {
	      int iphider=iphiderbin;
	      if (iphider>(1<<(nphiderbits_-1))) iphider-=(1<<nphiderbits_);
	      iphider=iphider<<(settings_->nbitsphiprojderL123()-nphiderbits_);
	      
	      double rproj=ir*krprojshiftdisk;
	      double phider=iphider*GlobalHistTruth::ITC_L1L2()->der_phiD_final.get_K();
	      double t=zmean[idisk]/rproj;
	      
	      if (isignbin) t=-t;
	  
	      double rinv=-phider*(2.0*t);

	      double bendproj=0.5*bend(rproj,rinv);

	    
	      int ibendproj=2.0*bendproj+15.5;
	      if (ibendproj<0) ibendproj=0;
	      if (ibendproj>31) ibendproj=31;
	      
	      bendtable[idisk].push_back(ibendproj);

	    }
	  }
	}
      }
    }

    

    return bendtable[diskindex][bendindex];

  }

  
    
private:

  int layer_;
  int disk_;
  bool barrel_;
  int nvm_;  //VMs in sector
  int nvmbits_; //# of bits for VMs in sector
  int nvmbins_; //VMs in in phi region

  int fact_;
  int icorrshift_;
  int icorzshift_;
  int phi0shift_;

  double phimin_;
  double phimax_;
  double phioffset_;

  unsigned int phimatchcut_[12];
  unsigned int zmatchcut_[12];

  unsigned int rphicutPS_[12];
  unsigned int rphicut2S_[12];
  unsigned int rcutPS_[12];
  unsigned int rcut2S_[12];

  double phifact_;
  double rzfact_;

  int nrbits_;
  int nphiderbits_;
  
  AllStubsMemory* allstubs_;
  vector<VMStubsMEMemory*> vmstubs_;
  vector<TrackletProjectionsMemory*> inputprojs_;

  int ialphafactinner_[10];
  int ialphafactouter_[10];

  //FIXME should index by iSeed
  vector<FullMatchMemory*> fullmatches_;

  //used in the layers
  vector<bool> table_;

  //used in the disks
  vector<bool> tablePS_;
  vector<bool> table2S_;

  unsigned int nMatchEngines_;  
  vector<MatchEngineUnit> matchengines_;

  CircularBuffer<ProjectionTemp> inputProjBuffer_;

  
};

#endif

//This class implements the tracklet processor
#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletProcessor_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletProcessor_h

#include "IMATH_TrackletCalculator.h"
#include "IMATH_TrackletCalculatorDisk.h"
#include "IMATH_TrackletCalculatorOverlap.h"

#include "TrackletCalculatorBase.h"
#include "VMStubsTEMemory.h"
#include "StubPairsMemory.h"
#include "TrackletProjectionsMemory.h"
#include "GlobalHistTruth.h"
#include "Util.h"


using namespace std;

class TrackletProcessor:public TrackletCalculatorBase{

public:

 TrackletProcessor(string name, const Settings* const settings, unsigned int iSector):
  TrackletCalculatorBase(name,settings,iSector){
    double dphi=2*M_PI/NSector;
    double dphiHG=0.5*settings_->dphisectorHG()-M_PI/NSector;
    phimin_=iSector_*dphi-dphiHG;
    phimax_=phimin_+dphi+2*dphiHG;
    phimin_-=M_PI/NSector;
    phimax_-=M_PI/NSector;
    phimin_=Util::phiRange(phimin_);
    phimax_=Util::phiRange(phimax_);
    if (phimin_>phimax_)  phimin_-=2*M_PI;
    phioffset_=phimin_;
    
    maxtracklet_=127;

    for(unsigned int ilayer=0;ilayer<6;ilayer++){
      vector<TrackletProjectionsMemory*> tmp(settings_->nallstubs(ilayer),0);
      trackletprojlayers_.push_back(tmp);
    }

    for(unsigned int idisk=0;idisk<5;idisk++){
      vector<TrackletProjectionsMemory*> tmp(settings_->nallstubs(idisk+6),0);
      trackletprojdisks_.push_back(tmp);
    }
          
    
    layer_=0;
    disk_=0;

    if (name_[3]=='L') {
      layer_=name_[4]-'0';
      if (name_[5]=='D') disk_=name_[6]-'0';    
    }
    if (name_[3]=='D') disk_=name_[4]-'0';    

    extra_=(layer_==2&&disk_==0);


    // set TC index
    if      (name_[7]=='A') iTC_ =0;
    else if (name_[7]=='B') iTC_ =1;
    else if (name_[7]=='C') iTC_ =2;
    else if (name_[7]=='D') iTC_ =3;
    else if (name_[7]=='E') iTC_ =4;
    else if (name_[7]=='F') iTC_ =5;
    else if (name_[7]=='G') iTC_ =6;
    else if (name_[7]=='H') iTC_ =7;
    else if (name_[7]=='I') iTC_ =8;
    else if (name_[7]=='J') iTC_ =9;
    else if (name_[7]=='K') iTC_ =10;
    else if (name_[7]=='L') iTC_ =11;
    else if (name_[7]=='M') iTC_ =12;
    else if (name_[7]=='N') iTC_ =13;
    else if (name_[7]=='O') iTC_ =14;
    
    assert(iTC_!=-1);
    
    if (name_.substr(3,4)=="L1L2") iSeed_ = 0;
    else if (name_.substr(3,4)=="L3L4") iSeed_ = 2;
    else if (name_.substr(3,4)=="L5L6") iSeed_ = 3;
    else if (name_.substr(3,4)=="D1D2") iSeed_ = 4;
    else if (name_.substr(3,4)=="D3D4") iSeed_ = 5;
    else if (name_.substr(3,4)=="D1L1") iSeed_ = 6;
    else if (name_.substr(3,4)=="D1L2") iSeed_ = 7;
    else if (name_.substr(3,4)=="L1D1") iSeed_ = 6;
    else if (name_.substr(3,4)=="L2D1") iSeed_ = 7;
    else if (name_.substr(3,4)=="L2L3") iSeed_ = 1;
    
    assert(iSeed_!=-1);

    TCIndex_ = (iSeed_<<4) + iTC_;
    assert(TCIndex_>=0 && TCIndex_<128);
   
    assert((layer_!=0)||(disk_!=0));

   
    if (iSeed_==0||iSeed_==1||iSeed_==2||iSeed_==3) {
      if (layer_==1) {
	rproj_[0]=settings_->rmean(2);
	rproj_[1]=settings_->rmean(3);
	rproj_[2]=settings_->rmean(4);
	rproj_[3]=settings_->rmean(5);
	lproj_[0]=3;
	lproj_[1]=4;
	lproj_[2]=5;
	lproj_[3]=6;
      }
      
      if (layer_==2) {
	rproj_[0]=settings_->rmean(0);
	rproj_[1]=settings_->rmean(3);
	rproj_[2]=settings_->rmean(4);
	rproj_[3]=settings_->rmean(5);
	lproj_[0]=1;
	lproj_[1]=4;
	lproj_[2]=5;
	lproj_[3]=6;
      }
      
      if (layer_==3) {
	rproj_[0]=settings_->rmean(0);
	rproj_[1]=settings_->rmean(1);
	rproj_[2]=settings_->rmean(4);
	rproj_[3]=settings_->rmean(5);
	lproj_[0]=1;
	lproj_[1]=2;
	lproj_[2]=5;
	lproj_[3]=6;
      }
      
      if (layer_==5) {
	rproj_[0]=settings_->rmean(0);
	rproj_[1]=settings_->rmean(1);
	rproj_[2]=settings_->rmean(2);
	rproj_[3]=settings_->rmean(3);
	lproj_[0]=1;
	lproj_[1]=2;
	lproj_[2]=3;
	lproj_[3]=4;
      }
    }
    
    if (iSeed_==4||iSeed_==5) {
      if (disk_==1) {
	zproj_[0]=settings_->zmean(2);
	zproj_[1]=settings_->zmean(3);
	zproj_[2]=settings_->zmean(4);
	dproj_[0]=3;
	dproj_[1]=4;
	dproj_[2]=5;
      }
      
      if (disk_==3) {
	zproj_[0]=settings_->zmean(0);
	zproj_[1]=settings_->zmean(1);
	zproj_[2]=settings_->zmean(4);
	dproj_[0]=1;
	dproj_[1]=2;
	dproj_[2]=5;
      }
    }
    
    
    if (iSeed_==6||iSeed_==7) {
      zprojoverlap_[0]=settings_->zmean(1);
      zprojoverlap_[1]=settings_->zmean(2);
      zprojoverlap_[2]=settings_->zmean(3);
      zprojoverlap_[3]=settings_->zmean(4);
    }
    
    if (settings_->usephicritapprox()) {
      double phicritFactor = 0.5 * settings_->rcrit() *GlobalHistTruth::ITC_L1L2()->rinv_final.get_K() /GlobalHistTruth::ITC_L1L2()->phi0_final.get_K();
      if (std::abs(phicritFactor - 2.) > 0.25)
        cout << "TrackletProcessor::TrackletProcessor phicrit approximation may be invalid! Please check." << endl;
    }
  }
  
  void addOutputProjection(TrackletProjectionsMemory* &outputProj, MemoryBase* memory){
      outputProj=dynamic_cast<TrackletProjectionsMemory*>(memory);
      assert(outputProj!=0);
  }
  
  void addOutput(MemoryBase* memory,string output){
    if (settings_->writetrace()) {
      cout << "In "<<name_<<" adding output to "<<memory->getName() << " to output "<<output<<endl;
    }
    if (output=="trackpar"){
      TrackletParametersMemory* tmp=dynamic_cast<TrackletParametersMemory*>(memory);
      assert(tmp!=0);
      trackletpars_=tmp;
      return;
    }

    if (output.substr(0,7)=="projout") {
      //output is on the form 'projoutL2PHIC' or 'projoutD3PHIB'
      TrackletProjectionsMemory* tmp=dynamic_cast<TrackletProjectionsMemory*>(memory);
      assert(tmp!=0);

      unsigned int layerdisk=output[8]-'1'; //layer or disk counting from 0
      unsigned int phiregion=output[12]-'A'; //phiregion counting from 0

      if (output[7]=='L') {
	assert(layerdisk<6);
	assert(phiregion<trackletprojlayers_[layerdisk].size());
	//check that phiregion not already initialized
	assert(trackletprojlayers_[layerdisk][phiregion]==0);
	trackletprojlayers_[layerdisk][phiregion]=tmp;
	return;
      }

      if (output[7]=='D') {
	assert(layerdisk<5);
	assert(phiregion<trackletprojdisks_[layerdisk].size());
	//check that phiregion not already initialized
	assert(trackletprojdisks_[layerdisk][phiregion]==0);
	trackletprojdisks_[layerdisk][phiregion]=tmp;
	return;
      }

    }

    cout << "Could not find output : "<<output<<endl;
    assert(0);


  }

  void addInput(MemoryBase* memory,string input){
    if (settings_->writetrace()) {
      cout << "In "<<name_<<" adding input from "<<memory->getName() << " to input "<<input<<endl;
    }

    if (input=="innervmstubin"){
      VMStubsTEMemory* tmp=dynamic_cast<VMStubsTEMemory*>(memory);
      assert(tmp!=0);
      innervmstubs_.push_back(tmp);
      setVMPhiBin();
      return;
    }
    if (input=="outervmstubin"){
      VMStubsTEMemory* tmp=dynamic_cast<VMStubsTEMemory*>(memory);
      assert(tmp!=0);
      outervmstubs_.push_back(tmp);
      setVMPhiBin();
      return;
    }
    if (input=="innerallstubin"){
      AllStubsMemory* tmp=dynamic_cast<AllStubsMemory*>(memory);
      assert(tmp!=0);
      innerallstubs_.push_back(tmp);
      return;
    }
    if (input=="outerallstubin"){
      AllStubsMemory* tmp=dynamic_cast<AllStubsMemory*>(memory);
      assert(tmp!=0);
      outerallstubs_.push_back(tmp);
      return;
    }
    assert(0);
  }

  void execute() {

    unsigned int countall=0;
    unsigned int countsel=0;

    unsigned int countteall=0;
    unsigned int counttepass=0;

    StubPairsMemory stubpairs("tmp",settings_,iSector_,0.0,1.0); //dummy arguments for now

    bool print=false;
    
    assert(innervmstubs_.size()==outervmstubs_.size());

    if (!settings_->useSeed(iSeed_)) return;
    
    for (unsigned int ivmmem=0;ivmmem<innervmstubs_.size();ivmmem++) {

      unsigned int innerphibin=innervmstubs_[ivmmem]->phibin();
      unsigned int outerphibin=outervmstubs_[ivmmem]->phibin();

      unsigned phiindex=32*innerphibin+outerphibin;
      
      
      //overlap seeding
      if (disk_==1 && (layer_==1 || layer_==2) ) {

	for(unsigned int i=0;i<innervmstubs_[ivmmem]->nVMStubs();i++){

	  VMStubTE innervmstub=innervmstubs_[ivmmem]->getVMStubTE(i); 
	
	  int lookupbits=innervmstub.vmbits().value();

	  int rdiffmax=(lookupbits>>7);	
	  int newbin=(lookupbits&127);
	  int bin=newbin/8;
	  
	  int rbinfirst=newbin&7;

	  int start=(bin>>1);
	  int last=start+(bin&1);
	  
	  for(int ibin=start;ibin<=last;ibin++) {
	    if (settings_->debugTracklet()) cout << getName() << " looking for matching stub in bin "<<ibin
			     <<" with "<<outervmstubs_[ivmmem]->nVMStubsBinned(ibin)<<" stubs"<<endl;
	    for(unsigned int j=0;j<outervmstubs_[ivmmem]->nVMStubsBinned(ibin);j++){
	      //if (countall>=settings_->maxStep("TE")) break;
	      countall++;
	      countteall++;
	      
	      VMStubTE outervmstub=outervmstubs_[ivmmem]->getVMStubTEBinned(ibin,j);
	      int rbin=(outervmstub.vmbits().value()&7);
	      if (start!=ibin) rbin+=8;
	      if ((rbin<rbinfirst)||(rbin-rbinfirst>rdiffmax)) {
		if (settings_->debugTracklet()) {
		  cout << getName() << " layer-disk stub pair rejected because rbin cut : "
		     <<rbin<<" "<<rbinfirst<<" "<<rdiffmax<<endl;
		}
		continue;
	      }

	      int ir=((start&3)<<3)+rbinfirst;
	    
	      assert(innerphibits_!=-1);
	      assert(outerphibits_!=-1);
	      
	      FPGAWord iphiinnerbin=innervmstub.finephi();
	      FPGAWord iphiouterbin=outervmstub.finephi();

	      assert(iphiouterbin==outervmstub.finephi());
	      
	      unsigned int index = (((iphiinnerbin.value()<<outerphibits_)+iphiouterbin.value())<<5)+ir;
	      
	      assert(index<phitable_[phiindex].size());
	      
	      if (!phitable_[phiindex][index]) {
		if (settings_->debugTracklet()) {
		  cout << "Stub pair rejected because of tracklet pt cut"<<endl;
		}
		continue;
	      }
	      
	      FPGAWord innerbend=innervmstub.bend();
	      FPGAWord outerbend=outervmstub.bend();

	      int ptinnerindex=(index<<innerbend.nbits())+innerbend.value();
	      int ptouterindex=(index<<outerbend.nbits())+outerbend.value();
	      
	      if (!(pttableinner_[phiindex][ptinnerindex]&&pttableouter_[phiindex][ptouterindex])) {
		if (settings_->debugTracklet()) {
		  cout << "Stub pair rejected because of stub pt cut bends : "
		       <<Stub::benddecode(innervmstub.bend().value(),innervmstub.isPSmodule())
		       <<" "
		       <<Stub::benddecode(outervmstub.bend().value(),outervmstub.isPSmodule())
		       <<endl;
		}		
		continue;
	      }

	      
	      if (settings_->debugTracklet()) cout << "Adding layer-disk pair in " <<getName()<<endl;
	      if (settings_->writeMonitorData("Seeds")) {
		ofstream fout("seeds.txt", ofstream::app);
		fout << __FILE__ << ":" << __LINE__ << " " << name_ << "_" << iSector_ << " " << iSeed_ << endl;
		fout.close();
	      }
	      stubpairs.addStubPair(innervmstub,outervmstub,0,innervmstubs_[ivmmem]->getName()+" "+outervmstubs_[ivmmem]->getName());
	      counttepass++;
	      countall++;
	    }
	  }
	}

      } else {



	for(unsigned int i=0;i<innervmstubs_[ivmmem]->nVMStubs();i++){
	  if (settings_->debugTracklet()) {
	    cout << "In "<<getName()<<" have inner stub"<<endl;
	  }

	  
	  if ((layer_==1 && disk_==0)||
	      (layer_==2 && disk_==0)||
	      (layer_==3 && disk_==0)||
	      (layer_==5 && disk_==0)) {	  

	    

	    VMStubTE innervmstub=innervmstubs_[ivmmem]->getVMStubTE(i);

	    
	    int lookupbits=(int)innervmstub.vmbits().value();
	    int zdiffmax=(lookupbits>>7);	
	    int newbin=(lookupbits&127);
	    int bin=newbin/8;
	    
	    int zbinfirst=newbin&7;
	    
	    int start=(bin>>1);
	    int last=start+(bin&1);
	    
	    if (print) {
	      cout << "start last : "<<start<<" "<<last<<endl;
	    }
	    
	    if (settings_->debugTracklet()) {
	      cout << "Will look in zbins "<<start<<" to "<<last<<endl;
	    }
	    for(int ibin=start;ibin<=last;ibin++) {
	      for(unsigned int j=0;j<outervmstubs_[ivmmem]->nVMStubsBinned(ibin);j++){
		if (settings_->debugTracklet()) {
		  cout << "In "<<getName()<<" have outer stub"<<endl;
		}
		
		//if (countall>=settings_->maxStep("TE")) break;
		countteall++;
		countall++;

		VMStubTE outervmstub=outervmstubs_[ivmmem]->getVMStubTEBinned(ibin,j);

		int zbin=(outervmstub.vmbits().value()&7);
		
		if (start!=ibin) zbin+=8;
		
		if (zbin<zbinfirst||zbin-zbinfirst>zdiffmax) {
		  if (settings_->debugTracklet()) {
		    cout << "Stubpair rejected because of wrong fine z"<<endl;
		  }
		  continue;
		}

		if (print) {
		  cout << "ibin j "<<ibin<<" "<<j<<endl;
		}
	      
		assert(innerphibits_!=-1);
		assert(outerphibits_!=-1);
	      
		FPGAWord iphiinnerbin=innervmstub.finephi();
		FPGAWord iphiouterbin=outervmstub.finephi();
		
		int index = (iphiinnerbin.value()<<outerphibits_)+iphiouterbin.value();
		
		assert(index<(int)phitable_[phiindex].size());		
		
		
		if (!phitable_[phiindex][index]) {
		  if (settings_->debugTracklet()) {
		    cout << "Stub pair rejected because of tracklet pt cut"<<endl;
		  }
		  continue;
		}
		
		FPGAWord innerbend=innervmstub.bend();
		FPGAWord outerbend=outervmstub.bend();

		
		int ptinnerindex=(index<<innerbend.nbits())+innerbend.value();
		int ptouterindex=(index<<outerbend.nbits())+outerbend.value();
		
		
		if (!(pttableinner_[phiindex][ptinnerindex]&&pttableouter_[phiindex][ptouterindex])) {
		  if (settings_->debugTracklet()) {
		    cout << "Stub pair rejected because of stub pt cut bends : "
			 <<Stub::benddecode(innervmstub.bend().value(),innervmstub.isPSmodule())
			 <<" "
			 <<Stub::benddecode(outervmstub.bend().value(),outervmstub.isPSmodule())
			 <<endl;
		  }		
		  continue;
		}
		
		if (settings_->debugTracklet()) cout << "Adding layer-layer pair in " <<getName()<<endl;
		if (settings_->writeMonitorData("Seeds")) {
		  ofstream fout("seeds.txt", ofstream::app);
		  fout << __FILE__ << ":" << __LINE__ << " " << name_ << "_" << iSector_ << " " << iSeed_ << endl;
		  fout.close();
		}
		stubpairs.addStubPair(innervmstub,outervmstub,0,innervmstubs_[ivmmem]->getName()+" "+outervmstubs_[ivmmem]->getName());
		counttepass++;	      
		countall++;
	      }
	    }
	    
	  } else if ((disk_==1 && layer_==0)||
		     (disk_==3 && layer_==0)) {
	    
	    if (settings_->debugTracklet()) cout << getName()<<"["<<iSector_<<"] Disk-disk pair" <<endl;

	    VMStubTE innervmstub=innervmstubs_[ivmmem]->getVMStubTE(i);

	    
	    int lookupbits=(int)innervmstub.vmbits().value();
	    bool negdisk=innervmstub.stub().first->disk().value()<0; //FIXME
	    int rdiffmax=(lookupbits>>6);	
	    int newbin=(lookupbits&63);
	    int bin=newbin/8;
	    
	    int rbinfirst=newbin&7;

	    int start=(bin>>1);
	    if (negdisk) start+=4;
	    int last=start+(bin&1);
	    for(int ibin=start;ibin<=last;ibin++) {
	      if (settings_->debugTracklet()) cout << getName() << " looking for matching stub in bin "<<ibin
			       <<" with "<<outervmstubs_[ivmmem]->nVMStubsBinned(ibin)<<" stubs"<<endl;
	      for(unsigned int j=0;j<outervmstubs_[ivmmem]->nVMStubsBinned(ibin);j++){
		//if (countall>=settings_->maxStep("TE")) break;
		countall++;
		countteall++;
		
		VMStubTE outervmstub=outervmstubs_[ivmmem]->getVMStubTEBinned(ibin,j);

		int vmbits=(int)outervmstub.vmbits().value();
		int rbin=(vmbits&7);
		if (start!=ibin) rbin+=8;
		if (rbin<rbinfirst) continue;
		if (rbin-rbinfirst>rdiffmax) continue;
		
	      
		unsigned int irouterbin=vmbits>>2;
		
		FPGAWord iphiinnerbin=innervmstub.finephi();
		FPGAWord iphiouterbin=outervmstub.finephi();
	      	      
		unsigned int index = (irouterbin<<(outerphibits_+innerphibits_))+(iphiinnerbin.value()<<outerphibits_)+iphiouterbin.value();
		
		assert(index<phitable_[phiindex].size());		
		if (!phitable_[phiindex][index]) {
		  if (settings_->debugTracklet()) {
		    cout << "Stub pair rejected because of tracklet pt cut"<<endl;
		  }
		  continue;
		}
		
		FPGAWord innerbend=innervmstub.bend();
		FPGAWord outerbend=outervmstub.bend();
		
		unsigned int ptinnerindex=(index<<innerbend.nbits())+innerbend.value();
		unsigned int ptouterindex=(index<<outerbend.nbits())+outerbend.value();
	      
		assert(ptinnerindex<pttableinner_[phiindex].size());
		assert(ptouterindex<pttableouter_[phiindex].size());
	      
		if (!(pttableinner_[phiindex][ptinnerindex]&&pttableouter_[phiindex][ptouterindex])) {
		  if (settings_->debugTracklet()) {
		    cout << "Stub pair rejected because of stub pt cut bends : "
			 <<Stub::benddecode(innervmstub.bend().value(),innervmstub.isPSmodule())
			 <<" "
			 <<Stub::benddecode(outervmstub.bend().value(),outervmstub.isPSmodule())
			 <<" pass : "<<pttableinner_[phiindex][ptinnerindex]<<" "<<pttableouter_[phiindex][ptouterindex]
			 <<endl;
		  }
		  continue;
		}

		if (settings_->debugTracklet()) cout << "Adding disk-disk pair in " <<getName()<<endl;
	      
		if (settings_->writeMonitorData("Seeds")) {
		  ofstream fout("seeds.txt", ofstream::app);
		  fout << __FILE__ << ":" << __LINE__ << " " << name_ << "_" << iSector_ << " " << iSeed_ << endl;
		  fout.close();
		}
		stubpairs.addStubPair(innervmstub,outervmstub,0,innervmstubs_[ivmmem]->getName()+" "+outervmstubs_[ivmmem]->getName());
		countall++;
		counttepass++;
	
	      }
	    }
	  }
	}
      
      }
    }

    if (settings_->writeMonitorData("TE")) {
      static ofstream out("trackletprocessor.txt");
      out << getName()<<" "<<countteall<<" "<<counttepass<<endl;
    }

    for(unsigned int i=0;i<stubpairs.nStubPairs();i++){

      if (trackletpars_->nTracklets()>=maxtracklet_) {
	cout << "Will break on too many tracklets in "<<getName()<<endl;
	break;
      }
      countall++;
      L1TStub* innerStub=stubpairs.getL1TStub1(i);
      Stub* innerFPGAStub=stubpairs.getFPGAStub1(i);
      
      L1TStub* outerStub=stubpairs.getL1TStub2(i);
      Stub* outerFPGAStub=stubpairs.getFPGAStub2(i);
      
      if (settings_->debugTracklet()) {
	cout << "TrackletProcessor execute "<<getName()<<"["<<iSector_<<"]"<<endl;
      }

      bool accept=false;
      
      if (innerFPGAStub->isBarrel()&&(getName()!="TC_D1L2A"&&getName()!="TC_D1L2B")){
	
	if (outerFPGAStub->isDisk()) {
	  
	  //overlap seeding                                              
	  accept = overlapSeeding(outerFPGAStub,outerStub,innerFPGAStub,innerStub);
	  
	} else {
	  
	  //barrel+barrel seeding	  
	    accept = barrelSeeding(innerFPGAStub,innerStub,outerFPGAStub,outerStub);
	}
	
      }  else {
	
	if (outerFPGAStub->isDisk()) {
	  
	  //disk+disk seeding
	  
	  accept = diskSeeding(innerFPGAStub,innerStub,outerFPGAStub,outerStub);
	  
	} else if (innerFPGAStub->isDisk()) {
	  
	  //layer+disk seeding
	  
	  accept = overlapSeeding(innerFPGAStub,innerStub,outerFPGAStub,outerStub);
	  
	} else {
	  
	  assert(0);
	  
	}
      }

      if (accept) countsel++;
      
      if (settings_->writeMonitorData("TP")) {
	static ofstream out("tc_seedpairs.txt");
	out << stubpairs.getTEDName(i)<<" "<<accept<<endl;
      }

      if (trackletpars_->nTracklets()>=maxtracklet_) {
	cout << "Will break on number of tracklets in "<<getName()<<endl;
	break;
      }
      
      if (countall>=settings_->maxStep("TP")) {
	if (settings_->debugTracklet()) cout << "Will break on MAXTC 1"<<endl;
	break;
      }
      if (settings_->debugTracklet()) {
	cout << "TrackletProcessor execute done"<<endl;
      }
      
    }
    if (countall>=settings_->maxStep("TP")) {
      if (settings_->debugTracklet()) cout << "Will break on MAXTC 2"<<endl;
      //break;
    }
      
    
    if (settings_->writeMonitorData("TP")) {
      static ofstream out("trackletcalculator.txt");
      out << getName()<<" "<<countall<<" "<<countsel<<endl;
    }
        
  }


  void setVMPhiBin() {
    
    if (innervmstubs_.size()!=outervmstubs_.size() ) return;

    for(unsigned int ivmmem=0;ivmmem<innervmstubs_.size();ivmmem++){
    
      unsigned int innerphibin=innervmstubs_[ivmmem]->phibin();
      unsigned int outerphibin=outervmstubs_[ivmmem]->phibin();
      
      unsigned phiindex=32*innerphibin+outerphibin;


      if (phitable_.find(phiindex)!=phitable_.end()) continue;
      
      innervmstubs_[ivmmem]->setother(outervmstubs_[ivmmem]);
      outervmstubs_[ivmmem]->setother(innervmstubs_[ivmmem]);


    
      if ((layer_==1 && disk_==0)||
	  (layer_==2 && disk_==0)||
	  (layer_==3 && disk_==0)||
	  (layer_==5 && disk_==0)){
      
	innerphibits_=settings_->nfinephi(0,iSeed_);
	outerphibits_=settings_->nfinephi(1,iSeed_);
	
	int innerphibins=(1<<innerphibits_);
	int outerphibins=(1<<outerphibits_);
	
	double innerphimin, innerphimax;
	innervmstubs_[ivmmem]->getPhiRange(innerphimin,innerphimax,iSeed_,0);
	double rinner=settings_->rmean(layer_-1);
	
	double outerphimin, outerphimax;
	outervmstubs_[ivmmem]->getPhiRange(outerphimin,outerphimax,iSeed_,1);
	double router=settings_->rmean(layer_);
	
	double phiinner[2];
	double phiouter[2];
	
	std::vector<bool> vmbendinner;
	std::vector<bool> vmbendouter;
	unsigned int nbins1=8;
	if (layer_>=4) nbins1=16;
	for (unsigned int i=0;i<nbins1;i++) {
	  vmbendinner.push_back(false);
	}
	
	unsigned int nbins2=8;
	if (layer_>=3) nbins2=16;
	for (unsigned int i=0;i<nbins2;i++) {
	  vmbendouter.push_back(false);
	}
	
	for (int iphiinnerbin=0;iphiinnerbin<innerphibins;iphiinnerbin++){
	  phiinner[0]=innerphimin+iphiinnerbin*(innerphimax-innerphimin)/innerphibins;
	  phiinner[1]=innerphimin+(iphiinnerbin+1)*(innerphimax-innerphimin)/innerphibins;
	  for (int iphiouterbin=0;iphiouterbin<outerphibins;iphiouterbin++){
	    phiouter[0]=outerphimin+iphiouterbin*(outerphimax-outerphimin)/outerphibins;
	    phiouter[1]=outerphimin+(iphiouterbin+1)*(outerphimax-outerphimin)/outerphibins;
	    
	    double bendinnermin=20.0;
	    double bendinnermax=-20.0;
	    double bendoutermin=20.0;
	    double bendoutermax=-20.0;
	    double rinvmin=1.0; 
	    for(int i1=0;i1<2;i1++) {
	      for(int i2=0;i2<2;i2++) {
		double rinv1=rinv(phiinner[i1],phiouter[i2],rinner,router);
		double abendinner=bend(rinner,rinv1); 
		double abendouter=bend(router,rinv1);
		if (abendinner<bendinnermin) bendinnermin=abendinner;
		if (abendinner>bendinnermax) bendinnermax=abendinner;
		if (abendouter<bendoutermin) bendoutermin=abendouter;
		if (abendouter>bendoutermax) bendoutermax=abendouter;
		if (std::abs(rinv1)<rinvmin) {
		  rinvmin=std::abs(rinv1);
		}
		
	      }
	    }
	    
	    phitable_[phiindex].push_back(rinvmin<rinvcutte);
	    
	    int nbins1=8;
	    if (layer_>=4) nbins1=16;
	    for(int ibend=0;ibend<nbins1;ibend++) {
	      double bend=Stub::benddecode(ibend,layer_<=3); 
	      
	      bool passinner=bend-bendinnermin>-bendcut&&bend-bendinnermax<bendcut;	    
	      if (passinner) vmbendinner[ibend]=true;
	      pttableinner_[phiindex].push_back(passinner);
	      
	    }
	    
	    int nbins2=8;
	    if (layer_>=3) nbins2=16;
	    for(int ibend=0;ibend<nbins2;ibend++) {
	      double bend=Stub::benddecode(ibend,layer_<=2); 
	      
	      bool passouter=bend-bendoutermin>-bendcut&&bend-bendoutermax<bendcut;
	      if (passouter) vmbendouter[ibend]=true;
	      pttableouter_[phiindex].push_back(passouter);
	      
	    }   
	  }
	}

	innervmstubs_[ivmmem]->setbendtable(settings_,vmbendinner);
	outervmstubs_[ivmmem]->setbendtable(settings_,vmbendouter);
      
	if (iSector_==0&&settings_->writeTable()) writeTETable();
      
      }

      if ((disk_==1 && layer_==0)||
	  (disk_==3 && layer_==0)){
	
	innerphibits_=settings_->nfinephi(0,iSeed_);
	outerphibits_=settings_->nfinephi(0,iSeed_);
	
	
	int outerrbits=3;
	
	int outerrbins=(1<<outerrbits);
	int innerphibins=(1<<innerphibits_);
	int outerphibins=(1<<outerphibits_);
	
	double innerphimin, innerphimax;
	innervmstubs_[ivmmem]->getPhiRange(innerphimin,innerphimax,iSeed_,0);
	
	double outerphimin, outerphimax;
	outervmstubs_[ivmmem]->getPhiRange(outerphimin,outerphimax,iSeed_,1);
	
	
	double phiinner[2];
	double phiouter[2];
	double router[2];
	
	std::vector<bool> vmbendinner;
	std::vector<bool> vmbendouter;
	
	for (unsigned int i=0;i<8;i++) {
	  vmbendinner.push_back(false);
	  vmbendouter.push_back(false);
	}
	

	for (int irouterbin=0;irouterbin<outerrbins;irouterbin++){
	  router[0]=rmindiskvm+irouterbin*(rmaxdiskvm-rmindiskvm)/outerrbins;
	  router[1]=rmindiskvm+(irouterbin+1)*(rmaxdiskvm-rmindiskvm)/outerrbins;
	  for (int iphiinnerbin=0;iphiinnerbin<innerphibins;iphiinnerbin++){
	    phiinner[0]=innerphimin+iphiinnerbin*(innerphimax-innerphimin)/innerphibins;
	    phiinner[1]=innerphimin+(iphiinnerbin+1)*(innerphimax-innerphimin)/innerphibins;
	    for (int iphiouterbin=0;iphiouterbin<outerphibins;iphiouterbin++){
	      phiouter[0]=outerphimin+iphiouterbin*(outerphimax-outerphimin)/outerphibins;
	      phiouter[1]=outerphimin+(iphiouterbin+1)*(outerphimax-outerphimin)/outerphibins;
	      
	      double bendinnermin=20.0;
	      double bendinnermax=-20.0;
	      double bendoutermin=20.0;
	      double bendoutermax=-20.0;
	      double rinvmin=1.0; 
	      double rinvmax=-1.0; 
	      for(int i1=0;i1<2;i1++) {
		for(int i2=0;i2<2;i2++) {
		  for(int i3=0;i3<2;i3++) {
		    double rinner=router[i3]*settings_->zmean(disk_-1)/settings_->zmean(disk_);
		    double rinv1=rinv(phiinner[i1],phiouter[i2],rinner,router[i3]);
		    double abendinner=bend(rinner,rinv1);
		    double abendouter=bend(router[i3],rinv1);
		    if (abendinner<bendinnermin) bendinnermin=abendinner;
		    if (abendinner>bendinnermax) bendinnermax=abendinner;
		    if (abendouter<bendoutermin) bendoutermin=abendouter;
		    if (abendouter>bendoutermax) bendoutermax=abendouter;
		    if (std::abs(rinv1)<rinvmin) {
		      rinvmin=std::abs(rinv1);
		    }
		    if (std::abs(rinv1)>rinvmax) {
		      rinvmax=std::abs(rinv1);
		    }
		  }
		}
	      }
	      
	      phitable_[phiindex].push_back(rinvmin<rinvcutte);


	      for(int ibend=0;ibend<8;ibend++) {
		double bend=Stub::benddecode(ibend,true); 
		
		bool passinner=bend-bendinnermin>-bendcutdisk&&bend-bendinnermax<bendcutdisk;	    
		if (passinner) vmbendinner[ibend]=true;
		pttableinner_[phiindex].push_back(passinner);
		
	      }
	      
	      for(int ibend=0;ibend<8;ibend++) {
		double bend=Stub::benddecode(ibend,true); 
		
		bool passouter=bend-bendoutermin>-bendcutdisk&&bend-bendoutermax<bendcutdisk;
		if (passouter) vmbendouter[ibend]=true;
		pttableouter_[phiindex].push_back(passouter);
		
	      }
	    
	    }
	  }
	}
	
	innervmstubs_[ivmmem]->setbendtable(settings_,vmbendinner);
	outervmstubs_[ivmmem]->setbendtable(settings_,vmbendouter);
	
	if (iSector_==0&&settings_->writeTable()) writeTETable();
	
      } else if (disk_==1 && (layer_==1 || layer_==2)) {
	
	innerphibits_=settings_->nfinephi(0,iSeed_);
	outerphibits_=settings_->nfinephi(1,iSeed_);
	unsigned int nrbits=5;
	
	int innerphibins=(1<<innerphibits_);
	int outerphibins=(1<<outerphibits_);
	
	double innerphimin, innerphimax;
	innervmstubs_[ivmmem]->getPhiRange(innerphimin,innerphimax,iSeed_,0);
	
	double outerphimin, outerphimax;
	outervmstubs_[ivmmem]->getPhiRange(outerphimin,outerphimax,iSeed_,1);
	
	double phiinner[2];
	double phiouter[2];
	double router[2];
	

	std::vector<bool> vmbendinner;
	std::vector<bool> vmbendouter;
	
	for (unsigned int i=0;i<8;i++) {
	  vmbendinner.push_back(false);
	  vmbendouter.push_back(false);
	}
	
	double dr=(rmaxdiskvm-rmindiskvm)/(1<<nrbits);

	for (int iphiinnerbin=0;iphiinnerbin<innerphibins;iphiinnerbin++){
	  phiinner[0]=innerphimin+iphiinnerbin*(innerphimax-innerphimin)/innerphibins;
	  phiinner[1]=innerphimin+(iphiinnerbin+1)*(innerphimax-innerphimin)/innerphibins;
	  for (int iphiouterbin=0;iphiouterbin<outerphibins;iphiouterbin++){
	    phiouter[0]=outerphimin+iphiouterbin*(outerphimax-outerphimin)/outerphibins;
	    phiouter[1]=outerphimin+(iphiouterbin+1)*(outerphimax-outerphimin)/outerphibins;
	    for (int irbin=0;irbin<(1<<nrbits);irbin++){
	      router[0]=rmindiskvm+dr*irbin;
	      router[1]=router[0]+dr; 
	      double bendinnermin=20.0;
	      double bendinnermax=-20.0;
	      double bendoutermin=20.0;
	      double bendoutermax=-20.0;
	      double rinvmin=1.0; 
	      for(int i1=0;i1<2;i1++) {
		for(int i2=0;i2<2;i2++) {
		  for(int i3=0;i3<2;i3++) {
		    double rinner=settings_->rmean(layer_-1);
		    double rinv1=rinv(phiinner[i1],phiouter[i2],rinner,router[i3]);
		    double abendinner=bend(rinner,rinv1);
		    double abendouter=bend(router[i3],rinv1);
		    if (abendinner<bendinnermin) bendinnermin=abendinner;
		    if (abendinner>bendinnermax) bendinnermax=abendinner;
		    if (abendouter<bendoutermin) bendoutermin=abendouter;
		    if (abendouter>bendoutermax) bendoutermax=abendouter;
		    if (std::abs(rinv1)<rinvmin) {
		      rinvmin=std::abs(rinv1);
		    }
		  }
		}
	      }
	    
	      phitable_[phiindex].push_back(rinvmin<rinvcutte);
	      
	      
	      for(int ibend=0;ibend<8;ibend++) {
		double bend=Stub::benddecode(ibend,true); 
		
		bool passinner=bend-bendinnermin>-bendcut&&bend-bendinnermax<bendcut;	    
		if (passinner) vmbendinner[ibend]=true;
		pttableinner_[phiindex].push_back(passinner);
		
	      }
	      
	      for(int ibend=0;ibend<8;ibend++) {
		double bend=Stub::benddecode(ibend,true); 
		
		bool passouter=bend-bendoutermin>-bendcut&&bend-bendoutermax<bendcut;
		if (passouter) vmbendouter[ibend]=true;
		pttableouter_[phiindex].push_back(passouter);
		
	      }
	    }
	  }
	}
	
    
	innervmstubs_[ivmmem]->setbendtable(settings_,vmbendinner);
	outervmstubs_[ivmmem]->setbendtable(settings_,vmbendouter);
	
	if (iSector_==0&&settings_->writeTable()) writeTETable();
	
	
      }
    } 
  }
    
  double rinv(double phi1, double phi2,double r1, double r2){
    
    if (r2<r1) { //can not form tracklet
      return 20.0; 
    }
    
    assert(r2>r1);

    double dphi=phi2-phi1;
    double dr=r2-r1;
    
    return 2.0*sin(dphi)/dr/sqrt(1.0+2*r1*r2*(1.0-cos(dphi))/(dr*dr));
    
  }

  double bend(double r, double rinv) {

    double dr=0.18;
    
    double delta=r*dr*0.5*rinv;

    double bend=delta/0.009;
    if (r<55.0) bend=delta/0.01;

    return bend;
    
  }


    void writeTETable() {

    ofstream outptcut;
    outptcut.open(getName()+"_ptcut.tab");
    outptcut << "{"<<endl;
    //for(unsigned int i=0;i<phitable_.size();i++){
    //  if (i!=0) outptcut<<","<<endl;
    //  outptcut << phitable_[i];
    //}
    outptcut <<endl<<"};"<<endl;
    outptcut.close();

    ofstream outstubptinnercut;
    outstubptinnercut.open(getName()+"_stubptinnercut.tab");
    outstubptinnercut << "{"<<endl;
    //for(unsigned int i=0;i<pttableinner_.size();i++){
    //  if (i!=0) outstubptinnercut<<","<<endl;
    //  outstubptinnercut << pttableinner_[i];
    //}
    outstubptinnercut <<endl<<"};"<<endl;
    outstubptinnercut.close();
    
    ofstream outstubptoutercut;
    outstubptoutercut.open(getName()+"_stubptoutercut.tab");
    outstubptoutercut << "{"<<endl;
    //for(unsigned int i=0;i<pttableouter_.size();i++){
    //  if (i!=0) outstubptoutercut<<","<<endl;
    //  outstubptoutercut << pttableouter_[i];
    //}
    outstubptoutercut <<endl<<"};"<<endl;
    outstubptoutercut.close();

    
  }
  

  
    
private:
  
  int iTC_;

  unsigned int maxtracklet_; //maximum numbor of tracklets that be stored
  

  vector<VMStubsTEMemory*> innervmstubs_;
  vector<VMStubsTEMemory*> outervmstubs_;
  
  vector<AllStubsMemory*> innerallstubs_;
  vector<AllStubsMemory*> outerallstubs_;


  bool extra_;

  
  map<unsigned int, vector<bool> > phitable_;
  map<unsigned int, vector<bool> > pttableinner_;
  map<unsigned int, vector<bool> > pttableouter_;
  
  int innerphibits_;
  int outerphibits_;


public:

};


#endif

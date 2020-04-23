//This class implementes the tracklet engine
#ifndef TRACKLETCALCULATOR_H
#define TRACKLETCALCULATOR_H

#include "IMATH_TrackletCalculator.h"
#include "IMATH_TrackletCalculatorDisk.h"
#include "IMATH_TrackletCalculatorOverlap.h"

#include "TrackletCalculatorBase.h"
#include "TrackletProjectionsMemory.h"
#include "Util.h"


using namespace std;

class TrackletCalculator:public TrackletCalculatorBase{

public:

 TrackletCalculator(string name, const Settings* const settings, unsigned int iSector):
  TrackletCalculatorBase(name, settings, iSector){
      double dphi=2*M_PI/NSector;
      double dphiHG=0.5*dphisectorHG-M_PI/NSector;
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
	vector<TrackletProjectionsMemory*> tmp(nallstubslayers[ilayer],0);
	trackletprojlayers_.push_back(tmp);
      }

      for(unsigned int idisk=0;idisk<5;idisk++){
	vector<TrackletProjectionsMemory*> tmp(nallstubsdisks[idisk],0);
	trackletprojdisks_.push_back(tmp);
      }
          
    
      layer_=0;
      disk_=0;
    
      if (name_[3]=='L') layer_=name_[4]-'0';    
      if (name_[3]=='D') disk_=name_[4]-'0';    
    

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
	  rproj_[0]=rmeanL3;
	  rproj_[1]=rmeanL4;
	  rproj_[2]=rmeanL5;
	  rproj_[3]=rmeanL6;
	  lproj_[0]=3;
	  lproj_[1]=4;
	  lproj_[2]=5;
	  lproj_[3]=6;
	}
	
	if (layer_==2) {
	  rproj_[0]=rmeanL1;
	  rproj_[1]=rmeanL4;
	  rproj_[2]=rmeanL5;
	  rproj_[3]=rmeanL6;
	  lproj_[0]=1;
	  lproj_[1]=4;
	  lproj_[2]=5;
	  lproj_[3]=6;
	}
      
      if (layer_==3) {
	rproj_[0]=rmeanL1;
	rproj_[1]=rmeanL2;
	rproj_[2]=rmeanL5;
	rproj_[3]=rmeanL6;
	lproj_[0]=1;
	lproj_[1]=2;
	lproj_[2]=5;
	lproj_[3]=6;
      }
      
      if (layer_==5) {
	rproj_[0]=rmeanL1;
	rproj_[1]=rmeanL2;
	rproj_[2]=rmeanL3;
	rproj_[3]=rmeanL4;
	lproj_[0]=1;
	lproj_[1]=2;
	lproj_[2]=3;
	lproj_[3]=4;
      }
    }
    
    if (iSeed_==4||iSeed_==5) {
      if (disk_==1) {
	zproj_[0]=zmeanD3;
	zproj_[1]=zmeanD4;
	zproj_[2]=zmeanD5;
	dproj_[0]=3;
	dproj_[1]=4;
	dproj_[2]=5;
      }
      
      if (disk_==3) {
	zproj_[0]=zmeanD1;
	zproj_[1]=zmeanD2;
	zproj_[2]=zmeanD5;
	dproj_[0]=1;
	dproj_[1]=2;
	dproj_[2]=5;
      }
    }
    
    
    if (iSeed_==6||iSeed_==7) {
      zprojoverlap_[0]=zmeanD2;
      zprojoverlap_[1]=zmeanD3;
      zprojoverlap_[2]=zmeanD4;
      zprojoverlap_[3]=zmeanD5;
    }
    
    if (usephicritapprox) {
      double phicritFactor = 0.5 * rcrit * TrackletCalculator::ITC_L1L2.rinv_final.get_K() / TrackletCalculator::ITC_L1L2.phi0_final.get_K();
      if (std::abs(phicritFactor - 2.) > 0.25)
        cout << "TrackletCalculator::TrackletCalculator phicrit approximation may be invalid! Please check." << endl;
    }
  }
  
  void addOutputProjection(TrackletProjectionsMemory* &outputProj, MemoryBase* memory){
      outputProj=dynamic_cast<TrackletProjectionsMemory*>(memory);
      assert(outputProj!=0);
  }
  
  void addOutput(MemoryBase* memory,string output){
    if (writetrace) {
      cout << "In "<<name_<<" adding output to "<<memory->getName()
	   << " to output "<<output<<endl;
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
    if (writetrace) {
      cout << "In "<<name_<<" adding input from "<<memory->getName()
	   << " to input "<<input<<endl;
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
    if (input.substr(0,8)=="stubpair"){
      StubPairsMemory* tmp=dynamic_cast<StubPairsMemory*>(memory);
      assert(tmp!=0);
      stubpairs_.push_back(tmp);
      return;
    }
    assert(0);
  }

  void execute() {

    unsigned int countall=0;
    unsigned int countsel=0;

    for(unsigned int l=0;l<stubpairs_.size();l++){
      if (trackletpars_->nTracklets()>=maxtracklet_) {
	cout << "Will break on too many tracklets in "<<getName()<<endl;
	break;
      }
      for(unsigned int i=0;i<stubpairs_[l]->nStubPairs();i++){

	countall++;
	L1TStub* innerStub=stubpairs_[l]->getL1TStub1(i);
	Stub* innerFPGAStub=stubpairs_[l]->getFPGAStub1(i);

	L1TStub* outerStub=stubpairs_[l]->getL1TStub2(i);
	Stub* outerFPGAStub=stubpairs_[l]->getFPGAStub2(i);

	if (debug1) {
	  cout << "TrackletCalculator execute "<<getName()<<"["<<iSector_<<"]"<<endl;
	}
	
	if (innerFPGAStub->isBarrel()&&(getName()!="TC_D1L2A"&&getName()!="TC_D1L2B")){


          if (outerFPGAStub->isDisk()) {

            //overlap seeding                                              
            bool accept = overlapSeeding(outerFPGAStub,outerStub,innerFPGAStub,innerStub);
            if (accept) countsel++;

          } else {
	  
	    //barrel+barrel seeding	  
	    bool accept = barrelSeeding(innerFPGAStub,innerStub,outerFPGAStub,outerStub);
	    
	    if (accept) countsel++;

	  }
	  
	}  else {

	  if (outerFPGAStub->isDisk()) {

	    //disk+disk seeding

	    bool accept = diskSeeding(innerFPGAStub,innerStub,outerFPGAStub,outerStub);

	    if (accept) countsel++;
	    

	  } else if (innerFPGAStub->isDisk()) {


	    //layer+disk seeding
	    
	    bool accept = overlapSeeding(innerFPGAStub,innerStub,outerFPGAStub,outerStub);

	    if (accept) countsel++;


	  } else {

	    assert(0);
	    
	  }
	}

	if (trackletpars_->nTracklets()>=maxtracklet_) {
	  cout << "Will break on number of tracklets in "<<getName()<<endl;
	  break;
	}
	
	if (countall>=MAXTC) {
	  if (debug1) cout << "Will break on MAXTC 1"<<endl;
	  break;
	}
	if (debug1) {
	  cout << "TrackletCalculator execute done"<<endl;
	}

      }
      if (countall>=MAXTC) {
	if (debug1) cout << "Will break on MAXTC 2"<<endl;
	break;
      }
    }

    if (writeTrackletCalculator) {
      static ofstream out("trackletcalculator.txt");
      out << getName()<<" "<<countall<<" "<<countsel<<endl;
    }

  }

private:
  
  int iTC_;

  unsigned int maxtracklet_; //maximum numbor of tracklets that be stored
  

  vector<AllStubsMemory*> innerallstubs_;
  vector<AllStubsMemory*> outerallstubs_;
  vector<StubPairsMemory*> stubpairs_;


public:
    
};



#endif

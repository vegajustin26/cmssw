#include "L1Trigger/TrackFindingTracklet/interface/MatchEngine.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"
#include "L1Trigger/TrackFindingTracklet/interface/Util.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std; 
using namespace Trklet; 

MatchEngine::MatchEngine(string name, const Settings* settings, Globals* global, unsigned int iSector):
  ProcessBase(name,settings,global,iSector){
    layer_=0;
    disk_=0;
    string subname=name.substr(3,2);
    if (subname=="L1") layer_=1;
    if (subname=="L2") layer_=2;
    if (subname=="L3") layer_=3;
    if (subname=="L4") layer_=4;
    if (subname=="L5") layer_=5;
    if (subname=="L6") layer_=6;
    if (subname=="D1") disk_=1;
    if (subname=="D2") disk_=2;
    if (subname=="D3") disk_=3;
    if (subname=="D4") disk_=4;
    if (subname=="D5") disk_=5;
    if (layer_==0&&disk_==0) {
      edm::LogPrint("Tracklet") << name<<" subname = "<<subname<<" "<<layer_<<" "<<disk_;
    }
    assert((layer_!=0)||(disk_!=0));

    if (layer_>0) {

      unsigned int nbits=3;
      if (layer_>=4) nbits=4;
      
      for(unsigned int irinv=0;irinv<32;irinv++){
	double rinv=(irinv-15.5)*(1<<(settings_->nbitsrinv()-5))*settings_->krinvpars();
	double projbend=bend(settings_->rmean(layer_-1),rinv);
	for(unsigned int ibend=0;ibend<(unsigned int)(1<<nbits);ibend++){
	  double stubbend=benddecode(ibend,layer_<=3);
	  bool pass=std::abs(stubbend-projbend)<settings_->bendcutme(layer_-1);
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
	  double stubbend=benddecode(ibend,true);
	  bool pass=std::abs(stubbend-projbend)<settings_->bendcutme(disk_+5);
	  tablePS_.push_back(pass);
	}
	for(unsigned int ibend=0;ibend<16;ibend++){
	  double stubbend=benddecode(ibend,false);
	  bool pass=std::abs(stubbend-projbend)<settings_->bendcutme(disk_+5);
	  table2S_.push_back(pass);
	}
      }
      
    }
    
}

void MatchEngine::addOutput(MemoryBase* memory,string output) {
  if (settings_->writetrace()) {
    edm::LogVerbatim("Tracklet") << "In "<<name_<<" adding output to "<<memory->getName() << " to output "<<output;
  }
  if (output=="matchout") {
    CandidateMatchMemory* tmp=dynamic_cast<CandidateMatchMemory*>(memory);
    assert(tmp!=0);
    candmatches_=tmp;
    return;
  }
  assert(0);
  
}

void MatchEngine::addInput(MemoryBase* memory,string input) {
  if (settings_->writetrace()) {
    edm::LogVerbatim("Tracklet") << "In "<<name_<<" adding input from "<<memory->getName() << " to input "<<input;
  }
  if (input=="vmstubin") {
    VMStubsMEMemory* tmp=dynamic_cast<VMStubsMEMemory*>(memory);
    assert(tmp!=0);
    vmstubs_=tmp;
    return;
  }
  if (input=="vmprojin") {
    VMProjectionsMemory* tmp=dynamic_cast<VMProjectionsMemory*>(memory);
    assert(tmp!=0);
    vmprojs_=tmp;
    return;
  }
  edm::LogPrint("Tracklet") << "Could not find input : "<<input;
  assert(0);
}

void MatchEngine::execute() {

    bool barrel=layer_>0;

    unsigned int countall=0;
    unsigned int countpass=0;
    
    constexpr unsigned int kNBitsBuffer=3;  

    int writeindex=0;
    int readindex=0;
    std::pair<int,int> projbuffer[1<<kNBitsBuffer]; //iproj zbin

    //The next projection to read, the number of projections and flag if we have more projections to read
    int iproj=0;
    int nproj=vmprojs_->nTracklets();
    bool moreproj=iproj<nproj;

    //Projection that is read from the buffer and compared to stubs  
    int rzbin=0;
    int projfinerz=0;
    int projfinerzadj=0;
    
    int projindex=0;
    int projrinv=0;
    bool isPSseed=false;
    
    //Number of stubs for current zbin and the stub being processed on this clock
    int nstubs=0;
    int istub=0;

    //Main processing loops starts here  
    for (unsigned int istep=0;istep<settings_->maxStep("ME");istep++) {

      countall++;
      
      int writeindexplus=(writeindex+1)%(1<<kNBitsBuffer);
      int writeindexplusplus=(writeindex+2)%(1<<kNBitsBuffer);

      //Determine if buffer is full - or near full as a projection
      //can point to two z bins we might fill two slots in the buffer
      bool bufferfull=(writeindexplus==readindex)||(writeindexplusplus==readindex);

      //Determin if buffer is empty
      bool buffernotempty=(writeindex!=readindex);

      //If we have more projections and the buffer is not full we read
      //next projection and put in buffer if there are stubs in the 
      //memory the projection points to

      if ((!moreproj)&&(!buffernotempty)) break;
		
      if (moreproj&&(!bufferfull)){
	Tracklet* proj=vmprojs_->getFPGATracklet(iproj);

	int iprojtmp=iproj;
	
	iproj++;
	moreproj=iproj<nproj;

	unsigned int rzfirst = barrel?proj->zbin1projvm(layer_):proj->rbin1projvm(disk_);
	unsigned int rzlast = rzfirst;
	bool second=(barrel?proj->zbin2projvm(layer_):proj->rbin2projvm(disk_))==1;
	if (second) rzlast += 1;

	//Check if there are stubs in the memory
	int nstubfirst=vmstubs_->nStubsBin(rzfirst);
	int nstublast=vmstubs_->nStubsBin(rzlast);
	bool savefirst=nstubfirst!=0;
	bool savelast=second&&(nstublast!=0);

	int writeindextmp=writeindex;
	int writeindextmpplus=(writeindex+1)%(1<<kNBitsBuffer);
	
	if (savefirst&&savelast) {
	  writeindex=writeindexplusplus;
	} else if (savefirst||savelast) {
	  writeindex=writeindexplus;
	}

	if (savefirst) { //FIXME code needs to be cleaner
	  std::pair<int,int> tmp(iprojtmp,rzfirst);
	  projbuffer[writeindextmp]=tmp;
	}
	if (savelast) {
	  std::pair<int,int> tmp(iprojtmp,rzlast+100); //FIXME hack to flag that this is second bin
	  if (savefirst) {
	    projbuffer[writeindextmpplus]=tmp;
	  } else {
	    projbuffer[writeindextmp]=tmp;
	  }
	}
      }
      

      //If the buffer is not empty we have a projection that we need to process.

      if (buffernotempty) {

	int istubtmp=istub;

	//New projection
	if (istub==0) {

	  projindex=projbuffer[readindex].first;
	  rzbin=projbuffer[readindex].second;
	  bool second=false;
	  if (rzbin>=100) {
	    rzbin-=100;
	    second=true;
	  }

	  Tracklet* proj=vmprojs_->getFPGATracklet(projindex);

	  nstubs=vmstubs_->nStubsBin(rzbin);

	  projfinerz = barrel?proj->finezvm(layer_):proj->finervm(disk_);

	  projrinv=barrel?(16+(((-2)*proj->fpgaphiprojder(layer_).value())>>(proj->fpgaphiprojder(layer_).nbits()-4))):proj->getBendIndex(disk_).value();
	  assert(projrinv>=0);
	  if (settings_->extended()&&projrinv==32){ //FIXME - should not need this
	    edm::LogPrint("Tracklet") << "Warning: projrinv:"<<projrinv;
	    projrinv=31;
	  }
	  assert(projrinv<32);

	  isPSseed=proj->PSseed()==1;
	  
	  //Calculate fine z position
	  if (second) {
	    projfinerzadj=projfinerz-8;
	  } else {
	    projfinerzadj=projfinerz;
	  }
	  if (nstubs==1) {
	    istub=0;
	    readindex=(readindex+1)%(1<<kNBitsBuffer);
	  } else {
	    istub++;
	  }
	} else {
	  //Check if last stub, if so, go to next buffer entry 
	  if (istub+1>=nstubs){
	    istub=0;
	    readindex=(readindex+1)%(1<<kNBitsBuffer);
	  } else {
	    istub++;
	  }
	}

	//Read vmstub memory and extract data fields
	//std::pair<Stub*,L1TStub*> stub=vmstubs_->getStubBin(rzbin,istubtmp);
	VMStubME vmstub=vmstubs_->getVMStubMEBin(rzbin,istubtmp);

	bool isPSmodule=vmstub.isPSmodule();
	
	int stubfinerz=vmstub.finerz().value();
	       
	int nbits=isPSmodule?3:4;

	//FIXME - not using finephi
	
	unsigned int index=(projrinv<<nbits)+vmstub.bend().value();

	//Check if stub z position consistent
	int idrz=stubfinerz-projfinerzadj;
	bool pass;
	
	if (barrel) {
	  if (isPSseed) {
	    pass=idrz>=-2&&idrz<=2;
	  } else {
	    pass=idrz>=-5&&idrz<=5;
	  }
	} else {
	  if (isPSmodule) {
	    pass=idrz>=-1&&idrz<=1;
	  } else {
	    pass=idrz>=-5&&idrz<=5;
	  }
	}

	//Check if stub bend and proj rinv consistent
	if (pass){
	  if (barrel?table_[index]:(isPSmodule?tablePS_[index]:table2S_[index])) {
	    Tracklet* proj=vmprojs_->getFPGATracklet(projindex);
	    std::pair<Tracklet*,int> tmp(proj,vmprojs_->getAllProjIndex(projindex));
            if (settings_->writeMonitorData("Seeds")) {
              ofstream fout("seeds.txt", ofstream::app);
              fout << __FILE__ << ":" << __LINE__ << " " << name_ << "_" << iSector_ << " " << proj->getISeed() << endl;
              fout.close();
            }
	    candmatches_->addMatch(tmp,vmstub.stub());
	    countpass++;
	  }
	}
      }
      
    }

    if (settings_->writeMonitorData("ME")) {
      globals_->ofstream("matchengine.txt")  << getName()<<" "<<countall<<" "<<countpass<<endl;
    }

}
 

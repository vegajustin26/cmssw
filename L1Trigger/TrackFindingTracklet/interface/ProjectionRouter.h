//This class implementes the projection router
#ifndef L1Trigger_TrackFindingTracklet_interface_ProjectionRouter_h
#define L1Trigger_TrackFindingTracklet_interface_ProjectionRouter_h

#include "ProcessBase.h"
#include "ProjectionRouterBendTable.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;

class ProjectionRouter:public ProcessBase{

public:

 ProjectionRouter(string name, const Settings* settings, Globals* global, unsigned int iSector):
   ProcessBase(name,settings,global,iSector){

    layerdisk_=initLayerDisk(3);

    vmprojs_.resize(settings_->nvmme(layerdisk_),0);

    nrbits_=5;
    nphiderbits_=6;
  }

  void addOutput(MemoryBase* memory,string output){
    if (settings_->writetrace()) {
      edm::LogVerbatim("Tracklet") << "In "<<name_<<" adding output to "<<memory->getName() << " to output "<<output;
    }
    if (output=="allprojout"){
      AllProjectionsMemory* tmp=dynamic_cast<AllProjectionsMemory*>(memory);
      assert(tmp!=0);
      allproj_=tmp;
      return;
    }
    
    unsigned int nproj=settings_->nallstubs(layerdisk_);
    unsigned int nprojvm=settings_->nvmme(layerdisk_);
    
    for (unsigned int iproj=0;iproj<nproj;iproj++) {
      for (unsigned int iprojvm=0;iprojvm<nprojvm;iprojvm++) {
	ostringstream oss;
	oss << "vmprojoutPHI"<<char(iproj+'A')<<iproj*nprojvm+iprojvm+1;
	string name=oss.str();
	if (output==name) {
	  VMProjectionsMemory* tmp=dynamic_cast<VMProjectionsMemory*>(memory);
	  assert(tmp!=0);
	  vmprojs_[iprojvm]=tmp;
	  return;
	}
      }
    }
    
    edm::LogPrint("Tracklet") << getName()<<" Did not find output : "<<output;
    assert(0);
  }

  void addInput(MemoryBase* memory,string input){
    if (settings_->writetrace()) {
      edm::LogVerbatim("Tracklet") << "In "<<name_<<" adding input from "<<memory->getName() << " to input "<<input;
    }
    if (input.substr(0,4)=="proj" && input.substr(input.size()-2,2)=="in"){
      TrackletProjectionsMemory* tmp=dynamic_cast<TrackletProjectionsMemory*>(memory);
      assert(tmp!=0);
      inputproj_.push_back(tmp);
      return;
    }
    edm::LogPrint("Tracklet") << "Could not find input : "<<input<<" in "<<getName();
    assert(0);
  }

  void execute() {


    if (globals_->projectionRouterBendTable()==0){
      ProjectionRouterBendTable* bendTablePtr=new ProjectionRouterBendTable();
      bendTablePtr->init(settings_,globals_, nrbits_, nphiderbits_);
      globals_->projectionRouterBendTable()=bendTablePtr;
    }

    unsigned int allprojcount=0;

    //These are just here to test that the order is correct. Does not affect
    //the actual execution
    
    int lastTCID=-1;
    
    for (unsigned int j=0;j<inputproj_.size();j++){

      for (unsigned int i=0;i<inputproj_[j]->nTracklets();i++){

	if (allprojcount>settings_->maxStep("PR")) continue;

	Tracklet* tracklet=inputproj_[j]->getFPGATracklet(i);

	FPGAWord fpgaphi;

	if (layerdisk_<6) {
	  fpgaphi=tracklet->fpgaphiproj(layerdisk_+1);
	} else {
	  int disk=layerdisk_-5;
	  fpgaphi=tracklet->fpgaphiprojdisk(disk);

	  //The next lines looks up the predicted bend based on:
	  // 1 - r projections
	  // 2 - phi derivative
	  // 3 - the sign - i.e. if track is forward or backward
	  int rindex=(tracklet->fpgarprojdisk(disk).value()>>(tracklet->fpgarprojdisk(disk).nbits()-nrbits_))&((1<<nrbits_)-1);
	  
	  int phiderindex=(tracklet->fpgaphiprojderdisk(disk).value()>>(tracklet->fpgaphiprojderdisk(disk).nbits()-nphiderbits_))&((1<<nphiderbits_)-1);
	  
	  int signindex=(tracklet->fpgarprojderdisk(disk).value()<0);
	  
	  int bendindex=(signindex<<(nphiderbits_+nrbits_))+(rindex<<(nphiderbits_))+phiderindex;
	  
	  int ibendproj=globals_->projectionRouterBendTable()->bendLoookup(disk-1,bendindex);

	  tracklet->setBendIndex(ibendproj,disk);
	}
	
	unsigned int iphivm=fpgaphi.bits(fpgaphi.nbits()-settings_->nbitsallstubs(layerdisk_)-settings_->nbitsvmme(layerdisk_),settings_->nbitsvmme(layerdisk_));
	  
	//This block of code just checks that the configuration is consistent
	if (lastTCID>=tracklet->TCID()) {
	  edm::LogPrint("Tracklet") << "Wrong TCID ordering for projections in "<<getName();
	} else {
	  lastTCID=tracklet->TCID();
	}
	
	allproj_->addTracklet(tracklet);

	vmprojs_[iphivm]->addTracklet(tracklet,allprojcount);

	allprojcount++;
	
      }
    }
    
    
    if (settings_->writeMonitorData("AP")) { 
      globals_->ofstream("allprojections.txt")  << getName() << " " << allproj_->nTracklets() << endl;
    } 
   

    if (settings_->writeMonitorData("VMP")) {
      ofstream& out=globals_->ofstream("chisq.txt"); 
      for (unsigned int i=0;i<8;i++) {
	if (vmprojs_[i]!=0) {
	  out << vmprojs_[i]->getName() << " " << vmprojs_[i]->nTracklets() << endl;
	}
      }
    }
  }

  
private:

  unsigned int layerdisk_;
  
  int nrbits_;
  int nphiderbits_;
  
  vector<TrackletProjectionsMemory*> inputproj_;

  AllProjectionsMemory* allproj_;
  std::vector<VMProjectionsMemory*> vmprojs_;

};

#endif

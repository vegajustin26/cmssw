//This class implementes the projection router
#ifndef L1Trigger_TrackFindingTracklet_interface_ProjectionRouter_h
#define L1Trigger_TrackFindingTracklet_interface_ProjectionRouter_h

#include "ProcessBase.h"

using namespace std;

class ProjectionRouter:public ProcessBase{

public:

 ProjectionRouter(string name, const Settings* settings, unsigned int iSector):
  ProcessBase(name,settings,iSector){

    layerdisk_=initLayerDisk(3);

    vmprojs_.resize(settings_->nvmme(layerdisk_),0);

    nrbits_=5;
    nphiderbits_=6;
  }

  void addOutput(MemoryBase* memory,string output){
    if (settings_->writetrace()) {
      cout << "In "<<name_<<" adding output to "<<memory->getName() << " to output "<<output<<endl;
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
    
    cout << getName()<<" Did not find output : "<<output<<endl;
    assert(0);
  }

  void addInput(MemoryBase* memory,string input){
    if (settings_->writetrace()) {
      cout << "In "<<name_<<" adding input from "<<memory->getName() << " to input "<<input<<endl;
    }
    if (input.substr(0,4)=="proj" && input.substr(input.size()-2,2)=="in"){
      TrackletProjectionsMemory* tmp=dynamic_cast<TrackletProjectionsMemory*>(memory);
      assert(tmp!=0);
      inputproj_.push_back(tmp);
      return;
    }
    cout << "Could not find input : "<<input<<" in "<<getName()<<endl;
    assert(0);
  }

  void execute() {


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
	  
	  int ibendproj=bendTable(disk-1,bendindex);

	  tracklet->setBendIndex(ibendproj,disk);
	}
	
	unsigned int iphivm=fpgaphi.bits(fpgaphi.nbits()-settings_->nbitsallstubs(layerdisk_)-settings_->nbitsvmme(layerdisk_),settings_->nbitsvmme(layerdisk_));
	  
	//This block of code just checks that the configuration is consistent
	if (lastTCID>=tracklet->TCID()) {
	  cout << "Wrong TCID ordering for projections in "<<getName()<<endl;
	} else {
	  lastTCID=tracklet->TCID();
	}
	
	allproj_->addTracklet(tracklet);

	vmprojs_[iphivm]->addTracklet(tracklet,allprojcount);

	allprojcount++;
	
      }
    }
    
    
    if (settings_->writeMonitorData("AP")) {
      static ofstream out("allprojections.txt"); 
      out << getName() << " " << allproj_->nTracklets() << endl;
    } 
   

    if (settings_->writeMonitorData("VMP")) {
      static ofstream out("vmprojections.txt");
      for (unsigned int i=0;i<8;i++) {
	if (vmprojs_[i]!=0) {
	  out << vmprojs_[i]->getName() << " " << vmprojs_[i]->nTracklets() << endl;
	}
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
	    ir=ir<<(nrbitsprojdisk-nrbits_);
	    for(unsigned int iphiderbin=0;iphiderbin<nphiderbins;iphiderbin++) {
	      int iphider=iphiderbin;
	      if (iphider>(1<<(nphiderbits_-1))) iphider-=(1<<nphiderbits_);
	      iphider=iphider<<(nbitsphiprojderL123-nphiderbits_);
	      
	      double rproj=ir*krprojshiftdisk;
	      double phider=iphider*settings_->ITC_L1L2()->der_phiD_final.get_K();
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

  unsigned int layerdisk_;
  
  int nrbits_;
  int nphiderbits_;
  
  vector<TrackletProjectionsMemory*> inputproj_;

  AllProjectionsMemory* allproj_;
  std::vector<VMProjectionsMemory*> vmprojs_;

};

#endif

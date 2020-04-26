//This class implements the tracklet engine
#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletEngine_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletEngine_h

#include "ProcessBase.h"
#include "Util.h"


using namespace std;

class TrackletEngine:public ProcessBase{

public:

 TrackletEngine(string name, const Settings* const settings,unsigned int iSector):
  ProcessBase(name,settings,iSector){

    stubpairs_=0;
    innervmstubs_=0;
    outervmstubs_=0;

    initLayerDisksandISeed(layerdisk1_,layerdisk2_,iSeed_);

    innerphibits_=settings->nfinephi(0,iSeed_);
    outerphibits_=settings->nfinephi(1,iSeed_);
  }

  void addOutput(MemoryBase* memory,string output){
    if (settings_->writetrace()) {
      cout << "In "<<name_<<" adding output to "<<memory->getName() << " to output "<<output<<endl;
    }
    if (output=="stubpairout") {
      StubPairsMemory* tmp=dynamic_cast<StubPairsMemory*>(memory);
      assert(tmp!=0);
      stubpairs_=tmp;
      return;
    }
    assert(0);
  }

  void addInput(MemoryBase* memory,string input){
    if (settings_->writetrace()) {
      cout << "In "<<name_<<" adding input from "<<memory->getName() << " to input "<<input<<endl;
    }
    if (input=="innervmstubin") {
      VMStubsTEMemory* tmp=dynamic_cast<VMStubsTEMemory*>(memory);
      assert(tmp!=0);
      innervmstubs_=tmp;
      setVMPhiBin();
      return;
    }
    if (input=="outervmstubin") {
      VMStubsTEMemory* tmp=dynamic_cast<VMStubsTEMemory*>(memory);
      assert(tmp!=0);
      outervmstubs_=tmp;
      setVMPhiBin();
      return;
    }
    cout << "Could not find input : "<<input<<endl;
    assert(0);
  }

  void execute() {

    if (!settings_->useSeed(iSeed_)) return;
    
    unsigned int countall=0;
    unsigned int countpass=0;

    assert(innervmstubs_!=0);
    assert(outervmstubs_!=0);

    for(unsigned int i=0;i<innervmstubs_->nVMStubs();i++){

      VMStubTE innervmstub=innervmstubs_->getVMStubTE(i);
      FPGAWord lookupbits=innervmstub.vmbits();

      unsigned int nbits=7;
      if (iSeed_==4 || iSeed_==5) nbits=6;
      int rzdiffmax=lookupbits.bits(nbits,lookupbits.nbits()-nbits);	
      int rzbinfirst=lookupbits.bits(0,3);
      int start=lookupbits.bits(4,nbits-4);
      int next=lookupbits.bits(3,1);

      if ((iSeed_==4 || iSeed_==5)&&innervmstub.stub().first->disk().value()<0) { //FIXME - need to store negative disk
	start+=4;
      }
      int last=start+next;

      for(int ibin=start;ibin<=last;ibin++) {
	
	for(unsigned int j=0;j<outervmstubs_->nVMStubsBinned(ibin);j++){
	  if (countall>=settings_->maxStep("TE")) break;
	  countall++;
	  VMStubTE outervmstub=outervmstubs_->getVMStubTEBinned(ibin,j);
	  
	  int rzbin=outervmstub.vmbits().bits(0,3);
	  
	  FPGAWord iphiinnerbin=innervmstub.finephi();
	  FPGAWord iphiouterbin=outervmstub.finephi();
	  
	  unsigned int index = (iphiinnerbin.value()<<outerphibits_)+iphiouterbin.value();
	  
	  if (iSeed_>=4){ //Also use r-position
	    int ir=((ibin&3)<<1)+(rzbin>>2);
	    index=(index<<3)+ir;
	  }
	  
	  if (start!=ibin) rzbin+=8;
	  if ((rzbin<rzbinfirst)||(rzbin-rzbinfirst>rzdiffmax)) {
	    continue;
	  }
	  
	  FPGAWord innerbend=innervmstub.bend();
	  FPGAWord outerbend=outervmstub.bend();
          
	  int ptinnerindex=(index<<innerbend.nbits())+innerbend.value();
	  int ptouterindex=(index<<outerbend.nbits())+outerbend.value();
	  
	  if (!(pttableinner_[ptinnerindex]&&pttableouter_[ptouterindex])) {
	    if (settings_->debugTracklet()) {
	      cout << "Stub pair rejected because of stub pt cut bends : "
		   <<Stub::benddecode(innervmstub.bend().value(),innervmstub.isPSmodule())
		   <<" "
		   <<Stub::benddecode(outervmstub.bend().value(),outervmstub.isPSmodule())
		   <<endl;
	    }		
	    continue;
	  }
	  
	  if (settings_->debugTracklet()) cout << "Adding stub pair in " <<getName()<<endl;

	  stubpairs_->addStubPair(innervmstub,outervmstub);
	  countpass++;
	}
      }
    }
	  
    if (settings_->writeMonitorData("TE")) {
      static ofstream out("trackletengine.txt");
      out << getName()<<" "<<countall<<" "<<countpass<<endl;
    }
    
  }

  void setVMPhiBin() {
    if (innervmstubs_==0 || outervmstubs_==0 ) return;

    innervmstubs_->setother(outervmstubs_);
    outervmstubs_->setother(innervmstubs_);

    int outerrbits=3;
    if (iSeed_<4) {
      outerrbits=0;
    }
    
    int outerrbins=(1<<outerrbits);
    int innerphibins=(1<<innerphibits_);
    int outerphibins=(1<<outerphibits_);
    
    double innerphimin, innerphimax;
    innervmstubs_->getPhiRange(innerphimin,innerphimax,iSeed_,0);

    double outerphimin, outerphimax;
    outervmstubs_->getPhiRange(outerphimin,outerphimax,iSeed_,1);

    double phiinner[2];
    double phiouter[2];
    double router[2];

    unsigned int nbendbitsinner=3;
    unsigned int nbendbitsouter=3;
    if (iSeed_==2) {
      nbendbitsouter=4;
    }
    if (iSeed_==3) {
      nbendbitsinner=4;
      nbendbitsouter=4;
    }
      
    std::vector<bool> vmbendinner((1<<nbendbitsinner),false);
    std::vector<bool> vmbendouter((1<<nbendbitsouter),false);
    
    for (int iphiinnerbin=0;iphiinnerbin<innerphibins;iphiinnerbin++){
      phiinner[0]=innerphimin+iphiinnerbin*(innerphimax-innerphimin)/innerphibins;
      phiinner[1]=innerphimin+(iphiinnerbin+1)*(innerphimax-innerphimin)/innerphibins;
      for (int iphiouterbin=0;iphiouterbin<outerphibins;iphiouterbin++){
	phiouter[0]=outerphimin+iphiouterbin*(outerphimax-outerphimin)/outerphibins;
	phiouter[1]=outerphimin+(iphiouterbin+1)*(outerphimax-outerphimin)/outerphibins;
	for (int irouterbin=0;irouterbin<outerrbins;irouterbin++){
	  if (iSeed_>=4) {
	    router[0]=rmindiskvm+irouterbin*(rmaxdiskvm-rmindiskvm)/outerrbins;
	    router[1]=rmindiskvm+(irouterbin+1)*(rmaxdiskvm-rmindiskvm)/outerrbins;
	  } else {
	    router[0]=settings_->rmean(layerdisk2_);
	    router[1]=settings_->rmean(layerdisk2_);
	  }
	  
	  double bendinnermin=20.0;
	  double bendinnermax=-20.0;
	  double bendoutermin=20.0;
	  double bendoutermax=-20.0;
	  double rinvmin=1.0; 
	  for(int i1=0;i1<2;i1++) {
	    for(int i2=0;i2<2;i2++) {
	      for(int i3=0;i3<2;i3++) {
		double rinner=0.0;
		if (iSeed_==4 || iSeed_==5) {
		  rinner=router[i3]*settings_->zmean(layerdisk1_-6)/settings_->zmean(layerdisk2_-6);
		} else {
		  rinner=settings_->rmean(layerdisk1_);
		}
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

	  bool passptcut=rinvmin<rinvcutte;
	  
	  for(int ibend=0;ibend<(1<<nbendbitsinner);ibend++) {
	    double bend=Stub::benddecode(ibend,nbendbitsinner==3); 
	    
	    bool passinner=bend-bendinnermin>-settings_->bendcutte(0,iSeed_)&&bend-bendinnermax<settings_->bendcutte(0,iSeed_); 
	    if (passinner) vmbendinner[ibend]=true;
	    pttableinner_.push_back(passinner&&passptcut);
	    
	  }
	    
	  for(int ibend=0;ibend<(1<<nbendbitsouter);ibend++) {
	    double bend=Stub::benddecode(ibend,nbendbitsouter==3); 
	    
	    bool passouter=bend-bendoutermin>-settings_->bendcutte(1,iSeed_)&&bend-bendoutermax<settings_->bendcutte(1,iSeed_); 
	    if (passouter) vmbendouter[ibend]=true;
	    pttableouter_.push_back(passouter&&passptcut);
	    
	  }
	  
	}
      }
    }
    
    innervmstubs_->setbendtable(settings_,vmbendinner);
    outervmstubs_->setbendtable(settings_,vmbendouter);
    
    if (iSector_==0&&settings_->writeTable()) writeTETable();
      
  }

  double rinv(double phi1, double phi2,double r1, double r2){

    if (r2<=r1) { //can not form tracklet
      return 20.0; 
    }

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

    ofstream outstubptinnercut;
    outstubptinnercut.open(getName()+"_stubptinnercut.tab");
    outstubptinnercut << "{"<<endl;
    for(unsigned int i=0;i<pttableinner_.size();i++){
      if (i!=0) outstubptinnercut<<","<<endl;
      outstubptinnercut << pttableinner_[i];
    }
    outstubptinnercut <<endl<<"};"<<endl;
    outstubptinnercut.close();
    
    ofstream outstubptoutercut;
    outstubptoutercut.open(getName()+"_stubptoutercut.tab");
    outstubptoutercut << "{"<<endl;
    for(unsigned int i=0;i<pttableouter_.size();i++){
      if (i!=0) outstubptoutercut<<","<<endl;
      outstubptoutercut << pttableouter_[i];
    }
    outstubptoutercut <<endl<<"};"<<endl;
    outstubptoutercut.close();
  }
  
private:

  //Which seed type and which layer/disk is used
  unsigned int iSeed_;
  unsigned int layerdisk1_;  //inner seeding layer
  unsigned int layerdisk2_;  //outer seeding layer

  //The input vmstubs memories
  VMStubsTEMemory* innervmstubs_;
  VMStubsTEMemory* outervmstubs_;

  //The output stub pair memory
  StubPairsMemory* stubpairs_;

  //The stub pt (bend) lookup table for the inner and outer stub
  vector<bool> pttableinner_;
  vector<bool> pttableouter_;

  //Number of phi bits used in the lookup table
  unsigned int innerphibits_;
  unsigned int outerphibits_;
  
};

#endif

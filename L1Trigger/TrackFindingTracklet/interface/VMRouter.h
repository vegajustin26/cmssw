//This class implementes the VM router
#ifndef VMROUTER_H
#define VMROUTER_H

#include "ProcessBase.h"
#include "VMStubTE.h"
#include "TETableOuter.h"
#include "TETableInner.h"
#include "TETableOuterDisk.h"
#include "TETableInnerDisk.h"
#include "TETableInnerOverlap.h"

using namespace std;

class VMRouter:public ProcessBase{

public:

 VMRouter(string name, const Settings* settings, unsigned int iSector):
  ProcessBase(name,settings,iSector){
      
    layerdisk_=initLayerDisk(4);
    initFineBinTable();

    vmstubsMEPHI_.resize(settings_->nvmme(layerdisk_),0);

    overlapbits_=7;
    nextrabits_=overlapbits_-(settings_->nbitsallstubs(layerdisk_)+settings_->nbitsvmme(layerdisk_));
    
  }
   
  void addOutput(MemoryBase* memory,string output){

    if (writetrace) {
      cout << "In "<<name_<<" adding output to "<<memory->getName()
	   << " to output "<<output<<endl;
    }

    if (output.substr(0,10)=="allstubout"){
      AllStubsMemory* tmp=dynamic_cast<AllStubsMemory*>(memory);
      assert(tmp!=0);
      allstubs_.push_back(tmp);
      return;
    }

    if (output.substr(0,12)=="vmstuboutPHI"){
      char seedtype=memory->getName().substr(11,1)[0];
      unsigned int pos=12;
      int vmbin=memory->getName().substr(pos,1)[0]-'0';
      pos++;
      if (pos<memory->getName().size()) {
	if (memory->getName().substr(pos,1)[0]!='n') {
	  vmbin=vmbin*10+memory->getName().substr(pos,1)[0]-'0';
	  pos++;
	}
      }

      int iseed=-1;
      unsigned int inner=1;
      if (memory->getName().substr(3,2)=="TE") {
	VMStubsTEMemory* tmp=dynamic_cast<VMStubsTEMemory*>(memory);
	assert(tmp!=0);
	if (seedtype<'I') {
	  if (layerdisk_==0||layerdisk_==1) iseed=0; 
	  if (layerdisk_==2||layerdisk_==3) iseed=2; 
	  if (layerdisk_==4||layerdisk_==5) iseed=3; 
	  if (layerdisk_==6||layerdisk_==7) iseed=4; 
	  if (layerdisk_==8||layerdisk_==9) iseed=5;
	  if (layerdisk_==0||layerdisk_==2||layerdisk_==4||layerdisk_==6||layerdisk_==8) inner=0;
	} else if (seedtype<'M') {
	  if (layerdisk_==1||layerdisk_==2) iseed=1; 
	  if (layerdisk_==1) inner=0;
	} else if (seedtype<='Z') {
	  if (layerdisk_==0||layerdisk_==6) iseed=6; 
	  if (layerdisk_==1||layerdisk_==6) iseed=7; 
	  if (layerdisk_==0||layerdisk_==1) inner=0;
	} else if (seedtype<'o' && seedtype>='a') {
	  if (layerdisk_==1||layerdisk_==2) iseed=10;
	  if (layerdisk_==1) inner=0;
	} else if (seedtype>'o' && seedtype<='z'){
	  if (layerdisk_==1) iseed=11;
	  if (layerdisk_==6) iseed=10;
	  inner=2; 
	} else {
	  assert(0);
	}
	assert(iseed!=-1);
	int seedindex=-1;
	for(unsigned int k=0;k<vmstubsTEPHI_.size();k++){
	  if (vmstubsTEPHI_[k].first.first==(unsigned int)iseed) {
	    seedindex=k;
	  }
	}
	if (seedindex==-1) {
	  seedindex=vmstubsTEPHI_.size();
	  vector<VMStubsTEMemory*> avectmp;
	  vector< vector<VMStubsTEMemory*> > vectmp(settings_->nvmte(inner,iseed),avectmp);
	  pair<unsigned int, unsigned int > tmppair(iseed,inner);
	  std::pair< std::pair<unsigned int, unsigned int>, vector< vector< VMStubsTEMemory* > > > atmp(tmppair,vectmp);
	  vmstubsTEPHI_.push_back(atmp);
	}
	vmstubsTEPHI_[seedindex].second[(vmbin-1)&(settings_->nvmte(inner,iseed)-1)].push_back(tmp);
	
      } else if (memory->getName().substr(3,2)=="ME") {
	VMStubsMEMemory* tmp=dynamic_cast<VMStubsMEMemory*>(memory);
	assert(tmp!=0);
	vmstubsMEPHI_[(vmbin-1)&(settings_->nvmme(layerdisk_)-1)]=tmp;
      } else {
	assert(0);
      }
      
      return;
    }
    
    cout << "Could not find : "<<output<<endl;
    assert(0);
  }

  void addInput(MemoryBase* memory,string input){
    if (writetrace) {
      cout << "In "<<name_<<" adding input from "<<memory->getName()
	   << " to input "<<input<<endl;
    }
    if (input=="stubin"){
      InputLinkMemory* tmp1=dynamic_cast<InputLinkMemory*>(memory);
      assert(tmp1!=0);
      if (tmp1!=0){
	stubinputs_.push_back(tmp1);
      }
      return;
    }
    cout << "Could not find input : "<<input<<endl;
    assert(0);
  }

  void execute() {

    unsigned int allStubCounter=0;

    //Loop over the input stubs
    for(unsigned int j=0;j<stubinputs_.size();j++){
      for(unsigned int i=0;i<stubinputs_[j]->nStubs();i++){
	if (allStubCounter>MAXVMROUTER) continue;
	if (allStubCounter>127) continue;
	std::pair<Stub*,L1TStub*> stub=stubinputs_[j]->getStub(i);

	bool negdisk=(stub.first->disk().value()<0);//FIXME how do we know that we have negative disk

	//use &127 to make sure we fit into the number of bits -
	//though we should have protected against overflows above
	FPGAWord allStubIndex(allStubCounter&127,7,true,__LINE__,__FILE__);
	
	stub.first->setAllStubIndex(allStubCounter);  //FIXME - should not be needed
	stub.second->setAllStubIndex(allStubCounter); //FIXME - should not be needed

	allStubCounter++;
	
	//Fill allstubs memories - in HLS this is the same write to multiple memories
	for (unsigned int l=0;l<allstubs_.size();l++){
	  allstubs_[l]->addStub(stub);
	}

	//Fill all the ME VM memories 

	FPGAWord iphi=stub.first->phicorr();
	unsigned int ivm=iphi.bits(iphi.nbits()-(settings_->nbitsallstubs(layerdisk_)+settings_->nbitsvmme(layerdisk_)),settings_->nbitsvmme(layerdisk_));
	unsigned int extrabits=iphi.bits(iphi.nbits()-overlapbits_,nextrabits_);

	unsigned int ivmPlus=ivm; 

	if (extrabits==((1U<<nextrabits_)-1)&&ivm!=((1U<<settings_->nbitsvmme(layerdisk_))-1)) ivmPlus++;
	unsigned int ivmMinus=ivm; 
	if (extrabits==0&&ivm!=0) ivmMinus--;

	//Calculate the z and r position for the vmstub FIXME should project to nominal layer/disk position
	int index=-1;
	if (layerdisk_>5) {
	  index=stub.first->r().value();
	  if (stub.first->isPSmodule()){
	    index=stub.first->r().value()>>(stub.first->r().nbits()-nbitsfinebintable_);
	  }
	} else {
	  //Take the top nbitsfinebintable_ bits of the z coordinate. The & is to handle the negative
	  //z values.
	  index=(stub.first->z().value()>>(stub.first->z().nbits()-nbitsfinebintable_))&((1<<nbitsfinebintable_)-1);
	}

	int rzfine=finebintable_[index];

	assert(rzfine>=0);

	int vmbin=rzfine>>3;
	if (negdisk) vmbin+=8;
	rzfine=rzfine&7; 
	
	if (layerdisk_>5) {    
	  stub.first->setfiner(rzfine);
	} else {
	  stub.first->setfinez(rzfine);
	}

	VMStubME vmstub(stub,
			stub.first->iphivmFineBins(5,3), //FIXME
			FPGAWord(rzfine,3,true,__LINE__,__FILE__),
			stub.first->bend(),
			allStubIndex);
			
	
	assert(vmstubsMEPHI_[ivmPlus]!=0);
	vmstubsMEPHI_[ivmPlus]->addStub(vmstub,vmbin);

	if (ivmMinus!=ivmPlus) {
	  assert(vmstubsMEPHI_[ivmMinus]!=0);
	  vmstubsMEPHI_[ivmMinus]->addStub(vmstub,vmbin);
	}

	//Fill the TE VM memories

	for (unsigned int i=0;i<vmstubsTEPHI_.size();i++) {

	  unsigned int iseed=vmstubsTEPHI_[i].first.first;
	  unsigned int inner=vmstubsTEPHI_[i].first.second;

	  if ((iseed==4||iseed==5||iseed==6||iseed==7)&&(!stub.first->isPSmodule())) continue;	  

	  FPGAWord binlookup=lookup(iseed,inner,stub.first->z(),stub.first->r(),negdisk,stub.first->isPSmodule());

	  if (binlookup.value()<0) continue;

	  unsigned int ivmte=iphi.bits(iphi.nbits()-(settings_->nbitsallstubs(layerdisk_)+settings_->nbitsvmte(inner,iseed)),
				       settings_->nbitsvmte(inner,iseed));

	  int bin=-1;
	  if (inner!=0) {
	    bin=binlookup.value()/8;
	    unsigned int tmp=binlookup.value()&7; //three bits in outer layers //FIXME
	    binlookup.set(tmp,3,true,__LINE__,__FILE__);
	  }

	  FPGAWord finephi=stub.first->iphivmFineBins(settings_->nphireg(inner,iseed),settings_->nfinephi(inner,iseed));

	  VMStubTE tmpstub(stub,
			   finephi,
			   stub.first->bend(),
			   binlookup, 
			   allStubIndex);

	  unsigned int nmem=vmstubsTEPHI_[i].second[ivmte].size();

	  assert(nmem>0);
	  
	  for (unsigned int l=0;l<nmem;l++){
	    if (settings_->debugTracklet()) {
	      cout << getName()<<" try adding stub to "<<vmstubsTEPHI_[i].second[ivmte][l]->getName()<<" inner="<<inner<<" bin="<<bin<<endl;
	    }
	    if (inner==0) {
	      vmstubsTEPHI_[i].second[ivmte][l]->addVMStub(tmpstub); 
	    } else {
	      vmstubsTEPHI_[i].second[ivmte][l]->addVMStub(tmpstub,bin); 
	    }
	  } 
	}	
      }
    }	  
  }

  void initFineBinTable(){

    nbitsfinebintable_=8;
    unsigned int nbins=1<<nbitsfinebintable_;

    if (layerdisk_<6) {
      finebintable_.resize(nbins,-1);
      
      for(unsigned int i=0;i<nbins;i++) {

	//awkward bit manipulations since the index is from a signed number...
	int index=(i+(1<<(nbitsfinebintable_-1)))&((1<<nbitsfinebintable_)-1);
	
	finebintable_[index]=(i>>(nbitsfinebintable_-6));
	
      }
    } else {
      for(unsigned int i=0;i<nbins;i++) {

	double rstub=kr*(i<<(nrbitsdisk-nbitsfinebintable_));;

	//use tabulated values for the 2S modules
	if (i<10) rstub=(layerdisk_<=7)?rDSSinner[i]:rDSSouter[i];
	  
	if (rstub<rmindiskvm) {
	  finebintable_.push_back(-1);
	} else {	
	  int rfine=64*(rstub-rmindiskvm)/(rmaxdisk-rmindiskvm);
	  finebintable_.push_back(rfine);
	}
      }
    }
    
    if (iSector_==0&&settings_->writeTable()) {

      ofstream outfinebin;
      outfinebin.open(getName()+"_finebin.tab");
      outfinebin << "{"<<endl;
      for(unsigned int i=0;i<finebintable_.size();i++) {
	if (i!=0) outfinebin<<","<<endl;
	outfinebin << finebintable_[i];
      }
      outfinebin <<endl<<"};"<<endl;
      outfinebin.close();
    }

    
  }

  FPGAWord lookup(unsigned int iseed, unsigned int inner,FPGAWord z, FPGAWord r, bool negdisk, bool isPSmodule){

    static bool first=true;
    static TETableBase* LUTs[3][12];
    
    if (first) {
      first=false;

      for (unsigned int i=0;i<3;i++){
	for (unsigned int iseed=0;iseed<12;iseed++) {
	  LUTs[i][iseed]=0;
	}
      }

      LUTs[0][0]=new TETableInner(settings_,1,2,-1,settings_->zbitstab(0,0),settings_->rbitstab(0,0));
      LUTs[0][1]=new TETableInner(settings_,2,3,-1,settings_->zbitstab(0,1),settings_->rbitstab(0,1));
      LUTs[0][2]=new TETableInner(settings_,3,4,2,settings_->zbitstab(0,2),settings_->rbitstab(0,2));
      LUTs[0][3]=new TETableInner(settings_,5,6,4,settings_->zbitstab(0,3),settings_->rbitstab(0,3));
      LUTs[0][4]=new TETableInnerDisk(settings_,1,2,-1,settings_->zbitstab(0,4),settings_->rbitstab(0,4));
      LUTs[0][5]=new TETableInnerDisk(settings_,3,4,-1,settings_->zbitstab(0,5),settings_->rbitstab(0,5));
      LUTs[0][6]=new TETableInnerOverlap(settings_,1,1,settings_->zbitstab(0,6),settings_->rbitstab(0,6));
      LUTs[0][7]=new TETableInnerOverlap(settings_,2,1,settings_->zbitstab(0,7),settings_->rbitstab(0,7));
      LUTs[0][10]=new TETableInner(settings_,2,3,1,settings_->zbitstab(0,10),settings_->rbitstab(0,10),true);

      
      LUTs[1][0]=new TETableOuter(settings_,2,settings_->zbitstab(1,0),settings_->rbitstab(1,0));
      LUTs[1][1]=new TETableOuter(settings_,3,settings_->zbitstab(1,1),settings_->rbitstab(1,1));
      LUTs[1][2]=new TETableOuter(settings_,4,settings_->zbitstab(1,2),settings_->rbitstab(1,2));
      LUTs[1][3]=new TETableOuter(settings_,6,settings_->zbitstab(1,3),settings_->rbitstab(1,3));
      LUTs[1][4]=new TETableOuterDisk(settings_,2,settings_->zbitstab(1,4),settings_->rbitstab(1,4));
      LUTs[1][5]=new TETableOuterDisk(settings_,4,settings_->zbitstab(1,5),settings_->rbitstab(1,5));
      LUTs[1][6]=new TETableOuterDisk(settings_,1,settings_->zbitstab(1,6),settings_->rbitstab(1,6));
      LUTs[1][7]=new TETableOuterDisk(settings_,1,settings_->zbitstab(1,7),settings_->rbitstab(1,7));
      LUTs[1][10]=new TETableOuter(settings_,3,settings_->zbitstab(1,10),settings_->rbitstab(1,10));

      LUTs[2][10]=new TETableOuterDisk(settings_,1,settings_->zbitstab(2,10),settings_->rbitstab(2,10));
      LUTs[2][11]=new TETableOuter(settings_,2,settings_->zbitstab(2,11),settings_->rbitstab(2,11));

    }

    if (iseed==10&&inner==2) {
      //return if radius to small (values <100 corresponds to the 2S modules)
      if (r.value()>100&&r.value()<rmindiskl3overlapvm/kr) return FPGAWord(-1,2,false,__LINE__,__FILE__);
      int bin=0;
      if(!isPSmodule) {
	bin = r.value(); // 0 to 9 //FIXME
	bin = bin >> 2; // 0 to 2
	bin += 1;
      }
      
      return FPGAWord(bin*8,settings_->lutwidthtabextended(inner,iseed),true,__LINE__,__FILE__);
    }
    
    assert(LUTs[inner][iseed]!=0);

    
    unsigned int zbits=settings_->zbitstab(inner,iseed);
    assert(zbits!=0);
    unsigned int rbits=settings_->rbitstab(inner,iseed);
    assert(rbits!=0);
    unsigned int lutwidth=settings_->lutwidthtab(inner,iseed);
    if (settings_->extended()){
      lutwidth=settings_->lutwidthtabextended(inner,iseed);
    }
    assert(lutwidth!=0);
    
    int zbin=(z.value()+(1<<(z.nbits()-1)))>>(z.nbits()-zbits);
    if (negdisk) zbin=7-zbin;//FIXME should not have hardcoded 7 here
    int rbin=(r.value()+(1<<(r.nbits()-1)))>>(r.nbits()-rbits);
    if (layerdisk_>=6) {
      rbin=r.value()>>(r.nbits()-rbits);
    } 

    int lutvalue=LUTs[inner][iseed]->lookup(zbin,rbin);

    if (lutvalue<0) {
      return FPGAWord(lutvalue,2,false,__LINE__,__FILE__);
    } else {
      return FPGAWord(lutvalue,lutwidth,true,__LINE__,__FILE__);
    }
  }
  

private:

  //0-5 are the layers and 6-10 are the disks
  unsigned int layerdisk_;

  //overlapbits_ is the top bits of phicorr used to add or subtract one to see if stub should be added to
  //two VMs. nextrabits_ is the number of bits beyond the bits for the phivm that is used by overlapbits_
  unsigned int overlapbits_;
  unsigned int nextrabits_;

  int nbitsfinebintable_;
  vector<int> finebintable_;

  //The input stub memories
  vector<InputLinkMemory*> stubinputs_;

  //The all stub memories
  vector<AllStubsMemory*> allstubs_;
  
  //The VM stubs memories used by the MEs
  vector<VMStubsMEMemory*> vmstubsMEPHI_;

  //The VM stubs memories used by the TEs
  // vmstubsTEPHI_[i].first.first is the seed number [0,11] for
  // vmstubsTEPHI_[i].first.second is the stub position in the seed 0 is inner 1 is outer 2 is
  //                               third stub in extended tracking
  // vmstubsTEPHI_[i].second[iVM][n] is the VMStubsTEMemory for iVM and the nth copy 
  vector<pair< pair<unsigned int, unsigned int >, vector< vector<VMStubsTEMemory*> > > > vmstubsTEPHI_;

  
};

#endif


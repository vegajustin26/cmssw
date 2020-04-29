// This class holds a list of stubs for an input link.
//This modules 'owns' the pointers to the stubs. All
//subsequent modules that handles stubs uses a pointer
//to the original stored here
 
#ifndef L1Trigger_TrackFindingTracklet_interface_InputLinkMemory_h
#define L1Trigger_TrackFindingTracklet_interface_InputLinkMemory_h

#include "L1TStub.h"
#include "Stub.h"
#include "MemoryBase.h"
#include "VMRouterPhiCorrTable.h"
#include "GlobalHistTruth.h"
#include <cmath>
#include <sstream>
#include <cctype>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;

class InputLinkMemory:public MemoryBase{

public:

 InputLinkMemory(string name, const Settings* const settings, unsigned int iSector, 
		 double, double):
  MemoryBase(name, settings, iSector){
    
    string subname=name.substr(5,7);
    phiregion_=subname[3]-'A';
    assert(phiregion_>=0&&phiregion_<8);
    
    layerdisk_=initLayerDisk(3);
  }

  bool addStub(const Settings* settings, GlobalHistTruth* globals, L1TStub& al1stub, Stub& stub, string dtc="") {
    
    if (layerdisk_<6&&globals->phiCorr(layerdisk_)==0) {
      globals->phiCorr(layerdisk_)=new VMRouterPhiCorrTable();
      int nbits=3;
      if (layerdisk_>=3) nbits=4;
      globals->phiCorr(layerdisk_)->init(settings,layerdisk_+1,nbits,3);
    }




    int stublayerdisk=stub.layer().value();
    if (stublayerdisk==-1) stublayerdisk=5+abs(stub.disk().value());
    assert(stublayerdisk>=0&&stublayerdisk<11);

    if (stublayerdisk!=layerdisk_) return false;
       
    if (layerdisk_<6) {
      FPGAWord r=stub.r();
      int bendbin=stub.bend().value();
      int rbin=(r.value()+(1<<(r.nbits()-1)))>>(r.nbits()-3);
      const VMRouterPhiCorrTable& phiCorrTable=*globals->phiCorr(layerdisk_);
      int iphicorr=phiCorrTable.getphiCorrValue(bendbin,rbin);
      stub.setPhiCorr(iphicorr);
    }
    
    int iphivmRaw=stub.iphivmRaw();
    int phibin=iphivmRaw/(32/settings_->nallstubs(layerdisk_));
    
    if (phibin!=phiregion_) return false;
    
    if (getName().substr(10,dtc.size())!=dtc) return false;
    
    string half=getName().substr(getName().size()-3,3);
    if (half[1]!='n') {
      half=getName().substr(getName().size()-1,1);
    }

    assert(half[0]=='A' || half[0]=='B');

    if (half[0]=='B' && iphivmRaw<=15) return false;
    if (half[0]=='A' && iphivmRaw>15) return false;
    
    if (settings_->debugTracklet()) {
      edm::LogVerbatim("Tracklet") << "Will add stub in "<<getName()<<" "<<"iphiwmRaw = "<<iphivmRaw<<" phi="<<al1stub.phi()<<" z="<<al1stub.z()<<" r="<<al1stub.r();
    }
    if (stubs_.size()<settings_->maxStep("Link")) {
      L1TStub* l1stub=new L1TStub(al1stub);
      Stub* stubptr=new Stub(stub);

      std::pair<Stub*,L1TStub*> tmp(stubptr,l1stub);
      stubs_.push_back(tmp);
    }
    return true;
  }

  unsigned int nStubs() const {return stubs_.size();}

  std::pair<Stub*,L1TStub*> getStub(unsigned int i) const {return stubs_[i];}

  void writeStubs(bool first) {

    openFile(first,"../data/MemPrints/InputStubs/InputStubs_");
    
    for (unsigned int j=0;j<stubs_.size();j++){
      string stub = stubs_[j].first->str();
      if (j<16) out_ <<"0";
      out_ << hex << j << dec;
      out_ << " " << stub <<" "<<Trklet::hexFormat(stub)<< endl;
    }
    out_.close();

  }


  void clean() {
    for(unsigned int i=0;i<stubs_.size();i++){
      delete stubs_[i].first;
      delete stubs_[i].second;
    }
    stubs_.clear();
  }
  
private:

  vector<std::pair<Stub*,L1TStub*> > stubs_;

  int phiregion_;
  
  int layerdisk_; //FIXME should be unsigned 

};

#endif

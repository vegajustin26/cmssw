//Base class for processing modules
#ifndef L1Trigger_TrackFindingTracklet_interface_ProcessBase_h
#define L1Trigger_TrackFindingTracklet_interface_ProcessBase_h

#include "Settings.h"
#include "Globals.h"
#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace Trklet;
using namespace std;

#include "Stub.h"

class ProcessBase{

public:

 ProcessBase(string name, const Settings* const settings, Globals* global, unsigned int iSector):
  settings_(settings), globals_(global) {
    name_=name;
    iSector_=iSector;
  }

  virtual ~ProcessBase() { } 

  // Add wire from pin "output" or "input" this proc module to memory instance "memory".

  virtual void addOutput(MemoryBase* memory,string output)=0;

  virtual void addInput(MemoryBase* memory,string input)=0;

  string getName() const {return name_;}

  unsigned int nbits(unsigned int power) {

    if (power==2) return 1;
    if (power==4) return 2;
    if (power==8) return 3;
    if (power==16) return 4;
    if (power==32) return 5;

    edm::LogVerbatim("Tracklet") << "nbits: power = "<<power;
    assert(0);

    return -1;
    
  }

  //method sets the layer and disk based on the name. pos is the position in the
  //memory name where the layer or disk is specified
  void initLayerDisk(unsigned int pos, int& layer, int& disk){

    string subname=name_.substr(pos,2);
    layer=0;
    disk=0;

    if (subname=="L1") layer=1;
    if (subname=="L2") layer=2;
    if (subname=="L3") layer=3;
    if (subname=="L4") layer=4;
    if (subname=="L5") layer=5;
    if (subname=="L6") layer=6;
    if (subname=="D1") disk=1;
    if (subname=="D2") disk=2;
    if (subname=="D3") disk=3;
    if (subname=="D4") disk=4;
    if (subname=="D5") disk=5;
    if (layer==0&&disk==0) {
      edm::LogPrint("Tracklet") << "Memoryname = "<<name_<<" subname = "<<subname<<" layer "<<layer<<" disk "<<disk;
    }
    assert((layer!=0)||(disk!=0));
  }

  void initLayerDisk(unsigned int pos, int& layer, int& disk, int& layerdisk){

    initLayerDisk(pos,layer,disk);

    layerdisk=layer-1;
    if (disk>0) layerdisk=5+disk;

  }

  unsigned int initLayerDisk(unsigned int pos){

    int layer,disk;
    initLayerDisk(pos,layer,disk);

    if (disk>0) return 5+disk;
    return layer-1;

  }

  
  //This function processes the name of a TE module to determine the layerdisks and
  //iseed
  void initLayerDisksandISeed(unsigned int& layerdisk1,
			      unsigned int& layerdisk2,
			      unsigned int& iSeed) {

    layerdisk1=99;
    layerdisk2=99;
    
    if (name_[3]=='L') {
      layerdisk1=name_[4]-'1';
    }
    if (name_[3]=='D') {
      layerdisk1=6+name_[4]-'1';
    }
    if (name_[11]=='L') {
      layerdisk2=name_[12]-'1';
    }
    if (name_[11]=='D') {
      layerdisk2=6+name_[12]-'1';
    }
    if (name_[12]=='L') {
      layerdisk2=name_[13]-'1';
    }
    if (name_[12]=='D') {
      layerdisk2=6+name_[13]-'1';
    }
    
    if (layerdisk1 == 0 && layerdisk2 == 1) iSeed = 0;
    else if (layerdisk1 == 1 && layerdisk2 == 2) iSeed = 1;
    else if (layerdisk1 == 2 && layerdisk2 == 3) iSeed = 2;
    else if (layerdisk1 == 4 && layerdisk2 == 5) iSeed = 3;
    else if (layerdisk1 == 6 && layerdisk2 == 7) iSeed = 4;
    else if (layerdisk1 == 8 && layerdisk2 == 9) iSeed = 5;
    else if (layerdisk1 == 0 && layerdisk2 == 6) iSeed = 6;
    else if (layerdisk1 == 1 && layerdisk2 == 6) iSeed = 7;
    else {
      assert(0);
    }
    
    

  }
  
  unsigned int getISeed(std::string name){
  
    //assumes here that namme is on the form XX_L1L2_XXX where L1L2 gives iSeed=0
    
    std::size_t pos = name.find("_");    
    std::string name1 = name.substr (pos+1); 

    pos = name1.find("_");    
    std::string name2 = name1.substr (0,pos); 

    if (name2=="L1L2") return 0;
    if (name2=="L2L3") return 1;
    if (name2=="L3L4") return 2;
    if (name2=="L5L6") return 3;
    if (name2=="D1D2") return 4;
    if (name2=="D3D4") return 5;
    if (name2=="L1D1") return 6;
    if (name2=="L2D1") return 7;
    
    if (name2=="L1L2XX") return 0;
    if (name2=="L2L3XX") return 1;
    if (name2=="L3L4XX") return 2;
    if (name2=="L5L6XX") return 3;
    if (name2=="D1D2XX") return 4;
    if (name2=="D3D4XX") return 5;
    if (name2=="L1D1XX") return 6;
    if (name2=="L2D1XX") return 7;
    if (name2=="L3L4L2") return 8;
    if (name2=="L5L6L4") return 9;
    if (name2=="L2L3D1") return 10;
    if (name2=="D1D2L2") return 11;
    
    edm::LogPrint("Tracklet") << getName()<<" name name1 name2 "<<name<<" - "<<name1<<" - "<<name2;
    assert(0);
    return 0;

  }

protected:

  string name_;
  unsigned int iSector_;

  const Settings* const settings_;
  Globals* globals_;
  
};

#endif

#ifndef SLHCEVENT_H
#define SLHCEVENT_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <map>
#include <ext/hash_set>
#include <math.h>
#include <assert.h>
#include "L1TStub.h"

using namespace std;

static double x_offset=0.199196*0.0;
static double y_offset=0.299922*0.0;


class L1SimTrack{

public:

  L1SimTrack() {
    eventid_=-1; 
    trackid_=-1;   
  }

  L1SimTrack(int eventid, int trackid, int type, double pt, double eta, double phi, 
	     double vx, double vy, double vz) {
    eventid_=eventid;
    trackid_=trackid;
    type_=type;
    pt_=pt;
    eta_=eta;
    phi_=phi;
    vx_=vx;
    vy_=vy;
    vz_=vz;
  }

  void write(ofstream& out){
    
    if (pt_ > -2.0) {
    out << "SimTrack: " 
	<< eventid_ << "\t" 
	<< trackid_ << "\t" 
	<< type_ << "\t" 
	<< pt_ << "\t" 
	<< eta_ << "\t" 
	<< phi_ << "\t" 
	<< vx_ << "\t" 
	<< vy_ << "\t" 
	<< vz_ << "\t" << endl; 
    }
	
  }
  void write(ostream& out){
    
    if (pt_ > -2) {
    out << "SimTrack: " 
	<< eventid_ << "\t" 
	<< trackid_ << "\t" 
	<< type_ << "\t" 
	<< pt_ << "\t" 
	<< eta_ << "\t" 
	<< phi_ << "\t" 
	<< vx_ << "\t" 
	<< vy_ << "\t" 
	<< vz_ << "\t" << endl; 
    }

  }
  
  int eventid() const { return eventid_; }
  int trackid() const { return trackid_; }
  int type() const { return type_; }
  double pt() const { return pt_; }
  double rinv() const { return charge()*0.01*0.3*3.8/pt_; }
  double eta() const { return eta_; }
  double phi() const { return phi_; }
  double vx() const { return vx_; }
  double vy() const { return vy_; }
  double vz() const { return vz_; }
  double dxy() const { return -vx() * sin(phi()) + vy() * cos(phi()); }
  double d0() const { return -dxy(); }
  int charge() const {
     if (type_==11) return -1;
     if (type_==13) return -1;
     if (type_==-211) return -1;
     if (type_==-321) return -1;
     if (type_==-2212) return -1;
     return 1;
  }
  
private:

  int eventid_;
  int trackid_;
  int type_;
  double pt_;
  double eta_;
  double phi_;
  double vx_;
  double vy_;
  double vz_;

};



class SLHCEvent{

public:


  SLHCEvent() {
    //empty constructor to be used with 'filler' functions
    eventnum_=0;
  }

  void setIPx(double x) { x_offset=x;}
  void setIPy(double y) { y_offset=y;}

  void setEventNum(int eventnum) { eventnum_=eventnum; }

  void addL1SimTrack(int eventid,int trackid,int type,double pt,double eta,double phi,
		     double vx,double vy,double vz){

    vx-=x_offset;
    vy-=y_offset;
    L1SimTrack simtrack(eventid,trackid,type,pt,eta,phi,vx,vy,vz);
    simtracks_.push_back(simtrack);

  }

  bool addStub(int layer,int ladder,int module, int strip, int eventid, vector<int> tps, 
              double pt,double bend,
              double x,double y,double z,
              vector<bool> innerStack,
              vector<int> irphi,
              vector<int> iz,
              vector<int> iladder,
              vector<int> imodule,
              int isPSmodule,
              int isFlipped){

    
    if (layer>999&&layer<1999&& z<0.0) {
      layer+=1000;
    }
    
    layer--;   
    x-=x_offset;
    y-=y_offset;

    
    L1TStub stub(eventid,tps,-1,-1,layer, ladder, module, strip, 
		 x, y, z, -1.0, -1.0, pt, bend, isPSmodule, isFlipped);

    for(unsigned int i=0;i<innerStack.size();i++){
      if (innerStack[i]) {
	stub.AddInnerDigi(iladder[i],imodule[i],irphi[i],iz[i]);
      }
      else {
	stub.AddOuterDigi(iladder[i],imodule[i],irphi[i],iz[i]);
      }
    }   

    stub.setiphi(stub.diphi());
    stub.setiz(stub.diz());

    stubs_.push_back(stub);
    return true;
  }

  L1TStub lastStub(){
    return stubs_.back();
  }

  SLHCEvent(istream& in) {

    string tmp;
    in >> tmp;
    while (tmp=="Map:") {
      in>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp;
      in>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp;
    }
    if (tmp=="EndMap") {
      in>>tmp;
    }
    if (tmp!="Event:") {
      cout << "Expected to read 'Event:' but found:"<<tmp<<endl;
      if (tmp=="") {
	cout << "WARNING: fewer events to process than specified!" << endl;
	return;
      }
      else {
	cout << "ERROR, aborting reading file" << endl;
	abort();
      }
    }
    in >> eventnum_;


    // read the SimTracks
    in >> tmp;
    while (tmp!="SimTrackEnd"){
      if (!(tmp=="SimTrack:"||tmp=="SimTrackEnd")) {
	cout << "Expected to read 'SimTrack:' or 'SimTrackEnd' but found:"
	     << tmp << endl;
	abort();
      }
      int eventid;
      int trackid;
      int type;
      string pt_str;
      string eta_str;
      string phi_str;
      string vx_str;
      string vy_str;
      string vz_str;
      double pt;
      double eta;
      double phi;
      double vx;
      double vy;
      double vz;
      in >> eventid >> trackid >> type >> pt_str >> eta_str >> phi_str >> vx_str >> vy_str >> vz_str;
      pt = strtod(pt_str.c_str(), NULL);
      eta = strtod(eta_str.c_str(), NULL);
      phi = strtod(phi_str.c_str(), NULL);
      vx = strtod(vx_str.c_str(), NULL);
      vy = strtod(vy_str.c_str(), NULL);
      vz = strtod(vz_str.c_str(), NULL);
      vx-=x_offset;
      vy-=y_offset;
      L1SimTrack simtrack(eventid,trackid,type,pt,eta,phi,vx,vy,vz);
      simtracks_.push_back(simtrack);
      in >> tmp;
    }

    int nlayer[11];
    for (int i=0;i<10;i++) {
      nlayer[i]=0;
    }

    int oldlayer=0;
    int oldladder=0;
    int oldmodule=0;
    int oldcbc=-1;
    int count=1;
    double oldz=-1000.0;
    
    //read stubs
    in >> tmp;
    while (tmp!="StubEnd"){

      if (!in.good()) {
	cout << "File not good"<<endl;
	abort();
      };
      if (!(tmp=="Stub:"||tmp=="StubEnd")) {
	cout << "Expected to read 'Stub:' or 'StubEnd' but found:"
	     << tmp << endl;
	abort();
      }
      int layer;
      int ladder;
      int module;
      int eventid;
      vector<int> tps;
      int strip;
      double pt;
      double x;
      double y;
      double z;
      double bend;
      int isPSmodule;
      int isFlipped;

      unsigned int ntps;

      in >> layer >> ladder >> module >> strip >> eventid >> pt >> x >> y >> z >> bend >> isPSmodule >> isFlipped >> ntps;

      for(unsigned int itps=0;itps<ntps;itps++){
	int tp;
	in >> tp;
	tps.push_back(tp);
      }

      if (layer>999&&layer<1999&& z<0.0) {
	layer+=1000;
      }

      int cbc=strip/126;
      if (layer>3&&layer==oldlayer&&ladder==oldladder&&module==oldmodule&&cbc==oldcbc&&std::abs(oldz-z)<1.0){
	count++;
      } else {
	oldlayer=layer;
	oldladder=ladder;
	oldmodule=module;
	oldcbc=cbc;
	oldz=z;
	count=1;
      }

      layer--;   
      x-=x_offset;
      y-=y_offset;

      if (layer < 10) nlayer[layer]++;

      L1TStub stub(eventid,tps,-1,-1,layer, ladder, module, strip, x, y, z, -1.0, -1.0, pt, bend, isPSmodule, isFlipped);

      in >> tmp;

      while (tmp=="InnerStackDigi:"||tmp=="OuterStackDigi:"){
	int irphi;
	int iz;
        int iladder;
        int imodule;
	in >> irphi;
	in >> iz;
	in >> iladder; 
        in >> imodule;
	if (tmp=="InnerStackDigi:") stub.AddInnerDigi(iladder,imodule,irphi,iz);
	if (tmp=="OuterStackDigi:") stub.AddOuterDigi(iladder,imodule,irphi,iz);
	in >> tmp;
      }   

      double t=std::abs(stub.z())/stub.r();
      double eta=asinh(t);
      
      if (std::abs(eta)<2.6&&count<=100) {
	stubs_.push_back(stub);
      }
      
    }
  }

  void write(ofstream& out){
    
    out << "Event: "<<eventnum_ << endl;
      
    for (unsigned int i=0; i<simtracks_.size(); i++) {
      simtracks_[i].write(out);
    }
    out << "SimTrackEnd" << endl;

    for (unsigned int i=0; i<stubs_.size(); i++) {
      stubs_[i].write(out);
    }
    out << "StubEnd" << endl;
    
  }

  void write(ostream& out){
    
    out << "Event: "<<eventnum_ << endl;
      
    for (unsigned int i=0; i<simtracks_.size(); i++) {
      simtracks_[i].write(out);
    }
    out << "SimTrackEnd" << endl;

    for (unsigned int i=0; i<stubs_.size(); i++) {
      stubs_[i].write(out);
    }
    out << "StubEnd" << endl;
    
  }

  unsigned int layersHit(int tpid, int &nlayers, int &ndisks){

    int l1=0;
    int l2=0;
    int l3=0;
    int l4=0;
    int l5=0;
    int l6=0;

    int d1=0;
    int d2=0;
    int d3=0;
    int d4=0;
    int d5=0;

    for (unsigned int istub=0; istub<stubs_.size(); istub++){
      if (stubs_[istub].tpmatch(tpid)){
	if (stubs_[istub].layer()==0) l1=1;
        if (stubs_[istub].layer()==1) l2=1;
        if (stubs_[istub].layer()==2) l3=1;
        if (stubs_[istub].layer()==3) l4=1;
	if (stubs_[istub].layer()==4) l5=1;
        if (stubs_[istub].layer()==5) l6=1;

        if (abs(stubs_[istub].disk())==1) d1=1;
        if (abs(stubs_[istub].disk())==2) d2=1;
        if (abs(stubs_[istub].disk())==3) d3=1;
        if (abs(stubs_[istub].disk())==4) d4=1;
        if (abs(stubs_[istub].disk())==5) d5=1;
      }

    }

    nlayers=l1+l2+l3+l4+l5+l6;
    ndisks=d1+d2+d3+d4+d5;

    return l1+2*l2+4*l3+8*l4+16*l5+32*l6+64*d1+128*d2+256*d3+512*d4+1024*d5;

  }


  
  int nstubs() { return stubs_.size(); }

  L1TStub stub(int i) { return stubs_[i]; }

  unsigned int nsimtracks() { return simtracks_.size(); }

  L1SimTrack simtrack(int i) { return simtracks_[i]; }

  int eventnum() const { return eventnum_; }

  int getSimtrackFromSimtrackid(int simtrackid, int eventid=0) const {
    for(unsigned int i=0;i<simtracks_.size();i++){
      if (simtracks_[i].trackid()==simtrackid && simtracks_[i].eventid()==eventid) return i;
    }
    return -1;
  }


private:

  int eventnum_;
  vector<L1SimTrack> simtracks_;
  vector<L1TStub> stubs_;


};

#endif




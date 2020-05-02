#ifndef L1Trigger_TrackFindingTracklet_interface_Util_h
#define L1Trigger_TrackFindingTracklet_interface_Util_h

#include "FWCore/MessageLogger/interface/MessageLogger.h"

namespace Trklet{

  //method return phi in the -pi to +pi range
  inline double phiRange(double phi) {
    //catch if phi is very out of range, not a number etc
    assert(std::abs(phi) < 100.0);
    while (phi < -M_PI)
      phi += 2 * M_PI;
    while (phi > M_PI)
      phi -= 2 * M_PI;
    return phi;
  }

  //Converts string in binary to hex (used in writing out memory content)
  inline std::string hexFormat(std::string binary){

    std::string tmp="";
  
    unsigned int value=0;
    unsigned int bits=0;
    
    for(unsigned int i=0;i<binary.size();i++) {
      unsigned int slot=binary.size()-i-1;
      if (!(binary[slot]=='0'||binary[slot]=='1')) continue;
      value+=((binary[slot]-'0')<<bits);
      bits++;
      if (bits==4||i==binary.size()-1) {
	assert(value<16);
	if (value==0) tmp+="0";
	if (value==1) tmp+="1";
	if (value==2) tmp+="2";
	if (value==3) tmp+="3";
	if (value==4) tmp+="4";
	if (value==5) tmp+="5";
	if (value==6) tmp+="6";
	if (value==7) tmp+="7";
	if (value==8) tmp+="8";
	if (value==9) tmp+="9";
	if (value==10) tmp+="A";
	if (value==11) tmp+="B";
	if (value==12) tmp+="C";
	if (value==13) tmp+="D";
	if (value==14) tmp+="E";
	if (value==15) tmp+="F";
	bits=0;
	value=0;
      }
    }

    std::string hexstring="0x";
    for (unsigned int i=0;i<tmp.size();i++){
      hexstring+=tmp[tmp.size()-i-1];
    }
    
    return hexstring;
     
  }

  //Should be optimized by layer - now first implementation to make sure it works OK
  inline int bendencode(double bend, bool isPS) {
  
    int ibend=2.0*bend;
  
    assert(std::abs(ibend-2.0*bend)<0.1);
  
    if (isPS) {
 
      if (ibend==0||ibend==1) return 0;
      if (ibend==2||ibend==3) return 1;
      if (ibend==4||ibend==5) return 2;
      if (ibend>=6) return 3;
      if (ibend==-1||ibend==-2) return 4;
      if (ibend==-3||ibend==-4) return 5;
      if (ibend==-5||ibend==-6) return 6;
      if (ibend<=-7) return 7;
    
      assert(0);
    
    }
  
  
    if (ibend==0||ibend==1) return 0;
    if (ibend==2||ibend==3) return 1;
    if (ibend==4||ibend==5) return 2;
    if (ibend==6||ibend==7) return 3;
    if (ibend==8||ibend==9) return 4;
    if (ibend==10||ibend==11) return 5;
    if (ibend==12||ibend==13) return 6;
    if (ibend>=14) return 7;
    if (ibend==-1||ibend==-2) return 8;
    if (ibend==-3||ibend==-4) return 9;
    if (ibend==-5||ibend==-6) return 10;
    if (ibend==-7||ibend==-8) return 11;
    if (ibend==-9||ibend==-10) return 12;
    if (ibend==-11||ibend==-12) return 13;
    if (ibend==-13||ibend==-14) return 14;
    if (ibend<=-15) return 15;
  
    edm::LogVerbatim("Tracklet") << "bend ibend : "<<bend<<" "<<ibend;
  
    assert(0);
  
  }
      
  //Should be optimized by layer - now first implementation to make sure it works OK
  inline double benddecode(int ibend, bool isPS) {
  
    if (isPS) {
    
      if (ibend==0) return 0.25;
      if (ibend==1) return 1.25;
      if (ibend==2) return 2.25;
      if (ibend==3) return 3.25;
      if (ibend==4) return -0.75;
      if (ibend==5) return -1.75;
      if (ibend==6) return -2.75;
      if (ibend==7) return -3.75;
      
      assert(0);
    }
    
    if (ibend==0) return 0.25;
    if (ibend==1) return 1.25;
    if (ibend==2) return 2.25;
    if (ibend==3) return 3.25;
    if (ibend==4) return 4.25;
    if (ibend==5) return 5.25;
    if (ibend==6) return 6.25;
    if (ibend==7) return 7.25;
    if (ibend==8) return -0.75;
    if (ibend==9) return -1.75;
    if (ibend==10) return -2.75;
    if (ibend==11) return -3.75;
    if (ibend==12) return -4.75;
    if (ibend==13) return -5.75;
    if (ibend==14) return -6.75;
    if (ibend==15) return -7.75;
    
    assert(0);
    
  }

  inline double bend(double r, double rinv) {
    
    double dr=0.18;
    
    double delta=r*dr*0.5*rinv;
    
    double bend=-delta/0.009;
    if (r<55.0) bend=-delta/0.01;
  
    return bend;
  }
  
  inline double rinv(double phi1, double phi2,double r1, double r2) {
    
    if (r2<=r1) { //can not form tracklet
      return 20.0; 
    }
    
    double dphi=phi2-phi1;
    double dr=r2-r1;
    
    return 2.0*sin(dphi)/dr/sqrt(1.0+2*r1*r2*(1.0-cos(dphi))/(dr*dr));
    
  }

  inline int round_int(double r) { return (r > 0.0) ? (r + 0.5) : (r - 0.5); }


  
};

#endif

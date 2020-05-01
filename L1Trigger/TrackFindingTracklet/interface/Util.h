#ifndef L1Trigger_TrackFindingTracklet_interface_Util_h
#define L1Trigger_TrackFindingTracklet_interface_Util_h

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


  
};

#endif

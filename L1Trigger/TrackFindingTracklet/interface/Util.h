#ifndef L1Trigger_TrackFindingTracklet_interface_Util_h
#define L1Trigger_TrackFindingTracklet_interface_Util_h

#include <sstream>
#include <cassert>
#include <cmath>
#include <algorithm>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

namespace trklet {

  //Converts string in binary to hex (used in writing out memory content)
  inline std::string hexFormat(const std::string& binary) {
    std::stringstream ss;

    unsigned int radix = 1, value = 0;
    for (int i = binary.length() - 1; i >= 0; i--) {
      if (binary.at(i) != '0' && binary.at(i) != '1')
        continue;
      value += (binary.at(i) - '0') * radix;
      if (radix == 8) {
        ss << std::hex << value;
        radix = 1;
        value = 0;
      } else
        radix <<= 1;
    }
    if (radix != 1)
      ss << std::hex << value;

    std::string str = ss.str() + "x0";
    std::reverse(str.begin(), str.end());
    return str;
  }

  inline double bendstrip(double r, double rinv, double stripPitch) {
    constexpr double dr = 0.18;
    double delta = r * dr * 0.5 * rinv;
    double bend = delta / stripPitch;
    return bend;
  }

  inline double rinv(double phi1, double phi2, double r1, double r2) {

    if (r2 <= r1) {  //FIXME can not form tracklet should not call function with r2<=r1
      return 20.0;
    }

    double dphi = phi2 - phi1;
    double dr = r2 - r1;

    return 2.0 * sin(dphi) / dr / sqrt(1.0 + 2 * r1 * r2 * (1.0 - cos(dphi)) / (dr * dr));
  }

  inline std::string convertHexToBin(const std::string& stubwordhex) {

    std::string stubwordbin="";

    for(char word:stubwordhex){
      std::string hexword="";
      if (word=='0') hexword="0000";
      if (word=='1') hexword="0001";
      if (word=='2') hexword="0010";
      if (word=='3') hexword="0011";
      if (word=='4') hexword="0100";
      if (word=='5') hexword="0101";
      if (word=='6') hexword="0110";
      if (word=='7') hexword="0111";
      if (word=='8') hexword="1000";
      if (word=='9') hexword="1001";
      if (word=='A') hexword="1010";
      if (word=='B') hexword="1011";
      if (word=='C') hexword="1100";
      if (word=='D') hexword="1101";
      if (word=='E') hexword="1110";
      if (word=='F') hexword="1111";
      if (hexword=="") {
	throw cms::Exception("Inconsistency") << __FILE__ << " " << __LINE__ << " hex string format invalid: " << stubwordhex;
      }
      stubwordbin+=hexword;
    }
    return stubwordbin;
  }


};  // namespace trklet
#endif

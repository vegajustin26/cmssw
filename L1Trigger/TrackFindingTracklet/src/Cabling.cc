#include "L1Trigger/TrackFindingTracklet/interface/Cabling.h"
#include "L1Trigger/TrackFindingTracklet/interface/DTCLink.h"
#include "L1Trigger/TrackFindingTracklet/interface/DTC.h"
#include "L1Trigger/TrackFindingTracklet/interface/Util.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;

Cabling::Cabling() {}

void Cabling::init(string dtcconfig, string moduleconfig) {
  ifstream indtc(dtcconfig.c_str());
  assert(indtc.good());
  
  string dtc;
  int isec;
  
  while (indtc.good()) {
    indtc >> dtc >> isec;
    
    if (!indtc.good())
      continue;
    
    if (dtcs.find(dtc) == dtcs.end()) {
      dtcs[dtc].init(dtc);
    }
    
    dtcs[dtc].addSec(isec);
    
    string dtcbase = dtc.substr(2, dtc.size() - 2);
    if (dtc[0] == 'n') {
      dtcbase = "neg_" + dtc.substr(6, dtc.size() - 6);
    }
    if (dtcranges.find(dtcbase) == dtcranges.end()) {
      dtcranges[dtcbase].init(dtcbase);
    }
  }
  
  ifstream inmodules(moduleconfig.c_str());
  
  int layer, ladder, module;
  
  while (inmodules.good()) {
    inmodules >> layer >> ladder >> module >> dtc;
    if (module > 300) {
      if (layer == 1) {
	module = (module - 300) + 12;
      }
      if (layer == 2) {
	module = (module - 300) + 12;
      }
      if (layer == 3) {
	module = (module - 300) + 12;
      }
      if (layer > 3) {
	module = (module - 300);
      }
      }
    if (module > 200) {
      module = (module - 200);
    }
    if (module > 100) {
      if (layer == 1) {
	module = (module - 100) + 19;
      }
      if (layer == 2) {
	module = (module - 100) + 23;
      }
      if (layer == 3) {
	module = (module - 100) + 27;
      }
    }
    if (!inmodules.good())
      break;
    modules[layer][ladder][module] = dtc;
  }
}

string Cabling::dtc(int layer, int ladder, int module) {
  std::map<int, std::map<int, std::map<int, string> > >::const_iterator it1 = modules.find(layer);
  assert(it1 != modules.end());
  std::map<int, std::map<int, string> >::const_iterator it2 = it1->second.find(ladder);
  assert(it2 != it1->second.end());
  std::map<int, string>::const_iterator it3 = it2->second.find(module);
  if (it3 == it2->second.end()) {
    edm::LogPrint("Tracklet") << "Could not add stub " << layer << " " << ladder << " " << module;
    assert(0);
  }
  string dtc = it3->second;
  return dtc;
}

void Cabling::addphi(string dtc, double phi, int layer, int module) {
  int layerdisk = layer - 1;
  
  if (layer > 1000)
    layerdisk = module + 5;
  
  assert(layerdisk >= 0);
  assert(layerdisk < 11);
  
  int isec = dtc[0] - '0';
  
  string dtcbase = dtc.substr(2, dtc.size() - 2);
  if (dtc[0] == 'n') {
    dtcbase = "neg_" + dtc.substr(6, dtc.size() - 6);
    isec = dtc[4] - '0';
  }
  
  double phisec = Trklet::phiRange(phi - isec * 2 * M_PI / 9.0);  //Nonant cabling hardcoded here
  
  assert(dtcranges.find(dtcbase) != dtcranges.end());
  
  dtcranges[dtcbase].addphi(phisec, layerdisk);
}

void Cabling::writephirange() {
  ofstream out("dtcphirange.txt");
  
  std::map<string, DTC>::const_iterator it = dtcranges.begin();
  for (; it != dtcranges.end(); ++it) {
    for (unsigned int i = 0; i < 11; i++) {
      double min = it->second.min(i);
      double max = it->second.max(i);
      if (min < max) {
	out << it->first << " " << i + 1 << " " << min << " " << max << endl;
        }
    }
  }
}

std::vector<string> Cabling::DTCs() const {
  std::vector<string> tmp;
  
  for (auto it = dtcs.begin(); it != dtcs.end(); ++it) {
    tmp.push_back(it->first);
  }
  
  return tmp;
}

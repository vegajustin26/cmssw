#include "L1Trigger/TrackFindingTracklet/interface/Cabling.h"
#include "L1Trigger/TrackFindingTracklet/interface/DTCLink.h"
#include "L1Trigger/TrackFindingTracklet/interface/DTC.h"
#include "L1Trigger/TrackFindingTracklet/interface/Util.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;
using namespace trklet;

Cabling::Cabling(string dtcconfig, string moduleconfig, const Settings *settings) {

  settings_=settings;
  
  ifstream indtc(dtcconfig.c_str());
  assert(indtc.good());

  string dtc;
  int isec;

  while (indtc.good()) {
    indtc >> dtc >> isec;

    if (!indtc.good())
      continue;

    if (dtcs_.find(dtc) == dtcs_.end()) {
      dtcs_[dtc].setName(dtc);
    }

    dtcs_[dtc].addSec(isec);

    string dtcbase = dtc.substr(2, dtc.size() - 2);
    if (dtc[0] == 'n') {
      dtcbase = "neg_" + dtc.substr(6, dtc.size() - 6);
    }
    if (dtcranges_.find(dtcbase) == dtcranges_.end()) {
      dtcranges_[dtcbase].setName(dtcbase);
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
    modules_[layer][ladder][module] = dtc;
  }
}

const string& Cabling::dtc(int layer, int ladder, int module) const {
  auto it1 = modules_.find(layer);
  assert(it1 != modules_.end());
  auto it2 = it1->second.find(ladder);
  assert(it2 != it1->second.end());
  auto it3 = it2->second.find(module);
  if (it3 == it2->second.end()) {
    edm::LogPrint("Tracklet") << "Could not find stub " << layer << " " << ladder << " " << module;
    assert(0);
  }
  return it3->second;
}

void Cabling::addphi(const string& dtc, double phi, int layer, int module) {
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

  double phisec = trklet::phiRange(phi - isec * settings_->dphisector());

  assert(dtcranges_.find(dtcbase) != dtcranges_.end());

  dtcranges_[dtcbase].addphi(phisec, layerdisk);
}

void Cabling::writephirange() const {
  ofstream out("dtcphirange.txt");

  for (auto&& it : dtcranges_) {
    for (unsigned int i = 0; i < N_LAYERDISK; i++) {
      double min = it.second.min(i);
      double max = it.second.max(i);
      if (min < max) {
        out << it.first << " " << i + 1 << " " << min << " " << max << endl;
      }
    }
  }
}

std::vector<string> Cabling::DTCs() const {
  std::vector<string> tmp;

  for (const auto& it : dtcs_) {
    tmp.push_back(it.first);
  }

  return tmp;
}

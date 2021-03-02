#include "L1Trigger/TrackFindingTracklet/interface/Sector.h"
#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"

#include "L1Trigger/TrackFindingTracklet/interface/InputLinkMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllInnerStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubsTEMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubsMEMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/StubPairsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/StubTripletsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletParametersMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletProjectionsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllProjectionsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMProjectionsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/CandidateMatchMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/FullMatchMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackFitMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/CleanTrackMemory.h"

#include "L1Trigger/TrackFindingTracklet/interface/VMRouterCM.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMRouter.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletEngine.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletEngineDisplaced.h"
#include "L1Trigger/TrackFindingTracklet/interface/TripletEngine.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletCalculator.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletProcessor.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletCalculatorDisplaced.h"
#include "L1Trigger/TrackFindingTracklet/interface/ProjectionRouter.h"
#include "L1Trigger/TrackFindingTracklet/interface/MatchEngine.h"
#include "L1Trigger/TrackFindingTracklet/interface/MatchCalculator.h"
#include "L1Trigger/TrackFindingTracklet/interface/MatchProcessor.h"
#include "L1Trigger/TrackFindingTracklet/interface/FitTrack.h"
#include "L1Trigger/TrackFindingTracklet/interface/PurgeDuplicate.h"
#include "L1Trigger/TrackFindingTracklet/interface/Util.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMRouterPhiCorrTable.h"

using namespace std;
using namespace trklet;

Sector::Sector(unsigned int i, Settings const& settings, Globals* globals) : settings_(settings), globals_(globals) {
  isector_ = i;
  double dphi = 2 * M_PI / N_SECTOR;
  double dphiHG = 0.5 * settings_.dphisectorHG() - M_PI / N_SECTOR;
  phimin_ = isector_ * dphi - dphiHG;
  phimax_ = phimin_ + dphi + 2 * dphiHG;
  phimin_ -= M_PI / N_SECTOR;
  phimax_ -= M_PI / N_SECTOR;
  phimin_ = reco::reduceRange(phimin_);
  phimax_ = reco::reduceRange(phimax_);
  if (phimin_ > phimax_) {
    phimin_ -= 2 * M_PI;
  }
}

Sector::~Sector() = default;

bool Sector::addStub(L1TStub stub, string dtc) {

  unsigned int layerdisk=stub.layerdisk();

  double stubphi = stub.phi();
  if (stubphi-phimin_>M_PI){
    stubphi-=2*M_PI;
  }

  if (layerdisk < N_LAYER && globals_->phiCorr(layerdisk) == nullptr) {
    globals_->phiCorr(layerdisk) = new VMRouterPhiCorrTable(settings_);

    int nbits = 3;
    if (layerdisk >= N_PSLAYER)
      nbits = 4;
    globals_->phiCorr(layerdisk)->init(layerdisk + 1, nbits, 3);
  }
  
  Stub fpgastub(stub, settings_, *globals_);

  if (layerdisk < N_LAYER) {
    FPGAWord r = fpgastub.r();
    int bendbin = fpgastub.bend().value();
    int rbin = (r.value() + (1 << (r.nbits() - 1))) >> (r.nbits() - 3);
    const VMRouterPhiCorrTable& phiCorrTable = *globals_->phiCorr(layerdisk);
    int iphicorr = phiCorrTable.getphiCorrValue(bendbin, rbin);
    fpgastub.setPhiCorr(iphicorr);
  }
  
  FPGAWord phi=fpgastub.phicorr();
  int ireg=phi.value()>>(phi.nbits()-settings_.nbitsallstubs(layerdisk));

  int nadd=0;
  for (unsigned int i = 0; i < IL_.size(); i++) {
    const string& name=IL_[i]->getName();
    if (name.find("_" + dtc) == string::npos) continue;
    if ( (name[3]=='L' && name[4]-'1'==(int)layerdisk) ||  (name[3]=='D' && name[4]-'1'==(int)layerdisk-6) ) {
      if (name[8]-'A'==ireg){
	IL_[i]->addStub(stub, fpgastub);
	nadd++;
      }
    }
  }

  assert(nadd==1);

  if (settings_.writeMem() && isector_ == (int)settings_.writememsect() ) {
    writeLink(fpgastub);
  }
  
  return true;
}

void Sector::writeLinkNewEvent(int event) {

  if (DTCLink_ofstreams_.size()==0) {
    //FIXME should be in settings
    static string dtcbasenames[12]={"PS10G_1","PS10G_2","PS10G_3","PS10G_4","PS_1","PS_2","2S_1","2S_2","2S_3","2S_4","2S_5","2S_6"};     
    for (const auto& dtcbasename: dtcbasenames) {
      for (int neg=0;neg<2;neg++) {
	string dtcnametmp=(neg==0?"":"neg")+dtcbasename;
	for (int ab=0;ab<2;ab++) {
	  string dtcname=dtcnametmp+(ab==0?"_A":"_B");	
	  string fname = "../data/MemPrints/InputStubs/Link_"+dtcname+"_"+to_string(isector_+1)+".dat";
	  ofstream* out = new ofstream;
	  out->open(fname.c_str());
	  DTCLink_ofstreams_[dtcname] = out;
	}
      }
    }
  }

  for(auto& out:DTCLink_ofstreams_){
    FPGAWord bx(event%8,3,true);
    (*out.second) << "BX "<<bx.str()<<" Event : "<< event << endl;
  }

}

void Sector::writeLink(const Stub& fpgastub) {

  //FIXME should be in settings
  static map<string, vector<int> > dtclayers{
    {"PS10G_1", {0,6,8,10}},
    {"PS10G_2", {0,7,9}},
    {"PS10G_3", {1,7}},
    {"PS10G_4", {6,8,10}},
    {"PS_1", {2,7}},
    {"PS_2", {2,9}},
    {"2S_1", {3,4}},
    {"2S_2", {4}},
    {"2S_3", {5}},
    {"2S_4", {5,8}},
    {"2S_5", {6,9}},
    {"2S_6", {7,10}}
  };
  
  string dtcname=fpgastub.l1tstub()->DTClink();
  int layerdisk=fpgastub.l1tstub()->layerdisk();

  int start=dtcname.substr(0,3)=="neg"?3:0;
  
  string dtcbase=dtcname.substr(start,dtcname.size()-2-start);

  vector<int> layers=dtclayers[dtcbase];

  int lcode=-1;
  for(unsigned int index=0;index<layers.size();index++) {
    if (layerdisk==layers[index]) {
      lcode=index;
    }
  }
  assert(lcode!=-1);
  
  assert(DTCLink_ofstreams_.find(dtcname) != DTCLink_ofstreams_.end());

  FPGAWord ldcode(lcode,2,true);

  string dataword = "1|" + ldcode.str() + "|" + fpgastub.str();

  (*DTCLink_ofstreams_.find(dtcname)->second) << dataword << " " << trklet::hexFormat(dataword) << endl;

    
}

void Sector::addMem(string memType, string memName) {
  if (memType == "InputLink:") {
    addMemToVec(IL_, memName, settings_, isector_, phimin_, phimax_);
  } else if (memType == "AllStubs:") {
    addMemToVec(AS_, memName, settings_, isector_);
  } else if (memType == "AllInnerStubs:") {
    addMemToVec(AIS_, memName, settings_, isector_);
  } else if (memType == "VMStubsTE:") {
    addMemToVec(VMSTE_, memName, settings_, isector_);
  } else if (memType == "VMStubsME:") {
    addMemToVec(VMSME_, memName, settings_, isector_);
  } else if (memType == "StubPairs:" || memType == "StubPairsDisplaced:") {
    addMemToVec(SP_, memName, settings_, isector_);
  } else if (memType == "StubTriplets:") {
    addMemToVec(ST_, memName, settings_, isector_);
  } else if (memType == "TrackletParameters:") {
    addMemToVec(TPAR_, memName, settings_, isector_);
  } else if (memType == "TrackletProjections:") {
    addMemToVec(TPROJ_, memName, settings_, isector_);
  } else if (memType == "AllProj:") {
    addMemToVec(AP_, memName, settings_, isector_);
  } else if (memType == "VMProjections:") {
    addMemToVec(VMPROJ_, memName, settings_, isector_);
  } else if (memType == "CandidateMatch:") {
    addMemToVec(CM_, memName, settings_, isector_);
  } else if (memType == "FullMatch:") {
    addMemToVec(FM_, memName, settings_, isector_);
  } else if (memType == "TrackFit:") {
    addMemToVec(TF_, memName, settings_, isector_, phimin_, phimax_);
  } else if (memType == "CleanTrack:") {
    addMemToVec(CT_, memName, settings_, isector_, phimin_, phimax_);
  } else {
    edm::LogPrint("Tracklet") << "Don't know of memory type: " << memType;
    exit(0);
  }
}

void Sector::addProc(string procType, string procName) {
  if (procType == "VMRouter:") {
    addProcToVec(VMR_, procName, settings_, globals_, isector_);
  } else if (procType == "VMRouterCM:") {
    addProcToVec(VMRCM_, procName, settings_, globals_, isector_);
  } else if (procType == "TrackletEngine:") {
    addProcToVec(TE_, procName, settings_, globals_, isector_);
  } else if (procType == "TrackletEngineDisplaced:") {
    addProcToVec(TED_, procName, settings_, globals_, isector_);
  } else if (procType == "TripletEngine:") {
    addProcToVec(TRE_, procName, settings_, globals_, isector_);
  } else if (procType == "TrackletCalculator:") {
    addProcToVec(TC_, procName, settings_, globals_, isector_);
  } else if (procType == "TrackletProcessor:") {
    addProcToVec(TP_, procName, settings_, globals_, isector_);
  } else if (procType == "TrackletCalculatorDisplaced:") {
    addProcToVec(TCD_, procName, settings_, globals_, isector_);
  } else if (procType == "ProjectionRouter:") {
    addProcToVec(PR_, procName, settings_, globals_, isector_);
  } else if (procType == "MatchEngine:") {
    addProcToVec(ME_, procName, settings_, globals_, isector_);
  } else if (procType == "MatchCalculator:" ||
             procType == "DiskMatchCalculator:") {  //TODO should not be used in configurations
    addProcToVec(MC_, procName, settings_, globals_, isector_);
  } else if (procType == "MatchProcessor:") {
    addProcToVec(MP_, procName, settings_, globals_, isector_);
  } else if (procType == "FitTrack:") {
    addProcToVec(FT_, procName, settings_, globals_, isector_);
  } else if (procType == "PurgeDuplicate:") {
    addProcToVec(PD_, procName, settings_, globals_, isector_);
  } else {
    edm::LogPrint("Tracklet") << "Don't know of processing type: " << procType;
    exit(0);
  }
}

void Sector::addWire(string mem, string procinfull, string procoutfull) {
  stringstream ss1(procinfull);
  string procin, output;
  getline(ss1, procin, '.');
  getline(ss1, output);

  stringstream ss2(procoutfull);
  string procout, input;
  getline(ss2, procout, '.');
  getline(ss2, input);

  MemoryBase* memory = getMem(mem);

  if (!procin.empty()) {
    ProcessBase* inProc = getProc(procin);
    inProc->addOutput(memory, output);
  }

  if (!procout.empty()) {
    ProcessBase* outProc = getProc(procout);
    outProc->addInput(memory, input);
  }
}

ProcessBase* Sector::getProc(string procName) {
  auto it = Processes_.find(procName);

  if (it != Processes_.end()) {
    return it->second;
  }
  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " Could not find process : " << procName << endl;
  return nullptr;
}

MemoryBase* Sector::getMem(string memName) {
  auto it = Memories_.find(memName);

  if (it != Memories_.end()) {
    return it->second;
  }
  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " Could not find memory : " << memName;
  return nullptr;
}

void Sector::writeInputStubs(bool first) {
  for (auto& i : IL_) {
    i->writeStubs(first);
  }
}

void Sector::writeVMSTE(bool first) {
  for (auto& i : VMSTE_) {
    i->writeStubs(first);
  }
}

void Sector::writeVMSME(bool first) {
  for (auto& i : VMSME_) {
    i->writeStubs(first);
  }
}

void Sector::writeAS(bool first) {
  for (auto& i : AS_) {
    i->writeStubs(first);
  }
}

void Sector::writeAIS(bool first) {
  for (auto& i : AIS_) {
    i->writeStubs(first);
  }
}

void Sector::writeSP(bool first) {
  for (auto& i : SP_) {
    i->writeSP(first);
  }
}

void Sector::writeST(bool first) {
  for (auto& i : ST_) {
    i->writeST(first);
  }
}

void Sector::writeTPAR(bool first) {
  for (auto& i : TPAR_) {
    i->writeTPAR(first);
  }
}

void Sector::writeTPROJ(bool first) {
  for (auto& i : TPROJ_) {
    i->writeTPROJ(first);
  }
}

void Sector::writeAP(bool first) {
  for (auto& i : AP_) {
    i->writeAP(first);
  }
}

void Sector::writeVMPROJ(bool first) {
  for (auto& i : VMPROJ_) {
    i->writeVMPROJ(first);
  }
}

void Sector::writeCM(bool first) {
  for (auto& i : CM_) {
    i->writeCM(first);
  }
}

void Sector::writeMC(bool first) {
  for (auto& i : FM_) {
    i->writeMC(first);
  }
}

void Sector::writeTF(bool first) {
  for (auto& i : TF_) {
    i->writeTF(first);
  }
}

void Sector::writeCT(bool first) {
  for (auto& i : CT_) {
    i->writeCT(first);
  }
}

void Sector::clean() {
  if (settings_.writeMonitorData("NMatches")) {
    int matchesL1 = 0;
    int matchesL3 = 0;
    int matchesL5 = 0;
    for (auto& i : TPAR_) {
      i->writeMatches(globals_, matchesL1, matchesL3, matchesL5);
    }
    globals_->ofstream("nmatchessector.txt") << matchesL1 << " " << matchesL3 << " " << matchesL5 << endl;
  }

  for (auto& mem : MemoriesV_) {
    mem->clean();
  }
}

void Sector::executeVMR() {
  if (settings_.writeMonitorData("IL")) {
    ofstream& out = globals_->ofstream("inputlink.txt");
    for (auto& i : IL_) {
      out << i->getName() << " " << i->nStubs() << endl;
    }
  }
  for (auto& i : VMR_) {
    i->execute();
  }
  for (auto& i : VMRCM_) {
    i->execute();
  }
}

void Sector::executeTE() {
  for (auto& i : TE_) {
    i->execute();
  }
}

void Sector::executeTED() {
  for (auto& i : TED_) {
    i->execute();
  }
}

void Sector::executeTRE() {
  for (auto& i : TRE_) {
    i->execute();
  }
}

void Sector::executeTP() {
  for (auto& i : TP_) {
    i->execute();
  }
}

void Sector::executeTC() {
  for (auto& i : TC_) {
    i->execute();
  }

  if (settings_.writeMonitorData("TrackProjOcc")) {
    ofstream& out = globals_->ofstream("trackprojocc.txt");
    for (auto& i : TPROJ_) {
      out << i->getName() << " " << i->nTracklets() << endl;
    }
  }
}

void Sector::executeTCD() {
  for (auto& i : TCD_) {
    i->execute();
  }
}

void Sector::executePR() {
  for (auto& i : PR_) {
    i->execute();
  }
}

void Sector::executeME() {
  for (auto& i : ME_) {
    i->execute();
  }
}

void Sector::executeMC() {
  for (auto& i : MC_) {
    i->execute();
  }
}

void Sector::executeMP() {
  for (auto& i : MP_) {
    i->execute();
  }
}

void Sector::executeFT() {
  for (auto& i : FT_) {
    i->execute();
  }
}

void Sector::executePD(std::vector<Track*>& tracks) {
  for (auto& i : PD_) {
    i->execute(tracks);
  }
}

std::vector<Tracklet*> Sector::getAllTracklets() const {
  std::vector<Tracklet*> tmp;
  for (auto& tpar : TPAR_) {
    for (unsigned int j = 0; j < tpar->nTracklets(); j++) {
      tmp.push_back(tpar->getTracklet(j));
    }
  }
  return tmp;
}

std::vector<const Stub*> Sector::getStubs() const {
  std::vector<const Stub*> tmp;

  for (auto& imem : IL_) {
    for (unsigned int istub = 0; istub < imem->nStubs(); istub++) {
      tmp.push_back(imem->getStub(istub));
    }
  }

  return tmp;
}

std::unordered_set<int> Sector::seedMatch(int itp) const {
  std::unordered_set<int> tmpSeeds;
  for (auto& i : TPAR_) {
    unsigned int nTracklet = i->nTracklets();
    for (unsigned int j = 0; j < nTracklet; j++) {
      if (i->getTracklet(j)->tpseed() == itp) {
        tmpSeeds.insert(i->getTracklet(j)->getISeed());
      }
    }
  }
  return tmpSeeds;
}

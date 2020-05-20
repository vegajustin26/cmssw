#include "L1Trigger/TrackFindingTracklet/interface/Sector.h"
#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"

#include "L1Trigger/TrackFindingTracklet/interface/InputLinkMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"
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

using namespace std;
using namespace trklet;

Sector::Sector(unsigned int i, const Settings* settings, Globals* globals) : settings_(settings), globals_(globals) {
  isector_ = i;
  double dphi = 2 * M_PI / N_SECTOR;
  double dphiHG = 0.5 * settings_->dphisectorHG() - M_PI / N_SECTOR;
  phimin_ = isector_ * dphi - dphiHG;
  phimax_ = phimin_ + dphi + 2 * dphiHG;
  phimin_ -= M_PI / N_SECTOR;
  phimax_ -= M_PI / N_SECTOR;
  phimin_ = phiRange(phimin_);
  phimax_ = phiRange(phimax_);
  if (phimin_ > phimax_)
    phimin_ -= 2 * M_PI;
}

Sector::~Sector() {
  for (auto& mem : MemoriesV_) {
    mem->clean();
  }
}

bool Sector::addStub(L1TStub stub, const string& dtc) {
  bool add = false;

  double phi = stub.phi();
  double dphi = 0.5 * settings_->dphisectorHG() - M_PI / N_SECTOR;

  std::map<string, std::vector<int> >& ILindex = globals_->ILindex();
  std::vector<int>& tmp = ILindex[dtc];
  if (tmp.empty()) {
    for (unsigned int i = 0; i < IL_.size(); i++) {
      if (IL_[i]->getName().find("_" + dtc) != string::npos) {
        tmp.push_back(i);
      }
    }
  }

  if (((phi > phimin_ - dphi) && (phi < phimax_ + dphi)) ||
      ((phi > 2 * M_PI + phimin_ - dphi) && (phi < 2 * M_PI + phimax_ + dphi))) {
    Stub fpgastub(stub, settings_, phimin_, phimax_);
    std::vector<int>& tmp = ILindex[dtc];
    assert(tmp.size() != 0);
    for (int i : tmp) {
      if (IL_[i]->addStub(settings_, globals_, stub, fpgastub, dtc))
        add = true;
    }
  }

  return add;
}

void Sector::addMem(const string& memType, const string& memName) {
  if (memType == "InputLink:") {
    addMemToVec(IL_, new InputLinkMemory(memName, settings_, isector_, phimin_, phimax_), memName);
  } else if (memType == "AllStubs:") {
    addMemToVec(AS_, new AllStubsMemory(memName, settings_, isector_), memName);
  } else if (memType == "VMStubsTE:") {
    addMemToVec(VMSTE_, new VMStubsTEMemory(memName, settings_, isector_), memName);
  } else if (memType == "VMStubsME:") {
    addMemToVec(VMSME_, new VMStubsMEMemory(memName, settings_, isector_), memName);
  } else if (memType == "StubPairs:" || memType == "StubPairsDisplaced:") {
    addMemToVec(SP_, new StubPairsMemory(memName, settings_, isector_), memName);
  } else if (memType == "StubTriplets:") {
    addMemToVec(ST_, new StubTripletsMemory(memName, settings_, isector_), memName);
  } else if (memType == "TrackletParameters:") {
    addMemToVec(TPAR_, new TrackletParametersMemory(memName, settings_, isector_), memName);
  } else if (memType == "TrackletProjections:") {
    addMemToVec(TPROJ_, new TrackletProjectionsMemory(memName, settings_, isector_), memName);
  } else if (memType == "AllProj:") {
    addMemToVec(AP_, new AllProjectionsMemory(memName, settings_, isector_), memName);
  } else if (memType == "VMProjections:") {
    addMemToVec(VMPROJ_, new VMProjectionsMemory(memName, settings_, isector_), memName);
  } else if (memType == "CandidateMatch:") {
    addMemToVec(CM_, new CandidateMatchMemory(memName, settings_, isector_), memName);
  } else if (memType == "FullMatch:") {
    addMemToVec(FM_, new FullMatchMemory(memName, settings_, isector_), memName);
  } else if (memType == "TrackFit:") {
    addMemToVec(TF_, new TrackFitMemory(memName, settings_, isector_, phimin_, phimax_), memName);
  } else if (memType == "CleanTrack:") {
    addMemToVec(CT_, new CleanTrackMemory(memName, settings_, isector_, phimin_, phimax_), memName);
  } else {
    edm::LogPrint("Tracklet") << "Don't know of memory type: " << memType;
    exit(0);
  }
}

void Sector::addProc(const string& procType, const string& procName) {
  if (procType == "VMRouter:") {
    addProcToVec(VMR_, new VMRouter(procName, settings_, globals_, isector_), procName);
  } else if (procType == "TrackletEngine:") {
    addProcToVec(TE_, new TrackletEngine(procName, settings_, globals_, isector_), procName);
  } else if (procType == "TrackletEngineDisplaced:") {
    addProcToVec(TED_, new TrackletEngineDisplaced(procName, settings_, globals_, isector_), procName);
  } else if (procType == "TripletEngine:") {
    addProcToVec(TRE_, new TripletEngine(procName, settings_, globals_, isector_), procName);
  } else if (procType == "TrackletCalculator:") {
    addProcToVec(TC_, new TrackletCalculator(procName, settings_, globals_, isector_), procName);
  } else if (procType == "TrackletProcessor:") {
    addProcToVec(TP_, new TrackletProcessor(procName, settings_, globals_, isector_), procName);
  } else if (procType == "TrackletCalculatorDisplaced:") {
    addProcToVec(TCD_, new TrackletCalculatorDisplaced(procName, settings_, globals_, isector_), procName);
  } else if (procType == "ProjectionRouter:") {
    addProcToVec(PR_, new ProjectionRouter(procName, settings_, globals_, isector_), procName);
  } else if (procType == "MatchEngine:") {
    addProcToVec(ME_, new MatchEngine(procName, settings_, globals_, isector_), procName);
  } else if (procType == "MatchCalculator:" ||
             procType == "DiskMatchCalculator:") {  //TODO should not be used in configurations
    addProcToVec(MC_, new MatchCalculator(procName, settings_, globals_, isector_), procName);
  } else if (procType == "MatchProcessor:") {
    addProcToVec(MP_, new MatchProcessor(procName, settings_, globals_, isector_), procName);
  } else if (procType == "FitTrack:") {
    addProcToVec(FT_, new FitTrack(procName, settings_, globals_, isector_), procName);
  } else if (procType == "PurgeDuplicate:") {
    addProcToVec(PD_, new PurgeDuplicate(procName, settings_, globals_, isector_), procName);
  } else {
    edm::LogPrint("Tracklet") << "Don't know of processing type: " << procType;
    exit(0);
  }
}

void Sector::addWire(string mem, const string& procinfull, const string& procoutfull) {
  stringstream ss1(procinfull);
  string procin, output;
  getline(ss1, procin, '.');
  getline(ss1, output);

  stringstream ss2(procoutfull);
  string procout, input;
  getline(ss2, procout, '.');
  getline(ss2, input);

  MemoryBase* memory = getMem(std::move(mem));

  if (procin != "") {
    ProcessBase* inProc = getProc(procin);
    inProc->addOutput(memory, output);
  }

  if (procout != "") {
    ProcessBase* outProc = getProc(procout);
    outProc->addInput(memory, input);
  }
}

ProcessBase* Sector::getProc(const string& procName) {
  auto it = Processes_.find(procName);

  if (it != Processes_.end()) {
    return it->second.get();
  }
  edm::LogPrint("Tracklet") << "Could not find process with name : " << procName;
  assert(0);
  return nullptr;
}

MemoryBase* Sector::getMem(const string& memName) {
  auto it = Memories_.find(memName);

  if (it != Memories_.end()) {
    return it->second.get();
  }
  edm::LogPrint("Tracklet") << "Could not find memory with name : " << memName;
  assert(0);
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
  if (settings_->writeMonitorData("NMatches")) {
    int matchesL1 = 0;
    int matchesL3 = 0;
    int matchesL5 = 0;
    for (auto& i : TPAR_) {
      i->writeMatches(globals_, matchesL1, matchesL3, matchesL5);
    }
    globals_->ofstream("nmatchessector.txt") << matchesL1 << " " << matchesL3 << " " << matchesL5 << endl;
  }

  for (auto& i : MemoriesV_) {
    i->clean();
  }
}

void Sector::executeVMR() {
  if (settings_->writeMonitorData("IL")) {
    ofstream& out = globals_->ofstream("inputlink.txt");
    for (auto& i : IL_) {
      out << i->getName() << " " << i->nStubs() << endl;
    }
  }

  for (auto& i : VMR_) {
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

  if (settings_->writeMonitorData("TrackProjOcc")) {
    ofstream& out = globals_->ofstream("trackprojocc.txt");
    for (unsigned int i = 0; i < TPROJ_.size(); i++) {
      out << TPROJ_[i]->getName() << " " << TPROJ_[i]->nTracklets() << endl;
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
  for (unsigned int i = 0; i < TPAR_.size(); i++) {
    for (unsigned int j = 0; j < TPAR_[i]->nTracklets(); j++) {
      tmp.push_back(TPAR_[i]->getTracklet(j));
    }
  }
  return tmp;
}

std::vector<const Stub*> Sector::getStubs() const {
  std::vector<const Stub*> tmp;

  for (auto imem : IL_) {
    for (unsigned int istub = 0; istub < imem->nStubs(); istub++) {
      tmp.push_back(imem->getStub(istub));
    }
  }
  return tmp;
}

std::unordered_set<int> Sector::seedMatch(int itp) const {
  std::unordered_set<int> tmpSeeds;
  for (unsigned int i = 0; i < TPAR_.size(); i++) {
    unsigned int nTracklet = TPAR_[i]->nTracklets();
    for (unsigned int j = 0; j < nTracklet; j++) {
      if (TPAR_[i]->getTracklet(j)->tpseed() == itp) {
        tmpSeeds.insert(TPAR_[i]->getTracklet(j)->getISeed());
      }
    }
  }
  return tmpSeeds;
}

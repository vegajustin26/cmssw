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
using namespace Trklet;

Sector::Sector(unsigned int i, const Settings* settings, Globals* globals) : settings_(settings), globals_(globals) {
  isector_ = i;
  double dphi = 2 * M_PI / settings_->NSector();
  double dphiHG = 0.5 * settings_->dphisectorHG() - M_PI / settings_->NSector();
  phimin_ = isector_ * dphi - dphiHG;
  phimax_ = phimin_ + dphi + 2 * dphiHG;
  phimin_ -= M_PI / settings_->NSector();
  phimax_ -= M_PI / settings_->NSector();
  phimin_ = Trklet::phiRange(phimin_);
  phimax_ = Trklet::phiRange(phimax_);
  if (phimin_ > phimax_)
    phimin_ -= 2 * M_PI;
}

Sector::~Sector() {
  for (const auto& p : Processes_) {
    ProcessBase* proc = p.second;
    delete proc;
  }
  for (const auto& m : Memories_) {
    MemoryBase* mem = m.second;
    delete mem;
  }
}

bool Sector::addStub(L1TStub stub, string dtc) {
  bool add = false;

  double phi = stub.phi();
  double dphi = 0.5 * settings_->dphisectorHG() - M_PI / settings_->NSector();

  std::map<string, std::vector<int> >& ILindex = globals_->ILindex();
  std::vector<int>& tmp = ILindex[dtc];
  if (tmp.size() == 0) {
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
    for (unsigned int i = 0; i < tmp.size(); i++) {
      if (IL_[tmp[i]]->addStub(settings_, globals_, stub, fpgastub, dtc))
        add = true;
    }
  }

  return add;
}

void Sector::addMem(string memType, string memName) {
  if (memType == "InputLink:") {
    IL_.push_back(new InputLinkMemory(memName, settings_, isector_, phimin_, phimax_));
    Memories_[memName] = IL_.back();
    MemoriesV_.push_back(IL_.back());
  } else if (memType == "AllStubs:") {
    AS_.push_back(new AllStubsMemory(memName, settings_, isector_));
    Memories_[memName] = AS_.back();
    MemoriesV_.push_back(AS_.back());
  } else if (memType == "VMStubsTE:") {
    VMSTE_.push_back(new VMStubsTEMemory(memName, settings_, isector_));
    Memories_[memName] = VMSTE_.back();
    MemoriesV_.push_back(VMSTE_.back());
  } else if (memType == "VMStubsME:") {
    VMSME_.push_back(new VMStubsMEMemory(memName, settings_, isector_));
    Memories_[memName] = VMSME_.back();
    MemoriesV_.push_back(VMSME_.back());
  } else if (memType == "StubPairs:" || memType == "StubPairsDisplaced:") {
    SP_.push_back(new StubPairsMemory(memName, settings_, isector_));
    Memories_[memName] = SP_.back();
    MemoriesV_.push_back(SP_.back());
  } else if (memType == "StubTriplets:") {
    ST_.push_back(new StubTripletsMemory(memName, settings_, isector_));
    Memories_[memName] = ST_.back();
    MemoriesV_.push_back(ST_.back());
  } else if (memType == "TrackletParameters:") {
    TPAR_.push_back(new TrackletParametersMemory(memName, settings_, isector_));
    Memories_[memName] = TPAR_.back();
    MemoriesV_.push_back(TPAR_.back());
  } else if (memType == "TrackletProjections:") {
    TPROJ_.push_back(new TrackletProjectionsMemory(memName, settings_, isector_));
    Memories_[memName] = TPROJ_.back();
    MemoriesV_.push_back(TPROJ_.back());
  } else if (memType == "AllProj:") {
    AP_.push_back(new AllProjectionsMemory(memName, settings_, isector_));
    Memories_[memName] = AP_.back();
    MemoriesV_.push_back(AP_.back());
  } else if (memType == "VMProjections:") {
    VMPROJ_.push_back(new VMProjectionsMemory(memName, settings_, isector_));
    Memories_[memName] = VMPROJ_.back();
    MemoriesV_.push_back(VMPROJ_.back());
  } else if (memType == "CandidateMatch:") {
    CM_.push_back(new CandidateMatchMemory(memName, settings_, isector_));
    Memories_[memName] = CM_.back();
    MemoriesV_.push_back(CM_.back());
  } else if (memType == "FullMatch:") {
    FM_.push_back(new FullMatchMemory(memName, settings_, isector_));
    Memories_[memName] = FM_.back();
    MemoriesV_.push_back(FM_.back());
  } else if (memType == "TrackFit:") {
    TF_.push_back(new TrackFitMemory(memName, settings_, isector_, phimin_, phimax_));
    Memories_[memName] = TF_.back();
    MemoriesV_.push_back(TF_.back());
  } else if (memType == "CleanTrack:") {
    CT_.push_back(new CleanTrackMemory(memName, settings_, isector_, phimin_, phimax_));
    Memories_[memName] = CT_.back();
    MemoriesV_.push_back(CT_.back());
  } else {
    edm::LogPrint("Tracklet") << "Don't know of memory type: " << memType;
    exit(0);
  }
}

void Sector::addProc(string procType, string procName) {
  if (procType == "VMRouter:") {
    VMR_.push_back(new VMRouter(procName, settings_, globals_, isector_));
    Processes_[procName] = VMR_.back();
  } else if (procType == "TrackletEngine:") {
    TE_.push_back(new TrackletEngine(procName, settings_, globals_, isector_));
    Processes_[procName] = TE_.back();
  } else if (procType == "TrackletEngineDisplaced:") {
    TED_.push_back(new TrackletEngineDisplaced(procName, settings_, globals_, isector_));
    Processes_[procName] = TED_.back();
  } else if (procType == "TripletEngine:") {
    TRE_.push_back(new TripletEngine(procName, settings_, globals_, isector_));
    Processes_[procName] = TRE_.back();
  } else if (procType == "TrackletCalculator:") {
    TC_.push_back(new TrackletCalculator(procName, settings_, globals_, isector_));
    Processes_[procName] = TC_.back();
  } else if (procType == "TrackletProcessor:") {
    TP_.push_back(new TrackletProcessor(procName, settings_, globals_, isector_));
    Processes_[procName] = TP_.back();
  } else if (procType == "TrackletCalculatorDisplaced:") {
    TCD_.push_back(new TrackletCalculatorDisplaced(procName, settings_, globals_, isector_));
    Processes_[procName] = TCD_.back();
  } else if (procType == "ProjectionRouter:") {
    PR_.push_back(new ProjectionRouter(procName, settings_, globals_, isector_));
    Processes_[procName] = PR_.back();
  } else if (procType == "MatchEngine:") {
    ME_.push_back(new MatchEngine(procName, settings_, globals_, isector_));
    Processes_[procName] = ME_.back();
  } else if (procType == "MatchCalculator:" ||
             procType == "DiskMatchCalculator:") {  //TODO should not be used in configurations
    MC_.push_back(new MatchCalculator(procName, settings_, globals_, isector_));
    Processes_[procName] = MC_.back();
  } else if (procType == "MatchProcessor:") {
    MP_.push_back(new MatchProcessor(procName, settings_, globals_, isector_));
    Processes_[procName] = MP_.back();
  } else if (procType == "FitTrack:") {
    FT_.push_back(new FitTrack(procName, settings_, globals_, isector_));
    Processes_[procName] = FT_.back();
  } else if (procType == "PurgeDuplicate:") {
    PD_.push_back(new PurgeDuplicate(procName, settings_, globals_, isector_));
    Processes_[procName] = PD_.back();
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

  if (procin != "") {
    ProcessBase* inProc = getProc(procin);
    inProc->addOutput(memory, output);
  }

  if (procout != "") {
    ProcessBase* outProc = getProc(procout);
    outProc->addInput(memory, input);
  }
}

ProcessBase* Sector::getProc(string procName) {
  map<string, ProcessBase*>::iterator it = Processes_.find(procName);

  if (it != Processes_.end()) {
    return it->second;
  }
  edm::LogPrint("Tracklet") << "Could not find process with name : " << procName;
  assert(0);
  return 0;
}

MemoryBase* Sector::getMem(string memName) {
  map<string, MemoryBase*>::iterator it = Memories_.find(memName);

  if (it != Memories_.end()) {
    return it->second;
  }
  edm::LogPrint("Tracklet") << "Could not find memory with name : " << memName;
  assert(0);
  return 0;
}

void Sector::writeInputStubs(bool first) {
  for (unsigned int i = 0; i < IL_.size(); i++) {
    IL_[i]->writeStubs(first);
  }
}

void Sector::writeVMSTE(bool first) {
  for (unsigned int i = 0; i < VMSTE_.size(); i++) {
    VMSTE_[i]->writeStubs(first);
  }
}

void Sector::writeVMSME(bool first) {
  for (unsigned int i = 0; i < VMSME_.size(); i++) {
    VMSME_[i]->writeStubs(first);
  }
}

void Sector::writeAS(bool first) {
  for (unsigned int i = 0; i < AS_.size(); i++) {
    AS_[i]->writeStubs(first);
  }
}

void Sector::writeSP(bool first) {
  for (unsigned int i = 0; i < SP_.size(); i++) {
    SP_[i]->writeSP(first);
  }
}

void Sector::writeST(bool first) {
  for (unsigned int i = 0; i < ST_.size(); i++) {
    ST_[i]->writeST(first);
  }
}

void Sector::writeTPAR(bool first) {
  for (unsigned int i = 0; i < TPAR_.size(); i++) {
    TPAR_[i]->writeTPAR(first);
  }
}

void Sector::writeTPROJ(bool first) {
  for (unsigned int i = 0; i < TPROJ_.size(); i++) {
    TPROJ_[i]->writeTPROJ(first);
  }
}

void Sector::writeAP(bool first) {
  for (unsigned int i = 0; i < AP_.size(); i++) {
    AP_[i]->writeAP(first);
  }
}

void Sector::writeVMPROJ(bool first) {
  for (unsigned int i = 0; i < VMPROJ_.size(); i++) {
    VMPROJ_[i]->writeVMPROJ(first);
  }
}

void Sector::writeCM(bool first) {
  for (unsigned int i = 0; i < CM_.size(); i++) {
    CM_[i]->writeCM(first);
  }
}

void Sector::writeMC(bool first) {
  for (unsigned int i = 0; i < FM_.size(); i++) {
    FM_[i]->writeMC(first);
  }
}

void Sector::writeTF(bool first) {
  for (unsigned int i = 0; i < TF_.size(); ++i) {
    TF_[i]->writeTF(first);
  }
}

void Sector::writeCT(bool first) {
  for (unsigned int i = 0; i < CT_.size(); ++i) {
    CT_[i]->writeCT(first);
  }
}

void Sector::clean() {
  if (settings_->writeMonitorData("NMatches")) {
    int matchesL1 = 0;
    int matchesL3 = 0;
    int matchesL5 = 0;
    for (unsigned int i = 0; i < TPAR_.size(); i++) {
      TPAR_[i]->writeMatches(globals_, matchesL1, matchesL3, matchesL5);
    }
    globals_->ofstream("nmatchessector.txt") << matchesL1 << " " << matchesL3 << " " << matchesL5 << endl;
  }

  for (unsigned int i = 0; i < MemoriesV_.size(); i++) {
    MemoriesV_[i]->clean();
  }
}

void Sector::executeVMR() {
  if (settings_->writeMonitorData("IL")) {
    ofstream& out = globals_->ofstream("inputlink.txt");
    for (unsigned int i = 0; i < IL_.size(); i++) {
      out << IL_[i]->getName() << " " << IL_[i]->nStubs() << endl;
    }
  }

  for (unsigned int i = 0; i < VMR_.size(); i++) {
    VMR_[i]->execute();
  }
}

void Sector::executeTE() {
  for (unsigned int i = 0; i < TE_.size(); i++) {
    TE_[i]->execute();
  }
}

void Sector::executeTED() {
  for (unsigned int i = 0; i < TED_.size(); i++) {
    TED_[i]->execute();
  }
}

void Sector::executeTRE() {
  for (unsigned int i = 0; i < TRE_.size(); i++) {
    TRE_[i]->execute();
  }
}

void Sector::executeTP() {
  for (unsigned int i = 0; i < TP_.size(); i++) {
    TP_[i]->execute();
  }
}

void Sector::executeTC() {
  for (unsigned int i = 0; i < TC_.size(); i++) {
    TC_[i]->execute();
  }

  if (settings_->writeMonitorData("TrackProjOcc")) {
    ofstream& out = globals_->ofstream("trackprojocc.txt");
    for (unsigned int i = 0; i < TPROJ_.size(); i++) {
      out << TPROJ_[i]->getName() << " " << TPROJ_[i]->nTracklets() << endl;
    }
  }
}

void Sector::executeTCD() {
  for (unsigned int i = 0; i < TCD_.size(); i++) {
    TCD_[i]->execute();
  }
}

void Sector::executePR() {
  for (unsigned int i = 0; i < PR_.size(); i++) {
    PR_[i]->execute();
  }
}

void Sector::executeME() {
  for (unsigned int i = 0; i < ME_.size(); i++) {
    ME_[i]->execute();
  }
}

void Sector::executeMC() {
  for (unsigned int i = 0; i < MC_.size(); i++) {
    MC_[i]->execute();
  }
}

void Sector::executeMP() {
  for (unsigned int i = 0; i < MP_.size(); i++) {
    MP_[i]->execute();
  }
}

void Sector::executeFT() {
  for (unsigned int i = 0; i < FT_.size(); i++) {
    FT_[i]->execute();
  }
}

void Sector::executePD(std::vector<Track*>& tracks) {
  for (unsigned int i = 0; i < PD_.size(); i++) {
    PD_[i]->execute(tracks);
  }
}

bool Sector::foundTrack(ofstream& outres, L1SimTrack simtrk) {
  bool match = false;
  for (unsigned int i = 0; i < TF_.size(); i++) {
    if (TF_[i]->foundTrack(outres, simtrk))
      match = true;
  }
  return match;
}

std::vector<Tracklet*> Sector::getAllTracklets() {
  std::vector<Tracklet*> tmp;
  for (unsigned int i = 0; i < TPAR_.size(); i++) {
    for (unsigned int j = 0; j < TPAR_[i]->nTracklets(); j++) {
      tmp.push_back(TPAR_[i]->getFPGATracklet(j));
    }
  }
  return tmp;
}

std::vector<std::pair<Stub*, L1TStub*> > Sector::getStubs() const {
  std::vector<std::pair<Stub*, L1TStub*> > tmp;

  for (unsigned int imem = 0; imem < IL_.size(); imem++) {
    for (unsigned int istub = 0; istub < IL_[imem]->nStubs(); istub++) {
      tmp.push_back(IL_[imem]->getStub(istub));
    }
  }

  return tmp;
}

std::set<int> Sector::seedMatch(int itp) {
  std::set<int> tmpSeeds;
  for (unsigned int i = 0; i < TPAR_.size(); i++) {
    unsigned int nTracklet = TPAR_[i]->nTracklets();
    for (unsigned int j = 0; j < nTracklet; j++) {
      if (TPAR_[i]->getFPGATracklet(j)->tpseed() == itp) {
        tmpSeeds.insert(TPAR_[i]->getFPGATracklet(j)->getISeed());
      }
    }
  }
  return tmpSeeds;
}

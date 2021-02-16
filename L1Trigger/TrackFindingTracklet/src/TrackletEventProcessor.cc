#include "L1Trigger/TrackFindingTracklet/interface/TrackletEventProcessor.h"
#include "L1Trigger/TrackFindingTracklet/interface/SLHCEvent.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"
#include "L1Trigger/TrackFindingTracklet/interface/SLHCEvent.h"
#include "L1Trigger/TrackFindingTracklet/interface/Sector.h"
#include "L1Trigger/TrackFindingTracklet/interface/HistBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/Track.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletConfigBuilder.h"
#include "L1Trigger/TrackFindingTracklet/interface/IMATH_TrackletCalculator.h"

#include "DataFormats/Math/interface/deltaPhi.h"

#include <iomanip>
#include <filesystem>

using namespace trklet;
using namespace std;

TrackletEventProcessor::TrackletEventProcessor() : settings_(nullptr) {}

TrackletEventProcessor::~TrackletEventProcessor() {
  if (settings_ && settings_->bookHistos()) {
    histbase_->close();
  }
}

void TrackletEventProcessor::init(Settings const& theSettings) {
  settings_ = &theSettings;

  globals_ = make_unique<Globals>(*settings_);

  //Verify consistency
  if (settings_->kphi0pars() != globals_->ITC_L1L2()->phi0_final.K()) {
    throw cms::Exception("Inconsistency") << "phi0 conversion parameter inconsistency\n";
  }

  if (settings_->krinvpars() != globals_->ITC_L1L2()->rinv_final.K()) {
    throw cms::Exception("Inconsistency") << "ring conversion parameter inconsistency\n";
  }

  if (settings_->ktpars() != globals_->ITC_L1L2()->t_final.K()) {
    throw cms::Exception("Inconsistency") << "t conversion parameter inconsistency\n";
  }

  if (settings_->kphider() != globals_->ITC_L1L2()->der_phiL_final.K()) {
    throw cms::Exception("Inconsistency") << "t conversion parameter inconsistency:"
      <<settings_->kphider()/globals_->ITC_L1L2()->der_phiL_final.K()<<"\n";
  }

  
  if (settings_->debugTracklet()) {
    edm::LogVerbatim("Tracklet") << "========================================================= \n"
                                 << "Conversion factors for global coordinates: \n"
                                 << "z    kz            = " << settings_->kz() << "\n"
                                 << "r    kr            = " << settings_->kr() << "\n"
                                 << "phi  kphi1         = " << settings_->kphi1() << "\n"
                                 << "========================================================= \n"
                                 << "Conversion factors for track(let) parameters: \n"
                                 << "rinv krinvpars     = " << settings_->krinvpars() << "\n"
                                 << "phi0 kphi0pars     = " << settings_->kphi0pars() << "\n"
                                 << "d0   kd0pars       = " << settings_->kd0pars() << "\n"
                                 << "t    ktpars        = " << settings_->ktpars() << "\n"
                                 << "z0   kz0pars       = " << settings_->kz0pars() << "\n"
                                 << "========================================================= \n"
                                 << "phi0bitshift = " << settings_->phi0bitshift() << "\n"
                                 << "d0bitshift   = ??? \n"
                                 << "=========================================================";
  }

  if (settings_->bookHistos()) {
    histbase_ = new HistBase;
    histbase_->open();
    histbase_->bookLayerResidual();
    histbase_->bookDiskResidual();
    histbase_->bookTrackletParams();
    histbase_->bookSeedEff();

    globals_->histograms() = histbase_;
  }

  // create the sector processors (1 sector processor = 1 board)
  sectors_.resize(N_SECTOR);

  for (unsigned int i = 0; i < N_SECTOR; i++) {
    sectors_[i] = make_unique<Sector>(i, *settings_, globals_.get());
  }

  if (settings_->extended()) {
 
    ifstream inmem(settings_->memoryModulesFile().c_str());
    assert(inmem.good());

    ifstream inproc(settings_->processingModulesFile().c_str());
    assert(inproc.good());

    ifstream inwire(settings_->wiresFile().c_str());
    assert(inwire.good());
    
    configure(inwire,inmem,inproc);

  } else {
    TrackletConfigBuilder config(*settings_);

    //Write configurations to file.
    if (settings_->writeConfig()){
      std::ofstream wires=openfile(settings_->tablePath(), "wires.dat", __FILE__, __LINE__);
      std::ofstream memories=openfile(settings_->tablePath(), "memories.dat", __FILE__, __LINE__);
      std::ofstream modules=openfile(settings_->tablePath(), "modules.dat", __FILE__, __LINE__);
      
      config.writeAll(wires,memories,modules);
    }
    
    std::stringstream wires;
    std::stringstream memories;
    std::stringstream modules;

    config.writeAll(wires,memories,modules);
    configure(wires,memories,modules);

  }

}


void TrackletEventProcessor::configure(istream& inwire, istream& inmem, istream& inproc) {

  // get the memory modules
  if (settings_->debugTracklet()) {
    edm::LogVerbatim("Tracklet") << "Will read memory modules";
  }

  while (inmem.good()) {
    string memType, memName, size;
    inmem >> memType >> memName >> size;
    if (!inmem.good())
      continue;
    if (settings_->writetrace()) {
      edm::LogVerbatim("Tracklet") << "Read memory: " << memType << " " << memName;
    }
    for (auto& sector : sectors_) {
      sector->addMem(memType, memName);
    }
  }

  // get the processing modules
  if (settings_->debugTracklet()) {
    edm::LogVerbatim("Tracklet") << "Will read processing modules";
  }


  while (inproc.good()) {
    string procType, procName;
    inproc >> procType >> procName;
    if (!inproc.good())
      continue;
    if (settings_->writetrace()) {
      edm::LogVerbatim("Tracklet") << "Read process: " << procType << " " << procName;
    }
    for (auto& sector : sectors_) {
      sector->addProc(procType, procName);
    }
  }

  // get the wiring information
  if (settings_->debugTracklet()) {
    edm::LogVerbatim("Tracklet") << "Will read wiring information";
  }


  while (inwire.good()) {
    string line;
    getline(inwire, line);
    if (!inwire.good())
      continue;
    if (settings_->writetrace()) {
      edm::LogVerbatim("Tracklet") << "Line : " << line;
    }
    stringstream ss(line);
    string mem, tmp1, procin, tmp2, procout;
    ss >> mem >> tmp1 >> procin;
    if (procin == "output=>") {
      procin = "";
      ss >> procout;
    } else {
      ss >> tmp2 >> procout;
    }

    for (auto& sector : sectors_) {
      sector->addWire(mem, procin, procout);
    }
  }

}


void TrackletEventProcessor::event(SLHCEvent& ev) {

  globals_->event() = &ev;

  tracks_.clear();

  eventnum_++;
  bool first = (eventnum_ == 1);

  cleanTimer_.start();
  for (unsigned int k = 0; k < N_SECTOR; k++) {
    sectors_[k]->clean();
  }
  cleanTimer_.stop();

  if (settings_->writeMem()) {      
    sectors_[settings_->writememsect()]->writeLinkNewEvent(eventnum_);
  }

  addStubTimer_.start();

  vector<int> layerstubs(N_LAYER+N_DISK,0);
  vector<int> layerstubssector(N_SECTOR*(N_LAYER+N_DISK),0);
  
  for (int j = 0; j < ev.nstubs(); j++) {
    L1TStub stub = ev.stub(j);

    string dtc=stub.DTClink();
    unsigned int isector = stub.region();

    layerstubs[stub.layerdisk()]++;
    layerstubssector[isector*(N_LAYER+N_DISK)+stub.layerdisk()]++;
    
    sectors_[isector]->addStub(stub, dtc);
  }

  if (settings_->writeMonitorData("StubsLayerSector")) {
    for (unsigned int index=0;index<layerstubssector.size();index++) {
      int layerdisk=index%(N_LAYER+N_DISK);
      int sector=index/(N_LAYER+N_DISK);
      globals_->ofstream("stubslayersector.txt") << layerdisk << " " << sector << " " << layerstubssector[index] << endl;
    }
  }

  if (settings_->writeMonitorData("StubsLayer")) {
    for (unsigned int layerdisk=0;layerdisk<layerstubs.size();layerdisk++) {
      globals_->ofstream("stubslayer.txt") << layerdisk << " " << layerstubs[layerdisk] << endl;
    }
  }

  addStubTimer_.stop();

  // ----------------------------------------------------------------------------------------
  // Now start the tracklet processing

  // VM router
  VMRouterTimer_.start();
  for (unsigned int k = 0; k < N_SECTOR; k++) {
    sectors_[k]->executeVMR();
    if (settings_->writeMem() && k == settings_->writememsect()) {
      sectors_[k]->writeInputStubs(first);
      sectors_[k]->writeVMSTE(first);
      sectors_[k]->writeVMSME(first);
      sectors_[k]->writeAS(first);
    }
  }
  VMRouterTimer_.stop();

  // tracklet engine
  TETimer_.start();
  for (unsigned int k = 0; k < N_SECTOR; k++) {
    sectors_[k]->executeTE();
  }
  TETimer_.stop();

  // tracklet engine displaced
  TEDTimer_.start();
  for (unsigned int k = 0; k < N_SECTOR; k++) {
    sectors_[k]->executeTED();
  }
  TEDTimer_.stop();

  // triplet engine
  TRETimer_.start();
  for (unsigned int k = 0; k < N_SECTOR; k++) {
    sectors_[k]->executeTRE();
    if (settings_->writeMem() && k == settings_->writememsect()) {
      sectors_[k]->writeST(first);
    }
  }
  TRETimer_.stop();

  // tracklet processor (alternative implementation to TE+TC)
  TPTimer_.start();
  for (unsigned int k = 0; k < N_SECTOR; k++) {
    sectors_[k]->executeTP();
  }
  TPTimer_.stop();

  for (unsigned int k = 0; k < N_SECTOR; k++) {
    if (settings_->writeMem() && k == settings_->writememsect()) {
      sectors_[k]->writeSP(first);
    }
  }

  // tracklet calculator
  TCTimer_.start();
  for (unsigned int k = 0; k < N_SECTOR; k++) {
    sectors_[k]->executeTC();
  }
  TCTimer_.stop();

  int nTP = globals_->event()->nsimtracks();
  for (int iTP = 0; iTP < nTP; iTP++) {
    L1SimTrack simtrk = globals_->event()->simtrack(iTP);
    if (simtrk.pt() < 2.0)
      continue;
    if (std::abs(simtrk.vz()) > 15.0)
      continue;
    if (hypot(simtrk.vx(), simtrk.vy()) > 0.1)
      continue;
    bool electron = (abs(simtrk.type()) == 11);
    bool muon = (abs(simtrk.type()) == 13);
    bool pion = (abs(simtrk.type()) == 211);
    bool kaon = (abs(simtrk.type()) == 321);
    bool proton = (abs(simtrk.type()) == 2212);
    if (!(electron || muon || pion || kaon || proton))
      continue;
    int nlayers = 0;
    int ndisks = 0;
    int simtrackid = simtrk.trackid();
    unsigned int hitmask = ev.layersHit(simtrackid, nlayers, ndisks);
    if (nlayers + ndisks < 4)
      continue;

    if (settings_->writeMonitorData("HitEff")) {
      static ofstream outhit("hiteff.txt");
      outhit << simtrk.eta() << " " << (hitmask & 1) << " " << (hitmask & 2) << " " << (hitmask & 4) << " "
             << (hitmask & 8) << " " << (hitmask & 16) << " " << (hitmask & 32) << " " << (hitmask & 64) << " "
             << (hitmask & 128) << " " << (hitmask & 256) << " " << (hitmask & 512) << " " << (hitmask & 1024) << endl;
    }

    std::unordered_set<int> matchseed;
    for (unsigned int k = 0; k < N_SECTOR; k++) {
      std::unordered_set<int> matchseedtmp = sectors_[k]->seedMatch(iTP);
      matchseed.insert(matchseedtmp.begin(), matchseedtmp.end());
    }
    if (settings_->bookHistos()) {
      for (int iseed = 0; iseed < 8; iseed++) {
        bool eff = matchseed.find(iseed) != matchseed.end();
        globals_->histograms()->fillSeedEff(iseed, simtrk.eta(), eff);
      }
    }
  }

  // tracklet calculator displaced
  TCDTimer_.start();
  for (unsigned int k = 0; k < N_SECTOR; k++) {
    sectors_[k]->executeTCD();
  }
  TCDTimer_.stop();

  for (unsigned int k = 0; k < N_SECTOR; k++) {
    if (settings_->writeMem() && k == settings_->writememsect()) {
      sectors_[k]->writeTPAR(first);
      sectors_[k]->writeTPROJ(first);
    }
  }

  // projection router
  PRTimer_.start();
  for (unsigned int k = 0; k < N_SECTOR; k++) {
    sectors_[k]->executePR();
    if (settings_->writeMem() && k == settings_->writememsect()) {
      sectors_[k]->writeVMPROJ(first);
      sectors_[k]->writeAP(first);
    }
  }
  PRTimer_.stop();

  // match engine
  METimer_.start();
  for (unsigned int k = 0; k < N_SECTOR; k++) {
    sectors_[k]->executeME();
    if (settings_->writeMem() && k == settings_->writememsect()) {
      sectors_[k]->writeCM(first);
    }
  }
  METimer_.stop();

  // match calculator
  MCTimer_.start();
  for (unsigned int k = 0; k < N_SECTOR; k++) {
    sectors_[k]->executeMC();
  }
  MCTimer_.stop();

  // match processor (alternative to ME+MC)
  MPTimer_.start();
  for (unsigned int k = 0; k < N_SECTOR; k++) {
    sectors_[k]->executeMP();
  }
  MPTimer_.stop();

  for (unsigned int k = 0; k < N_SECTOR; k++) {
    if (settings_->writeMem() && k == settings_->writememsect()) {
      sectors_[k]->writeMC(first);
    }
  }

  // fit track
  FTTimer_.start();
  for (unsigned int k = 0; k < N_SECTOR; k++) {
    sectors_[k]->executeFT();
#ifndef USEHYBRID  //don't try to print these memories if running hybrid
    if ((settings_->writeMem() || settings_->writeMonitorData("IFit")) && k == settings_->writememsect()) {
      sectors_[k]->writeTF(first);
    }
#endif
  }
  FTTimer_.stop();

  // purge duplicate
  PDTimer_.start();
  for (unsigned int k = 0; k < N_SECTOR; k++) {
    sectors_[k]->executePD(tracks_);
#ifndef USEHYBRID  //don't try to print these memories if running hybrid
    if (((settings_->writeMem() || settings_->writeMonitorData("IFit")) && k == settings_->writememsect()) ||
        settings_->writeMonitorData("CT")) {
      sectors_[k]->writeCT(first);
    }
#endif
  }
  PDTimer_.stop();

}

void TrackletEventProcessor::printSummary() {
  //if (settings_->writeMonitorData("Cabling")) {
  //  cabling_->writephirange();
  //}

  if (settings_->bookHistos()) {
    globals_->histograms()->close();
  }

  edm::LogVerbatim("Tracklet") << "Process             Times called   Average time (ms)      Total time (s)";
  edm::LogVerbatim("Tracklet") << "Cleaning              " << setw(10) << cleanTimer_.ntimes() << setw(20)
                               << setprecision(3) << cleanTimer_.avgtime() * 1000.0 << setw(20) << setprecision(3)
                               << cleanTimer_.tottime();
  edm::LogVerbatim("Tracklet") << "Add Stubs             " << setw(10) << addStubTimer_.ntimes() << setw(20)
                               << setprecision(3) << addStubTimer_.avgtime() * 1000.0 << setw(20) << setprecision(3)
                               << addStubTimer_.tottime();
  edm::LogVerbatim("Tracklet") << "VMRouter              " << setw(10) << VMRouterTimer_.ntimes() << setw(20)
                               << setprecision(3) << VMRouterTimer_.avgtime() * 1000.0 << setw(20) << setprecision(3)
                               << VMRouterTimer_.tottime();
  if (settings_->combined()) {
    edm::LogVerbatim("Tracklet") << "TrackletProcessor     " << setw(10) << TPTimer_.ntimes() << setw(20)
                                 << setprecision(3) << TPTimer_.avgtime() * 1000.0 << setw(20) << setprecision(3)
                                 << TPTimer_.tottime();
    edm::LogVerbatim("Tracklet") << "MatchProcessor        " << setw(10) << MPTimer_.ntimes() << setw(20)
                                 << setprecision(3) << MPTimer_.avgtime() * 1000.0 << setw(20) << setprecision(3)
                                 << MPTimer_.tottime();
  } else {
    edm::LogVerbatim("Tracklet") << "TrackletEngine        " << setw(10) << TETimer_.ntimes() << setw(20)
                                 << setprecision(3) << TETimer_.avgtime() * 1000.0 << setw(20) << setprecision(3)
                                 << TETimer_.tottime();
    edm::LogVerbatim("Tracklet") << "TrackletEngineDisplaced" << setw(10) << TEDTimer_.ntimes() << setw(20)
                                 << setprecision(3) << TEDTimer_.avgtime() * 1000.0 << setw(20) << setprecision(3)
                                 << TEDTimer_.tottime();
    edm::LogVerbatim("Tracklet") << "TripletEngine         " << setw(10) << TRETimer_.ntimes() << setw(20)
                                 << setprecision(3) << TRETimer_.avgtime() * 1000.0 << setw(20) << setprecision(3)
                                 << TRETimer_.tottime();
    edm::LogVerbatim("Tracklet") << "TrackletCalculator    " << setw(10) << TCTimer_.ntimes() << setw(20)
                                 << setprecision(3) << TCTimer_.avgtime() * 1000.0 << setw(20) << setprecision(3)
                                 << TCTimer_.tottime();
    edm::LogVerbatim("Tracklet") << "TrackletCalculatorDisplaced" << setw(10) << TCDTimer_.ntimes() << setw(20)
                                 << setprecision(3) << TCDTimer_.avgtime() * 1000.0 << setw(20) << setprecision(3)
                                 << TCDTimer_.tottime();
    edm::LogVerbatim("Tracklet") << "ProjectionRouter      " << setw(10) << PRTimer_.ntimes() << setw(20)
                                 << setprecision(3) << PRTimer_.avgtime() * 1000.0 << setw(20) << setprecision(3)
                                 << PRTimer_.tottime();
    edm::LogVerbatim("Tracklet") << "MatchEngine           " << setw(10) << METimer_.ntimes() << setw(20)
                                 << setprecision(3) << METimer_.avgtime() * 1000.0 << setw(20) << setprecision(3)
                                 << METimer_.tottime();
    edm::LogVerbatim("Tracklet") << "MatchCalculator       " << setw(10) << MCTimer_.ntimes() << setw(20)
                                 << setprecision(3) << MCTimer_.avgtime() * 1000.0 << setw(20) << setprecision(3)
                                 << MCTimer_.tottime();
  }
  edm::LogVerbatim("Tracklet") << "FitTrack              " << setw(10) << FTTimer_.ntimes() << setw(20)
                               << setprecision(3) << FTTimer_.avgtime() * 1000.0 << setw(20) << setprecision(3)
                               << FTTimer_.tottime();
  edm::LogVerbatim("Tracklet") << "PurgeDuplicate        " << setw(10) << PDTimer_.ntimes() << setw(20)
                               << setprecision(3) << PDTimer_.avgtime() * 1000.0 << setw(20) << setprecision(3)
                               << PDTimer_.tottime();
}

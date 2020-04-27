// ROOT includes
#include "TMath.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TLatex.h>

#include <iostream>

#include "../interface/IMATH_TrackletCalculator.h"
#include "../interface/IMATH_TrackletCalculatorDisk.h"
#include "../interface/IMATH_TrackletCalculatorOverlap.h"

#include "../interface/slhcevent.h"

#include "../interface/Sector.h"
#include "../interface/Cabling.h"
#include "../interface/FPGAWord.h"
#include "../interface/CPUTimer.h"
#include "../interface/StubVariance.h"
#include "../interface/Settings.h"

#include "../interface/GlobalHistTruth.h"
#include "../interface/HistImp.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#ifdef IMATH_ROOT
TFile* var_base::h_file_=0;
bool   var_base::use_root = false;
#endif

#include <iomanip>
#include <fstream>
#include <string>

//Uncomment if you want root output
//#define USEROOT


// Include file to define ROOT-Tree
// --------------------------------
#ifdef USEROOT
#include "FPGAEvent.h"
#endif
// --------------------------------

int main(const int argc, const char** argv)
{

  Trklet::Settings settings;

  // tracklet calculators 
  IMATH_TrackletCalculator* ITC_L1L2 = new IMATH_TrackletCalculator(&settings,1,2);
  IMATH_TrackletCalculator* ITC_L2L3 = new IMATH_TrackletCalculator(&settings,2,3);
  IMATH_TrackletCalculator* ITC_L3L4 = new IMATH_TrackletCalculator(&settings,3,4);
  IMATH_TrackletCalculator* ITC_L5L6 = new IMATH_TrackletCalculator(&settings,5,6);
  
  IMATH_TrackletCalculatorDisk* ITC_F1F2 = new IMATH_TrackletCalculatorDisk(&settings,1,2);
  IMATH_TrackletCalculatorDisk* ITC_F3F4 = new IMATH_TrackletCalculatorDisk(&settings,3,4);
  IMATH_TrackletCalculatorDisk* ITC_B1B2 = new IMATH_TrackletCalculatorDisk(&settings,-1,-2);
  IMATH_TrackletCalculatorDisk* ITC_B3B4 = new IMATH_TrackletCalculatorDisk(&settings,-3,-4);
  
  IMATH_TrackletCalculatorOverlap* ITC_L1F1 = new IMATH_TrackletCalculatorOverlap(&settings,1,1);
  IMATH_TrackletCalculatorOverlap* ITC_L2F1 = new IMATH_TrackletCalculatorOverlap(&settings,2,1);
  IMATH_TrackletCalculatorOverlap* ITC_L1B1 = new IMATH_TrackletCalculatorOverlap(&settings,1,-1);
  IMATH_TrackletCalculatorOverlap* ITC_L2B1 = new IMATH_TrackletCalculatorOverlap(&settings,2,-1);

  GlobalHistTruth::ITC_L1L2()=ITC_L1L2;
  GlobalHistTruth::ITC_L2L3()=ITC_L2L3;
  GlobalHistTruth::ITC_L3L4()=ITC_L3L4;
  GlobalHistTruth::ITC_L5L6()=ITC_L5L6;

  GlobalHistTruth::ITC_F1F2()=ITC_F1F2;
  GlobalHistTruth::ITC_F3F4()=ITC_F3F4;
  GlobalHistTruth::ITC_B1B2()=ITC_B1B2;
  GlobalHistTruth::ITC_B3B4()=ITC_B3B4;

  GlobalHistTruth::ITC_L1F1()=ITC_L1F1;
  GlobalHistTruth::ITC_L2F1()=ITC_L2F1;
  GlobalHistTruth::ITC_L1B1()=ITC_L1B1;
  GlobalHistTruth::ITC_L2B1()=ITC_L2B1;

  settings.krinvpars() = ITC_L1L2->rinv_final.get_K();
  settings.kphi0pars() = ITC_L1L2->phi0_final.get_K();
  settings.kd0pars()   = settings.kd0();
  settings.ktpars()    = ITC_L1L2->t_final.get_K();
  settings.kz0pars()   = ITC_L1L2->z0_final.get_K();
  settings.kphiproj123() =ITC_L1L2->phi0_final.get_K()*4;
  settings.kzproj()=settings.kz();
  settings.kphider()=ITC_L1L2->rinv_final.get_K()*(1<<settings.phiderbitshift());
  settings.kzder()=ITC_L1L2->t_final.get_K()*(1<<settings.zderbitshift());
  settings.krprojshiftdisk() = ITC_L1L2->rD_0_final.get_K();
  settings.kphiprojdisk()=ITC_L1L2->phi0_final.get_K()*4.0;
  settings.krdisk() = settings.kr();
  settings.kzpars() = settings.kz();  


  edm::LogVerbatim("Tracklet") << "=========================================================";
  edm::LogVerbatim("Tracklet") << "Conversion factors for global coordinates:";
  edm::LogVerbatim("Tracklet") << "z    kz            = "<< settings.kz() ;
  edm::LogVerbatim("Tracklet") << "r    kr            = "<< settings.kr() ;
  edm::LogVerbatim("Tracklet") << "phi  kphi1         = "<< settings.kphi1() ;
  edm::LogVerbatim("Tracklet") << "=========================================================";
  edm::LogVerbatim("Tracklet") << "Conversion factors for track(let) parameters:";
  edm::LogVerbatim("Tracklet") << "rinv krinvpars     = "<< settings.krinvpars() ;
  edm::LogVerbatim("Tracklet") << "phi0 kphi0pars     = "<< settings.kphi0pars() ;
  edm::LogVerbatim("Tracklet") << "d0   kd0pars       = "<< settings.kd0pars() ;
  edm::LogVerbatim("Tracklet") << "t    ktpars        = "<< settings.ktpars() ;
  edm::LogVerbatim("Tracklet") << "z0   kz0pars       = "<< settings.kzpars() ;
  edm::LogVerbatim("Tracklet") << "=========================================================";
  edm::LogVerbatim("Tracklet") << "phi0bitshift = "<<settings.phi0bitshift();
  edm::LogVerbatim("Tracklet") << "d0bitshift   = "<<"???";
  edm::LogVerbatim("Tracklet");
  edm::LogVerbatim("Tracklet") << "=========================================================";

#include "../plugins/WriteInvTables.icc"
#include "../plugins/WriteDesign.icc"
  
  using namespace std;
  if (argc<4)
    edm::LogVerbatim("Tracklet") << "Need to specify the input ascii file and the number of events to run on and if you want to filter on MC truth";

  HistImp* histimp=new HistImp;
  histimp->init();
  histimp->bookLayerResidual();
  histimp->bookDiskResidual();
  histimp->bookTrackletParams();
  histimp->bookSeedEff();
  
  GlobalHistTruth::histograms()=histimp;
  
  int nevents = atoi(argv[2]);

  int selectmu = atoi(argv[3]);

  assert((selectmu==0)||(selectmu==1));

  ifstream infile;
  istream* in = &cin;
  if(strcmp(argv[1],"stdin")){
    infile.open(argv[1]);
    in = &infile;
  }

  ofstream outres;
  if (settings.writeMonitorData("ResEff")) outres.open("trackres.txt");

  ofstream outeff;
  if (settings.writeMonitorData("ResEff")) outeff.open("trackeff.txt");

  ofstream outpars;
  if (settings.writeMonitorData("Pars")) outpars.open("trackpars.txt");


//Open file to hold ROOT-Tree
// --------------------------
#ifdef USEROOT
  TFile  *hfile = new TFile("myTest.root","RECREATE","Simple ROOT Ntuple"); 
  TTree *trackTree = new TTree("FPGAEvent","L1Track Tree");
  FPGAEvent *fpgaEvent = new FPGAEvent;  
  fpgaEvent->reset();
  trackTree->Branch("Event",&fpgaEvent);
#endif
// --------------------------



// Define Sectors (boards)	 
  Sector** sectors=new Sector*[settings.NSector()];

  Cabling cabling;

  cabling.init("../data/calcNumDTCLinks.txt","../data/modules_T5v3_27SP_nonant_tracklet.dat");


  
  for (unsigned int i=0;i<settings.NSector();i++) {
    sectors[i]=new Sector(i,&settings);
  }  


  edm::LogVerbatim("Tracklet") << "Will read memory modules file";

  string memfile="../data/memorymodules_"+settings.geomext()+".dat";
  ifstream inmem(memfile.c_str());
  assert(inmem.good());

  while (inmem.good()){
    string memType, memName, size;
    inmem >>memType>>memName>>size;
    if (!inmem.good()) continue;
    if (settings.writetrace()) {
      edm::LogVerbatim("Tracklet") << "Read memory: "<<memType<<" "<<memName;
    }
    for (unsigned int i=0;i<settings.NSector();i++) {
      sectors[i]->addMem(memType,memName);
    }
    
  }


  edm::LogVerbatim("Tracklet") << "Will read processing modules file";

  string procfile="../data/processingmodules_"+settings.geomext()+".dat";
  ifstream inproc(procfile.c_str());
  assert(inproc.good());

  while (inproc.good()){
    string procType, procName;
    inproc >>procType>>procName;
    if (!inproc.good()) continue;
    if (settings.writetrace()) {
      edm::LogVerbatim("Tracklet") << "Read process: "<<procType<<" "<<procName;
    }
    for (unsigned int i=0;i<settings.NSector();i++) {
      sectors[i]->addProc(procType,procName);
    }
    
  }


  edm::LogVerbatim("Tracklet") << "Will read wiring information";

  string wirefile="../data/wires_"+settings.geomext()+".dat";
  ifstream inwire(wirefile.c_str());
  assert(inwire.good());


  while (inwire.good()){
    string line;
    getline(inwire,line);
    if (!inwire.good()) continue;
    if (settings.writetrace()) {
      edm::LogVerbatim("Tracklet") << "Line : "<<line;
    }
    stringstream ss(line);
    string mem,tmp1,procin,tmp2,procout;
    ss>>mem>>tmp1>>procin;
    if (procin=="output=>") {
      procin="";
      ss>>procout;
    }
    else{
      ss>>tmp2>>procout;
    }

    for (unsigned int i=0;i<settings.NSector();i++) {
      sectors[i]->addWire(mem,procin,procout);
    }
  
  }

  std::map<string,vector<int> > dtclayerdisk;

  ifstream indtc("../data/dtclinklayerdisk.dat");
  assert(indtc.good());
  string dtc;
  indtc >> dtc;
  while (indtc.good()){
    vector<int> tmp;
    dtclayerdisk[dtc]=tmp;
    int layerdisk;
    indtc >> layerdisk;
    while (layerdisk>0) {
      dtclayerdisk[dtc].push_back(layerdisk);
      indtc >> layerdisk;
    }
    indtc >> dtc;
  }
  
  ofstream skimout;
  if (settings.skimfile() != "") skimout.open(settings.skimfile().c_str());

  CPUTimer readTimer;
  CPUTimer cleanTimer;
  CPUTimer addStubTimer;
  CPUTimer VMRouterTimer;  
  CPUTimer TETimer;
  CPUTimer TEDTimer;
  CPUTimer TRETimer;
  CPUTimer TCTimer;
  CPUTimer TCDTimer;
  CPUTimer PRTimer;
  CPUTimer METimer;
  CPUTimer MCTimer;
  CPUTimer MPTimer;
  CPUTimer FTTimer;
  CPUTimer PDTimer;

  if (settings.writeMonitorData("Seeds")) {
    ofstream fout("seeds.txt", ofstream::out);
    fout.close();
  }

  bool first=true;

  for (int eventnum=0;eventnum<nevents&&!in->eof();eventnum++){
    
    readTimer.start();
    SLHCEvent ev(*in);
    readTimer.stop();

    GlobalHistTruth::event()=&ev;

    
    L1SimTrack simtrk;


// setup ROOT Tree and Add Monte Carlo tracks to the ROOT-Tree Event
// -----------------------------------------------------------------
#ifdef USEROOT
    fpgaEvent->reset();
    fpgaEvent->nevt = eventnum;
    for(int nst=0; nst<ev.nsimtracks(); nst++) {
      simtrk = ev.simtrack(nst);
      FPGAEventMCTrack *mcTrack = new FPGAEventMCTrack(simtrk.type(),simtrk.pt(),simtrk.eta(),simtrk.phi(),simtrk.vx(),simtrk.vy(),simtrk.vz());
      fpgaEvent->mcTracks.push_back(*mcTrack);
    }
#endif
// ------------------------------------------------------------------	 

    
    if (selectmu==1) {

      if (ev.nsimtracks()==0) {
      	eventnum--;
	continue;
      }

      
      simtrk=ev.simtrack(0);

      if (settings.debugTracklet()) {
	edm::LogVerbatim("Tracklet") <<"nstub simtrkid pt phi eta t vz:"<<ev.nstubs()<<" "<<simtrk.trackid()<<" "<<simtrk.pt()<<" "<<simtrk.phi()<<" "
	     <<simtrk.eta()<<" "
	     <<sinh(simtrk.eta())<<" "
	     <<simtrk.vz();
      }
      
      bool good=(fabs(simtrk.pt())>2.0
		 //&&fabs(simtrk.pt())<3.0
                 &&fabs(simtrk.vz())<15.0
		 //&&simtrk.eta()>0.0
		 //&&fabs(fabs(simtrk.eta())-0.0)>1.6
		 //&&fabs(fabs(simtrk.eta())-0.0)<1.9
		 &&fabs(fabs(simtrk.eta())-0.0)<2.4
		 //&&fabs(fabs(simtrk.eta())-0.0)>1.8
		 //&&fabs(fabs(simtrk.phi())-0.05)<0.05
		 //&&fabs(simtrk.eta()-0.0)>1.6
		 //&&fabs(simtrk.eta()-1.5)<0.2
		 //&&fabs(phisector-0.61)<0.03
		 //&&fabs(fabs(simtrk.eta())-1.4)<0.1
                 //&&fabs(simtrk.d0())<1.0
		 );


      if (!good) {
	eventnum--;
	continue;
      }

      if (settings.skimfile() != "") ev.write(skimout);
 
    } 

    if (settings.writeMonitorData("Seeds")) {
      ofstream fout("seeds.txt", ofstream::app);
      fout << "======== Event " << eventnum << " ========" << endl;
      for(unsigned nst=0; nst<ev.nsimtracks(); nst++) {
        const L1SimTrack &simtrk = ev.simtrack(nst);
        fout << "SimTrk " << simtrk.pt() << " " << simtrk.eta() << " " << simtrk.phi() << " " << simtrk.d0() << " ";

        vector<string> hitPattern;
        for(int i=0; i<ev.nstubs(); i++) {
          const L1TStub stub = ev.stub(i);
          if (!stub.tpmatch(simtrk.trackid()))
            continue;
          if (stub.layer() < 999) {
            switch (stub.layer()) {
              case 0: hitPattern.push_back("L1"); break;
              case 1: hitPattern.push_back("L2"); break;
              case 2: hitPattern.push_back("L3"); break;
              case 3: hitPattern.push_back("L4"); break;
              case 4: hitPattern.push_back("L5"); break;
              case 5: hitPattern.push_back("L6"); break;
              default: edm::LogVerbatim("Tracklet") << "Stub layer: " << stub.layer(); assert(0);
            }
          }
          else {
            string d = (stub.isPSmodule() ? "D" : "d");
            switch (abs(stub.disk())) {
              case 1: hitPattern.push_back(d+"1"); break;
              case 2: hitPattern.push_back(d+"2"); break;
              case 3: hitPattern.push_back(d+"3"); break;
              case 4: hitPattern.push_back(d+"4"); break;
              case 5: hitPattern.push_back(d+"5"); break;
              default: edm::LogVerbatim("Tracklet") << "Stub disk: " << stub.disk(); assert(0);
            }
          }
        }
        bool (*compare)(const string &, const string &) = [](const string &a, const string &b) -> bool {
          if (a.at(0) == 'L' && b.at(0) == 'D')
            return true;
          else if (a.at(0) == 'D' && b.at(0) == 'L')
            return false;
          else
            return a.at(1) < b.at(1);
        };
        sort(hitPattern.begin(), hitPattern.end(), compare);
        hitPattern.erase(unique(hitPattern.begin(), hitPattern.end()), hitPattern.end());
        for (const auto &stub : hitPattern)
          fout << stub;
        if (hitPattern.empty())
          fout << "XX";
        fout << endl;
      }
      fout.close();
    }

    if (settings.writeMonitorData("Variance")) {
      StubVariance variance(ev);
    }

    edm::LogVerbatim("Tracklet") <<"Process event: "<<eventnum<<" with "<<ev.nstubs()<<" stubs and "<<ev.nsimtracks()<<" simtracks";

    std::vector<Track*> tracks;

    int nlayershit=0;


    
// Processesing done in FPGA.icc
#include "../plugins/FPGA.icc"

// Block for producing ROOT-Tree
// ------------------------------
#ifdef USEROOT
#include "FPGATree.icc"
#endif
// ------------------------------

    if (settings.writeMonitorData("ResEff")) {
      outeff << simtrk.pt()*simtrk.trackid()/fabs(simtrk.trackid())<<" "<<simtrk.eta()
	     <<" "<<simtrk.phi();
      if (match) outeff << " 1"<<endl;
       else outeff << " 0"<<endl;
    }

    if (settings.writeMonitorData("MatchEff")) { 
      static ofstream out("matcheff.txt");
      int nsim=0;
      for(unsigned int isimtrack=0;isimtrack<ev.nsimtracks();isimtrack++){
        L1SimTrack simtrack=ev.simtrack(isimtrack);
        if (simtrack.pt()<2.0) continue;
        if (fabs(simtrack.eta())>2.4) continue;
        if (fabs(simtrack.vz())>15.0) continue;
        if (hypot(simtrack.vx(),simtrack.vy())>0.1) continue;
        bool electron=(abs(simtrack.type())==11);
        bool muon=(abs(simtrack.type())==13);
        bool pion=(abs(simtrack.type())==211);
        bool kaon=(abs(simtrack.type())==321);
        bool proton=(abs(simtrack.type())==2212);
        if (!(electron||muon||pion||kaon||proton)) continue;
        int nlayers=0;
        int ndisks=0;
        int simeventid=simtrack.eventid();
        int simtrackid=simtrack.trackid();
        ev.layersHit(simtrackid,nlayers,ndisks);
        if (nlayers+ndisks<4) continue;
	nsim++;
	for (int seed=-1;seed<8;seed++){
	  bool eff=false;
	  bool effloose=false;
	  //int layerdisk=0;
	  int itrackmatch=-1;
	  for (unsigned int itrack=0;itrack<tracks.size();itrack++) {
	    std::vector<L1TStub*> stubs=tracks[itrack]->stubs();
	    if (seed==-1) {
	      if (tracks[itrack]->duplicate()) continue;
	    } else {
	      if (seed!=tracks[itrack]->seed()) continue;
	    }
	  
	    unsigned int nmatch=0;
	    //layerdisk=0;
	    for(unsigned int istub=0;istub<stubs.size();istub++){
	      if (stubs[istub]->tpmatch(simtrackid)) {
		nmatch++;
	      } else {
		if (stubs[istub]->layer()<999) {
		  //layerdisk=stubs[istub]->layer()+1;
		} else {
		  //layerdisk=-abs(stubs[istub]->disk());
		}
	      }
	    }

	    if (nmatch==stubs.size()) {
	      eff=true;
	      itrackmatch=itrack;
	    }
	    if (nmatch>=stubs.size()-1) {
	      effloose=true;
	      if (!eff) itrackmatch=itrack;
	    }

	  }
	  double dpt=-999;
	  double dphi=-999;
	  double deta=-999;
	  double dz0=-999;
	  int q=1;
	  if (simtrack.type()==11||simtrack.type()==13||
	      simtrack.type()==-211||simtrack.type()==-321||simtrack.type()==-2212){
	    q=-1;
	  }

	  if (itrackmatch>=0) {
	    dpt=tracks[itrackmatch]->pt(&settings)-q*simtrack.pt();
	    dphi=tracks[itrackmatch]->phi0(&settings)-simtrack.phi();
	    if (dphi>M_PI) dphi-=2*M_PI;
	    if (dphi<-M_PI) dphi+=2*M_PI;
	    deta=tracks[itrackmatch]->eta(&settings)-simtrack.eta();
	    dz0=tracks[itrackmatch]->z0(&settings)-simtrack.vz();
	  }
	
	  out <<eventnum<<" "<<simeventid<<" "<<seed<<" "<<simtrackid<<" "<<simtrack.type()<<" "
	      <<simtrack.pt()<<" "<<simtrack.eta()<<" "<<simtrack.phi()<<" "
	      <<simtrack.vx()<<" "<<simtrack.vy()<<" "<<simtrack.vz()<<" "
	      <<eff<<" "<<effloose<<" "
	      <<dpt<<" "<<dphi<<" "<<deta<<" "<<dz0
	      <<endl;
	}
      }
    }

    
    
    // Clean up 

    //edm::LogVerbatim("Tracklet") << "Duplicates : ";
    int ntrack=0;
    for(unsigned int l=0;l<tracks.size();l++) {
      //  edm::LogVerbatim("Tracklet") <<tracks[l].duplicate()<<" ";
      if (settings.writeMonitorData("Pars")) {
	double phi=tracks[l]->iphi0()*settings.kphi0pars()+tracks[l]->sector()*2*M_PI/settings.NSector();
	if (phi>M_PI) phi-=2*M_PI;
	double phisec=phi-2*M_PI;
	while (phisec<0.0) phisec+=2*M_PI/settings.NSector();
	outpars  <<tracks[l]->duplicate()<<" "<<asinh(tracks[l]->it()*settings.ktpars())<<" "
		 <<phi<<" "<<tracks[l]->iz0()*settings.kz()<<" "<<phisec/(2*M_PI/settings.NSector())<<" "
		 <<tracks[l]->irinv()*settings.krinvpars();
      }   	
      if (!tracks[l]->duplicate()) {
	//edm::LogVerbatim("Tracklet") << "FPGA Track pt, eta, phi, z0, chi2 = " 
	//   << tracks[l]->pt() << " " << tracks[l]->eta() << " " << tracks[l]->phi0() << " " << tracks[l]->z0() << " " << tracks[l]->chisq() 
	//   << " seed " << tracks[l]->seed() << " duplicate " << tracks[l]->duplicate();
	//edm::LogVerbatim("Tracklet") << " ---------- not duplicate";
	//edm::LogVerbatim("Tracklet") << "tapprox "<<tracks[l].eta();
	ntrack++;
	//edm::LogVerbatim("Tracklet") << "eta = "<<tracks[l].eta();
      }
    }
    
    edm::LogVerbatim("Tracklet") << "Number layers/disks hit = "<<nlayershit<<" number of found tracks : "<<tracks.size()<<" unique "<<ntrack;


  // dump what was found   
//	 printf("Track Parameters: \n");
//	 for(std::vector<Track*>::iterator trk=tracks.begin(); trk!=tracks.end(); trk++){
//	   printf("irinv = %i \n", (*trk)->irinv() );
//		 printf("iphi0 = %i \n", (*trk)->iphi0() );
//		 printf("iz0   = %i \n", (*trk)->iz0() );
//		 printf("it    = %i \n", (*trk)->it() );
//		 printf("stubID=");
//		 std::map<int, int> stubs = (*trk)->stubID();
//		 for(std::map<int, int>::iterator sb=stubs.begin(); sb!=stubs.end(); sb++) printf(" %i -- %i ",sb->first,sb->second);
//		 printf("\n");
//		 printf("dup   = %i\n \n", (*trk)->duplicate());
//		 printf("chisq   = %f\n \n", (*trk)->chisq());
//	 } 




    first=false;

  }

  if (settings.writeMonitorData("Cabling")) {
    cabling.writephirange();
  }
  
  edm::LogVerbatim("Tracklet") << "Process             Times called   Average time (ms)      Total time (s)";
  edm::LogVerbatim("Tracklet") << "Reading               "
       <<setw(10)<<readTimer.ntimes()
       <<setw(20)<<setprecision(3)<<readTimer.avgtime()*1000.0
       <<setw(20)<<setprecision(3)<<readTimer.tottime();
  edm::LogVerbatim("Tracklet") << "Cleaning              "
       <<setw(10)<<cleanTimer.ntimes()
       <<setw(20)<<setprecision(3)<<cleanTimer.avgtime()*1000.0
       <<setw(20)<<setprecision(3)<<cleanTimer.tottime();
  edm::LogVerbatim("Tracklet") << "Add Stubs             "
       <<setw(10)<<addStubTimer.ntimes()
       <<setw(20)<<setprecision(3)<<addStubTimer.avgtime()*1000.0
       <<setw(20)<<setprecision(3)<<addStubTimer.tottime();
  edm::LogVerbatim("Tracklet") << "VMRouter              "
       <<setw(10)<<VMRouterTimer.ntimes()
       <<setw(20)<<setprecision(3)<<VMRouterTimer.avgtime()*1000.0
       <<setw(20)<<setprecision(3)<<VMRouterTimer.tottime();
  edm::LogVerbatim("Tracklet") << "TrackletEngine        "
       <<setw(10)<<TETimer.ntimes()
       <<setw(20)<<setprecision(3)<<TETimer.avgtime()*1000.0
       <<setw(20)<<setprecision(3)<<TETimer.tottime();
  edm::LogVerbatim("Tracklet") << "TrackletEngineDisplaced"
       <<setw(10)<<TEDTimer.ntimes()
       <<setw(20)<<setprecision(3)<<TEDTimer.avgtime()*1000.0
       <<setw(20)<<setprecision(3)<<TEDTimer.tottime();
  edm::LogVerbatim("Tracklet") << "TripletEngine         "
       <<setw(10)<<TRETimer.ntimes()
       <<setw(20)<<setprecision(3)<<TRETimer.avgtime()*1000.0
       <<setw(20)<<setprecision(3)<<TRETimer.tottime();
  edm::LogVerbatim("Tracklet") << "TrackletCalculator    "
       <<setw(10)<<TCTimer.ntimes()
       <<setw(20)<<setprecision(3)<<TCTimer.avgtime()*1000.0
       <<setw(20)<<setprecision(3)<<TCTimer.tottime();
  edm::LogVerbatim("Tracklet") << "TrackletCalculatorDisplaced"
       <<setw(10)<<TCDTimer.ntimes()
       <<setw(20)<<setprecision(3)<<TCDTimer.avgtime()*1000.0
       <<setw(20)<<setprecision(3)<<TCDTimer.tottime();
  edm::LogVerbatim("Tracklet") << "ProjectionRouter      "
       <<setw(10)<<PRTimer.ntimes()
       <<setw(20)<<setprecision(3)<<PRTimer.avgtime()*1000.0
       <<setw(20)<<setprecision(3)<<PRTimer.tottime();
  edm::LogVerbatim("Tracklet") << "MatchEngine           "
       <<setw(10)<<METimer.ntimes()
       <<setw(20)<<setprecision(3)<<METimer.avgtime()*1000.0
       <<setw(20)<<setprecision(3)<<METimer.tottime();
  edm::LogVerbatim("Tracklet") << "MatchCalculator       "
       <<setw(10)<<MCTimer.ntimes()
       <<setw(20)<<setprecision(3)<<MCTimer.avgtime()*1000.0
       <<setw(20)<<setprecision(3)<<MCTimer.tottime();
  edm::LogVerbatim("Tracklet") << "MatchProcessor        "
       <<setw(10)<<MPTimer.ntimes()
       <<setw(20)<<setprecision(3)<<MPTimer.avgtime()*1000.0
       <<setw(20)<<setprecision(3)<<MPTimer.tottime();
  edm::LogVerbatim("Tracklet") << "FitTrack              "
       <<setw(10)<<FTTimer.ntimes()
       <<setw(20)<<setprecision(3)<<FTTimer.avgtime()*1000.0
       <<setw(20)<<setprecision(3)<<FTTimer.tottime();
  edm::LogVerbatim("Tracklet") << "PurgeDuplicate        "
       <<setw(10)<<PDTimer.ntimes()
       <<setw(20)<<setprecision(3)<<PDTimer.avgtime()*1000.0
       <<setw(20)<<setprecision(3)<<PDTimer.tottime();


  if (settings.skimfile()!="") skimout.close();

  histimp->close();
  
// Write and Close ROOT-Tree  
// -------------------------
#ifdef USEROOT
   hfile->Write();
   hfile->Close();
#endif
// --------------------------  

  for (unsigned int i=0;i<settings.NSector();i++)
    delete sectors[i];

} 

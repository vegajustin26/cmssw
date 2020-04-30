#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletEventProcessor_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletEventProcessor_h

#include "../interface/GlobalHistTruth.h"
#include "../interface/HistImp.h"

namespace Trklet {
  class TrackletEventProcessor {

  public:
    
    TrackletEventProcessor() {}

    ~TrackletEventProcessor() {
      for (unsigned int i=0;i<settings_->NSector();i++)
	delete sectors_[i];

      delete globals_;
      
    }
    
    void init(const Settings* theSettings) {
      settings_=theSettings;

      
      globals_ = new GlobalHistTruth(settings_);
      
      settings_->krinvpars() = globals_->ITC_L1L2()->rinv_final.get_K();
      settings_->kphi0pars() = globals_->ITC_L1L2()->phi0_final.get_K();
      settings_->kd0pars()   = settings_->kd0();
      settings_->ktpars()    = globals_->ITC_L1L2()->t_final.get_K();
      settings_->kz0pars()   = globals_->ITC_L1L2()->z0_final.get_K();
      settings_->kphiproj123() =globals_->ITC_L1L2()->phi0_final.get_K()*4;
      settings_->kzproj()=settings_->kz();
      settings_->kphider()=globals_->ITC_L1L2()->rinv_final.get_K()*(1<<settings_->phiderbitshift());
      settings_->kzder()=globals_->ITC_L1L2()->t_final.get_K()*(1<<settings_->zderbitshift());
      settings_->krprojshiftdisk() = globals_->ITC_L1L2()->rD_0_final.get_K();
      settings_->kphiprojdisk()=globals_->ITC_L1L2()->phi0_final.get_K()*4.0;
      settings_->krdisk() = settings_->kr();
      settings_->kzpars() = settings_->kz();  

      edm::LogVerbatim("Tracklet") << "=========================================================";
      edm::LogVerbatim("Tracklet") << "Conversion factors for global coordinates:";
      edm::LogVerbatim("Tracklet") << "z    kz            = "<< settings_->kz() ;
      edm::LogVerbatim("Tracklet") << "r    kr            = "<< settings_->kr() ;
      edm::LogVerbatim("Tracklet") << "phi  kphi1         = "<< settings_->kphi1() ;
      edm::LogVerbatim("Tracklet") << "=========================================================";
      edm::LogVerbatim("Tracklet") << "Conversion factors for track(let) parameters:";
      edm::LogVerbatim("Tracklet") << "rinv krinvpars     = "<< settings_->krinvpars() ;
      edm::LogVerbatim("Tracklet") << "phi0 kphi0pars     = "<< settings_->kphi0pars() ;
      edm::LogVerbatim("Tracklet") << "d0   kd0pars       = "<< settings_->kd0pars() ;
      edm::LogVerbatim("Tracklet") << "t    ktpars        = "<< settings_->ktpars() ;
      edm::LogVerbatim("Tracklet") << "z0   kz0pars       = "<< settings_->kzpars() ;
      edm::LogVerbatim("Tracklet") << "=========================================================";
      edm::LogVerbatim("Tracklet") << "phi0bitshift = "<<settings_->phi0bitshift();
      edm::LogVerbatim("Tracklet") << "d0bitshift   = "<<"???";
      edm::LogVerbatim("Tracklet") << "=========================================================";

      const Settings& settings=*settings_;
      GlobalHistTruth* globals=globals_;
#include "../plugins/WriteInvTables.icc"
#include "../plugins/WriteDesign.icc"


      HistImp* histimp;
      if (settings_->bookHistos()) {
	histimp = new HistImp;
	histimp->init();
	histimp->bookLayerResidual();
	histimp->bookDiskResidual();
	histimp->bookTrackletParams();
	histimp->bookSeedEff();
  
	globals_->histograms()=histimp;
      }
  

      
      // Crate the sectors (boards)	 
      sectors_=new Sector*[settings_->NSector()];

      for (unsigned int i=0;i<settings_->NSector();i++) {
	sectors_[i]=new Sector(i,settings_,globals_);
      }  


      edm::LogVerbatim("Tracklet") << "Will read memory modules file";

      string memfile="../data/memorymodules_"+settings_->geomext()+".dat";
      ifstream inmem(memfile.c_str());
      assert(inmem.good());

      while (inmem.good()){
	string memType, memName, size;
	inmem >>memType>>memName>>size;
	if (!inmem.good()) continue;
	if (settings_->writetrace()) {
	  edm::LogVerbatim("Tracklet") << "Read memory: "<<memType<<" "<<memName;
	}
	for (unsigned int i=0;i<settings_->NSector();i++) {
	  sectors_[i]->addMem(memType,memName);
	}
	
      }
      
      
      edm::LogVerbatim("Tracklet") << "Will read processing modules file";
      
      string procfile="../data/processingmodules_"+settings_->geomext()+".dat";
      ifstream inproc(procfile.c_str());
      assert(inproc.good());
      
      while (inproc.good()){
	string procType, procName;
	inproc >>procType>>procName;
	if (!inproc.good()) continue;
	if (settings_->writetrace()) {
	  edm::LogVerbatim("Tracklet") << "Read process: "<<procType<<" "<<procName;
	}
	for (unsigned int i=0;i<settings_->NSector();i++) {
	  sectors_[i]->addProc(procType,procName);
	}
	
      }
      
      
      edm::LogVerbatim("Tracklet") << "Will read wiring information";
      
      string wirefile="../data/wires_"+settings_->geomext()+".dat";
      ifstream inwire(wirefile.c_str());
      assert(inwire.good());
      
      
      while (inwire.good()){
	string line;
	getline(inwire,line);
	if (!inwire.good()) continue;
	if (settings_->writetrace()) {
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

	for (unsigned int i=0;i<settings_->NSector();i++) {
	  sectors_[i]->addWire(mem,procin,procout);
	}
  
      }

      ifstream indtc("../data/dtclinklayerdisk.dat");
      assert(indtc.good());
      string dtc;
      indtc >> dtc;
      while (indtc.good()){
	vector<int> tmp;
	dtclayerdisk_[dtc]=tmp;
	int layerdisk;
	indtc >> layerdisk;
	while (layerdisk>0) {
	  dtclayerdisk_[dtc].push_back(layerdisk);
	  indtc >> layerdisk;
	}
	indtc >> dtc;
      }
      
      cabling_.init("../data/calcNumDTCLinks.txt","../data/modules_T5v3_27SP_nonant_tracklet.dat");

      
    }

    void event(SLHCEvent& ev) {

      globals_->event()=&ev;

      tracks_.clear();
      
      if (settings_->writeMonitorData("Variance")) {
	StubVariance variance(ev,globals_);
      }

      eventnum_++;
      bool first=(eventnum_==1);
      
      cleanTimer_.start();
      for (unsigned int k=0;k<settings_->NSector();k++) {
	sectors_[k]->clean();
      }
      cleanTimer_.stop();

      bool hitlayer[6];
      bool hitdisk[5];
      int stublayer[6];
      int stublayer1[6][settings_->NSector()];
      int stubdisk1[5][settings_->NSector()];
      for (unsigned int ll=0;ll<6;ll++){
	hitlayer[ll]=false;
	stublayer[ll]=0;
	for (unsigned int jj=0;jj<settings_->NSector();jj++){
	  stublayer1[ll][jj]=0;
	}
      }
      for (unsigned int ll=0;ll<5;ll++){
	hitdisk[ll]=false;
	for (unsigned int jj=0;jj<settings_->NSector();jj++){
	  stubdisk1[ll][jj]=0;
	}
      }	
      
      
      int stubcount[6][24*settings_->NSector()];	
      for (unsigned int ll=0;ll<24*settings_->NSector();ll++){
	stubcount[0][ll]=0;
	stubcount[1][ll]=0;
	stubcount[2][ll]=0;
	stubcount[3][ll]=0;
	stubcount[4][ll]=0;
	stubcount[5][ll]=0;
      }    
      
      addStubTimer_.start();
      
      for (int j=0;j<ev.nstubs();j++){
	
	L1TStub stub=ev.stub(j);
	
	if (settings_->debugTracklet()) edm::LogVerbatim("Tracklet") << "Stub: layer="<<stub.layer()+1
								   <<" disk="<<stub.disk()  
								   <<" phi="<<stub.phi()
								   <<" r="<<stub.r()
								   <<" z="<<stub.z();
	
	
	double phi=stub.phi();
	phi+=0.5*settings_->dphisectorHG();
	
	if (phi<0.0) phi+=2*M_PI;
	unsigned int isector=settings_->NSector()*phi/(2*M_PI);
	assert(isector<settings_->NSector());
      
	

	if (stub.layer()<7) {
	  
	  hitlayer[stub.layer()]=true;
	  
	  stub.lorentzcor(-40.0/10000.0);
	  
	  double phi=stub.phi();
	  if (phi<0.0) phi+=2*M_PI;
	  unsigned int iphi=24*settings_->NSector()*phi/(2*M_PI);
	  assert(iphi<24*settings_->NSector());
	  double max=115.0;
	  if (stub.layer()==0) max=70.0;
	  if (std::abs(stub.z())<max) stubcount[stub.layer()][iphi]++;
	  unsigned int isector=iphi/24;
	  assert(isector<settings_->NSector());
	  stublayer1[stub.layer()][isector]++;
	  stublayer[stub.layer()]++;
	  
	} else {
	  stubdisk1[abs(stub.disk())-1][isector]++;
	  hitdisk[abs(stub.disk())-1]=true;
	}
	

	int layer=stub.layer()+1;
	int ladder=stub.ladder(); 
	int module=stub.module();
	
	string dtc=cabling_.dtc(layer,ladder,module);
	string dtcbase=dtc.substr(2,dtc.size()-2);
	if (dtc[0]=='n') {
	  dtcbase=dtc.substr(0,4)+dtc.substr(6,dtc.size()-6);
	}     
	
	cabling_.addphi(dtc,stub.phi(),layer, module);
	
	for (unsigned int k=0;k<settings_->NSector();k++) {
	  int diff=k-isector;
	  int nSector=settings_->NSector();
	  if (diff>nSector/2) diff-=settings_->NSector();
	  if (diff<-nSector/2) diff+=settings_->NSector();
	  if (abs(diff)>1) continue;
	  double phiminsect=k*2*M_PI/settings_->NSector()-0.5*(settings_->dphisectorHG()-2*M_PI/settings_->NSector())-M_PI/settings_->NSector();
	  double dphi=stub.phi()-phiminsect;
	  if (dphi>M_PI) dphi-=2*M_PI;
	  while (dphi<0.0) dphi+=2*M_PI;
	  if (dphi>settings_->dphisectorHG()) continue;
	  bool add=sectors_[k]->addStub(stub,dtcbase);
	  
	  static std::map<string,ofstream*> dtcstubs;

	  if (settings_->writeMem()) {
	    vector<string>  dtcs=cabling_.DTCs();
	    for(auto it=dtcs.begin();it!=dtcs.end();++it){
	      string dtc=*it;
	      string dtcbase=dtc.substr(2,dtc.size()-2);
	      if (dtc[0]=='n') {
		dtcbase=dtc.substr(0,4)+dtc.substr(6,dtc.size()-6);
	      }     
	      
	      string fname="../data/MemPrints/InputStubs/Link_";
	      fname+=dtcbase;
	      if (dtcstubs.find(dtcbase+"A")!=dtcstubs.end()) continue;
	      fname+="_A.dat";
	      ofstream* out=new ofstream;
	      out->open(fname.c_str());
	      dtcstubs[dtcbase+"A"]=out;
	      
	      fname="../data/MemPrints/InputStubs/Link_";
	      fname+=dtcbase;
	      if (dtcstubs.find(dtcbase+"B")!=dtcstubs.end()) continue;
	      fname+="_B.dat";
	      out=new ofstream;
	      out->open(fname.c_str());
	      dtcstubs[dtcbase+"B"]=out;
	    }
	    
	    static int oldevent=-1;
	    if (eventnum_!=oldevent) {
	      oldevent=eventnum_;
	      for(auto it=dtcstubs.begin();it!=dtcstubs.end();++it) {
		FPGAWord tmp;
	      tmp.set(eventnum_%8,3);
	      (*(it->second)) << "BX "<<tmp.str()<<" Event : "<<eventnum_+1<<endl;
	      }	       		
	    }
	  }
	  
	  if(add&&settings_->writeMem()&&k==settings_->writememsect()) {
	    Stub fpgastub(stub,settings_,sectors_[k]->phimin(),sectors_[k]->phimax());
	    FPGAWord phi=fpgastub.phi();
	    int topbit=phi.value()>>(phi.nbits()-1);
	    std::vector<int> tmp=dtclayerdisk_[dtcbase];
	    int layerdisk=stub.layer()+1;
	    if (layerdisk>999){
	      layerdisk=10+abs(stub.disk());
	    }
	    int layerdiskcode=-1;
	    for(unsigned int i=0;i<tmp.size();i++){
	      if (tmp[i]==layerdisk) layerdiskcode=i;
	    }
	    if (layerdiskcode==-1) {
	      edm::LogVerbatim("Tracklet") << "dtcbase layerdisk layer disk : "<<dtcbase<<" "<<layerdisk<<" "<<stub.layer()+1<<" "<<stub.disk();
	    }
	    assert(layerdiskcode>=0);
	    assert(layerdiskcode<4);
	    FPGAWord ldcode;
	    ldcode.set(layerdiskcode,2);
	    string dataword=ldcode.str()+"|"+fpgastub.str();
	  if (topbit==0) {
	    (*dtcstubs[dtcbase+"A"]) << dataword<<" "<<Trklet::hexFormat(dataword)<<endl;
	  } else {
	    (*dtcstubs[dtcbase+"B"]) << dataword<<" "<<Trklet::hexFormat(dataword)<<endl;
          }
	  }
	  
	}
	
	

	
	
      }

      /*
      for(int ii=0;ii<6;ii++){
	if (hitlayer[ii]) nlayershit++;
	if (ii<5) {
	  if (hitdisk[ii]) nlayershit++;
	}
      } 
      */  

      if (settings_->writeMem()) {
	for (unsigned int k=0;k<settings_->NSector();k++) {
	  if(k==settings_->writememsect())
	    sectors_[k]->writeInputStubs(first);
	}
      }
      
      
      if (settings_->writeMonitorData("StubsLayer")) {
	static ofstream out("stubslayer.txt");
	out <<stublayer[0]<<" "<<stublayer[1]<<" "<<stublayer[2]<<" "
	    <<stublayer[3]<<" "<<stublayer[4]<<" "<<stublayer[5]<<endl;
      }     


      if (settings_->writeMonitorData("StubsLayerSector")) {
	static ofstream out("stubslayerpersector.txt");
	for(unsigned int jj=0;jj<settings_->NSector();jj++){
	  out <<stublayer1[0][jj]<<" "<<stublayer1[1][jj]<<" "
	      <<stublayer1[2][jj]<<" "
	      <<stublayer1[3][jj]<<" "<<stublayer1[4][jj]<<" "
	      <<stublayer1[5][jj]<<endl; 
	}
	static ofstream out1("stubsdiskpersector.txt");
	for(unsigned int jj=0;jj<settings_->NSector();jj++){
	  out1 <<stubdisk1[0][jj]<<" "<<stubdisk1[1][jj]<<" "
	       <<stubdisk1[2][jj]<<" "
	       <<stubdisk1[3][jj]<<" "<<stubdisk1[4][jj]<<endl; 
	}
      }     
      
      
      addStubTimer_.stop();


      //Now start processing
      
      VMRouterTimer_.start();
      for (unsigned int k=0;k<settings_->NSector();k++) {
	sectors_[k]->executeVMR();	 
	if(settings_->writeMem()&&k==settings_->writememsect()) {
	  sectors_[k]->writeInputStubs(first);	 
	  sectors_[k]->writeVMSTE(first);	 
	  sectors_[k]->writeVMSME(first);	 
	  sectors_[k]->writeAS(first);	 
	}      
      }
      VMRouterTimer_.stop();
      
      TETimer_.start();
      for (unsigned int k=0;k<settings_->NSector();k++) {
	sectors_[k]->executeTE();	
      }
      TETimer_.stop();
      
      TEDTimer_.start();
      for (unsigned int k=0;k<settings_->NSector();k++) {
	sectors_[k]->executeTED();	
      }
      TEDTimer_.stop();
      
      TRETimer_.start();
      for (unsigned int k=0;k<settings_->NSector();k++) {
	sectors_[k]->executeTRE();	
	if(settings_->writeMem()&&k==settings_->writememsect()){
	  sectors_[k]->writeST(first);
	} 
      }
      TRETimer_.stop();

      for (unsigned int k=0;k<settings_->NSector();k++) {
	sectors_[k]->executeTP();	 
	if(settings_->writeMem()&&k==settings_->writememsect()){
	  sectors_[k]->writeTPAR(first);
	} 
      }
      
      for (unsigned int k=0;k<settings_->NSector();k++) {
	if(settings_->writeMem()&&k==settings_->writememsect()){
	  sectors_[k]->writeSP(first);
	} 
      }
      

      TCTimer_.start();
      for (unsigned int k=0;k<settings_->NSector();k++) {
	sectors_[k]->executeTC();	 
	if(settings_->writeMem()&&k==settings_->writememsect()){
	  sectors_[k]->writeTPAR(first);
	} 
      }
      TCTimer_.stop();
      
      int nTP=globals_->event()->nsimtracks();
      for (int iTP=0;iTP<nTP;iTP++){
	L1SimTrack simtrk=globals_->event()->simtrack(iTP);
	if (simtrk.pt()<2.0) continue;
	if (std::abs(simtrk.vz())>15.0) continue;
	if (hypot(simtrk.vx(),simtrk.vy())>0.1) continue;
	bool electron=(abs(simtrk.type())==11);
	bool muon=(abs(simtrk.type())==13);
	bool pion=(abs(simtrk.type())==211);
	bool kaon=(abs(simtrk.type())==321);
	bool proton=(abs(simtrk.type())==2212);
	if (!(electron||muon||pion||kaon||proton)) continue;
	int nlayers=0;
	int ndisks=0;
	int simtrackid=simtrk.trackid();
	unsigned int hitmask=ev.layersHit(simtrackid,nlayers,ndisks);
	if (nlayers+ndisks<4) continue;
	
	if (settings_->writeMonitorData("HitEff")) {
	  static ofstream outhit("hiteff.txt");
	  outhit << simtrk.eta()<<" "
		 <<(hitmask&1) << " " << (hitmask&2) << " "
		 << (hitmask&4) << " " << (hitmask&8) << " "	  
		 << (hitmask&16) << " " << (hitmask&32)<<" "
		 << (hitmask&64) << " " << (hitmask&128)<<" "
		 << (hitmask&256) << " " << (hitmask&512)<<" "
		 << (hitmask&1024) << endl;
	}
	
	
	std::set<int> matchseed;
	for (unsigned int k=0;k<settings_->NSector();k++) {
	  std::set<int> matchseedtmp=sectors_[k]->seedMatch(iTP);
	  matchseed.insert(matchseedtmp.begin(),matchseedtmp.end());
	}
	if (settings_->bookHistos()) {
	  for(int iseed=0;iseed<8;iseed++){
	    bool eff=matchseed.find(iseed)!=matchseed.end();
	    globals_->histograms()->fillSeedEff(iseed,simtrk.eta(),eff);
	  }
	}
      } 
      
      TCDTimer_.start();
      for (unsigned int k=0;k<settings_->NSector();k++) {
	sectors_[k]->executeTCD();	 
	if(settings_->writeMem()&&k==settings_->writememsect()){
	  sectors_[k]->writeTPAR(first);
	  sectors_[k]->writeTPROJ(first);
	} 
      }
      TCDTimer_.stop();
      
      
      PRTimer_.start();
      for (unsigned int k=0;k<settings_->NSector();k++) {
	sectors_[k]->executePR();	
	if(settings_->writeMem()&&k==settings_->writememsect()){
	  sectors_[k]->writeVMPROJ(first);
	  sectors_[k]->writeAP(first);
	}
      }
      PRTimer_.stop();
      
      METimer_.start();
      for (unsigned int k=0;k<settings_->NSector();k++) {
	sectors_[k]->executeME();	
	if(settings_->writeMem()&&k==settings_->writememsect()){
	  sectors_[k]->writeCM(first);
	} 
      }
      METimer_.stop();
      
      MCTimer_.start();
      for (unsigned int k=0;k<settings_->NSector();k++) {
	sectors_[k]->executeMC();
      }
      MCTimer_.stop();
      
      
      MPTimer_.start();
      for (unsigned int k=0;k<settings_->NSector();k++) {
	sectors_[k]->executeMP();
      }
      MPTimer_.stop();

      for (unsigned int k=0;k<settings_->NSector();k++) {
	if(settings_->writeMem()&&k==settings_->writememsect()){
	  sectors_[k]->writeMC(first);
	}
      }
      
      
      
      FTTimer_.start();
      for (unsigned int k=0;k<settings_->NSector();k++) {
	sectors_[k]->executeFT();	 
	if((settings_->writeMem()||settings_->writeMonitorData("IFit"))&&k==settings_->writememsect()){
	  sectors_[k]->writeTF(first);
	}
      }
      FTTimer_.stop();
      
      PDTimer_.start();
      for (unsigned int k=0;k<settings_->NSector();k++) {
	sectors_[k]->executePD(tracks_);	  
	if(((settings_->writeMem()||settings_->writeMonitorData("IFit"))&&k==settings_->writememsect()) || settings_->writeMonitorData("CT")){
	  sectors_[k]->writeCT(first);
	}
      }
      PDTimer_.stop();
      
    }

    void printSummary() {
      
      if (settings_->writeMonitorData("Cabling")) {
	cabling_.writephirange();
      }

      if (settings_->bookHistos()) {
	globals_->histograms()->close();
      }

    
      edm::LogVerbatim("Tracklet") << "Process             Times called   Average time (ms)      Total time (s)";
      edm::LogVerbatim("Tracklet") << "Cleaning              "
				   <<setw(10)<<cleanTimer_.ntimes()
				   <<setw(20)<<setprecision(3)<<cleanTimer_.avgtime()*1000.0
				   <<setw(20)<<setprecision(3)<<cleanTimer_.tottime();
      edm::LogVerbatim("Tracklet") << "Add Stubs             "
				   <<setw(10)<<addStubTimer_.ntimes()
				   <<setw(20)<<setprecision(3)<<addStubTimer_.avgtime()*1000.0
				   <<setw(20)<<setprecision(3)<<addStubTimer_.tottime();
      edm::LogVerbatim("Tracklet") << "VMRouter              "
			       <<setw(10)<<VMRouterTimer_.ntimes()
				   <<setw(20)<<setprecision(3)<<VMRouterTimer_.avgtime()*1000.0
			       <<setw(20)<<setprecision(3)<<VMRouterTimer_.tottime();
      edm::LogVerbatim("Tracklet") << "TrackletEngine        "
				   <<setw(10)<<TETimer_.ntimes()
				   <<setw(20)<<setprecision(3)<<TETimer_.avgtime()*1000.0
				   <<setw(20)<<setprecision(3)<<TETimer_.tottime();
      edm::LogVerbatim("Tracklet") << "TrackletEngineDisplaced"
				   <<setw(10)<<TEDTimer_.ntimes()
				   <<setw(20)<<setprecision(3)<<TEDTimer_.avgtime()*1000.0
				   <<setw(20)<<setprecision(3)<<TEDTimer_.tottime();
      edm::LogVerbatim("Tracklet") << "TripletEngine         "
				   <<setw(10)<<TRETimer_.ntimes()
				   <<setw(20)<<setprecision(3)<<TRETimer_.avgtime()*1000.0
				   <<setw(20)<<setprecision(3)<<TRETimer_.tottime();
      edm::LogVerbatim("Tracklet") << "TrackletCalculator    "
				   <<setw(10)<<TCTimer_.ntimes()
				   <<setw(20)<<setprecision(3)<<TCTimer_.avgtime()*1000.0
				   <<setw(20)<<setprecision(3)<<TCTimer_.tottime();
      edm::LogVerbatim("Tracklet") << "TrackletCalculatorDisplaced"
				   <<setw(10)<<TCDTimer_.ntimes()
				   <<setw(20)<<setprecision(3)<<TCDTimer_.avgtime()*1000.0
				   <<setw(20)<<setprecision(3)<<TCDTimer_.tottime();
      edm::LogVerbatim("Tracklet") << "ProjectionRouter      "
				   <<setw(10)<<PRTimer_.ntimes()
				   <<setw(20)<<setprecision(3)<<PRTimer_.avgtime()*1000.0
				   <<setw(20)<<setprecision(3)<<PRTimer_.tottime();
      edm::LogVerbatim("Tracklet") << "MatchEngine           "
				   <<setw(10)<<METimer_.ntimes()
				   <<setw(20)<<setprecision(3)<<METimer_.avgtime()*1000.0
				   <<setw(20)<<setprecision(3)<<METimer_.tottime();
      edm::LogVerbatim("Tracklet") << "MatchCalculator       "
				   <<setw(10)<<MCTimer_.ntimes()
				   <<setw(20)<<setprecision(3)<<MCTimer_.avgtime()*1000.0
				   <<setw(20)<<setprecision(3)<<MCTimer_.tottime();
      edm::LogVerbatim("Tracklet") << "MatchProcessor        "
				   <<setw(10)<<MPTimer_.ntimes()
				   <<setw(20)<<setprecision(3)<<MPTimer_.avgtime()*1000.0
				   <<setw(20)<<setprecision(3)<<MPTimer_.tottime();
      edm::LogVerbatim("Tracklet") << "FitTrack              "
				   <<setw(10)<<FTTimer_.ntimes()
				   <<setw(20)<<setprecision(3)<<FTTimer_.avgtime()*1000.0
				   <<setw(20)<<setprecision(3)<<FTTimer_.tottime();
      edm::LogVerbatim("Tracklet") << "PurgeDuplicate        "
				   <<setw(10)<<PDTimer_.ntimes()
				   <<setw(20)<<setprecision(3)<<PDTimer_.avgtime()*1000.0
				   <<setw(20)<<setprecision(3)<<PDTimer_.tottime();
      
    }
      
    std::vector<Track*>& tracks() {return tracks_;}

  private:

    const Settings* settings_{0};

    GlobalHistTruth* globals_{};
    
    Sector** sectors_{};

    int eventnum_={0};
    
    Cabling cabling_;
    
    CPUTimer cleanTimer_;
    CPUTimer addStubTimer_;
    CPUTimer VMRouterTimer_;  
    CPUTimer TETimer_;
    CPUTimer TEDTimer_;
    CPUTimer TRETimer_;
    CPUTimer TCTimer_;
    CPUTimer TCDTimer_;
    CPUTimer PRTimer_;
    CPUTimer METimer_;
    CPUTimer MCTimer_;
    CPUTimer MPTimer_;
    CPUTimer FTTimer_;
    CPUTimer PDTimer_;

    std::vector<Track*> tracks_;

    std::map<string,vector<int> > dtclayerdisk_;
    
  };
};
      
#endif
      

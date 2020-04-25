#ifndef L1Trigger_TrackFindingTracklet_interface_Settings_h
#define L1Trigger_TrackFindingTracklet_interface_Settings_h

#include <array>
#include <set>
#include <map>

namespace Trklet{

  class Settings{
    
  public:

    Settings(){

      //Uncomment to run the hybrid algorithm
      //#ifdef CMSSW_GIT_HASH
      //#define USEHYBRID
      //#endif

      geomTkTDR_=false;
      nzbitsdisk_=7;

      NSector_=9;
      rcrit_=55.0;

      rinvmax_=0.01*0.3*3.8/2.0; //0.01 to convert to cm-1 -  FIXME should not have all these hardcoded numbers
      
      nzbitsstub_={{12,12,12,8,8,8,7,7,7,7,7}};
      nphibitsstub_={{14,14,14,17,17,17,14,14,14,14,14}};
      nrbitsstub_={{7,7,7,7,7,7,12,12,12,12,12}};

      nrbitsprojderdisk_=9;
      nbitsphiprojderL123_=8+2;
      nbitsphiprojderL456_=8+2;
      nbitszprojderL123_=8+2;
      nbitszprojderL456_=7+2;

      
      useseeding_={{0,1,2,3,4,5,6,7,8,9,10,11}};

      nbitsvmte_[0]={{2,2,2,2,2,2,1,1,2,2,3,2}};
      nbitsvmte_[1]={{3,2,3,3,2,2,2,2,3,3,2,2}};
      nbitsvmte_[2]={{0,0,0,0,0,0,0,0,0,0,2,1}};

      nbitsvmme_={{2,3,3,3,3,3,3,2,2,2,2}};
      
      nbitsallstubs_={{3,2,2,2,2,2,2,2,2,2,2}}; 

      // geometry related 
      zlength_=120.0;
      rmaxdisk_=120.0;

      drmax_ = rmaxdisk_/32.0;
      dzmax_ = zlength_/32.0;

      rmean_={{geomTkTDR_?(rmaxdisk_*858)/4096:(rmaxdisk_*851)/4096,  //FIXME - can not depend on geomTkTDR
	      geomTkTDR_?(rmaxdisk_*1279)/4096:(rmaxdisk_*1269)/4096,
	      geomTkTDR_?(rmaxdisk_*1795)/4096:(rmaxdisk_*1784)/4096,
	      geomTkTDR_?(rmaxdisk_*2347)/4096:(rmaxdisk_*2347)/4096,
	      geomTkTDR_?(rmaxdisk_*2937)/4096:(rmaxdisk_*2936)/4096,
	      geomTkTDR_?(rmaxdisk_*3783)/4096:(rmaxdisk_*3697)/4096}};

      zmean_={{(zlength_*2239)/2048,
	      (zlength_*2645)/2048,
	      (zlength_*3163)/2048,
	      (zlength_*3782)/2048,
	      (zlength_*4523)/2048}};

      // stub bend cuts
      bendcutte_[0]={{1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25}};   //inner
      bendcutte_[1]={{1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25}};   //outer

      nfinephi_[0]={{2,2,2,2,2,2,2,2,2,2,2,2}};     //inner
      nfinephi_[1]={{3,3,3,3,3,3,3,3,3,3,3,3}};     //outer
      nfinephi_[2]={{0,0,0,0,0,0,0,0,0,0,3,3}};   //outermost

      //These are the number of bits used for the VM regions in the TE by seedindex
      nphireg_[0]={{5,4,4,4,4,4,4,3,4,4,5,4}};     //inner
      nphireg_[1]={{5,4,5,5,4,4,4,4,4,4,4,4}};    //outer
      nphireg_[2]={{0,0,0,0,0,0,0,0,0,0,4,4}};   //outermost


      zbitstab_[0]={{7,7,7,7,3,3,7,7,0,0,7,0}};
      zbitstab_[1]={{7,7,7,7,3,3,3,3,0,0,7,0}};
      zbitstab_[2]={{0,0,0,0,0,0,0,0,0,0,3,7}};

      rbitstab_[0]={{4,4,4,4,8,8,3,3,0,0,4,0}};
      rbitstab_[1]={{4,4,4,4,7,7,7,7,0,0,4,0}};
      rbitstab_[2]={{0,0,0,0,0,0,0,0,0,0,7,4}};

      lutwidthtab_[0]={{10,11,11,11,11,11,11,11, 0, 0,11, 0}};
      lutwidthtab_[1]={{ 6, 6, 6, 6,10,10,10,10, 0, 0, 6, 0}};
      lutwidthtab_[2]={{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 6}};

      lutwidthtabextended_[0]={{11,11,21,21,21,21,11,11, 0, 0,21, 0}};
      lutwidthtabextended_[1]={{ 6, 6, 6, 6,10,10,10,10, 0, 0, 6, 0}};
      lutwidthtabextended_[2]={{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 6}};

      //projection layers by seed index. For each seeding index (row) the list of layers that we consider projections to
      projlayers_[0]={{3, 4, 5, 6}};  //0 L1L2
      projlayers_[1]={{1, 4, 5, 6}};  //1 L2L3
      projlayers_[2]={{1, 2, 5, 6}};  //2 L3L4
      projlayers_[3]={{1, 2, 3, 4}};  //3 L5L6    
      projlayers_[4]={{1, 2}};        //4 D1D2    
      projlayers_[5]={{1}};           //5 D3D4    
      projlayers_[6]={{}};            //6 L1D1    
      projlayers_[7]={{1}};           //7 L2D1    
      projlayers_[8]={{1, 5, 6}};     //8 L2L3L4
      projlayers_[9]={{1, 2, 3}};     //9 L4L5L6
      projlayers_[10]={{1}};          //10 L2L3D1
      projlayers_[11]={{1}};            //11 D1D2L2

      //projection disks by seed index. For each seeding index (row) the list of diks that we consider projections to
      projdisks_[0]={{1, 2, 3, 4}};    //0 L1L2
      projdisks_[1]={{1, 2, 3, 4}};    //1 L2L3    
      projdisks_[2]={{1, 2}};          //2 L3L4
      projdisks_[3]={{}};              //3 L5L6
      projdisks_[4]={{3, 4, 5}};       //4 D1D2    
      projdisks_[5]={{1, 2, 5}};       //5 D3D4    
      projdisks_[6]={{2, 3, 4, 5}};    //6 L1D1    
      projdisks_[7]={{2, 3, 4}};       //7 L2D1    
      projdisks_[8]={{1, 2}};          //8 L2L3L4
      projdisks_[9]={{}};              //9 L4L5L6
      projdisks_[10]={{2, 3, 4}};       //10 L2L3D1
      projdisks_[11]={{3, 4}};           //11 D1D2L2

      //rphi cuts for layers - the column is the seedindex
      rphimatchcut_[0]={{0.0, 0.1, 0.07, 0.08, 0.07, 0.05, 0.0, 0.05, 0.08, 0.15, 0.125, 0.15}};   //Layer 1
      rphimatchcut_[1]={{0.0, 0.0, 0.06, 0.08, 0.05, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0}};   //Layer 2
      rphimatchcut_[2]={{0.1, 0.0, 0.0, 0.08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.08, 0.0, 0.0}};  //Layer 3
      rphimatchcut_[3]={{0.19, 0.19, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};   //Layer 4
      rphimatchcut_[4]={{0.4, 0.4, 0.08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.08, 0.0, 0.0, 0.0}};   //Layer 5
      rphimatchcut_[5]={{0.5, 0.0, 0.19, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0}};   //Layer 6
      
      //z cuts for layers - the column is the seedindex
      zmatchcut_[0]={{0.0, 0.7, 5.5, 15.0, 1.5, 2.0, 0.0, 1.5, 1.0, 8.0, 1.0, 1.5}};   //Layer 1
      zmatchcut_[1]={{0.0, 0.0, 3.5, 15.0, 1.25, 0.0, 0.0, 0.0, 0.0, 7.0, 0.0, 0.0}};   //Layer 2
      zmatchcut_[2]={{0.7, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0}};   //Layer 3
      zmatchcut_[3]={{3.0, 3.0, 0.0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};  //Layer 4
      zmatchcut_[4]={{3.0, 3.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.5, 0.0, 0.0, 0.0}};   //Layer 5
      zmatchcut_[5]={{4.0, 0.0, 9.5, 0.0, 0.0, 0.0, 0.0, 0.0, 4.5, 0.0, 0.0, 0.0}};   //Layer 6
      
      //rphi cuts for PS modules in disks - the column is the seedindex
      rphicutPS_[0]={{0.2, 0.2, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};   //disk 1
      rphicutPS_[1]={{0.2, 0.2, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.0, 0.0, 0.15, 0.0}};   //disk 2
      rphicutPS_[2]={{0.25, 0.2, 0.0, 0.0, 0.15, 0.0, 0.2, 0.15, 0.0, 0.0, 0.0, 0.2}};   //disk 3
      rphicutPS_[3]={{0.5, 0.2, 0.0, 0.0, 0.2, 0.0, 0.3, 0.5, 0.0, 0.0, 0.0, 0.0}};   //disk 4
      rphicutPS_[4]={{0.0, 0.0, 0.0, 0.0, 0.25, 0.1, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0}};   //disk 5

      //r cuts for PS modules in disks - the column is the seedindex
      rcutPS_[0]={{0.5, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};   //disk 1
      rcutPS_[1]={{0.5, 0.5, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.0, 0.0, 0.5, 0.0}};   //disk 2
      rcutPS_[2]={{0.5, 0.5, 0.0, 0.0, 0.5, 0.0, 0.6, 0.8, 0.0, 0.0, 0.0, 0.4}};   //disk 3
      rcutPS_[3]={{0.5, 0.5, 0.0, 0.0, 0.8, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0}};   //disk 4
      rcutPS_[4]={{0.0, 0.0, 0.0, 0.0, 1.0, 0.5, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0}};   //disk 5
      
      //rphi cuts for 2S modules in disks = the column is the seedindex
      rphicut2S_[0]={{0.5, 0.5, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0}};   //disk 1
      rphicut2S_[1]={{0.5, 0.5, 0.8, 0.0, 0.0, 0.0, 0.5, 0.15, 0.3, 0.0, 0.68, 0.0}};   //disk 2
      rphicut2S_[2]={{0.5, 0.5, 0.0, 0.0, 0.15, 0.0, 0.2, 0.25, 0.0, 0.0, 0.8, 0.1}};   //disk 3
      rphicut2S_[3]={{0.5, 0.5, 0.0, 0.0, 0.2, 0.0, 0.25, 0.5, 0.0, 0.0, 0.6, 0.4}};   //disk 4
      rphicut2S_[4]={{0.0, 0.0, 0.0, 0.0, 0.4, 0.2, 0.4, 0.0, 0.0, 0.0, 0.0, 0.8}};   //disk 5
      
      //r cuts for 2S modules in disks -the column is the seedindex
      rcut2S_[0]={{3.8, 3.8, 3.8, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0}};   //disk 1
      rcut2S_[1]={{3.8, 3.8, 3.8, 0.0, 0.0, 0.0, 3.8, 3.4, 3.0, 0.0, 3.0, 0.0}};   //disk 2
      rcut2S_[2]={{3.6, 3.8, 0.0, 0.0, 3.6, 0.0, 3.6, 3.8, 0.0, 0.0, 3.8, 3.0}};  //disk 3
      rcut2S_[3]={{3.6, 3.8, 0.0, 0.0, 3.6, 0.0, 3.5, 3.8, 0.0, 0.0, 3.0, 3.0}};   //disk 4
      rcut2S_[4]={{0.0, 0.0, 0.0, 0.0, 3.6, 3.4, 3.7, 0.0, 0.0, 0.0, 0.0, 3.0}};   //disk 5

      maxstepoffset_=10000;
      maxstep_={{"Link",108},
		{"MC",108},
		{"ME",108},
		{"MP",108},
		{"PR",108},
		{"TC",108},
		{"TE",108},
		{"TP",108},
		{"TRE",108},
		{"VMR",108}};

      writeMonitorData_={{"IL",false},
			 {"TE",false},
			 {"CT",false},
			 {"HitPattern",false},
			 {"ChiSq",false},
			 {"Seeds",false},
			 {"FT",false},
			 {"Residuals",false},
			 {"MC",false},
			 {"ME",false},
			 {"AP",false},
			 {"VMP",false},
			 {"NMatches",false},
			 {"TrackProjOcc",false},
			 {"TC",false},
			 {"Pars",false},
			 {"TPars",false},
			 {"TPD",false},
			 {"TrackletPars",false},
			 {"TED",false},
			 {"TP",false},
			 {"TRE",false},
			 {"VMR",false},
			 {"Variance",false},
			 {"StubsLayer",false},
			 {"StubsLayerSector",false},
			 {"ResEff",false},
			 {"HitEff",false},
			 {"MatchEff",false},
			 {"Cabling",false},
			 {"IFit",false},
			 {"AS",false}};


      double rDSSinner_mod1 = geomTkTDR_?69.2345:68.9391;
      double rDSSinner_mod2 = geomTkTDR_?80.0056:78.7750;
      double rDSSinner_mod3 = geomTkTDR_?87.3444:85.4550;
      double rDSSinner_mod4 = geomTkTDR_?98.2515:96.3150;
      double rDSSinner_mod5 = geomTkTDR_?104.9750:102.3160;
      
      double rDSSouter_mod1 = geomTkTDR_?67.6317:66.4903;
      double rDSSouter_mod2 = geomTkTDR_?78.1300:76.7750;
      double rDSSouter_mod3 = geomTkTDR_?86.4293:84.4562;
      double rDSSouter_mod4 = geomTkTDR_?97.1316:94.9920;
      double rDSSouter_mod5 = geomTkTDR_?104.9750:102.3160;

      double halfstrip = 2.5; //we want the center of the two strip positions in a module, not just the center of a module 

      rDSSinner_ = {{rDSSinner_mod1-halfstrip, rDSSinner_mod1+halfstrip,  //FIXME can not depend on geomTkTDR
		    rDSSinner_mod2-halfstrip, rDSSinner_mod2+halfstrip,
		    rDSSinner_mod3-halfstrip, rDSSinner_mod3+halfstrip,
		    rDSSinner_mod4-halfstrip, rDSSinner_mod4+halfstrip,
		    rDSSinner_mod5-halfstrip, rDSSinner_mod5+halfstrip}};

      rDSSouter_ = {{rDSSouter_mod1-halfstrip, rDSSouter_mod1+halfstrip,
		    rDSSouter_mod2-halfstrip, rDSSouter_mod2+halfstrip,
		    rDSSouter_mod3-halfstrip, rDSSouter_mod3+halfstrip, 
		    rDSSouter_mod4-halfstrip, rDSSouter_mod4+halfstrip,
		    rDSSouter_mod5-halfstrip, rDSSouter_mod5+halfstrip}};

      // various printouts for debugging and warnings 
      printDebugKF_=false; // if true print lots of debugging statements related to the KF fit
      debugTracklet_=false; //Print detailed debug information about tracklet tracking
      writetrace_=false; //Print out details about parsing configuration files
      
      warnNoMem_=false;  //If true will print out warnings about missing projection memories
      warnNoDer_=false;  //If true will print out warnings about missing track fit derivatives

      writeMem_=false;    //If true will print out content of memories to files
      writeTable_=false;  //IF true will print out content of LUTs to files

      // Write various lookup tables and autogenerated code (from iMath)
      writeVerilog_=false;     //Write out auto-generated Verilog mudules used by TCs
      writeHLS_=false;         //Write out auto-generated HLS mudules used by TCs
      writeInvTable_=false;    //Write out tables of drinv and invt in tracklet calculator for Verilog module
      writeHLSInvTable_=false; //Write out tables of drinv and invt in tracklet calculator for HLS module

      writememsect_=3;         //writemem only for this sector (note that the files will have _4 extension)
      
      writeTripletTables_=false; //Train and write the TED and TRE tables. N.B.: the tables
                                 //cannot be applied while they are being trained, i.e.,
                                 //this flag effectively turns off the cuts in
                                 //TrackletEngineDisplaced and TripletEngine

      writeoutReal_ = false; 

      
      bookHistos_=false; //set to true/false to turn on/off histogram booking internal to the tracking (class "HistImp")
      
      // pt constants
      ptcut_ = 1.91; //Minimum pt
      rinvcut_ = 0.01*0.3*3.8/ptcut_; //0.01 to convert to cm-1

      // Parameters for bit sizes
      alphashift_ = 12;  
      nbitsalpha_ = 4;  //bits used to store alpha
      alphaBitsTable_ = 2; //For number of bits in track derivative table
      nrinvBitsTable_ = 3; //number of bits for tabulating rinv dependence
      
      MEBinsBits_ = 3;
      MEBins_ = (1<<MEBinsBits_);
      MEBinsDisks_ = 8; //on each side



      
      // Options for chisq fit
      useMSFit_ = false;
      exactderivatives_ = false;  //for both the integer and float
      exactderivativesforfloating_ = true; //only for the floating point
      useapprox_ = true; //use approximate postion based on integer representation for floating point
      usephicritapprox_ = false; //use floating point approximate version of phicrit cut if true

      // Duplicate Removal
      // "merge" (hybrid dup removal)
      // "ichi" (pairwise, keep track with best ichisq), "nstub" (pairwise, keep track with more stubs)
      // "grid" (TMTT-like removal), "" (no removal)
      minIndStubs_ = 3; // not used with merge removal
      removalType_ = "ichi";
      doKF_ = false;
      
#ifdef USEHYBRID
      removalType_ = "merge";
      // "CompareBest" (recommended) Compares only the best stub in each track for each region (best = smallest phi residual)
      // and will merge the two tracks if stubs are shared in three or more regions
      // "CompareAll" Compares all stubs in a region, looking for matches, and will merge the two tracks if stubs are shared in three or more regions
      mergeComparison_ = "CompareBest";
      doKF_ = true; 
#endif
      fakefit_ = false; // if true, run a dummy fit, producing TTracks directly from output of tracklet pattern reco stage

      
      // configurable options
      nHelixPar_ = 4; // 4 or 5 param helix fit
      extended_ = false; // turn on displaced tracking
      
      fitpatternfile_ = "../data/fitpattern.txt"; // list of the different hit patterns for fits (read in through python config for CMSSW running)
      skimfile_ = ""; // if this string is non-empty, write ascii file with processed events
      
    }

    bool geomTkTDR() const {return geomTkTDR_;}
    unsigned int nzbitsdisk() const {return nzbitsdisk_;}
    unsigned int nzbitsstub(unsigned int layerdisk) const {return nzbitsstub_[layerdisk];}
    unsigned int nphibitsstub(unsigned int layerdisk) const {return nphibitsstub_[layerdisk];}
    unsigned int nrbitsstub(unsigned int layerdisk) const {return nrbitsstub_[layerdisk];}

    unsigned int nrbitsprojderdisk() const {return nrbitsprojderdisk_;}
    unsigned int nbitsphiprojderL123() const {return nbitsphiprojderL123_;}
    unsigned int nbitsphiprojderL456() const {return nbitsphiprojderL456_;}
    unsigned int nbitszprojderL123() const {return nbitszprojderL123_;}
    unsigned int nbitszprojderL456() const {return nbitszprojderL456_;}

    bool useSeed(unsigned int iSeed) const {return useseeding_.find(iSeed)!=useseeding_.end();}
    unsigned int nbitsvmte(unsigned int inner, unsigned int iSeed) const {return nbitsvmte_[inner][iSeed];}
    unsigned int nvmte(unsigned int inner, unsigned int iSeed) const {return (1<<nbitsvmte_[inner][iSeed]);}

    unsigned int nbitsvmme(unsigned int layerdisk) const {return nbitsvmme_[layerdisk];}
    unsigned int nvmme(unsigned int layerdisk) const {return (1<<nbitsvmme_[layerdisk]);}

    unsigned int nbitsallstubs(unsigned int layerdisk) const {return nbitsallstubs_[layerdisk];}
    unsigned int nallstubs(unsigned int layerdisk) const {return (1<<nbitsallstubs_[layerdisk]);}

    bool writeMonitorData(std::string module) const {
      if (writeMonitorData_.find(module)==writeMonitorData_.end()){
	std::cout << "Settings::writeMonitorData module = "<<module<<" not known" << std::endl;
	assert(0);
      }
      return writeMonitorData_.at(module);
    }

    unsigned int maxStep(std::string module) const {
      if (maxstep_.find(module)==maxstep_.end()){
	std::cout << "Settings::maxStep module = "<<module<<" not known" << std::endl;
	assert(0);
      }
      return maxstep_.at(module)+maxstepoffset_;;
    }
    

    double zlength() const { return zlength_;}
    double rmaxdisk() const { return rmaxdisk_;}

    double drmax() const {return drmax_;}
    double dzmax() const {return dzmax_;}


    
    double bendcutte(unsigned int inner, unsigned int iSeed) const {return bendcutte_[inner][iSeed];}
    double nfinephi(unsigned int inner, unsigned int iSeed) const {return nfinephi_[inner][iSeed];}
    double nphireg(unsigned int inner, unsigned int iSeed) const {return nphireg_[inner][iSeed];}
    double zbitstab(unsigned int inner, unsigned int iSeed) const {return zbitstab_[inner][iSeed];}
    double rbitstab(unsigned int inner, unsigned int iSeed) const {return rbitstab_[inner][iSeed];}
    double lutwidthtab(unsigned int inner, unsigned int iSeed) const {return lutwidthtab_[inner][iSeed];}
    double lutwidthtabextended(unsigned int inner, unsigned int iSeed) const {return lutwidthtabextended_[inner][iSeed];}

    double projlayers(unsigned int iSeed, unsigned int i) const {return projlayers_[iSeed][i];}
    double projdisks(unsigned int iSeed, unsigned int i) const {return projdisks_[iSeed][i];}
    double rphimatchcut(unsigned int iSeed, unsigned int ilayer) const {return rphimatchcut_[ilayer][iSeed];}
    double zmatchcut(unsigned int iSeed, unsigned int ilayer) const {return zmatchcut_[ilayer][iSeed];}
    double rphicutPS(unsigned int iSeed, unsigned int idisk) const {return rphicutPS_[idisk][iSeed];}
    double rcutPS(unsigned int iSeed, unsigned int idisk) const {return rcutPS_[idisk][iSeed];}
    double rphicut2S(unsigned int iSeed, unsigned int idisk) const {return rphicut2S_[idisk][iSeed];}
    double rcut2S(unsigned int iSeed, unsigned int idisk) const {return rcut2S_[idisk][iSeed];}

    double rmean(unsigned int iLayer) const {return rmean_[iLayer];}
    double zmean(unsigned int iDisk) const {return zmean_[iDisk];}

    double rDSSinner(unsigned int iBin) const {return rDSSinner_[iBin];} 
    double rDSSouter(unsigned int iBin) const {return rDSSouter_[iBin];} 

    
    bool printDebugKF() const {return printDebugKF_;}
    bool debugTracklet() const {return debugTracklet_;}
    bool writetrace() const {return writetrace_;}
    
    bool warnNoMem() const {return warnNoMem_;}
    bool warnNoDer() const {return warnNoDer_;}

    bool writeMem() const {return writeMem_;}
    bool writeTable() const {return writeTable_;}

    bool writeVerilog() const {return writeVerilog_;}
    bool writeHLS() const {return writeHLS_;}
    bool writeInvTable() const {return writeInvTable_;}
    bool writeHLSInvTable() const {return writeHLSInvTable_;}

    unsigned int writememsect() const {return writememsect_;}

    bool writeTripletTables() const {return writeTripletTables_;}

    bool writeoutReal() const {return writeoutReal_;}
 
    bool bookHistos() const {return bookHistos_;}
    
    double ptcut() const {return ptcut_;}
    double rinvcut() const {return rinvcut_;}

    int alphashift() const {return alphashift_;}
    int nbitsalpha() const {return nbitsalpha_;}
    int alphaBitsTable() const {return alphaBitsTable_;}
    int nrinvBitsTable() const {return nrinvBitsTable_;}
      
    unsigned int MEBinsBits() const {return MEBinsBits_;}
    unsigned int MEBins() const {return MEBins_;}
    unsigned int MEBinsDisks() const {return MEBinsDisks_;}

    std::string geomext() const {return extended_?"hourglassExtended":"hourglass";}  

    bool useMSFit() const {return useMSFit_;}
    bool exactderivatives() const {return exactderivatives_;}
    bool exactderivativesforfloating() const {return exactderivativesforfloating_;}
    bool useapprox() const {return useapprox_;}
    bool usephicritapprox() const {return usephicritapprox_;}
    
    unsigned int minIndStubs() const {return minIndStubs_;}
    std::string removalType() const {return removalType_;}
    std::string mergeComparison() const {return mergeComparison_;}
    bool doKF() const {return doKF_;}
    bool fakefit() const {return fakefit_;}

    
    // configurable 
    unsigned int nHelixPar() const {return nHelixPar_;}
    void setNHelixPar(unsigned int nHelixPar) {nHelixPar_ = nHelixPar;}

    bool extended() const {return extended_;}
    void setExtended(bool extended) {extended_ = extended;}
    
    std::string fitpatternfile() const {return fitpatternfile_;}
    void setFitpatternfile(std::string fitpatternfile) {fitpatternfile_ = fitpatternfile;}

    std::string skimfile() const {return skimfile_;}
    void setSkimfile(std::string skimfile) {skimfile_ = skimfile;}

    double dphisectorHG() const { return 2*M_PI/NSector_+2*fmax(std::abs(asin(0.5*rinvmax_*rmean(0))-asin(0.5*rinvmax_*rcrit_)),
								std::abs(asin(0.5*rinvmax_*rmean(5))-asin(0.5*rinvmax_*rcrit_)));}
    
    double phicritmin() const { return 0.5*dphisectorHG()-M_PI/NSector_; }
    double phicritmax() const { return dphisectorHG()-0.5*dphisectorHG()+M_PI/NSector_; }

    
  private:

    bool geomTkTDR_;

    unsigned int NSector_;

    double rcrit_;
    double rinvmax_;
    
    std::array<double,6> rmean_;
    std::array<double,5> zmean_;
    
    unsigned int nzbitsdisk_;

    std::array<unsigned int,11> nzbitsstub_;
    std::array<unsigned int,11> nphibitsstub_;
    std::array<unsigned int,11> nrbitsstub_;

    unsigned int nrbitsprojderdisk_;
    unsigned int nbitsphiprojderL123_;
    unsigned int nbitsphiprojderL456_;
    unsigned int nbitszprojderL123_;
    unsigned int nbitszprojderL456_;

    std::set<unsigned int> useseeding_;

    std::array<unsigned int,11> nbitsallstubs_;
    std::array<unsigned int,11> nbitsvmme_;
    std::array<std::array<unsigned int, 12>, 3> nbitsvmte_;

    std::array<std::array<double, 8>, 2> bendcutte_;

    double zlength_;
    double rmaxdisk_;

    double drmax_;
    double dzmax_;

    std::array<std::array<unsigned int, 12>, 3> nfinephi_;
    std::array<std::array<unsigned int, 12>, 3> nphireg_;
    std::array<std::array<unsigned int, 12>, 3> zbitstab_;
    std::array<std::array<unsigned int, 12>, 3> rbitstab_;
    std::array<std::array<unsigned int, 12>, 3> lutwidthtab_;
    std::array<std::array<unsigned int, 12>, 3> lutwidthtabextended_;

    std::array<std::array<unsigned int, 4>, 12> projlayers_;
    std::array<std::array<unsigned int, 5>, 12> projdisks_;
    std::array<std::array<double, 12>, 6> rphimatchcut_;
    std::array<std::array<double, 12>, 6> zmatchcut_;
    std::array<std::array<double, 12>, 5> rphicutPS_;
    std::array<std::array<double, 12>, 5> rcutPS_;
    std::array<std::array<double, 12>, 5> rphicut2S_;
    std::array<std::array<double, 12>, 5> rcut2S_;  

    unsigned int maxstepoffset_;
    std::map<std::string,unsigned int> maxstep_;
    std::map<std::string,bool> writeMonitorData_;

    std::array<double,10> rDSSinner_; 
    std::array<double,10> rDSSouter_;
			       
    bool printDebugKF_;
    bool debugTracklet_;
    bool writetrace_;
    
    bool warnNoMem_;
    bool warnNoDer_;

    bool writeMem_;
    bool writeTable_;

    bool writeVerilog_;
    bool writeHLS_;
    bool writeInvTable_;
    bool writeHLSInvTable_;

    unsigned int writememsect_;

    bool writeTripletTables_;

    bool writeoutReal_;
    
    bool bookHistos_;
	  
    double ptcut_;
    double rinvcut_;

    int alphashift_;
    int nbitsalpha_;
    int alphaBitsTable_;
    int nrinvBitsTable_;
    
    unsigned int MEBinsBits_;
    unsigned int MEBins_;
    unsigned int MEBinsDisks_;

    /*
    mutable IMATH_TrackletCalculator* ITC_L1L2_;
    mutable IMATH_TrackletCalculator* ITC_L2L3_;
    mutable IMATH_TrackletCalculator* ITC_L3L4_;
    mutable IMATH_TrackletCalculator* ITC_L5L6_;
    
    mutable IMATH_TrackletCalculatorDisk* ITC_F1F2_;
    mutable IMATH_TrackletCalculatorDisk* ITC_F3F4_;
    mutable IMATH_TrackletCalculatorDisk* ITC_B1B2_;
    mutable IMATH_TrackletCalculatorDisk* ITC_B3B4_;
    
    mutable IMATH_TrackletCalculatorOverlap* ITC_L1F1_;
    mutable IMATH_TrackletCalculatorOverlap* ITC_L2F1_;
    mutable IMATH_TrackletCalculatorOverlap* ITC_L1B1_;
    mutable IMATH_TrackletCalculatorOverlap* ITC_L2B1_;
    */

    bool useMSFit_;
    bool exactderivatives_;
    bool exactderivativesforfloating_;
    bool useapprox_;
    bool usephicritapprox_;

    unsigned int minIndStubs_;
    std::string removalType_;
    std::string mergeComparison_;
    bool doKF_;
    bool fakefit_;

    unsigned int nHelixPar_;
    bool extended_;

    std::string fitpatternfile_;
    std::string skimfile_;

  };
}

#endif

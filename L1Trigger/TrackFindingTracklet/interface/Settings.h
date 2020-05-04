#ifndef L1Trigger_TrackFindingTracklet_interface_Settings_h
#define L1Trigger_TrackFindingTracklet_interface_Settings_h

#include <iostream>
#include <string>
#include <array>
#include <set>
#include <map>
#include <cassert>
#include <cmath>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

namespace Trklet {

  class Settings {
  public:
    Settings() {
      //Uncomment to run the hybrid algorithm
      #ifdef CMSSW_GIT_HASH
      #define USEHYBRID
      #endif
    }

    // processing & memory modules, wiring, etc.
    std::string DTCLinkFile() const { return DTCLinkFile_; }
    std::string moduleCablingFile() const { return moduleCablingFile_; }
    std::string DTCLinkLayerDiskFile() const { return DTCLinkLayerDiskFile_; }
    std::string fitPatternFile() const { return fitPatternFile_; }
    std::string processingModulesFile() const { return processingModulesFile_; }
    std::string memoryModulesFile() const { return memoryModulesFile_; }
    std::string wiresFile() const { return wiresFile_; }

    void setDTCLinkFile(std::string DTCLinkFileName) { DTCLinkFile_ = DTCLinkFileName; }
    void setModuleCablingFile(std::string moduleCablingFileName) { moduleCablingFile_ = moduleCablingFileName; }
    void setDTCLinkLayerDiskFile(std::string DTCLinkLayerDiskFileName) {
      DTCLinkLayerDiskFile_ = DTCLinkLayerDiskFileName;
    }
    void setFitPatternFile(std::string fitPatternFileName) { fitPatternFile_ = fitPatternFileName; }
    void setProcessingModulesFile(std::string processingModulesFileName) {
      processingModulesFile_ = processingModulesFileName;
    }
    void setMemoryModulesFile(std::string memoryModulesFileName) { memoryModulesFile_ = memoryModulesFileName; }
    void setWiresFile(std::string wiresFileName) { wiresFile_ = wiresFileName; }

    unsigned int nzbitsstub(unsigned int layerdisk) const { return nzbitsstub_[layerdisk]; }
    unsigned int nphibitsstub(unsigned int layerdisk) const { return nphibitsstub_[layerdisk]; }
    unsigned int nrbitsstub(unsigned int layerdisk) const { return nrbitsstub_[layerdisk]; }

    unsigned int nrbitsprojderdisk() const { return nrbitsprojderdisk_; }
    unsigned int nbitsphiprojderL123() const { return nbitsphiprojderL123_; }
    unsigned int nbitsphiprojderL456() const { return nbitsphiprojderL456_; }
    unsigned int nbitszprojderL123() const { return nbitszprojderL123_; }
    unsigned int nbitszprojderL456() const { return nbitszprojderL456_; }

    bool useSeed(unsigned int iSeed) const { return useseeding_.find(iSeed) != useseeding_.end(); }
    unsigned int nbitsvmte(unsigned int inner, unsigned int iSeed) const { return nbitsvmte_[inner][iSeed]; }
    unsigned int nvmte(unsigned int inner, unsigned int iSeed) const { return (1 << nbitsvmte_[inner][iSeed]); }

    unsigned int nbitsvmme(unsigned int layerdisk) const { return nbitsvmme_[layerdisk]; }
    unsigned int nvmme(unsigned int layerdisk) const { return (1 << nbitsvmme_[layerdisk]); }

    unsigned int nbitsallstubs(unsigned int layerdisk) const { return nbitsallstubs_[layerdisk]; }
    unsigned int nallstubs(unsigned int layerdisk) const { return (1 << nbitsallstubs_[layerdisk]); }

    bool writeMonitorData(std::string module) const {
      if (writeMonitorData_.find(module) == writeMonitorData_.end()) {
        edm::LogPrint("Tracklet") << "Settings::writeMonitorData module = " << module << " not known";
        assert(0);
      }
      return writeMonitorData_.at(module);
    }

    unsigned int maxStep(std::string module) const {
      if (maxstep_.find(module) == maxstep_.end()) {
        edm::LogPrint("Tracklet") << "Settings::maxStep module = " << module << " not known";
        assert(0);
      }
      return maxstep_.at(module) + maxstepoffset_;
      ;
    }

    double zlength() const { return zlength_; }
    double rmaxdisk() const { return rmaxdisk_; }

    double drmax() const { return rmaxdisk_ / 32.0; }
    double dzmax() const { return zlength_ / 32.0; }

    double half2SmoduleWidth() const { return half2SmoduleWidth_; }

    double bendcutte(unsigned int inner, unsigned int iSeed) const { return bendcutte_[inner][iSeed]; }
    double bendcutme(unsigned int layerdisk) const { return bendcutme_[layerdisk]; }
    double nfinephi(unsigned int inner, unsigned int iSeed) const { return nfinephi_[inner][iSeed]; }
    double nphireg(unsigned int inner, unsigned int iSeed) const { return nphireg_[inner][iSeed]; }
    double zbitstab(unsigned int inner, unsigned int iSeed) const { return zbitstab_[inner][iSeed]; }
    double rbitstab(unsigned int inner, unsigned int iSeed) const { return rbitstab_[inner][iSeed]; }
    double lutwidthtab(unsigned int inner, unsigned int iSeed) const { return lutwidthtab_[inner][iSeed]; }
    double lutwidthtabextended(unsigned int inner, unsigned int iSeed) const {
      return lutwidthtabextended_[inner][iSeed];
    }

    double projlayers(unsigned int iSeed, unsigned int i) const { return projlayers_[iSeed][i]; }
    double projdisks(unsigned int iSeed, unsigned int i) const { return projdisks_[iSeed][i]; }
    double rphimatchcut(unsigned int iSeed, unsigned int ilayer) const { return rphimatchcut_[ilayer][iSeed]; }
    double zmatchcut(unsigned int iSeed, unsigned int ilayer) const { return zmatchcut_[ilayer][iSeed]; }
    double rphicutPS(unsigned int iSeed, unsigned int idisk) const { return rphicutPS_[idisk][iSeed]; }
    double rcutPS(unsigned int iSeed, unsigned int idisk) const { return rcutPS_[idisk][iSeed]; }
    double rphicut2S(unsigned int iSeed, unsigned int idisk) const { return rphicut2S_[idisk][iSeed]; }
    double rcut2S(unsigned int iSeed, unsigned int idisk) const { return rcut2S_[idisk][iSeed]; }

    double rmean(unsigned int iLayer) const { return irmean_[iLayer] * rmaxdisk_ / 4096; }
    double rmax(unsigned int iLayer) const { return rmean(iLayer) + drmax(); }
    double rmin(unsigned int iLayer) const { return rmean(iLayer) - drmax(); }
    double zmean(unsigned int iDisk) const { return izmean_[iDisk] * zlength_ / 2048; }
    double zmax(unsigned int iDisk) const { return zmean(iDisk) + dzmax(); }
    double zmin(unsigned int iDisk) const { return zmean(iDisk) - dzmax(); }

    double rDSSinner(unsigned int iBin) const {
      return rDSSinner_mod_[iBin / 2] + halfstrip_ * ((iBin % 2 == 0) ? -1 : 1);
    }
    double rDSSouter(unsigned int iBin) const {
      return rDSSouter_mod_[iBin / 2] + halfstrip_ * ((iBin % 2 == 0) ? -1 : 1);
    }

    bool printDebugKF() const { return printDebugKF_; }
    bool debugTracklet() const { return debugTracklet_; }
    bool writetrace() const { return writetrace_; }

    bool warnNoMem() const { return warnNoMem_; }
    bool warnNoDer() const { return warnNoDer_; }

    bool writeMem() const { return writeMem_; }
    bool writeTable() const { return writeTable_; }

    bool writeVerilog() const { return writeVerilog_; }
    bool writeHLS() const { return writeHLS_; }
    bool writeInvTable() const { return writeInvTable_; }
    bool writeHLSInvTable() const { return writeHLSInvTable_; }

    unsigned int writememsect() const { return writememsect_; }

    bool writeTripletTables() const { return writeTripletTables_; }

    bool writeoutReal() const { return writeoutReal_; }

    bool bookHistos() const { return bookHistos_; }

    double ptcut() const { return ptcut_; }
    double rinvcut() const { return 0.01 * c_ * bfield_ / ptcut_; }  //0.01 to convert to cm-1

    double c() const { return c_; }
    double bfield() const { return bfield_; }

    double rinvmax() const { return 0.01 * c_ * bfield_ / ptmin_; }

    int alphashift() const { return alphashift_; }
    int nbitsalpha() const { return nbitsalpha_; }
    int alphaBitsTable() const { return alphaBitsTable_; }
    int nrinvBitsTable() const { return nrinvBitsTable_; }

    unsigned int MEBinsBits() const { return MEBinsBits_; }
    unsigned int MEBins() const { return 1u << MEBinsBits_; }
    unsigned int MEBinsDisks() const { return MEBinsDisks_; }

    std::string geomext() const { return extended_ ? "hourglassExtended" : "hourglass"; }

    bool useMSFit() const { return useMSFit_; }
    bool exactderivatives() const { return exactderivatives_; }
    bool exactderivativesforfloating() const { return exactderivativesforfloating_; }
    bool useapprox() const { return useapprox_; }
    bool usephicritapprox() const { return usephicritapprox_; }

    unsigned int minIndStubs() const { return minIndStubs_; }
    std::string removalType() const { return removalType_; }
    std::string mergeComparison() const { return mergeComparison_; }
    bool doKF() const { return doKF_; }
    bool fakefit() const { return fakefit_; }

    // configurable
    unsigned int nHelixPar() const { return nHelixPar_; }
    void setNHelixPar(unsigned int nHelixPar) { nHelixPar_ = nHelixPar; }

    bool extended() const { return extended_; }
    void setExtended(bool extended) { extended_ = extended; }

    std::string skimfile() const { return skimfile_; }
    void setSkimfile(std::string skimfile) { skimfile_ = skimfile; }

    double dphisectorHG() const {
      return 2 * M_PI / NSector_ +
             2 * fmax(std::abs(asin(0.5 * rinvmax() * rmean(0)) - asin(0.5 * rinvmax() * rcrit_)),
                      std::abs(asin(0.5 * rinvmax() * rmean(5)) - asin(0.5 * rinvmax() * rcrit_)));
    }

    double rcrit() const { return rcrit_; }

    double dphisector() const { return 2 * M_PI / NSector_; }

    unsigned int NSector() const { return NSector_; }

    double phicritmin() const { return 0.5 * dphisectorHG() - M_PI / NSector_; }
    double phicritmax() const { return dphisectorHG() - 0.5 * dphisectorHG() + M_PI / NSector_; }

    double phicritminmc() const { return phicritmin() - dphicritmc_; }
    double phicritmaxmc() const { return phicritmax() + dphicritmc_; }

    double kphi() const { return dphisectorHG() / (1 << nphibitsstub(0)); }
    double kphi1() const { return dphisectorHG() / (1 << nphibitsstub(5)); }

    double kz() const { return 2 * zlength_ / (1 << nzbitsstub_[0]); }
    double kr() const { return rmaxdisk_ / (1 << nrbitsstub_[6]); }

    double maxrinv() const { return maxrinv_; }
    double maxd0() const { return maxd0_; }
    unsigned int nbitsd0() const { return nbitsd0_; }

    double kd0() const { return 2 * maxd0_ / (1 << nbitsd0_); }

    double rinvcutte() const { return 0.01 * c_ * bfield_ / ptcutte_; }  //0.01 to convert to cm-1

    double rmindiskvm() const { return rmindiskvm_; }
    double rmaxdiskvm() const { return rmaxdiskvm_; }

    double rmaxdiskl1overlapvm() const { return rmaxdiskl1overlapvm_; }
    double rmindiskl2overlapvm() const { return rmindiskl2overlapvm_; }
    double rmindiskl3overlapvm() const { return rmindiskl3overlapvm_; }

    double z0cut() const { return z0cut_; }

    // Obsolete - only used in TrackletCalculatorDisplaced (Ryd - 2020-01-16)
    int iphicritminmc() const { return iphicritminmc_; }
    int iphicritmaxmc() const { return iphicritmaxmc_; }

    unsigned int NLONGVMBITS() const { return NLONGVMBITS_; }
    unsigned int NLONGVMBINS() const { return (1 << NLONGVMBITS_); }

    unsigned int ntrackletmax() const { return ntrackletmax_; }

    //Bits used to store track parameter in tracklet
    int nbitsrinv() const { return nbitsrinv_; }
    int nbitsphi0() const { return nbitsphi0_; }
    int nbitst() const { return nbitst_; }
    int nbitsz0() const { return nbitsz0_; }

    //track and tracklet parameters
    int rinv_shift() const { return rinv_shift_; }
    int phi0_shift() const { return phi0_shift_; }
    int t_shift() const { return t_shift_; }
    int z0_shift() const { return z0_shift_; }

    //projections are coarsened from global to stub precision

    //projection to R parameters
    int SS_phiL_shift() const { return SS_phiL_shift_; }
    int PS_zL_shift() const { return PS_zL_shift_; }

    int SS_phiderL_shift() const { return SS_phiderL_shift_; }
    int PS_zderL_shift() const { return PS_zderL_shift_; }
    int SS_zderL_shift() const { return SS_zderL_shift_; }

    //projection to Z parameters
    int SS_phiD_shift() const { return SS_phiD_shift_; }
    int PS_rD_shift() const { return PS_rD_shift_; }

    int SS_phiderD_shift() const { return SS_phiderD_shift_; }
    int PS_rderD_shift() const { return PS_rderD_shift_; }

    //numbers needed for matches & fit, unclear what they are.
    int phi0bitshift() const { return phi0bitshift_; }
    int phiderbitshift() const { return phiderbitshift_; }
    int zderbitshift() const { return zderbitshift_; }

    int phiresidbits() const { return phiresidbits_; }
    int zresidbits() const { return zresidbits_; }
    int rresidbits() const { return rresidbits_; }

    //Trackfit
    int fitrinvbitshift() const { return fitrinvbitshift_; }
    int fitphi0bitshift() const { return fitphi0bitshift_; }
    int fittbitshift() const { return fittbitshift_; }
    int fitz0bitshift() const { return fitz0bitshift_; }

    //r correction bits
    int rcorrbits() const { return rcorrbits_; }

    int chisqphifactbits() const { return chisqphifactbits_; }
    int chisqzfactbits() const { return chisqzfactbits_; }

    //should not return reference...
    double& krinvpars() const { return krinvpars_; }
    double& kphi0pars() const { return kphi0pars_; }
    double& kd0pars() const { return kd0pars_; }
    double& ktpars() const { return ktpars_; }
    double& kz0pars() const { return kz0pars_; }
    double& kphiproj123() const { return kphiproj123_; }
    double& kzproj() const { return kzproj_; }
    double& kphider() const { return kphider_; }
    double& kzder() const { return kzder_; }
    double& krprojshiftdisk() const { return krprojshiftdisk_; }
    double& kphiprojdisk() const { return kphiprojdisk_; }
    double& krdisk() const { return krdisk_; }
    double& kzpars() const { return kzpars_; }

  private:
    std::string DTCLinkFile_;
    std::string moduleCablingFile_;
    std::string DTCLinkLayerDiskFile_;
    std::string fitPatternFile_;
    std::string processingModulesFile_;
    std::string memoryModulesFile_;
    std::string wiresFile_;

    unsigned int NSector_{9};

    double rcrit_{55.0};

    double dphicritmc_{0.005};

    std::array<unsigned int, 6> irmean_{{851, 1269, 1784, 2347, 2936, 3697}};
    std::array<unsigned int, 5> izmean_{{2239, 2645, 3163, 3782, 4523}};

    std::array<unsigned int, 11> nzbitsstub_{{12, 12, 12, 8, 8, 8, 7, 7, 7, 7, 7}};
    std::array<unsigned int, 11> nphibitsstub_{{14, 14, 14, 17, 17, 17, 14, 14, 14, 14, 14}};
    std::array<unsigned int, 11> nrbitsstub_{{7, 7, 7, 7, 7, 7, 12, 12, 12, 12, 12}};

    unsigned int nrbitsprojderdisk_{9};
    unsigned int nbitsphiprojderL123_{10};
    unsigned int nbitsphiprojderL456_{10};
    unsigned int nbitszprojderL123_{10};
    unsigned int nbitszprojderL456_{9};

    std::set<unsigned int> useseeding_{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

    std::array<unsigned int, 11> nbitsallstubs_{{3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}};
    std::array<unsigned int, 11> nbitsvmme_{{2, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2}};
    std::array<std::array<unsigned int, 12>, 3> nbitsvmte_{{{{2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 3, 2}},
                                                            {{3, 2, 3, 3, 2, 2, 2, 2, 3, 3, 2, 2}},
                                                            {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1}}}};

    std::array<std::array<double, 8>, 2> bendcutte_{{{{1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25}},    //inner
                                                     {{1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25}}}};  //outer

    std::array<double, 11> bendcutme_{{2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.5, 1.5, 1.5, 1.5, 1.5}};

    double rmindiskvm_{22.5};
    double rmaxdiskvm_{67.0};

    double rmaxdiskl1overlapvm_{45.0};
    double rmindiskl2overlapvm_{40.0};
    double rmindiskl3overlapvm_{50.0};

    double z0cut_{15.0};

    // Obsolete - only used in TrackletCalculatorDisplaced (Ryd - 2020-01-16)
    int iphicritminmc_{9253};
    int iphicritmaxmc_{56269};

    unsigned int NLONGVMBITS_{3};

    double zlength_{120.0};
    double rmaxdisk_{120.0};

    double half2SmoduleWidth_{4.57};

    double maxrinv_{0.006};
    double maxd0_{10.0};

    unsigned int nbitsd0_{13};

    double ptmin_{2.0};  //minumim pt for tracks

    double ptcutte_{1.8};  //Minimum pt in TE

    unsigned int ntrackletmax_{127};  //maximum number of tracklets that can be stored

    //Bits used to store track parameter in tracklet
    int nbitsrinv_{14};
    int nbitsphi0_{18};
    int nbitst_{14};
    int nbitsz0_{10};

    //track and tracklet parameters
    int rinv_shift_{-8};  // Krinv = 2^shift * Kphi/Kr
    int phi0_shift_{1};   // Kphi0 = 2^shift * Kphi
    int t_shift_{-10};    // Kt    = 2^shift * Kz/Kr
    int z0_shift_{0};     // Kz0   = 2^shift * kz

    //projections are coarsened from global to stub precision

    //projection to R parameters
    int SS_phiL_shift_{0};
    int PS_zL_shift_{0};  // z projections have global precision in ITC

    int SS_phiderL_shift_{-5};
    int PS_zderL_shift_{-7};  // Kderz = 2^shift * Kz/Kr
    int SS_zderL_shift_{-7};

    //projection to Z parameters
    int SS_phiD_shift_{3};
    int PS_rD_shift_{1};  // a bug?! coarser by a factor of two then stubs??

    int SS_phiderD_shift_{-4};
    int PS_rderD_shift_{-6};  //Kderrdisk = 2^shift * Kr/Kz

    //numbers needed for matches & fit, unclear what they are.
    int phi0bitshift_{1};
    int phiderbitshift_{7};
    int zderbitshift_{6};

    int phiresidbits_{12};
    int zresidbits_{9};
    int rresidbits_{7};

    //Trackfit
    int fitrinvbitshift_{9};  //6 OK?
    int fitphi0bitshift_{6};  //4 OK?
    int fittbitshift_{10};    //4 OK? //lower number gives rounding problems
    int fitz0bitshift_{8};    //6 OK?

    //r correction bits
    int rcorrbits_{6};

    int chisqphifactbits_{14};
    int chisqzfactbits_{14};

    std::array<std::array<unsigned int, 12>, 3> nfinephi_{{{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}},    //inner
                                                           {{3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}},    //outer
                                                           {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3}}}};  //outermost

    //These are the number of bits used for the VM regions in the TE by seedindex
    std::array<std::array<unsigned int, 12>, 3> nphireg_{{{{5, 4, 4, 4, 4, 4, 4, 3, 4, 4, 5, 4}},    //inner
                                                          {{5, 4, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4}},    //outer
                                                          {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4}}}};  //outermost

    std::array<std::array<unsigned int, 12>, 3> zbitstab_{{{{7, 7, 7, 7, 3, 3, 7, 7, 0, 0, 7, 0}},
                                                           {{7, 7, 7, 7, 3, 3, 3, 3, 0, 0, 7, 0}},
                                                           {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 7}}}};

    std::array<std::array<unsigned int, 12>, 3> rbitstab_{{{{4, 4, 4, 4, 8, 8, 3, 3, 0, 0, 4, 0}},
                                                           {{4, 4, 4, 4, 7, 7, 7, 7, 0, 0, 4, 0}},
                                                           {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 4}}}};

    std::array<std::array<unsigned int, 12>, 3> lutwidthtab_{{{{10, 11, 11, 11, 11, 11, 11, 11, 0, 0, 11, 0}},
                                                              {{6, 6, 6, 6, 10, 10, 10, 10, 0, 0, 6, 0}},
                                                              {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 6}}}};

    std::array<std::array<unsigned int, 12>, 3> lutwidthtabextended_{{{{11, 11, 21, 21, 21, 21, 11, 11, 0, 0, 21, 0}},
                                                                      {{6, 6, 6, 6, 10, 10, 10, 10, 0, 0, 6, 0}},
                                                                      {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 6}}}};

    //projection layers by seed index. For each seeding index (row) the list of layers that we consider projections to
    std::array<std::array<unsigned int, 4>, 12> projlayers_{{{{3, 4, 5, 6}},  //0 L1L2
                                                             {{1, 4, 5, 6}},  //1 L2L3
                                                             {{1, 2, 5, 6}},  //2 L3L4
                                                             {{1, 2, 3, 4}},  //3 L5L6
                                                             {{1, 2}},        //4 D1D2
                                                             {{1}},           //5 D3D4
                                                             {{}},            //6 L1D1
                                                             {{1}},           //7 L2D1
                                                             {{1, 5, 6}},     //8 L2L3L4
                                                             {{1, 2, 3}},     //9 L4L5L6
                                                             {{1}},           //10 L2L3D1
                                                             {{1}}}};         //11 D1D2L2

    //projection disks by seed index. For each seeding index (row) the list of diks that we consider projections to
    std::array<std::array<unsigned int, 5>, 12> projdisks_{{{{1, 2, 3, 4}},  //0 L1L2
                                                            {{1, 2, 3, 4}},  //1 L2L3
                                                            {{1, 2}},        //2 L3L4
                                                            {{}},            //3 L5L6
                                                            {{3, 4, 5}},     //4 D1D2
                                                            {{1, 2, 5}},     //5 D3D4
                                                            {{2, 3, 4, 5}},  //6 L1D1
                                                            {{2, 3, 4}},     //7 L2D1
                                                            {{1, 2}},        //8 L2L3L4
                                                            {{}},            //9 L4L5L6
                                                            {{2, 3, 4}},     //10 L2L3D1
                                                            {{3, 4}}}};      //11 D1D2L2

    //rphi cuts for layers - the column is the seedindex
    std::array<std::array<double, 12>, 6> rphimatchcut_{
        {{{0.0, 0.1, 0.07, 0.08, 0.07, 0.05, 0.0, 0.05, 0.08, 0.15, 0.125, 0.15}},  //Layer 1
         {{0.0, 0.0, 0.06, 0.08, 0.05, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0}},         //Layer 2
         {{0.1, 0.0, 0.0, 0.08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.08, 0.0, 0.0}},          //Layer 3
         {{0.19, 0.19, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},         //Layer 4
         {{0.4, 0.4, 0.08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.08, 0.0, 0.0, 0.0}},          //Layer 5
         {{0.5, 0.0, 0.19, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0}}}};         //Layer 6

    //z cuts for layers - the column is the seedindex
    std::array<std::array<double, 12>, 6> zmatchcut_{
        {{{0.0, 0.7, 5.5, 15.0, 1.5, 2.0, 0.0, 1.5, 1.0, 8.0, 1.0, 1.5}},   //Layer 1
         {{0.0, 0.0, 3.5, 15.0, 1.25, 0.0, 0.0, 0.0, 0.0, 7.0, 0.0, 0.0}},  //Layer 2
         {{0.7, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0}},    //Layer 3
         {{3.0, 3.0, 0.0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},    //Layer 4
         {{3.0, 3.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.5, 0.0, 0.0, 0.0}},    //Layer 5
         {{4.0, 0.0, 9.5, 0.0, 0.0, 0.0, 0.0, 0.0, 4.5, 0.0, 0.0, 0.0}}}};  //Layer 6

    //rphi cuts for PS modules in disks - the column is the seedindex
    std::array<std::array<double, 12>, 5> rphicutPS_{
        {{{0.2, 0.2, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},     //disk 1
         {{0.2, 0.2, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.0, 0.0, 0.15, 0.0}},    //disk 2
         {{0.25, 0.2, 0.0, 0.0, 0.15, 0.0, 0.2, 0.15, 0.0, 0.0, 0.0, 0.2}},  //disk 3
         {{0.5, 0.2, 0.0, 0.0, 0.2, 0.0, 0.3, 0.5, 0.0, 0.0, 0.0, 0.0}},     //disk 4
         {{0.0, 0.0, 0.0, 0.0, 0.25, 0.1, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0}}}};  //disk 5

    //r cuts for PS modules in disks - the column is the seedindex
    std::array<std::array<double, 12>, 5> rcutPS_{
        {{{0.5, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},    //disk 1
         {{0.5, 0.5, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.0, 0.0, 0.5, 0.0}},    //disk 2
         {{0.5, 0.5, 0.0, 0.0, 0.5, 0.0, 0.6, 0.8, 0.0, 0.0, 0.0, 0.4}},    //disk 3
         {{0.5, 0.5, 0.0, 0.0, 0.8, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0}},    //disk 4
         {{0.0, 0.0, 0.0, 0.0, 1.0, 0.5, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0}}}};  //disk 5

    //rphi cuts for 2S modules in disks = the column is the seedindex
    std::array<std::array<double, 12>, 5> rphicut2S_{
        {{{0.5, 0.5, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0}},    //disk 1
         {{0.5, 0.5, 0.8, 0.0, 0.0, 0.0, 0.5, 0.15, 0.3, 0.0, 0.68, 0.0}},  //disk 2
         {{0.5, 0.5, 0.0, 0.0, 0.15, 0.0, 0.2, 0.25, 0.0, 0.0, 0.8, 0.1}},  //disk 3
         {{0.5, 0.5, 0.0, 0.0, 0.2, 0.0, 0.25, 0.5, 0.0, 0.0, 0.6, 0.4}},   //disk 4
         {{0.0, 0.0, 0.0, 0.0, 0.4, 0.2, 0.4, 0.0, 0.0, 0.0, 0.0, 0.8}}}};  //disk 5

    //r cuts for 2S modules in disks -the column is the seedindex
    std::array<std::array<double, 12>, 5> rcut2S_{
        {{{3.8, 3.8, 3.8, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0}},    //disk 1
         {{3.8, 3.8, 3.8, 0.0, 0.0, 0.0, 3.8, 3.4, 3.0, 0.0, 3.0, 0.0}},    //disk 2
         {{3.6, 3.8, 0.0, 0.0, 3.6, 0.0, 3.6, 3.8, 0.0, 0.0, 3.8, 3.0}},    //disk 3
         {{3.6, 3.8, 0.0, 0.0, 3.6, 0.0, 3.5, 3.8, 0.0, 0.0, 3.0, 3.0}},    //disk 4
         {{0.0, 0.0, 0.0, 0.0, 3.6, 3.4, 3.7, 0.0, 0.0, 0.0, 0.0, 3.0}}}};  //disk 5

    unsigned int maxstepoffset_{10000};

    std::map<std::string, unsigned int> maxstep_{{"Link", 108},
                                                 {"MC", 108},
                                                 {"ME", 108},
                                                 {"MP", 108},
                                                 {"PR", 108},
                                                 {"TC", 108},
                                                 {"TE", 108},
                                                 {"TP", 108},
                                                 {"TRE", 108},
                                                 {"VMR", 108}};

    std::map<std::string, bool> writeMonitorData_{{"IL", false},           {"TE", false},
                                                  {"CT", false},           {"HitPattern", false},
                                                  {"ChiSq", false},        {"Seeds", false},
                                                  {"FT", false},           {"Residuals", false},
                                                  {"MC", false},           {"ME", false},
                                                  {"AP", false},           {"VMP", false},
                                                  {"NMatches", false},     {"TrackProjOcc", false},
                                                  {"TC", false},           {"Pars", false},
                                                  {"TPars", false},        {"TPD", false},
                                                  {"TrackletPars", false}, {"TED", false},
                                                  {"TP", false},           {"TRE", false},
                                                  {"VMR", false},          {"Variance", false},
                                                  {"StubsLayer", false},   {"StubsLayerSector", false},
                                                  {"ResEff", false},       {"HitEff", false},
                                                  {"MatchEff", false},     {"Cabling", false},
                                                  {"IFit", false},         {"AS", false}};

    std::array<double, 5> rDSSinner_mod_{{68.9391, 78.7750, 85.4550, 96.3150, 102.3160}};

    std::array<double, 5> rDSSouter_mod_{{66.4903, 76.7750, 84.4562, 94.9920, 102.3160}};

    double halfstrip_{
        2.5};  //we want the center of the two strip positions in a module, not just the center of a module

    // various printouts for debugging and warnings
    bool printDebugKF_{false};   // if true print lots of debugging statements related to the KF fit
    bool debugTracklet_{false};  //Print detailed debug information about tracklet tracking
    bool writetrace_{false};     //Print out details about parsing configuration files

    bool warnNoMem_{false};  //If true will print out warnings about missing projection memories
    bool warnNoDer_{false};  //If true will print out warnings about missing track fit derivatives

    bool writeMem_{false};    //If true will print out content of memories to files
    bool writeTable_{false};  //IF true will print out content of LUTs to files

    // Write various lookup tables and autogenerated code (from iMath)
    bool writeVerilog_{false};      //Write out auto-generated Verilog mudules used by TCs
    bool writeHLS_{false};          //Write out auto-generated HLS mudules used by TCs
    bool writeInvTable_{false};     //Write out tables of drinv and invt in tracklet calculator for Verilog module
    bool writeHLSInvTable_{false};  //Write out tables of drinv and invt in tracklet calculator for HLS module

    unsigned int writememsect_{3};  //writemem only for this sector (note that the files will have _4 extension)

    bool writeTripletTables_{false};  //Train and write the TED and TRE tables. N.B.: the tables
                                      //cannot be applied while they are being trained, i.e.,
                                      //this flag effectively turns off the cuts in
                                      //TrackletEngineDisplaced and TripletEngine

    bool writeoutReal_{false};

    bool bookHistos_{
        false};  //set to true/false to turn on/off histogram booking internal to the tracking (class "HistImp")

    // pt constants
    double ptcut_{1.91};  //Minimum pt cut

    // Parameters for bit sizes
    int alphashift_{12};
    int nbitsalpha_{4};      //bits used to store alpha
    int alphaBitsTable_{2};  //For number of bits in track derivative table
    int nrinvBitsTable_{3};  //number of bits for tabulating rinv dependence

    unsigned int MEBinsBits_{3};
    unsigned int MEBinsDisks_{8};  //on each side

    // Options for chisq fit
    bool useMSFit_{false};  //for both the integer and float
    bool exactderivatives_{false};
    bool exactderivativesforfloating_{true};  //only for the floating point
    bool useapprox_{true};          //use approximate postion based on integer representation for floating point
    bool usephicritapprox_{false};  //use floating point approximate version of phicrit cut if true

    // Duplicate Removal
    // "merge" (hybrid dup removal)
    // "ichi" (pairwise, keep track with best ichisq), "nstub" (pairwise, keep track with more stubs)
    // "grid" (TMTT-like removal), "" (no removal)
    unsigned int minIndStubs_{3};  // not used with merge removal

#ifdef USEHYBRID
    std::string removalType_{"merge"};
    // "CompareBest" (recommended) Compares only the best stub in each track for each region (best = smallest phi residual)
    // and will merge the two tracks if stubs are shared in three or more regions
    // "CompareAll" Compares all stubs in a region, looking for matches, and will merge the two tracks if stubs are shared in three or more regions
    std::string mergeComparison_{"CompareBest"};
    bool doKF_{true};
#endif

#ifndef USEHYBRID
    bool doKF_{false};
    std::string removalType_{"ichi"};
    std::string mergeComparison_{""};
#endif

    bool fakefit_{
        false};  // if true, run a dummy fit, producing TTracks directly from output of tracklet pattern reco stage

    unsigned int nHelixPar_{4};  // 4 or 5 param helix fit
    bool extended_{false};       // turn on displaced tracking

    std::string skimfile_{""};  //if not empty events will be written out in ascii format to this file

    //constants derivative from the above - TODO should be calculated in Settings.h, not set externally
    //then we can remove the mutable

    double bfield_{3.8};  //B-field in T
    double c_{0.3};       //speed of light m/ns

    mutable double krinvpars_;
    mutable double kphi0pars_;
    mutable double kd0pars_;
    mutable double ktpars_;
    mutable double kz0pars_;
    mutable double kphiproj123_;
    mutable double kzproj_;
    mutable double kphider_;
    mutable double kzder_;
    mutable double krprojshiftdisk_;
    mutable double kphiprojdisk_;
    mutable double krdisk_;
    mutable double kzpars_;
  };
}  // namespace Trklet

#endif

#include "L1Trigger/TrackFindingTracklet/interface/MatchProcessor.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"
#include "L1Trigger/TrackFindingTracklet/interface/Util.h"
#include "L1Trigger/TrackFindingTracklet/interface/ProjectionRouterBendTable.h"
#include "L1Trigger/TrackFindingTracklet/interface/HistBase.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include <filesystem>

using namespace std;
using namespace trklet;

MatchProcessor::MatchProcessor(string name, Settings const& settings, Globals* global, unsigned int iSector)
    : ProcessBase(name, settings, global, iSector), fullmatches_(12), inputProjBuffer_(3) {

  phiregion_ = name[8] - 'A';

  layerdisk_ = initLayerDisk(3);

  barrel_ = layerdisk_ < N_LAYER;
  
  phishift_ = settings_.nphibitsstub(N_LAYER-1)-settings_.nphibitsstub(layerdisk_);
  dzshift_ = settings_.nzbitsstub(0) - settings_.nzbitsstub(layerdisk_);

  if (barrel_) {
    icorrshift_=ilog2(settings_.kphi(layerdisk_)/(settings_.krbarrel()*settings_.kphider()));
    icorzshift_=ilog2(settings_.kz(layerdisk_)/(settings_.krbarrel()*settings_.kzder()));
  } else {
    icorrshift_=ilog2(settings_.kphi(layerdisk_)/(settings_.kz()*settings_.kphiderdisk()));
    icorzshift_=ilog2(settings_.krprojshiftdisk()/(settings_.kz()*settings_.krder()));
  }
   
  nrbits_ = 5;
  nphiderbits_ = 6;

  nrinv_=NRINVBITS;
  double rinvhalf=0.5*((1<<nrinv_)-1);

  for (unsigned int iSeed = 0; iSeed < 12; iSeed++) {
    if (layerdisk_ < N_LAYER) {
      phimatchcut_[iSeed] =
          settings_.rphimatchcut(iSeed, layerdisk_) / (settings_.kphi1() * settings_.rmean(layerdisk_));
      zmatchcut_[iSeed] = settings_.zmatchcut(iSeed, layerdisk_) / settings_.kz();
    } else {
      rphicutPS_[iSeed] = settings_.rphicutPS(iSeed, layerdisk_ - N_LAYER) / (settings_.kphi() * settings_.kr());
      rphicut2S_[iSeed] = settings_.rphicut2S(iSeed, layerdisk_ - N_LAYER) / (settings_.kphi() * settings_.kr());
      rcut2S_[iSeed] = settings_.rcut2S(iSeed, layerdisk_ - N_LAYER) / settings_.krprojshiftdisk();
      rcutPS_[iSeed] = settings_.rcutPS(iSeed, layerdisk_ - N_LAYER) / settings_.krprojshiftdisk();
    }
  }

  if (iSector_ == 0 && barrel_ && settings_.writeTable()) {

    ofstream outphicut=openfile(settings_.tablePath(), getName() + "_phicut.tab", __FILE__, __LINE__);

    outphicut << "{" << endl;
    for (unsigned int seedindex = 0; seedindex < 12; seedindex++) {
      if (seedindex != 0)
        outphicut << "," << endl;
      outphicut << phimatchcut_[seedindex];
    }
    outphicut << endl << "};" << endl;
    outphicut.close();

    ofstream outzcut=openfile(settings_.tablePath(), getName() + "_zcut.tab", __FILE__, __LINE__);

    outzcut << "{" << endl;
    for (unsigned int seedindex = 0; seedindex < N_SEED; seedindex++) {
      if (seedindex != 0)
        outzcut << "," << endl;
      outzcut << zmatchcut_[seedindex];
    }
    outzcut << endl << "};" << endl;
    outzcut.close();
  }

  if (barrel_) {
    unsigned int nbits = (layerdisk_ < N_PSLAYER) ? N_BENDBITS_PS : N_BENDBITS_2S;
    
    for (unsigned int irinv = 0; irinv < (1u<<nrinv_); irinv++) {
      double rinv = (irinv - rinvhalf) * (1 << (settings_.nbitsrinv() - nrinv_)) * settings_.krinvpars();
      double stripPitch = settings_.stripPitch(layerdisk_<N_PSLAYER);
      double projbend = bendstrip(settings_.rmean(layerdisk_), rinv, stripPitch);
      for (unsigned int ibend = 0; ibend < (unsigned int)(1 << nbits); ibend++) {
        double stubbend = settings_.benddecode(ibend, layerdisk_, layerdisk_ < (int)N_PSLAYER);
        bool pass = std::abs(stubbend - projbend) < settings_.bendcutme(ibend, layerdisk_, layerdisk_ < (int)N_PSLAYER);
        table_.push_back(pass);
      }
    }
  } else {
    table_.resize(32*16*2,false);
    for (unsigned int iprojbend = 0; iprojbend < 32; iprojbend++) {
      double projbend = 0.5 * (iprojbend - rinvhalf);
      for (unsigned int ibend = 0; ibend < 8; ibend++) {
        double stubbend = settings_.benddecode(ibend, layerdisk_, true);
        bool pass = std::abs(stubbend - projbend) < settings_.bendcutme(ibend, layerdisk_, true);
	table_[512+8*iprojbend+ibend]=pass;
      }
      for (unsigned int ibend = 0; ibend < 16; ibend++) {
        double stubbend = settings_.benddecode(ibend, layerdisk_, false);
        bool pass = std::abs(stubbend - projbend) < settings_.bendcutme(ibend, layerdisk_, false);
	table_[16*iprojbend+ibend]=pass;
      }
    }
  }
  
  if (settings_.writeTable()) {

    char layerdisk = barrel_ ? '0' + layerdisk_ + 1 : '0' + layerdisk_ - N_LAYER + 1;
    string fname = barrel_ ? "METable_L" : "METable_D";
    fname += layerdisk;
    fname += ".tab";

    ofstream out=openfile(settings_.tablePath(), fname, __FILE__, __LINE__);
    out << "{" << endl;
    for (unsigned int i = 0; i < table_.size(); i++) {
      if (i != 0) {
	out << "," << endl;
      }
      out << table_[i];
    }
    out << "};" << endl;
    out.close();
  }

  
  for (unsigned int i = 0; i < N_DSS_MOD * 2; i++) {
    ialphafactinner_[i] = (1 << settings_.alphashift()) * settings_.krprojshiftdisk() * settings_.half2SmoduleWidth() /
                          (1 << (settings_.nbitsalpha() - 1)) / (settings_.rDSSinner(i) * settings_.rDSSinner(i)) /
                          settings_.kphi();
    ialphafactouter_[i] = (1 << settings_.alphashift()) * settings_.krprojshiftdisk() * settings_.half2SmoduleWidth() /
                          (1 << (settings_.nbitsalpha() - 1)) / (settings_.rDSSouter(i) * settings_.rDSSouter(i)) /
                          settings_.kphi();
  }

  nvm_ = settings_.nvmme(layerdisk_) * settings_.nallstubs(layerdisk_);
  nvmbins_ = settings_.nvmme(layerdisk_);
  nvmbits_ = settings_.nbitsvmme(layerdisk_)+settings_.nbitsallstubs(layerdisk_);

  nMatchEngines_ = 4;
  for (unsigned int iME = 0; iME < nMatchEngines_; iME++) {
    MatchEngineUnit tmpME(barrel_, layerdisk_, table_);
    matchengines_.push_back(tmpME);
  }

  if (globals_->projectionRouterBendTable() == nullptr) { 
    auto* bendTablePtr = new ProjectionRouterBendTable();
    bendTablePtr->init(settings_, globals_, nrbits_, nphiderbits_);
    globals_->projectionRouterBendTable() = bendTablePtr;
  }
  
}

void MatchProcessor::addOutput(MemoryBase* memory, string output) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding output to " << memory->getName() << " to output "
                                 << output;
  }
  if (output.find("matchout") != std::string::npos) {
    auto* tmp = dynamic_cast<FullMatchMemory*>(memory);
    assert(tmp != nullptr);
    unsigned int iSeed = getISeed(tmp->getName());
    assert(iSeed < fullmatches_.size());
    assert(fullmatches_[iSeed] == nullptr);
    fullmatches_[iSeed] = tmp;
    return;
  }
  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " could not find output: " << output;
}

void MatchProcessor::addInput(MemoryBase* memory, string input) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding input from " << memory->getName() << " to input "
                                 << input;
  }
  if (input == "allstubin") {
    auto* tmp = dynamic_cast<AllStubsMemory*>(memory);
    assert(tmp != nullptr);
    allstubs_ = tmp;
    return;
  }
  if (input == "vmstubin") {
    auto* tmp = dynamic_cast<VMStubsMEMemory*>(memory);
    assert(tmp != nullptr);
    vmstubs_.push_back(tmp);  //to allow more than one stub in?  vmstubs_=tmp;
    return;
  }
  if (input == "projin") {
    auto* tmp = dynamic_cast<TrackletProjectionsMemory*>(memory);
    assert(tmp != nullptr);
    inputprojs_.push_back(tmp);
    return;
  }
  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " could not find input: " << input;
}

void MatchProcessor::execute() {
  assert(vmstubs_.size() == 1);

  /*
    The code is organized in three 'steps' corresponding to the PR, ME, and MC functions. The output from
    the PR step is buffered in a 'circular' buffer, and similarly the ME output is put in a circular buffer. 
    
    The implementation is done in steps, emulating what can be done in firmware. One each step we do:
    
    1) A projection is read and if there is space it is insert into the inputProjBuffer_
    
    2) Process next match in the ME - if there is an idle ME the next projection is inserted
    
    3) Readout match from ME and send to match calculator
    
  */

  //bool print = getName()=="MP_L3PHIC" && iSector_==3;

  Tracklet* oldTracklet = nullptr;

  unsigned int countme = 0;
  unsigned int countall = 0;
  unsigned int countsel = 0;
  unsigned int countinputproj = 0;

  unsigned int iprojmem = 0;
  unsigned int iproj = 0;

  inputProjBuffer_.reset();

  for (unsigned int i = 0; i < inputprojs_.size(); i++) {
    countinputproj += inputprojs_[i]->nTracklets();
  }

  for (unsigned int iME = 0; iME < nMatchEngines_; iME++) {
    matchengines_[iME].reset();
  }

  for (unsigned int istep = 0; istep < settings_.maxStep("MP") ; istep++) {
    bool projdone = false;
    bool medone = true;
    //Step 1
    //First step here checks if we have more input projections to put into
    //the input puffer for projections

    if (istep < settings_.maxStep("MP")) {
      if (iprojmem < inputprojs_.size()) {
        TrackletProjectionsMemory* projMem = inputprojs_[iprojmem];
        if (projMem->nTracklets() == 0) {
          iprojmem++;
        } else if (iproj < projMem->nTracklets()) {
          if (!inputProjBuffer_.almostfull()) {
            if (settings_.debugTracklet()) {
              edm::LogVerbatim("Tracklet") << getName() << " have projection in memory : " << projMem->getName();
            }

            Tracklet* proj = projMem->getTracklet(iproj);
            FPGAWord fpgaphi = barrel_ ? proj->layerProj(layerdisk_+1).fpgaphiproj() : proj->fpgaphiprojdisk(layerdisk_-N_LAYER+1);

            unsigned int iphi = (fpgaphi.value() >> (fpgaphi.nbits() - nvmbits_)) & (nvmbins_ - 1);

	    int nextrabits = 2;
	    int overlapbits = nvmbits_+nextrabits;
	      
	    unsigned int extrabits = fpgaphi.bits(fpgaphi.nbits() - overlapbits, nextrabits);

	    unsigned int ivmPlus = iphi;

	    int plusShift=0;
	    int negShift=0;
	    
	    
	    if (extrabits == ((1U << nextrabits) - 1) && iphi != ((1U << settings_.nbitsvmme(layerdisk_)) - 1)) {
	      plusShift = 1;
	      ivmPlus++;
	    }
	    unsigned int ivmMinus = iphi;
	    if (extrabits == 0 && iphi != 0){
	      negShift = -1;
	      ivmMinus--;
	    }

            int projrinv = -1;
            if (barrel_) {
              projrinv = 16 + (proj->fpgarinv().value() >> (proj->fpgarinv().nbits() - nrinv_));
            } else {
              //The next lines looks up the predicted bend based on:
              // 1 - r projections
              // 2 - phi derivative
              // 3 - the sign - i.e. if track is forward or backward
              int rindex = (proj->fpgarprojdisk(layerdisk_-N_LAYER+1).value() >> (proj->fpgarprojdisk(layerdisk_-N_LAYER+1).nbits() - nrbits_)) &
                           ((1 << nrbits_) - 1);

              int phiderindex = (proj->fpgaphiprojderdisk(layerdisk_-N_LAYER+1).value() >>
                                 (proj->fpgaphiprojderdisk(layerdisk_-N_LAYER+1).nbits() - nphiderbits_)) &
                                ((1 << nphiderbits_) - 1);
	      
              int signindex = proj->fpgarprojderdisk(layerdisk_-N_LAYER+1).value() < 0;

              int bendindex = (signindex << (nphiderbits_ + nrbits_)) + (rindex << (nphiderbits_)) + phiderindex;

              projrinv = globals_->projectionRouterBendTable()->bendLoookup(layerdisk_-N_LAYER, bendindex);

              proj->setBendIndex(projrinv, layerdisk_-N_LAYER+1);
            }
            assert(projrinv >= 0);

            unsigned int slot = barrel_ ? proj->layerProj(layerdisk_+1).fpgazbin1projvm().value() : proj->rbin1projvm(layerdisk_-N_LAYER+1);
            bool second = (barrel_ ? proj->layerProj(layerdisk_+1).fpgazbin2projvm().value() : proj->rbin2projvm(layerdisk_-N_LAYER+1));

            unsigned int projfinephi = (fpgaphi.value() >> (fpgaphi.nbits() - (nvmbits_ + 3))) & 7;
            int projfinerz = barrel_ ? proj->layerProj(layerdisk_+1).fpgafinezvm().value() : proj->finervm(layerdisk_-N_LAYER+1);

            bool isPSseed = proj->PSseed();

            int nbins = 8;
            if (layerdisk_ >= 6)
              nbins = 16;

            VMStubsMEMemory* stubmem = vmstubs_[0];
            bool usefirstPlus = stubmem->nStubsBin(ivmPlus * nbins + slot) != 0;
            bool usesecondPlus = (second && (stubmem->nStubsBin(ivmPlus * nbins + slot + 1) != 0));
            bool usefirstMinus = stubmem->nStubsBin(ivmMinus * nbins + slot) != 0;
            bool usesecondMinus = (second && (stubmem->nStubsBin(ivmMinus * nbins + slot + 1) != 0));

            if (usefirstPlus) {
              ProjectionTemp tmpProj(proj, slot, projrinv, projfinerz, projfinephi, ivmPlus, plusShift, usesecondPlus, isPSseed);
              inputProjBuffer_.store(tmpProj);
            } else if (usesecondPlus) {
              ProjectionTemp tmpProj(proj, slot + 1, projrinv, projfinerz - 8, projfinephi, ivmPlus, plusShift, false, isPSseed);
              inputProjBuffer_.store(tmpProj);
            }
	    if (ivmPlus!=ivmMinus) {
	      if (usefirstMinus) {
		ProjectionTemp tmpProj(proj, slot, projrinv, projfinerz, projfinephi, ivmMinus, negShift, usesecondMinus, isPSseed);
		inputProjBuffer_.store(tmpProj);
	      } else if (usesecondMinus) {
		ProjectionTemp tmpProj(proj, slot + 1, projrinv, projfinerz - 8, projfinephi, ivmMinus, negShift, false, isPSseed);
		inputProjBuffer_.store(tmpProj);
	      }
	    }
            iproj++;
            if (iproj == projMem->nTracklets()) {
              iproj = 0;
              iprojmem++;
            }
          }
        }
      } else {
        projdone = true;
      }
    }

    //Step 2
    //Check if we have ME that can process projection

    bool addedProjection = false;
    for (unsigned int iME = 0; iME < nMatchEngines_; iME++) {
      if (!matchengines_[iME].idle())
	countme++;
      matchengines_[iME].step();
      //if match engine empty and we have queued projections add to match engine
      if ((!addedProjection) && matchengines_[iME].idle() && (!inputProjBuffer_.empty())) {
	ProjectionTemp tmpProj = inputProjBuffer_.read();
	VMStubsMEMemory* stubmem = vmstubs_[0];
	
	if (settings_.debugTracklet()) {
	  edm::LogVerbatim("Tracklet") << getName() << " adding projection to match engine";
	}
	
	int nbins = 8;
	if (layerdisk_ >= 6)
	  nbins = 16;
	
	matchengines_[iME].init(stubmem,
				tmpProj.iphi() * nbins + tmpProj.slot(),
				tmpProj.projrinv(),
				tmpProj.projfinerz(),
				tmpProj.projfinephi(),
				tmpProj.shift(),
				tmpProj.usesecond(),
				tmpProj.isPSseed(),
				tmpProj.proj());
	addedProjection = true;
      }
    }

    //Step 3
    //Check if we have candidate match to process
    
    unsigned int iMEbest = nMatchEngines_;
    int bestTCID = -1;
    bool bestInPipeline = false;
    for (unsigned int iME = 0; iME < nMatchEngines_; iME++) {
      bool empty = matchengines_[iME].empty();
      medone = medone && (empty && matchengines_[iME].idle());
      if (empty && matchengines_[iME].idle())
	continue;
      int currentTCID = empty ? matchengines_[iME].currentProj()->TCID() : matchengines_[iME].peek().first->TCID();
      if ((iMEbest == nMatchEngines_) || (currentTCID < bestTCID)) {
	iMEbest = iME;
	bestTCID = currentTCID;
	bestInPipeline = empty;
      }
    }
    
    if (iMEbest != nMatchEngines_ && (!bestInPipeline)) {
      std::pair<Tracklet*, const Stub*> candmatch = matchengines_[iMEbest].read();
      
      const Stub* fpgastub = candmatch.second;
      Tracklet* tracklet = candmatch.first;
      
      if (oldTracklet != nullptr) {
	//allow equal here since we can have more than one cadidate match per tracklet projection
	assert(oldTracklet->TCID() <= tracklet->TCID());
      }
      oldTracklet = tracklet;
      
      bool match = matchCalculator(tracklet, fpgastub);
      
      if (settings_.debugTracklet() && match) {
	edm::LogVerbatim("Tracklet") << getName() << " have match";
      }
      
      countall++;
      if (match)
	countsel++;
    }
  
    if ((projdone && medone) || (istep == settings_.maxStep("MP") - 1)) {
      if (settings_.writeMonitorData("MP")) {
	globals_->ofstream("matchprocessor.txt") << getName() << " " << istep << " " << countall << " " << countsel << " "
						 << countme << " " << countinputproj << endl;
      }
      break;
    }
  }
  
  if (settings_.writeMonitorData("MC")) {
    globals_->ofstream("matchcalculator.txt") << getName() << " " << countall << " " << countsel << endl;
  }
}

bool MatchProcessor::matchCalculator(Tracklet* tracklet, const Stub* fpgastub) {
  const L1TStub* stub = fpgastub->l1tstub();

  if (layerdisk_ < N_LAYER) {
    const LayerProjection& layerProj = tracklet->layerProj(layerdisk_+1);
    int ir = fpgastub->r().value();
    int iphi = layerProj.fpgaphiproj().value();
    int icorr = (ir * layerProj.fpgaphiprojder().value()) >> icorrshift_;
    iphi += icorr;

    int iz = layerProj.fpgazproj().value();
    int izcor = (ir * layerProj.fpgazprojder().value() + (1 << (icorzshift_ - 1))) >> icorzshift_;
    iz += izcor;

    int ideltaz = fpgastub->z().value() - iz;
    int ideltaphi = (fpgastub->phi().value() - iphi) << phishift_;
    

    //Floating point calculations

    double phi = stub->phi();
    double r = stub->r();
    double z = stub->z();

    if (settings_.useapprox()) {
      double dphi = reco::reduceRange(phi - fpgastub->phiapprox(phimin_, phimax_));
      assert(std::abs(dphi) < 0.001);
      phi = fpgastub->phiapprox(phimin_, phimax_);
      z = fpgastub->zapprox();
      r = fpgastub->rapprox();
    }

    if (phi < 0)
      phi += 2 * M_PI;
    phi -= phimin_;

    double dr = r - settings_.rmean(layerdisk_);
    assert(std::abs(dr) < settings_.drmax());

    double dphi = reco::reduceRange(phi - (layerProj.phiproj() + dr * layerProj.phiprojder()));

    double dz = z - (layerProj.zproj() + dr * layerProj.zprojder());

    double dphiapprox =
        reco::reduceRange(phi - (layerProj.phiprojapprox() + dr * layerProj.phiprojderapprox()));

    double dzapprox = z - (layerProj.zprojapprox() + dr * layerProj.zprojderapprox());

    int seedindex = tracklet->getISeed();

    assert(phimatchcut_[seedindex] > 0);
    assert(zmatchcut_[seedindex] > 0);

    if (settings_.bookHistos()) {
      bool truthmatch = tracklet->stubtruthmatch(stub);

      HistBase* hists = globals_->histograms();
      hists->FillLayerResidual(layerdisk_+1,
                               seedindex,
                               dphiapprox * settings_.rmean(layerdisk_),
                               ideltaphi * settings_.kphi1() * settings_.rmean(layerdisk_),
                               (ideltaz << dzshift_ ) * settings_.kz(),
                               dz,
                               truthmatch);
    }

    if (settings_.writeMonitorData("Residuals")) {
      double pt = 0.01 * settings_.c() * settings_.bfield() / std::abs(tracklet->rinv());

      globals_->ofstream("layerresiduals.txt")
          << layerdisk_+1 << " " << seedindex << " " << pt << " "
          << ideltaphi * settings_.kphi1() * settings_.rmean(layerdisk_) << " "
          << dphiapprox * settings_.rmean(layerdisk_) << " "
          << phimatchcut_[seedindex] * settings_.kphi1() * settings_.rmean(layerdisk_) << "   "
          << (ideltaz << dzshift_) * settings_.kz() << " " << dz << " " << zmatchcut_[seedindex] * settings_.kz() << endl;
    }

    bool imatch = ((unsigned int)std::abs(ideltaphi) <= phimatchcut_[seedindex]) &&
      ((unsigned int)std::abs(ideltaz << dzshift_) <= zmatchcut_[seedindex]);

    if (settings_.debugTracklet()) {
      edm::LogVerbatim("Tracklet") << getName() << " imatch = " << imatch << " ideltaphi cut " << ideltaphi << " "
                                   << phimatchcut_[seedindex] << " ideltaz<<dzshift cut " << (ideltaz << dzshift_) << " "
                                   << zmatchcut_[seedindex];
    }

    //This would catch significant consistency problems in the configuration - helps to debug if there are problems.
    if (std::abs(dphi) > 0.5*settings_.dphisectorHG() || std::abs(dphiapprox) > 0.5*settings_.dphisectorHG()) {
      throw cms::Exception("LogicError") << "WARNING dphi and/or dphiapprox too large : "
					 << dphi << " " << dphiapprox << endl;
    }

    if (imatch) {
      tracklet->addMatch(layerdisk_+1,
                         ideltaphi,
                         ideltaz,
                         dphi,
                         dz,
                         dphiapprox,
                         dzapprox,
                         (phiregion_ << 7) + fpgastub->stubindex().value(),
                         stub->r(),
                         fpgastub);

      if (settings_.debugTracklet()) {
        edm::LogVerbatim("Tracklet") << "Accepted full match in layer " << getName() << " " << tracklet << " "
                                     << iSector_;
      }

      int iSeed = tracklet->getISeed();
      assert(fullmatches_[iSeed] != nullptr);
      fullmatches_[iSeed]->addMatch(tracklet, fpgastub);

      return true;
    } else {
      return false;
    }
  } else {  //disk matches

    //check that stubs and projections in same half of detector
    assert(stub->z() * tracklet->t() > 0.0);

    int sign = (tracklet->t() > 0.0) ? 1 : -1;
    int disk = sign * (layerdisk_-N_LAYER+1);
    assert(disk != 0);

    //Perform integer calculations here

    int iz = fpgastub->z().value();

    int iphi = tracklet->fpgaphiprojdisk(disk).value();
    int iphicorr = (iz * tracklet->fpgaphiprojderdisk(disk).value()) >> icorrshift_;
    iphi += iphicorr;

    int ir = tracklet->fpgarprojdisk(disk).value();
    int ircorr = (iz * tracklet->fpgarprojderdisk(disk).value()) >> icorzshift_;
    ir += ircorr;

    int ideltaphi = fpgastub->phi().value() - iphi;

    int irstub = fpgastub->r().value();
    int ialphafact = 0;
    if (!stub->isPSmodule()) {
      assert(irstub < (int)N_DSS_MOD * 2);
      if (layerdisk_-N_LAYER <= 1) {
        ialphafact = ialphafactinner_[irstub];
        irstub = settings_.rDSSinner(irstub) / settings_.kr();
      } else {
        ialphafact = ialphafactouter_[irstub];
        irstub = settings_.rDSSouter(irstub) / settings_.kr();
      }
    }

    int ideltar = (irstub * settings_.kr()) / settings_.krprojshiftdisk() - ir;

    if (!stub->isPSmodule()) {
      int ialphanew = fpgastub->alphanew().value();
      int iphialphacor = ((ideltar * ialphanew * ialphafact) >> settings_.alphashift());
      ideltaphi += iphialphacor;
    }

    //Perform floating point calculations here

    double phi = stub->phi();
    double z = stub->z();
    double r = stub->r();

    if (settings_.useapprox()) {
      double dphi = reco::reduceRange(phi - fpgastub->phiapprox(phimin_, phimax_));
      assert(std::abs(dphi) < 0.001);
      phi = fpgastub->phiapprox(phimin_, phimax_);
      z = fpgastub->zapprox();
      r = fpgastub->rapprox();
    }

    if (phi < 0)
      phi += 2 * M_PI;
    phi -= phimin_;

    double dz = z - sign * settings_.zmean(layerdisk_ - N_LAYER);

    if (std::abs(dz) > settings_.dzmax()) {
      edm::LogProblem("Tracklet") << __FILE__ << ":" << __LINE__ << " " << name_ << "_" << iSector_ << " "
                                  << tracklet->getISeed();
      edm::LogProblem("Tracklet") << "stub " << stub->z() << " disk " << disk << " " << dz;
      assert(std::abs(dz) < settings_.dzmax());
    }

    double phiproj = tracklet->phiprojdisk(disk) + dz * tracklet->phiprojderdisk(disk);
    double rproj = tracklet->rprojdisk(disk) + dz * tracklet->rprojderdisk(disk);
    double deltar = r - rproj;

    double dr = stub->r() - rproj;
    double drapprox = stub->r() - (tracklet->rprojapproxdisk(disk) + dz * tracklet->rprojderapproxdisk(disk));

    double dphi = reco::reduceRange(phi - phiproj);
    double dphiapprox =
        reco::reduceRange(phi - (tracklet->phiprojapproxdisk(disk) + dz * tracklet->phiprojderapproxdisk(disk)));

    double drphi = dphi * stub->r();
    double drphiapprox = dphiapprox * stub->r();

    if (!stub->isPSmodule()) {
      double alphanorm = stub->alphanorm();
      dphi += dr * alphanorm * settings_.half2SmoduleWidth() / stub->r2();
      ;
      dphiapprox += drapprox * alphanorm * settings_.half2SmoduleWidth() / stub->r2();

      drphi += dr * alphanorm * settings_.half2SmoduleWidth() / stub->r();
      drphiapprox += dr * alphanorm * settings_.half2SmoduleWidth() / stub->r();
    }

    int seedindex = tracklet->getISeed();

    int idrphicut = rphicutPS_[seedindex];
    int idrcut = rcutPS_[seedindex];
    if (!stub->isPSmodule()) {
      idrphicut = rphicut2S_[seedindex];
      idrcut = rcut2S_[seedindex];
    }

    double drphicut = idrphicut * settings_.kphi() * settings_.kr();
    double drcut = idrcut * settings_.krprojshiftdisk();

    if (settings_.writeMonitorData("Residuals")) {
      double pt = 0.01 * settings_.c() * settings_.bfield() / std::abs(tracklet->rinv());

      globals_->ofstream("diskresiduals.txt")
          << layerdisk_-N_LAYER+1 << " " << stub->isPSmodule() << " " << tracklet->layer() << " " << abs(tracklet->disk()) << " " << pt
          << " " << ideltaphi * settings_.kphi() * stub->r() << " " << drphiapprox << " " << drphicut << " "
          << ideltar * settings_.krprojshiftdisk() << " " << deltar << " " << drcut << " " << endl;
    }

    bool match = (std::abs(drphi) < drphicut) && (std::abs(deltar) < drcut);
    bool imatch = (std::abs(ideltaphi * irstub) < idrphicut) && (std::abs(ideltar) < idrcut);

    if (settings_.debugTracklet()) {
      edm::LogVerbatim("Tracklet") << "imatch match disk: " << imatch << " " << match << " " << std::abs(ideltaphi)
                                   << " " << drphicut / (settings_.kphi() * stub->r()) << " " << std::abs(ideltar)
                                   << " " << drcut / settings_.krprojshiftdisk() << " r = " << stub->r();
    }

    if (imatch) {
      if (settings_.debugTracklet()) {
        edm::LogVerbatim("Tracklet") << "MatchCalculator found match in disk " << getName();
      }

      if (std::abs(dphi) >= 0.25) {
        edm::LogPrint("Tracklet") << "dphi " << dphi << " ISeed " << tracklet->getISeed();
      }
      assert(std::abs(dphi) < 0.25);
      assert(std::abs(dphiapprox) < 0.25);

      tracklet->addMatchDisk(disk,
                             ideltaphi,
                             ideltar,
                             drphi / stub->r(),
                             dr,
                             drphiapprox / stub->r(),
                             drapprox,
                             stub->alpha(settings_.stripPitch(stub->isPSmodule())),
                             (phiregion_ << 7) + fpgastub->stubindex().value(),
                             stub->z(),
                             fpgastub);
      if (settings_.debugTracklet()) {
        edm::LogVerbatim("Tracklet") << "Accepted full match in disk " << getName() << " " << tracklet << " "
                                     << iSector_;
      }

      int iSeed = tracklet->getISeed();
      assert(fullmatches_[iSeed] != nullptr);
      fullmatches_[iSeed]->addMatch(tracklet, fpgastub);

      return true;
    } else {
      return false;
    }
  }
}

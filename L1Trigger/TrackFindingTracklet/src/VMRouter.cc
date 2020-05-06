#include "L1Trigger/TrackFindingTracklet/interface/VMRouter.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"
#include "L1Trigger/TrackFindingTracklet/interface/TETableOuter.h"
#include "L1Trigger/TrackFindingTracklet/interface/TETableInner.h"
#include "L1Trigger/TrackFindingTracklet/interface/TETableOuterDisk.h"
#include "L1Trigger/TrackFindingTracklet/interface/TETableInnerDisk.h"
#include "L1Trigger/TrackFindingTracklet/interface/TETableInnerOverlap.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubTE.h"
#include "L1Trigger/TrackFindingTracklet/interface/InputLinkMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubsMEMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubsTEMemory.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;
using namespace Trklet;

VMRouter::VMRouter(string name, const Settings* settings, Globals* global, unsigned int iSector)
    : ProcessBase(name, settings, global, iSector) {
  layerdisk_ = initLayerDisk(4);
  initFineBinTable();

  vmstubsMEPHI_.resize(settings_->nvmme(layerdisk_), 0);

  overlapbits_ = 7;
  nextrabits_ = overlapbits_ - (settings_->nbitsallstubs(layerdisk_) + settings_->nbitsvmme(layerdisk_));
}

void VMRouter::addOutput(MemoryBase* memory, string output) {
  if (settings_->writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding output to " << memory->getName() << " to output "
                                 << output;
  }

  if (output.substr(0, 10) == "allstubout") {
    AllStubsMemory* tmp = dynamic_cast<AllStubsMemory*>(memory);
    assert(tmp != 0);
    allstubs_.push_back(tmp);
    return;
  }

  if (output.substr(0, 12) == "vmstuboutPHI") {
    char seedtype = memory->getName().substr(11, 1)[0];
    unsigned int pos = 12;
    int vmbin = memory->getName().substr(pos, 1)[0] - '0';
    pos++;
    if (pos < memory->getName().size()) {
      if (memory->getName().substr(pos, 1)[0] != 'n') {
        vmbin = vmbin * 10 + memory->getName().substr(pos, 1)[0] - '0';
        pos++;
      }
    }

    int iseed = -1;
    unsigned int inner = 1;
    if (memory->getName().substr(3, 2) == "TE") {
      VMStubsTEMemory* tmp = dynamic_cast<VMStubsTEMemory*>(memory);
      assert(tmp != 0);
      if (seedtype < 'I') {
        if (layerdisk_ == 0 || layerdisk_ == 1)
          iseed = 0;
        if (layerdisk_ == 2 || layerdisk_ == 3)
          iseed = 2;
        if (layerdisk_ == 4 || layerdisk_ == 5)
          iseed = 3;
        if (layerdisk_ == 6 || layerdisk_ == 7)
          iseed = 4;
        if (layerdisk_ == 8 || layerdisk_ == 9)
          iseed = 5;
        if (layerdisk_ == 0 || layerdisk_ == 2 || layerdisk_ == 4 || layerdisk_ == 6 || layerdisk_ == 8)
          inner = 0;
      } else if (seedtype < 'M') {
        if (layerdisk_ == 1 || layerdisk_ == 2)
          iseed = 1;
        if (layerdisk_ == 1)
          inner = 0;
      } else if (seedtype <= 'Z') {
        if (layerdisk_ == 0 || layerdisk_ == 6)
          iseed = 6;
        if (layerdisk_ == 1 || layerdisk_ == 6)
          iseed = 7;
        if (layerdisk_ == 0 || layerdisk_ == 1)
          inner = 0;
      } else if (seedtype < 'o' && seedtype >= 'a') {
        if (layerdisk_ == 1 || layerdisk_ == 2)
          iseed = 10;
        if (layerdisk_ == 1)
          inner = 0;
      } else if (seedtype > 'o' && seedtype <= 'z') {
        if (layerdisk_ == 1)
          iseed = 11;
        if (layerdisk_ == 6)
          iseed = 10;
        inner = 2;
      } else {
        assert(0);
      }
      assert(iseed != -1);
      int seedindex = -1;
      for (unsigned int k = 0; k < vmstubsTEPHI_.size(); k++) {
        if (vmstubsTEPHI_[k].first.first == (unsigned int)iseed) {
          seedindex = k;
        }
      }
      if (seedindex == -1) {
        seedindex = vmstubsTEPHI_.size();
        vector<VMStubsTEMemory*> avectmp;
        vector<vector<VMStubsTEMemory*> > vectmp(settings_->nvmte(inner, iseed), avectmp);
        pair<unsigned int, unsigned int> tmppair(iseed, inner);
        std::pair<std::pair<unsigned int, unsigned int>, vector<vector<VMStubsTEMemory*> > > atmp(tmppair, vectmp);
        vmstubsTEPHI_.push_back(atmp);
      }
      vmstubsTEPHI_[seedindex].second[(vmbin - 1) & (settings_->nvmte(inner, iseed) - 1)].push_back(tmp);

    } else if (memory->getName().substr(3, 2) == "ME") {
      VMStubsMEMemory* tmp = dynamic_cast<VMStubsMEMemory*>(memory);
      assert(tmp != 0);
      vmstubsMEPHI_[(vmbin - 1) & (settings_->nvmme(layerdisk_) - 1)] = tmp;
    } else {
      assert(0);
    }

    return;
  }

  edm::LogPrint("Tracklet") << "Could not find : " << output;
  assert(0);
}

void VMRouter::addInput(MemoryBase* memory, string input) {
  if (settings_->writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding input from " << memory->getName() << " to input "
                                 << input;
  }
  if (input == "stubin") {
    InputLinkMemory* tmp1 = dynamic_cast<InputLinkMemory*>(memory);
    assert(tmp1 != 0);
    if (tmp1 != 0) {
      stubinputs_.push_back(tmp1);
    }
    return;
  }
  edm::LogPrint("Tracklet") << "Could not find input : " << input;
  assert(0);
}

void VMRouter::execute() {
  unsigned int allStubCounter = 0;

  //Loop over the input stubs
  for (unsigned int j = 0; j < stubinputs_.size(); j++) {
    for (unsigned int i = 0; i < stubinputs_[j]->nStubs(); i++) {
      if (allStubCounter > settings_->maxStep("VMR"))
        continue;
      if (allStubCounter > 127)
        continue;
      std::pair<Stub*, L1TStub*> stub = stubinputs_[j]->getStub(i);

      bool negdisk =
          (stub.first->disk().value() <
           0);  //Note - this information is not part of the stub - but rather from which input memory we are reading

      //use &127 to make sure we fit into the number of bits -
      //though we should have protected against overflows above
      FPGAWord allStubIndex(allStubCounter & 127, 7, true, __LINE__, __FILE__);

      stub.first->setAllStubIndex(
          allStubCounter);  //TODO - should not be needed - but need to migrate some other pieces of code before removing
      stub.second->setAllStubIndex(
          allStubCounter);  //TODO - should not be needed - but need to migrate some other pieces of code before removing

      allStubCounter++;

      //Fill allstubs memories - in HLS this is the same write to multiple memories
      for (unsigned int l = 0; l < allstubs_.size(); l++) {
        allstubs_[l]->addStub(stub);
      }

      //Fill all the ME VM memories

      FPGAWord iphi = stub.first->phicorr();
      unsigned int ivm =
          iphi.bits(iphi.nbits() - (settings_->nbitsallstubs(layerdisk_) + settings_->nbitsvmme(layerdisk_)),
                    settings_->nbitsvmme(layerdisk_));
      unsigned int extrabits = iphi.bits(iphi.nbits() - overlapbits_, nextrabits_);

      unsigned int ivmPlus = ivm;

      if (extrabits == ((1U << nextrabits_) - 1) && ivm != ((1U << settings_->nbitsvmme(layerdisk_)) - 1))
        ivmPlus++;
      unsigned int ivmMinus = ivm;
      if (extrabits == 0 && ivm != 0)
        ivmMinus--;

      //Calculate the z and r position for the vmstub

      //Take the top nbitszfinebintable_ bits of the z coordinate
      int indexz = (((1 << (stub.first->z().nbits() - 1)) + stub.first->z().value()) >>
                    (stub.first->z().nbits() - nbitszfinebintable_));
      int indexr = -1;
      if (layerdisk_ > 5) {
        if (negdisk) {
          indexz = (1 << nbitszfinebintable_) - indexz;
        }
        indexr = stub.first->r().value();
        if (stub.first->isPSmodule()) {
          indexr = stub.first->r().value() >> (stub.first->r().nbits() - nbitsrfinebintable_);
        }
      } else {
        //Take the top nbitsfinebintable_ bits of the z coordinate. The & is to handle the negative z values.
        indexr = (((1 << (stub.first->r().nbits() - 1)) + stub.first->r().value()) >>
                  (stub.first->r().nbits() - nbitsrfinebintable_));
      }

      assert(indexz >= 0);
      assert(indexr >= 0);
      assert(indexz < (1 << nbitszfinebintable_));
      assert(indexr < (1 << nbitsrfinebintable_));
      unsigned int index = (indexz << nbitsrfinebintable_) + indexr;
      assert(index < finebintable_.size());

      int rzfine = finebintable_[index];

      assert(rzfine >= 0);

      int vmbin = rzfine >> 3;
      if (negdisk)
        vmbin += 8;
      rzfine = rzfine & 7;

      if (layerdisk_ > 5) {
        stub.first->setfiner(rzfine);
      } else {
        stub.first->setfinez(rzfine);
      }

      VMStubME vmstub(stub,
                      stub.first->iphivmFineBins(
                          iphi.nbits() - (settings_->nbitsallstubs(layerdisk_) + settings_->nbitsvmme(layerdisk_)),
                          settings_->nbitsvmme(layerdisk_)),
                      FPGAWord(rzfine, 3, true, __LINE__, __FILE__),
                      stub.first->bend(),
                      allStubIndex);

      assert(vmstubsMEPHI_[ivmPlus] != 0);
      vmstubsMEPHI_[ivmPlus]->addStub(vmstub, vmbin);

      if (ivmMinus != ivmPlus) {
        assert(vmstubsMEPHI_[ivmMinus] != 0);
        vmstubsMEPHI_[ivmMinus]->addStub(vmstub, vmbin);
      }

      //Fill the TE VM memories

      for (unsigned int i = 0; i < vmstubsTEPHI_.size(); i++) {
        unsigned int iseed = vmstubsTEPHI_[i].first.first;
        unsigned int inner = vmstubsTEPHI_[i].first.second;

        if ((iseed == 4 || iseed == 5 || iseed == 6 || iseed == 7) && (!stub.first->isPSmodule()))
          continue;

        FPGAWord binlookup = lookup(iseed, inner, stub.first->z(), stub.first->r(), negdisk, stub.first->isPSmodule());

        if (binlookup.value() < 0)
          continue;

        unsigned int ivmte =
            iphi.bits(iphi.nbits() - (settings_->nbitsallstubs(layerdisk_) + settings_->nbitsvmte(inner, iseed)),
                      settings_->nbitsvmte(inner, iseed));

        int bin = -1;
        if (inner != 0) {
          bin = binlookup.value() / 8;
          unsigned int tmp = binlookup.value() & 7;  //three bits in outer layers - this could be coded cleaner...
          binlookup.set(tmp, 3, true, __LINE__, __FILE__);
        }

        FPGAWord finephi =
            stub.first->iphivmFineBins(settings_->nphireg(inner, iseed), settings_->nfinephi(inner, iseed));

        VMStubTE tmpstub(stub, finephi, stub.first->bend(), binlookup, allStubIndex);

        unsigned int nmem = vmstubsTEPHI_[i].second[ivmte].size();

        assert(nmem > 0);

        for (unsigned int l = 0; l < nmem; l++) {
          if (settings_->debugTracklet()) {
            edm::LogVerbatim("Tracklet") << getName() << " try adding stub to "
                                         << vmstubsTEPHI_[i].second[ivmte][l]->getName() << " inner=" << inner
                                         << " bin=" << bin;
          }
          if (inner == 0) {
            vmstubsTEPHI_[i].second[ivmte][l]->addVMStub(tmpstub);
          } else {
            vmstubsTEPHI_[i].second[ivmte][l]->addVMStub(tmpstub, bin);
          }
        }
      }
    }
  }
}

void VMRouter::initFineBinTable() {
  if (layerdisk_ < 6) {
    nbitszfinebintable_ = 7;
    nbitsrfinebintable_ = 4;
    //unsigned int nbins = 1 << (nbitszfinebintable_+nbitsrfinebintable_);
    //finebintable_.resize(nbins, -1);

    double dr = 2 * settings_->drmax() / (1 << nbitsrfinebintable_);
    double dz = 2 * settings_->zlength() / (1 << nbitszfinebintable_);

    double rmin = settings_->rmean(layerdisk_) - settings_->drmax();
    double zmin = -settings_->zlength();
    int NBINS = settings_->NLONGVMBINS() * settings_->NLONGVMBINS();

    for (int izbin = 0; izbin < (1 << nbitszfinebintable_); izbin++) {
      for (int irbin = 0; irbin < (1 << nbitsrfinebintable_); irbin++) {
        double z = zmin + (izbin + 0.5) * dz;
        double r = rmin + (irbin + 0.5) * dr;

        double zproj = z * settings_->rmean(layerdisk_) / r;

        int zbin = NBINS * (zproj + settings_->zlength()) / (2 * settings_->zlength());

        if (zbin < 0)
          zbin = 0;
        if (zbin >= NBINS)
          zbin = NBINS - 1;
        finebintable_.push_back(zbin);
      }
    }
  } else {
    nbitszfinebintable_ = 3;
    nbitsrfinebintable_ = 7;
    //unsigned int nbins = 1 << (nbitszfinebintable_+nbitsrfinebintable_);

    double rmin = 0.0;
    double dr = settings_->rmaxdisk() / (1 << nbitsrfinebintable_);
    ;

    double zmin = settings_->zmean(layerdisk_ - 6) - settings_->dzmax();
    double dz = 2 * settings_->dzmax() / (1 << nbitszfinebintable_);
    ;

    int NBINS = settings_->NLONGVMBINS() * settings_->NLONGVMBINS();

    for (int izbin = 0; izbin < (1 << nbitszfinebintable_); izbin++) {
      for (int irbin = 0; irbin < (1 << nbitsrfinebintable_); irbin++) {
        double r = rmin + (irbin + 0.5) * dr;
        double z = zmin + (izbin + 0.5) * dz;

        if (irbin < 10)  //special case for the tabulated radii
          r = (layerdisk_ <= 7) ? settings_->rDSSinner(irbin) : settings_->rDSSouter(irbin);

        double rproj = r * settings_->zmean(layerdisk_ - 6) / z;

        int rbin = NBINS * (rproj - settings_->rmindiskvm()) / (settings_->rmaxdisk() - settings_->rmindiskvm());

        if (rbin < 0)
          rbin = 0;
        if (rbin >= NBINS)
          rbin = NBINS - 1;

        finebintable_.push_back(rbin);
      }
    }
  }

  if (iSector_==0&&settings_->writeTable()) {
    
    ofstream outfinebin;
    outfinebin.open(getName()+"_finebin.tab");
    outfinebin << "{"<<endl;
    for(unsigned int i=0;i<finebintable_.size();i++) {
      if (i!=0) outfinebin<<","<<endl;
      outfinebin << finebintable_[i];
    }
    outfinebin <<endl<<"};"<<endl;
    outfinebin.close();
  }
} 


FPGAWord VMRouter::lookup(
    unsigned int iseed, unsigned int inner, FPGAWord z, FPGAWord r, bool negdisk, bool isPSmodule) {
  if (globals_->teTable(0, 0) == 0) {
    globals_->teTable(0, 0) =
        new TETableInner(settings_, 1, 2, -1, settings_->zbitstab(0, 0), settings_->rbitstab(0, 0));
    globals_->teTable(0, 1) =
        new TETableInner(settings_, 2, 3, -1, settings_->zbitstab(0, 1), settings_->rbitstab(0, 1));
    globals_->teTable(0, 2) =
        new TETableInner(settings_, 3, 4, 2, settings_->zbitstab(0, 2), settings_->rbitstab(0, 2));
    globals_->teTable(0, 3) =
        new TETableInner(settings_, 5, 6, 4, settings_->zbitstab(0, 3), settings_->rbitstab(0, 3));
    globals_->teTable(0, 4) =
        new TETableInnerDisk(settings_, 1, 2, -1, settings_->zbitstab(0, 4), settings_->rbitstab(0, 4));
    globals_->teTable(0, 5) =
        new TETableInnerDisk(settings_, 3, 4, -1, settings_->zbitstab(0, 5), settings_->rbitstab(0, 5));
    globals_->teTable(0, 6) =
        new TETableInnerOverlap(settings_, 1, 1, settings_->zbitstab(0, 6), settings_->rbitstab(0, 6));
    globals_->teTable(0, 7) =
        new TETableInnerOverlap(settings_, 2, 1, settings_->zbitstab(0, 7), settings_->rbitstab(0, 7));
    globals_->teTable(0, 10) =
        new TETableInner(settings_, 2, 3, 1, settings_->zbitstab(0, 10), settings_->rbitstab(0, 10), true);

    globals_->teTable(1, 0) = new TETableOuter(settings_, 2, settings_->zbitstab(1, 0), settings_->rbitstab(1, 0));
    globals_->teTable(1, 1) = new TETableOuter(settings_, 3, settings_->zbitstab(1, 1), settings_->rbitstab(1, 1));
    globals_->teTable(1, 2) = new TETableOuter(settings_, 4, settings_->zbitstab(1, 2), settings_->rbitstab(1, 2));
    globals_->teTable(1, 3) = new TETableOuter(settings_, 6, settings_->zbitstab(1, 3), settings_->rbitstab(1, 3));
    globals_->teTable(1, 4) = new TETableOuterDisk(settings_, 2, settings_->zbitstab(1, 4), settings_->rbitstab(1, 4));
    globals_->teTable(1, 5) = new TETableOuterDisk(settings_, 4, settings_->zbitstab(1, 5), settings_->rbitstab(1, 5));
    globals_->teTable(1, 6) = new TETableOuterDisk(settings_, 1, settings_->zbitstab(1, 6), settings_->rbitstab(1, 6));
    globals_->teTable(1, 7) = new TETableOuterDisk(settings_, 1, settings_->zbitstab(1, 7), settings_->rbitstab(1, 7));
    globals_->teTable(1, 10) = new TETableOuter(settings_, 3, settings_->zbitstab(1, 10), settings_->rbitstab(1, 10));

    globals_->teTable(2, 10) =
        new TETableOuterDisk(settings_, 1, settings_->zbitstab(2, 10), settings_->rbitstab(2, 10));
    globals_->teTable(2, 11) = new TETableOuter(settings_, 2, settings_->zbitstab(2, 11), settings_->rbitstab(2, 11));
  }

  if (iseed == 10 && inner == 2) {
    //return if radius to small (values <100 corresponds to the 2S modules)
    if (r.value() > 100 && r.value() < settings_->rmindiskl3overlapvm() / settings_->kr())
      return FPGAWord(-1, 2, false, __LINE__, __FILE__);
    int bin = 0;
    if (!isPSmodule) {
      bin = r.value();  // 0 to 9 for the ring index
      bin = bin >> 2;   // 0 to 2
      bin += 1;
    }

    return FPGAWord(bin * 8, settings_->lutwidthtabextended(inner, iseed), true, __LINE__, __FILE__);
  }

  assert(globals_->teTable(inner, iseed) != 0);

  unsigned int zbits = settings_->zbitstab(inner, iseed);
  assert(zbits != 0);
  unsigned int rbits = settings_->rbitstab(inner, iseed);
  assert(rbits != 0);
  unsigned int lutwidth = settings_->lutwidthtab(inner, iseed);
  if (settings_->extended()) {
    lutwidth = settings_->lutwidthtabextended(inner, iseed);
  }
  assert(lutwidth != 0);

  int zbin = (z.value() + (1 << (z.nbits() - 1))) >> (z.nbits() - zbits);
  if (negdisk)
    zbin = (1 << zbits) - 1 - zbin;
  int rbin = (r.value() + (1 << (r.nbits() - 1))) >> (r.nbits() - rbits);
  if (layerdisk_ >= 6) {
    rbin = r.value() >> (r.nbits() - rbits);
  }

  int lutvalue = globals_->teTable(inner, iseed)->lookup(zbin, rbin);

  if (lutvalue < 0) {
    return FPGAWord(lutvalue, 2, false, __LINE__, __FILE__);
  } else {
    return FPGAWord(lutvalue, lutwidth, true, __LINE__, __FILE__);
  }
}

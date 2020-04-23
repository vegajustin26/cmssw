#ifndef VMSTUBME_H
#define VMSTUBME_H

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>

#include "FPGAWord.h"
#include "Stub.h"
#include "L1TStub.h"

using namespace std;

class VMStubME {
public:
  VMStubME() {}

  VMStubME(std::pair<Stub*, L1TStub*> stub, FPGAWord finephi, FPGAWord finerz, FPGAWord bend, FPGAWord allstubindex) {
    stub_ = stub;
    finephi_ = finephi;
    finerz_ = finerz;
    bend_ = bend;
    allStubIndex_ = allstubindex;
  }

  ~VMStubME() {}

  FPGAWord finephi() const { return finephi_; }

  FPGAWord finerz() const { return finerz_; }

  FPGAWord bend() const { return bend_; }

  std::pair<Stub*, L1TStub*> stub() const { return stub_; }

  bool isPSmodule() const { return stub_.first->isPSmodule(); }

  FPGAWord stubindex() const { return allStubIndex_; }

  //return binary string for memory printout
  std::string str() const {
    string stub = allStubIndex_.str();
    stub += "|";
    stub += bend_.str();
    stub += "|";
    stub += finephi_.str();
    stub += "|";
    stub += finerz_.str();

    return stub;
  }

private:
  FPGAWord allStubIndex_;
  FPGAWord finephi_;
  FPGAWord finerz_;
  FPGAWord bend_;
  std::pair<Stub*, L1TStub*> stub_;
};

#endif

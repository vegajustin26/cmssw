//Base class for processing modules
#ifndef PROCESSBASE_H
#define PROCESSBASE_H

using namespace std;

#include "Stub.h"

class ProcessBase {
public:
  ProcessBase(string name, unsigned int iSector) {
    name_ = name;
    iSector_ = iSector;
    extended_ = hourglassExtended;
    nHelixPar_ = nHelixPar;
  }

  virtual ~ProcessBase() {}

  // Add wire from pin "output" or "input" this proc module to memory instance "memory".

  virtual void addOutput(MemoryBase* memory, string output) = 0;

  virtual void addInput(MemoryBase* memory, string input) = 0;

  string getName() const { return name_; }

  unsigned int nbits(unsigned int power) {
    if (power == 2)
      return 1;
    if (power == 4)
      return 2;
    if (power == 8)
      return 3;
    if (power == 16)
      return 4;
    if (power == 32)
      return 5;

    cout << "nbits: power = " << power << endl;
    assert(0);

    return -1;
  }

  //method sets the layer and disk based on the name. pos is the position in the
  //memory name where the layer or disk is specified
  void initLayerDisk(unsigned int pos, int& layer, int& disk) {
    string subname = name_.substr(pos, 2);
    layer = 0;
    disk = 0;

    if (subname == "L1")
      layer = 1;
    if (subname == "L2")
      layer = 2;
    if (subname == "L3")
      layer = 3;
    if (subname == "L4")
      layer = 4;
    if (subname == "L5")
      layer = 5;
    if (subname == "L6")
      layer = 6;
    if (subname == "D1")
      disk = 1;
    if (subname == "D2")
      disk = 2;
    if (subname == "D3")
      disk = 3;
    if (subname == "D4")
      disk = 4;
    if (subname == "D5")
      disk = 5;
    if (layer == 0 && disk == 0) {
      cout << "Memoryname = " << name_ << " subname = " << subname << " layer " << layer << " disk " << disk << endl;
    }
    assert((layer != 0) || (disk != 0));
  }

  void initLayerDisk(unsigned int pos, int& layer, int& disk, int& layerdisk) {
    initLayerDisk(pos, layer, disk);

    layerdisk = layer - 1;
    if (disk > 0)
      layerdisk = 5 + disk;
  }

  unsigned int initLayerDisk(unsigned int pos) {
    int layer, disk;
    initLayerDisk(pos, layer, disk);

    if (disk > 0)
      return 5 + disk;
    return layer - 1;
  }

  unsigned int getISeed(std::string name) {
    //assumes here that namme is on the form XX_L1L2_XXX where L1L2 gives iSeed=0

    std::size_t pos = name.find("_");
    std::string name1 = name.substr(pos + 1);

    pos = name1.find("_");
    std::string name2 = name1.substr(0, pos);

    if (name2 == "L1L2")
      return 0;
    if (name2 == "L2L3")
      return 1;
    if (name2 == "L3L4")
      return 2;
    if (name2 == "L5L6")
      return 3;
    if (name2 == "D1D2")
      return 4;
    if (name2 == "D3D4")
      return 5;
    if (name2 == "L1D1")
      return 6;
    if (name2 == "L2D1")
      return 7;

    if (name2 == "L1L2XX")
      return 0;
    if (name2 == "L2L3XX")
      return 1;
    if (name2 == "L3L4XX")
      return 2;
    if (name2 == "L5L6XX")
      return 3;
    if (name2 == "D1D2XX")
      return 4;
    if (name2 == "D3D4XX")
      return 5;
    if (name2 == "L1D1XX")
      return 6;
    if (name2 == "L2D1XX")
      return 7;
    if (name2 == "L3L4L2")
      return 8;
    if (name2 == "L5L6L4")
      return 9;
    if (name2 == "L2L3D1")
      return 10;
    if (name2 == "D1D2L2")
      return 11;

    cout << getName() << " name name1 name2 " << name << " - " << name1 << " - " << name2 << endl;
    assert(0);
    return 0;
  }

protected:
  string name_;
  unsigned int iSector_;
  bool extended_;
  unsigned int nHelixPar_;
};

#endif

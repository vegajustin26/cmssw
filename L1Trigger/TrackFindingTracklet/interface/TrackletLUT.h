#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletLUT_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletLUT_h

#include <string>
#include <vector>

namespace trklet {

  class Settings;

  class TrackletLUT {
    public:

    TrackletLUT(const Settings& settings);

    TrackletLUT& operator=(const TrackletLUT& other) {

      name_ = other.name_;
      table_ = other.table_;
      nbits_ = other.nbits_;
      positive_ = other.positive_;
      
      return *this;
    }
  
    ~TrackletLUT() = default;

    void initteptlut(bool fillInner, bool fillTEMem, unsigned int iSeed, unsigned int layerdisk1, unsigned int layerdisk2,
		     unsigned int innerphibits, unsigned int outerphibits,
		     double innerphimin, double innerphimax, double outerphimin, double outerphimax);

    
    void initProjectionBend(double k_phider, unsigned int idisk, unsigned int nrbits, unsigned int nphiderbits);
    
    void initBendMatch(unsigned int layerdisk);
    
    enum VMRTableType { me, disk, inner, inneroverlap, innerthird };

    void initVMRTable(unsigned int layerdisk, VMRTableType type);

    void initPhiCorrTable(unsigned int layerdisk, unsigned int rbits);

    void writeTable() const;
    
    int lookup(unsigned int index) const;

  private:

    int getphiCorrValue(unsigned int layerdisk, unsigned int ibend, unsigned int irbin,
			double rmean, double dr, double drmax) const;

    int getVMRLookup(unsigned int layerdisk, double z, double r, double dz, double dr, int iseed = -1) const;
    
    const Settings& settings_;
    
    std::string name_;

    std::vector<int> table_;

    unsigned int nbits_;
    
    bool positive_;
    
  };
};  // namespace trklet
#endif

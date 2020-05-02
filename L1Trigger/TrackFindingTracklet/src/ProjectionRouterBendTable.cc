#include "L1Trigger/TrackFindingTracklet/interface/ProjectionRouterBendTable.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"

using namespace Trklet;

void ProjectionRouterBendTable::init(const Settings* settings, Globals* globals, unsigned int nrbits, unsigned int nphiderbits) {

    for (unsigned int idisk=0;idisk<5;idisk++) {

      unsigned int nsignbins=2;
      unsigned int nrbins=1<<(nrbits);
      unsigned int nphiderbins=1<<(nphiderbits);
      
      for(unsigned int isignbin=0;isignbin<nsignbins;isignbin++) {
	for(unsigned int irbin=0;irbin<nrbins;irbin++) {
	  int ir=irbin;
	  if (ir>(1<<(nrbits-1))) ir-=(1<<nrbits);
	  ir=ir<<(settings->nrbitsstub(6)-nrbits);
	  for(unsigned int iphiderbin=0;iphiderbin<nphiderbins;iphiderbin++) {
	    int iphider=iphiderbin;
	    if (iphider>(1<<(nphiderbits-1))) iphider-=(1<<nphiderbits);
	    iphider=iphider<<(settings->nbitsphiprojderL123()-nphiderbits);
	      
	    double rproj=ir*settings->krprojshiftdisk();
	    double phider=iphider*globals->ITC_L1L2()->der_phiD_final.get_K();
	    double t=settings->zmean(idisk)/rproj;
	    
	    if (isignbin) t=-t;
	    
	    double rinv=-phider*(2.0*t);
	    
	    double bendproj=0.5*bend(rproj,rinv);
	    
	    
	    int ibendproj=2.0*bendproj+15.5;
	    if (ibendproj<0) ibendproj=0;
	    if (ibendproj>31) ibendproj=31;
	    
	    bendtable_[idisk].push_back(ibendproj);
	    
	  }
	}
      }
    }
}
   
int ProjectionRouterBendTable::bendLoookup(int diskindex,int bendindex) {
  return bendtable_[diskindex][bendindex];
}

double ProjectionRouterBendTable::bend(double r, double rinv) {
  
  double dr=0.18;
  double delta=r*dr*0.5*rinv;
  
  double bend=-delta/0.009;
  if (r<55.0) bend=-delta/0.01;
  
  return bend;
}

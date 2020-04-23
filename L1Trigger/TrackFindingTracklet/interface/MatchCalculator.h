//This class implementes the match calculator
#ifndef MATCHCALCULATOR_H
#define MATCHCALCULATOR_H

#include "ProcessBase.h"
#include "GlobalHistTruth.h"
#include "Util.h"

using namespace std;

class MatchCalculator:public ProcessBase{

public:

 MatchCalculator(string name, const Settings* settings, unsigned int iSector):
  ProcessBase(name,settings,iSector){
    
    double dphi=2*M_PI/NSector;
    double dphiHG=0.5*dphisectorHG-M_PI/NSector;
    phimin_=iSector_*dphi-dphiHG;
    phimax_=phimin_+dphi+2*dphiHG;
    phimin_-=M_PI/NSector;
    phimax_-=M_PI/NSector;
    if (phimin_>M_PI) {
      phimin_-=2*M_PI;
      phimax_-=2*M_PI;
    }
    phioffset_=phimin_;

    layerdisk_=initLayerDisk(3);

    fullMatches_.resize(12,0);
    
    //FIXME should sort out constants here
    icorrshift_=7+idrinvbits+phi0bitshift-rinvbitshift-phiderbitshift;
    if (layerdisk_<3) {
      icorzshift_=-1-PS_zderL_shift;
    } else {
      icorzshift_=-1-SS_zderL_shift;      
    }
    phi0shift_=3;
    fact_=1;
    if (layerdisk_>=3&&layerdisk_<6) {
      fact_=(1<<(nbitszprojL123-nbitszprojL456));
      icorrshift_-=(10-nbitsrL456);
      icorzshift_+=(nbitszprojL123-nbitszprojL456+nbitsrL456-nbitsrL123);
      phi0shift_=0;
    }

    
    for(unsigned int iSeed=0;iSeed<12;iSeed++) {
      if (layerdisk_<6) {
	phimatchcut_[iSeed]=rphimatchcut[layerdisk_][iSeed]/(kphi1*rmean[layerdisk_]);
	zmatchcut_[iSeed]=zmatchcut[layerdisk_][iSeed]/kz;
      } else {
	rphicutPS_[iSeed]=rphicutPS[layerdisk_-6][iSeed]/(kphiproj123*kr);
	rphicut2S_[iSeed]=rphicut2S[layerdisk_-6][iSeed]/(kphiproj123*kr);
	rcut2S_[iSeed]=rcut2S[layerdisk_-6][iSeed]/krprojshiftdisk;
	rcutPS_[iSeed]=rcutPS[layerdisk_-6][iSeed]/krprojshiftdisk;
      }
    }

    if (iSector_==0&&layerdisk_<6&&settings_->writeTable()) {

      ofstream outphicut;
      outphicut.open(getName()+"_phicut.tab");
      outphicut << "{"<<endl;
      for(unsigned int seedindex=0;seedindex<12;seedindex++){
	if (seedindex!=0) outphicut<<","<<endl;
	outphicut << phimatchcut_[seedindex];
      }
      outphicut <<endl<<"};"<<endl;
      outphicut.close();

      ofstream outzcut;
      outzcut.open(getName()+"_zcut.tab");
      outzcut << "{"<<endl;
      for(unsigned int seedindex=0;seedindex<12;seedindex++){
	if (seedindex!=0) outzcut<<","<<endl;
	outzcut << zmatchcut_[seedindex];
      }
      outzcut <<endl<<"};"<<endl;
      outzcut.close();
    }

    if (iSector_==0&&layerdisk_>=6&&settings_->writeTable()) {

      ofstream outphicut;
      outphicut.open(getName()+"_PSphicut.tab");
      outphicut << "{"<<endl;
      for(unsigned int seedindex=0;seedindex<12;seedindex++){
	if (seedindex!=0) outphicut<<","<<endl;
	outphicut << rphicutPS_[seedindex];
      }
      outphicut <<endl<<"};"<<endl;
      outphicut.close();

      outphicut.open(getName()+"_2Sphicut.tab");
      outphicut << "{"<<endl;
      for(unsigned int seedindex=0;seedindex<12;seedindex++){
	if (seedindex!=0) outphicut<<","<<endl;
	outphicut << rphicut2S_[seedindex];
      }
      outphicut <<endl<<"};"<<endl;
      outphicut.close();

      ofstream outzcut;
      outzcut.open(getName()+"_PSrcut.tab");
      outzcut << "{"<<endl;
      for(unsigned int seedindex=0;seedindex<12;seedindex++){
	if (seedindex!=0) outzcut<<","<<endl;
	outzcut << rcutPS_[seedindex];
      }
      outzcut <<endl<<"};"<<endl;
      outzcut.close();

      outzcut.open(getName()+"_2Srcut.tab");
      outzcut << "{"<<endl;
      for(unsigned int seedindex=0;seedindex<12;seedindex++){
	if (seedindex!=0) outzcut<<","<<endl;
	outzcut << rcut2S_[seedindex];
      }
      outzcut <<endl<<"};"<<endl;
      outzcut.close();
    }
      
    for (unsigned int i=0; i<10; i++) {
      ialphafactinner_[i]= 
	(1<<alphashift)*krprojshiftdisk*half2SmoduleWidth/(1<<(nbitsalpha-1))/(rDSSinner[i]*rDSSinner[i])/kphiproj123;
      ialphafactouter_[i]= 
	(1<<alphashift)*krprojshiftdisk*half2SmoduleWidth/(1<<(nbitsalpha-1))/(rDSSouter[i]*rDSSouter[i])/kphiproj123;
    }
    
  }

  void addOutput(MemoryBase* memory,string output){
    if (writetrace) {
      cout << "In "<<name_<<" adding output to "<<memory->getName()
	   << " to output "<<output<<endl;
    }
    if (output.substr(0,8)=="matchout"){
      FullMatchMemory* tmp=dynamic_cast<FullMatchMemory*>(memory);
      assert(tmp!=0);
      unsigned int iSeed=getISeed(memory->getName());
      fullMatches_[iSeed]=tmp;
      return;
    }
    cout << "Count not fined output = "<<output<<endl;
    assert(0);
  }

  void addInput(MemoryBase* memory,string input){
    if (writetrace) {
      cout << "In "<<name_<<" adding input from "<<memory->getName()
	   << " to input "<<input<<endl;
    }
    if (input=="allstubin"){
      AllStubsMemory* tmp=dynamic_cast<AllStubsMemory*>(memory);
      assert(tmp!=0);
      allstubs_=tmp;
      return;
    }
    if (input=="allprojin"){
      AllProjectionsMemory* tmp=dynamic_cast<AllProjectionsMemory*>(memory);
      assert(tmp!=0);
      allprojs_=tmp;
      return;
    }
    if (input.substr(0,5)=="match" && input.substr(input.size()-2,2)=="in"){
      CandidateMatchMemory* tmp=dynamic_cast<CandidateMatchMemory*>(memory);
      assert(tmp!=0);
      matches_.push_back(tmp);
      return;
    }
    cout << getName()<<" Could not find input "<<input<<endl;
    assert(0);
  }

  void execute() {

    unsigned int countall=0;
    unsigned int countsel=0;

    Tracklet* oldTracklet=0;

    std::vector<std::pair<std::pair<Tracklet*,int>,std::pair<Stub*,L1TStub*> > > mergedMatches=mergeMatches(matches_);
    
    for(unsigned int j=0;j<mergedMatches.size();j++){
	
      if (settings_->debugTracklet()&&j==0) {
        cout << getName() <<" has "<<mergedMatches.size()<<" candidate matches"<<endl;
      }
      
      countall++;
	
      L1TStub* stub=mergedMatches[j].second.second;
      Stub* fpgastub=mergedMatches[j].second.first;
      Tracklet* tracklet=mergedMatches[j].first.first;

      //check that the matches are orderd correctly
      //allow equal here since we can have more than one cadidate match per tracklet projection
      if (oldTracklet!=0) {
	assert(oldTracklet->TCID()<=tracklet->TCID());
      }
      oldTracklet=tracklet;
      
      if (layerdisk_<6) {
	  
	//Integer calculation
	
	int ir=fpgastub->r().value();
       	int iphi=tracklet->fpgaphiproj(layerdisk_+1).value();
	int icorr=(ir*tracklet->fpgaphiprojder(layerdisk_+1).value())>>icorrshift_;	
	iphi+=icorr;
	
	int iz=tracklet->fpgazproj(layerdisk_+1).value();
	int izcor=(ir*tracklet->fpgazprojder(layerdisk_+1).value()+(1<<(icorzshift_-1)))>>icorzshift_;
	iz+=izcor;	

	int ideltaz=fpgastub->z().value()-iz;
	int ideltaphi=(fpgastub->phi().value()<<phi0shift_)-(iphi<<(phi0bitshift-1+phi0shift_)); 


	//Floating point calculations

	double phi=stub->phi()-phioffset_;
	double r=stub->r();
	double z=stub->z();

	
	if (useapprox) {
	  double dphi=Util::phiRange(phi-fpgastub->phiapprox(0.0,0.0));
	  assert(std::abs(dphi)<0.001);
	  phi=fpgastub->phiapprox(0.0,0.0);
	  z=fpgastub->zapprox();
	  r=fpgastub->rapprox();
	}
	
	if (phi<0) phi+=2*M_PI;

	double dr=r-tracklet->rproj(layerdisk_+1);
	assert(std::abs(dr)<drmax);

	double dphi=Util::phiRange(phi-(tracklet->phiproj(layerdisk_+1)+dr*tracklet->phiprojder(layerdisk_+1)));
	
       	double dz=z-(tracklet->zproj(layerdisk_+1)+dr*tracklet->zprojder(layerdisk_+1));
	
	double dphiapprox=Util::phiRange(phi-(tracklet->phiprojapprox(layerdisk_+1)+
						  dr*tracklet->phiprojderapprox(layerdisk_+1)));
	    
	double dzapprox=z-(tracklet->zprojapprox(layerdisk_+1)+
				   dr*tracklet->zprojderapprox(layerdisk_+1));
	
	int seedindex=tracklet->getISeed();

	assert(phimatchcut_[seedindex]>0);
	assert(zmatchcut_[seedindex]>0);

	bool truthmatch=tracklet->stubtruthmatch(stub);

	HistBase* hists=GlobalHistTruth::histograms();
	hists->FillLayerResidual(layerdisk_+1, seedindex,
				 dphiapprox*rmean[layerdisk_],
				 ideltaphi*kphi1*rmean[layerdisk_],
				 ideltaz*fact_*kz, dz,truthmatch);
    

	if (std::abs(dphi)>0.2 || std::abs(dphiapprox)>0.2 ) {
	  cout << "WARNING dphi and/or dphiapprox too large : "
	       <<dphi<<" "<<dphiapprox<<endl;
	}
	
	assert(std::abs(dphi)<0.2);
	assert(std::abs(dphiapprox)<0.2);


	
	if (writeResiduals) {
	  static ofstream out("layerresiduals.txt");
	  
	  double pt=0.003*3.8/std::abs(tracklet->rinv());
	  
	  out << layerdisk_+1<<" "<<seedindex<<" "<<pt<<" "<<ideltaphi*kphi1*rmean[layerdisk_]
	      <<" "<<dphiapprox*rmean[layerdisk_]
	      <<" "<<phimatchcut_[seedindex]*kphi1*rmean[layerdisk_]
	      <<"   "<<ideltaz*fact_*kz<<" "<<dz<<" "<<zmatchcut_[seedindex]*kz<<endl;	  
	}

	
	bool imatch=(std::abs(ideltaphi)<=(int)phimatchcut_[seedindex])&&(std::abs(ideltaz*fact_)<=(int)zmatchcut_[seedindex]);

	if (settings_->debugTracklet()) {
	  cout << getName()<<" imatch = "<<imatch<<" ideltaphi cut "<<ideltaphi
	       <<" "<<phimatchcut_[seedindex]
	       <<" ideltaz*fact cut "<<ideltaz*fact_<<" "<<zmatchcut_[seedindex]<<endl;
	}

	
	if (imatch) {
	  
	  countsel++;
	  
	  tracklet->addMatch(layerdisk_+1,ideltaphi,ideltaz,
			     dphi,dz,dphiapprox,dzapprox,
			     (fpgastub->phiregion().value()<<7)+fpgastub->stubindex().value(),
			     stub->r(),mergedMatches[j].second);
	  

	  if (settings_->debugTracklet()) {
	    cout << "Accepted full match in layer " <<getName()
		 << " "<<tracklet
		 << " "<<iSector_<<endl;	   
	  }
	  
	  fullMatches_[seedindex]->addMatch(tracklet,mergedMatches[j].second);

	}
      } else {  //disk matches
	
	
	//check that stubs and projections in same half of detector
	assert(stub->z()*tracklet->t()>0.0);

	int sign=(tracklet->t()>0.0)?1:-1;
	int disk=sign*(layerdisk_-5);
	assert(disk!=0);
	  
	//Perform integer calculations here

	int iz=fpgastub->z().value();
	int iphi=tracklet->fpgaphiprojdisk(disk).value();

	int shifttmp=t2bits+tbitshift+phi0bitshift+2-rinvbitshiftdisk-phiderdiskbitshift-PS_phiderD_shift;
	assert(shifttmp>=0);
	int iphicorr=(iz*tracklet->fpgaphiprojderdisk(disk).value())>>shifttmp;
	
	iphi+=iphicorr;

	int ir=tracklet->fpgarprojdisk(disk).value();
	
	  
	int shifttmp2=rprojdiskbitshift+t3shift-rderdiskbitshift;
	
	assert(shifttmp2>=0);
	int ircorr=(iz*tracklet->fpgarprojderdisk(disk).value())>>shifttmp2;
										
	ir+=ircorr;

	int ideltaphi=fpgastub->phi().value()*kphi/kphiproj123-iphi; 


	int irstub = fpgastub->r().value();
	int ialphafact=0;
	if(!stub->isPSmodule()){
	  assert(irstub<10);
	  if (abs(disk)<=2) {
	    ialphafact = ialphafactinner_[irstub];
	    irstub = rDSSinner[irstub]/kr;
	  } else {
	    ialphafact = ialphafactouter_[irstub];
	    irstub = rDSSouter[irstub]/kr;
	  }
	}

	//FIXME stub and projection r should not use different # bits...
	int ideltar=(irstub>>1)-ir;
	
	if (!stub->isPSmodule()) {	  
	  int ialphanew=fpgastub->alphanew().value();
	  int iphialphacor=((ideltar*ialphanew*ialphafact)>>alphashift);
	  ideltaphi+=iphialphacor;
	}
	
	//Perform floating point calculations here
	
	double phi=stub->phi()-phioffset_;
	double z=stub->z();
	double r=stub->r();
		
	if (useapprox) {
	  double dphi=Util::phiRange(phi-fpgastub->phiapprox(0.0,0.0));
	  assert(std::abs(dphi)<0.001);
	  phi=fpgastub->phiapprox(0.0,0.0);
	  z=fpgastub->zapprox();
	  r=fpgastub->rapprox();
	}

	if (phi<0) phi+=2*M_PI;

	double dz=z-sign*zmean[layerdisk_-6];
	
	if(std::abs(dz) > dzmax){
	  cout << __FILE__ << ":" << __LINE__ << " " << name_ << "_" << iSector_ << " " << tracklet->getISeed() << endl;
	  cout << "stub "<<stub->z() <<" disk "<<disk<<" "<<dz<<endl;
	  assert(std::abs(dz)<dzmax);
	}	
		  

	double phiproj=tracklet->phiprojdisk(disk)+dz*tracklet->phiprojderdisk(disk);

	double rproj=tracklet->rprojdisk(disk)+dz*tracklet->rprojderdisk(disk);
		
	double deltar=r-rproj;
	
	double dr=stub->r()-rproj;
	
	double dphi=Util::phiRange(phi-phiproj);

	double dphiapprox=Util::phiRange(phi-(tracklet->phiprojapproxdisk(disk)+
						  dz*tracklet->phiprojderapproxdisk(disk)));

	double drapprox=stub->r()-(tracklet->rprojapproxdisk(disk)+
				   dz*tracklet->rprojderapproxdisk(disk));
	
	double drphi=dphi*stub->r();
	double drphiapprox=dphiapprox*stub->r();


	
	if (!stub->isPSmodule()) {
	  double alphanew=stub->alphanew();
	  dphi+=dr*alphanew*half2SmoduleWidth/stub->r2();;
	  dphiapprox+=drapprox*alphanew*half2SmoduleWidth/stub->r2();
	  
	  drphi+=dr*alphanew*half2SmoduleWidth/stub->r();
	  drphiapprox+=dr*alphanew*half2SmoduleWidth/stub->r();
	}


	int seedindex=tracklet->getISeed();
	
	int idrphicut=rphicutPS_[seedindex];
	int idrcut=rcutPS_[seedindex]; 
	if (!stub->isPSmodule()) {
	  idrphicut=rphicut2S_[seedindex];
	  idrcut=rcut2S_[seedindex]; 
	}

	double drphicut=idrphicut*kphiproj123*kr;
	double drcut=idrcut*krprojshiftdisk;

	
	if (writeResiduals) {
	  static ofstream out("diskresiduals.txt");
	  
	  double pt=0.003*3.8/std::abs(tracklet->rinv());
	  
	  out << disk<<" "<<stub->isPSmodule()<<" "<<tracklet->layer()<<" "
	      <<abs(tracklet->disk())<<" "<<pt<<" "
	      <<ideltaphi*kphiproj123*stub->r()<<" "<<drphiapprox<<" "
	      <<drphicut<<" "
	      <<ideltar*krprojshiftdisk<<" "<<deltar<<" "
	      <<drcut<<" "
	      <<endl;	  
	}

	
	bool match=(std::abs(drphi)<drphicut)&&(std::abs(deltar)<drcut);
	
	bool imatch=(std::abs(ideltaphi*irstub)<idrphicut)&&(std::abs(ideltar)<idrcut);


	if (settings_->debugTracklet()) {
	  cout << "imatch match disk: "<<imatch<<" "<<match<<" "
	       <<std::abs(ideltaphi)<<" "<<drphicut/(kphiproj123*stub->r())<<" "
	       <<std::abs(ideltar)<<" "<<drcut/krprojshiftdisk<<" r = "<<stub->r()<<endl;
	}
		
	  
	if (imatch) {
	  
	  std::pair<Stub*,L1TStub*> tmp(fpgastub,stub);
	    
	  countsel++;
	  
	  if (settings_->debugTracklet()) {
	    cout << "MatchCalculator found match in disk "<<getName()<<endl;
	  }

          if(std::abs(dphi)>=0.25){
            cout<<"dphi "<<dphi<<"\n";
            cout<<"Seed / ISeed "<<tracklet->getISeed()<<"\n";
          }
          assert(std::abs(dphi)<0.25);
          assert(std::abs(dphiapprox)<0.25);

	  tracklet->addMatchDisk(disk,ideltaphi,ideltar,
				 drphi/stub->r(),dr,drphiapprox/stub->r(),drapprox,
				 stub->alpha(),
				 (fpgastub->phiregion().value()<<7)+fpgastub->stubindex().value(),
				 stub->z(),tmp);
	  if (settings_->debugTracklet()) {
	    cout << "Accepted full match in disk " <<getName()
		 << " "<<tracklet
		 << " "<<iSector_<<endl;	   
	  }

	  fullMatches_[seedindex]->addMatch(tracklet,mergedMatches[j].second);
	}
      }
      if (countall>=MAXMC) break;
    }


    if (writeMatchCalculator) {
      static ofstream out("matchcalculator.txt");
      out << getName()<<" "<<countall<<" "<<countsel<<endl;
    }

    
  }

  
  std::vector<std::pair<std::pair<Tracklet*,int>,std::pair<Stub*,L1TStub*> > > mergeMatches(vector<CandidateMatchMemory*>& candmatch) {
    
    std::vector<std::pair<std::pair<Tracklet*,int>,std::pair<Stub*,L1TStub*> > >  tmp;
    
    
    std::vector<unsigned int> indexArray;
    indexArray.reserve(candmatch.size());
    for (unsigned int i=0;i<candmatch.size();i++) {
      indexArray.push_back(0);
    }

    int bestIndex=-1;
    do {
      int bestSector=100;
      int bestTCID=(1<<16);
      bestIndex=-1;
      for (unsigned int i=0;i<candmatch.size();i++) {
	if (indexArray[i]>=candmatch[i]->nMatches()) {
	  // skip as we were at the end
	  continue;
	}
	int TCID=candmatch[i]->getFPGATracklet(indexArray[i])->TCID();
	int dSector=0;
	if (dSector>2) dSector-=NSector;
	if (dSector<-2) dSector+=NSector;
	assert(abs(dSector)<2);
	if (dSector==-1) dSector=2;
	if (dSector<bestSector) {
	  bestSector=dSector;
	  bestTCID=TCID;
	  bestIndex=i;
	}
	if (dSector==bestSector) {
	  if (TCID<bestTCID) {
	    bestTCID=TCID;
	    bestIndex=i;
	  }
	}
      }
      if (bestIndex!=-1) {
	tmp.push_back(candmatch[bestIndex]->getMatch(indexArray[bestIndex]));
	indexArray[bestIndex]++;
      }
    } while (bestIndex!=-1);
    
    if (layerdisk_<6) {
      
      int lastTCID=-1;
      bool error=false;
      
      //Allow equal TCIDs since we can have multiple candidate matches
      for(unsigned int i=1;i<tmp.size();i++){
	if (lastTCID>tmp[i].first.first->TCID()) {
	  cout << "Wrong TCID ordering for projections in "<<getName()<<" last "<<lastTCID<<" "<<tmp[i].first.first->TCID()<<endl;
	  error=true;
	} else {
	  lastTCID=tmp[i].first.first->TCID();
	}
      }
      
      if (error) {
	for(unsigned int i=1;i<tmp.size();i++){
	  cout << "Wrong order for in "<<getName()<<" "<<i<<" "<<tmp[i].first.first<<" "<<tmp[i].first.first->TCID()<<endl;
	}
      }
      
    }
    
    for (unsigned int i=0;i<tmp.size();i++) {
      if (i>0) {
	//This allows for equal TCIDs. This means that we can e.g. have a track seeded
	//in L1L2 that projects to both L3 and D4. The algorithm will pick up the first hit and
	//drop the second

	assert(tmp[i-1].first.first->TCID()<=tmp[i].first.first->TCID());

      }
    }
    
    return tmp;
  }

    
private:

  unsigned int layerdisk_;

  int fact_;
  int icorrshift_;
  int icorzshift_;
  int phi0shift_;
  double phimin_;
  double phimax_;
  double phioffset_;

  unsigned int phimatchcut_[12];
  unsigned int zmatchcut_[12];
  unsigned int rphicutPS_[12];
  unsigned int rphicut2S_[12];
  unsigned int rcutPS_[12];
  unsigned int rcut2S_[12];

  int ialphafactinner_[10];
  int ialphafactouter_[10];
  
  AllStubsMemory* allstubs_;
  AllProjectionsMemory* allprojs_;

  vector<CandidateMatchMemory*> matches_;

  vector<FullMatchMemory*> fullMatches_;

};

#endif

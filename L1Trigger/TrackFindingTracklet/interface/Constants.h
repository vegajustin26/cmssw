#ifndef L1Trigger_TrackFindingTracklet_interface_Constants_h
#define L1Trigger_TrackFindingTracklet_interface_Constants_h

static const bool geomTkTDR=false; // false => newest T14 tracker, true => "TDR" (T5/T6 tracker, D21/D11/D17 CMS geometries)

//Bits used to store track parameter in tracklet
static const int nbitsrinv=14;
static const int nbitsphi0=18;
static const int nbitsd0=13;
static const int nbitst=14;
static const int nbitsz0=10;

//track and tracklet parameters
const int rinv_shift = -8;  // Krinv = 2^shift * Kphi/Kr
const int phi0_shift = 1;   // Kphi0 = 2^shift * Kphi
const int t_shift    = -10; // Kt    = 2^shift * Kz/Kr
const int z0_shift   = 0;   // Kz0   = 2^shift * kz

//projections are coarsened from global to stub precision  

//projection to R parameters
const int PS_phiL_shift = 0;   // phi projections have global precision in ITC
const int SS_phiL_shift = 0;   
const int PS_zL_shift   = 0;   // z projections have global precision in ITC
const int SS_zL_shift   = 0;

const int PS_phiderL_shift = -5;   // Kderphi = 2^shift * Kphi/Kr
const int SS_phiderL_shift = -5; 
const int PS_zderL_shift   = -7;  // Kderz = 2^shift * Kz/Kr
const int SS_zderL_shift   = -7;  
  
//projection to Z parameters
const int PS_phiD_shift = 3;   
const int SS_phiD_shift = 3;   
const int PS_rD_shift   = 1;   // a bug?! coarser by a factor of two then stubs??
const int SS_rD_shift   = 1;

const int PS_phiderD_shift = -4; //Kderphidisk = 2^shift * Kphi/Kz
const int SS_phiderD_shift = -4; 
const int PS_rderD_shift   = -6;  //Kderrdisk = 2^shift * Kr/Kz
const int SS_rderD_shift   = -6;  

//numbers needed for matches & fit, unclear what they are.
static const int idrinvbits=19;
static const int phi0bitshift=1;
static const int rinvbitshift=13;
static const int tbitshift=9;
static const int z0bitshift=0;
static const int phiderbitshift=7;
static const int zderbitshift=6;
static const int t2bits=23;
static const int t3shift=8;
static const int rinvbitshiftdisk=13; 
static const int rprojdiskbitshift=6;
static const int phiderdiskbitshift=20;
static const int rderdiskbitshift=7;


static const int phiresidbits=12; 
static const int zresidbits=9;
static const int rresidbits=7;

//Trackfit
static const int fitrinvbitshift=9;  //6 OK?
static const int fitphi0bitshift=6;  //4 OK?
static const int fittbitshift=10;     //4 OK? //lower number gives rounding problems
static const int fitz0bitshift=8;    //6 OK?

//r correction bits
static const int rcorrbits=6;

static const int chisqphifactbits=14;
static const int chisqzfactbits=14;


//
//Geometry 
//

//These define the length scale for both r and z
static const double zlength=120.0;
static const double rmaxdisk=120.0;


// these assume either "TDR" tracker geometry (T5 or T6), or otherwise most recent T14 tracker 
// T5: http://cms-tklayout.web.cern.ch/cms-tklayout/layouts/recent-layouts/OT616_200_IT404/layout.html
// T14: http://cms-tklayout.web.cern.ch/cms-tklayout/layouts/recent-layouts/OT616_200_IT404/layout.html

static const double rmeanL1=geomTkTDR?(rmaxdisk*858)/4096:(rmaxdisk*851)/4096;
static const double rmeanL2=geomTkTDR?(rmaxdisk*1279)/4096:(rmaxdisk*1269)/4096;
static const double rmeanL3=geomTkTDR?(rmaxdisk*1795)/4096:(rmaxdisk*1784)/4096;
static const double rmeanL4=geomTkTDR?(rmaxdisk*2347)/4096:(rmaxdisk*2347)/4096;
static const double rmeanL5=geomTkTDR?(rmaxdisk*2937)/4096:(rmaxdisk*2936)/4096;
static const double rmeanL6=geomTkTDR?(rmaxdisk*3783)/4096:(rmaxdisk*3697)/4096;

static const double zmeanD1=(zlength*2239)/2048;
static const double zmeanD2=(zlength*2645)/2048;
static const double zmeanD3=(zlength*3163)/2048;
static const double zmeanD4=(zlength*3782)/2048;
static const double zmeanD5=(zlength*4523)/2048;


static const double rmindiskvm=22.5;
static const double rmaxdiskvm=67.0;

static const double rmaxdiskl1overlapvm=45.0;
static const double rmindiskl2overlapvm=40.0;
static const double rmindiskl3overlapvm=50.0;

static const double half2SmoduleWidth=4.57;


static const double drmax=rmaxdisk/32.0;
static const double dzmax=zlength/32.0;

static const double rmean[6]={rmeanL1,rmeanL2,rmeanL3,rmeanL4,rmeanL5,rmeanL6};

static const double zmean[5]={zmeanD1,zmeanD2,zmeanD3,zmeanD4,zmeanD5};

static const double rcrit=55.0;

static const double rmaxL6=rmeanL6+drmax; 


static const double zmaxD5=zmeanD5+dzmax; 

static const double two_pi=2*M_PI;

static const double ptcut=1.91; //Minimum pt // ONLY USED IN IMATH STUFF 
static const double rinvcut=0.01*0.3*3.8/ptcut; //0.01 to convert to cm-1  // ONLY USED IN IMATH STUFF 

static const double ptcutte=1.8; //Minimum pt in TE
static const double rinvcutte=0.01*0.3*3.8/ptcutte; //0.01 to convert to cm-1 in TE
static const double bendcut=1.25;  //Obsolete should be removed
static const double bendcutdisk=1.25;   //Obsolete should be removed
static const double z0cut=15.0;
static const double mecut=2.0;
static const double mecutdisk=1.5;

static const double bendcut_[2][8]={ {1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25},     //inner
				     {1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25} };     //outer

static const unsigned int NSector=9; 

static const double rinvmax=0.01*0.3*3.8/2.0; //0.01 to convert to cm-1

static const double dphisectorHG=2*M_PI/NSector+2*fmax(std::abs(asin(0.5*rinvmax*rmean[0])-asin(0.5*rinvmax*rcrit)),
						std::abs(asin(0.5*rinvmax*rmean[5])-asin(0.5*rinvmax*rcrit)));

static const double phicritmin=0.5*dphisectorHG-M_PI/NSector;
static const double phicritmax=dphisectorHG-0.5*dphisectorHG+M_PI/NSector;

static const double dphicritmc=0.005; //lose for MC
static const double phicritminmc=phicritmin-dphicritmc;
static const double phicritmaxmc=phicritmax+dphicritmc;

// Obsolete - only used in TrackletCalculatorDisplaced (Ryd - 2020-01-16)
static const int iphicritminmc=9253;
static const int iphicritmaxmc=56269;

static const unsigned int NLONGVMBITS=3; 

static const unsigned int NLONGVMBINS=(1<<NLONGVMBITS);

static const double dphisector=2*M_PI/NSector;


//Minimal ranges for track parameters
static const double maxrinv=0.006;
static const double maxd0=10.;


//These are constants defining global coordinate system

static const double kphi=dphisectorHG/(1<<14); //nbitsphistubL123
static const double kphi1=dphisectorHG/(1<<17); //nbitsphistubL456
static const double kz=2*zlength/(1<<12);   //nbitszL123
static const double kr=rmaxdisk/(1<<12);  //nrbitsdisk
static const double kd0 = 2*maxd0/(1<<nbitsd0);


//constants derivative from the above
static double krinvpars, kphi0pars, kd0pars, ktpars, kz0pars;
static double kphiproj123, kphiproj456, kzproj, kphider, kzder;
static double krprojshiftdisk, kphiprojdisk,krprojderdisk;
static double krdisk,krprojderdiskshift, kzpars;

#endif


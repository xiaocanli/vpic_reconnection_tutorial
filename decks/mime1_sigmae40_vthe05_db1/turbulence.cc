///////////////////////////////////////////////////////////////////////////////
//
//  Turbulence Problem
//  The boundary condition can be open boundary or periodic.
//
//  Use head version of VPIC.
//
///////////////////////////////////////////////////////////////////////////////
/* #include <H5Part.h> */
#include <math.h>
#include <list>
#include <iterator>
#include "vpic/dumpmacros.h"
#include "injection.hh"   //  Subroutine to compute re-injection velocity
#include "tracer.hh"
#include "species_advance/species_advance_efield.h"
#include "hdf5.h"

#define NUMFOLD (rank()/32)

// structure to hold the data for energy diagnostics
struct edata {
  species_id sp_id;       /* species id */
  double       vth;       /* thermal energy */
  char  fname[256];       /* file to save data */
};

// electric and magnetic fields
typedef struct emf {
  float ex, ey, ez;     // Electric field and div E error
  float cbx, cby, cbz;  // Magnetic field and div B error
} emf_t;

// Whether to simulation turbulence
#define TURBULENCE_SIMULATION

// Whether to use a Harris sheet. If not defined, a forcefree sheet is used
/* #define HARRIS_SHEET */

// Whether only electrons carry current (for forcefree sheet only)
/* #define ELECTRONS_CARRRY_CURRENT */

// Whether to use local perturbation. If not defined, global perturbation is used
/* #define LOCAL_PERTURBATION */

// Whether to use open boundary conditions
/* #define OPEN_BOUNDARY_CONDITIONS */

// Whether to dump three successive time steps of field and hydro data
/* #define DUMP_THREE_SUCCESSIVE_STEPS */

// Whether to dump hydro quantities for different energy band
/* #define ENERGY_BAND_DIAGNOSTICS */

// Whether to perform diagnostics on turbulence mixing
/* #define TURBULENCE_MIXING_DIAGNOSTICS */

// Whether to perturb initially using large-scale waves
// You still need to fine-tune the wave by yourself
/* #define WAVE_PERTURBATION */

// Whether to use relativistic plasma
/* #define RELATIVISTIC_PLASMA */

// Whether to have tracer particles advanced with part of electric field
/* #define SPECIES_ADVANCED_WITH_PART_EFIELD */

// Whether to have tracer particles advanced without perpendicular electric field
/* #define SPECIES_ADVANCED_WITHOUT_EPERP */

// Whether to check whether particles crossing regions with E>B
/* #define CHECK_PARTICLE_CROSSING_EGTB */

// Whether to enable radiation reaction force on the particles
/* #define ENABLE_RADIATION_REACTION */

// Whether to use HDF5 format for dumping fields and hydro
#define DUMP_WITH_HDF5

#ifdef DUMP_WITH_HDF5
// Deck only works if VPIC was build with HDF support. Check for that:
#ifndef VPIC_ENABLE_HDF5
#error "VPIC_ENABLE_HDF5" is required
#endif
#endif

// naming convention for the dump files
#define HYDRO_FILE_FORMAT "hydro/%d/T.%d/%s.%d.%d"
#define SPEC_FILE_FORMAT "hydro/%d/T.%d/spectrum-%s.%d.%d"

// directory on scratch space for dumping data
#define DUMP_DIR_FORMAT "./%s"

// array access
#define LOCAL_CELL_ID(x,y,z) VOXEL(x, y, z, grid->nx, grid->ny, grid->nz)
#define HYDRO(x,y,z) hi[LOCAL_CELL_ID(x,y,z)]
#define HYDRO_TOT(x,y,z) htot[LOCAL_CELL_ID(x,y,z)]

// Vadim's in-line average
#define ALLOCATE(A,LEN,TYPE)                                    \
  if ( !((A)=(TYPE *)malloc((size_t)(LEN)*sizeof(TYPE))) )      \
    ERROR(("Cannot allocate."));

void
checkpt_subdir( const char * fbase, int tag )
{
  char fname[256];
  if( !fbase ) ERROR(( "NULL filename base" ));
  sprintf( fname, "%s/%i/restore.0.%i", fbase, tag, world_rank );
  if( world_rank==0 ) log_printf( "*** Checkpointing to \"%s\"\n", fbase );
  checkpt_objects( fname );
}

begin_globals {

  int restart_interval;
  int energies_interval;
  int spectrum_interval;
  int fields_interval;
  int ehydro_interval;
  int Hhydro_interval;
  int eparticle_interval;
  int Hparticle_interval;
  int quota_check_interval;  //  How frequently to check if quote exceeded

  int rtoggle;               // enables save of last two restart dumps
  double quota_sec;          // Run quota in seconds
  double b0;                 // B0
  double bg;                 // Guide field
  double v_A;                // Alfven speed
  double c0;                 // light speed
  double topology_x;         // domain topology
  double topology_y;
  double topology_z;

  // Variables for new output format
  DumpParameters fdParams;
  std::vector<DumpParameters *> outputParams;
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
  DumpParameters eTopdParams;
  DumpParameters eBotdParams;
  DumpParameters HTopdParams;
  DumpParameters HBotdParams;
  edata edeTop;
  edata edeBot;
  edata edHTop;
  edata edHBot;
#else
  DumpParameters hedParams;
  DumpParameters hHdParams;
  edata ede;            // electron species information
  edata edH;            // ion species information
#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
  edata ede_wo_epara;
  edata edH_wo_epara;
  edata ede_wo_eparay;
  edata edH_wo_eparay;
  edata ede_wo_egtb;
  edata edH_wo_egtb;
#ifdef SPECIES_ADVANCED_WITHOUT_EPERP
  edata ede_wo_eperp;
  edata edH_wo_eperp;
#endif
#endif
#ifdef CHECK_PARTICLE_CROSSING_EGTB
  edata ede_egtb;
  edata edH_egtb;
#endif
#endif

  // particle spectrum and energy band diagnostics
  double emax_band;     // maximum energy for energy band diagnostics
  double emin_band;     // minimum energy for energy band diagnostics
  int nbands;           // # of energy bands
  double emax_spect;    // maximum energy for energy spectrum diagnostics
  double emin_spect;    // minimum energy for energy spectrum diagnostics
  int nbins;            // # of energy bins for energy spectrum diagnostics
  int nx_zone;          // # of cells / zone along x for energy spectrum diagnostics
  int ny_zone;          // # of cells / zone along y for energy spectrum diagnostics
  int nz_zone;          // # of cells / zone along z for energy spectrum diagnostics

  int restart_step;     // time step for restart dump

  int nsp;              // number of species
  int nsp_efield;       // number of species for electric field diagnostics

  // // time-averaging diagnostic
  // int dis_nav;         // number of steps to average over
  // int dis_interval;    // number of steps between outputs.
  // int dis_iter;        // iteration count. 0 means we are not averaging at the moment
  // int dis_begin_int;   // the first time step of the interval
  // char fields_dir[128];
  // char hydro_dir[128];
  // char restart_avg_dir[128];

  // particle tracking
  int tracer_interval;         // tracer info is saved or dumped every tracer_interval steps
  int tracer_pass1_interval;   // tracer interval for the 1st run. A multiple of tracer_interval
  int tracer_pass2_interval;   // tracer interval for the re-run. A multiple of tracer_interval
  int tracer_file_interval;    // interval when multiple frames of tracers are saved in the same file
  int emf_at_tracer;           // 0 or 1, electric and magnetic fields at tracer
  int hydro_at_tracer;         // 0 or 1, hydro fields at tracer
  int ve_at_tracer;            // 0 or 1, electron bulk velocity at tracer
  int num_tracer_fields_add;   // additional number of tracer fields
  int particle_tracing;
  int particle_select;
  int particle_tracing_start;
  int dump_traj_directly;
  species_t *tracers_list;
#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
  species_t *species_list_efield;
#endif
#ifdef CHECK_PARTICLE_CROSSING_EGTB
  species_t *species_list_egtb;
#endif
  int tag;
  double mi_me;
  int Ntracer;

  int stride_particle_dump;  // stride for particle dump

#ifdef OPEN_BOUNDARY_CONDITIONS
  // Variables for Open BC Model
  double nb;        // Background density
  double vth[2];    // Thermal velocity of Harris components
  double vthb[2];   // Thermal velocity of background components
  double q[2];      // Species charge
  double L_de;      // Initial Harris sheet thickness
  double uf[2];     // Initial Fluid Drift in Harris
  double nfac;      // Normalization factor to convert particles per cell into density
  double rin[3];    // Relaxation parameter for inflow boundary moments
  double rout[3];   // Relaxation parameter for outlfow boundary moments
  double sort[2];   // Intervals where we know particles are sorted
  double edrive;    // Drive field for inflow boundary
  double tdrive;
  int left,right,top,bottom;  // Keep track of boundary domains
  double *nbot, *ubot, *pbot, *bbot, *fbot;           // Moments for bottom injectors
  double *ntop, *utop, *ptop, *btop, *ftop;           // Moments for top injectors
  double *nleft, *uleft, *pleft, *bleft, *fleft;      // Moments for left injectors
  double *nright, *uright, *pright, *bright, *fright; // Moments for right injectors
#endif
};


begin_initialization {

  // use natural PIC units
  double ec   = 1;        // Charge normalization
  double me   = 1;        // Mass normalization
  double c    = 1;        // Speed of light
  double de   = 1;        // Length normalization (electron inertial length)
  double eps0 = 1;        // Permittivity of space

  double cfl_req   = 0.9; // How close to Courant should we try to run
  double wpedt_max = 0.3; // How big a timestep is allowed if Courant is not too restrictive
  int rng_seed     = 1;   // Random number seed increment

  // particle tracking
  int particle_tracing = 1; // 0: notracing, 1: forward tracing 2: tracing from particle files
  int particle_select = 50000; // track one every particle_select particles
  int particle_tracing_start = 0; // the time step that particle tracking is triggered
                                  // this should be set to 0 for Pass1 and 2
  int nframes_per_tracer_file = 1000; // number of frames of tracers saved in one tracer file
  int dump_traj_directly = 0;     // dump particle trajectories in 1st pass
  int emf_at_tracer = 1;          // electric and magnetic fields at tracer
  int hydro_at_tracer = 0;        // hydro fields at tracer
  int ve_at_tracer = 0;           // electron bulk velocity
  int num_emf = 0;                // number of electric and magnetic field, change between passes
  int num_hydro = 0;              // number of hydro fields, change between passes
  if (emf_at_tracer) num_emf = 6; // Make sure the sum of these two == TRACER_NUM_ADDED_FIELDS
  if (hydro_at_tracer) {
    num_hydro = 5; // single fluid velocity, electron and ion number densities
    if (ve_at_tracer) num_hydro += 3;
  } else {
    ve_at_tracer = 0;
  }

#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
  int particle_select_part_efield = 20;
#endif
#ifdef CHECK_PARTICLE_CROSSING_EGTB
  int particle_select_egtb = 20;
#endif

  // Physics parameters
  double mi_me   = 1.0;  // Ion mass / electron mass
  double L_di;            // Half thickness of the current sheet

#ifdef RELATIVISTIC_PLASMA
  double sigma = mi_me * 0.1; // Magnetization factor (2*magnetic energy)/(electron rest energy)
  // Make sure that current sheet is thick enough to support the current density
#ifdef ELECTRONS_CARRRY_CURRENT
  L_di = 1.2 * sqrt(sigma) / sqrt(mi_me);
#else
  L_di = 0.6 * sqrt(sigma) / sqrt(mi_me);
#endif

#else
  L_di = 1.0;               // Sheet thickness / ion inertial length
#endif

  double Ti_Te   = 1.0;   // Ion temperature / electron temperature
#ifdef HARRIS_SHEET
  double nb_n0   = 0.23;    // background plasma density
  double Tbe_Te  = 0.76;    // Ratio of background electron temperature to Harris electron temperature
  double Tbi_Ti  = 0.76;    // Ratio of background ion temperature to Harris ion temperature
#else
  double nb_n0   = 1.0;   // background plasma density
  double Tbe_Te  = 1.0;   // Ratio of background electron temperature to Harris electron temperature
  double Tbi_Ti  = 1.0;   // Ratio of background ion temperature to Harris ion temperature
#ifdef RELATIVISTIC_PLASMA
  double Te      = 1.0;   // Electron temperature
#else
  double vthe    = 0.5;   // Electron thermal speed over c
#endif
#endif // #ifdef HARRIS_SHEET

#ifdef RELATIVISTIC_PLASMA
  double wpe_wce = 1/sqrt(sigma); // electron plasma freq / electron cyclotron freq
  double factor = 6.0;            // for generating Maxwell–Jüttner distribution
#else
  double wpe_wce = 1.0/sqrt(40.0); // electron plasma freq / electron cyclotron freq
#endif
  double bg = 0.0;        // guide field strength
  double theta   = 0.0;   // B0 = Bx
  double taui    = 100.0/(wpe_wce*mi_me); // simulation wci's to run

  double quota   = 47.5;          // run quota in hours
  double quota_sec = quota*3600;  // Run quota in seconds

  double pi = 3.1415927;
  double cs   = cos(theta/180.0*pi);
  double sn   = sin(theta/180.0*pi);

  //derived quantities
  double mi = me*mi_me;                 // Ion mass
#ifdef HARRIS_SHEET
  double Te = me*c*c/(2*eps0*wpe_wce*wpe_wce*(1+Ti_Te)); // Electron temperature
  double Ti = Te*Ti_Te;                                  // Ion temperature
  double vthe = sqrt(Te/me);                             // Electron thermal velocity
  double vthi = sqrt(Ti/mi);                             // Ion thermal velocity
  double vtheb = sqrt(Tbe_Te) * vthe;                    // Background electron thermal velocity
  double vthib = sqrt(Tbi_Ti) * vthi;                    // Background ion thermal velocity
#else
#ifdef RELATIVISTIC_PLASMA
  double vthe = sqrt(Te/me);      // Electron momentum
#else
  double Te   = me * vthe * vthe; // Electron temperature
#endif
  double vthi = vthe*sqrt(Ti_Te/mi_me); // Ion thermal velocity
  double Ti = Te * Ti_Te;               // Ion temperature
  double vtheb = vthe;                  // Background electron thermal velocity
  double vthib = vthi;                  // Background ion thermal velocity
#endif

  double wci  = 1.0/(mi_me*wpe_wce);    // Ion cyclotron frequency
  double wce  = wci*mi_me;              // Electron cyclotron frequency
  double wpe  = wce*wpe_wce;            // electron plasma frequency
  double wpi  = wpe/sqrt(mi_me);        // ion plasma frequency
  double di   = c/wpi;                  // ion inertial length
  double L    = L_di*di;                // Sheet thickness in c/wpe
  double rhoi_L = sqrt(Ti_Te/(1.0+Ti_Te))/L_di;
  double v_A  = (wci/wpi)/sqrt(nb_n0);  // based on nb

#ifdef OPEN_BOUNDARY_CONDITIONS
  // Parameters for Open BC model
  //double rin[3] =  {0.005, 0.005, 0.00}; //  Relaxation - density, velocity + particle flux, pressure tensor
  //double rout[3] = {0.01, 0.01, 0.00};   //  Relaxation - density, velocity + particle flux, pressure tensor
  //double edrive = 0.04*v_A/(wpe_wce);    //  Setting edrive = 0 will give undriven limit
  double rin[3] =  {0.000, 0.06, 0.00}; //  Relaxation - density, velocity + particle flux, pressure tensor
  double rout[3] = {0.002, 0.002, 0.002};   //  Relaxation - density, velocity + particle flux, pressure tensor
  double edrive = 0.08*v_A/(wpe_wce);    //  Setting edrive = 0 will give undriven limit
  //double edrive = 0.007;
  double tdrive = 600.0;                 //  Time scale to turn on driving electric field
#endif

  // Numerical parameters
  double ion_sort_interval = 25;        // Injector moments are also updated at this internal
  double electron_sort_interval = 25;   // Injector moments are also updated at this internal
  double nppc  =  100;                  // Average number of macro particle per cell per species

  double nx = 1024;
  double ny = 1024;
  double nz = 1;

  double Lx  = 256.0/sqrt(mi_me)*di;   // size of box in x dimension
  double Ly  = ny * Lx / nx;           // size of box in y dimension
  double Lz  = nz * Lx / nx;           // size of box in z dimension

  double topology_x = 8;  // Number of domains in x, y, and z
  double topology_y = 8;
  double topology_z = 1;

  /* double nx = 128; */
  /* double ny = 128; */
  /* double nz = 128; */

  /* double Lx  = 64.0/sqrt(mi_me)*di;   // size of box in x dimension */
  /* double Ly  = ny * Lx / nx;           // size of box in y dimension */
  /* double Lz  = nz * Lx / nx;           // size of box in z dimension */

  /* double topology_x = 8;  // Number of domains in x, y, and z */
  /* double topology_y = 8; */
  /* double topology_z = 1; */

  double hx = Lx/nx;
  double hy = Ly/ny;
  double hz = Lz/nz;

  double b0 = me*c*wce/ec;                // Asymptotic magnetic field strength
  double n0 = me*eps0*wpe*wpe/(ec*ec);    // Peak electron (ion) density
  double Ne = nppc*nx*ny*nz;              // total macro electrons in box
#ifdef HARRIS_SHEET
  double Npe_sheet     = 2*n0*Lx*Ly*L*tanh(0.5*Lz/L); // number of physical electrons in sheet
  double Npe_back      = nb_n0*n0*Ly*Lz*Lx;           // Number of physical electrons in background
  double Npe           = Npe_sheet + Npe_back;
  double Ne_sheet      = Ne*Npe_sheet/Npe;
  double Ne_back       = Ne*Npe_back/Npe;
  Ne_sheet  = trunc_granular(Ne_sheet,nproc()); // Make it divisible by number of processors
  Ne_back  = trunc_granular(Ne_back,nproc());   // Make it divisible by number of processors
  Ne = Ne_sheet + Ne_back;
  double weight_s = Npe_sheet/Ne_sheet;    // particle weight in the sheet
  double weight_b = Npe_back/Ne_back;      // particle weight in the background
  double qe_s = -ec*weight_s;              // Charge per macro electron in the sheet
  double qe_b = -ec*weight_b;              // Charge per macro electron in the background
  double qi_s =  ec*weight_s;              // Charge per macro ion in the sheet
  double qi_b =  ec*weight_b;              // Charge per macro ion in the background

  double vdri = 2*c*Ti/(ec*b0*L);          // Ion drift velocity
  double vdre = -vdri/(Ti_Te);             // electron drift velocity
  double gdri = 1/sqrt(1-vdri*vdri/(c*c)); // gamma of ion drift frame
  double gdre = 1/sqrt(1-vdre*vdre/(c*c)); // gamma of electron drift frame
  double udri = vdri*gdri;                 // 4-velocity of ion drift frame
  double udre = vdre*gdre;                 // 4-velocity of electron drift frame
  double tanhf = tanh(0.5*Lz/L);
#else
  double Np = n0*Lx*Ly*Lz;                 // total number of physical electrons
  Ne = trunc_granular(Ne,nproc());         // Make it divisible by number of processors
  double weight = Np/Ne;                   // particle weight
  double qe = -ec*weight;                  // Charge per macro electron
  double qi =  ec*weight;                  // Charge per macro ion
  double qe_s = qe;                        // Charge per macro electron in the sheet
  double qe_b = qe;                        // Charge per macro electron in the background
  double qi_s = qi;                        // Charge per macro ion in the sheet
  double qi_b = qi;                        // Charge per macro ion in the background
#endif
  double nfac = qi_s/(hx*hy*hz);           // Convert density to particles per cell
  double Ntracer = Ne / particle_select;   // Number of particle tracers for each species
  Ntracer = trunc_granular(Ntracer, nproc());
  /* double Lpert = 1.6*Lx;                       // wavelength of perturbation */
  double Lpert = Lx;                       // wavelength of perturbation
  double dbz = 0.02*b0;                    // Perturbation in Bz rel. to Bo (Only change here)
  double dbx = -dbz*Lpert/(2.0*Lz);        // Set Bx perturbation so that div(B) = 0

  // parameters for turbulence simulations
#ifdef TURBULENCE_SIMULATION
  double db_rms = b0;                   // Turbulence amplitude db/b0
  int nmodes = 8;                       // Number of modes
  int nmodes_para = 2;                  // 3D only: Number of modes parallel to the mean field along z
  double db_amp;                        // Wave amplitude for each mode
  if (nz > 1) { // 3D
    db_amp = db_rms / (nmodes * 2) / sqrt(2 * nmodes_para); // both forward and backward
  } else {
    db_amp = db_rms / (nmodes * 2); // both forward and backward
  }
#endif

#ifdef ENABLE_RADIATION_REACTION
  double re_de0 = 5.3068E-19;  // classical radius of electrons normalized by de if plasma number density n = 1cm^-3
  double gamma_e   = 1.0E4;    // electron Lorentz factor that we are interested in
  double erho_rad  = 0.0;      // energy density of the radiation field
  double tauc_taua_targeted = 200.0;
  double tau_acc  = 1.0 / (0.0005 * b0); // acceleration time scale, assuming 1st-order Fermi
  double tau_cool = 1.5 / (gamma_e * b0 * b0 * re_de0); // radiation cooling time scale
  double tauc_taua = tau_cool / tau_acc;
  double re_de = tauc_taua * re_de0 / tauc_taua_targeted;
  double gamma_min = 1E0;      // minimum Lorentz factor to turn on radiation reaction
#endif

  // Determine the time step
  double dg = courant_length(Lx,Ly,Lz,nx,ny,nz); // courant length
  double dt = cfl_req*dg/c;                      // courant limited time step
  if( wpe*dt>wpedt_max) dt=wpedt_max/wpe;        // override timestep if plasma frequency limited

  int restart_interval = int(10000.0/(wpe*dt));
  /* int restart_interval = 50;                     // for testing */
  int energies_interval = 100;
  /* int energies_interval = 1;                     // for testing */
  int interval = int(20.0/(wpe*dt));
  /* int interval = 50;                             // for testing */
  /* int tracer_interval = int(0.25/(wpe*dt)); */
  int tracer_interval = 1;                      // for testing
  if (tracer_interval == 0) tracer_interval = 1;
  int tracer_pass1_interval = tracer_interval;
  int tracer_pass2_interval = tracer_interval;
  int spectrum_interval  = interval;
  int fields_interval    = interval;
  int ehydro_interval    = interval;
  int Hhydro_interval    = interval;
  int eparticle_interval = 200000*interval;
  int Hparticle_interval = 200000*interval;
  int quota_check_interval = 100;

  // particle spectrum and energy band diagnostics
  double emax_band = 120.0;     // maximum energy for energy band diagnostics
  double emin_band = 0.0;       // minimum energy for energy band diagnostics
  int nbands = 0;               // # of energy bands
  double emax_spect = 1.0E6;    // maximum energy for energy spectrum diagnostics
  double emin_spect = 1.0E-4;   // minimum energy for energy spectrum diagnostics
  int nbins = 1000;             // # of energy bins for energy spectrum diagnostics
  int nx_zone = 8;             // # of cells / zone along x for energy spectrum diagnostics
                                // MAKE SURE that (nx / topology_x) is divisible by nx_zone
  int ny_zone = 8;             // # of cells / zone along y for energy spectrum diagnostics
                                // MAKE SURE that (ny / topology_y) is divisible by ny_zone
  int nz_zone = 128;           // # of cells / zone along z for energy spectrum diagnostics
                                // MAKE SURE that (nz / topology_z) is divisible by nz_zone

  // particle dump
  int stride_particle_dump = 40; // stride for particle dump

  // Determine which domains area along the boundaries - Use macro from
  // grid/partition.c

# define RANK_TO_INDEX(rank,ix,iy,iz) BEGIN_PRIMITIVE {                 \
    int _ix, _iy, _iz;                                                  \
    _ix  = (rank);                /* ix = ix+gpx*( iy+gpy*iz ) */       \
    _iy  = _ix/int(topology_x);   /* iy = iy+gpy*iz */                  \
    _ix -= _iy*int(topology_x);   /* ix = ix */                         \
    _iz  = _iy/int(topology_y);   /* iz = iz */                         \
    _iy -= _iz*int(topology_y);   /* iy = iy */                         \
    (ix) = _ix;                                                         \
    (iy) = _iy;                                                         \
    (iz) = _iz;                                                         \
  } END_PRIMITIVE

  int ix, iy, iz ;
  RANK_TO_INDEX( int(rank()), ix, iy, iz );
#ifdef OPEN_BOUNDARY_CONDITIONS
  int left = 0, right = 0, top = 0, bottom = 0;
  if (ix == 0) left=1;
  if (ix == topology_x-1) right=1;
  if (iz == 0) bottom=1;
  if (iz == topology_z-1) top=1;
#endif

  ///////////////////////////////////////////////
  // Setup high level simulation parameters
  sim_log("Setting up high-level simulation parameters. ");
  num_step             = int(taui/(wci*dt));
  /* num_step             = 100;  // for testing */
  status_interval      = 200;
  sync_shared_interval = status_interval/2;
  clean_div_e_interval = status_interval/2;
  clean_div_b_interval = status_interval/2;
  field_interval = fields_interval;
  hydro_interval = hydro_interval;

  global->mi_me                = mi_me;
  global->restart_interval     = restart_interval;
  global->energies_interval    = energies_interval;
  global->spectrum_interval    = spectrum_interval;
  global->fields_interval      = fields_interval;
  global->ehydro_interval      = ehydro_interval;
  global->Hhydro_interval      = Hhydro_interval;
  global->eparticle_interval   = eparticle_interval;
  global->Hparticle_interval   = Hparticle_interval;
  global->quota_check_interval = quota_check_interval;
  global->quota_sec            = quota_sec;
  global->rtoggle              = 0;
  global->restart_step         = 0;

  global->b0  = b0;
  global->bg  = bg;
  global->v_A = v_A;
  global->c0 = c;

  global->topology_x = topology_x;
  global->topology_y = topology_y;
  global->topology_z = topology_z;

  // particle tracking
  global->particle_tracing      = particle_tracing;
  global->tracer_interval       = tracer_interval;
  global->tracer_file_interval  = nframes_per_tracer_file * tracer_interval;
  global->tracer_pass1_interval = tracer_pass1_interval;
  global->tracer_pass2_interval = tracer_pass2_interval;
  global->Ntracer = int(Ntracer);
  global->dump_traj_directly = dump_traj_directly;
  global->emf_at_tracer   = emf_at_tracer;
  global->hydro_at_tracer = hydro_at_tracer;
  global->ve_at_tracer = ve_at_tracer;
  global->num_tracer_fields_add = num_emf + num_hydro;

#ifdef OPEN_BOUNDARY_CONDITIONS
  global->edrive  = edrive;
  global->tdrive  = tdrive;
  global->left = left;
  global->right = right;
  global->top = top;
  global->bottom = bottom;
  global->nb  = nb_n0;
  global->rin[0]  = rin[0];
  global->rin[1]  = rin[1];
  global->rin[2]  = rin[2];
  global->rout[0] = rout[0];
  global->rout[1] = rout[1];
  global->rout[2] = rout[2];
  global->vth[0]  = sqrt(2)*vthe;
  global->vth[1]  = sqrt(2)*vthi;
  global->vthb[0] = sqrt(2)*vtheb;
  global->vthb[1] = sqrt(2)*vthib;
#ifdef HARRIS_SHEET
  global->q[0]  = weight_s;
  global->q[1]  = weight_s;
  global->uf[0] = udre;
  global->uf[1] = udri;
#else
  global->q[0]  = weight;
  global->q[1]  = weight;
  global->uf[0] = 0.0;
  global->uf[1] = 0.0;
#endif
  global->sort[0] = electron_sort_interval;
  global->sort[1] = ion_sort_interval;
  global->nfac = nfac;
  global->L_de = L;
#endif

  // particle dump
  global->stride_particle_dump = stride_particle_dump;

  /////////////////////////////////////////////////////////////////////////////
  // Setup the grid

  // Setup basic grid parameters
  grid->dx = hx;
  grid->dy = hy;
  grid->dz = hz;
  grid->dt = dt;
  grid->cvac = c;
  grid->eps0 = eps0;

#ifdef ENABLE_RADIATION_REACTION
  grid->re = re_de;
  grid->erho_rad  = erho_rad;
  grid->gamma_min = gamma_min;
#endif

  // Define the grid

  define_periodic_grid(  0, -0.5*Ly, -0.5*Lz,     // Low corner
                        Lx,  0.5*Ly,  0.5*Lz,     // High corner
                        nx, ny, nz,               // Resolution
                        topology_x, topology_y, topology_z); // Topology

  // ***** Set Field Boundary Conditions *****

#ifndef TURBULENCE_SIMULATION
  sim_log("Conducting fields on Z-boundaries");
  if ( iz==0 )            set_domain_field_bc( BOUNDARY(0,0,-1), pec_fields );
  if ( iz==topology_z-1 ) set_domain_field_bc( BOUNDARY( 0,0,1), pec_fields );

#ifdef OPEN_BOUNDARY_CONDITIONS
  sim_log("Absorb fields on X-boundaries");
  if ( ix==0 )            set_domain_field_bc( BOUNDARY(-1,0,0), absorb_fields );
  if ( ix==topology_x-1 ) set_domain_field_bc( BOUNDARY( 1,0,0), absorb_fields );

  // ***** Set Particle Boundary Conditions *****
  sim_log("Absorb particles on X & Z-boundaries");
  if ( iz==0 )            set_domain_particle_bc( BOUNDARY(0,0,-1), absorb_particles );
  if ( iz==topology_z-1 ) set_domain_particle_bc( BOUNDARY(0,0,1), absorb_particles );
  if ( ix==0 )            set_domain_particle_bc( BOUNDARY(-1,0,0), absorb_particles );
  if ( ix==topology_x-1 ) set_domain_particle_bc( BOUNDARY(1,0,0), absorb_particles );
#else
  sim_log("Reflective particles on Z-boundaries");
  if ( iz==0 )            set_domain_particle_bc( BOUNDARY(0,0,-1), reflect_particles );
  if ( iz==topology_z-1 ) set_domain_particle_bc( BOUNDARY(0,0,1), reflect_particles );
#endif
#endif // #ifndef TURBULENCE_SIMULATION

  /////////////////////////////////////////////////////////////////////////////
  // Setup materials

  sim_log("Setting up materials. ");

  define_material( "vacuum", 1 );
#ifdef OPEN_BOUNDARY_CONDITIONS
  material_t * resistive = define_material( "resistive",1,1,10);
#endif

  // Note: define_material defaults to isotropic materials with mu=1,sigma=0
  // Tensor electronic, magnetic and conductive materials are supported
  // though. See "shapes" for how to define them and assign them to regions.
  // Also, space is initially filled with the first material defined.

  /////////////////////////////////////////////////////////////////////////////
  // Finalize Field Advance

  sim_log("Finalizing Field Advance");

  define_field_array(NULL);

#ifdef OPEN_BOUNDARY_CONDITIONS
  // Define resistive layer surrounding boundary --> set thickness=0 to eliminate this feature
  double thickness = 1;
//#define resistive_layer (x < hx*thickness || x > Lx-hx*thickness || z <-Lz/2+hz*thickness  || z > Lz/2-hz*thickness )
#define resistive_layer (x < hx*thickness || x > Lx-hx*thickness )

  if (thickness > 0) set_region_material(resistive_layer, resistive, resistive);
#endif

  /////////////////////////////////////////////////////////////////////////////
  // Setup the species

  sim_log("Setting up species. ");
  double nparticles_per_proc = 13*Ne/nproc();
  double nmovers_per_proc    = 0.3*nparticles_per_proc;
  int ntracers_per_proc      = 10*nparticles_per_proc/particle_select;

  sim_log( "num_particles_per_proc = "<<nparticles_per_proc );
  sim_log( "num_tracer_particles_per_proc = "<<ntracers_per_proc );

#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
  int nparticles_part_efield_per_proc = nparticles_per_proc/particle_select_part_efield;
  double nmovers_part_efield_per_proc = 0.1*nparticles_part_efield_per_proc;
  sim_log( "nparticles_part_efield_per_proc = "<<nparticles_part_efield_per_proc );
#endif

#ifdef CHECK_PARTICLE_CROSSING_EGTB
  int nparticles_egtb_per_proc = nparticles_per_proc/particle_select_egtb;
  double nmovers_egtb_per_proc = 0.1*nparticles_egtb_per_proc;
  sim_log( "nparticles_egtb_per_proc = "<<nparticles_egtb_per_proc );
#endif

  double sort_method = 1;   //  0=in place and 1=out of place
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
  species_t *electronTop = define_species("electronTop",-ec, me, nparticles_per_proc,
          nmovers_per_proc, electron_sort_interval, sort_method);
  species_t *electronBot = define_species("electronBot",-ec, me, nparticles_per_proc,
          nmovers_per_proc, electron_sort_interval, sort_method);
  species_t *ionTop = define_species("ionTop", ec, mi, nparticles_per_proc,
          nmovers_per_proc, ion_sort_interval, sort_method);
  species_t *ionBot = define_species("ionBot", ec, mi, nparticles_per_proc,
          nmovers_per_proc, ion_sort_interval, sort_method);
  global->nsp = 4;
#else
  species_t *electron = define_species("electron",-ec, me, nparticles_per_proc,
          nmovers_per_proc, electron_sort_interval, sort_method);
  species_t *ion = define_species("ion", ec, mi, nparticles_per_proc,
          nmovers_per_proc, ion_sort_interval, sort_method);
  species_t *s, *s0;
  global->nsp = 2;
  global->nsp_efield = 0;
#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
  sim_log("Setting up species advanced with part of electric field");
  species_t *electron_wo_epara = define_species("electron_wo_epara",-ec, me,
      nparticles_part_efield_per_proc, nmovers_part_efield_per_proc, electron_sort_interval, sort_method);
  species_t *ion_wo_epara = define_species("ion_wo_epara", ec, mi,
      nparticles_part_efield_per_proc, nmovers_part_efield_per_proc, ion_sort_interval, sort_method);
  species_t *electron_wo_eparay = define_species("electron_wo_eparay",-ec, me,
      nparticles_part_efield_per_proc, nmovers_part_efield_per_proc, electron_sort_interval, sort_method);
  species_t *ion_wo_eparay = define_species("ion_wo_eparay", ec, mi,
      nparticles_part_efield_per_proc, nmovers_part_efield_per_proc, ion_sort_interval, sort_method);
  species_t *electron_wo_egtb = define_species("electron_wo_egtb",-ec, me,
      nparticles_part_efield_per_proc, nmovers_part_efield_per_proc, electron_sort_interval, sort_method);
  species_t *ion_wo_egtb = define_species("ion_wo_egtb", ec, mi,
      nparticles_part_efield_per_proc, nmovers_part_efield_per_proc, ion_sort_interval, sort_method);
  global->nsp_efield += 6;
#ifdef SPECIES_ADVANCED_WITHOUT_EPERP
  species_t *electron_wo_eperp = define_species("electron_wo_eperp",-ec, me,
      nparticles_part_efield_per_proc, nmovers_part_efield_per_proc, electron_sort_interval, sort_method);
  species_t *ion_wo_eperp = define_species("ion_wo_eperp", ec, mi,
      nparticles_part_efield_per_proc, nmovers_part_efield_per_proc, ion_sort_interval, sort_method);
  global->nsp_efield += 2;
#endif
  s = species_list;
  global->species_list_efield = species_list;
  for(int i = 1; i < global->nsp_efield; ++i) { // advance global->nsp_efield elements down species_list
    s0 = s;
    s = s->next;
    UNREGISTER_OBJECT( s0 );
  }
  // set species_list to point to the first physical particle species
  species_list = s->next;
  UNREGISTER_OBJECT( s );
  s->next = NULL ;
#endif
#ifdef CHECK_PARTICLE_CROSSING_EGTB
  sim_log("Setting up species that will check whether particles pass E>B regions");
  species_t *electron_egtb = define_species("electron_egtb",-ec, me,
      nparticles_egtb_per_proc, nmovers_egtb_per_proc, electron_sort_interval, sort_method);
  species_t *ion_egtb = define_species("ion_egtb", ec, mi,
      nparticles_egtb_per_proc, nmovers_egtb_per_proc, ion_sort_interval, sort_method);
  s = species_list;
  global->species_list_egtb = species_list;
  for(int i = 1; i < 2 ; ++i) { // advance 2 elements down species_list
    s0 = s;
    s = s->next;
    UNREGISTER_OBJECT( s0 );
  }
  // set species_list to point to the first physical particle species
  species_list = s->next;
  UNREGISTER_OBJECT( s );
  s->next = NULL ;
#endif
#endif // #ifdef TURBULENCE_MIXING_DIAGNOSTICS

  // particle tracer species
  sim_log("Setting up tracer electrons.");
  sort_method = 0;   //  0=in place and 1=out of place
  species_t *e_tracer = define_species("electron_tracer",-ec, me, ntracers_per_proc,
          nmovers_per_proc, electron_sort_interval, sort_method);
  sim_log("Setting up tracer ions.");
  species_t *H_tracer = define_species("H_tracer", ec, mi, ntracers_per_proc,
          nmovers_per_proc, ion_sort_interval, sort_method);

  hijack_tracers(2);

  /////////////////////////////////////////////////////////////////////////////
  // Log diagnostic information about this simulation

  sim_log( "***********************************************" );
  sim_log("* Topology:                       " << topology_x
    << " " << topology_y << " " << topology_z);
#ifdef RELATIVISTIC_PLASMA
  sim_log ( "sigma   = " << sigma );
#endif
  sim_log ( "L_di   = " << L_di );
  sim_log ( "L/de   = " << L/de);
  sim_log ( "rhoi/L = " << rhoi_L );
  sim_log ( "Ti/Te = " << Ti_Te ) ;
  sim_log ( "Tbe/Te = " << Tbe_Te ) ;
  sim_log ( "Tbi/Ti = " << Tbi_Ti ) ;
  sim_log ( "Te = " << Te ) ;
  sim_log ( "Ti = " << Ti ) ;
  sim_log ( "nb/n0 = " << nb_n0 ) ;
  sim_log ( "wpe/wce = " << wpe_wce );
#ifdef RELATIVISTIC_PLASMA
  sim_log ( "factor   = " << factor );
#endif
  sim_log ( "mi/me = " << mi_me );
  sim_log ( "theta = " << theta );
  sim_log ( "taui = " << taui );
  sim_log ( "num_step = " << num_step );
  sim_log ( "Lx/di = " << Lx/di );
  sim_log ( "Lx/de = " << Lx/de );
  sim_log ( "Ly/di = " << Ly/di );
  sim_log ( "Ly/de = " << Ly/de );
  sim_log ( "Lz/di = " << Lz/di );
  sim_log ( "Lz/de = " << Lz/de );
  sim_log ( "nx = " << nx );
  sim_log ( "ny = " << ny );
  sim_log ( "nz = " << nz );
  sim_log ( "courant = " << c*dt/dg );
  sim_log ( "nproc = " << float(nproc())  );
  sim_log ( "nppc = " << nppc );
  sim_log ( "b0 = " << b0 );
  sim_log ( "v_A (based on nb) = " << v_A );
  sim_log ( "di = " << di );
  sim_log ( "Ne = " << Ne );
  sim_log ( "total # of particles = " << 2*Ne );
#ifdef HARRIS_SHEET
  sim_log ( "Npe_sheet = " << Npe_sheet );
  sim_log ( "Npe_back = " << Npe_back );
  sim_log ( "Npe = " << Npe );
  sim_log ( "Ne_sheet = " << Ne_sheet );
  sim_log ( "Ne_back = " << Ne_back );
  sim_log ( "weight_s = " << weight_s );
  sim_log ( "weight_b = " << weight_b );
  sim_log ( "qi_s = " << qi_s );
  sim_log ( "qi_s = " << qi_s );
  sim_log ( "qi_b = " << qi_b );
  sim_log ( "qe_s = " << qe_s );
  sim_log ( "qe_b = " << qe_b );
  sim_log ( "vdri = " << vdri );
  sim_log ( "vdre = " << vdre );
  sim_log ( "gdri = " << gdri );
  sim_log ( "gdre = " << gdre );
  sim_log ( "udri = " << udri );
  sim_log ( "udre = " << udre );
  sim_log ( "tanhf = " << tanhf );
#else
  sim_log ( "qi = " << qi );
  sim_log ( "qe = " << qe );
#endif
  sim_log ( "nfac = " << nfac );
  sim_log ( "dt*wpe = " << wpe*dt );
  sim_log ( "dt*wce = " << wce*dt );
  sim_log ( "dt*wci = " << wci*dt );
  sim_log ( "energies_interval = " << energies_interval );
  sim_log ( "dx/de = " << Lx/(de*nx) );
  sim_log ( "dy/de = " << Ly/(de*ny) );
  sim_log ( "dz/de = " << Lz/(de*nz) );
  sim_log ( "dx/rhoe = " << (Lx/nx)/(vthe/wce)  );
  sim_log ( "L/debye = " << L/(vthe/wpe)  );
  sim_log ( "dx/rhoi = " << (Lx/nx)/(vthi/wci) );
  sim_log ( "dx/rhoe = " << (Lx/nx)/(vthe/wce) );
  sim_log ( "dx/debye = " << (Lx/nx)/(vthe/wpe)  );
  sim_log ( "n0 = " << n0 );
  sim_log ( "vthi/c = " << vthi/c );
  sim_log ( "vthe/c = " << vthe/c );
  sim_log ( "vthib/c = " << vthib/c );
  sim_log ( "vtheb/c = " << vtheb/c );
  sim_log ( "restart_interval = "      << restart_interval );
  sim_log ( "spectrum_interval = "     << spectrum_interval );
  sim_log ( "fields_interval = "       << fields_interval );
  sim_log ( "ehydro_interval = "       << ehydro_interval );
  sim_log ( "Hhydro_interval = "       << Hhydro_interval );
  sim_log ( "eparticle_interval = "    << eparticle_interval );
  sim_log ( "Hparticle_interval = "    << Hparticle_interval );
  sim_log ( "quota_check_interval = "  << quota_check_interval );
  sim_log ( "particle_tracing = "      << particle_tracing );
  sim_log ( "tracer_interval = "       << tracer_interval );
  sim_log ( "tracer_file_interval = "  << global->tracer_file_interval );
  sim_log ( "tracer_pass1_interval = " << tracer_pass1_interval );
  sim_log ( "tracer_pass2_interval = " << tracer_pass2_interval );
  sim_log ( "Ntracer = "               << global->Ntracer );
  sim_log ( "emf_at_tracer = "         << emf_at_tracer );
  sim_log ( "hydro_at_tracer = "       << hydro_at_tracer );
  sim_log ( "ve_at_tracer = "          << ve_at_tracer );
  sim_log ( "dump_traj_directly = "    << dump_traj_directly );
  sim_log ( "num_tracer_fields_add = " << global->num_tracer_fields_add );
  sim_log ( "emax_band = " << emax_band );
  sim_log ( "emin_band = " << emin_band );
  sim_log ( "nbands = " << nbands );
  sim_log ( "emax_spect = " << emax_spect );
  sim_log ( "emin_spect = " << emin_spect );
  sim_log ( "nbins = " << nbins );
  sim_log ( "nx_zone = " << nx_zone );
  sim_log ( "ny_zone = " << ny_zone );
  sim_log ( "nz_zone = " << nz_zone );
  sim_log ( "stride_particle_dump = " << stride_particle_dump );

  /////////////////////////////////////////////////////////////////////////////
  // Dump simulation information to file "info"
  if (rank() == 0 ) {

    FileIO fp_info;

    if ( ! (fp_info.open("info", io_write)==ok) ) ERROR(("Cannot open file."));

    fp_info.print( "           ***** Simulation parameters ***** \n");
#ifdef RELATIVISTIC_PLASMA
    fp_info.print( " sigma = %e\n", sigma);
#endif
    fp_info.print( " L/di = %e\n", L_di);
    fp_info.print( " L/de = %e\n", L/de);
    fp_info.print( " rhoi/L = %e\n", rhoi_L);
    fp_info.print( " Ti/Te = %e\n", Ti_Te );
    fp_info.print( " Tbe/Te = %e\n", Tbe_Te ) ;
    fp_info.print( " Tbi/Ti = %e\n", Tbi_Ti ) ;
    fp_info.print( " Te = %e\n", Te ) ;
    fp_info.print( " Ti = %e\n", Ti ) ;
    fp_info.print( " nb/n0 = %e\n", nb_n0 );
    fp_info.print( " wpe/wce = %e\n", wpe_wce );
#ifdef RELATIVISTIC_PLASMA
    fp_info.print( " factor = %e\n", factor );
#endif
    fp_info.print( " mi/me = %e\n", mi_me );
    fp_info.print( " theta = %e\n", theta );
    fp_info.print( " taui = %e\n", taui );
    fp_info.print( " num_step = %i\n", num_step );
    fp_info.print( " Lx/de = %e\n", Lx/de );
    fp_info.print( " Ly/de = %e\n", Ly/de );
    fp_info.print( " Lz/de = %e\n", Lz/de );
    fp_info.print( " Lx/di = %e\n", Lx/di );
    fp_info.print( " Ly/di = %e\n", Ly/di );
    fp_info.print( " Lz/di = %e\n", Lz/di );
    fp_info.print( " nx = %e\n", nx );
    fp_info.print( " ny = %e\n", ny );
    fp_info.print( " nz = %e\n", nz );
    fp_info.print( " courant = %e\n", c*dt/dg );
    fp_info.print( " nproc = %i\n", nproc() );
    fp_info.print( " nppc = %e\n", nppc );
    fp_info.print( " b0 = %e\n", b0 );
    fp_info.print( " v_A (based on nb) = %e\n", v_A );
    fp_info.print( " di = %e\n", di );
    fp_info.print( " Ne = %e\n", Ne );
    fp_info.print( " total # of particles = %e\n", 2*Ne );
#ifdef HARRIS_SHEET
    fp_info.print( " Npe_sheet = %e\n", Npe_sheet );
    fp_info.print( " Npe_back = %e\n", Npe_back );
    fp_info.print( " Npe = %e\n", Npe );
    fp_info.print( " Ne_sheet = %e\n", Ne_sheet );
    fp_info.print( " Ne_back = %e\n", Ne_back );
    fp_info.print( " weight_s = %e\n", weight_s );
    fp_info.print( " weight_b = %e\n", weight_b );
    fp_info.print( " qi_s = %e\n", qi_s );
    fp_info.print( " qi_b = %e\n", qi_b );
    fp_info.print( " qe_s = %e\n", qe_s );
    fp_info.print( " qe_b = %e\n", qe_b );
    fp_info.print( " vdri = %e\n", vdri );
    fp_info.print( " vdre = %e\n", vdre );
    fp_info.print( " gdri = %e\n", gdri );
    fp_info.print( " gdre = %e\n", gdre );
    fp_info.print( " udri = %e\n", udri );
    fp_info.print( " udre = %e\n", udre );
    fp_info.print( " tanhf = %e\n", tanhf );
#else
    fp_info.print( " qi = %e\n", qi );
    fp_info.print( " qe = %e\n", qe );
#endif
    fp_info.print( " nfac = %e\n", nfac );
    fp_info.print( " dt*wpe = %e\n", wpe*dt );
    fp_info.print( " dt*wce = %e\n", wce*dt );
    fp_info.print( " dt*wci = %e\n", wci*dt );
    fp_info.print( " energies_interval = %i\n", energies_interval);
    fp_info.print( " dx/de = %e\n", Lx/(de*nx) );
    fp_info.print( " dy/de = %e\n", Ly/(de*ny) );
    fp_info.print( " dz/de = %e\n", Lz/(de*nz) );
    fp_info.print( " dx/rhoe = %e\n", (Lx/nx)/(vthe/wce)  );
    fp_info.print( " L/debye = %e\n", L/(vthe/wpe) );
    fp_info.print( " dx/rhoi = %e\n", (Lx/nx)/(vthi/wci) );
    fp_info.print( " dx/rhoe = %e\n", (Lx/nx)/(vthe/wce) );
    fp_info.print( " dx/debye = %e\n", (Lx/nx)/(vthe/wpe) );
    fp_info.print( " n0 = %e\n", n0 );
    fp_info.print( " vthi/c = %e\n", vthi/c );
    fp_info.print( " vthe/c = %e\n", vthe/c );
    fp_info.print( " vthib/c = %e\n", vthib/c );
    fp_info.print( " vtheb/c = %e\n", vtheb/c );
    fp_info.print( " restart_interval = %i\n", restart_interval );
    fp_info.print( " spectrum_interval = %i\n", spectrum_interval );
    fp_info.print( " fields_interval = %i\n", fields_interval );
    fp_info.print( " ehydro_interval = %i\n", ehydro_interval );
    fp_info.print( " Hhydro_interval = %i\n", Hhydro_interval );
    fp_info.print( " eparticle_interval = %i\n", eparticle_interval );
    fp_info.print( " Hparticle_interval = %i\n", Hparticle_interval );
    fp_info.print( " quota_check_interval = %i\n", quota_check_interval );
    fp_info.print( " particle_tracing = %i\n", particle_tracing );
    fp_info.print( " tracer_interval = %i\n", tracer_interval );
    fp_info.print( " tracer_file_interval = %i\n", global->tracer_file_interval );
    fp_info.print( " tracer_pass1_interval = %i\n", tracer_pass1_interval );
    fp_info.print( " tracer_pass2_interval = %i\n", tracer_pass2_interval );
    fp_info.print( " Ntracer = %i\n", global->Ntracer );
    fp_info.print( " emf_at_tracer = %i\n", emf_at_tracer );
    fp_info.print( " hydro_at_tracer = %i\n", hydro_at_tracer );
    fp_info.print( " ve_at_tracer = %i\n", ve_at_tracer );
    fp_info.print( " dump_traj_directly = %i\n", dump_traj_directly );
    fp_info.print( " num_tracer_fields_add = %i\n", global->num_tracer_fields_add );
    fp_info.print( " emax_band = %e\n", emax_band );
    fp_info.print( " emin_band = %e\n", emin_band );
    fp_info.print( " nbands = %i\n", nbands );
    fp_info.print( " emax_spect = %e\n", emax_spect );
    fp_info.print( " emin_spect = %e\n", emin_spect );
    fp_info.print( " nbins = %i\n", nbins );
    fp_info.print( " nx_zone = %i\n", nx_zone );
    fp_info.print( " ny_zone = %i\n", ny_zone );
    fp_info.print( " nz_zone = %i\n", nz_zone );
    fp_info.print( " stride_particle_dump = %i\n", stride_particle_dump );
    fp_info.print( " ***************************\n");
    fp_info.close();


    // for the parallelized translate.f90 written by Vadim
    // write binary info file

    if ( ! (fp_info.open("info.bin", io_write)==ok) )
      ERROR(("Cannot open file."));

    fp_info.write(&topology_x, 1 );
    fp_info.write(&topology_y, 1 );
    fp_info.write(&topology_z, 1 );

    fp_info.write(&Lx, 1 );
    fp_info.write(&Ly, 1 );
    fp_info.write(&Lz, 1 );

    fp_info.write(&nx, 1 );
    fp_info.write(&ny, 1 );
    fp_info.write(&nz, 1 );

    fp_info.write(&dt, 1 );

    fp_info.write(&mi_me, 1 );
    fp_info.write(&wpe_wce, 1 );
    fp_info.write(&vthe, 1 );
    fp_info.write(&vthi, 1 );

    fp_info.close();
  }

  /////////////////////////////////////////////////////////////////////////////
  // Load fields


  // Define some function to load profiles

  sim_log( "Loading fields" );


#ifdef TURBULENCE_SIMULATION
  // the setup follows Comisso & Sironi 19: 10.3847/1538-4357/ab4c33
  set_region_field( everywhere, 0, 0, 0, 0, 0, b0);
  seed_entropy( rank() );
  const double _x0 = grid->x0, _y0 = grid->y0, _z0 = grid->z0;
  const double _dx = grid->dx, _dy = grid->dy, _dz = grid->dz;
  const double _c  = grid->cvac;
  const int    _nx = grid->nx, _ny = grid->ny, _nz = grid->nz;
  double *phi_array;

  // set the phases of the modes
  int nphases;
  if (nz > 1) { // 3D
    nphases = nmodes * nmodes * nmodes_para * 4;
  } else {
    nphases = nmodes * nmodes * 2;
  }
  phi_array = new double[nphases];
  if (rank() == 0) {
    for (int m=0; m<nmodes; m++) {
      for (int n=0; n<nmodes; n++) {
        if (nz > 1) { // 3D
          for (int l=0; l<nmodes_para; l++) {
            int imode = 4 * ((m * nmodes + n) * nmodes_para + l); 
            for (int im=0; im<4; im++) {
              phi_array[imode+im] = uniform(rng(0), 0, 2 * M_PI);
            }
          }
        } else { // 2D
          int imode = 2 * (m * nmodes + n); 
          phi_array[imode] = uniform(rng(0), 0, 2 * M_PI);
          phi_array[imode+1] = uniform(rng(0), 0, 2 * M_PI);
        }
      }
    }
  }
  mp_barrier(); // Just to be safe
  MPI_Bcast(phi_array, nphases, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for( int _k=0; _k<_nz+2; _k++ ) {
    const double _ze = _z0 + _dz*_k;
    const double _zc = _z0 + _dz*(_k-0.5);
    for( int _j=0; _j<_ny+2; _j++ ) {
      const double _ye = _y0 + _dy*_j;
      const double _yc = _y0 + _dy*(_j-0.5);
      field_t *_f = &field(0,_j,_k);
      for( int _i=0; _i<_nx+2; _i++ ) {
        const double _xe = _x0 + _dx*_i;
        const double _xc = _x0 + _dx*(_i-0.5);
        for (int m=1; m<nmodes+1; m++) {
          double kx = 2 * M_PI * m / Lx;
          for (int n=1; n<nmodes+1; n++) {
            double ky = 2 * M_PI * n / Ly;
            double kx_kxy = kx / sqrt(kx*kx + ky*ky);
            double ky_kxy = ky / sqrt(kx*kx + ky*ky);
            double _phi;
            if (nz == 1) { // 2D simulation
              // factor 2 means that forward and backward modes have the same amplitude but opposite phases
              int imode = 2 * ((m-1)*nmodes + n - 1);
              _phi = phi_array[imode];
              _f->cbx -= 2 * _c * db_amp * ky_kxy * sin(kx*_xe + ky*_yc + _phi);
              _f->cby += 2 * _c * db_amp * kx_kxy * sin(kx*_xc + ky*_ye + _phi);
              _phi = phi_array[imode+1];
              _f->cbx -= 2 * _c * db_amp * ky_kxy * sin(-kx*_xe + ky*_yc + _phi);
              _f->cby -= 2 * _c * db_amp * kx_kxy * sin(-kx*_xc + ky*_ye + _phi);
            } else { // 3D simulation
              for (int l=1; l<nmodes_para+1; l++) {
                double kz = 2 * M_PI * l / Lz;
                int imode = 4 * (((m - 1) * nmodes + n - 1) * nmodes_para + l - 1); 
                _phi = phi_array[imode];
                _f->cbx -= 2 * _c * db_amp * ky_kxy * sin(kx*_xe + ky*_yc + kz*_zc + _phi);
                _f->cby += 2 * _c * db_amp * kx_kxy * sin(kx*_xc + ky*_ye + kz*_zc + _phi);
                _phi = phi_array[imode+1];
                _f->cbx -= 2 * _c * db_amp * ky_kxy * sin(-kx*_xe + ky*_yc + kz*_zc + _phi);
                _f->cby -= 2 * _c * db_amp * kx_kxy * sin(-kx*_xc + ky*_ye + kz*_zc + _phi);
                _phi = phi_array[imode+2];
                _f->cbx += 2 * _c * db_amp * ky_kxy * sin(kx*_xe - ky*_yc + kz*_zc + _phi);
                _f->cby += 2 * _c * db_amp * kx_kxy * sin(kx*_xc - ky*_ye + kz*_zc + _phi);
                _phi = phi_array[imode+3];
                _f->cbx += 2 * _c * db_amp * ky_kxy * sin(-kx*_xe - ky*_yc + kz*_zc + _phi);
                _f->cby -= 2 * _c * db_amp * kx_kxy * sin(-kx*_xc - ky*_ye + kz*_zc + _phi);
              } // l
            }
          } // n
        } // m
        _f++;
      } // _i
    } // _j
  } // _k
  delete [] phi_array;
#else //#ifdef TURBULENCE_SIMULATION

#ifdef LOCAL_PERTURBATION
#define R2 (((x-0.5*Lx)*(x-0.5*Lx) + z*z)/(L*L))
#define PERT  (0.2*tanh(R2)/cosh(R2))
#define DBX (-z*PERT)
#define DBZ ((x-0.5*Lx)*PERT)
#else
  // modified perturbation
#define DBX (dbx*cos(2.0*pi*(x-0.5*Lx)/Lpert)*sin(pi*z/Lz))
#define DBZ (dbz*cos(pi*z/Lz)*sin(2.0*pi*(x-0.5*Lx)/Lpert))
#endif

#define BX (b0*tanh(z/L))

#ifdef HARRIS_SHEET
#define BY (bg*b0)
#else
#define BY sqrt(b0*b0 + bg*bg*b0*b0 - BX*BX)
#endif // ifdef HARRIS_SHEET

  // Define initial perturbations to B to drive turbulence
  // Fundamental mode numbers
  double kx = 2.0*pi/Lx;
  double ky = 2.0*pi/Ly;
  double kz = pi/Lz;

  // l --> mode number in X
  // m --> mode number in Y
  // n --> mode number in Z

  //define DBX(a,n,m,phi) -a*((n*kz)/(m*kx))*sin(n*kz*z)*cos(m*kx*x + phi)

#define DBY_WAVE(a,l,n,phi) (a*cos(l*kx*x+phi)*cos(n*kz*z))
#define DBZ_WAVE(a,l,m,phi) (a*cos(l*kx*x)*sin(m*ky*y+phi))

#define BYWAVE DBY_WAVE(0.02,2,1,0.0)  + DBY_WAVE(0.03,3,2,0.2)  + DBY_WAVE(0.04,4,1,-0.5) + DBY_WAVE(0.05,5,3,0.6) + DBY_WAVE(0.06,6,5,-0.8) + DBY_WAVE(0.07,7,5,0.8)
#define BZWAVE DBZ_WAVE(0.025,2,1,0.5) + DBZ_WAVE(0.02,3,2,-0.2) + DBZ_WAVE(0.045,4,3,-0.3) + DBZ_WAVE(0.04,5,4,0.3) + DBZ_WAVE(0.05,6,5,0.8) + DBZ_WAVE(0.075,7,6,0.8)

  // Note: everywhere is a region that encompasses the entire simulation
  // In general, regions are specified as logical equations (i.e. x>0 && x+y<2)

#ifdef WAVE_PERTURBATION
  set_region_field(everywhere, 0, 0, 0, (BX+DBX)*cs+(BY+BYWAVE)*sn,
                   -(BX+DBX)*sn+(BY+BYWAVE)*cs, DBZ + BZWAVE);
#else
  set_region_field( everywhere, 0, 0, 0, (BX+DBX)*cs+BY*sn, -(BX+DBX)*sn+BY*cs, DBZ);
#endif

#endif //#ifdef TURBULENCE_SIMULATION

  /////////////////////////////////////////////////////////////////////////////
  // LOAD PARTICLES

  // particle tracking parameters
  int itracer    = 0;   // tracer index
  int iparticle  = 0;   // particle index

  sim_log( "Loading particles" );

  seed_entropy( rank() );  //Generators desynchronized
  double xmin = grid->x0 , xmax = grid->x0+(grid->dx)*(grid->nx);
  double ymin = grid->y0 , ymax = grid->y0+(grid->dy)*(grid->ny);
  double zmin = grid->z0 , zmax = grid->z0+(grid->dz)*(grid->nz);

#ifdef TURBULENCE_SIMULATION

  repeat ( Ne/nproc() ) {
    double x, y, z, ux, uy, uz, upa1, upe1, uz1, gu1;

    x = uniform(rng(0), xmin, xmax);
    y = uniform(rng(0), ymin, ymax);
    z = uniform(rng(0), zmin, zmax);

    // inject_particles() will return an error for particles no on this
    // node and will not inject particle locally

#ifdef RELATIVISTIC_PLASMA
    double fs0, f, fs, u2, gg, ue_range, ui_range;
    fs0 = 1.01;
    ue_range = factor*sqrt((Te/me) * (Te/me + 2));
    f = 1.0;
    fs = 0.0;
    do {
      upa1 = uniform(rng(0), -ue_range, ue_range);
      upe1 = uniform(rng(0), -ue_range, ue_range);
      uz1 = uniform(rng(0), -ue_range, ue_range);
      u2 = upa1*upa1 + upe1*upe1 + uz1*uz1;
      gg = sqrt(1.0 + u2);
      /* fs = exp(-u2/(1.0+gg)*me/Te); */
      fs = exp(-(gg-1)*me/Te);
      f  = uniform(rng(0), 0, 1) * fs0;
    } while(f > fs);
#else
    upa1 = normal(rng(0), 0, vthe);
    upe1 = normal(rng(0), 0, vthe);
    uz1 = normal(rng(0), 0, vthe);
#endif

    ux = upa1;
    uy = upe1;
    uz = uz1;

    inject_particle(electron, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
    ++iparticle;
#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
    if (iparticle%particle_select_part_efield == 0) {
      inject_particle(electron_wo_epara, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
      inject_particle(electron_wo_eparay, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
      inject_particle(electron_wo_egtb, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
#ifdef SPECIES_ADVANCED_WITHOUT_EPERP
      inject_particle(electron_wo_eperp, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
#endif
    }
#endif
#ifdef CHECK_PARTICLE_CROSSING_EGTB
    if (iparticle%particle_select_egtb == 0) {
      inject_particle(electron_egtb, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
    }
#endif

#ifdef RELATIVISTIC_PLASMA
    ui_range = factor*sqrt((Ti/mi) * (Ti/mi + 2));
    do {
      upa1 = uniform(rng(0), -ui_range, ui_range);
      upe1 = uniform(rng(0), -ui_range, ui_range);
      uz1 = uniform(rng(0), -ui_range, ui_range);
      u2 = upa1*upa1 + upe1*upe1 + uz1*uz1;
      gg = sqrt(1.0 + u2);
      /* fs = exp(-u2/(1.0+gg)*mi/Ti); */
      fs = exp(-(gg-1)*mi/Ti);
      f  = uniform(rng(0), 0, 1) * fs0;
    } while(f > fs);
#else
    upa1 = normal(rng(0), 0, vthi);
    upe1 = normal(rng(0), 0, vthi);
    uz1 = normal(rng(0), 0, vthi);
#endif

    ux = upa1;
    uy = upe1;
    uz = uz1;

    inject_particle(ion, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
    if (iparticle%particle_select_part_efield == 0) {
      inject_particle(ion_wo_epara, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
      inject_particle(ion_wo_eparay, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
      inject_particle(ion_wo_egtb, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
#ifdef SPECIES_ADVANCED_WITHOUT_EPERP
      inject_particle(ion_wo_eperp, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
#endif
    }
#endif
#ifdef CHECK_PARTICLE_CROSSING_EGTB
    if (iparticle%particle_select_egtb == 0) {
      inject_particle(ion_egtb, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
    }
#endif

    if (particle_tracing == 1) { // only tag particles in the 1st pass
      if (iparticle%particle_select == 0) {
        itracer++;
        /* int tag = ((((int)rank())<<19) | (itracer & 0x7ffff)); // 13 bits (8192) for rank and 19 bits (~520k) */
        int tag = ((((int)rank())<<13) | (itracer & 0x1fff)); // 19 bits (520k) for rank and 13 bits (8192)
        tag_tracer( (electron->p + electron->np-1), e_tracer, tag );
        tag_tracer( (ion->p      + ion->np-1),      H_tracer, tag );
      }
    }
  }

#else // ifdef TURBULENCE_SIMULATION

#ifdef HARRIS_SHEET
  // Harris current sheet
  sim_log( "-> Main Harris Sheet" );

  repeat ( Ne_sheet/nproc() ) {
    double x, y, z, ux, uy, uz, d0 ;

    do {
      z = L*atanh(uniform(rng(0),-1,1)*tanhf);
    } while( z<= zmin || z>=zmax );
    x = uniform(rng(0), xmin, xmax);
    y = uniform(rng(0), ymin, ymax);

    // inject_particles() will return an error for particles no on this
    // node and will not inject particle locally

    ux = normal(rng(0), 0, vthe);
    uy = normal(rng(0), 0, vthe);
    uz = normal(rng(0), 0, vthe);
    d0 = gdre*uy + sqrt(ux*ux + uy*uy + uz*uz + 1)*udre;
    uy = d0*cs - ux*sn;
    ux = d0*sn + ux*cs;

#ifdef TURBULENCE_MIXING_DIAGNOSTICS
    if (z > 0) {
      inject_particle(electronTop, x, y, z, ux, uy, uz, weight_s, 0, 0 );
    } else {
      inject_particle(electronBot, x, y, z, ux, uy, uz, weight_s, 0, 0 );
    }
    ++iparticle;
#else
    inject_particle(electron, x, y, z, ux, uy, uz, weight_s, 0, 0 );
    ++iparticle;
#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
    if (iparticle%particle_select_part_efield == 0) {
      inject_particle(electron_wo_epara, x, y, z, ux, uy, uz, weight_s, 0, 0 );
      inject_particle(electron_wo_eparay, x, y, z, ux, uy, uz, weight_s, 0, 0 );
      inject_particle(electron_wo_egtb, x, y, z, ux, uy, uz, weight_s, 0, 0 );
#ifdef SPECIES_ADVANCED_WITHOUT_EPERP
      inject_particle(electron_wo_eperp, x, y, z, ux, uy, uz, weight_s, 0, 0 );
#endif
    }
#endif
#ifdef CHECK_PARTICLE_CROSSING_EGTB
    if (iparticle%particle_select_egtb == 0) {
      inject_particle(electron_egtb, x, y, z, ux, uy, uz, weight_s, 0, 0 );
    }
#endif
#endif

    ux = normal(rng(0), 0, vthi);
    uy = normal(rng(0), 0, vthi);
    uz = normal(rng(0), 0, vthi);
    d0 = gdri*uy + sqrt(ux*ux + uy*uy + uz*uz + 1)*udri;
    uy = d0*cs - ux*sn;
    ux = d0*sn + ux*cs;

#ifdef TURBULENCE_MIXING_DIAGNOSTICS
    if (z > 0) {
      inject_particle(ionTop, x, y, z, ux, uy, uz, weight_s, 0, 0 );
    } else {
      inject_particle(ionBot, x, y, z, ux, uy, uz, weight_s, 0, 0 );
    }
#else
    inject_particle(ion, x, y, z, ux, uy, uz, weight_s, 0, 0 );
#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
    if (iparticle%particle_select_part_efield == 0) {
      inject_particle(ion_wo_epara, x, y, z, ux, uy, uz, weight_s, 0, 0 );
      inject_particle(ion_wo_eparay, x, y, z, ux, uy, uz, weight_s, 0, 0 );
      inject_particle(ion_wo_egtb, x, y, z, ux, uy, uz, weight_s, 0, 0 );
#ifdef SPECIES_ADVANCED_WITHOUT_EPERP
      inject_particle(ion_wo_eperp, x, y, z, ux, uy, uz, weight_s, 0, 0 );
#endif
    }
#endif
#ifdef CHECK_PARTICLE_CROSSING_EGTB
    if (iparticle%particle_select_egtb == 0) {
      inject_particle(ion_egtb, x, y, z, ux, uy, uz, weight_s, 0, 0 );
    }
#endif
#endif

    // tag and inject tracers
    if (particle_tracing == 1) { // only tag particles in the 1st pass
      if (iparticle%particle_select == 0) {
        itracer++;
        /* int tag = ((((int)rank())<<19) | (itracer & 0x7ffff)); // 13 bits (8192) for rank and 19 bits (~520k) */
        int tag = ((((int)rank())<<13) | (itracer & 0x1fff)); // 19 bits (520k) for rank and 13 bits (8192)
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
        if (z > 0) {
          tag_tracer( (electronTop->p + electronTop->np-1), e_tracer, tag );
          tag_tracer( (ionTop->p      + ionTop->np-1),      H_tracer, tag );
        } else {
          tag_tracer( (electronBot->p + electronBot->np-1), e_tracer, tag );
          tag_tracer( (ionBot->p      + ionBot->np-1),      H_tracer, tag );
        }
#else
        tag_tracer( (electron->p + electron->np-1), e_tracer, tag );
        tag_tracer( (ion->p      + ion->np-1),      H_tracer, tag );
#endif
      }
    }
  }

  sim_log( "-> Background Population" );

  repeat ( Ne_back/nproc() ) {
    double x, y, z, ux, uy, uz;

    x = uniform(rng(0), xmin, xmax);
    y = uniform(rng(0), ymin, ymax);
    z = uniform(rng(0), zmin, zmax);

    ux = normal(rng(0), 0, vtheb);
    uy = normal(rng(0), 0, vtheb);
    uz = normal(rng(0), 0, vtheb);

#ifdef TURBULENCE_MIXING_DIAGNOSTICS
    if (z > 0) {
      inject_particle(electronTop, x, y, z, ux, uy, uz, weight_b, 0, 0 );
    } else {
      inject_particle(electronBot, x, y, z, ux, uy, uz, weight_b, 0, 0 );
    }
    ++iparticle;
#else
    inject_particle(electron, x, y, z, ux, uy, uz, weight_b, 0, 0 );
    ++iparticle;
#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
    if (iparticle%particle_select_part_efield == 0) {
      inject_particle(electron_wo_epara, x, y, z, ux, uy, uz, weight_b, 0, 0 );
      inject_particle(electron_wo_eparay, x, y, z, ux, uy, uz, weight_b, 0, 0 );
      inject_particle(electron_wo_egtb, x, y, z, ux, uy, uz, weight_b, 0, 0 );
#ifdef SPECIES_ADVANCED_WITHOUT_EPERP
      inject_particle(electron_wo_eperp, x, y, z, ux, uy, uz, weight_b, 0, 0 );
#endif
    }
#endif
#ifdef CHECK_PARTICLE_CROSSING_EGTB
    if (iparticle%particle_select_egtb == 0) {
      inject_particle(electron_egtb, x, y, z, ux, uy, uz, weight_b, 0, 0 );
    }
#endif
#endif

    ux = normal(rng(0), 0, vthib);
    uy = normal(rng(0), 0, vthib);
    uz = normal(rng(0), 0, vthib);

#ifdef TURBULENCE_MIXING_DIAGNOSTICS
    if (z > 0) {
      inject_particle(ionTop, x, y, z, ux, uy, uz, weight_b, 0, 0 );
    } else {
      inject_particle(ionBot, x, y, z, ux, uy, uz, weight_b, 0, 0 );
    }
#else
    inject_particle(ion, x, y, z, ux, uy, uz, weight_b, 0, 0 );
#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
    if (iparticle%particle_select_part_efield == 0) {
      inject_particle(ion_wo_epara, x, y, z, ux, uy, uz, weight_b, 0, 0 );
      inject_particle(ion_wo_eparay, x, y, z, ux, uy, uz, weight_b, 0, 0 );
      inject_particle(ion_wo_egtb, x, y, z, ux, uy, uz, weight_b, 0, 0 );
#ifdef SPECIES_ADVANCED_WITHOUT_EPERP
      inject_particle(ion_wo_eperp, x, y, z, ux, uy, uz, weight_b, 0, 0 );
#endif
    }
#endif
#ifdef CHECK_PARTICLE_CROSSING_EGTB
    if (iparticle%particle_select_egtb == 0) {
      inject_particle(ion_egtb, x, y, z, ux, uy, uz, weight_b, 0, 0 );
    }
#endif
#endif

    // tag and inject tracers
    if (particle_tracing == 1) { // only tag particles in the 1st pass
      if (iparticle%particle_select == 0) {
        itracer++;
        /* int tag = ((((int)rank())<<19) | (itracer & 0x7ffff)); // 13 bits (8192) for rank and 19 bits (~520k) */
        int tag = ((((int)rank())<<13) | (itracer & 0x1fff)); // 19 bits (520k) for rank and 13 bits (8192)
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
        if (z > 0) {
          tag_tracer( (electronTop->p + electronTop->np-1), e_tracer, tag );
          tag_tracer( (ionTop->p      + ionTop->np-1),      H_tracer, tag );
        } else {
          tag_tracer( (electronBot->p + electronBot->np-1), e_tracer, tag );
          tag_tracer( (ionBot->p      + ionBot->np-1),      H_tracer, tag );
        }
#else
        tag_tracer( (electron->p + electron->np-1), e_tracer, tag );
        tag_tracer( (ion->p      + ion->np-1),      H_tracer, tag );
#endif
      }
    }
  }

#else // Forcefree current sheet

#ifdef ELECTRONS_CARRRY_CURRENT
#define VDY -(b0/L)/(cosh(z/L)*cosh(z/L))  // Assumes electrons carry currents
#else
#define VDY -0.5*(b0/L)/(cosh(z/L)*cosh(z/L))  // Assumes both electrons and ions carry currents
#endif

#define VDX VDY*BX/BY
#define VD sqrt(VDX*VDX+VDY*VDY)
#define GVD 1./sqrt(1.-VD*VD/(c*c))

  sim_log( "-> Force Free Sheet" );

  repeat ( Ne/nproc() ) {
    double x, y, z, ux, uy, uz, upa1, upe1, uz1, gu1;

    x = uniform(rng(0), xmin, xmax);
    y = uniform(rng(0), ymin, ymax);
    z = uniform(rng(0), zmin, zmax);

    // inject_particles() will return an error for particles no on this
    // node and will not inject particle locally

    // Load electrons as drifting Maxwellian with velocity specified to be
    // consistent with B field

#ifdef RELATIVISTIC_PLASMA
    double fs0, f, fs, u2, gg, ue_range, ui_range;
    fs0 = 1.01;
    ue_range = factor*sqrt((Te/me) * (Te/me + 2));
    f = 1.0;
    fs = 0.0;
    do {
      upa1 = uniform(rng(0), -ue_range, ue_range);
      upe1 = uniform(rng(0), -ue_range, ue_range);
      uz1 = uniform(rng(0), -ue_range, ue_range);
      u2 = upa1*upa1 + upe1*upe1 + uz1*uz1;
      gg = sqrt(1.0 + u2);
      /* fs = exp(-u2/(1.0+gg)*me/Te); */
      fs = exp(-(gg-1)*me/Te);
      f  = uniform(rng(0), 0, 1) * fs0;
    } while(f > fs);
#else
    upa1 = normal(rng(0), 0, vthe);
    upe1 = normal(rng(0), 0, vthe);
    uz1 = normal(rng(0), 0, vthe);
#endif

    gu1 = sqrt(1.0+upa1*upa1+upe1*upe1+uz1*uz1);
    ux = (GVD*upa1*VDX/VD - upe1*VDY/VD) + GVD*VDX*gu1;
    uy = (GVD*upa1*VDY/VD + upe1*VDX/VD) + GVD*VDY*gu1;
    uz = uz1;

#ifdef TURBULENCE_MIXING_DIAGNOSTICS
    if (z > 0) {
      inject_particle(electronTop, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
    } else {
      inject_particle(electronBot, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
    }
    ++iparticle;
#else
    inject_particle(electron, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
    ++iparticle;
#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
    if (iparticle%particle_select_part_efield == 0) {
      inject_particle(electron_wo_epara, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
      inject_particle(electron_wo_eparay, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
      inject_particle(electron_wo_egtb, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
#ifdef SPECIES_ADVANCED_WITHOUT_EPERP
      inject_particle(electron_wo_eperp, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
#endif
    }
#endif
#ifdef CHECK_PARTICLE_CROSSING_EGTB
    if (iparticle%particle_select_egtb == 0) {
      inject_particle(electron_egtb, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
    }
#endif
#endif

#ifdef RELATIVISTIC_PLASMA
    ui_range = factor*sqrt((Ti/mi) * (Ti/mi + 2));
    do {
      upa1 = uniform(rng(0), -ui_range, ui_range);
      upe1 = uniform(rng(0), -ui_range, ui_range);
      uz1 = uniform(rng(0), -ui_range, ui_range);
      u2 = upa1*upa1 + upe1*upe1 + uz1*uz1;
      gg = sqrt(1.0 + u2);
      /* fs = exp(-u2/(1.0+gg)*mi/Ti); */
      fs = exp(-(gg-1)*mi/Ti);
      f  = uniform(rng(0), 0, 1) * fs0;
    } while(f > fs);
#else
    upa1 = normal(rng(0), 0, vthi);
    upe1 = normal(rng(0), 0, vthi);
    uz1 = normal(rng(0), 0, vthi);
#endif

#ifdef ELECTRONS_CARRRY_CURRENT
    // Assume ions have no drift
    ux = upa1;
    uy = upe1;
    uz = uz1;
#else
    // Assume ions have the same drift as electrons
    gu1 = sqrt(1.0+upa1*upa1+upe1*upe1+uz1*uz1);
    ux = (-GVD*upa1*VDX/VD + upe1*VDY/VD) - GVD*VDX*gu1;
    uy = (-GVD*upa1*VDY/VD - upe1*VDX/VD) - GVD*VDY*gu1;
    uz = uz1;
#endif

#ifdef TURBULENCE_MIXING_DIAGNOSTICS
    if (z > 0) {
      inject_particle(ionTop, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
    } else {
      inject_particle(ionBot, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
    }
#else
    inject_particle(ion, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
    if (iparticle%particle_select_part_efield == 0) {
      inject_particle(ion_wo_epara, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
      inject_particle(ion_wo_eparay, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
      inject_particle(ion_wo_egtb, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
#ifdef SPECIES_ADVANCED_WITHOUT_EPERP
      inject_particle(ion_wo_eperp, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
#endif
    }
#endif
#ifdef CHECK_PARTICLE_CROSSING_EGTB
    if (iparticle%particle_select_egtb == 0) {
      inject_particle(ion_egtb, x, y, z, ux*cs+uy*sn, -ux*sn+uy*cs, uz, weight, 0, 0);
    }
#endif
#endif

    if (particle_tracing == 1) { // only tag particles in the 1st pass
      if (iparticle%particle_select == 0) {
        itracer++;
        /* int tag = ((((int)rank())<<19) | (itracer & 0x7ffff)); // 13 bits (8192) for rank and 19 bits (~520k) */
        int tag = ((((int)rank())<<13) | (itracer & 0x1fff)); // 19 bits (520k) for rank and 13 bits (8192)
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
        if (z > 0) {
          tag_tracer( (electronTop->p + electronTop->np-1), e_tracer, tag );
          tag_tracer( (ionTop->p      + ionTop->np-1),      H_tracer, tag );
        } else {
          tag_tracer( (electronBot->p + electronBot->np-1), e_tracer, tag );
          tag_tracer( (ionBot->p      + ionBot->np-1),      H_tracer, tag );
        }
#else
        tag_tracer( (electron->p + electron->np-1), e_tracer, tag );
        tag_tracer( (ion->p      + ion->np-1),      H_tracer, tag );
#endif
      }
    }
  }
#endif // ifdef HARRIS_SHEET
#endif // ifdef TURBULENCE_SIMULATION

  // particle tracking
  if ( particle_tracing==2 ) { // load tracers for the 2nd pass

    union tag {
      int   tag;
      float q;
    };

    // BJA - note that this block assumes no additional user data written
    //       on the first pass, so 7 fields per particle. Since there's a
    //       postprocessing step between first and second runs, this seems
    //       adequate.

    char fname[256];
    float *ftracer_electron, *ftracer_H;
    int ntracer_electron, ntracer_H;

    // BJA - Note that we are changing this to put 0 for rank here. Each proc should
    //       attempt to write the entire tracer array and rely on the particle
    //       injection semantics to choose the right processor to write this
    //       upon - FIXME: Need to check that this works for the "inject_particle_raw"
    //       injection; if these semantics aren't there, we have to put into the loops
    //       below checks vs. global x,y,z range for domain and particle position

    sprintf( fname, "./tracer_select/T.%i/%s.%i",
             int(particle_tracing_start), "electron_tracer", 0 );
    read_tracer( ftracer_electron, ntracer_electron, fname );
    sprintf( fname, "./tracer_select/T.%i/%s.%i",
             int(particle_tracing_start), "H_tracer", 0 );
    read_tracer( ftracer_H, ntracer_H, fname );

    int       i, ix, iy, iz;
    float     x, y, z, q, ux, uy, uz;
    union tag tagp;
    i = 0;
    repeat( ntracer_electron ) {
      particle_t p[1];
      q      = ftracer_electron[7*i];
      tagp.q = q;
      x      = ftracer_electron[7*i+1];
      y      = ftracer_electron[7*i+2];
      z      = ftracer_electron[7*i+3];
      ix     = int( (x-grid->x0)/grid->dx ) + 1;
      iy     = int( (y-grid->y0)/grid->dy ) + 1;
      iz     = int( (z-grid->z0)/grid->dz ) + 1;
      p->i   = LOCAL_CELL_ID(ix,iy,iz);
      p->dx  = 2*( x - ( (ix-0.5)*grid->dx + grid->x0 ) ) / grid->dx;
      p->dy  = 2*( y - ( (iy-0.5)*grid->dy + grid->y0 ) ) / grid->dy;
      p->dz  = 2*( z - ( (iz-0.5)*grid->dz + grid->z0 ) ) / grid->dz;
      ux     = ftracer_electron[7*i+4];
      uy     = ftracer_electron[7*i+5];
      uz     = ftracer_electron[7*i+6];
      p->ux  = ux;
      p->uy  = uy;
      p->uz  = uz;

      tag_tracer( p, e_tracer, tagp.tag );

      // DEBUG
      sim_log_local( q << " q "<<tagp.q );
      sim_log_local( rank()<< " e " << i << " tagp.tag = "<<tagp.tag );
      sim_log_local( "x, y, z = "<<x<<" "<<y<<" "<<z<<" ux, uy, uz = "<<ux<<" "<<uy<<" "<<uz );

      i++;
    } // repeat


    i = 0;
    repeat( ntracer_H ) {
      particle_t p[1];
      q      = ftracer_H[7*i];
      tagp.q = q;
      x      = ftracer_H[7*i+1];
      y      = ftracer_H[7*i+2];
      z      = ftracer_H[7*i+3];
      ix     = int( (x-grid->x0)/grid->dx ) + 1;
      iy     = int( (y-grid->y0)/grid->dy ) + 1;
      iz     = int( (z-grid->z0)/grid->dz ) + 1;
      p->i   = LOCAL_CELL_ID(ix,iy,iz);
      p->dx  = 2*( x - ( (ix-0.5)*grid->dx + grid->x0 ) ) / grid->dx ;
      p->dy  = 2*( y - ( (iy-0.5)*grid->dy + grid->y0 ) ) / grid->dy ;
      p->dz  = 2*( z - ( (iz-0.5)*grid->dz + grid->z0 ) ) / grid->dz ;
      ux     = ftracer_H[7*i+4];
      uy     = ftracer_H[7*i+5];
      uz     = ftracer_H[7*i+6];
      p->ux  = ux;
      p->uy  = uy;
      p->uz  = uz;

      tag_tracer( p, H_tracer, tagp.tag );
      i++;

      // DEBUG
      sim_log_local( "H " << i << " tagp.tag = "<<tagp.tag );

    } // repeat

    // DEBUG
    if ( step()==0 ) {
      sim_log("e_tracer->np = "<<e_tracer->np);
      sim_log("H_tracer->np= "<<H_tracer->np);
    }

    delete [] ftracer_electron;
    delete [] ftracer_H;
  } // if particle_tracing==2

  sim_log( "Finished loading particles" );

  /*--------------------------------------------------------------------------
   * New dump definition
   *------------------------------------------------------------------------*/

  /*--------------------------------------------------------------------------
   * Set data output format
   *
   * This option allows the user to specify the data format for an output
   * dump.  Legal settings are 'band' and 'band_interleave'.  Band-interleave
   * format is the native storage format for data in VPIC.  For field data,
   * this looks something like:
   *
   *   ex0 ey0 ez0 div_e_err0 cbx0 ... ex1 ey1 ez1 div_e_err1 cbx1 ...
   *
   * Banded data format stores all data of a particular state variable as a
   * contiguous array, and is easier for ParaView to process efficiently.
   * Banded data looks like:
   *
   *   ex0 ex1 ex2 ... exN ey0 ey1 ey2 ...
   *
   *------------------------------------------------------------------------*/

  global->fdParams.format = band;

  sim_log ( "Fields output format = band" );

#ifdef TURBULENCE_MIXING_DIAGNOSTICS
  global->eTopdParams.format = band;
  global->eBotdParams.format = band;
#else
  global->hedParams.format = band;
#endif

  sim_log ( "Electron species output format = band" );

#ifdef TURBULENCE_MIXING_DIAGNOSTICS
  global->HTopdParams.format = band;
  global->HBotdParams.format = band;
#else
  global->hHdParams.format = band;
#endif

  sim_log ( "Ion species output format = band" );

  /*--------------------------------------------------------------------------
   * Set stride
   *
   * This option allows data down-sampling at output.  Data are down-sampled
   * in each dimension by the stride specified for that dimension.  For
   * example, to down-sample the x-dimension of the field data by a factor
   * of 2, i.e., half as many data will be output, select:
   *
   *   global->fdParams.stride_x = 2;
   *
   * The following 2-D example shows down-sampling of a 7x7 grid (nx = 7,
   * ny = 7.  With ghost-cell padding the actual extents of the grid are 9x9.
   * Setting the strides in x and y to equal 2 results in an output grid of
   * nx = 4, ny = 4, with actual extents 6x6.
   *
   * G G G G G G G G G
   * G X X X X X X X G
   * G X X X X X X X G         G G G G G G
   * G X X X X X X X G         G X X X X G
   * G X X X X X X X G   ==>   G X X X X G
   * G X X X X X X X G         G X X X X G
   * G X X X X X X X G         G X X X X G
   * G X X X X X X X G         G G G G G G
   * G G G G G G G G G
   *
   * Note that grid extents in each dimension must be evenly divisible by
   * the stride for that dimension:
   *
   *   nx = 150;
   *   global->fdParams.stride_x = 10; // legal -> 150/10 = 15
   *
   *   global->fdParams.stride_x = 8; // illegal!!! -> 150/8 = 18.75
   *------------------------------------------------------------------------*/

  // relative path to fields data from global header
#ifdef DUMP_WITH_HDF5
  sprintf(global->fdParams.baseDir, "field_hdf5");
  dump_mkdir(global->fdParams.baseDir);
#else
  dump_mkdir("fields");
  sprintf(global->fdParams.baseDir, "fields/%d",NUMFOLD);
  dump_mkdir(global->fdParams.baseDir);
#endif

  // base file name for fields output
  sprintf(global->fdParams.baseFileName, "fields");

  global->fdParams.stride_x = 1;
  global->fdParams.stride_y = 1;
  global->fdParams.stride_z = 1;

  // add field parameters to list
  global->outputParams.push_back(&global->fdParams);

  sim_log ( "Fields x-stride " << global->fdParams.stride_x );
  sim_log ( "Fields y-stride " << global->fdParams.stride_y );
  sim_log ( "Fields z-stride " << global->fdParams.stride_z );

  // relative path to electron species data from global header
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
#ifdef DUMP_WITH_HDF5
  sprintf(global->eTopdParams.baseDir, "hydro_hdf5");
  sprintf(global->eBotdParams.baseDir, "hydro_hdf5");
  dump_mkdir(global->eTopdParams.baseDir);
#else
  dump_mkdir("hydro");
  sprintf(global->eTopdParams.baseDir, "hydro/%d",NUMFOLD);
  sprintf(global->eBotdParams.baseDir, "hydro/%d",NUMFOLD);
  dump_mkdir(global->eTopdParams.baseDir);
#endif

  // base file name for fields output
  sprintf(global->eTopdParams.baseFileName, "eTophydro");
  sprintf(global->eBotdParams.baseFileName, "eBothydro");

  global->eTopdParams.stride_x = 1;
  global->eTopdParams.stride_y = 1;
  global->eTopdParams.stride_z = 1;
  global->eBotdParams.stride_x = 1;
  global->eBotdParams.stride_y = 1;
  global->eBotdParams.stride_z = 1;

  global->outputParams.push_back(&global->eTopdParams);
  global->outputParams.push_back(&global->eBotdParams);

  sim_log ( "Electron species x-stride " << global->eTopdParams.stride_x );
  sim_log ( "Electron species y-stride " << global->eTopdParams.stride_y );
  sim_log ( "Electron species z-stride " << global->eTopdParams.stride_z );

#ifdef DUMP_WITH_HDF5
  sprintf(global->HTopdParams.baseDir, "hydro_hdf5");
  sprintf(global->HBotdParams.baseDir, "hydro_hdf5");
  dump_mkdir(global->HTopdParams.baseDir);
#else
  sprintf(global->HTopdParams.baseDir, "hydro/%d",NUMFOLD);
  sprintf(global->HBotdParams.baseDir, "hydro/%d",NUMFOLD);
  dump_mkdir(global->HTopdParams.baseDir);
#endif

  // base file name for fields output
  sprintf(global->HTopdParams.baseFileName, "HTophydro");
  sprintf(global->HBotdParams.baseFileName, "HBothydro");

  global->HTopdParams.stride_x = 1;
  global->HTopdParams.stride_y = 1;
  global->HTopdParams.stride_z = 1;
  global->HBotdParams.stride_x = 1;
  global->HBotdParams.stride_y = 1;
  global->HBotdParams.stride_z = 1;

  global->outputParams.push_back(&global->HTopdParams);
  global->outputParams.push_back(&global->HBotdParams);

  sim_log ( "Ion species x-stride " << global->HTopdParams.stride_x );
  sim_log ( "Ion species y-stride " << global->HTopdParams.stride_y );
  sim_log ( "Ion species z-stride " << global->HTopdParams.stride_z );
#else // #ifdef TURBULENCE_MIXING_DIAGNOSTICS
#ifdef DUMP_WITH_HDF5
  sprintf(global->hedParams.baseDir, "hydro_hdf5");
  dump_mkdir(global->hedParams.baseDir);
#else
  sprintf(global->hedParams.baseDir, "hydro/%d",NUMFOLD);
  dump_mkdir(global->hedParams.baseDir);
#endif

  // base file name for fields output
  sprintf(global->hedParams.baseFileName, "ehydro");

  global->hedParams.stride_x = 1;
  global->hedParams.stride_y = 1;
  global->hedParams.stride_z = 1;

  // add electron species parameters to list
  global->outputParams.push_back(&global->hedParams);

  sim_log ( "Electron species x-stride " << global->hedParams.stride_x );
  sim_log ( "Electron species y-stride " << global->hedParams.stride_y );
  sim_log ( "Electron species z-stride " << global->hedParams.stride_z );

  // relative path to electron species data from global header
  dump_mkdir("hydro");
#ifdef DUMP_WITH_HDF5
  sprintf(global->hHdParams.baseDir, "hydro_hdf5");
  dump_mkdir(global->hHdParams.baseDir);
#else
  sprintf(global->hHdParams.baseDir, "hydro/%d",NUMFOLD);
  dump_mkdir(global->hHdParams.baseDir);
#endif

  // base file name for fields output
  sprintf(global->hHdParams.baseFileName, "Hhydro");

  global->hHdParams.stride_x = 1;
  global->hHdParams.stride_y = 1;
  global->hHdParams.stride_z = 1;

  sim_log ( "Ion species x-stride " << global->hHdParams.stride_x );
  sim_log ( "Ion species y-stride " << global->hHdParams.stride_y );
  sim_log ( "Ion species z-stride " << global->hHdParams.stride_z );

  // add electron species parameters to list
  global->outputParams.push_back(&global->hHdParams);
#endif // #ifdef TURBULENCE_MIXING_DIAGNOSTICS

  /*--------------------------------------------------------------------------
   * Set output fields
   *
   * It is now possible to select which state-variables are output on a
   * per-dump basis.  Variables are selected by passing an or-list of
   * state-variables by name.  For example, to only output the x-component
   * of the electric field and the y-component of the magnetic field, the
   * user would call output_variables like:
   *
   *   global->fdParams.output_variables( ex | cby );
   *
   * NOTE: OUTPUT VARIABLES ARE ONLY USED FOR THE BANDED FORMAT.  IF THE
   * FORMAT IS BAND-INTERLEAVE, ALL VARIABLES ARE OUTPUT AND CALLS TO
   * 'output_variables' WILL HAVE NO EFFECT.
   *
   * ALSO: DEFAULT OUTPUT IS NONE!  THIS IS DUE TO THE WAY THAT VPIC
   * HANDLES GLOBAL VARIABLES IN THE INPUT DECK AND IS UNAVOIDABLE.
   *
   * For convenience, the output variable 'all' is defined:
   *
   *   global->fdParams.output_variables( all );
   *------------------------------------------------------------------------*/
  /* CUT AND PASTE AS A STARTING POINT
   * REMEMBER TO ADD APPROPRIATE GLOBAL DUMPPARAMETERS VARIABLE
   *
   * XL: This only works when dump fields and hydro in binary format.

   output_variables( all );

   output_variables( electric | div_e_err | magnetic | div_b_err |
                     tca      | rhob      | current  | rhof |
                     emat     | nmat      | fmat     | cmat );

   output_variables( current_density  | charge_density |
                     momentum_density | ke_density     | stress_tensor );
   */


  global->fdParams.output_variables( electric | magnetic );
  //global->hedParams.output_variables( current_density | charge_density
  //                    | stress_tensor );
  //global->hHdParams.output_variables( current_density | charge_density
  //                    | stress_tensor );

#ifdef TURBULENCE_MIXING_DIAGNOSTICS
  global->eTopdParams.output_variables( all );
  global->eBotdParams.output_variables( all );
  global->HTopdParams.output_variables( all );
  global->HBotdParams.output_variables( all );
#else
  //global->fdParams.output_variables( all );
  global->hedParams.output_variables( all );
  global->hHdParams.output_variables( all );
#endif

  /*--------------------------------------------------------------------------
   * Convenience functions for simlog output
   *------------------------------------------------------------------------*/

  char varlist[512];
  create_field_list(varlist, global->fdParams);

  sim_log ( "Fields variable list: " << varlist );

#ifdef TURBULENCE_MIXING_DIAGNOSTICS
  create_hydro_list(varlist, global->eTopdParams);
  sim_log ( "Electron top species variable list: " << varlist );

  create_hydro_list(varlist, global->eBotdParams);
  sim_log ( "Electron bot species variable list: " << varlist );

  create_hydro_list(varlist, global->HTopdParams);
  sim_log ( "Ion top species variable list: " << varlist );

  create_hydro_list(varlist, global->HBotdParams);
  sim_log ( "Ion bot species variable list: " << varlist );
#else
  create_hydro_list(varlist, global->hedParams);
  sim_log ( "Electron species variable list: " << varlist );

  create_hydro_list(varlist, global->hHdParams);
  sim_log ( "Ion species variable list: " << varlist );
#endif

  /* ---------------------------------------------
     now add parameters for the energy diagnostics
     --------------------------------------------- */

#ifdef TURBULENCE_MIXING_DIAGNOSTICS
  global->edeTop.sp_id = electronTop->id;
  global->edeTop.vth = sqrt(2.0)*vthe;
  sprintf(global->edeTop.fname,global->eTopdParams.baseFileName);

  global->edeBot.sp_id = electronBot->id;
  global->edeBot.vth = sqrt(2.0)*vthe;
  sprintf(global->edeBot.fname,global->eBotdParams.baseFileName);

  global->edHTop.sp_id = ionTop->id;
  global->edHTop.vth = sqrt(2.0)*vthi;
  sprintf(global->edHTop.fname, global->HTopdParams.baseFileName);

  global->edHBot.sp_id = ionBot->id;
  global->edHBot.vth = sqrt(2.0)*vthi;
  sprintf(global->edHBot.fname, global->HBotdParams.baseFileName);
#else
  global->ede.sp_id = electron->id;
  global->ede.vth = sqrt(2.0)*vthe;
  sprintf(global->ede.fname, global->hedParams.baseFileName);

  global->edH.sp_id = ion->id;
  global->edH.vth = sqrt(2.0)*vthi;
  sprintf(global->edH.fname, global->hHdParams.baseFileName);
#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
  global->ede_wo_epara.sp_id = electron_wo_epara->id;
  global->ede_wo_epara.vth = sqrt(2.0)*vthe;
  sprintf(global->ede_wo_epara.fname, "e_wo_epara_hydro");

  global->edH_wo_epara.sp_id = ion_wo_epara->id;
  global->edH_wo_epara.vth = sqrt(2.0)*vthi;
  sprintf(global->edH_wo_epara.fname, "H_wo_epara_hydro");

  global->ede_wo_eparay.sp_id = electron_wo_eparay->id;
  global->ede_wo_eparay.vth = sqrt(2.0)*vthe;
  sprintf(global->ede_wo_eparay.fname, "e_wo_eparay_hydro");

  global->edH_wo_eparay.sp_id = ion_wo_eparay->id;
  global->edH_wo_eparay.vth = sqrt(2.0)*vthi;
  sprintf(global->edH_wo_eparay.fname, "H_wo_eparay_hydro");

  global->ede_wo_egtb.sp_id = electron_wo_egtb->id;
  global->ede_wo_egtb.vth = sqrt(2.0)*vthe;
  sprintf(global->ede_wo_egtb.fname, "e_wo_egtb_hydro");

  global->edH_wo_egtb.sp_id = ion_wo_egtb->id;
  global->edH_wo_egtb.vth = sqrt(2.0)*vthi;
  sprintf(global->edH_wo_egtb.fname, "H_wo_egtb_hydro");

#ifdef SPECIES_ADVANCED_WITHOUT_EPERP
  global->ede_wo_eperp.sp_id = electron_wo_eperp->id;
  global->ede_wo_eperp.vth = sqrt(2.0)*vthe;
  sprintf(global->ede_wo_eperp.fname, "e_wo_eperp_hydro");

  global->edH_wo_eperp.sp_id = ion_wo_eperp->id;
  global->edH_wo_eperp.vth = sqrt(2.0)*vthi;
  sprintf(global->edH_wo_eperp.fname, "H_wo_eperp_hydro");
#endif
#endif
#ifdef CHECK_PARTICLE_CROSSING_EGTB
  global->ede_egtb.sp_id = electron_egtb->id;
  global->ede_egtb.vth = sqrt(2.0)*vthe;
  sprintf(global->ede_egtb.fname, "e_egtb_hydro");

  global->edH_egtb.sp_id = ion_egtb->id;
  global->edH_egtb.vth = sqrt(2.0)*vthi;
  sprintf(global->edH_egtb.fname, "H_egtb_hydro");
#endif
#endif // #ifdef TURBULENCE_MIXING_DIAGNOSTICS

  global->emax_band  = emax_band;
  global->emin_band  = emin_band;
  global->nbands     = nbands;
  global->emax_spect = emax_spect;
  global->emin_spect = emin_spect;
  global->nbins      = nbins;
  global->nx_zone    = nx_zone;
  global->ny_zone    = ny_zone;
  global->nz_zone    = nz_zone;

  sim_log("*** Finished with user-specified initialization ***");


  // Upon completion of the initialization, the following occurs:
  // - The synchronization error (tang E, norm B) is computed between domains
  //   and tang E / norm B are synchronized by averaging where discrepancies
  //   are encountered.
  // - The initial divergence error of the magnetic field is computed and
  //   one pass of cleaning is done (for good measure)
  // - The bound charge density necessary to give the simulation an initially
  //   clean divergence e is computed.
  // - The particle momentum is uncentered from u_0 to u_{-1/2}
  // - The user diagnostics are called on the initial state
  // - The physics loop is started
  //
  // The physics loop consists of:
  // - Advance particles from x_0,u_{-1/2} to x_1,u_{1/2}
  // - User particle injection at x_{1-age}, u_{1/2} (use inject_particles)
  // - User current injection (adjust field(x,y,z).jfx, jfy, jfz)
  // - Advance B from B_0 to B_{1/2}
  // - Advance E from E_0 to E_1
  // - User field injection to E_1 (adjust field(x,y,z).ex,ey,ez,cbx,cby,cbz)
  // - Advance B from B_{1/2} to B_1
  // - (periodically) Divergence clean electric field
  // - (periodically) Divergence clean magnetic field
  // - (periodically) Synchronize shared tang e and norm b
  // - Increment the time step
  // - Call user diagnostics
  // - (periodically) Print a status message

} //begin_initialization

#define should_dump(x) \
  (global->x##_interval>0 && remainder(step(), global->x##_interval) == 0)

#define should_dump_pre(x) \
  (global->x##_interval>0 && remainder(step(), global->x##_interval) == -1)

#define should_dump_post(x) \
  (global->x##_interval>0 && remainder(step(), global->x##_interval) == 1)

begin_diagnostics {
  /*--------------------------------------------------------------------------
   * NOTE: YOU CANNOT DIRECTLY USE C FILE DESCRIPTORS OR SYSTEM CALLS ANYMORE
   *
   * To create a new directory, use:
   *
   *   dump_mkdir("full-path-to-directory/directoryname")
   *
   * To open a file, use: FileIO class
   *
   * Example for file creation and use:
   *
   *   // declare file and open for writing
   *   // possible modes are: io_write, io_read, io_append,
   *   // io_read_write, io_write_read, io_append_read
   *   FileIO fileIO;
   *   FileIOStatus status;
   *   status= fileIO.open("full-path-to-file/filename", io_write);
   *
   *   // formatted ASCII  output
   *   fileIO.print("format string", varg1, varg2, ...);
   *
   *   // binary output
   *   // Write n elements from array data to file.
   *   // T is the type, e.g., if T=double
   *   // fileIO.write(double * data, size_t n);
   *   // All basic types are supported.
   *   fileIO.write(T * data, size_t n);
   *
   *   // close file
   *   fileIO.close();
   *------------------------------------------------------------------------*/

  /*--------------------------------------------------------------------------
   * Data output directories
   * WARNING: The directory list passed to "global_header" must be
   * consistent with the actual directories where fields and species are
   * output using "field_dump" and "hydro_dump".
   *
   * DIRECTORY PATHES SHOULD BE RELATIVE TO
   * THE LOCATION OF THE GLOBAL HEADER!!!
   *------------------------------------------------------------------------*/


  /*--------------------------------------------------------------------------
   * Normal rundata dump
   *------------------------------------------------------------------------*/
  if(step()==0) {
    dump_mkdir("data");
    dump_mkdir("fields");
    dump_mkdir("hydro");
    dump_mkdir("particle");
    dump_mkdir("rundata");
    dump_mkdir("restore0");
    dump_mkdir("restore1");  // 1st backup
    dump_mkdir("restore2");  // 2nd backup
    dump_mkdir("tracer");
    dump_mkdir("tracer/tracer1");
    dump_mkdir("tracer/tracer2");
    dump_mkdir("tracer/traj1");
    dump_mkdir("tracer/traj2");

    // Make subfolders for restart
    char restorefold[128];
    sprintf(restorefold, "restore0/%i", NUMFOLD);
    dump_mkdir(restorefold);
    sprintf(restorefold, "restore1/%i", NUMFOLD);
    dump_mkdir(restorefold);
    sprintf(restorefold, "restore2/%i", NUMFOLD);
    dump_mkdir(restorefold);

    // And rundata
    char rundatafold[128];
    char rundatafile[128];
    sprintf(rundatafold, "rundata/%i", NUMFOLD);
    dump_mkdir(rundatafold);
    sprintf(rundatafile, "rundata/%i/grid", NUMFOLD);
    dump_grid(rundatafile);
    dump_materials("rundata/materials");
    dump_species("rundata/species");
    global_header("global", global->outputParams);
  } // if

  /*--------------------------------------------------------------------------
   * particle tracking
   * -----------------------------------------------------------------------*/

  // Set TRACER_ACCUM_HYDRO to 1 if we need to accumulate hydro moments before
  // writing trajectory output. Since this involves a pass through all the particles
  // in the system as well as synchronization (which hits MPI), don't do this step
  // unless we must.
#undef  TRACER_DO_ACCUM_HYDRO
#define TRACER_DO_ACCUM_HYDRO (0)       //CHANGE BETWEEN PASSES

//  // Setup data needed for hydro output
//# ifdef TRACER_DO_ACCUM_HYDRO
//    TRACER_HYDRO_SETUP( e, "electron" )
//    TRACER_HYDRO_SETUP( H, "H"       )
//# endif

  // Be careful! This number should be set correctly
#undef  TRACER_NUM_ADDED_FIELDS
#define TRACER_NUM_ADDED_FIELDS (11)       //CHANGE BETWEEN PASSES

#undef CALC_TRACER_USER_DEFINED_DATA
#define CALC_TRACER_USER_DEFINED_DATA                           \
  if ( global->emf_at_tracer ) {                                \
    CALC_EMFS_AT_TRACER;                                        \
  }                                                             \
  if ( global->hydro_at_tracer ) {                              \
      CALC_HYDRO_FIELDS_AT_TRACER;                              \
  }

  // We assume hydro fields are alway behind electric and magnetic fields
#undef TRACER_USER_DEFINED_DATA
#define TRACER_USER_DEFINED_DATA                                \
  if ( global->emf_at_tracer ) {                                \
    pout[index + 6 + 1]  = ex;                                  \
    pout[index + 6 + 2]  = ey;                                  \
    pout[index + 6 + 3]  = ez;                                  \
    pout[index + 6 + 4]  = bx;                                  \
    pout[index + 6 + 5]  = by;                                  \
    pout[index + 6 + 6]  = bz;                                  \
    if ( global->hydro_at_tracer ) {                            \
      pout[index + 12 + 1] = vx;                                \
      pout[index + 12 + 2] = vy;                                \
      pout[index + 12 + 3] = vz;                                \
      pout[index + 12 + 4] = ne;                                \
      pout[index + 12 + 5] = ni;                                \
      if ( global->ve_at_tracer ) {                             \
        pout[index + 12 + 6] = vex;                             \
        pout[index + 12 + 7] = vey;                             \
        pout[index + 12 + 8] = vez;                             \
      }                                                         \
    }                                                           \
  } else {                                                      \
    if ( global->hydro_at_tracer ) {                            \
      pout[index + 6 + 1] = vx;                                 \
      pout[index + 6 + 2] = vy;                                 \
      pout[index + 6 + 3] = vz;                                 \
      pout[index + 6 + 4] = ne;                                 \
      pout[index + 6 + 5] = ni;                                 \
      if ( global->ve_at_tracer ) {                             \
        pout[index + 6 + 6] = vex;                              \
        pout[index + 6 + 7] = vey;                              \
        pout[index + 6 + 8] = vez;                              \
      }                                                         \
    }                                                           \
  }

  // Hydro fields at tracer positions
  static hydro_array_t * hydro_tot_array;
  static hydro_t * ALIGNED(128) htot;
  static hydro_t * ALIGNED(128) hi;
  static hydro_t * RESTRICT ALIGNED(16) htot0;
  static hydro_t * RESTRICT ALIGNED(16) h0;
  const int nx = grid->nx;
  const int ny = grid->ny;
  const int nz = grid->nz;
  int frame;

  const int tracer_ratio1 = global->tracer_pass1_interval / global->tracer_interval;
  const int tracer_ratio2 = global->tracer_pass2_interval / global->tracer_interval;

  // initialize buffered tracer data
  if ( step() == 0 || (step()>1 && step()==global->restart_step+1) ) {
    if ( global->particle_tracing==1 && tracer_ratio1 > 1 ) {
      init_buffered_tracers(tracer_ratio1);
    } else if ( global->particle_tracing==2 && tracer_ratio2 > 1 ){
      init_buffered_tracers(tracer_ratio2);
    }
    if ( global->particle_tracing > 0 && (step()>1 &&
         step()==global->restart_step+1 && (tracer_ratio1 > 1 || tracer_ratio2 > 1)) ) {
      read_buffered_tracer_restart(global->rtoggle);
    }
  }

  // Accumulate hydro
  if ( global->particle_tracing > 0 && global->hydro_at_tracer ) {
    if ( step() == 0 || (step()>1 && step()==global->restart_step+1) ) {
      hydro_tot_array = new_hydro_array(grid);
      UNREGISTER_OBJECT(hydro_tot_array);
    }
  }

  if( global->particle_tracing > 0 && should_dump(tracer) ) {
    if ( global->hydro_at_tracer ) {  // accumulate hydro at tracer positions
      int x, y, z;
      float rho_tot;
      clear_hydro_array(hydro_tot_array);
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
      species_t * sp = find_species_name("electronTop", species_list);
      accumulate_hydro_p(hydro_tot_array, sp, interpolator_array);
      synchronize_hydro_array(hydro_tot_array);
      sp = find_species_name("electronBot", species_list);
      accumulate_hydro_p(hydro_array, sp, interpolator_array);
      synchronize_hydro_array(hydro_array);
      htot = hydro_tot_array->h;
      hi   = hydro_array->h;
      for (z = 1; z <= nz + 1; z++) {
        for (y = 1; y <= ny + 1; y++) {
          htot0 = &HYDRO_TOT(1, y, z);
          h0    = &HYDRO(1, y, z);
          for (x = 1; x <= nx + 1; x++) {
            htot0->rho += h0->rho;
            htot0->jx += h0->jx;
            htot0->jy += h0->jy;
            htot0->jz += h0->jz;
            // we use txx, tyy, and tzz as electron bulk velocity
            if (global->ve_at_tracer && fabs(htot0->rho) > 0) {
              htot0->txx = htot0->jx / htot0->rho;
              htot0->tyy = htot0->jy / htot0->rho;
              htot0->tzz = htot0->jz / htot0->rho;
            }
            htot0++;
            h0++;
          }
        }
      }

      sp = find_species_name("ionTop", species_list);
      clear_hydro_array(hydro_array);
      accumulate_hydro_p(hydro_array, sp, interpolator_array);
      synchronize_hydro_array(hydro_array);
      htot = hydro_tot_array->h;
      hi   = hydro_array->h;
      for (z = 1; z <= nz + 1; z++) {
        for (y = 1; y <= ny + 1; y++) {
          htot0 = &HYDRO_TOT(1, y, z);
          h0    = &HYDRO(1, y, z);
          for (x = 1; x <= nx + 1; x++) {
            // ne*ve + ni_top * vi_top * mi
            htot0->jx = -htot0->jx + h0->jx*global->mi_me;
            htot0->jy = -htot0->jy + h0->jy*global->mi_me;
            htot0->jz = -htot0->jz + h0->jz*global->mi_me;
            htot0->rho = fabs(htot0->rho); // Electron number density
            htot0->px = h0->rho;           // Ion number density at Top
            htot0++;
            h0++;
          }
        }
      }

      sp = find_species_name("ionBot", species_list);
      clear_hydro_array(hydro_array);
      accumulate_hydro_p(hydro_array, sp, interpolator_array);
      synchronize_hydro_array(hydro_array);

      htot = hydro_tot_array->h;
      hi   = hydro_array->h;
      for (z = 1; z <= nz + 1; z++) {
        for (y = 1; y <= ny + 1; y++) {
          htot0 = &HYDRO_TOT(1, y, z);
          h0    = &HYDRO(1, y, z);
          for (x = 1; x <= nx + 1; x++) {
            // px here is the ion number density at Top
            rho_tot = htot0->rho + (htot0->px + h0->rho) * global->mi_me;
            // jx, jy, jz are actually vx, vy, vz now
            htot0->jx = (htot0->jx + h0->jx*global->mi_me) / rho_tot;
            htot0->jy = (htot0->jy + h0->jy*global->mi_me) / rho_tot;
            htot0->jz = (htot0->jz + h0->jz*global->mi_me) / rho_tot;
            htot0->px += h0->rho; // Ion number density including both Top and Bottom
            htot0++;
            h0++;
          }
        }
      }
#else
      species_t * sp = find_species_name("electron", species_list);
      accumulate_hydro_p(hydro_tot_array, sp, interpolator_array);
      synchronize_hydro_array(hydro_tot_array);
      sp = find_species_name("ion", species_list);
      clear_hydro_array(hydro_array);
      accumulate_hydro_p(hydro_array, sp, interpolator_array);
      synchronize_hydro_array(hydro_array);
      htot = hydro_tot_array->h;
      hi   = hydro_array->h;
      for (z = 1; z <= nz + 1; z++) {
        for (y = 1; y <= ny + 1; y++) {
          htot0 = &HYDRO_TOT(1, y, z);
          h0    = &HYDRO(1, y, z);
          for (x = 1; x <= nx + 1; x++) {
            // we use txx, tyy, and tzz as electron bulk velocity
            if (global->ve_at_tracer && fabs(htot0->rho) > 0) {
              htot0->txx = htot0->jx / htot0->rho;
              htot0->tyy = htot0->jy / htot0->rho;
              htot0->tzz = htot0->jz / htot0->rho;
            }
            // Assuming electron has -1 charge, ion has +1 charge
            rho_tot = fabs(htot0->rho) + h0->rho * global->mi_me;
            // jx, jy, jz are actually vx, vy, vz for single fluid now
            htot0->jx = (-htot0->jx + h0->jx*global->mi_me) / rho_tot;
            htot0->jy = (-htot0->jy + h0->jy*global->mi_me) / rho_tot;
            htot0->jz = (-htot0->jz + h0->jz*global->mi_me) / rho_tot;
            htot0->rho = fabs(htot0->rho); // Electron number density
            htot0->px = h0->rho;           // Ion number density
            htot0++;
            h0++;
          }
        }
      }
#endif
    } // if global->hydro_at_tracer

    // Buffer tracer data for latter dump
    if (global->particle_tracing==1 && tracer_ratio1 > 1) {
      frame = ((step() % global->tracer_pass1_interval)-1) / global->tracer_interval;
      if (frame < 0) frame = 0;
      buffer_tracers(tracer_ratio1, frame);
    } else if (global->particle_tracing==2 && tracer_ratio2 > 1) {
      frame = ((step() % global->tracer_pass2_interval)-1) / global->tracer_interval;
      if (frame < 0) frame = 0;
      buffer_tracers(tracer_ratio2, frame);
    }

  } // if should_dump(tracer)

  if ( global->particle_tracing==1 ) {           // First pass
    if ( should_dump(tracer_pass1) || step() == num_step) {
      //if ( TRACER_DO_ACCUM_HYDRO ) {
      //  // accumulate electron hydro
      //  TRACER_ACCUM_HYDRO( e );
      //  // accumulate H hydro
      //  TRACER_ACCUM_HYDRO( H );
      //} // if
      if (global->dump_traj_directly) {
        dump_traj("tracer/traj1");
      } else {
        if (tracer_ratio1 == 1) { // tracer data is not buffered
          //dump_tracers("tracer/tracer1");
/* #include "dumptracer_h5part.cc" */
/* #include "dumptracer_hdf5.cc" */
#include "dumptracer_hdf5_single.cc"
        } else {
          dump_buffered_tracer(tracer_ratio1, "tracer/tracer1");
          clear_buffered_tracers(tracer_ratio1);
        }
      }
    } // if
  } else if ( global->particle_tracing==2 ) {    // Second pass
    if ( should_dump(tracer_pass2) || step() == num_step) {
      //if ( TRACER_DO_ACCUM_HYDRO ) {
      //  // accumulate electron hydro
      //  TRACER_ACCUM_HYDRO( e );
      //  // accumulate H hydro
      //  TRACER_ACCUM_HYDRO( H );
      //}  // if
      if (global->dump_traj_directly) {
        dump_traj("tracer/traj2");
      } else {
        if (tracer_ratio2 == 1) { // tracer data is not buffered
          //dump_tracers("tracer/tracer2");
/* #include "dumptracer_h5part.cc" */
/* #include "dumptracer_hdf5.cc" */
#include "dumptracer_hdf5_single.cc"
        } else {
          dump_buffered_tracer(tracer_ratio2, "tracer/tracer2");
          clear_buffered_tracers(tracer_ratio2);
        }
      }
    }  // if
  } // if global->particle_tracing

  /*--------------------------------------------------------------------------
   * Normal rundata energies dump
   *------------------------------------------------------------------------*/
  if(should_dump(energies)) {
    dump_energies("rundata/energies", step() == 0 ? 0 : 1);
  } // if

  /*--------------------------------------------------------------------------
   * Field data output
   *------------------------------------------------------------------------*/

#ifdef DUMP_WITH_HDF5
  field_dump_flag.disableE();
  field_dump_flag.disableCB();
  field_dump_flag.disableTCA();
  field_dump_flag.disableJF();
  field_dump_flag.disableEMAT();
  field_dump_flag.disableFMAT();
  field_dump_flag.ex = true;
  field_dump_flag.ey = true;
  field_dump_flag.ez = true;
  field_dump_flag.cbx = true;
  field_dump_flag.cby = true;
  field_dump_flag.cbz = true;

#ifdef DUMP_THREE_SUCCESSIVE_STEPS
  if(should_dump_pre(fields)) {
    double time_to_dump_fields = uptime();
    dump_fields_hdf5(global->fdParams.baseDir, 0);
    time_to_dump_fields = uptime() - time_to_dump_fields;
    sim_log("Time in dumping fields: "<< time_to_dump_fields << " s");
  }
  if(should_dump_post(fields)) {
    double time_to_dump_fields = uptime();
    dump_fields_hdf5(global->fdParams.baseDir, 0);
    time_to_dump_fields = uptime() - time_to_dump_fields;
    sim_log("Time in dumping fields: "<< time_to_dump_fields << " s");
  }
#endif
  if(step() == 1 || should_dump(fields)) {
    double time_to_dump_fields = uptime();
    dump_fields_hdf5(global->fdParams.baseDir, 0);
    time_to_dump_fields = uptime() - time_to_dump_fields;
    sim_log("Time in dumping fields: "<< time_to_dump_fields << " s");
  }

#else // #ifdef DUMP_WITH_HDF5

#ifdef DUMP_THREE_SUCCESSIVE_STEPS
  if(should_dump_pre(fields)) {
    double time_to_dump_fields = uptime();
    field_dump(global->fdParams);
    time_to_dump_fields = uptime() - time_to_dump_fields;
    sim_log("Time in dumping fields: "<< time_to_dump_fields << " s");
  }
  if(should_dump_post(fields)) {
    double time_to_dump_fields = uptime();
    field_dump(global->fdParams);
    time_to_dump_fields = uptime() - time_to_dump_fields;
    sim_log("Time in dumping fields: "<< time_to_dump_fields << " s");
  }
#endif
  if(step() == 1 || should_dump(fields)) {
    double time_to_dump_fields = uptime();
    field_dump(global->fdParams);
    time_to_dump_fields = uptime() - time_to_dump_fields;
    sim_log("Time in dumping fields: "<< time_to_dump_fields << " s");
  }

#endif // #ifdef DUMP_WITH_HDF5


#ifdef DUMP_WITH_HDF5
  hydro_dump_flag.resetToDefaults();

  /*--------------------------------------------------------------------------
   * Electron species output
   *------------------------------------------------------------------------*/

#ifdef TURBULENCE_MIXING_DIAGNOSTICS

#ifdef DUMP_THREE_SUCCESSIVE_STEPS
  if(should_dump_pre(ehydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("electronTop", global->eTopdParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump_post(ehydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("electronTop", global->eTopdParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump_pre(ehydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("electronBot", global->eBotdParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump_post(ehydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("electronBot", global->eBotdParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
#endif
  if(should_dump(ehydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("electronTop", global->eTopdParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump(ehydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("electronBot", global->eBotdParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }

#else // #ifdef TURBULENCE_MIXING_DIAGNOSTICS

#ifdef DUMP_THREE_SUCCESSIVE_STEPS
  if(should_dump_pre(ehydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("electron", global->hedParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump_post(ehydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("electron", global->hedParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
#endif
  if(should_dump(ehydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("electron", global->hedParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
#endif

  /*--------------------------------------------------------------------------
   * Ion species output
   *------------------------------------------------------------------------*/

#ifdef TURBULENCE_MIXING_DIAGNOSTICS

#ifdef DUMP_THREE_SUCCESSIVE_STEPS
  if(should_dump_pre(Hhydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("ionTop", global->HTopdParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump_post(Hhydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("ionTop", global->HTopdParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump_pre(Hhydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("ionBot", global->HBotdParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump_post(Hhydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("ionBot", global->HBotdParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
#endif
  if(should_dump(Hhydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("ionTop", global->HTopdParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump(Hhydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("ionBot", global->HBotdParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }

#else // #ifdef TURBULENCE_MIXING_DIAGNOSTICS

#ifdef DUMP_THREE_SUCCESSIVE_STEPS
  if(should_dump_pre(Hhydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("ion", global->hHdParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump_post(Hhydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("ion", global->hHdParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
#endif
  if(should_dump(Hhydro)) {
    double time_to_dump_hydro = uptime();
    dump_hydro_hdf5("ion", global->hHdParams.baseDir, 0);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
#endif

#else // #ifdef DUMP_WITH_HDF5

#ifdef TURBULENCE_MIXING_DIAGNOSTICS

#ifdef DUMP_THREE_SUCCESSIVE_STEPS
  if(should_dump_pre(ehydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("electronTop", global->eTopdParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump_post(ehydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("electronTop", global->eTopdParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump_pre(ehydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("electronBot", global->eBotdParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump_post(ehydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("electronBot", global->eBotdParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
#endif
  if(should_dump(ehydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("electronTop", global->eTopdParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump(ehydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("electronBot", global->eBotdParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }

#else // #ifdef TURBULENCE_MIXING_DIAGNOSTICS

#ifdef DUMP_THREE_SUCCESSIVE_STEPS
  if(should_dump_pre(ehydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("electron", global->hedParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump_post(ehydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("electron", global->hedParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
#endif
  if(should_dump(ehydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("electron", global->hedParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }

#endif // #ifdef TURBULENCE_MIXING_DIAGNOSTICS

  /*--------------------------------------------------------------------------
   * Ion species output
   *------------------------------------------------------------------------*/

#ifdef TURBULENCE_MIXING_DIAGNOSTICS

#ifdef DUMP_THREE_SUCCESSIVE_STEPS
  if(should_dump_pre(Hhydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("ionTop", global->HTopdParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump_post(Hhydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("ionTop", global->HTopdParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump_pre(Hhydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("ionBot", global->HBotdParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump_post(Hhydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("ionBot", global->HBotdParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
#endif
  if(should_dump(Hhydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("ionTop", global->HTopdParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump(Hhydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("ionBot", global->HBotdParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }

#else // #ifdef TURBULENCE_MIXING_DIAGNOSTICS

#ifdef DUMP_THREE_SUCCESSIVE_STEPS
  if(should_dump_pre(Hhydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("ion", global->hHdParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
  if(should_dump_post(Hhydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("ion", global->hHdParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }
#endif
  if(should_dump(Hhydro)) {
    double time_to_dump_hydro = uptime();
    hydro_dump("ion", global->hHdParams);
    time_to_dump_hydro = uptime() - time_to_dump_hydro;
    sim_log("Time in dumping hydro: "<< time_to_dump_hydro << " s");
  }

#endif // #ifdef TURBULENCE_MIXING_DIAGNOSTICS

#endif // #ifdef DUMP_WITH_HDF5

  /*--------------------------------------------------------------------------
   * Time-averaging field and hydro data output
   *------------------------------------------------------------------------*/
//#include "time_average_master.cc"

  /*--------------------------------------------------------------------------
   * Energy Spectrum Output
   *------------------------------------------------------------------------*/

#include "energy_local.cc"   // Subroutine to compute energy spectrum diagnostic

  //Vadim:
  //#include "dissipation.cxx"
  //#include "Ohms_exp_all_v2.cxx"

  /*--------------------------------------------------------------------------
   * Restart dump
   *------------------------------------------------------------------------*/

  if ( step()>0 && should_dump(restart) ) {
    static const char * restart_fbase[2] = { "restore1", "restore2" };
    double dumpstart = uptime();

    if ( global->rtoggle == 2 ) global->rtoggle = 0;

    //checkpt_subdir( restart_fbase[global->rtoggle], NUMFOLD );
    checkpt_subdir( "restore1", NUMFOLD );
    //checkpt_subdir( restart_fbase[global->rtoggle], 0 );

    // particle tracking
    if (global->particle_tracing > 0) {
      //dump_tracer_restart(global->rtoggle);
      dump_tracer_restart(1);
    }

    // dump buffered data for restart
    if (global->particle_tracing==1 && tracer_ratio1 > 1) {
      dump_buffered_tracer_restart(tracer_ratio1, global->rtoggle);
    } else if (global->particle_tracing==2 && tracer_ratio2 > 1) {
      dump_buffered_tracer_restart(tracer_ratio2, global->rtoggle);
    }

    double dumpelapsed = uptime() - dumpstart;

    sim_log("Restart duration "<<dumpelapsed);

    global->rtoggle^=1;

    global->restart_step = step();
    free_buffered_tracers();

    // free allocated memory for new hydro array
    if ( global->particle_tracing > 0 && global->hydro_at_tracer ) {
      FREE_ALIGNED( hydro_tot_array->h );
      FREE( hydro_tot_array );
    }

    if (rank() == 0) {
      FileIO fp_restart_info;
      if ( ! (fp_restart_info.open("latest_restart", io_write)==ok) )
        ERROR(("Cannot open file."));
      fp_restart_info.print("restart");
      fp_restart_info.print(restart_fbase[global->rtoggle]);
      fp_restart_info.close();
    }
  } // if

  // Shut down simulation if wall clock time exceeds global->quota_sec.
  // Note that the mp_elapsed() is guaranteed to return the same value for all
  // processors (elapsed time on proc #0), and therefore the abort should
  // be synchronized across processors.

  if ( step()>0 && global->quota_check_interval && (step()%global->quota_check_interval)==0 ) {
    if ( uptime() > global->quota_sec ) {

      global->rtoggle = 2;

      //checkpt_subdir( "restore0", 0 );
      //checkpt_subdir( "restore0", NUMFOLD );
      checkpt_subdir( "restore1", NUMFOLD );

      // particle tracking
      if (global->particle_tracing > 0) {
        //dump_tracer_restart(global->rtoggle);
        dump_tracer_restart(1);
      }

      // dump buffered tracer
#include "dumptracer_hdf5_single.cc"

      // dump buffered data for restart
      if (global->particle_tracing==1 && tracer_ratio1 > 1) {
        dump_buffered_tracer_restart(tracer_ratio1, global->rtoggle);
      } else if (global->particle_tracing==2 && tracer_ratio2 > 1) {
        dump_buffered_tracer_restart(tracer_ratio2, global->rtoggle);
      }

      sim_log( "Restart dump restart completed." );
      sim_log( "Allowed runtime exceeded for this job.  Terminating." );
      mp_barrier(); // Just to be safe

      global->restart_step = step();
      free_buffered_tracers();

      // free allocated memory for new hydro array
      if ( global->particle_tracing > 0 && global->hydro_at_tracer ) {
        FREE_ALIGNED( hydro_tot_array->h );
        FREE( hydro_tot_array );
      }

      if (rank() == 0) {
        FileIO fp_restart_info;
        if ( ! (fp_restart_info.open("latest_restart", io_write)==ok) )
          ERROR(("Cannot open file."));
        fp_restart_info.print("restart restore0");
        fp_restart_info.close();
      }
      exit(0);
    } // if
  } // if

  // The end of the simulation
  if ( step() == num_step ) {
    // free allocated memory for new hydro array
    if ( global->particle_tracing > 0 && global->hydro_at_tracer ) {
      FREE_ALIGNED( hydro_tot_array->h );
      FREE( hydro_tot_array );
    }
  }

  /*--------------------------------------------------------------------------
   * Particle dump
   *------------------------------------------------------------------------*/
  char subdir[36];
  if ( should_dump(eparticle) && step() !=0 &&
       step() > 0*(global->fields_interval)  ) {
    sprintf(subdir,"particle/T.%d",step());
    dump_mkdir(subdir);
    /* sprintf(subdir,"particle/T.%d/eparticle",step()); */
    /* dump_particles("electron",subdir); */
    /* sprintf(subdir,"particle/T.%d/hparticle",step()); */
    /* dump_particles("ion",subdir); */
    /* #include "dump_with_h5part.cc" */
    /* #include "dump_with_hdf5.cc" */
  }

  if ( should_dump(Hparticle) && step() !=0 &&
       step() > 0*(global->fields_interval)  ) {
    sprintf(subdir,"particle/T.%d",step());
    dump_mkdir(subdir);
    /* sprintf(subdir,"particle/T.%d/hparticle",step()); */
    /* dump_particles("ion",subdir); */
  }

} // end diagnostics

begin_particle_injection {
  // particle tracking
  // Note: read_tracer_restart is called if needed in advance_tracers()

  // if ( global->particle_tracing > 0 ) advance_tracers(global->rtoggle);
  if ( global->particle_tracing > 0 ) advance_tracers(1);

  static int a_initted_efield = 0;
  static accumulator_array_t *dummy_a_efield;
  species_t *s;
#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
  if (step() == 0) {
    s = global->species_list_efield;
    while( s ){
      uncenter_p( s, interpolator_array ); // necessary, or particles with be numerically heated
      s = s->next;
    }
  }
  if (a_initted_efield == 0){
    dummy_a_efield = new_accumulator_array(grid);
    UNREGISTER_OBJECT(dummy_a_efield);
    a_initted_efield = 1;
  }
  s = global->species_list_efield;
  while( s ){
    if ((s->sort_interval>0) && (step()%s->sort_interval==0)) sort_p(s);
    s = s->next;
  }
  s = global->species_list_efield;
  for (int isp = 0; isp < global->nsp_efield; isp++) {   // loop over species
    if (isp < global->nsp_efield-6) {        // remove perpendicular electric field but keep ExB drift
      advance_p_efield(s, dummy_a_efield, interpolator_array, 3);
    } else if (isp < global->nsp_efield-4) { // remove electric field in regions with |E| > |B|
      advance_p_efield(s, dummy_a_efield, interpolator_array, 2);
    } else if (isp < global->nsp_efield-2) { // remove the y-component of parallel electric field
      advance_p_efield(s, dummy_a_efield, interpolator_array, 1);
    } else if (isp < global->nsp_efield) { // remove parallel electric field
      advance_p_efield(s, dummy_a_efield, interpolator_array, 0);
    } else {
      advance_p_efield(s, dummy_a_efield, interpolator_array, 100);
    }
    s = s->next;
  }
  s = global->species_list_efield;
  for ( int npass=0; npass<num_comm_round; ++npass ) {
    boundary_p(particle_bc_list, s, field_array, dummy_a_efield);
  }
  s = global->species_list_efield;
#endif

#ifdef CHECK_PARTICLE_CROSSING_EGTB
  if (step() == 0) {
    s = global->species_list_egtb;
    while( s ){
      uncenter_p( s, interpolator_array ); // necessary, or particles with be numerically heated
      s = s->next;
    }
  }
  // advance particles
  if (a_initted_efield == 0){
    dummy_a_efield = new_accumulator_array(grid);
    UNREGISTER_OBJECT(dummy_a_efield);
    a_initted_efield = 1;
  }
  s = global->species_list_egtb;
  while( s ){
    if ((s->sort_interval>0) && (step()%s->sort_interval==0)) sort_p(s);
    s = s->next;
  }
  s = global->species_list_egtb;
  while( s ){
    advance_p(s, dummy_a_efield, interpolator_array);
    s = s->next ;
  }
  s = global->species_list_egtb;
  for ( int npass=0; npass<num_comm_round; ++npass ) {
    boundary_p(particle_bc_list, s, field_array, dummy_a_efield);
  }
  s = global->species_list_egtb;

  // check whether particles cross E>B regions
  float dx0, dy0, dz0;
  float ex, ey, ez, bx, by, bz;
  float etot2, btot2;
  int ii;
  particle_t * ALIGNED(32) p;
  const interpolator_t * ALIGNED(16) f0=interpolator_array->i;
  const interpolator_t * ALIGNED(16) f;
  while( s ){
    p = s->p;
    for (std::size_t i = 0; i < s->np; ++i ) { // loop over particles
      dx0 = p->dx;
      dy0 = p->dy;
      dz0 = p->dz;
      ii  = p->i;
      f   = f0 + ii;
      ex  = f->ex + dy0*f->dexdy + dz0*(f->dexdz+dy0*f->d2exdydz);
      ey  = f->ey + dz0*f->deydz + dx0*(f->deydx+dz0*f->d2eydzdx);
      ez  = f->ez + dx0*f->dezdx + dy0*(f->dezdy+dx0*f->d2ezdxdy);
      bx  = f->cbx + dx0*f->dcbxdx;
      by  = f->cby + dy0*f->dcbydy;
      bz  = f->cbz + dz0*f->dcbzdz;
      etot2 = ex * ex + ey * ey + ez * ez;
      btot2 = bx * bx + by * by + bz * bz;
      if (etot2 > btot2) p->w = 100;  // tag the particle with a large weight
      p++;
    }
    s = s->next;
  } /* while */
#endif

#ifdef OPEN_BOUNDARY_CONDITIONS
#include "open_boundary.cc"
#endif // OPEN_BOUNDARY_CONDITIONS

} // end particle injection

begin_current_injection {
} // end current injection

begin_field_injection {
#ifdef OPEN_BOUNDARY_CONDITIONS
  const int nx=grid->nx;
  const int ny=grid->ny;
  const int nz=grid->nz;
  double t=grid->dt*step();
  double tau = global->tdrive;
  int x,y,z;

// There macros are from local.c to apply boundary conditions

#define XYZ_LOOP(xl,xh,yl,yh,zl,zh) \
  for( z=zl; z<=zh; z++ )           \
    for( y=yl; y<=yh; y++ )         \
      for( x=xl; x<=xh; x++ )

#define xy_EDGE_LOOP(z) XYZ_LOOP(1,nx,1,ny+1,z,z)
#define yx_EDGE_LOOP(z) XYZ_LOOP(1,nx+1,1,ny,z,z)

  // Top Boundary
  if (global->top) {
    yx_EDGE_LOOP(nz+1) field(x,y,z).ey = (global->edrive)*(1-exp(-t/tau));
    xy_EDGE_LOOP(nz+1) field(x,y,z).ex = -(global->edrive)*(1-exp(-t/tau))*global->bg;
  }

  // Bottom Boundary
  if (global->bottom) {
    yx_EDGE_LOOP(1) field(x,y,z).ey = (global->edrive)*(1-exp(-t/tau));
    xy_EDGE_LOOP(1) field(x,y,z).ex = (global->edrive)*(1-exp(-t/tau))*global->bg;
  }
#endif // whether to use open boundary conditions
} // end field injection

begin_particle_collisions {
} // end particle collisions

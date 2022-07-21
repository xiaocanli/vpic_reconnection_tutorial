/* Modified tracer header file.
 *
 * Changes to legacy version:
 * - Now a proper header file with safeties for multiple inclusion
 * - Added user-configurable number of fields and a facility for
 *   adding user-defined fields to tracer particle and trajectory
 *   dumps
 * - Added status checks on files to make more graceful error
 *   reporting when read/write files can't be opened
 * - Removed first-pass buffered write macros
 * - Added loop to advance_tracers boundary_p calls to reflect
 *   num_comm_round semantics
 * - Added hydro accumulation macros + switch to turn these on
 *   as necessary
 * - Added undefs of several of the helper macros local to this
 *   file.
 * - Tracers changed to dump checkpoint/restart into
 *   restart/restart[tag] where tag = 0, 1, 2
 * - We now use fseek (-like) file writing of trajectory files
 *   to ensure we are checkpoint/restart safe without introducing
 *   glitches into the traj files
 *
 * FIXME:
 * - Add typesafing of macro parameters for safety. (Obviously,
 *   this will be limited - as C/C++ macros are just a lexical
 *   substitution, any exceptions would occur at runtime.)
 *
 * Mods added by Brian Albright, XTD-PRI, 4 Aug. 2016
 *--------------------------------------------------------------*/

// Old comments:
/* This file contains all the necessary routines for tagging
 * and advancing tracer particles in VPIC. The methedology here
 * is to copy all tagged particles into a new tracer species.
 * This tracer species does not back react on the fields/plasma
 * so that the evolution of this tracer species is identical
 * to that of the actual simulation particles. The reason for
 * copying is that we can then use the q subfield in the
 * particle structure as a tag for the tracers, but by altering
 * q we cannot accuratly accumulate moments and thus there
 * can be no back-reaction and so we need to copy. The particle
 * push is unaffected since this depends on q/m which is
 * unmodified and specified external to the particle structure.
 *
 *
 *********************************************************************
 * Example of how these routines should be called form the input deck:
 *********************************************************************
 * begin_globals{
 *      ...
 *      species_t * tracers_list ;
 *      ...
 * }
 *
 * begin_initilization{
 *  ...
 *      define_species(,...);       // Because I'm lazy, tracer species
 *      ...                             // must be defined after all other
 *      ...                 // species !!
 *      tag_tracer(...);        // Call to tag/create a new tracer particle
 *      ...
 *      hijack_tracers(...);        // MUST BE CALLED AFTER DEFINING TRACERS
 *  ...
 * }
 *
 * begin_diagnostics{
 *      ...
 *      tag_tracer(...);                // Call to tag/create a new tracer particle
 *      ...
 *      dump_tracers(...);              // Call to dump tracers when required
 *      ...
 *      dump_tracer_restart();          // Call to dump tracers to a restart file...
 *      ...
 * }
 *
 * begin_particle_injection{
 *      advance_tracers();              // MUST BE CALLED AT EVERY TIME STEP
 *      ...
 *  }
 *
 */

//--------------------------------------------------------------
#ifndef __TRACER_INTERP_HXX__
#define __TRACER_INTERP_HXX__


//--------------------------------------------------------------
// General purpose allocation macro
#ifndef ALLOCATE
#define ALLOCATE(A,LEN,TYPE)                                                  \
  if ( !((A)=(TYPE *)malloc((size_t)(LEN)*sizeof(TYPE))) )                    \
    ERROR(("Cannot allocate."));
#endif // ALLOCATE

//--------------------------------------------------------------
// Turn on some useful output for debugging
#ifndef TRACER_VERBOSE
#define TRACER_VERBOSE 0
#endif // TRACER_VERBOSE

//--------------------------------------------------------------
// Default number of fields written per particle

#ifndef TRACER_NUM_FIELDS_BASE
#define TRACER_NUM_FIELDS_BASE (7)
#endif // TRACER_NUM_FIELDS_BASE

//--------------------------------------------------------------
// Allow user to override using the following syntax:
//
// #undef  TRACER_NUM_ADDED_FIELDS
// #define TRACER_NUM_ADDED_FIELDS (N)
//
// where N >= 0

#ifndef TRACER_NUM_ADDED_FIELDS
#define TRACER_NUM_ADDED_FIELDS (0)
#endif // TRACER_NUM_ADDED_FIELDS

//--------------------------------------------------------------
// User can add additional fields to write per particle. This is
// potentially useful if data such as, e.g., kinetic energy, Lorentz
// factor, or hydro dump data are desired. Syntax is as follows:
//
// #undef CALC_TRACER_USER_DEFINED_DATA
// #define CALC_TRACER_USER_DEFINED_DATA                    \
//     ex = ...                                             \
//     ey = ...                                             \
//
// #undef TRACER_USER_DEFINED_DATA
// #define TRACER_USER_DEFINED_DATA \
//         pout[index + 12 + 1]                       = ... ; \
//         pout[index + 12 + 2]                       = ... ; \
//         ...
//         pout[index + 12 + TRACER_NUM_ADDED_FIELDS] = ... ;
//

#ifndef CALC_TRACER_USER_DEFINED_DATA
#define CALC_TRACER_USER_DEFINED_DATA
#endif // CALC_TRACER_USER_DEFINED_DATA

#ifndef TRACER_USER_DEFINED_DATA
#define TRACER_USER_DEFINED_DATA
#endif // TRACER_USER_DEFINED_DATA

//--------------------------------------------------------------
// Redefine TRACER_DO_ACCUM_HYDRO to 1 in the input deck if we need to accumulate
// hydro moments before writing trajectory output. Since this involves a pass
// through all the particles in the system as well as synchronization (which
// hits MPI), don't do this step unless we need these data!

#ifndef TRACER_DO_ACCUM_HYDRO
#define TRACER_DO_ACCUM_HYDRO 0
#endif // #ifndef TRACER_ACCUM_HYDRO

//--------------------------------------------------------------
// Declare variables and allocate space for hydro arrays. TAG
// identifies which hydro species type is used; ID is the
// species name (type const char[]).  Note: because namespace needs i
// to be declared outside this block, there's no BEGIN_PRIMITIVE

#ifndef TRACER_HYDRO_SETUP
#define TRACER_HYDRO_SETUP(TAG,ID)                                            \
    static int TAG##_initted = 0;                                             \
    static species_t * TAG##_species;                                         \
    static hydro_t * ALIGNED(128) TAG##_hydro;                                \
    hydro_array_t  * TAG##_hydro_array;                                       \
    if ( ! TAG##_initted ) {                                                  \
      TAG##_initted     = 1;                                                  \
      TAG##_species     = find_species(ID);                                   \
      TAG##_hydro_array = new_hydro_array( grid );                            \
      UNREGISTER_OBJECT( TAG##_hydro_array );                                 \
    } // if initted
#endif // TRACER_HYDRO_SETUP

//--------------------------------------------------------------
// Macro to accumulate hydro data. TAG identifies which
// hydro moments before writing trajectory output. Since this involves a pass
// through all the particles in the system as well as synchronization (which
// hits MPI), don't do this step unless we need these data!

#ifndef TRACER_ACCUM_HYDRO
#define TRACER_ACCUM_HYDRO(TAG)                                               \
    BEGIN_PRIMITIVE {                                                         \
      clear_hydro_array( TAG##_hydro_array );                                 \
      accumulate_hydro_p( TAG##_hydro_array,                                  \
                          TAG##_species,                                      \
                          interpolator_array );                               \
      synchronize_hydro_array( TAG##_hydro_array );                           \
    } END_PRIMITIVE
#endif // TRACER_ACCUM_HYDRO

//--------------------------------------------------------------
// Here, p is the particle to copy and tag and tracer is the species
// into which we will inject the tracer. If only one species is
// defined then tracer = tracers_list. tag shold be a unique
// 4 byte long identifier for this tracer; it is the only way to
// distinguish tracers! that
//
// Warning:  inject_particle_raw doesn't check for free storage
// space, so we're silently assuming all nodes have enough room to
// hold all tracers
//

#define tag_tracer(p, tracer, tag) BEGIN_PRIMITIVE{                           \
    float q = * reinterpret_cast<float *>( &tag ) ;                           \
    double xmin = grid->x0;                                                   \
    double xmax = grid->x0+grid->nx*grid->dx;                                 \
    double ymin = grid->y0;                                                   \
    double ymax = grid->y0+grid->ny*grid->dy;                                 \
    double zmin = grid->z0;                                                   \
    double zmax = grid->z0+grid->nz*grid->dz;                                 \
                                                                              \
    if (    x >= xmin && x <= xmax                                            \
         && y >= ymin && y <= ymax                                            \
         && z >= zmin && z <= zmax ) {                                        \
      inject_particle_raw(tracer, p->dx, p->dy, p->dz, p->i,                  \
                          p->ux, p->uy, p->uz, q);                            \
    } /* if */                                                                \
} END_PRIMITIVE

//--------------------------------------------------------------
// BJA - This makes use of VPIC species semantics. Provided tracers
// are defined last so they're at the head of the species_list linked
// list, we can snip them from the list and evolve them separately
// so they are not accumulated by the field solve. However, this
// puts the onus on the user to handle checkpoint/restart of
// species in the input deck.
//
// FIXME: It would be cleaner and less error prone if this
//        this functionality were enabled in the VPIC source.
//

#define hijack_tracers( num_tracer_species )                                  \
  BEGIN_PRIMITIVE{                                                            \
    species_t *s         = species_list ;                                     \
    global->tracers_list = species_list ;                                     \
    species_t *s0;                                                            \
                                                                              \
    /* advance num_tracer_species elements down species_list */               \
    for(int i = 1; i < num_tracer_species ; ++i) {                            \
      s0 = s;                                                                 \
      s = s->next;                                                            \
      UNREGISTER_OBJECT( s0 );                                                \
    }                                                                         \
                                                                              \
    /* set species_list to point to the first physical */                     \
    /* particle species */                                                    \
    species_list = s->next;                                                   \
    UNREGISTER_OBJECT( s );                                                   \
    s->next      = NULL ;                                                     \
  } END_PRIMITIVE

//--------------------------------------------------------------
// advance_p takes care of all local moves, but to resolve cross-domain
// moves and boundary interactions we need to call boundary_p as well.
//
// Note: We've configured this to pass through the boundary handler
//       num_comm_round times instead of assuming three passes maximum.
//       This is done in preparation for being able to process tracers
//       through randomizing boundary handlers such as maxwellian_reflux
//
// BJA - Note: randomizing boundary handlers (e.g., maxwellian_reflux)
//       cannot be used with tracer particles presently since tracers'
//       trajectories will be randomized upon interaction with the
//       boundary in a non-deterministic fashion. FIXME: in the
//       new head version, have each particle generate its own
//       private RNG sequence
//

#define advance_tracers(RESTART_TAG) BEGIN_PRIMITIVE{                         \
    static int a_initted = 0;                                                 \
    static accumulator_array_t *dummy_a;                                      \
    if (a_initted == 0){                                                      \
      dummy_a = new_accumulator_array(grid);                                  \
      UNREGISTER_OBJECT(dummy_a);                                             \
      if(step()) read_tracer_restart(RESTART_TAG);                            \
      a_initted = 1;                                                          \
    } /* if */                                                                \
                                                                              \
    species_t *s = global->tracers_list ;                                     \
    while( s ){                                                               \
      if ((s->sort_interval>0) && (step()%s->sort_interval==0)) sort_p(s);    \
      s = s->next ;                                                           \
    } /* while */                                                             \
    s = global->tracers_list ;                                                \
    while( s ){                                                               \
      advance_p(s, dummy_a, interpolator_array);                              \
      s = s->next ;                                                           \
    } /* while */                                                             \
    s = global->tracers_list ;                                                \
    for ( int npass=0; npass<num_comm_round; ++npass ) {                      \
      boundary_p(particle_bc_list, s, field_array, dummy_a);                  \
    } /* for */                                                               \
    s = global->tracers_list ;                                                \
} END_PRIMITIVE

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
template <typename T> int is_negative(T val) {
    return (val < T(0));
}

//--------------------------------------------------------------
// Electric and magnetic fields at tracer positions
#ifndef CALC_EMFS_AT_TRACER
#define CALC_EMFS_AT_TRACER                                                   \
    ex  = f->ex + dy0*f->dexdy + dz0*(f->dexdz+dy0*f->d2exdydz);              \
    ey  = f->ey + dz0*f->deydz + dx0*(f->deydx+dz0*f->d2eydzdx);              \
    ez  = f->ez + dx0*f->dezdx + dy0*(f->dezdy+dx0*f->d2ezdxdy);              \
    bx  = f->cbx + dx0*f->dcbxdx;                                             \
    by  = f->cby + dy0*f->dcbydy;                                             \
    bz  = f->cbz + dz0*f->dcbzdz
#endif

//--------------------------------------------------------------
// Single fluid velocity fields and number densities
// NOTE that we reuse jx, jy, jz for velocities and rho for electron
// number density, and px for ion number density
#ifndef VELS_RHOS
#define VELS_RHOS(wn, i, j, k)                                                \
    hy   = &HYDRO_TOT(i, j, k);                                               \
    vx  += wn*hy->jx;                                                         \
    vy  += wn*hy->jy;                                                         \
    vz  += wn*hy->jz;                                                         \
    vex += wn*hy->txx;                                                        \
    vey += wn*hy->tyy;                                                        \
    vez += wn*hy->tzz;                                                        \
    ne  += wn*hy->rho;                                                        \
    ni  += wn*hy->px
#endif

//--------------------------------------------------------------
// Hydro fields at tracer positions
#ifndef CALC_HYDRO_FIELDS_AT_TRACER
#define CALC_HYDRO_FIELDS_AT_TRACER                                           \
    /* Compute the trilinear coefficients */                                  \
    dx=dx0; dy=dy0; dz=dz0;                                                   \
    w0  = r8V;      /* w0 = 1/(8V) = w/8               */                     \
    dx *= w0;       /* dx = wx                         */                     \
    w1  = w0+dx;    /* w1 = w/8 + wx/8 = (w/8)(1+x)    */                     \
    w0 -= dx;       /* w0 = w/8 - wx/8 = (w/8)(1-x)    */                     \
    w3  = 1+dy;     /* w3 = 1+y                        */                     \
    w2  = w0*w3;    /* w2 = (w/8)(1-x)(1+y)            */                     \
    w3 *= w1;       /* w3 = (w/8)(1+x)(1+y)            */                     \
    dy  = 1-dy;     /* dy = 1-y                        */                     \
    w0 *= dy;       /* w0 = (w/8)(1-x)(1-y)            */                     \
    w1 *= dy;       /* w1 = (w/8)(1+x)(1-y)            */                     \
    w7  = 1+dz;     /* w7 = 1+z                        */                     \
    w4  = w0*w7;    /* w4 = (w/8)(1-x)(1-y)(1+z) *Done */                     \
    w5  = w1*w7;    /* w5 = (w/8)(1+x)(1-y)(1+z) *Done */                     \
    w6  = w2*w7;    /* w6 = (w/8)(1-x)(1+y)(1+z) *Done */                     \
    w7 *= w3;       /* w7 = (w/8)(1+x)(1+y)(1+z) *Done */                     \
    dz  = 1-dz;     /* dz = 1-z                        */                     \
    w0 *= dz;       /* w0 = (w/8)(1-x)(1-y)(1-z) *Done */                     \
    w1 *= dz;       /* w1 = (w/8)(1+x)(1-y)(1-z) *Done */                     \
    w2 *= dz;       /* w2 = (w/8)(1-x)(1+y)(1-z) *Done */                     \
    w3 *= dz;       /* w3 = (w/8)(1+x)(1+y)(1-z) *Done */                     \
    iz = ii / nx2ny2;                                                         \
    iy = (ii % nx2ny2) / nx2;                                                 \
    ix = ii % nx2;                                                            \
    vx = 0.0; vy = 0.0; vz = 0.0;                                             \
    vex = 0.0; vey = 0.0; vez = 0.0;                                          \
    ne = 0.0; ni = 0.0;                                                       \
    VELS_RHOS(w0, ix,   iy,   iz);                                            \
    VELS_RHOS(w1, ix+1, iy,   iz);                                            \
    VELS_RHOS(w2, ix,   iy+1, iz);                                            \
    VELS_RHOS(w3, ix+1, iy+1, iz);                                            \
    VELS_RHOS(w4, ix,   iy,   iz+1);                                          \
    VELS_RHOS(w5, ix+1, iy,   iz+1);                                          \
    VELS_RHOS(w6, ix,   iy+1, iz+1);                                          \
    VELS_RHOS(w7, ix+1, iy+1, iz+1)
#endif

//--------------------------------------------------------------
// Symbols used in macros below - Note: need to undef at the end
// of include file so we don't pollute namespace

#define nxg_tracer (grid->nx + 2)
#define nyg_tracer (grid->ny + 2)
#define nzg_tracer (grid->nz + 2)
#define i0_tracer (ii%nxg_tracer)
#define j0_tracer ((ii/nxg_tracer)%nyg_tracer)
#define k0_tracer (ii/(nxg_tracer*nyg_tracer))
#define tracer_x ((i0_tracer + (dx0-1)*0.5) * grid->dx + grid->x0)
#define tracer_y ((j0_tracer + (dy0-1)*0.5) * grid->dy + grid->y0)
#define tracer_z ((k0_tracer + (dz0-1)*0.5) * grid->dz + grid->z0)

//--------------------------------------------------------------
// This dump routine converts particle positions to global
// rather than local coordinates. User defined custom data
// fields are enabled as described above.
//
// Note: extra field added for number of fields per particle

#define dump_tracers(fbase)                                                   \
  BEGIN_PRIMITIVE {                                                           \
    char dname[256], fname[256] ;                                             \
    float ex, ey, ez, bx, by, bz;                                             \
    float dx0, dy0, dz0;                                                      \
    float ux, uy, uz;                                                         \
    int ii, j, nvar, index;                                                   \
    float dx, dy, dz;                                                         \
    float vx, vy, vz;                                                         \
    float ne, ni;                                                             \
    float w0, w1, w2, w3, w4, w5, w6, w7;                                     \
    int ix, iy, iz;                                                           \
    species_t *s = global->tracers_list ;                                     \
    const int nx2 = grid->nx + 2;                                             \
    const int nx2ny2 = (grid->ny+2) * nx2;                                    \
    const particle_t     * ALIGNED(32) p;                                     \
    const interpolator_t * ALIGNED(16) f0=interpolator_array->i;              \
    const interpolator_t * ALIGNED(16) f;                                     \
    const grid_t * g = grid;                                                  \
    const hydro_t * ALIGNED(32) hy;                                           \
    const float r8V = 0.125;                                                  \
    float *pout;                                                              \
    FileIO       fh;                                                          \
    FileIOStatus status;                                                      \
    sprintf(dname, "%s/T.%d", fbase, step() );                                \
    dump_mkdir(dname);                                                        \
                                                                              \
    nvar = TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS;                  \
                                                                              \
    while( s ) {                                                              \
      if ( s->np > 0 ) {                                                      \
        pout = new float[s->np*nvar];                                         \
        for ( p=s->p, j=0; j<s->np; j++, p++ ) {                              \
          index = j*nvar;                                                     \
          dx0 = p->dx;                                                        \
          dy0 = p->dy;                                                        \
          dz0 = p->dz;                                                        \
          ii  = p->i;                                                         \
          ux  = p->ux;                                                        \
          uy  = p->uy;                                                        \
          uz  = p->uz;                                                        \
          f   = f0 + ii;                                                      \
          CALC_TRACER_USER_DEFINED_DATA;                                      \
          pout[index + 0]  = (float) p->w;                                    \
          pout[index + 1]  = (float) tracer_x ;                               \
          pout[index + 2]  = (float) tracer_y ;                               \
          pout[index + 3]  = (float) tracer_z ;                               \
          pout[index + 4]  = ux;                                              \
          pout[index + 5]  = uy;                                              \
          pout[index + 6]  = uz;                                              \
          TRACER_USER_DEFINED_DATA;                                           \
        }                                                                     \
        sprintf(fname, "%s/%s.%d", dname , s->name, (int)rank());             \
        status = fh.open(fname, io_write);                                    \
        if ( status == fail ) ERROR(("Could not open file %s", fname));       \
        WRITE(int,    s->np ,    fh ) ;                                       \
        WRITE(int,    nvar,      fh ) ;                                       \
        WRITE(float,  s->q ,     fh ) ;                                       \
        WRITE(float,  s->m ,     fh ) ;                                       \
        WRITE(double, step()*grid->dt , fh ) ;                                \
        WRITE(float,  1.0 ,      fh ) ;                                       \
        WRITE(double, 1.0 ,      fh ) ;                                       \
        fh.write(pout, nvar*s->np);                                           \
        fh.close();                                                           \
        delete [] pout;                                                       \
      }                                                                       \
      s = s->next;                                                            \
    }                                                                         \
  } END_PRIMITIVE

//--------------------------------------------------------------
// In input deck, need to execute this dump_tracer_restart macro
// every time we do checkpoint/restart.

#define dump_tracer_restart(RESTART_TAG)                                      \
  BEGIN_PRIMITIVE {                                                           \
    species_t *s = global->tracers_list ;                                     \
    char fname[256] ;                                                         \
    int numspecies = 0;                                                       \
    FileIO       f;                                                           \
    FileIOStatus status;                                                      \
    sprintf(fname, "restore%d/%d/tracer%d.%d", RESTART_TAG,                   \
            NUMFOLD, global->particle_tracing, (int)rank() ) ;                \
    status = f.open(fname, io_write);                                         \
    if ( status == fail ) ERROR(("Could not open file %s", fname));           \
                                                                              \
    while( s ){                                                               \
        ++numspecies;                                                         \
        s = s->next;                                                          \
    }                                                                         \
    WRITE(int, numspecies, f);                                                \
    s = global->tracers_list;                                                 \
    while( s ){                                                               \
        WRITE_STRING(s->name, f);                                             \
        WRITE(species_id, s->id, f);                                          \
        WRITE(float, s->q, f);                                                \
        WRITE(float, s->m, f);                                                \
        WRITE(int, s->max_np, f);                                             \
        WRITE(int, s->max_nm, f);                                             \
        WRITE(int, s->sort_interval, f);                                      \
        WRITE(int, s->sort_out_of_place, f);                                  \
        WRITE(int, s->np, f);                                                 \
        if (s->np > 0 ) f.write(s->p, s->np);                                 \
        s = s->next;                                                          \
    }                                                                         \
    f.close();                                                                \
  } END_PRIMITIVE

//--------------------------------------------------------------
// In input deck, need to execute this read_tracer_restart macro
// before advancing tracers in begin_particle_injection block
//
// FIXME: Need to test that tag = 2 (code runs to completion &
//        writes restart dumps) restarts properly

#define read_tracer_restart(RESTART_TAG)                                      \
  BEGIN_PRIMITIVE {                                                           \
    char fname[256] ;                                                         \
    FileIO f;                                                                 \
    FileIOStatus status;                                                      \
    int namelen, maxnp, maxnm, sort ;                                         \
    int sorttype, np, numspecies, i ;                                         \
    float q, m;                                                               \
    species_id id;                                                            \
    char *name;                                                               \
    species_t * s;                                                            \
                                                                              \
    global->tracers_list = NULL ;                                             \
    sprintf(fname, "restore%d/%d/tracer%d.%d", RESTART_TAG,                   \
            NUMFOLD, global->particle_tracing, (int)rank() ) ;                \
    status = f.open(fname, io_read);                                          \
    if ( status == fail ) ERROR(("Could not open file %s", fname));           \
                                                                              \
    READ(int, numspecies, f);                                                 \
    for( i=0 ; i < numspecies ; ++i){                                         \
        READ(int, namelen, f) ;                                               \
        name = (char *)malloc(namelen+1) ;                                    \
        f.read(name, namelen) ;                                               \
        name[namelen] = '\0' ;                                                \
        READ(species_id,    id,   f);                                         \
        READ(float,     q,   f);                                              \
        READ(float,     m,   f);                                              \
        READ(int,       maxnp,    f);                                         \
        READ(int,       maxnm,    f);                                         \
        READ(int,       sort,     f);                                         \
        READ(int,       sorttype, f);                                         \
        READ(int,       np,   f);                                             \
        s = append_species( species(name, q, m, maxnp, maxnm, sort, sorttype, \
                    grid), &(global->tracers_list) );                         \
        s->id = id ;                                                          \
        s->np = np ;                                                          \
        if ( np > 0 ) f.read(s->p, np);                                       \
    }                                                                         \
    f.close();                                                                \
  } END_PRIMITIVE

//--------------------------------------------------------------
// dump tracer by particle trajectory - requires global->tracer2_interval
// to be defined
//
// FIXME: buffered write
//

/* #   undef ACCUM_HYDRO */
// dump tracer by particle trajectory
#define dump_traj(fbase) BEGIN_PRIMITIVE{                                     \
    char dname[256], fname[256] ;                                             \
    char sp_name[16], tracer_name[16];                                        \
    char *pch;                                                                \
    species_t *s = global->tracers_list ;                                     \
    float ex, ey, ez, bx, by, bz;                                             \
    float dx0, dy0, dz0;                                                      \
    float ux, uy, uz, q;                                                      \
    int ii, n, nvar, index, tag;                                              \
    float dx, dy, dz;                                                         \
    float vx, vy, vz;                                                         \
    float vex, vey, vez;                                                      \
    float ne, ni;                                                             \
    float w0, w1, w2, w3, w4, w5, w6, w7;                                     \
    int ix, iy, iz;                                                           \
    const int nx2 = grid->nx + 2;                                             \
    const int nx2ny2 = (grid->ny+2) * nx2;                                    \
    const particle_t     * ALIGNED(32) p;                                     \
    const particle_t     * ALIGNED(32) p0;                                    \
    const interpolator_t * ALIGNED(16) f0=interpolator_array->i;              \
    const interpolator_t * ALIGNED(16) f;                                     \
    const grid_t  * g = grid;                                                 \
    const hydro_t * ALIGNED(32) hy;                                           \
    const float r8V = 0.125;                                                  \
    float *pout;                                                              \
    FileIO fh;                                                                \
    FileIOStatus status;                                                      \
                                                                              \
    sprintf(dname, "%s", fbase );                                             \
    dump_mkdir(dname);                                                        \
    nvar = TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS;                  \
    index = 0;                                                                \
    pout = new float[nvar];                                                   \
    while( s ){                                                               \
      n = s->np;                                                              \
      if ( n > 0 ){                                                           \
        p0 = s->p;                                                            \
        for ( p=p0; n; n--, p++ ){                                            \
          dx0 = p->dx;                                                        \
          dy0 = p->dy;                                                        \
          dz0 = p->dz;                                                        \
          ii = p->i;                                                          \
          ux = p->ux;                                                         \
          uy = p->uy;                                                         \
          uz = p->uz;                                                         \
          q = p->w;                                                           \
          int tag = *reinterpret_cast<int*>(&q);                              \
          f = f0 + ii;                                                        \
          if (tag != 0) {                                                     \
            CALC_TRACER_USER_DEFINED_DATA;                                    \
            sprintf(fname, "%s/%s.%i", dname , s->name, tag);                 \
            status = fh.open(fname,io_append);                                \
            if ( status == fail ) ERROR(("Could not open file %s", fname));   \
            pout[index + 0] = step()*grid->dt ;                               \
            pout[index + 1] = (float) tracer_x ;                              \
            pout[index + 2] = (float) tracer_y ;                              \
            pout[index + 3] = (float) tracer_z ;                              \
            pout[index + 4] = ux;                                             \
            pout[index + 5] = uy;                                             \
            pout[index + 6] = uz;                                             \
            TRACER_USER_DEFINED_DATA;                                         \
            fh.write(pout,nvar);                                              \
            fh.close();                                                       \
          }                                                                   \
        }                                                                     \
      }                                                                       \
      s = s->next;                                                            \
    }                                                                         \
    delete [] pout;                                                           \
} END_PRIMITIVE

//--------------------------------------------------------------

#define read_tracer(arr,npart,fname)                                          \
  BEGIN_PRIMITIVE {                                                           \
    float        qm, const1;                                                  \
    double       timestep, const2;                                            \
    FileIO       fp;                                                          \
    FileIOStatus status;                                                      \
                                                                              \
    status = fp.open( fname, io_read );                                       \
    if ( status == fail )                                                     \
      ERROR(("Could not open file %s", fname));                               \
    fp.read( &npart,    1 );                                                  \
    fp.read( &qm,       1 );                                                  \
    fp.read( &timestep, 1 );                                                  \
    fp.read( &const1,   1 );                                                  \
    fp.read( &const2,   1 );                                                  \
    if ( TRACER_VERBOSE != 0 ) {                                              \
      sim_log_local("npart:                        "<<npart);                 \
      sim_log_local("qm:                           "<<qm);                    \
      sim_log_local("time:                         "<<timestep);              \
      sim_log_local("sanity check (should be 1 1): "<<const1<<" "<<const2);   \
    }                                                                         \
    if ( const1 != 1.0 || const2 != 1.0 ) {                                   \
      ERROR(("Mangled header: const1 = %f, const2 = %f\n", const1, const2));  \
      exit(1);                                                                \
    }                                                                         \
    ALLOCATE( arr, npart*7, float );                                          \
    if ( npart > 0 ) {                                                        \
      fp.read( arr, npart*7 );                                                \
      sim_log_local( npart << " loaded from " << fname );                     \
      if ( TRACER_VERBOSE != 0 ) {                                            \
        for (int tmpval=0; tmpval<npart*7; ++tmpval ) {                       \
          sim_log(tmpval<<" "<<arr[tmpval]);                                  \
        }                                                                     \
      }                                                                       \
    }                                                                         \
  } END_PRIMITIVE

static float  *pout;
static int    *np_tracer;
static double *t_tracer;
static int    *ntraj_point_tracer;

//--------------------------------------------------------------
// This routine initializes buffered tracer information.

#define init_buffered_tracers(nframes)                                        \
  BEGIN_PRIMITIVE{                                                            \
    static int numspecies = 0;                                                \
    static int psize = 0;                                                     \
    static int ntracer = global->Ntracer / nproc();                           \
    static int nvar = TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS;       \
    numspecies = num_species(global->tracers_list);                           \
    psize = 2*numspecies*ntracer*nvar*nframes;                                \
    pout      = new float[psize];                                             \
    np_tracer = new int[nframes*numspecies];                                  \
    t_tracer  = new double[nframes];                                          \
    ntraj_point_tracer = new int[numspecies];                                 \
    clear_buffered_tracers(nframes);                                          \
} END_PRIMITIVE

//--------------------------------------------------------------
// This routine clears buffered tracer information.

#define clear_buffered_tracers(nframes)                                       \
  BEGIN_PRIMITIVE{                                                            \
    static int numspecies = 0;                                                \
    static int psize = 0;                                                     \
    static int ntracer = global->Ntracer / nproc();                           \
    static int nvar = TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS;       \
    numspecies = num_species(global->tracers_list);                           \
    psize = 2*numspecies*ntracer*nvar*nframes;                                \
    memset(pout,      0, sizeof(float) * psize);                              \
    memset(np_tracer, 0, sizeof(int) * nframes*numspecies);                   \
    memset(t_tracer,  0, sizeof(double) * nframes);                           \
    memset(ntraj_point_tracer, 0, sizeof(int) * numspecies);                  \
} END_PRIMITIVE

//--------------------------------------------------------------
// This routine frees buffered tracer information.

#define free_buffered_tracers()                                               \
  BEGIN_PRIMITIVE{                                                            \
    delete [] ntraj_point_tracer;                                             \
    delete [] pout;                                                           \
    delete [] np_tracer;                                                      \
    delete [] t_tracer;                                                       \
} END_PRIMITIVE

//--------------------------------------------------------------
// This routine buffers tracer information. Theses tracers will be
// dumped every few steps. This reduces I/O operations.

#define buffer_tracers(nframes, frame)                                        \
  BEGIN_PRIMITIVE{                                                            \
    const particle_t     * ALIGNED(32) p;                                     \
    const particle_t     * ALIGNED(32) p0;                                    \
    const interpolator_t * ALIGNED(16) f;                                     \
    const interpolator_t * ALIGNED(16) f0=interpolator_array->i;              \
    const grid_t * g = grid;                                                  \
    const int nx2 = grid->nx + 2;                                             \
    const int nx2ny2 = (grid->ny+2) * nx2;                                    \
    const hydro_t * ALIGNED(32) hy;                                           \
    const float r8V = 0.125;                                                  \
    species_t *s = global->tracers_list ;                                     \
    static int numspecies = 0;                                                \
    static int ntracer = global->Ntracer / nproc();                           \
    int ii, j, n, nvar, isp, index;                                           \
    float ex, ey, ez, bx, by, bz;                                             \
    float dx0, dy0, dz0;                                                      \
    float ux, uy, uz, q;                                                      \
    float dx, dy, dz;                                                         \
    float vx, vy, vz;                                                         \
    float vex, vey, vez;                                                      \
    float ne, ni;                                                             \
    float w0, w1, w2, w3, w4, w5, w6, w7;                                     \
    int ix, iy, iz;                                                           \
                                                                              \
    numspecies = num_species(s);                                              \
    nvar = TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS;                  \
                                                                              \
    t_tracer[frame] = step() * grid->dt;                                      \
    index = 2 * numspecies * ntracer * nvar * frame;                          \
    isp = 0;                                                                  \
    while( s ){                                                               \
      n = s->np;                                                              \
      np_tracer[isp*nframes + frame] = n;                                     \
      ntraj_point_tracer[isp] += n;                                           \
      if ( n > 0 ){                                                           \
        p0 = s->p;                                                            \
        j = 0;                                                                \
        for ( p=p0; n; n--, j++, p++ ){                                       \
          dx0 = p->dx;                                                        \
          dy0 = p->dy;                                                        \
          dz0 = p->dz;                                                        \
          ii  = p->i;                                                         \
          ux  = p->ux;                                                        \
          uy  = p->uy;                                                        \
          uz  = p->uz;                                                        \
          f   = f0 + ii;                                                      \
          index += nvar;                                                      \
          CALC_TRACER_USER_DEFINED_DATA;                                      \
          pout[index + 0]  = (float) p->w;                                    \
          pout[index + 1]  = (float) tracer_x ;                               \
          pout[index + 2]  = (float) tracer_y ;                               \
          pout[index + 3]  = (float) tracer_z ;                               \
          pout[index + 4]  = ux;                                              \
          pout[index + 5]  = uy;                                              \
          pout[index + 6]  = uz;                                              \
          TRACER_USER_DEFINED_DATA;                                           \
        }                                                                     \
      }                                                                       \
      isp++;                                                                  \
      s = s->next;                                                            \
    }                                                                         \
} END_PRIMITIVE

//--------------------------------------------------------------
// This routine dumps buffered tracer information.

#define dump_buffered_tracer(nframes, fbase)                                  \
  BEGIN_PRIMITIVE{                                                            \
    static int psize = 0;                                                     \
    static int ntracer = global->Ntracer / nproc();                           \
    static int nvar = TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS;       \
    psize = 2*ntracer*nvar*nframes;                                           \
    char dname[256], fname[256];                                              \
    FileIO       fh;                                                          \
    FileIOStatus status;                                                      \
    species_t *s = global->tracers_list ;                                     \
    sprintf(dname, "%s/T.%d", fbase, step() );                                \
    dump_mkdir(dname);                                                        \
    int isp = 0;                                                              \
    while ( s ) {                                                             \
        sprintf(fname, "%s/%s_%s.%d", dname , s->name, "buffer", (int)rank());\
        status = fh.open(fname, io_write);                                    \
        if ( status == fail ) ERROR(("Could not open file %s", fname));       \
        WRITE(int,    nvar, fh);                                              \
        WRITE(float,  s->q, fh);                                              \
        WRITE(float,  s->m, fh);                                              \
        WRITE(int,    nvar, fh);                                              \
        WRITE(int, nframes, fh);                                              \
        WRITE(float,  1.0 , fh);                                              \
        WRITE(double, 1.0 , fh);                                              \
        fh.write(np_tracer + isp*nframes, nframes);                           \
        fh.write(t_tracer, nframes);                                          \
        fh.write(pout + psize*isp, ntraj_point_tracer[isp]*nvar);             \
        fh.close();                                                           \
        s = s->next;                                                          \
        isp++;                                                                \
    }                                                                         \
} END_PRIMITIVE

//--------------------------------------------------------------
// This routine dumps buffered tracer information for restart.

#define dump_buffered_tracer_restart(nframes, RESTART_TAG)                    \
  BEGIN_PRIMITIVE{                                                            \
    static int psize = 0;                                                     \
    static int ntracer = global->Ntracer / nproc();                           \
    static int nvar = TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS;       \
    species_t *s = global->tracers_list ;                                     \
    int numspecies = 0;                                                       \
    char fname[256];                                                          \
    FileIO fh;                                                                \
    FileIOStatus status;                                                      \
    numspecies = num_species(s);                                              \
    psize = 2*numspecies*ntracer*nvar*nframes;                                \
    sprintf(fname, "restart/restart%d/buffered_tracer%d.%d", RESTART_TAG,     \
            global->particle_tracing, (int)rank() ) ;                         \
    status = fh.open(fname, io_write);                                        \
    if ( status == fail ) ERROR(("Could not open file %s", fname));           \
    WRITE(int, numspecies, fh);                                               \
    WRITE(int, nframes, fh);                                                  \
    WRITE(int, psize, fh);                                                    \
    fh.write(np_tracer, nframes*numspecies);                                  \
    fh.write(t_tracer, nframes);                                              \
    fh.write(ntraj_point_tracer, numspecies);                                 \
    fh.write(pout, psize);                                                    \
    fh.close();                                                               \
} END_PRIMITIVE

//--------------------------------------------------------------
// This routine reads buffered tracer information at restart.

#define read_buffered_tracer_restart(RESTART_TAG)                             \
  BEGIN_PRIMITIVE{                                                            \
    int numspecies, nframes, psize;                                           \
    FileIO fh;                                                                \
    FileIOStatus status;                                                      \
    char fname[256];                                                          \
    sprintf(fname, "restart/restart%d/buffered_tracer%d.%d", RESTART_TAG,     \
            global->particle_tracing, (int)rank() ) ;                         \
    status = fh.open(fname, io_read);                                         \
    if ( status == fail ) ERROR(("Could not open file %s", fname));           \
    READ(int, numspecies, fh);                                                \
    READ(int, nframes, fh);                                                   \
    READ(int, psize, fh);                                                     \
    fh.read(np_tracer, nframes*numspecies);                                   \
    fh.read(t_tracer, nframes);                                               \
    fh.read(ntraj_point_tracer, numspecies);                                  \
    fh.read(pout, psize);                                                     \
    fh.close();                                                               \
} END_PRIMITIVE


/****************************************************************
 ******************* Copied from dumpmacros.h *******************
 ***************************************************************/

#ifndef WRITE
#define WRITE(type,value,fileIO) BEGIN_PRIMITIVE {                            \
    type __WRITE_tmp = (type)(value);                                         \
    fileIO.write( &__WRITE_tmp, 1 );                                          \
} END_PRIMITIVE
#endif

// Note: strlen does not include the terminating NULL

#ifndef WRITE_STRING
#define WRITE_STRING(string,fileIO) BEGIN_PRIMITIVE {                         \
    int __WRITE_STRING_len = 0;                                               \
    if( string!=NULL ) __WRITE_STRING_len = strlen(string);                   \
    fileIO.write( &__WRITE_STRING_len, 1 );                                   \
    if( __WRITE_STRING_len>0 ) {                                              \
        fileIO.write( string, __WRITE_STRING_len );                           \
    }                                                                         \
} END_PRIMITIVE
#endif

#ifndef READ
#define READ(type,value,fileIO) BEGIN_PRIMITIVE {                             \
    type __READ_tmp;                                                          \
    fileIO.read(&__READ_tmp, 1 );                                             \
    (value) = __READ_tmp;                                                     \
} END_PRIMITIVE
#endif

#endif // __TRACER_HXX__

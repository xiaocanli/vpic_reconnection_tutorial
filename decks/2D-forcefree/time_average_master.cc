/*
 * This piece of code computes time-averaged quantities
 * Assumptions:
 *  1. The stride is 1 for all directions.
 *  2. global->dis_nav is smaller than global->dis_interval/3.
 *  3. We only need to dump electric field and magnetic field.
 *  4. dis_nav is odd.
 *
 * Global variables (initialized here):
 *  global->dis_nav: number of steps to average over
 *  global->dis_interval: number of steps between outputs
 *  global->dis_iter: iteration count. 0 means we are not averaging at the moment
 *  global->dis_begin_int: the first time step of the interval
 *
 * The adjustable parameters of these, dis_interval and dis_nav, are now set by
 * the variables AVG_SPACING and AVG_TOTAL_STEPS in the beginning so they can
 * be more easily found and edited.
 *
 * The data are saved in the three intervals:
 *  1. j*dis_interval - dis_nav*3/2 <= step() <= j*dis_interval - dis_nav/2 - 1
 *  2. j*dis_interval - dis_nav/2 <= step() <= j*dis_interval + dis_nav/2
 *  3. j*dis_interval + dis_nav/2 + 1 <= step() <= j*dis_interval + dis_nav*3/2
 * where j is an integer. The output is assigned time index corresponding to
 * the center of the interval.
 *
 * The code is restart-aware.
 *
 */

{ // start the diagnostics

  static int dis_initialized = 0;   // the flag that signals initialization

  FileIO fileIO;
  FileIOStatus status;
  char fname[256];
  char fdir[256];
  int dim[3];

  const int nx = grid->nx;
  const int ny = grid->ny;
  const int nz = grid->nz;

  //PARAMETERS FOR AVERAGING INTERVAL
  const int AVG_SPACING = 1; //Spacing between averaging. 1 should mean every time the
  //program writes non-averaged data, 2 every other, etc. Better to be odd.
  const int AVG_TOTAL_STEPS = 21; //Total number of steps to average over (NOT dt's)

#ifdef TURBULENCE_MIXING_DIAGNOSTICS
  static hydro_t * ALIGNED(128) hydro_avg_etop;
  static hydro_t * ALIGNED(128) hydro_avg_ebot;
  static hydro_t * ALIGNED(128) hydro_avg_htop;
  static hydro_t * ALIGNED(128) hydro_avg_hbot;
#else
  static hydro_t * ALIGNED(128) hydro_avg_e;
  static hydro_t * ALIGNED(128) hydro_avg_h;
#endif
  static hydro_t * ALIGNED(128) hi;
  static hydro_t * RESTRICT ALIGNED(16) h0;
  static hydro_t * RESTRICT ALIGNED(16) havg;

  static emf_t * ALIGNED(128) emf_avg;
  static field_t * ALIGNED(128) fi;
  static field_t * RESTRICT ALIGNED(16) f0;
  static emf_t * RESTRICT ALIGNED(16) favg;

  if (!dis_initialized) { // initialization upon start or re-start
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
    MALLOC_ALIGNED( hydro_avg_etop, grid->nv, 128 );
    MALLOC_ALIGNED( hydro_avg_ebot, grid->nv, 128 );
    MALLOC_ALIGNED( hydro_avg_htop, grid->nv, 128 );
    MALLOC_ALIGNED( hydro_avg_hbot, grid->nv, 128 );
#else
    MALLOC_ALIGNED( hydro_avg_e, grid->nv, 128 );
    MALLOC_ALIGNED( hydro_avg_h, grid->nv, 128 );
#endif
    MALLOC_ALIGNED( emf_avg, grid->nv, 128 );

    if (step() == 0) {
#ifdef DUMP_WITH_HDF5
      sprintf(global->fields_dir, "fields-avg-hdf5");
      sprintf(global->hydro_dir, "hydro-avg-hdf5");
      sprintf(global->restart_avg_dir, "restart-avg-hdf5");
#else
      dump_mkdir("fields-avg");
      dump_mkdir("hydro-avg");
      dump_mkdir("restart-avg");
      sprintf(global->fields_dir, "fields-avg/%d", NUMFOLD);
      sprintf(global->hydro_dir, "hydro-avg/%d", NUMFOLD);
      sprintf(global->restart_avg_dir, "restart-avg/%d", NUMFOLD);
#endif
      dump_mkdir(global->fields_dir);
      dump_mkdir(global->hydro_dir);
      dump_mkdir(global->restart_avg_dir);

      global->dis_iter = 0;      // initialize iteration count
      global->dis_begin_int = 0; // initialize the start of the averaging interval
      global->dis_nav = AVG_TOTAL_STEPS;
      global->dis_interval = AVG_SPACING*global->fields_interval;
      global->dis_begin_int = global->dis_interval - (global->dis_nav*3)/2;  // Skip t=0 averaging
    }

    // do we need to restart
    if (global->dis_iter>0) {
#ifdef DUMP_WITH_HDF5
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
      READ_TIME_AVERAGE_HYDRO_HDF5("electronTop", global->restart_avg_dir,
          hydro_avg_etop, step() - 1)
      READ_TIME_AVERAGE_HYDRO_HDF5("electronBot", global->restart_avg_dir,
          hydro_avg_ebot, step() - 1)
      READ_TIME_AVERAGE_HYDRO_HDF5("ionTop", global->restart_avg_dir,
          hydro_avg_htop, step() - 1)
      READ_TIME_AVERAGE_HYDRO_HDF5("ionBot", global->restart_avg_dir,
          hydro_avg_hbot, step() - 1)
#else
      READ_TIME_AVERAGE_HYDRO_HDF5("electron", global->restart_avg_dir,
          hydro_avg_e, step() - 1)
      READ_TIME_AVERAGE_HYDRO_HDF5("ion", global->restart_avg_dir,
          hydro_avg_h, step() - 1)
#endif
      READ_TIME_AVERAGE_FIELDS_HDF5(global->restart_avg_dir, emf_avg, step() - 1)
#else
      sprintf(fname, "%s/%s.%d", global->restart_avg_dir, "restart-avg", rank());
      status = fileIO.open(fname, io_read);
      if(status == fail) ERROR(("Failed opening file: %s", fname));
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
      fileIO.read(hydro_avg_etop, (nx+2)*(ny+2)*(nz+2));
      fileIO.read(hydro_avg_ebot, (nx+2)*(ny+2)*(nz+2));
      fileIO.read(hydro_avg_htop, (nx+2)*(ny+2)*(nz+2));
      fileIO.read(hydro_avg_hbot, (nx+2)*(ny+2)*(nz+2));
#else
      fileIO.read(hydro_avg_e, (nx+2)*(ny+2)*(nz+2));
      fileIO.read(hydro_avg_h, (nx+2)*(ny+2)*(nz+2));
#endif
      fileIO.read(emf_avg, (nx+2)*(ny+2)*(nz+2));
      if( fileIO.close() ) ERROR(( "File close failed on restart dump of averaged data!!!" ));
#endif
    } // end of file read

    dis_initialized = 1;

  } // end of initialization

  // This way of computing remainder has issues for some reason - Bill
  // int r=remainder(step, global->dis_interval);

  int rstep = (step() + (global->dis_nav*3)/2)%global->dis_interval;

  species_t *sp;

  // check that we are inside the averaging interval
  if (( rstep >=0) && (rstep < (global->dis_nav * 3) && step()>= global->dis_begin_int )) {
    if (global->dis_iter == 0) { //  this is the first iteration
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
      CLEAR(hydro_avg_etop, grid->nv);
      CLEAR(hydro_avg_ebot, grid->nv);
      CLEAR(hydro_avg_htop, grid->nv);
      CLEAR(hydro_avg_hbot, grid->nv);
#else
      CLEAR(hydro_avg_e, grid->nv);
      CLEAR(hydro_avg_h, grid->nv);
#endif
      CLEAR(emf_avg, grid->nv);
      global->dis_begin_int = step();    // mark the beginning of the interval
    }

#ifdef TURBULENCE_MIXING_DIAGNOSTICS
    sp = find_species_name("electronTop", species_list);
    clear_hydro_array( hydro_array );
    accumulate_hydro_p( hydro_array, sp, interpolator_array );
    synchronize_hydro_array( hydro_array );
    hi = hydro_array->h;
    int x, y, z;
    LOOP_YZ_FACE {
      havg = &HYDRO_AVG_ETOP(0, y, z);
      h0 = &HYDRO_STEP(0, y, z);
      HYDRO_ALONG_X;
    }
    sp = find_species_name("electronBot", species_list);
    clear_hydro_array( hydro_array );
    accumulate_hydro_p( hydro_array, sp, interpolator_array );
    synchronize_hydro_array( hydro_array );
    hi = hydro_array->h;
    LOOP_YZ_FACE {
      havg = &HYDRO_AVG_EBOT(0, y, z);
      h0 = &HYDRO_STEP(0, y, z);
      HYDRO_ALONG_X;
    }
#else
    sp = find_species_name("electron", species_list);
    clear_hydro_array( hydro_array );
    accumulate_hydro_p( hydro_array, sp, interpolator_array );
    synchronize_hydro_array( hydro_array );
    hi = hydro_array->h;
    int x, y, z;
    LOOP_YZ_FACE {
      havg = &HYDRO_AVG_E(0, y, z);
      h0 = &HYDRO_STEP(0, y, z);
      HYDRO_ALONG_X;
    }
#endif

#ifdef TURBULENCE_MIXING_DIAGNOSTICS
    sp = find_species_name("ionTop", species_list);
    clear_hydro_array( hydro_array );
    accumulate_hydro_p( hydro_array, sp, interpolator_array );
    synchronize_hydro_array( hydro_array );
    hi = hydro_array->h;
    LOOP_YZ_FACE {
      havg = &HYDRO_AVG_HTOP(0, y, z);
      h0 = &HYDRO_STEP(0, y, z);
      HYDRO_ALONG_X;
    }
    sp = find_species_name("ionBot", species_list);
    clear_hydro_array( hydro_array );
    accumulate_hydro_p( hydro_array, sp, interpolator_array );
    synchronize_hydro_array( hydro_array );
    hi = hydro_array->h;
    LOOP_YZ_FACE {
      havg = &HYDRO_AVG_HBOT(0, y, z);
      h0 = &HYDRO_STEP(0, y, z);
      HYDRO_ALONG_X;
    }
#else
    sp = find_species_name("ion", species_list);
    clear_hydro_array( hydro_array );
    accumulate_hydro_p( hydro_array, sp, interpolator_array );
    synchronize_hydro_array( hydro_array );
    hi = hydro_array->h;
    LOOP_YZ_FACE {
      havg = &HYDRO_AVG_H(0, y, z);
      h0 = &HYDRO_STEP(0, y, z);
      HYDRO_ALONG_X;
    }
#endif

    fi = field_array->f;
    LOOP_YZ_FACE {
      favg = &EMF_AVG(0, y, z);
      f0 = &FIELD_STEP(0, y, z);
      EMF_ALONG_X;
    }

    global->dis_iter++;  // iteration count
    global->dis_iter = global->dis_iter%global->dis_nav; // reset

    // write files if it's the last iteration OR we need to save restart
    bool write_restart1 = (step()>0 && should_dump(restart));
    bool write_restart2 = (step()>0 && (global->quota_check_interval &&
                                        (step()%global->quota_check_interval)==0 &&
                                        uptime() > global->quota_sec));
    bool write_restart = (write_restart1 || write_restart2);
    static float isteps;
    if (write_restart) { // write the restart file
#ifdef DUMP_WITH_HDF5
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
      DUMP_TIME_AVERAGE_HYDRO_HDF5("electronTop", global->restart_avg_dir,
          hydro_avg_etop, step())
      DUMP_TIME_AVERAGE_HYDRO_HDF5("electronBot", global->restart_avg_dir,
          hydro_avg_ebot, step())
      DUMP_TIME_AVERAGE_HYDRO_HDF5("ionTop", global->restart_avg_dir,
          hydro_avg_htop, step())
      DUMP_TIME_AVERAGE_HYDRO_HDF5("ionBot", global->restart_avg_dir,
          hydro_avg_hbot, step())
#else
      DUMP_TIME_AVERAGE_HYDRO_HDF5("electron", global->restart_avg_dir,
          hydro_avg_e, step())
      DUMP_TIME_AVERAGE_HYDRO_HDF5("ion", global->restart_avg_dir,
          hydro_avg_h, step())
#endif
      DUMP_TIME_AVERAGE_FIELDS_HDF5(global->restart_avg_dir, emf_avg, step())
#else
      sprintf(fname, "%s/%s.%d", global->restart_avg_dir, "restart-avg", rank());
      status = fileIO.open(fname, io_write);
      if(status == fail) ERROR(("Failed opening file: %s", fname));
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
      fileIO.write(hydro_avg_etop, (nx+2)*(ny+2)*(nz+2));
      fileIO.write(hydro_avg_ebot, (nx+2)*(ny+2)*(nz+2));
      fileIO.write(hydro_avg_htop, (nx+2)*(ny+2)*(nz+2));
      fileIO.write(hydro_avg_hbot, (nx+2)*(ny+2)*(nz+2));
#else
      fileIO.write(hydro_avg_e, (nx+2)*(ny+2)*(nz+2));
      fileIO.write(hydro_avg_h, (nx+2)*(ny+2)*(nz+2));
#endif
      fileIO.write(emf_avg, (nx+2)*(ny+2)*(nz+2));
      if( fileIO.close() ) ERROR(( "File close failed on restart dump of averaged data!!!" ));
#endif

      if (write_restart2) { // release the resources
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
        if (!hydro_avg_etop) FREE_ALIGNED( hydro_avg_etop );
        if (!hydro_avg_ebot) FREE_ALIGNED( hydro_avg_ebot );
        if (!hydro_avg_htop) FREE_ALIGNED( hydro_avg_htop );
        if (!hydro_avg_hbot) FREE_ALIGNED( hydro_avg_hbot );
#else
        if (!hydro_avg_e) FREE_ALIGNED( hydro_avg_e );
        if (!hydro_avg_h) FREE_ALIGNED( hydro_avg_h );
#endif
        if (!emf_avg)     FREE_ALIGNED( emf_avg );
      }
    } else { // regular write
      if (global->dis_iter == 0) {
        isteps = 1.0 / global->dis_nav;
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
        LOOP_YZ_FACE {
          havg = &HYDRO_AVG_ETOP(0, y, z);
          AVG_HYDRO;
        }
        LOOP_YZ_FACE {
          havg = &HYDRO_AVG_EBOT(0, y, z);
          AVG_HYDRO;
        }
        LOOP_YZ_FACE {
          havg = &HYDRO_AVG_HTOP(0, y, z);
          AVG_HYDRO;
        }
        LOOP_YZ_FACE {
          havg = &HYDRO_AVG_HBOT(0, y, z);
          AVG_HYDRO;
        }
#else
        LOOP_YZ_FACE {
          havg = &HYDRO_AVG_E(0, y, z);
          AVG_HYDRO;
        }
        LOOP_YZ_FACE {
          havg = &HYDRO_AVG_H(0, y, z);
          AVG_HYDRO;
        }
#endif
        LOOP_YZ_FACE {
          favg = &EMF_AVG(0, y, z);
          AVG_EMF;
        }

        int step_dump = step() - global->dis_nav/2;

#ifdef DUMP_WITH_HDF5
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
        DUMP_TIME_AVERAGE_HYDRO_HDF5("electronTop", global->hydro_dir, hydro_avg_etop, step_dump)
        DUMP_TIME_AVERAGE_HYDRO_HDF5("electronBot", global->hydro_dir, hydro_avg_ebot, step_dump)
        DUMP_TIME_AVERAGE_HYDRO_HDF5("ionTop", global->hydro_dir, hydro_avg_htop, step_dump)
        DUMP_TIME_AVERAGE_HYDRO_HDF5("ionBot", global->hydro_dir, hydro_avg_hbot, step_dump)
#else
        DUMP_TIME_AVERAGE_HYDRO_HDF5("electron", global->hydro_dir, hydro_avg_e, step_dump)
        DUMP_TIME_AVERAGE_HYDRO_HDF5("ion", global->hydro_dir, hydro_avg_h, step_dump)
#endif
        DUMP_TIME_AVERAGE_FIELDS_HDF5(global->fields_dir, emf_avg, step_dump)
#else
        size_t numvars;
        size_t * varlist;
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
        DUMP_HYDRO_AVG("etophydro", "electronTop", e, HYDRO_AVG_ETOP);
        DUMP_HYDRO_AVG("ebothydro", "electronBot", e, HYDRO_AVG_EBOT);
        DUMP_HYDRO_AVG("Htophydro", "ionTop", H, HYDRO_AVG_HTOP);
        DUMP_HYDRO_AVG("Hbothydro", "ionBot", H, HYDRO_AVG_HBOT);
#else
        DUMP_HYDRO_AVG("ehydro", "electron", e, HYDRO_AVG_E);
        DUMP_HYDRO_AVG("Hhydro", "ion", H, HYDRO_AVG_H);
#endif
        DUMP_EMF_AVG("fields", EMF_AVG);
#endif

        sim_log("Time-averaging diagnostic");

        if (rank() == 0) {
          MESSAGE(("Time Averaging ---> start=%d   end=%d  steps=%d\n",
                   global->dis_begin_int,step(),1+step()-global->dis_begin_int));
        }
      }
    }
  } // end of inside the averaging interval

  // The end of the simulation
  if ( step() == num_step ) {
      // release the resources
#ifdef TURBULENCE_MIXING_DIAGNOSTICS
      if (!hydro_avg_etop) FREE_ALIGNED( hydro_avg_etop );
      if (!hydro_avg_ebot) FREE_ALIGNED( hydro_avg_ebot );
      if (!hydro_avg_htop) FREE_ALIGNED( hydro_avg_htop );
      if (!hydro_avg_hbot) FREE_ALIGNED( hydro_avg_hbot );
#else
      if (!hydro_avg_e) FREE_ALIGNED( hydro_avg_e );
      if (!hydro_avg_h) FREE_ALIGNED( hydro_avg_h );
#endif
      if (!emf_avg)     FREE_ALIGNED( emf_avg );
  }
} // end of time-averaging diagnostic

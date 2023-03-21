/*
 * This piece of code computes time-averaged quantities
 * Assumptions:
 *  1. The stride is 1 for all directions.
 *  2. global->dis_nav is smaller than global->dis_iter*2/3.
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

#define HYDRO_STEP(x,y,z) hi[LOCAL_CELL_ID(x,y,z)]
#define FIELD_STEP(x,y,z) fi[LOCAL_CELL_ID(x,y,z)]
#define HYDRO_AVG_E(x,y,z) hydro_avg_e[LOCAL_CELL_ID(x,y,z)]
#define HYDRO_AVG_H(x,y,z) hydro_avg_h[LOCAL_CELL_ID(x,y,z)]
#define EMF_AVG(x,y,z) emf_avg[LOCAL_CELL_ID(x,y,z)]

// loop over yz face
#define LOOP_YZ_FACE                        \
    for (z = 0; z < nz + 2; z++)            \
        for (y = 0; y < ny + 2; y++)

#define HYDRO_ALONG_X                       \
    for (x = 0; x < nx + 2; x++) {          \
        havg->jx += h0->jx;                 \
        havg->jy += h0->jy;                 \
        havg->jz += h0->jz;                 \
        havg->rho += h0->rho;               \
        havg->px += h0->px;                 \
        havg->py += h0->py;                 \
        havg->pz += h0->pz;                 \
        havg->ke += h0->ke;                 \
        havg->txx += h0->txx;               \
        havg->tyy += h0->tyy;               \
        havg->tzz += h0->tzz;               \
        havg->tyz += h0->tyz;               \
        havg->tzx += h0->tzx;               \
        havg->txy += h0->txy;               \
        havg++;                             \
        h0++;                               \
    }

#define EMF_ALONG_X                         \
    for (x = 0; x < nx + 2; x++) {          \
        favg->ex += f0->ex;                 \
        favg->ey += f0->ey;                 \
        favg->ez += f0->ez;                 \
        favg->cbx += f0->cbx;               \
        favg->cby += f0->cby;               \
        favg->cbz += f0->cbz;               \
        favg++;                             \
        f0++;                               \
    }

#define AVG_HYDRO                           \
    for (x = 0; x < nx + 2; x++) {          \
        havg->jx *= isteps;                 \
        havg->jy *= isteps;                 \
        havg->jz *= isteps;                 \
        havg->rho *= isteps;                \
        havg->px *= isteps;                 \
        havg->py *= isteps;                 \
        havg->pz *= isteps;                 \
        havg->ke *= isteps;                 \
        havg->txx *= isteps;                \
        havg->tyy *= isteps;                \
        havg->tzz *= isteps;                \
        havg->tyz *= isteps;                \
        havg->tzx *= isteps;                \
        havg->txy *= isteps;                \
        havg++;                             \
    }

#define AVG_EMF                             \
    for (x = 0; x < nx + 2; x++) {          \
        favg->ex *= isteps;                 \
        favg->ey *= isteps;                 \
        favg->ez *= isteps;                 \
        favg->cbx *= isteps;                \
        favg->cby *= isteps;                \
        favg->cbz *= isteps;                \
        favg++;                             \
    }

#define DUMP_HYDRO_AVG(filename, species, sname, ha)                            \
    sprintf(fdir, "%s/T.%d", global->hydro_dir, step_dump);                     \
    dump_mkdir(fdir);                                                           \
    sprintf(fname, "%s/%s.%ld.%d", fdir, filename, step_dump, rank());          \
    status = fileIO.open(fname, io_write);                                      \
    if(status == fail) ERROR(("Failed opening file: %s", fname));               \
    sp = find_species_name(species, species_list);                              \
    /* Time step in the header is not correct. */                               \
    /* It should not be used when doing analysis */                             \
    WRITE_HEADER_V0(2, sp->id, sp->q/sp->m, fileIO);                            \
    dim[0] = nx + 2;                                                            \
    dim[1] = ny + 2;                                                            \
    dim[2] = nz + 2;                                                            \
    WRITE_ARRAY_HEADER(hydro_array->h, 3, dim, fileIO);                         \
    /* Create a variable list of hydro values to output */                      \
    numvars = std::min(global->h##sname##dParams.output_vars.bitsum(),          \
                              total_hydro_variables);                           \
    varlist = new size_t[numvars];                                              \
    for(size_t i(0), c(0); i<total_hydro_variables; i++)                        \
        if( global->h##sname##dParams.output_vars.bitset(i) ) varlist[c++] = i; \
                                                                                \
    for(size_t v(0); v<numvars; v++)                                            \
    for(size_t k(0); k<nz+2; k++)                                               \
    for(size_t j(0); j<ny+2; j++)                                               \
    for(size_t i(0); i<nx+2; i++) {                                             \
        const uint32_t * href = reinterpret_cast<uint32_t *>(&ha(i,j,k));       \
        fileIO.write(&href[varlist[v]], 1);                                     \
    }                                                                           \
                                                                                \
    if( fileIO.close() ) ERROR(( "File close failed on hydro dump!!!" ));       \
    delete[] varlist;

#define DUMP_EMF_AVG(filename, fa)                                              \
    sprintf(fdir, "%s/T.%d", global->fields_dir, step_dump);                    \
    dump_mkdir(fdir);                                                           \
    sprintf(fname, "%s/%s.%ld.%d", fdir, filename, step_dump, rank());          \
    status = fileIO.open(fname, io_write);                                      \
    if(status == fail) ERROR(("Failed opening file: %s", fname));               \
    /* Time step in the header is not correct. */                               \
    /* It should not be used when doing analysis */                             \
    WRITE_HEADER_V0(1, -1, 0, fileIO);                                          \
    dim[0] = nx + 2;                                                            \
    dim[1] = ny + 2;                                                            \
    dim[2] = nz + 2;                                                            \
    WRITE_ARRAY_HEADER(field_array->f, 3, dim, fileIO);                         \
    for(size_t v(0); v<numvars; v++)                                            \
    for(size_t k(0); k<nz+2; k++)                                               \
    for(size_t j(0); j<ny+2; j++)                                               \
    for(size_t i(0); i<nx+2; i++) {                                             \
        const uint32_t * fref = reinterpret_cast<uint32_t *>(&fa(i,j,k));       \
        fileIO.write(&fref[v], 1);                                              \
    }                                                                           \
                                                                                \
    if( fileIO.close() ) ERROR(( "File close failed on field dump!!!" ));


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

    static hydro_t * ALIGNED(128) hydro_avg_e;
    static hydro_t * ALIGNED(128) hydro_avg_h;
    static hydro_t * ALIGNED(128) hi;
    static hydro_t * RESTRICT ALIGNED(16) h0;
    static hydro_t * RESTRICT ALIGNED(16) havg;

    static emf_t * ALIGNED(128) emf_avg;
    static field_t * ALIGNED(128) fi;
    static field_t * RESTRICT ALIGNED(16) f0;
    static emf_t * RESTRICT ALIGNED(16) favg;

    if (!dis_initialized) { // initialization upon start or re-start
        MALLOC_ALIGNED( hydro_avg_e, grid->nv, 128 );
        MALLOC_ALIGNED( hydro_avg_h, grid->nv, 128 );
        MALLOC_ALIGNED( emf_avg, grid->nv, 128 );

        if (step() == 0) {
            dump_mkdir("fields-avg");
            dump_mkdir("hydro-avg");
            dump_mkdir("restart-avg");
            sprintf(global->fields_dir, "fields-avg/%d", NUMFOLD);
            sprintf(global->hydro_dir, "hydro-avg/%d", NUMFOLD);
            sprintf(global->restart_avg_dir, "restart-avg/%d", NUMFOLD);
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
            sprintf(fname, "%s/%s.%d", global->restart_avg_dir, "restart-avg", rank());
            status = fileIO.open(fname, io_read);
            if(status == fail) ERROR(("Failed opening file: %s", fname));
            fileIO.read(hydro_avg_e, (nx+2)*(ny+2)*(nz+2));
            fileIO.read(hydro_avg_h, (nx+2)*(ny+2)*(nz+2));
            fileIO.read(emf_avg, (nx+2)*(ny+2)*(nz+2));
            if( fileIO.close() ) ERROR(( "File close failed on restart dump of averaged data!!!" ));
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
            CLEAR(hydro_avg_e, grid->nv);
            CLEAR(hydro_avg_h, grid->nv);
            CLEAR(emf_avg, grid->nv);
            global->dis_begin_int = step();    // mark the beginning of the interval
        }

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
            sprintf(fname, "%s/%s.%d", global->restart_avg_dir, "restart-avg", rank());
            status = fileIO.open(fname, io_write);
            if(status == fail) ERROR(("Failed opening file: %s", fname));
            fileIO.write(hydro_avg_e, (nx+2)*(ny+2)*(nz+2));
            fileIO.write(hydro_avg_h, (nx+2)*(ny+2)*(nz+2));
            fileIO.write(emf_avg, (nx+2)*(ny+2)*(nz+2));
            if( fileIO.close() ) ERROR(( "File close failed on restart dump of averaged data!!!" ));

            if (write_restart2) { // release the resources
                if (!hydro_avg_e) FREE_ALIGNED( hydro_avg_e );
                if (!hydro_avg_h) FREE_ALIGNED( hydro_avg_h );
                if (!emf_avg)     FREE_ALIGNED( emf_avg );
            }
        } else { // regular write
            if (global->dis_iter == 0) {
                isteps = 1.0 / global->dis_nav;
                LOOP_YZ_FACE {
                    havg = &HYDRO_AVG_E(0, y, z);
                    AVG_HYDRO;
                }
                LOOP_YZ_FACE {
                    havg = &HYDRO_AVG_H(0, y, z);
                    AVG_HYDRO;
                }
                LOOP_YZ_FACE {
                    favg = &EMF_AVG(0, y, z);
                    AVG_EMF;
                }

                int step_dump = step() - global->dis_nav/2;
                size_t numvars;
                size_t * varlist;
                DUMP_HYDRO_AVG("ehydro", "electron", e, HYDRO_AVG_E);
                DUMP_HYDRO_AVG("Hhydro", "ion", H, HYDRO_AVG_H);
                DUMP_EMF_AVG("fields", EMF_AVG);

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
        if (!hydro_avg_e) FREE_ALIGNED( hydro_avg_e );
        if (!hydro_avg_h) FREE_ALIGNED( hydro_avg_h );
        if (!emf_avg)     FREE_ALIGNED( emf_avg );
    }
} // end of time-averaging diagnostic

#undef HYDRO_STEP
#undef HYDRO_AVG_E
#undef HYDRO_AVG_H
#undef FIELD_STEP
#undef EMF_AVG
#undef LOOP_YZ_FACE
#undef HYDRO_ALONG_X
#undef EMF_ALONG_X
#undef AVG_HYDRO
#undef AVG_EMF
#undef DUMP_HYDRO_AVG
#undef DUMP_EMF_AVG

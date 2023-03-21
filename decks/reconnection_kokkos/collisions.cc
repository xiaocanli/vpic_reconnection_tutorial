/* Datatypes required by collision algorithm. */

/* Key datatype, used for shuffling particles */
typedef struct _pkey_t {
  int i;           /* Particle array index */
  float r;         /* Random variable used to shuffle+sort keys */
} pkey_t;

typedef struct _cdata_t {
  species_t *s;    /* Particle species */
  float     cvar;  /* Base variance of polar angle--see T&A paper */
  float     m;     /* Particle mass / electron mass */
  float     z;     /* Particle core charge */
} cdata_t;

/* Heapsort of key array; modified from Numerical Recipes in C
   2nd Ed., Press et al.  Note: N.R.C. arrays all are indexed
   starting with 1, so KARR needs to be offset by one element. */

#define COPY_KEY(K1,K2)   (K1).i=(K2).i, (K1).r=(K2).r
#define HEAPSORT_KEYS(N,KARR)                   \
  if ( (N)>1 )                                  \
    do {                                        \
      unsigned long n;                          \
      pkey_t *ka;                               \
      unsigned long i,ir,j,l;                   \
      pkey_t rka;                               \
                                                \
      n=(N);                                    \
      ka=(KARR);                                \
      l=(n >> 1)+1;                             \
      ir=n;                                     \
      for (;;) {                                \
        if ( l>1 ) {                            \
          --l;                                  \
          COPY_KEY( rka, ka[l] );               \
        } else {                                \
          COPY_KEY( rka, ka[ir] );              \
          COPY_KEY( ka[ir], ka[1] );            \
          if ( --ir==1 ) {                      \
            COPY_KEY( ka[1], rka );             \
            break;                              \
          }                                     \
        }                                       \
        i=l;                                    \
        j=l+l;                                  \
        while ( j<=ir ) {                       \
          if ( j<ir && ka[j].r<ka[j+1].r ) j++; \
          if ( rka.r<ka[j].r ) {                \
            COPY_KEY( ka[i], ka[j] );           \
            i=j;                                \
            j <<= 1;                            \
          } else break;                         \
        }                                       \
        COPY_KEY( ka[i], rka );                 \
      }                                         \
    } while (0)

/* Takizuka and Abe model for 2-particle elastic scattering from JCP 25, 205
   (1977).
   I, J:    particle indices
   P1, P2:  particle array heads
   M1, M2:  relative masses of particles
   var:     variance of collision operator */

#define SCATTER_PARTICLES(I,J,P1,P2,M1,M2,SQRT_VAR)                     \
  do {                                                                  \
    float ux, uy, uz, uperp, u, dux, duy, duz, mratio1, mratio2;        \
    float sin_theta, one_m_cos_theta, sin_phi, cos_phi, delta, phi;     \
                                                                        \
    mratio1=(float)M2/(float)(M1+M2);                                   \
    mratio2=(float)M1/(float)(M1+M2);                                   \
    ux=(P1)[I].ux-(P2)[J].ux;                                           \
    uy=(P1)[I].uy-(P2)[J].uy;                                           \
    uz=(P1)[I].uz-(P2)[J].uz;                                           \
    uperp=sqrt(ux*ux+uy*uy);                                            \
    u=sqrt(uperp*uperp+uz*uz);                                          \
    delta=normal( rng(0), 0, SQRT_VAR)*pow(u,-1.5);                     \
    sin_theta=2*delta/(1+delta*delta);                                  \
    one_m_cos_theta=sin_theta*delta;                                    \
    phi=2*M_PI*uniform( rng(0), 0, 1 );                                 \
    sin_phi=sin(phi);                                                   \
    cos_phi=cos(phi);                                                   \
    if ( uperp>0 ) {  /* General case */                                \
      dux=((ux*uz*sin_theta*cos_phi)-uy*u*sin_theta*sin_phi)/uperp      \
        - ux*one_m_cos_theta;                                           \
      duy=((uy*uz*sin_theta*cos_phi)+ux*u*sin_theta*sin_phi)/uperp      \
        - uy*one_m_cos_theta;                                           \
      duz=-uperp*sin_theta*cos_phi-uz*one_m_cos_theta;                  \
    } else { /* Handle purely z-directed difference vectors separately */ \
      dux=u*sin_theta*cos_phi;                                          \
      duy=u*sin_theta*sin_phi;                                          \
      duz=-u*one_m_cos_theta;                                           \
    }                                                                   \
    (P1)[I].ux+=mratio1*dux;                                            \
    (P1)[I].uy+=mratio1*duy;                                            \
    (P1)[I].uz+=mratio1*duz;                                            \
    (P2)[J].ux-=mratio2*dux;                                            \
    (P2)[J].uy-=mratio2*duy;                                            \
    (P2)[J].uz-=mratio2*duz;                                            \
  } while (0);


begin_particle_collisions {

  // Put collision model here since particle injection occurs after the particle
  // sort (collision models require sorted particle array).

  // Eventually this needs to be integrated into the VPIC source tree, but for
  // now leave this in the input deck.

  int do_collisions = global->ee_collisions && global->ei_collisions
    && global->ii_collisions ;

  if ( do_collisions &&  step() > 0 && step()%global->tstep_coll==0 ) {
    species_t *s1, *s2;
    int icell, i, j, npar1, npar2, isp1, isp2;
    static int initted=0, nsp;
    static size_t karr_max_length;
    static cdata_t cd[2 /* =nsp */];
    static pkey_t *karr1, *karr2;
    float var, sqrt_var, sqrt_half_var, density;
    particle_t *phead1, *phead2;

    if ( initted==0 ) {

      MESSAGE(("Initializing collisional operator\n"));

      karr_max_length = (size_t)global->nppc_max;
      karr1=(pkey_t *)malloc( karr_max_length*sizeof(pkey_t) );
      if ( !karr1 )
        ERROR(("Could not allocate key array in collision handler."));
      karr2=(pkey_t *)malloc( karr_max_length*sizeof(pkey_t) );
      if ( !karr2 )
        ERROR(("Could not allocate key array in collision handler."));

      isp1=0;
      cd[isp1].s    = find_species( "electron" );
      cd[isp1].m    = 1;
      cd[isp1].z    = -1;
      cd[isp1].cvar = // Fix Brian's stupid bug!
        global->cvar*pow(cd[isp1].z,4)/pow(cd[isp1].m,2)*4;
      MESSAGE(("cvar = %e\n", cd[isp1].cvar));

      ++isp1;
      cd[isp1].s    = find_species( "ion" );
      cd[isp1].m    = global->mi_me;  /* in m_e */
      cd[isp1].z    = global->Z;     /* in e */
      cd[isp1].cvar = // Fix Brian's stupid bug!
        global->cvar*pow(cd[isp1].z,4)/pow(cd[isp1].m,2)*4;
      MESSAGE(("cvar = %e\n", cd[isp1].cvar));

      nsp=isp1+1;
      initted=1;
    }

    // Self-scattering among species
    for ( isp1=0; isp1<nsp; ++isp1 ) {
      // MESSAGE(("Self-Scatter"));
      s1=cd[isp1].s;

      for ( icell=((grid->nx+2)*(grid->ny+2)*(grid->nz+2))-2;
            icell>=0; --icell ) {
        npar1  = s1->partition[icell+1]-s1->partition[icell];
        density = npar1*(global->nfac);
        //if ( density > 0.98 )
        //  MESSAGE(("npar1, density = %i  %e\n",npar1,density));

        var=cd[isp1].cvar*density;
        sqrt_var=sqrt(var);
        sqrt_half_var=sqrt(var/2);

        phead1 = s1->p;
        for ( i=0; i<npar1; ++i ) {
          karr1[i].i=i+s1->partition[icell];
          karr1[i].r=uniform( rng(0), 0, 1 );
        }
        HEAPSORT_KEYS( npar1, karr1-1 );
        if ( npar1>4 ) {
          for ( i=0; i<npar1-4; i+=2 )
            SCATTER_PARTICLES( karr1[i  ].i, karr1[i+1].i, phead1, phead1,
                               1.0, 1.0, sqrt_var );
          if ( i==npar1-4 ) {
            SCATTER_PARTICLES( karr1[i  ].i, karr1[i+1].i, phead1, phead1,
                               1.0, 1.0, sqrt_var );
            SCATTER_PARTICLES( karr1[i+2].i, karr1[i+3].i, phead1, phead1,
                               1.0, 1.0, sqrt_var  );
          } else {
            SCATTER_PARTICLES( karr1[i  ].i, karr1[i+1].i, phead1, phead1,
                               1.0, 1.0, sqrt_half_var );
            SCATTER_PARTICLES( karr1[i  ].i, karr1[i+2].i, phead1, phead1,
                               1.0, 1.0, sqrt_half_var );
            SCATTER_PARTICLES( karr1[i+1].i, karr1[i+2].i, phead1, phead1,
                               1.0, 1.0, sqrt_half_var );
          }
        }
      }
    }

    // Cross-species scattering
    if ( 1 ) {
      //   MESSAGE(("Cross-Scatter"));
      for ( isp1=1; isp1<nsp; ++isp1 ) {
#if 1
        s1=cd[isp1].s;
        for ( isp2=0; isp2<isp1; ++isp2 ) {
          float m1=cd[isp1].m, m2=cd[isp2].m;
          float rmass=m1*m2/(m1+m2);

          s2=cd[isp2].s;

          for ( icell=((grid->nx+2)*(grid->ny+2)*(grid->nz+2))-2; icell>=0;
                --icell ) {
            npar1  = s1->partition[icell+1]-s1->partition[icell];
            npar2  = s2->partition[icell+1]-s2->partition[icell];

            density = npar1*(global->nfac);
            var=density*(global->cvar)*pow(cd[isp1].z,2)
              *pow(cd[isp2].z,2)/(rmass*rmass);
            sqrt_var=sqrt(var);

            phead1 = s1->p;
            phead2 = s2->p;

            //if ( npar1>(int)karr_max_length || npar2>(int)karr_max_length )
            // MESSAGE(("Error: too many keys needed!  %d %d (%d)\n", npar1,
            // npar2, (int)karr_max_length ));

            /* Generate key arrays and sort */
            for ( i=0; i<npar1; ++i ) {
              karr1[i].i=i+s1->partition[icell];
              karr1[i].r=uniform( rng(0), 0, 1 );
            }
            HEAPSORT_KEYS( npar1, karr1-1 );
            for ( i=0; i<npar2; ++i ) {
              karr2[i].i=i+s2->partition[icell];
              karr2[i].r=uniform( rng(0), 0, 1 );
            }
            HEAPSORT_KEYS( npar2, karr2-1 );

            // Call this with npar1 >= npar2
#if 1
#define CROSS_COLLISION_LOOP(NPAR1,NPAR2,P1,P2,K1,K2,M1,M2,SQRT_VAR)    \
            do {                                                        \
              if ( NPAR1>0 && NPAR2>0 ) {                               \
                int ii, iimax, index1=-1, index2=-1;                    \
                                                                        \
                ii=(NPAR1)/(NPAR2);                                     \
                iimax=(NPAR1)-(NPAR2)*ii;                               \
                for ( i=0; i<iimax; ++i) {                              \
                  ++index2;                                             \
                  for ( j=0; j<ii+1; ++j ) {                            \
                    ++index1;                                           \
                    SCATTER_PARTICLES( K1[index1].i, K2[index2].i, P1, P2, \
                                       M1, M2, (SQRT_VAR)/sqrt(ii+1) ); \
                  }                                                     \
                }                                                       \
                for ( ; i<(NPAR2); ++i ) {                              \
                  ++index2;                                             \
                  for ( j=0; j<ii; ++j ) {                              \
                    ++index1;                                           \
                    SCATTER_PARTICLES( K1[index1].i, K2[index2].i, P1, P2, \
                                       M1, M2, (SQRT_VAR)/sqrt(ii) );   \
                  }                                                     \
                }                                                       \
              }                                                         \
            } while (0)

            if ( npar1>npar2 )
              CROSS_COLLISION_LOOP( npar1, npar2, phead1,
                                    phead2, karr1, karr2, m1, m2, sqrt_var );
            else
              CROSS_COLLISION_LOOP( npar2, npar1, phead2, phead1, karr2,
                                    karr1, m2, m1, sqrt_var );
#endif

#if 0
            // De-macro-ize the above to check where we die.
            if ( npar1>npar2 ) {

              if ( npar1>0 && npar2>0 ) {
                int ii, iimax, index1=-1, index2=-1;

                ii=(npar1)/(npar2);
                iimax=(npar1)-(npar2)*ii;
                for ( i=0; i<iimax; ++i) {
                  ++index2;
                  for ( j=0; j<ii+1; ++j ) {
                    ++index1;
                    SCATTER_PARTICLES( karr1[index1].i, karr2[index2].i,
                                       phead1, phead2, m1, m2,
                                       (sqrt_var)/sqrt(ii+1) );
                  }
                }
                for ( ; i<(npar2); ++i ) {
                  ++index2;
                  for ( j=0; j<ii; ++j ) {
                    ++index1;
                    SCATTER_PARTICLES( karr1[index1].i, karr2[index2].i,
                                       phead1, phead2, m1, m2,
                                       (sqrt_var)/sqrt(ii) );
                  }
                }
              }

            } else {

              if ( npar2>0 && npar1>0 ) {
                int ii, iimax, index1=-1, index2=-1;

                ii=(npar2)/(npar1);
                iimax=(npar2)-(npar1)*ii;
                for ( i=0; i<iimax; ++i) {
                  ++index2;
                  for ( j=0; j<ii+1; ++j ) {
                    ++index1;
                    SCATTER_PARTICLES( karr2[index1].i, karr1[index2].i,
                                       phead2, phead1, m2, m1,
                                       (sqrt_var)/sqrt(ii+1) );
                  }
                }
                for ( ; i<(npar1); ++i ) {
                  ++index2;
                  for ( j=0; j<ii; ++j ) {
                    ++index1;
                    SCATTER_PARTICLES( karr2[index1].i, karr1[index2].i,
                                       phead2, phead1, m2, m1,
                                       (sqrt_var)/sqrt(ii) );
                  }
                }
              }
            } // end if
#endif
          }
        }
#endif
      }
    }
  }
}

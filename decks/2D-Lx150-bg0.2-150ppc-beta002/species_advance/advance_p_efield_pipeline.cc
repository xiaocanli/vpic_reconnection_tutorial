// FIXME: PARTICLE MOVERS NEED TO BE OVERALLOCATED IN STRUCTORS TO
// ACCOUNT FOR SPLITTING THE MOVER ARRAY BETWEEN HOST AND PIPELINES

/* #define IN_spa */

#define HAS_V4_PIPELINE
#define HAS_V8_PIPELINE
#define HAS_V16_PIPELINE

#include "spa_private_efield.h"

#include "util/pipelines/pipelines_exec.h"

//----------------------------------------------------------------------------//
// Reference implementation for an advance_p pipeline function which does not
// make use of explicit calls to vector intrinsic functions.
//----------------------------------------------------------------------------//

void
advance_p_efield_pipeline_scalar( advance_p_efield_pipeline_args_t * args,
                                  int pipeline_rank,
                                  int n_pipeline )
{
  particle_t           * ALIGNED(128) p0 = args->p0;
  accumulator_t        * ALIGNED(128) a0 = args->a0;
  const interpolator_t * ALIGNED(128) f0 = args->f0;
  const grid_t *                      g  = args->g;

  particle_t           * ALIGNED(32)  p;
  particle_mover_t     * ALIGNED(16)  pm;
  const interpolator_t * ALIGNED(16)  f;
  float                * ALIGNED(16)  a;

  const int etype            = args->etype; // electric field type

  const float qdt_2mc        = args->qdt_2mc;
  const float cdt_dx         = args->cdt_dx;
  const float cdt_dy         = args->cdt_dy;
  const float cdt_dz         = args->cdt_dz;
  const float qsp            = args->qsp;
  const float msp            = args->msp;
  const float gamma_sironi   = args->gamma_sironi;
  const float one            = 1.0;
  const float one_third      = 1.0/3.0;
  const float two_fifteenths = 2.0/15.0;

  float dx, dy, dz, ux, uy, uz, q;
  float hax, hay, haz, cbx, cby, cbz;
  float v0, v1, v2, v3, v4, v5;
  float ex, ey, ez;
  float eparax, eparay, eparaz, ie2b2, w2, wtmp;
  float edotb, btot2, etot2, ib2;
  float gamma_scale;
  int   ii;

  int itmp, n, nm, max_nm;

  DECLARE_ALIGNED_ARRAY( particle_mover_t, 16, local_pm, 1 );

  // Determine which quads of particles quads this pipeline processes.

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, itmp, n );

  p = args->p0 + itmp;

  // Determine which movers are reserved for this pipeline.
  // Movers (16 bytes) should be reserved for pipelines in at least
  // multiples of 8 such that the set of particle movers reserved for
  // a pipeline is 128-byte aligned and a multiple of 128-byte in
  // size.  The host is guaranteed to get enough movers to process its
  // particles with this allocation.

  max_nm = args->max_nm - ( args->np&15 );

  if ( max_nm < 0 ) max_nm = 0;

  DISTRIBUTE( max_nm, 8, pipeline_rank, n_pipeline, itmp, max_nm );

  if ( pipeline_rank == n_pipeline ) max_nm = args->max_nm - itmp;

  pm   = args->pm + itmp;
  nm   = 0;
  itmp = 0;

  // Determine which accumulator array to use
  // The host gets the first accumulator array.

  if ( pipeline_rank != n_pipeline )
    a0 += ( 1 + pipeline_rank ) *
          POW2_CEIL( (args->nx+2)*(args->ny+2)*(args->nz+2), 2 );

  // Process particles for this pipeline.

  for( ; n; n--, p++ )
  {
    dx   = p->dx;                             // Load position
    dy   = p->dy;
    dz   = p->dz;
    ii   = p->i;

    f    = f0 + ii;                           // Interpolate E

    ex = (f->ex + dy*f->dexdy) + dz*(f->dexdz + dy*f->d2exdydz);
    ey = (f->ey + dz*f->deydz) + dx*(f->deydx + dz*f->d2eydzdx);
    ez = (f->ez + dx*f->dezdx) + dy*(f->dezdy + dx*f->d2ezdxdy);

    cbx  = f->cbx + dx*f->dcbxdx;             // Interpolate B
    cby  = f->cby + dy*f->dcbydy;
    cbz  = f->cbz + dz*f->dcbzdz;

    ux   = p->ux;                             // Load momentum
    uy   = p->uy;
    uz   = p->uz;
    q    = p->w;

    btot2 = cbx*cbx + cby*cby + cbz*cbz;
    ib2 = 1.0 / btot2;
    edotb = ex * cbx + ey * cby + ez * cbz;
    eparax = edotb * cbx * ib2;
    eparay = edotb * cby * ib2;
    eparaz = edotb * cbz * ib2;

    if (etype == 0) { // remove all parallel electric field
      hax  = qdt_2mc * (ex - eparax);
      hay  = qdt_2mc * (ey - eparay);
      haz  = qdt_2mc * (ez - eparaz);
    } else if (etype == 1) { // remove the y-component of parallel electric field
      hax  = qdt_2mc * ex;
      hay  = qdt_2mc * (ey - eparay);
      haz  = qdt_2mc * ez;
    } else if (etype == 2 || etype == 3 || etype == 4) { // in regions where |E| > |B|
      // 2: remove electric field
      // 3: or set gamma to a few (Sironi's method)
      etot2 = ex*ex + ey*ey + ez*ez;
      if (etot2 < btot2) {
        hax  = qdt_2mc * ex;
        hay  = qdt_2mc * ey;
        haz  = qdt_2mc * ez;
      } else {
        hax  = 0.0;
        hay  = 0.0;
        haz  = 0.0;
        if (etype == 3 || (etype == 4 && p->w < 10)) { // set gamma to a few (Sironi's method)
          v0 = one / sqrtf( one + ( ux*ux+ ( uy*uy + uz*uz ) ) );
          gamma_scale = gamma_sironi * v0 / msp;
          ux *= gamma_scale;
          uy *= gamma_scale;
          uz *= gamma_scale;
        }
        p->w = 100; // use w to mark these particles
      }
    } else if (etype == 5) { // remove perpendicular electric field but keep ExB drift
      hax  = qdt_2mc * eparax;
      hay  = qdt_2mc * eparay;
      haz  = qdt_2mc * eparaz;
    } else {
      hax  = qdt_2mc * ex;
      hay  = qdt_2mc * ey;
      haz  = qdt_2mc * ez;
    }

    ux  += hax;                               // Half advance E
    uy  += hay;
    uz  += haz;

    v0   = qdt_2mc / sqrtf( one + ( ux*ux + ( uy*uy + uz*uz ) ) );

                                              // Boris - scalars
    v1   = cbx*cbx + ( cby*cby + cbz*cbz );
    v2   = ( v0*v0 ) * v1;
    v3   = v0 * ( one + v2 * ( one_third + v2 * two_fifteenths ) );
    v4   = v3 / ( one + v1 * ( v3 * v3 ) );
    v4  += v4;

    v0   = ux + v3*( uy*cbz - uz*cby );       // Boris - uprime
    v1   = uy + v3*( uz*cbx - ux*cbz );
    v2   = uz + v3*( ux*cby - uy*cbx );

    ux  += v4*( v1*cbz - v2*cby );            // Boris - rotation
    uy  += v4*( v2*cbx - v0*cbz );
    uz  += v4*( v0*cby - v1*cbx );

    ux  += hax;                               // Half advance E
    uy  += hay;
    uz  += haz;

    p->ux = ux;                               // Store momentum
    p->uy = uy;
    p->uz = uz;

    v0   = one / sqrtf( one + ( ux*ux+ ( uy*uy + uz*uz ) ) );
                                              // Get norm displacement
                                              //
    ux  *= cdt_dx;
    uy  *= cdt_dy;
    uz  *= cdt_dz;

    ux  *= v0;
    uy  *= v0;
    uz  *= v0;

    if (etype == 5) { // remove perpendicular electric field but keep ExB drift
      /* eparax = (ey * cbz - ez * cby) * ib2; // reuse epara as ExB drift */
      /* eparay = (ez * cbx - ex * cbz) * ib2; */
      /* eparaz = (ex * cby - ey * cbx) * ib2; */
      // https://aip.scitation.org/doi/10.1063/1.873736 for calculating ExB drift
      w2 = ex*ex + ey*ey + ez*ez + cbx*cbx + cby*cby + cbz*cbz;
      if (w2 > 0.0) {
        ie2b2 = 1.0/w2;
      } else {
        ie2b2 = 0.0;
      }
      eparax = (ey * cbz - ez * cby) * ie2b2; // reuse epara as ExB/(E^2+B^2)
      eparay = (ez * cbx - ex * cbz) * ie2b2;
      eparaz = (ex * cby - ey * cbx) * ie2b2;
      w2 = eparax*eparax + eparay*eparay + eparaz*eparaz;
      if (w2 < 0.25) {
        wtmp = sqrt(1.0-4*w2);
      } else {
        wtmp = 0.0;
      }
      if (w2 > 0.0) {
        w2 = 0.5 * (1.0 - wtmp) / w2; // reuse w2
      } else {
        w2 = 0.0;
      }
      eparax *= w2; // reuse epara as ExB drift
      eparay *= w2;
      eparaz *= w2;
      ux += eparax * cdt_dx;
      uy += eparay * cdt_dy;
      uz += eparaz * cdt_dz;
    }

    v0   = dx + ux;                           // Streak midpoint (inbnds)
    v1   = dy + uy;
    v2   = dz + uz;

    v3   = v0 + ux;                           // New position
    v4   = v1 + uy;
    v5   = v2 + uz;

    // FIXME-KJB: COULD SHORT CIRCUIT ACCUMULATION IN THE CASE WHERE QSP==0!
    if (  v3 <= one &&  v4 <= one &&  v5 <= one &&   // Check if inbnds
         -v3 <= one && -v4 <= one && -v5 <= one )
    {
      // Common case (inbnds).  Note: accumulator values are 4 times
      // the total physical charge that passed through the appropriate
      // current quadrant in a time-step.

      q *= qsp;

      p->dx = v3;                             // Store new position
      p->dy = v4;
      p->dz = v5;

      dx = v0;                                // Streak midpoint
      dy = v1;
      dz = v2;

      v5 = q*ux*uy*uz*one_third;              // Compute correction

      a  = (float *)( a0 + ii );              // Get accumulator

#     define ACCUMULATE_J(X,Y,Z,offset)                                 \
      v4  = q*u##X;   /* v2 = q ux                            */        \
      v1  = v4*d##Y;  /* v1 = q ux dy                         */        \
      v0  = v4-v1;    /* v0 = q ux (1-dy)                     */        \
      v1 += v4;       /* v1 = q ux (1+dy)                     */        \
      v4  = one+d##Z; /* v4 = 1+dz                            */        \
      v2  = v0*v4;    /* v2 = q ux (1-dy)(1+dz)               */        \
      v3  = v1*v4;    /* v3 = q ux (1+dy)(1+dz)               */        \
      v4  = one-d##Z; /* v4 = 1-dz                            */        \
      v0 *= v4;       /* v0 = q ux (1-dy)(1-dz)               */        \
      v1 *= v4;       /* v1 = q ux (1+dy)(1-dz)               */        \
      v0 += v5;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */        \
      v1 -= v5;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */        \
      v2 -= v5;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */        \
      v3 += v5;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */        \
      a[offset+0] += v0;                                                \
      a[offset+1] += v1;                                                \
      a[offset+2] += v2;                                                \
      a[offset+3] += v3

      ACCUMULATE_J( x, y, z, 0 );
      ACCUMULATE_J( y, z, x, 4 );
      ACCUMULATE_J( z, x, y, 8 );

#     undef ACCUMULATE_J
    }

    else                                        // Unlikely
    {
      local_pm->dispx = ux;
      local_pm->dispy = uy;
      local_pm->dispz = uz;

      local_pm->i     = p - p0;

      if ( move_p( p0, local_pm, a0, g, qsp ) ) // Unlikely
      {
        if ( nm < max_nm )
        {
          pm[nm++] = local_pm[0];
        }

        else
        {
          itmp++;                               // Unlikely
        }
      }
    }
  }

  args->seg[pipeline_rank].pm        = pm;
  args->seg[pipeline_rank].max_nm    = max_nm;
  args->seg[pipeline_rank].nm        = nm;
  args->seg[pipeline_rank].n_ignored = itmp;
}

//----------------------------------------------------------------------------//
// Top level function to select and call the proper advance_p pipeline
// function.
//----------------------------------------------------------------------------//

// advance particles without parallel electric field
void
advance_p_efield_pipeline( species_t * RESTRICT sp,
                           accumulator_array_t * RESTRICT aa,
                           const interpolator_array_t * RESTRICT ia,
                           const int efield_type,
                           const float gamma_sironi )
{
  DECLARE_ALIGNED_ARRAY( advance_p_efield_pipeline_args_t, 128, args, 1 );

  DECLARE_ALIGNED_ARRAY( particle_mover_seg_t, 128, seg, MAX_PIPELINE + 1 );

  int rank;

  if ( !sp || !aa || !ia || sp->g != aa->g || sp->g != ia->g )
  {
    ERROR( ( "Bad args" ) );
  }

  args->p0      = sp->p;
  args->pm      = sp->pm;
  args->a0      = aa->a;
  args->f0      = ia->i;
  args->seg     = seg;
  args->g       = sp->g;

  args->qdt_2mc = (sp->q*sp->g->dt)/(2*sp->m*sp->g->cvac);
  args->cdt_dx  = sp->g->cvac*sp->g->dt*sp->g->rdx;
  args->cdt_dy  = sp->g->cvac*sp->g->dt*sp->g->rdy;
  args->cdt_dz  = sp->g->cvac*sp->g->dt*sp->g->rdz;
  args->qsp     = sp->q;
  args->msp     = sp->m;
  args->gamma_sironi = gamma_sironi;

  args->np      = sp->np;
  args->max_nm  = sp->max_nm;
  args->nx      = sp->g->nx;
  args->ny      = sp->g->ny;
  args->nz      = sp->g->nz;
  args->etype   = efield_type;

  // Have the host processor do the last incomplete bundle if necessary.
  // Note: This is overlapped with the pipelined processing.  As such,
  // it uses an entire accumulator.  Reserving an entire accumulator
  // for the host processor to handle at most 15 particles is wasteful
  // of memory.  It is anticipated that it may be useful at some point
  // in the future have pipelines accumulating currents while the host
  // processor is doing other more substantive work (e.g. accumulating
  // currents from particles received from neighboring nodes).
  // However, it is worth reconsidering this at some point in the
  // future.

  EXEC_PIPELINES( advance_p_efield, args, 0 );

  WAIT_PIPELINES();

  // FIXME: HIDEOUS HACK UNTIL BETTER PARTICLE MOVER SEMANTICS
  // INSTALLED FOR DEALING WITH PIPELINES.  COMPACT THE PARTICLE
  // MOVERS TO ELIMINATE HOLES FROM THE PIPELINING.

  sp->nm = 0;
  for( rank = 0; rank <= N_PIPELINE; rank++ )
  {
    if ( args->seg[rank].n_ignored )
    {
      WARNING( ( "Pipeline %i ran out of storage for %i movers",
                 rank, args->seg[rank].n_ignored ) );
    }

    if ( sp->pm + sp->nm != args->seg[rank].pm )
    {
      MOVE( sp->pm + sp->nm, args->seg[rank].pm, args->seg[rank].nm );
    }

    sp->nm += args->seg[rank].nm;
  }
}

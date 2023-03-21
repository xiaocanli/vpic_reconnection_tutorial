#ifndef _spa_private_efield_h_
#define _spa_private_efield_h_

#include "species_advance_efield.h"

///////////////////////////////////////////////////////////////////////////////
// advance_p_pipeline interface

typedef struct particle_mover_seg
{
  MEM_PTR( particle_mover_t, 16 ) pm; // First mover in segment
  int max_nm;                         // Maximum number of movers
  int nm;                             // Number of movers used
  int n_ignored;                      // Number of movers ignored

  PAD_STRUCT( SIZEOF_MEM_PTR+3*sizeof(int) )

} particle_mover_seg_t;

typedef struct advance_p_efield_pipeline_args
{
  MEM_PTR( particle_t,           128 ) p0;       // Particle array
  MEM_PTR( particle_mover_t,     128 ) pm;       // Particle mover array
  MEM_PTR( accumulator_t,        128 ) a0;       // Accumulator arrays
  MEM_PTR( const interpolator_t, 128 ) f0;       // Interpolator array
  MEM_PTR( particle_mover_seg_t, 128 ) seg;      // Dest for return values
  MEM_PTR( const grid_t,         1   ) g;        // Local domain grid params

  float                                qdt_2mc;  // Particle/field coupling
  float                                cdt_dx;   // x-space/time coupling
  float                                cdt_dy;   // y-space/time coupling
  float                                cdt_dz;   // z-space/time coupling
  float                                qsp;      // Species particle charge

  int                                  np;       // Number of particles
  int                                  max_nm;   // Number of movers
  int                                  nx;       // x-mesh resolution
  int                                  ny;       // y-mesh resolution
  int                                  nz;       // z-mesh resolution
  int                                  etype;    // electric field type

  PAD_STRUCT( 6*SIZEOF_MEM_PTR + 5*sizeof(float) + 6*sizeof(int) )

} advance_p_efield_pipeline_args_t;

void
advance_p_efield_pipeline_scalar( advance_p_efield_pipeline_args_t * args,
                                  int pipeline_rank,
                                  int n_pipeline );
void
advance_p_efield_pipeline_v4( advance_p_efield_pipeline_args_t * args,
                              int pipeline_rank,
                              int n_pipeline );

void
advance_p_efield_pipeline_v8( advance_p_efield_pipeline_args_t * args,
                              int pipeline_rank,
                              int n_pipeline );

#endif // _spa_private_efield_h_

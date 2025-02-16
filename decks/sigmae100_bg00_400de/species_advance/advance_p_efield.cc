/* #define IN_spa */

#include "species_advance_efield.h"

//----------------------------------------------------------------------------//
// Top level function to select and call particle advance function using the
// desired particle advance abstraction. Currently, the only abstraction
// available is the pipeline abstraction. Particles advanced by this function
// experience only part of electric field and have no feedback to the fields.
// efield_type=0: without parallel electric field
// efield_type=1: without the y-component of parallel electric field
// efield_type=2: without electric field in regions where |E|>|B|
//----------------------------------------------------------------------------//

void
advance_p_efield( species_t * RESTRICT sp,
                  accumulator_array_t * RESTRICT aa,
                  const interpolator_array_t * RESTRICT ia,
                  const int efield_type,
                  const float gamma_sironi )
{
  // Once more options are available, this should be conditionally executed
  // based on user choice.
  advance_p_efield_pipeline( sp, aa, ia, efield_type, gamma_sironi );
}

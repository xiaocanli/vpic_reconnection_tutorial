#ifndef _species_advance_efield_h_
#define _species_advance_efield_h_

#include "sf_interface/sf_interface.h"
//----------------------------------------------------------------------------//
// Choose between using AoSoA or AoS data layout for the particles.
//----------------------------------------------------------------------------//

#include "species_advance/species_advance_aos.h"
#include "species_advance/species_advance.h"

//----------------------------------------------------------------------------//
// Declare methods.
//----------------------------------------------------------------------------//

BEGIN_C_DECLS

// In advance_p_efield.cc

void
advance_p_efield( species_t * RESTRICT sp,
                  accumulator_array_t * RESTRICT aa,
                  const interpolator_array_t * RESTRICT ia,
                  const int efield_type,
                  const float gamma_sironi);

void
advance_p_efield_pipeline( species_t * RESTRICT sp,
                           accumulator_array_t * RESTRICT aa,
                           const interpolator_array_t * RESTRICT ia,
                           const int efield_type,
                           const float gamma_sironi);

END_C_DECLS

#endif // _species_advance_efield_h_

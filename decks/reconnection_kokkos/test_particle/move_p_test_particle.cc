#define IN_spa
#include "species_advance/species_advance.h"
#include "vpic/kokkos_helpers.h"

int
move_p_test_particle( particle_t       * ALIGNED(128) p0,
                      particle_mover_t * ALIGNED(16)  pm,
                      const grid_t     *              g,
                      const float                     qsp ) {

  float cx = 0.25 * g->rdy * g->rdz / g->dt;
  float cy = 0.25 * g->rdz * g->rdx / g->dt;
  float cz = 0.25 * g->rdx * g->rdy / g->dt;

  float s_midx, s_midy, s_midz;
  float s_dispx, s_dispy, s_dispz;
  float s_dir[3];
  float v0, v1, v2, v3;
  int axis, face;
  int64_t neighbor;
  //float *a;
  particle_t * ALIGNED(32) p = p0 + pm->i;

  for(;;) {
    s_midx = p->dx;
    s_midy = p->dy;
    s_midz = p->dz;

    s_dispx = pm->dispx;
    s_dispy = pm->dispy;
    s_dispz = pm->dispz;

    s_dir[0] = (s_dispx>0) ? 1 : -1;
    s_dir[1] = (s_dispy>0) ? 1 : -1;
    s_dir[2] = (s_dispz>0) ? 1 : -1;
    
    // Compute the twice the fractional distance to each potential
    // streak/cell face intersection.
    v0 = (s_dispx==0) ? 3.4e38 : (s_dir[0]-s_midx)/s_dispx;
    v1 = (s_dispy==0) ? 3.4e38 : (s_dir[1]-s_midy)/s_dispy;
    v2 = (s_dispz==0) ? 3.4e38 : (s_dir[2]-s_midz)/s_dispz;

    // Determine the fractional length and axis of current streak. The
    // streak ends on either the first face intersected by the
    // particle track or at the end of the particle track.
    // 
    //   axis 0,1 or 2 ... streak ends on a x,y or z-face respectively
    //   axis 3        ... streak ends at end of the particle track
    /**/      v3=2,  axis=3;
    if(v0<v3) v3=v0, axis=0;
    if(v1<v3) v3=v1, axis=1;
    if(v2<v3) v3=v2, axis=2;
    v3 *= 0.5;

    // Compute the midpoint and the normalized displacement of the streak
    s_dispx *= v3;
    s_dispy *= v3;
    s_dispz *= v3;
    s_midx += s_dispx;
    s_midy += s_dispy;
    s_midz += s_dispz;

    // Compute the remaining particle displacment
    pm->dispx -= s_dispx;
    pm->dispy -= s_dispy;
    pm->dispz -= s_dispz;

    // Compute the new particle offset
    p->dx += s_dispx+s_dispx;
    p->dy += s_dispy+s_dispy;
    p->dz += s_dispz+s_dispz;

    // If an end streak, return success (should be ~50% of the time)

    if( axis==3 ) break;

    // Determine if the particle crossed into a local cell or if it
    // hit a boundary and convert the coordinate system accordingly.
    // Note: Crossing into a local cell should happen ~50% of the
    // time; hitting a boundary is usually a rare event.  Note: the
    // entry / exit coordinate for the particle is guaranteed to be
    // +/-1 _exactly_ for the particle.

    v0 = s_dir[axis];
    (&(p->dx))[axis] = v0; // Avoid roundoff fiascos--put the particle
                           // _exactly_ on the boundary.
    face = axis; if( v0>0 ) face += 3;
    neighbor = g->neighbor[ 6*p->i + face ];
    
    if( UNLIKELY( neighbor==reflect_particles ) ) {
      // Hit a reflecting boundary condition.  Reflect the particle
      // momentum and remaining displacement and keep moving the
      // particle.
      (&(p->ux    ))[axis] = -(&(p->ux    ))[axis];
      (&(pm->dispx))[axis] = -(&(pm->dispx))[axis];
      continue;
    }

    if( UNLIKELY( neighbor<g->rangel || neighbor>g->rangeh ) ) {
      // Cannot handle the boundary condition here.  Save the updated
      // particle position, face it hit and update the remaining
      // displacement in the particle mover.
      p->i = 8*p->i + face;
      return 1; // Return "mover still in use"
    }

    // Crossed into a normal voxel.  Update the voxel index, convert the
    // particle coordinate system and keep moving the particle.
    
    p->i = neighbor - g->rangel; // Compute local index of neighbor
    /**/                         // Note: neighbor - g->rangel < 2^31 / 6
    (&(p->dx))[axis] = -v0;      // Convert coordinate system
  }

  return 0; // Return "mover not in use"
}

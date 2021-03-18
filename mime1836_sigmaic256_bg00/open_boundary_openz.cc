{
  int inject;
  double x, y, z, age, flux, vtherm, vd;
  double uv[3];
  double zcell;
  const int nsp=global->nsp;
  const int nx=grid->nx;
  const int ny=grid->ny;
  const int nz=grid->nz;
  const double rin[3]={global->rin[0],global->rin[1],global->rin[2]};
  const double rout[3]={global->rout[0],global->rout[1],global->rout[2]};
  const double sqpi =1.772453850905516;
  const double dt=grid->dt;
  const double hx=grid->dx;
  const double hy=grid->dy;
  const double hz=grid->dz;
  const double nb=global->nb;
  const double nfac=global->nfac;
  const double muR=global->muR;

  // Initialize the injectors on the first call

  static int initted=0;
  if ( !initted ) {
    initted=1;

    if (rank() == 0) MESSAGE(("----------------Initializing the Particle Injectors-----------------"));

    // Intialize injectors for Harris Sheet with a uniform background

    if (global->right) {
      DEFINE_INJECTOR(right,ny,nz);

      // double muR =3 ;

      if (step() ==0) {
        for ( int n=1; n<=nsp; n++ ) {
          double cn = (uf(2)/vth(2))/(vth(n)/vth(2));
          for ( int k=1;k<=nz; k++ ) {
            for ( int j=1;j<=ny; j++ ) {
              bright(n,k,j) = 0;
              zcell = (grid->z0 + k*hz-hz/2)/(global->L_de);
              nright(n,k,j) = (nb + muR*muR/(cosh(muR*zcell)*cosh(muR*zcell)))/nfac;
              uright(1,n,k,j) = 0;
              uright(2,n,k,j) = uf(n)/(1+nb/muR/muR*cosh(muR*zcell)*cosh(muR*zcell));
              uright(3,n,k,j) = 0;
              pright(1,2,n,k,j)=pright(2,1,n,k,j)=pright(1,3,n,k,j)=0;
              pright(3,1,n,k,j)=pright(2,3,n,k,j)=pright(3,2,n,k,j)=0;
              pright(1,1,n,k,j) = (nb*vthb(n)*vthb(n) +
                  vth(n)*vth(n)*muR*muR/(cosh(muR*zcell)*cosh(muR*zcell)))/(2*nfac);
              pright(2,2,n,k,j) = (nb*vthb(n)*vthb(n) +
                  vth(n)*vth(n)*(muR*muR/(cosh(muR*zcell)*cosh(muR*zcell)) +
                    2*nb*cn*cn/(1+nb/muR/muR*cosh(muR*zcell)*cosh(muR*zcell))))/(2*nfac);
              pright(3,3,n,k,j) = pright(1,1,n,k,j);
            }
          }
        }  // end for
      } else {
        READ_INJECTOR(right,ny,nz,0);
      }
    } //end right boundary

    if (global->left) {

      DEFINE_INJECTOR(left,ny,nz);

      if (step()==0) {
        for ( int n=1; n<=nsp; n++ ) {
          double cn = (uf(2)/vth(2))/(vth(n)/vth(2));
          for ( int k=1;k<nz+1; k++ ) {
            for ( int j=1;j<=ny; j++ ) {
              bleft(n,k,j) = 0;
              zcell = (grid->z0 + k*hz-hz/2)/(global->L_de);
              nleft(n,k,j) = (nb + 1/(cosh(zcell)*cosh(zcell)))/nfac;
              uleft(1,n,k,j) = 0;
              uleft(2,n,k,j) = uf(n)/(1+nb*cosh(zcell)*cosh(zcell));
              uleft(3,n,k,j) = 0;
              pleft(1,2,n,k,j)=pleft(2,1,n,k,j)=pleft(1,3,n,k,j)=0;
              pleft(3,1,n,k,j)=pleft(2,3,n,k,j)=pleft(3,2,n,k,j)=0;
              pleft(1,1,n,k,j) = (nb*vthb(n)*vthb(n) +
                  vth(n)*vth(n)/(cosh(zcell)*cosh(zcell)))/(2*nfac);
              pleft(2,2,n,k,j) = (nb*vthb(n)*vthb(n) +
                  vth(n)*vth(n)*(1/(cosh(zcell)*cosh(zcell)) +
                    2*nb*cn*cn/(1+nb*cosh(zcell)*cosh(zcell))))/(2*nfac);
              pleft(3,3,n,k,j) = pleft(1,1,n,k,j);
            }
          }
        } // end for
      } else {
        READ_INJECTOR(left,ny,nz,0);
      }
    } // end left boundary

    if (global->top) {

      DEFINE_INJECTOR(top,ny,nx);

      if (step()==0) {
        for ( int n=1; n<=nsp; n++ ) {
          for ( int i=1;i<=nx; i++ ) {
            for ( int j=1;j<=ny; j++ ) {
              btop(n,i,j) = 0;
              ntop(n,i,j) = nb/nfac;
              //ftop(n,i,j) = ntop(n,i,j)*vthb(n)/(2*hz*sqpi);
              utop(1,n,i,j) = 0;
              utop(2,n,i,j) = 0;
              utop(3,n,i,j) = 0;
              ptop(1,2,n,i,j)=ptop(2,1,n,i,j)=ptop(1,3,n,i,j)=0;
              ptop(3,1,n,i,j)=ptop(2,3,n,i,j)=ptop(3,2,n,i,j)=0;
              ptop(1,1,n,i,j) = ntop(n,i,j)*vthb(n)*vthb(n)/2;;
              ptop(2,2,n,i,j) = ntop(n,i,j)*vthb(n)*vthb(n)/2;;
              ptop(3,3,n,i,j) = ntop(n,i,j)*vthb(n)*vthb(n)/2;;
            }
          }
        } // end for
      } else {
        READ_INJECTOR(top,ny,nx,0);
      }

    } // end top boundary

    if (global->bottom) {
      DEFINE_INJECTOR(bot,ny,nx);
      if (step() ==0) {
        for ( int n=1; n<=nsp; n++ ) {
          for ( int i=1;i<=nx; i++ ) {
            for ( int j=1;j<=ny; j++ ) {
              bbot(n,i,j) = 0;
              nbot(n,i,j) = nb/nfac;
              //fbot(n,i,j) = nbot(n,i,j)*vthb(n)/(2*hz*sqpi);
              ubot(1,n,i,j) = 0.0;
              ubot(2,n,i,j) = 0.0;
              ubot(3,n,i,j) = 0.0;
              pbot(1,2,n,i,j)=pbot(2,1,n,i,j)=pbot(1,3,n,i,j)=0;
              pbot(3,1,n,i,j)=pbot(2,3,n,i,j)=pbot(3,2,n,i,j)=0;
              pbot(1,1,n,i,j) = nbot(n,i,j)*vthb(n)*vthb(n)/2;
              pbot(2,2,n,i,j) = nbot(n,i,j)*vthb(n)*vthb(n)/2;
              pbot(3,3,n,i,j) = nbot(n,i,j)*vthb(n)*vthb(n)/2;
            }
          }
        } // end for
      } else {
        READ_INJECTOR(bot,ny,nx,0);
      }
    }  // end bottom boundary

    if (rank() == 0) MESSAGE(("-------------------------------------------------------------------"));

  }// End of Intialization

  //  Inject particles on Left Boundary

  if (global->left) {
    double t=grid->dt*step();
    for ( int n=1; n<=nsp; n++ ) {
      species_t * species = find_species_id(n-1,species_list );
      for ( int k=1;k<=nz; k++ ) {
        for ( int j=1;j<=ny; j++ ) {
          // bleft(n,k,j) = bleft(n,k,j) + dt*fleft(n,k,j);
          vtherm = sqrt(2.0*pleft(1,1,n,k,j)/nleft(n,k,j));
          double vexb_x, vexb_y, vexb_z;
          if (step() <= global->nsteps_relax) {
            vexb_x = uleft(1,n,k,j);
            vexb_y = uleft(2,n,k,j);
            vexb_z = uleft(3,n,k,j);
          } else {
            // https://aip.scitation.org/doi/10.1063/1.873736 for calculating ExB drift
            double ycell = grid->y0 + (j-0.5)*grid->dy; // ycell ranges from -Ly/2 to Ly/2
            double zcell = grid->z0 + (k-1.0)*grid->dz; // zcell ranges from -Lz/2 to Lz/2
            /* double ex = field(1,j,k).ex; */
            /* double ez = field(1,j,k).ez; */
            double ex = 0.0;
            double ez = 0.0;
            double ey = global->eyd0 * pow(cosh(ycell/global->ldy), -2);
            ey *= pow(cosh((zcell-global->zd)/global->ldz), -2) +
                  pow(cosh((zcell+global->zd)/global->ldz), -2);
            ey *= 1.0 - exp(-(t-global->trelax)/global->tdrive_lu);
            double e2 = ex*ex + ey*ey + ez*ez;
            double bx = field(1,j,k).cbx;
            double by = field(1,j,k).cby;
            double bz = field(1,j,k).cbz;
            double b2 = bx*bx + by*by + bz*bz;
            double ieb2 = 1.0 / (e2 + b2);
            double wx = (ey*bz - ez*by) * ieb2;
            double wy = (ez*bx - ex*bx) * ieb2;
            double wz = (ex*by - ey*bz) * ieb2;
            double w2 = wx*wx + wy*wy + wz*wz;
            double wtmp = 0.5 * (1.0 - sqrt(1.0-4.0*w2)) / w2;
            vexb_x = wx * wtmp;
            vexb_y = wy * wtmp;
            vexb_z = wz * wtmp;
          }
          vd = uleft(1,n,k,j)/vtherm;
          bleft(n,k,j) += dt*nleft(n,k,j)*vtherm*(exp(-vd*vd)/sqpi+vd*(erf(vd)+1))/(2*hx);
          inject = (int) bleft(n,k,j);
          bleft(n,k,j) = bleft(n,k,j) - (double) inject;
          double uflow[3] = {uleft(1,n,k,j),uleft(2,n,k,j),uleft(3,n,k,j)};
          /* double uflow[3] = {vexb_x,vexb_y,vexb_z}; */
          /* if (fabs(global->eyd0) > 0.0) { */
          /*   uflow[0] = vexb_x; */
          /*   uflow[2] = vexb_z; */
          /* } */
          double press[9] = {pleft(1,1,n,k,j),pleft(1,2,n,k,j),pleft(1,3,n,k,j),
                             pleft(2,1,n,k,j),pleft(2,2,n,k,j),pleft(2,3,n,k,j),
                             pleft(3,1,n,k,j),pleft(3,2,n,k,j),pleft(3,3,n,k,j)};
          repeat(inject) {
            compute_injection(uv,nleft(n,k,j),uflow,press,1,2,3,rng(0));
            x = grid->x0;
            y = grid->y0 + hy*(j-1) + hy*uniform(rng(0),0,1);
            z = grid->z0 + hz*(k-1) + hz*uniform(rng(0),0,1);
            age = 0;
            inject_particle(species, x, y, z, uv[0], uv[1], uv[2], q(n), age, 1 );
          }
        }
      }
    }
  } // end left injector

    //  Inject particles on Right Boundary

  if (global->right) {
    for ( int n=1; n<=nsp; n++ ) {
      species_t * species = find_species_id(n-1,species_list );
      for ( int k=1;k<=nz; k++ ) {
        for ( int j=1;j<=ny; j++ ) {
          //bright(n,k,j) = bright(n,k,j) + dt*fright(n,k,j);
          vtherm = sqrt(2.0*pright(1,1,n,k,j)/nright(n,k,j));
          vd = -uright(1,n,k,j)/vtherm;
          bright(n,k,j) += dt*nright(n,k,j)*vtherm*(exp(-vd*vd)/sqpi+vd*(erf(vd)+1))/(2*hx);

          inject = (int) bright(n,k,j);
          bright(n,k,j) = bright(n,k,j) - (double) inject;
          double uflow[3] = {uright(1,n,k,j),uright(2,n,k,j),uright(3,n,k,j)};
          double press[9] = {pright(1,1,n,k,j),pright(1,2,n,k,j),pright(1,3,n,k,j),
                             pright(2,1,n,k,j),pright(2,2,n,k,j),pright(2,3,n,k,j),
                             pright(3,1,n,k,j),pright(3,2,n,k,j),pright(3,3,n,k,j)};
          repeat(inject) {
            compute_injection(uv,nright(n,k,j),uflow,press,-1,2,3,rng(0));
            x = grid->x1;
            y = grid->y0 + hy*(j-1) + hy*uniform(rng(0),0,1);
            z = grid->z0 + hz*(k-1) + hz*uniform(rng(0),0,1);
            age = 0;
            inject_particle(species, x, y, z, uv[0], uv[1], uv[2], q(n), age, 1 );
          }
        }
      }
    }
  } // end right injector

  //  Inject particles on Top Boundary

  if (global->top) {
    for ( int n=1; n<=nsp; n++ ) {
      species_t * species = find_species_id(n-1,species_list );
      for ( int i=1;i<=nx; i++ ) {
        for ( int j=1;j<=ny; j++ ) {
          vtherm = sqrt(2.0*ptop(3,3,n,i,j)/ntop(n,i,j));
          vd = -utop(3,n,i,j)/vtherm;
          btop(n,i,j) += dt*ntop(n,i,j)*vtherm*(exp(-vd*vd)/sqpi+vd*(erf(vd)+1))/(2*hz);

          inject = (int) btop(n,i,j);
          btop(n,i,j) = btop(n,i,j)- (double) inject;
          double uflow[3] = {0,0,utop(3,n,i,j)};

          double press[9] = {ptop(1,1,n,i,j),ptop(1,2,n,i,j),ptop(1,3,n,i,j),
                             ptop(2,1,n,i,j),ptop(2,2,n,i,j),ptop(2,3,n,i,j),
                             ptop(3,1,n,i,j),ptop(3,2,n,i,j),ptop(3,3,n,i,j)};
          repeat(inject) {
            compute_injection(uv,ntop(n,i,j),uflow,press,-3,2,1,rng(0));
            double gam = sqrt(1+uv[0]*uv[0]+uv[1]*uv[1]+uv[2]*uv[2]);
            double vz = uv[2]/gam;

            x = grid->x0 + hx*(i-1) + hx*uniform(rng(0),0,1);
            y = grid->y0 + hy*(j-1) + hy*uniform(rng(0),0,1);
            //z = grid->z1;
            z = grid->z1 + 0.0*vz*dt;
            age=0;
            inject_particle(species, x, y, z, uv[0], uv[1], uv[2], q(n), age, 1 );
          }
        }
      }
    }
  }  // end top injector

  //  Inject particles on Bottom Boundary

  if (global->bottom) {
    for ( int n=1; n<=nsp; n++ ) {
      species_t * species = find_species_id(n-1,species_list );
      for ( int i=1;i<=nx; i++ ) {
        for ( int j=1;j<=ny; j++ ) {
          vtherm = sqrt(2.0*pbot(3,3,n,i,j)/nbot(n,i,j));
          vd = ubot(3,n,i,j)/vtherm;
          bbot(n,i,j) += dt*nbot(n,i,j)*vtherm*(exp(-vd*vd)/sqpi+vd*(erf(vd)+1))/(2*hz);

          inject = (int) bbot(n,i,j);
          bbot(n,i,j) = bbot(n,i,j)- (double) inject;
          double uflow[3] = {0,0,ubot(3,n,i,j)};
          double press[9] = {pbot(1,1,n,i,j),pbot(1,2,n,i,j),pbot(1,3,n,i,j),
                             pbot(2,1,n,i,j),pbot(2,2,n,i,j),pbot(2,3,n,i,j),
                             pbot(3,1,n,i,j),pbot(3,2,n,i,j),pbot(3,3,n,i,j)};
          repeat(inject) {
            compute_injection(uv,nbot(n,i,j),uflow,press,3,2,1,rng(0));
            double gam = sqrt(1+uv[0]*uv[0]+uv[1]*uv[1]+uv[2]*uv[2]);
            double vz = uv[2]/gam;
            x = grid->x0 + hx*(i-1) + hx*uniform(rng(0),0,1);
            y = grid->y0 + hy*(j-1) + hy*uniform(rng(0),0,1);
            //z = grid->z0;
            z = grid->z0 + 0.0*dt*vz;
            age = 0;
            inject_particle(species, x, y, z, uv[0], uv[1], uv[2], q(n), age, 1 );
          }
        }
      }
    }
  } // end bottom injector

  //  *******  Update the injector moments at every sort interval *********

  double v[3];
  double u[3];
  double p[9];

#define icell(i,j,k) INDEX_FORTRAN_3(i,j,k,0,nx+1,0,ny+1,0,nz+1)
#define v(i) v[INDEX_FORTRAN_1(i,1,3)]
#define u(i) u[INDEX_FORTRAN_1(i,1,3)]
#define p(i,j) p[INDEX_FORTRAN_2(i,j,1,3,1,3)]

  // Parameters for measuring moments  - BE CAREFUL - don't make too big or we will go off the node

  int noff = 1; // Offset from edge - to measure moments
                // noff = 0  --> start with cell directly on boundary

  int nav = 4;   // How many cells to include in the "inward" direction in the averging
  int navin = 4; // How many cells to include in the "inward" direction in the averging

  // Right boundary Moments

  if (global->right) {
    for ( int n=1; n<=nsp; n++ ) {
      species_t * species = find_species_id(n-1,species_list );
      particle_t * part;
      if (remainder(step(), global->sort[n-1]) == 0) {
        double npart;
        for ( int k=1;k<=nz; k++ ) {
          for ( int j=1;j<=ny; j++ ) {
            npart = 0;
            flux = 0;
            u[0] = u[1] = u[2] = 0;
            p[0] = p[1] = p[2] = p[3]= p[4] = p[5] = p[6]= p[7] = p[8] = 0;
            for ( int i=nx-noff; i>nx-nav-noff; i-- ){
              int nstart = species->partition[icell(i,j,k)];
              int nstop  = species->partition[icell(i,j,k)+1];
              int ncell  = nstop - nstart;
              npart = npart + ncell;
              for (int np=nstart; np<nstop ; np++) {
                part=&species->p[np];
                double gamma = sqrt(1.0+part->ux*part->ux +
                    part->uy*part->uy + part->uz*part->uz);
                v(1) = part->ux;
                v(2) = part->uy;
                v(3) = part->uz;
                if (v(1) < 0) flux = flux - v(1)/gamma;
                for ( int a=1;a<=3; a++ ) {
                  u(a) = u(a) + v(a);
                  for ( int b=1;b<=3; b++ ) {
                    p(a,b) = p(a,b) + v(a)*v(b);
                  }
                }
              } // end particle loop for single cell
            } // end cells included for these moments

            if ( npart > 0 ) {
              nright(n,k,j) = (1.0-rout[0])*nright(n,k,j) + rout[0]*npart/nav;
              for ( int a=1; a<=3; a++ ) {
                uright(a,n,k,j) = (1-rout[1])*uright(a,n,k,j) + rout[1]*u(a)/npart;
                for ( int b=1;b<=3; b++ ) {
                  p(a,b) = (p(a,b) - u(a)*u(b)/npart)/nav;
                  pright(a,b,n,k,j) = (1-rout[2])*pright(a,b,n,k,j) + rout[2]*p(a,b);
                }
              }
            }

          }
        }
      }
    }
  }  // end right moment update

  // Left boundary Moments

  if (global->left) {
    for ( int n=1; n<=nsp; n++ ) {
      species_t * species = find_species_id(n-1,species_list );
      particle_t * part;
      if (remainder(step(), global->sort[n-1]) == 0) {
        double npart;
        for ( int k=1;k<=nz; k++ ) {
          for ( int j=1;j<=ny; j++ ) {
            npart = 0;
            flux = 0;
            u[0] = u[1] = u[2] = 0;
            p[0] = p[1] = p[2] = p[3]= p[4] = p[5] = p[6]= p[7] = p[8] = 0;
            for ( int i=1+noff;i<=nav+noff; i++ ) {
              int nstart = species->partition[icell(i,j,k)];
              int nstop  = species->partition[icell(i,j,k)+1];
              int ncell  = nstop - nstart;
              npart = npart + ncell;
              for (int np=nstart; np < nstop ; np++) {
                part=&species->p[np];
                double gamma = sqrt(1.0+part->ux*part->ux +
                    part->uy*part->uy + part->uz*part->uz);
                v(1) = part->ux;
                v(2) = part->uy;
                v(3) = part->uz;
                if (v(1) > 0) flux = flux + v(1)/gamma;
                for ( int a=1; a<=3; a++ ) {
                  u(a) = u(a) + v(a);
                  for ( int b=1;b<=3; b++ ) {
                    p(a,b) = p(a,b) + v(a)*v(b);
                  }
                }
              } // end particle loop for single cell
            } // end cells included for these moments

            if ( npart > 0 ) {
              nleft(n,k,j) = (1.0-rout[0])*nleft(n,k,j) + rout[0]*npart/nav;
              for ( int a=1; a<=3; a++ ) {
                uleft(a,n,k,j) = (1-rout[1])*uleft(a,n,k,j) + rout[1]*u(a)/npart;
                for ( int b=1;b<=3; b++ ) {
                  p(a,b) = (p(a,b) - u(a)*u(b)/npart)/nav;
                  pleft(a,b,n,k,j) = (1-rout[2])*pleft(a,b,n,k,j) + rout[2]*p(a,b);
                }
              }
            }

          }
        }
      }
    }
  }  // end left moment update

  // Top boundary Moments

  if (global->top) {
    for ( int n=1; n<=nsp; n++ ) {
      species_t * species = find_species_id(n-1,species_list );
      particle_t * part;
      if (remainder(step(), global->sort[n-1]) == 0) {
        double npart;
        for ( int i=1;i<=nx; i++ ) {
          for ( int j=1;j<=ny; j++ ) {
            npart = 0;
            flux = 0;
            u[0] = u[1] = u[2] = 0;
            p[0] = p[1] = p[2] = p[3]= p[4] = p[5] = p[6]= p[7] = p[8] = 0;
            for ( int k=nz-noff; k>nz-navin-noff; k-- ) {
              int nstart = species->partition[icell(i,j,k)];
              int nstop  = species->partition[icell(i,j,k)+1];
              int ncell  = nstop - nstart;
              npart = npart + ncell;
              for (int np=nstart; np < nstop ; np++) {
                part=&species->p[np];
                double gamma = sqrt(1.0+part->ux*part->ux +
                    part->uy*part->uy + part->uz*part->uz);
                v(1) = part->ux;
                v(2) = part->uy;
                v(3) = part->uz;
                if (v(3) < 0) flux = flux - v(3)/gamma;
                for ( int a=1; a<=3; a++ ) {
                  u(a) = u(a) + v(a);
                  for ( int b=1;b<=3; b++ ) {
                    p(a,b) = p(a,b) + v(a)*v(b);
                  }
                }
              } // end particle loop for single cell
            } // end cells included for these moments

            if ( npart > 0 ) {
              ntop(n,i,j) = (1-rin[0])*ntop(n,i,j) + rin[0]*npart/navin;
              for ( int a=1;a<=3; a++ ) {
                utop(a,n,i,j) = (1-rin[1])*utop(a,n,i,j) + rin[1]*u(a)/npart;
                for ( int b=1;b<=3; b++ ) {
                  p(a,b) = (p(a,b) - v(a)*v(b)/npart)/navin;
                  ptop(a,b,n,i,j) = (1-rin[2])*ptop(a,b,n,i,j) + rin[2]*p(a,b);
                }
              }
            }
          }

        }
      }
    }
  }  // end bottom moment update

  if (global->bottom) {
    for ( int n=1; n<=nsp; n++ ) {
      species_t * species = find_species_id(n-1,species_list );
      particle_t * part;
      if (remainder(step(), global->sort[n-1]) == 0) {
        double npart;
        for ( int i=1;i<=nx; i++ ) {
          for ( int j=1;j<=ny; j++ ) {
            npart = 0;
            flux = 0;
            u[0] = u[1] = u[2] = 0;
            p[0] = p[1] = p[2] = p[3]= p[4] = p[5] = p[6]= p[7] = p[8] = 0;
            for ( int k=1+noff;k<=navin+noff; k++ ) {
              int nstart = species->partition[icell(i,j,k)];
              int nstop  = species->partition[icell(i,j,k)+1];
              int ncell  = nstop - nstart;
              npart = npart + ncell;
              for (int np=nstart; np < nstop ; np++) {
                part=&species->p[np];
                double gamma = sqrt(1.0+part->ux*part->ux +
                    part->uy*part->uy + part->uz*part->uz);
                v(1) = part->ux;
                v(2) = part->uy;
                v(3) = part->uz;
                if (v(3) > 0) flux = flux + v(3)/gamma;
                for ( int a=1; a<=3; a++ ) {
                  u(a) = u(a) + v(a);
                  for ( int b=1;b<=3; b++ ) {
                    p(a,b) = p(a,b) + v(a)*v(b);
                  }
                }
              } // end particle loop for single cell
            } // end cells included for these moments

            if ( npart > 0) {
              nbot(n,i,j) = (1-rin[0])*nbot(n,i,j) + rin[0]*npart/navin;
              for ( int a=1;a<=3; a++ ) {
                ubot(a,n,i,j) = (1-rin[1])*ubot(a,n,i,j) + rin[1]*u(a)/npart;
                for ( int b=1;b<=3; b++ ) {
                  p(a,b) = (p(a,b) - v(a)*v(b)/npart)/navin;
                  pbot(a,b,n,i,j) = (1-rin[2])*pbot(a,b,n,i,j) + rin[2]*p(a,b);
                }
              }
            }

          }
        }
      }
    }
  }  // end bottom moment update

  //  Periodically save injector moments on outflow boundaries
  //  Only do this on the outflow boundaries, and only if we have
  //  a single domain on each boundary - otherwise would have to combine data files.

  if ( global->topology_y == 1 && global->topology_z ==1 ) {

    // How often to write moments to file -
    int nskip = 10;

    if (global->left) {
      int j = 1;
      for ( int n=1; n<=nsp; n++ ) {
        if (remainder(step(), nskip*global->sort[n-1]) == 0) {
          char buffer[20];
          sprintf(buffer, "injectors/left%i.dat", n);
          FileIO fileIO;
          FileIOStatus status;
          status= fileIO.open(buffer, io_append);
          for ( int k=1;k<=nz; k++ ) {
            zcell = (grid->z0 + k*hz-hz/2)/(global->L_de);
            fileIO.print("%6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g \n",
                zcell, nleft(n,k,j), uleft(1,n,k,j), uleft(2,n,k,j), uleft(3,n,k,j),
                pleft(1,1,n,k,j), pleft(2,2,n,k,j), pleft(3,3,n,k,j),
                pleft(1,2,n,k,j), pleft(1,3,n,k,j), pleft(2,3,n,k,j));
          } //end for
          fileIO.print("  \n \n");
          fileIO.close();
        }
      }
    }  //end left output

    if (global->right) {
      int j = 1;
      for ( int n=1; n<=nsp; n++ ) {
        if (remainder(step(), nskip*global->sort[n-1]) == 0) {
          char buffer[20];
          sprintf(buffer, "injectors/right%i.dat", n);
          FileIO fileIO;
          FileIOStatus status;
          status= fileIO.open(buffer, io_append);
          for ( int k=1;k<=nz; k++ ) {
            zcell = (grid->z0 + k*hz-hz/2)/(global->L_de);
            fileIO.print("%6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g \n",
                zcell, nright(n,k,j), uright(1,n,k,j), uright(2,n,k,j), uright(3,n,k,j),
                pright(1,1,n,k,j), pright(2,2,n,k,j), pright(3,3,n,k,j),
                pright(1,2,n,k,j), pright(1,3,n,k,j), pright(2,3,n,k,j));
          } //end for
          fileIO.print("  \n \n");
          fileIO.close();
        }
      }
    }  //end right output
  } //end output for injector moments

}

//  Use Fortran Indexing of arrays
#define vel(n)    vel[INDEX_FORTRAN_1(n,1,3)]
#define uv(n)     uv[INDEX_FORTRAN_1(n,1,3)]
#define p(a,b)    p[INDEX_FORTRAN_2(a,b,1,3,1,3)]
#define beta(a,b) beta[INDEX_FORTRAN_2(a,b,1,3,1,3)]

void compute_injection(double uv[],    // relativistic momentum of injected particle
                       double ne,      // density of boundary cell
                       double vel[],   // fluid momentum of boundary cell
                       double p[],     // relativistic stress tensor of boundary cell
                       int i,          // inward normal direction for injector
                       int j,          // j x k = n
                       int k,          //
                       rng_t *rng) {

  // Declarations
  const double sqpi =1.772453850905516;
  const int itmax = 100;
  const double convergence=1e-7;
  double beta[9];
  double r0,arg,vold,vn,vnold,vmid,fmid,error,f,dfdx;
  int n,d;

  // Normal direction

  n = abs(i);
  d = i/n;

//  Compute inverse of pressure tensor

  double det=p(1,1)*p(2,2)*p(3,3)+2*p(1,2)*p(2,3)*p(1,3)-p(1,1)*p(2,3)*p(2,3)-p(2,2)*p(1,3)*p(1,3)-p(3,3)*p(1,2)*p(1,2);
  beta(1,1) = (p(2,2)*p(3,3)-p(2,3)*p(2,3))/det;
  beta(1,2) = (p(2,3)*p(1,3)-p(1,2)*p(3,3))/det;
  beta(2,1) = beta(1,2);
  beta(1,3) = (p(1,2)*p(2,3) - p(2,2)*p(1,3))/det;
  beta(3,1) = beta(1,3);
  beta(2,2) = (p(1,1)*p(3,3)-p(1,3)*p(1,3))/det;
  beta(2,3) = (p(1,2)*p(1,3)-p(1,1)*p(2,3))/det;
  beta(3,2) = beta(2,3);
  beta(3,3) = (p(1,1)*p(2,2)-p(1,2)*p(1,2))/det;

// Check if det is less than zero

  if (det < 0) {
    ERROR((" *** Determinant is Negative  det= %e",det));
    det=p(1,1)*p(2,2)*p(3,3);
    beta(1,1) = 1/p(1,1);
    beta(1,2) = 0;
    beta(2,1) = 0;
    beta(1,3) = 0;
    beta(3,1) = 0;
    beta(2,2) = 1/p(2,2);
    beta(2,3) = 0;
    beta(3,2) = 0;
    beta(3,3) = 1/p(3,3);
  }

  // Convert beta matrix to correct units  (velocity**2/c**2)

  for ( int a=1;a<=3; a++ ) {
    for ( int b=1;b<=3; b++ ) {
      beta(a,b) = beta(a,b)*ne/2;
    }
  }

  // Specify fluid drift velocity in the normal direction

  double vtherm = sqrt(2*p(n,n)/ne);
  double vdrift = d*vel(n)/vtherm;

  // Now compute the normal velocity for this particle

 reject:

  // Bracket the root over a specified range - if the root is outside this range
  // then pick a new random number

  double v1 = 0;
  double v2 = 1/vtherm;    // Pick initial range to be some fraction of c
  double f1=1;
  double f2=1;
  while (f1*f2 > 0) {
    r0 = frand(rng);
    arg=r0*(exp(-vdrift*vdrift)+sqpi*vdrift*(1+erf(vdrift)))-exp(-vdrift*vdrift)-sqpi*vdrift*erf(vdrift);
    f1=sqpi*vdrift*erf(v1-vdrift)-exp(-(v1-vdrift)*(v1-vdrift))-arg;
    f2=sqpi*vdrift*erf(v2-vdrift)-exp(-(v2-vdrift)*(v2-vdrift))-arg;
  }

  // First attempt Newton solve - but fall back to bisection if needed

  vn=sqrt(-log(1-r0)) ; //  ! Initial guess (analytic solution for vdrift=0)

  // Now use Newton's

  vnold=0;
  error=1;
  int it=0;
  while (error > convergence && it < itmax && vn > v1 && vn < v2) {
    f=sqpi*vdrift*erf(vn-vdrift)-exp(-(vn-vdrift)*(vn-vdrift))-arg;
    dfdx=2*vn*exp(-(vn-vdrift)*(vn-vdrift));
    vn=vn-f/dfdx;
    error = fabs(vn-vnold);
    vnold=vn;
    it=it+1;
  }

  // Fall back to bisection if we did not converge

  if (error > convergence) {
    it = 0;
    vmid = (v1+v2)/2;
    while( error > convergence && it < itmax) {
      it = it+1;
      vold = vmid;
      fmid=sqpi*vdrift*erf(vmid-vdrift)-exp(-(vmid-vdrift)*(vmid-vdrift))-arg;
      if ( fmid*f2 <= 0 )
  v1 = vmid;
      else {
  v2 = vmid;
  f2 = fmid;
      }
      vmid = (v1+v2)/2;
      error = fabs(vmid-vold);
    }
    vn=vmid;
  }

  if (it > itmax){
    MESSAGE((" Convergence Problem in Injector --> error=%",error));
  }

  // error=sqpi*vdrift*erf(vn-vdrift)-exp(-(vn-vdrift)*(vn-vdrift))-arg;
  //if (error >= convergence) MESSAGE(("Oh Shit!  How did this happen -> %g", error));

  // Transfer solution to correct normal velocity component

  uv(n) = d*vtherm*vn;

  // Now compute the two transverse components of momentum

  det=(beta(j,j)*beta(k,k)-beta(j,k)*beta(j,k));
  double xj=drandn(rng)/sqrt(2);
  double xk=drandn(rng)/sqrt(2);
  uv(j) = vel(j) + xj*sqrt(beta(k,k)/det)+(uv(n)-vel(n))*p(n,j)/p(n,n);
  uv(k) = vel(k) + (xk*sqrt(beta(k,k))-(uv(n)-vel(n))*beta(n,k)-(uv(j)-vel(j))*beta(j,k))/beta(k,k);

  // To make this relativistically correct - use rejection to decide if we keep th

  double gamma = sqrt(1+uv(1)*uv(1)+uv(2)*uv(2)+uv(3)*uv(3));
  double test=drand(rng);
  if (test > 1/gamma)  goto reject;

  //  Now return relativistic momentum

} //end of function

#undef vel
#undef uv
#undef p
#undef beta

//  Define other useful macros for open boundary model

//  Safe allocate

#define ALLOCATE(A,LEN,TYPE)                                             \
  if ( !((A)=(TYPE *)malloc((size_t)(LEN)*sizeof(TYPE))) ) ERROR(("Cannot allocate."));

#define DUMP_INJECTOR(side,dim1,dim2,set)                          \
  {                                                                \
    FileIO fileIO;                                                 \
    FileIOStatus status;                                           \
    char filename[20];                                             \
    sprintf(filename, "restart" #set "/" #side ".%d",int(rank())); \
    status= fileIO.open(filename, io_write);                       \
    size_t len = dim1*dim2*nsp;                                    \
    fileIO.write(global->n##side,len);                             \
    fileIO.write(global->b##side,len);                             \
    fileIO.write(global->f##side,len);                             \
    fileIO.write(global->u##side,3*len);                           \
    fileIO.write(global->p##side,9*len);                           \
    fileIO.close();                                                \
  }

#define READ_INJECTOR(side,dim1,dim2,set)                          \
  {                                                                \
    FileIO fileIO;                                                 \
    FileIOStatus status;                                           \
    char filename[20];                                             \
    sprintf(filename, "restart" #set "/" #side ".%d",int(rank())); \
    status= fileIO.open(filename, io_read);                        \
    size_t len = dim1*dim2*nsp;                                    \
    fileIO.read(global->n##side,len);                              \
    fileIO.read(global->b##side,len);                              \
    fileIO.read(global->f##side,len);                              \
    fileIO.read(global->u##side,3*len);                            \
    fileIO.read(global->p##side,9*len);                            \
    fileIO.close();                                                \
  }

#define DEFINE_INJECTOR(side,dim1,dim2)                            \
  size_t len = dim1*dim2*nsp;                                      \
  ALLOCATE(global->f##side,len,double);                            \
  ALLOCATE(global->n##side,len,double);                            \
  ALLOCATE(global->b##side,len,double);                            \
  ALLOCATE(global->u##side,3*len,double)                           \
  ALLOCATE(global->p##side,9*len,double);


# define DUMP_INJECTORS(set)                                       \
  if (global->left) DUMP_INJECTOR(left,ny,nz,set);                 \
  if (global->right) DUMP_INJECTOR(right,ny,nz,set);               \
  if (global->top) DUMP_INJECTOR(top,ny,nx,set);                   \
  if (global->bottom) DUMP_INJECTOR(bot,ny,nx,set);

  // Define Fortran style indexing of arrays

#define INDEX_FORTRAN_4(a,b,c,d,al,ah,bl,bh,cl,ch,dl,dh)  \
  ((a)-(al) + ((ah)-(al)+1)*((b)-(bl) +  ((bh)-(bl)+1)*(((c)-(cl)) + ((ch)-(cl)+1)*((d)-(dl)))))

#define INDEX_FORTRAN_5(a,b,c,d,e,al,ah,bl,bh,cl,ch,dl,dh,el,eh)    \
  ((a)-(al) + ((ah)-(al)+1)*((b)-(bl) +  ((bh)-(bl)+1)*(((c)-(cl)) + ((ch)-(cl)+1)*(((d)-(dl)) + ((dh)-(dl)+1)*((e)-(el))))))

#define bbot(n,i,j) global->bbot[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nx,1,ny)]
#define fbot(n,i,j) global->fbot[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nx,1,ny)]
#define nbot(n,i,j) global->nbot[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nx,1,ny)]
#define ubot(a,n,i,j) global->ubot[INDEX_FORTRAN_4(a,n,i,j,1,3,1,nsp,1,nx,1,ny)]
#define pbot(a,b,n,i,j) global->pbot[INDEX_FORTRAN_5(a,b,n,i,j,1,3,1,3,1,nsp,1,nx,1,ny)]

#define btop(n,i,j) global->btop[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nx,1,ny)]
#define ftop(n,i,j) global->ftop[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nx,1,ny)]
#define ntop(n,i,j) global->ntop[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nx,1,ny)]
#define utop(a,n,i,j) global->utop[INDEX_FORTRAN_4(a,n,i,j,1,3,1,nsp,1,nx,1,ny)]
#define ptop(a,b,n,i,j) global->ptop[INDEX_FORTRAN_5(a,b,n,i,j,1,3,1,3,1,nsp,1,nx,1,ny)]

#define bleft(n,i,j) global->bleft[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nz,1,ny)]
#define fleft(n,i,j) global->fleft[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nz,1,ny)]
#define nleft(n,i,j) global->nleft[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nz,1,ny)]
#define uleft(a,n,i,j) global->uleft[INDEX_FORTRAN_4(a,n,i,j,1,3,1,nsp,1,nz,1,ny)]
#define pleft(a,b,n,i,j) global->pleft[INDEX_FORTRAN_5(a,b,n,i,j,1,3,1,3,1,nsp,1,nz,1,ny)]

#define bright(n,i,j) global->bright[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nz,1,ny)]
#define fright(n,i,j) global->fright[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nz,1,ny)]
#define nright(n,i,j) global->nright[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nz,1,ny)]
#define uright(a,n,i,j) global->uright[INDEX_FORTRAN_4(a,n,i,j,1,3,1,nsp,1,nz,1,ny)]
#define pright(a,b,n,i,j) global->pright[INDEX_FORTRAN_5(a,b,n,i,j,1,3,1,3,1,nsp,1,nz,1,ny)]

#define uf(n)   global->uf[INDEX_FORTRAN_1(n,1,nsp)]
#define vth(n) global->vth[INDEX_FORTRAN_1(n,1,nsp)]
#define vthb(n) global->vthb[INDEX_FORTRAN_1(n,1,nsp)]
#define q(n) global->q[INDEX_FORTRAN_1(n,1,nsp)]

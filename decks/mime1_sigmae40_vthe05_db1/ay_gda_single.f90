!---------------------------------------------------------------------------------------
!  This program computes the flux surfaces - Ay, by reading the the inplane magnetic field
!  components from data/bx.gda and data/bz.gda, taking the curl, and solving the needed
!  Poisson equation
!
!  This version is for a single sheet - periodic in x, with conducting BC in z
!---------------------------------------------------------------------------------------
 

module fft
  implicit none
  integer(kind=8) forward,inverse
  real(kind=8) norm,kx,pi
  real(kind=8), allocatable, dimension(:) :: kx2,in,out
  INTEGER FFTW_PATIENT
  PARAMETER (FFTW_PATIENT=32)
  INTEGER FFTW_R2HC
  PARAMETER (FFTW_R2HC=0)
  INTEGER FFTW_HC2R
  PARAMETER (FFTW_HC2R=1)
  INTEGER FFTW_RODFT00
  PARAMETER (FFTW_RODFT00=7)
  INTEGER FFTW_RODFT01
  PARAMETER (FFTW_RODFT01=8)
  INTEGER FFTW_RODFT10
  PARAMETER (FFTW_RODFT10=9)
  INTEGER FFTW_REDFT01
  PARAMETER (FFTW_REDFT01=4)
  INTEGER FFTW_REDFT10
  PARAMETER (FFTW_REDFT10=5)
end module fft


program translate
  use fft
  implicit none
  integer input_record,output_record,input_error,output_error,record_length,it1,it2,i,j,it,interval
  integer(kind=4) :: nx,ny,nz
  real(kind=4)diff,xmax,ymax,zmax,time
  real(kind=4), allocatable, dimension(:,:) :: bx,bz,Ay,curlb
  real(kind=8) dx,dz
  character *6 cit

! Range of time slices to convert.  Make big and let error terminate program
  interval = 1571
  it1=0
  it2=60*interval
  output_record = 0

! Open data information file to find out the size of the arrays

  open(unit=10,file='data/info',access='stream',status='old',form='unformatted')

!  Read data from info  file

  read(10)nx,ny,nz
  read(10)xmax,ymax,zmax 

! Cell size

  dx = xmax/real(nx)
  dz = zmax/real(nz)

! Echo this information

  print *,"---------------------------------------------------"
  print *,"xmax=",xmax,"    zmax=",zmax
  print *,"nx=",nx,"   nz=",nz
  print *,"dx=",dx,"   dz=",dz
  print *,"---------------------------------------------------"

! Initialize FFT library
    
  allocate(in(nx))
  allocate(out(nx))
  norm = 1.0d0/sqrt(real(nx))
  call dfftw_plan_r2r_1d(forward,nx,in,out,FFTW_R2HC,FFTW_PATIENT)
  call dfftw_plan_r2r_1d(inverse,nx,in,out,FFTW_HC2R,FFTW_PATIENT)

! Wavectors in x direction ( Half-Complex Format )
! See page 10 of FFTW-3.0 manual for description

  allocate(kx2(nx))
  pi = acos(-1.0d0)
  kx2(1)=0.0d0
  do j=1,nx/2 
     kx = (2.0d0*pi*real(j)/xmax) 
     diff = sin(kx*dx/2.0d0)/(kx*dx/2.0d0)
     kx2(j+1) = (kx*diff)**2 
  enddo
  do j=1,(nx+1)/2-1
     kx = (2.0d0*pi*real(j)/xmax) 
     diff = sin(kx*dx/2.0d0)/(kx*dx/2.0d0)
     kx2(nx+1-j) = (kx*diff)**2 
  enddo

! Allocate storage space for fields and moments

  allocate(bx(nx,nz))
  allocate(bz(nx,nz))
  allocate(Ay(nx,nz))
  allocate(curlb(nx,nz))

! Set record length for gda files

! Use this for Bill's gda format with extra time and it records
!  inquire(iolength=record_length)bx,time,it

! Use this if you are only saving the matrix in the gda file
  inquire(iolength=record_length)bx

! *** WARNING - also make the read statments below consistent ***

  print *," Setting record length=",record_length

! Open all of the direct access binary files (Bill's preferred data format)
! Loop over time slices

  ! do input_record = it1,it2,interval
  do input_record = it1,it2
     write(cit,'(I0)') input_record

     ! open(unit=20,file='data/bx_'//cit//'.gda',access='direct',&
     !    recl=record_length,status='unknown',form='unformatted',action='read') 
     ! open(unit=30,file='data/bz_'//cit//'.gda',access='direct',&
     !    recl=record_length,status='unknown',form='unformatted',action='read')    

     open(unit=20,file='data/bx.gda',access='direct',&
        recl=record_length,status='unknown',form='unformatted',action='read')     
     open(unit=30,file='data/bz.gda',access='direct',&
        recl=record_length,status='unknown',form='unformatted',action='read')    

     print *,"Reading record=",input_record

! Read magnetic field data

     ! read(20,rec=input_record+1)bx,time,it
     ! read(30,rec=input_record+1)bz,time,it
     read(20, rec=input_record+1)bx
     read(30, rec=input_record+1)bz
  
! Compute Curl of magnetic field

     print *,"Computing curl"

        do i=2,nx-1
           do j=2,nz-1
              curlb(i,j) = (bx(i,j+1)-bx(i,j-1))/(2.0*dz) - (bz(i+1,j)-bz(i-1,j))/(2.0*dx)
!              curlb(i,j) = (bx(i,k,j+1)-bx(i,k,j-1))/(2.0*dz) - (bz(i+1,k,j)-bz(i-1,k,j))/(2.0*dx)
           enddo
        enddo
        curlb(1,:) = curlb(2,:)
        curlb(nx,:) = curlb(nx-1,:)
        curlb(:,1) = curlb(:,nz-1)
        curlb(:,nz) = curlb(:,2)

! Solve Poisson equation

        print *,"Solving Poisson"
       call poisson(dx,nx,dz,nz,curlb,Ay,bx)

! Save flux surfaces

       print *,"Saving Ay"
       output_record = output_record +1
       open(unit=40,file='data/Ay.gda',access='direct', &
        recl=record_length,status='unknown',form='unformatted',action='write')     
       !write(40,rec=output_record)Ay,time,it
       write(40,rec=output_record)Ay
       close(40)

    close(20)
    close(30)


  enddo

! Close files


end program translate

subroutine poisson(dx,nx,dz,nz,rho,phi,bx)
  use fft
  implicit none
  integer nx,nz,i,j
  real(kind=4) phi(nx,nz),rho(nx,nz),bx(nx,nz)
  real(kind=8) dx,dz,dx2,dz2,ak(nx,nz)
  real(kind=8) d(nz),g(nz),w(nz),rhs(nz),cycle(nx,nz)
  real(kind=8) error,maxerror,delsq,bsum

! Do internal calculation in double precision to avoid round-off on large grids

! Constants

dz2 = dz**2
dx2 = dx**2

! fft in x-direction
  do j=1,nz
     in(:) = -dz2*real(rho(:,j),kind=8)
     call dfftw_execute(forward)
     ak(:,j) = norm*out(:)
  enddo

! Now do a tri-diagonal solve across the z-direction

do i=1,nx
   d(:) = -2.0d0 - dz2*kx2(i)
   d(1) = -3.0d0 - dz2*kx2(i)  ! phi=0 at z=0
   d(nz) = -3.0d0 - dz2*kx2(i) ! phi = 0 at z=zmax - fix this below
   w(1) = 1.0d0/d(1)
   g(1) = ak(i,1)/d(1)
   do j=2,nz
      w(j) = 1.0d0/(d(j)-w(j-1))
      g(j) = (ak(i,j)-g(j-1))*w(j)
   enddo
   ak(i,nz) = g(nz)
   do j=nz-1,1,-1
      ak(i,j) = g(j) - w(j)*ak(i,j+1)
   enddo
enddo

! Now take inverse FFT

do j=1,nz
   in(:) = ak(:,j)
   call dfftw_execute(inverse)
   ak(:,j) = norm*out(:)
enddo

! Now add correction for oblique or asymmetric layers

  bsum = 0.0d0
  do i=1,nx
     do j=1,nz
        bsum = bsum + bx(i,j)
     enddo
  enddo
  bsum = -bsum/real(nz*nx)
  
  do i=1,nx
     do j=1,nz
        ak(i,j) = ak(i,j) + bsum*real(j-nz/2)*dz
     enddo
  enddo

! Do error check on solution

maxerror = 0.0
do i=2,nx-1
   do j=2,nz-1
      delsq = (ak(i+1,j)-2.0*ak(i,j)+ak(i-1,j))/dx2 + (ak(i,j+1)-2.0*ak(i,j)+ak(i,j-1))/dz2
      error = abs(real(rho(i,j),kind=8) + delsq)
      if (error > maxerror) maxerror=error
   enddo
enddo
print *,"Error in poisson solve = ", maxerror

! Return solution in single precision

  phi(:,:) = real(ak(:,:),kind=4)

  return
end subroutine poisson





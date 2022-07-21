!---------------------------------------------------------
!
!  Parallel Diagnostic for Computing Turbelent Spectra 
!
!  This version uses P3FFT package which peforms 
!  3D parallel Fourier transforms.  Note that the decomposion
!  must be in the y and z direction (not in x)
!--------------------------------------------------------

module MPI
  include "mpif.h"                                                                                         
  integer myid,numprocs,ierr                                                                               
  integer master  
  integer nfiles
  parameter(nfiles=1)
  integer sizes(3), subsizes(3), starts(3)
  integer fileinfo, ierror, fh(nfiles), filetype, status(MPI_STATUS_SIZE), output_format
  integer(kind=MPI_OFFSET_KIND) :: disp, offset
  parameter(master=0)                                                                                      
end module MPI

program turbulence
  use mpi
  use p3dfft
  implicit none
  character(50) fname,sname,time
  character(4), allocatable, dimension(:) ::  bname
  real(kind=4) xmax,ymax,zmax,pi,a0,a1,a2,alpha,window,dkx,dky,dkz,kperp,dk,kmax
  parameter(pi=3.141592654)
  integer(kind=4)nxt,nyt,nzt,ntot,i,j,k,n,i0,j0,k0,num,nbase,ip
  real(kind=4), allocatable, dimension(:,:) ::  lper,perp
  real(kind=4), allocatable, dimension(:) ::  lpar
  real(kind=4), allocatable, dimension(:,:) :: para,eperp
  logical found

! p3fft stuff

  integer istart(3),iend(3),isize(3)
  integer fstart(3),fend(3),fsize(3)
  logical overwrite
  integer(kind=4) topology(2)
  integer, parameter :: mytype=4
  real(mytype), allocatable, dimension(:,:,:) :: field
  complex(mytype), allocatable, dimension(:,:,:) :: local                                              

! Initialize MPI 

  call MPI_INIT(ierr)                                                                                      
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)                                                           
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)         

! Displacement and offset

  disp = 0 
  offset = 0

! Time slice

  time='14400'

! Variables to transform

  nbase = 3
  allocate(bname(nbase))
  fname="b"
  bname(1) = 'bx'
  bname(2) = 'by'
  bname(3) = 'bz'

!  fname="U"
!  bname(1) = 'uix'
!  bname(2) = 'uiy'
!  bname(3) = 'uiz'

! Read array sizes from file

  if (myid == master) then
     open(unit=10,file='data/info',status='unknown',form='unformatted')
     read(10)nxt,nyt,nzt
     read(10)xmax,ymax,zmax
     print *,"nxt=",nxt,"  xmax=",xmax
     print *,"nyt=",nyt,"  ymax=",ymax
     print *,"nzt=",nzt,"  zmax=",zmax
  endif
  call MPI_BCAST(nxt,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nyt,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nzt,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(xmax,1,MPI_REAL,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ymax,1,MPI_REAL,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(zmax,1,MPI_REAL,master,MPI_COMM_WORLD,ierr)

! Alternatively - can specify these manually

!  xmax = 720.0
!  ymax = 720.0
!  zmax = 240.0
!  nxt = 3072
!  nyt = 3072
!  nzt = 1024

! Topology for parallel processing

  topology(1)=8  ! y-direction
  topology(2)=8  ! z-direction

! Is the requested topology consistent ? 

  if ( topology(1)*topology(2) .ne. numprocs) then
     print *, " ** Inconsistent Topology --> Terminating ***"
     call MPI_FINALIZE(ierr)
     stop
  endif

! Initialize P3DFFT

 if (myid == master ) print *,"Initializing P3DFFT"
 overwrite=.False.

! Syntax for 2.6.1 of P3DFFT

 call p3dfft_setup(topology,nxt,nyt,nzt,MPI_COMM_WORLD)
! Syntax for 2.4 of P3DFFT
! call p3dfft_setup(topology,nxt,nyt,nzt,overwrite)

 call p3dfft_get_dims(istart,iend,isize,1)
 call p3dfft_get_dims(fstart,fend,fsize,2)

 if (myid == master) then    
    print *,"------------------------------------------------------"
    print *,"Decomposition of real data"
    print *,"istart=",istart(1)," iend=",iend(1)," isize=",isize(1)
    print *,"istart=",istart(2)," iend=",iend(2)," isize=",isize(2)
    print *,"istart=",istart(3)," iend=",iend(3)," isize=",isize(3)
    print *
    print *,"Decomposition of FFT data"
    print *,"fstart=",fstart(1)," fend=",fend(1)," fsize=",fsize(1)
    print *,"fstart=",fstart(2)," fend=",fend(2)," fsize=",fsize(2)
    print *,"fstart=",fstart(3)," fend=",fend(3)," fsize=",fsize(3)
    print *,"------------------------------------------------------"
 endif

 call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! if (myid == numprocs-1) then   
!    print *,"------------------ LAST NODE ---------------------------" 
!    print *,"istart=",istart(1)," iend=",iend(1)," isize=",isize(1)
!    print *,"istart=",istart(2)," iend=",iend(2)," isize=",isize(2)
!    print *,"istart=",istart(3)," iend=",iend(3)," isize=",isize(3)
!    print *
!    print *,"fstart=",fstart(1)," fend=",fend(1)," fsize=",fsize(1)
!    print *,"fstart=",fstart(2)," fend=",fend(2)," fsize=",fsize(2)
!    print *,"fstart=",fstart(3)," fend=",fend(3)," fsize=",fsize(3)
!    print *,"------------------------------------------------------" 
! endif


! Allocate storage space for calculations

 if (myid .eq. master) print *, "Allocating storage arrays"
 allocate(field(isize(1),isize(2),isize(3)))
 allocate(local(fsize(1),fsize(2),fsize(3)))

! allocate(lpar(nxt/2+1))
! allocate(para(nbase,nxt/2+1))
! allocate(lper(nyt,nzt))
! allocate(perp(nyt,nzt))

 allocate(lpar(nyt))
 allocate(para(nbase,nyt))
 allocate(lper(nxt/2+1,nzt))
 allocate(perp(nxt/2+1,nzt))

! Initialize variables

  para=0.0
  perp=0.0
  lpar=0.0
  lper=0.0

! Size of the global matrix

  sizes(1) = nxt
  sizes(2) = nyt
  sizes(3) = nzt

! size of the chunck seen by each process

  subsizes(1) = isize(1)
  subsizes(2) = isize(2)
  subsizes(3) = isize(3)

! where each chunck starts

  starts(1) = istart(1)-1
  starts(2) = istart(2)-1
  starts(3) = istart(3)-1

! Grid for reduced kperp calculation

  dkx = 2.0*pi/xmax
  dky = 2.0*pi/ymax
  dkz = 2.0*pi/zmax
  num = nzt/2
  kmax = pi*real(nzt)/zmax
  dk = kmax/real(num)          
  allocate(eperp(nbase,num))
  eperp=0.0

! Loop over base names

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if (myid .eq.0) print *,"Start of main loop"

  do n=1,nbase

!  Open the data file
     write(sname,"(A,A,A,A,A)")"data/",trim(bname(n)),"_",trim(time),".gda"
     inquire(file=trim(sname),exist=found)
     if (.not. found) then
        print *,trim(sname),"  does not exist --> Terminating"
        call MPI_FINALIZE(ierr)
        stop
     endif
     if (myid .eq. 0 ) print *,"Reading file -->",trim(sname)
     call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts, MPI_ORDER_FORTRAN, MPI_REAL4, filetype, ierror)
     call MPI_TYPE_COMMIT(filetype, ierror)
     call MPI_INFO_CREATE(fileinfo,ierror)
     call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sname), MPI_MODE_RDONLY, fileinfo, fh, ierror)
     call MPI_FILE_SET_VIEW(fh,disp, MPI_REAL4, filetype, 'native', MPI_INFO_NULL, ierror)

! Read data and close file

     ntot = isize(1)*isize(2)*isize(3)
     call MPI_FILE_READ_AT_ALL(fh, offset, field, ntot, MPI_REAL4, status, ierror)
     call MPI_FILE_CLOSE(fh, ierror)
     if (myid .eq. 0 ) print *,"Finished Reading data"

! Subtract off equilibrium By field

     if ( n== 2) then
        if (myid .eq. 0 ) print *,"Remove equilibrium By"     
        field = field - 0.46
     endif

! Apply Blackman window in z-direction
   
        if (myid .eq. 0 ) print *,"Apply Blackman window in z-direction"     
        alpha=0.16
        a0 = (1.0-alpha)/2.0
        a1=0.5
        a2=alpha/2.0
        do k=1,isize(3)
           j = k+istart(3)-1
           window = a0 -a1*cos(2.0*pi*j/(nzt-1))+a2*cos(4*pi*j/(nzt-1))
           field(:,:,k) = window*field(:,:,k) 
        enddo
     
! Take 3D transform

     if (myid .eq. 0 ) print *,"Taking 3D Fourier Transform"

!  This is the syntax for version 2.6.1 of P3FFT
     call p3dfft_ftran_r2c(field,local,"fff")

!  This is the syntax for version 2.4 of P3FFT
!     call p3dfft_ftran_r2c(field,local)


     local = local/real(nxt*nyt*nzt)

! Compute reduced power spectrums on each node

     lpar(:) = 0.0
     lper(:,:) = 0.0
     do i=1,fsize(1)
        do j=1,fsize(2)
           do k=1,fsize(3)
              i0 = fstart(1)+i-1
              j0 = fstart(2)+j-1
              k0 = fstart(3)+k-1
              lpar(j0) = lpar(j0) + abs(local(i,j,k))**2
              lper(i0,k0) = lper(i0,k0) + abs(local(i,j,k))**2
           enddo
        enddo
     enddo
     
! Sum across nodes and output global power spectrum on master

     call MPI_REDUCE(lpar,para(n,:),nyt,MPI_REAL,MPI_SUM,master,MPI_COMM_WORLD,ierr)
     call MPI_REDUCE(lper,perp,(nxt/2+1)*nzt,MPI_REAL,MPI_SUM,master,MPI_COMM_WORLD,ierr)

! Compute reduced perp spectrum

     if (myid .eq. master) then

        do i = 1,nxt/2+1
           do k = 1,nzt/2
              kperp =sqrt((dkx*real(i-1))**2 + (dkz*real(k-1))**2)
              ip = int(kperp/dk)+1
              if (ip .le. num ) eperp(n,ip) = eperp(n,ip) + (perp(i,k) + perp(i,nzt-k+1))
           enddo
        enddo
     endif

! End of loop over base quantities

  enddo

  if (myid .eq. 0) then

     write(sname,"(A,A,A)")trim(fname),trim(time),".par"
     open(unit=50,file=trim(sname),status='unknown')
     do i=1,nyt/2
        write(50,66)2.0*pi*real(i)/ymax,(para(1,i)+para(1,nyt+1-i)),(para(2,i)+para(2,nyt+1-i)),  &
             (para(3,i)+para(3,nyt+1-i)), &
             (para(1,i)+para(2,i)+para(3,i)+para(1,nyt+1-i)+para(2,nyt+1-i)+para(3,nyt+1-i))
     enddo
     close(50)


     write(sname,"(A,A,A)")trim(fname),trim(time),".per"
     open(unit=50,file=trim(sname),status='unknown')
     do i=1,num
        kperp = dk*real(i)
        write(50,66)kperp,(eperp(1,i)),(eperp(2,i)),(eperp(3,i)),(eperp(1,i)+eperp(2,i)+eperp(3,i))
     enddo
     close(50)

  endif

! Format statements

 66 format(5(1x,1pe13.5))

! End of Program

  call MPI_FINALIZE(ierr)

end program turbulence


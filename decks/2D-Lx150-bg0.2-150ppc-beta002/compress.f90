!---------------------------------------------------------
!
!  Parallel Diagnostic for Computing Turbelent Spectra 
!
!  This version decomposes in the z-direction across the layer
!  and each processors does FFT in x and y directions and then
!  we finally average across the z-direction
!--------------------------------------------------------

module MPI
  include "mpif.h"                                                                                         
  integer myid,numprocs,ierr                                                                               
  integer master  
  integer nfiles, nbands
  parameter(nfiles=1)
  integer insize(3), insub(3), instart(3),inbuf
  integer outsize(3), outsub(3), outstart(3),outbuf
  integer fileinfo, ierror, fh(nfiles), filetype, status(MPI_STATUS_SIZE)
  integer(kind=MPI_OFFSET_KIND) :: disp, offset
  parameter(master=0)                 
  parameter(disp =0,offset=0)                                                                     
end module MPI

program compress
  use mpi
  implicit none
  character(50) fname,sname
  character(10), allocatable, dimension (:) ::  base
  real(kind=4) xmax,ymax,zmax,c
  integer(kind=4)nx,ny,nz,nxt,nyt,nzt,tx,ty,tz,ix,iy,iz,i,j,k,l
  integer(kind=4)m,n,nbase,tstart,delta,nb,time,it0,nslice,it,tslice
  real(kind=4), allocatable, dimension(:,:,:) :: in,out
!  integer(kind=4), allocatable, dimension(:) :: tslice
  logical found

! Initialize MPI 

  call MPI_INIT(ierr)                                                                                      
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)                                                           
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)                                                       

! Topology for parallel processing

!  tx = 20
!  ty = 20
!  tz = 4
 
  tx = 1
  ty = 1
  tz = 32

! Choose number of base names 

  nbase = 31
  allocate(base(nbase))

! Choose starting time interval

  tstart = 1
  delta = 2217
  nslice=40

! Altertnativel - pick time slices directly
! Useful if not all slices are equally space in time

!  nslice=10
!  allocate(tslice(nslice))
!  tslice(1) = 2550
!  tslice(2) = 3060
!  tslice(3) = 3570
!  tslice(4) = 4080
!  tslice(5) = 4590
!  tslice(6) = 5100
!  tslice(7) = 5610
!  tslice(8) = 6120
!  tslice(9) = 6630
!  tslice(10) = 7140

! Choose base names

  base(1) = 'ne'
  base(2) = 'ni'
  base(3) = 'bx'
  base(4) = 'by'
  base(5) = 'bz'
  base(6) = 'ex'
  base(7) = 'ey'
  base(8) = 'ez'
  base(9) = 'jx'
  base(10) = 'jy'
  base(11) = 'jz'
  base(12) = 'absJ'
  base(13) = 'uex'
  base(14) = 'uey'
  base(15) = 'uez'
  base(16) = 'uix'
  base(17) = 'uiy'
  base(18) = 'uiz'
  base(19) = 'e-mix1'
  base(20) = 'i-mix1'
  base(21) = 'eEB01'
  base(22) = 'eEB02'
  base(23) = 'eEB03'
  base(24) = 'eEB04'
  base(25) = 'eEB05'
  base(26) = 'vex'
  base(27) = 'vey'
  base(28) = 'vez'
  base(29) = 'vix'
  base(30) = 'viy'
  base(31) = 'viz'

! Read array sizes from file

  if (myid == master) then
     ! open(unit=10,file='data/info',status='unknown',form='unformatted')
     open(unit=10, file="data/info", access='stream', status='unknown', &
          form='unformatted', action='read')
     read(10)nxt,nyt,nzt
     read(10)xmax,ymax,zmax
     print *,"nxt=",nxt,"  xmax=",xmax
     print *,"nyt=",nyt,"  ymax=",ymax
     print *,"nzt=",nzt,"  zmax=",zmax
     close(10)
     ! open(unit=10,file='data-smooth/info',status='unknown',form='unformatted')
     open(unit=10, file="data-smooth/info", access='stream', status='unknown', &
          form='unformatted', action='write')
     write(10)nxt/2,nyt/2,nzt/2
     print *,"nx=",nxt/2,"  xmax=",xmax
     print *,"ny=",nyt/2,"  ymax=",ymax
     print *,"nz=",nzt/2,"  zmax=",zmax
     write(10)xmax,ymax,zmax
     close(10)
  endif
  call MPI_BCAST(nxt,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nyt,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nzt,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(xmax,1,MPI_REAL,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ymax,1,MPI_REAL,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(zmax,1,MPI_REAL,master,MPI_COMM_WORLD,ierr)

! Is the requested topology consistent ? 

  if ( tx*ty*tz .ne. numprocs) then
     print *, " ** Inconsistent Topology --> Terminating ***"
     call MPI_FINALIZE(ierr)
     stop
  endif
  if ( mod(nxt,tx) .ne. 0 .or. mod(nyt,ty) .ne. 0 .or. mod(nzt,tz) .ne. 0 ) then
     print *, " ** Number of cells (nxy,nyt,nzt) need to divide evenly with topology --> Terminating"
     call MPI_FINALIZE(ierr)
     stop
  endif

! Domain decomposition

  nx = nxt/tx
  ny = nyt/ty
  nz = nzt/tz
  call rank_to_index(myid,ix,iy,iz,tx,ty,tz) 

! Size of the global input data

  insize(1) = nxt
  insize(2) = nyt
  insize(3) = nzt

! Size of the strided output data

  outsize(1) = nxt/2
  outsize(2) = nyt/2
  outsize(3) = nzt/2

! size of the chunck seen by each process

  insub(1) = nx+1
  insub(2) = ny+1
  insub(3) = nz+1
  if (ix .eq. tx-1) insub(1) = nx
  if (iy .eq. ty-1) insub(2) = ny
  if (iz .eq. tz-1) insub(3) = nz
  inbuf = insub(1)*insub(2)*insub(3)

  outsub(1) = nx/2
  outsub(2) = ny/2
  outsub(3) = nz/2
  outbuf = outsub(1)*outsub(2)*outsub(3)

! where each chunck starts

  instart(1) = ix*nx
  instart(2) = iy*ny
  instart(3) = iz*nz
  outstart(1) = ix*nx/2
  outstart(2) = iy*ny/2
  outstart(3) = iz*nz/2

! Allocate storage 

  allocate(in(nx+1,ny+1,nz+1))
  allocate(out(outsub(1),outsub(2),outsub(3)))

! Loop over base names

!  do nb = 20,nbase
  do nb = 26, 31

     found = .true.

! Loop over time steps

!     if (nb.eq.2) it0=78
     do it = 0, 40

        time=it*delta
        
        write(sname,"(A,A,A,I0,A)")"data/",trim(base(nb)),"_",time,".gda"
        inquire(file=trim(sname),exist=found)
        if (.not. found) then
           print *,trim(sname),"  does not exist --> Skipping"
           exit
!     call MPI_FINALIZE(ierr)
!     stop
        endif
        if (myid .eq. 0 ) print *,"Reading file -->",trim(sname)
        call MPI_TYPE_CREATE_SUBARRAY(3,insize,insub,instart, MPI_ORDER_FORTRAN, MPI_REAL4, filetype, ierror)
        call MPI_TYPE_COMMIT(filetype, ierror)
        call MPI_INFO_CREATE(fileinfo,ierror)
        call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sname), MPI_MODE_RDONLY, fileinfo, fh, ierror)
        call MPI_FILE_SET_VIEW(fh,disp, MPI_REAL4, filetype, 'native', MPI_INFO_NULL, ierror)
        
! Read data and close file

        call MPI_FILE_READ_AT_ALL(fh, offset, in(1:insub(1),1:insub(2),1:insub(3)), inbuf, MPI_REAL4, status, ierror)
        call MPI_FILE_CLOSE(fh, ierror)
        if (myid .eq. 0 ) print *,"Finished Reading Data"

! Stride through input data, apply digital filter and reduce matrix size by 2 in each direction

        if (myid .eq. 0 ) print *,"Apply digital filter and stride by 2x in each direction"

! For MPI domains on the maximum edge, our smoothing operation would need a ghost cell
! which we don't really have, so I am just copying the edge into ghost
  
        if (ix .eq. tx-1) in(nx+1,:,:) = in(nx,:,:)
        if (iy .eq. ty-1) in(:,ny+1,:) = in(:,ny,:)
        if (iz .eq. tz-1) in(:,:,nz+1) = in(:,:,nz)

        do i=2,insub(1),2
           do j=2,insub(2),2
              do k=2,insub(3),2
                 out(i/2,j/2,k/2) =  0.0
                 do l=-1,1
                    do m=-1,1
                       do n=-1,1

! This should give generalization of (1,2,1) template from 1D to 3D data

                          c=8.0/(2.0**(abs(l)+abs(m)+abs(n)))
!                          if (myid .eq. 0 .and. i .eq. 4 .and. j .eq. 4 .and. k .eq. 8) print *,l,m,n,c
                          out(i/2,j/2,k/2) =  out(i/2,j/2,k/2) + c*in(i+l,j+m,k+n)
                       enddo
                    enddo
                 enddo
                 out(i/2,j/2,k/2) =  out(i/2,j/2,k/2)/64.0
              enddo
           enddo
        enddo

! Save smoothed data

        write(sname,"(A,A,A,I0,A)")"data-smooth/",trim(base(nb)),"_",time,".gda"
        if (myid .eq. 0 ) print *,"Writing file -->",trim(sname)
        call MPI_TYPE_CREATE_SUBARRAY(3,outsize,outsub,outstart, MPI_ORDER_FORTRAN, MPI_REAL4, filetype, ierror)
        call MPI_TYPE_COMMIT(filetype,ierror)
        call MPI_INFO_CREATE(fileinfo,ierror)
        call MPI_INFO_SET(fileinfo,"romio_cb_write","enable",ierror)
        call MPI_INFO_SET(fileinfo,"romio_ds_write","disable",ierror)
        call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sname), MPI_MODE_RDWR+MPI_MODE_CREATE, fileinfo, fh, ierror)
        call MPI_FILE_SET_VIEW(fh,disp, MPI_REAL4, filetype, 'native', MPI_INFO_NULL, ierror)
        call MPI_FILE_WRITE_ALL(fh, out, outbuf, MPI_REAL4, status, ierror)
        call MPI_FILE_CLOSE(fh, ierror)
        if (myid .eq. 0 ) print *,"Finished Wrting Data"
!        time = time + delta
 666 continue
     enddo

  enddo

  deallocate(in, out)
  
! End of Program

  call MPI_FINALIZE(ierr)

end program compress

subroutine rank_to_index(rank,ix,iy,iz,topology_x,topology_y,topology_z) 
implicit none
integer iix, iiy, iiz, rank,ix,iy,iz,topology_x,topology_y,topology_z

iix  = rank
iiy  = iix/topology_x
iix  = iix - iiy*topology_x
iiz  = iiy/topology_y 
iiy  = iiy - iiz*topology_y

ix = iix
iy = iiy
iz = iiz

end subroutine rank_to_index





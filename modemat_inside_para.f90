program modecoupmat_inside

!  compute the convolution matrix due to a mask inside a square patch
!   inf: input file of shear mask field [1:ni,1:nj] (real)
!   inf2: input file of halo mask field [1:ni,1:nj] (real)
!   outf: output field of mode coupling matrix of EE, BB, EB power
!          [1:npmax*3,1:npmax*3] (real)

  use fieldinfo_
  use sub_
  implicit none

  real :: lside_amin
  integer :: i,j,i1,j1,id,jd,id1,jd1,im,jm,nn(2),np1,np2,k,dj,ifac
  integer, parameter :: ikmax=2
  real, pointer :: mask(:,:),mask2(:,:)
  real ::  num_part(npmax),num(npmax),fmask(ni*2,nj),fmask2(ni*2,nj)
  real ::  l,pbl,qbl,powm,ntotpix
  double precision :: lside(2)
  real ::  labs,labs2,g1,g2,cosj,cosj2,sinj,sinj2,lx,ly,l2x,l2y
  double precision ::  mmat_part(npmax,npmax,ikmax),mmat(npmax,npmax,ikmax),dmmat(ikmax)
  double precision :: mmat3(npmax3,npmax3)
  character(150) :: inf,inf2,cform
  character(150) :: outf
  character(20) :: clside_amin

  include 'mpif.h'
  integer :: rank,size,ierror,status(MPI_STATUS_SIZE)
  integer :: jstart,jend

  call getarg(1,inf)
  call getarg(2,inf2)
  call getarg(3,outf)

  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

  print*,'node',rank
  write(*,*) inf,inf2,outf

  call binninginit
  call readmask(inf,mask)
  call readmask(inf2,mask2)

  nn(1)=ni
  nn(2)=nj
  ntotpix=float(ni)*float(nj)

  mmat3=0.d0

  fmask=0.
  do j=1,nj
     do i=1,ni
        fmask(2*i-1,j)=mask(i,j)/ntotpix
     enddo
  enddo
  deallocate(mask)
  call fourn(fmask,nn,2,-1)

  fmask2=0.
  do j=1,nj
     do i=1,ni
        fmask2(2*i-1,j)=mask2(i,j)/ntotpix
     enddo
  enddo
  deallocate(mask2)
  call fourn(fmask2,nn,2,-1)

  lside(1)=rad1
  lside(2)=rad2
  
  call jstartend2(size,rank,jstart,jend,lside)
  write(*,*) rank,jstart,jend
  
  num_part=0.
  mmat_part=0.d0
  
  if (rank==0) then
     write(*,*) 'start calculation of mode coupling matrix &
          due to the simulated mask'
  endif
  
  dj=jend-jstart+1
  do j=jstart,jend
     if (mod(j,16)==0) write(*,*) 'rank:',rank, int(100.*(j-jstart)/dj),'% finish'
     j1 = j-1
     if (j1>nn(2)/2)  j1 = j1 - nn(2)
     do i=1,nn(1)/2+1
        i1 = i-1
        
        lx=float(i1)/lside(1)*(2*pi)
        ly=float(j1)/lside(2)*(2*pi)
        g1=sqrt(float(i1*i1+j1*j1))
        labs=sqrt(lx*lx+ly*ly)
        if (labs<llowbin(1)) cycle
        if (labs>=llowbin(npmax+1)) cycle
        pbl=labs*(labs+1)/(2.*pi)
        
        np1=1
        do while (labs>=llowbin(np1))
           np1=np1+1
        enddo
        np1=np1-1
        
        if ((i1==0).or.(i1==nn(1)/2)) then
           ifac=1
        else
           ! symmetry between i and nn(1)+2-i
           ifac=2
        endif
        
        num_part(np1)=num_part(np1)+ifac
        
        do jd=1,nn(2)
           jd1 = jd-1
           if (jd1>nn(2)/2)  jd1 = jd1 - nn(2)
           do id=1,nn(1)
              id1 = id-1
              if (id1>nn(1)/2)  id1 = id1 - nn(1)
              
              l2x=float(id1)/lside(1)*2.*pi
              l2y=float(jd1)/lside(2)*2.*pi
              g2=sqrt(float(id1*id1+jd1*jd1))
              labs2=sqrt(l2x*l2x+l2y*l2y)
              if (labs2<llowbin(1)) cycle
              if (labs2>=llowbin(npmax+1)) cycle
              qbl=(2.*pi)/labs2/(labs2+1)
              
              np2=1
              do while (labs2>=llowbin(np2))
                 np2=np2+1
              enddo
              np2=np2-1
              
              im=i1-id1
              jm=j1-jd1
!
!!! comment out if aliasing effect is included
!
!              if (abs(im)>nn(1)/2) cycle
!              if (abs(jm)>nn(2)/2) cycle
              
              im=mod(im+nn(1),nn(1))+1
              jm=mod(jm+nn(2),nn(2))+1
              
              powm=(fmask(2*im-1,jm)*fmask2(2*im-1,jm)+fmask(2*im,jm)*fmask2(2*im,jm))
              
              cosj=(i1*id1+j1*jd1)/g1/g2
              cosj2=2.*cosj*cosj-1.
              sinj=(i1*jd1-j1*id1)/g1/g2
              sinj2=2.*sinj*cosj
              
              dmmat(1)=pbl*qbl*powm*cosj2*cosj2
              dmmat(2)=pbl*qbl*powm*sinj2*sinj2
              
              do k=1,ikmax
                 mmat_part(np1,np2,k)=mmat_part(np1,np2,k)+ifac*dmmat(k)
              enddo
           enddo
        enddo
     enddo
  enddo
  
  call MPI_REDUCE(num_part,num,npmax,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  do k=1,ikmax
     do np2=1,npmax
        call MPI_REDUCE(mmat_part(:,np2,k),mmat(:,np2,k),npmax,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     enddo
  enddo
  
  if (rank==0) then
     
     do np2=1,npmax
        do np1=1,npmax
           do k=1,ikmax
              if (num(np1)>0.) then
                 mmat(np1,np2,k)=mmat(np1,np2,k)/num(np1)
              else
!  k=2 is off-diagonal components
                 if (np2==np1.and.k/=2) then
                    mmat(np1,np2,k)=1.d0
                 else
                    mmat(np1,np2,k)=0.d0
                 endif
              endif
           enddo
        enddo
     enddo
     
     do np2=1,npmax
        do np1=1,npmax
           mmat3(np1,np2)=mmat(np1,np2,1)
           mmat3(np1+npmax,np2+npmax)=mmat(np1,np2,1)
           mmat3(np1+npmax,np2)=mmat(np1,np2,2)
           mmat3(np1,np2+npmax)=mmat(np1,np2,2)
           mmat3(np1+npmax2,np2+npmax2)=mmat(np1,np2,1)-mmat(np1,np2,2)
        enddo
     enddo
     
  endif

  if (rank==0) then
     
     cform='(2(i4,1x),(1pe15.8,1x))'

     open(2,file=outf,status='unknown')
     do np2=1,npmax3
        do np1=1,npmax3
           write(2,cform) np1,np2,mmat3(np1,np2)
        enddo
        write(2,*)
     enddo
     close(2)
     
  endif

  write(*,*) 'finish: rank',rank
  call MPI_FINALIZE(ierror)

contains

  subroutine readmask(inf,mask)
    
    implicit none

    real, pointer :: mask(:,:)
    character :: inf*150

    allocate(mask(ni,nj))
    open(1,file=inf,status='old',form='unformatted')
    read(1) mask
    close(1)

  end subroutine readmask

  subroutine jstartend2(size,rank,jstart,jend,lside)

    implicit none

    integer :: i,j,i1,j1
    integer, intent(in) :: size,rank
    integer, intent(out) :: jstart,jend
    integer :: nmode,neach
    double precision :: lx,ly,labs
    double precision, intent(in) :: lside(2)

    nmode=0
    do j=1,nj
       j1=j-1
       if (j1>nj/2)  j1 = j1 - nj
       do i=1,ni/2+1
          i1=i-1
          lx=i1*(2*pi/lside(1))
          ly=j1*(2*pi/lside(2))
          labs=sqrt(lx*lx+ly*ly)
          if (labs<llowbin(1)) cycle
          if (labs>=llowbin(npmax+1)) cycle
          nmode=nmode+1
       enddo
    enddo
    neach=nmode/size

    nmode=0
    do j=1,nj
       j1=j-1
       if (j1>nj/2)  j1 = j1 - nj
       do i=1,ni/2+1
          i1=i-1
          lx=i1*(2*pi/lside(1))
          ly=j1*(2*pi/lside(2))
          labs=sqrt(lx*lx+ly*ly)
          if (labs<llowbin(1)) cycle
          if (labs>=llowbin(npmax+1)) cycle
          nmode=nmode+1
       enddo
       if (nmode<=neach*rank) jstart=j+1
       if (nmode<=neach*(rank+1)) jend=j
    enddo

    if (rank==0) jstart=1
    if (rank==size-1) jend=nj

  end subroutine jstartend2

end program modecoupmat_inside

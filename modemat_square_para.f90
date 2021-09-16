program modecoupmat_square

  use fieldinfo_
  use sub_
  implicit none

  integer :: np,np2,i,j,i1,j1,im,dj
  integer, parameter :: immax=4
  real ::  lx,ly,labs,num_part(npmax),num(npmax),fac
  double precision :: dmat(npmax,immax),mat_part(npmax,npmax,immax),mat(npmax,npmax,immax)
  double precision :: mat3(npmax3,npmax3),lside(2)
  character ::  outf*150,cform*150

  include 'mpif.h'
  integer :: rank,size,ierror,status(MPI_STATUS_SIZE)
  integer :: jstart,jend

  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

  call getarg(1,outf)
  call binninginit

  if (rank==0) write(*,*) 'Start computing the convolution matrix due to a square patch of sky'

  mat3=0.d0

  lside(1)=rad1
  lside(2)=rad2
     
  call jstartend(size,rank,jstart,jend,lside)
  write(*,*) rank,jstart,jend
  
  mat_part=0.d0
  num_part=0.
  dj=jend-jstart+1
  do j=jstart,jend
     if (mod(j,4)==0) write(*,*) 'rank:',rank, int(100.*(j-jstart)/dj),'% finish'
     j1 = j-1
     do i=j,ni/2+1
        i1 = i-1
        
        lx=float(i1)/lside(1)*(2.*pi)
        ly=float(j1)/lside(2)*(2.*pi)
        labs=sqrt(lx*lx+ly*ly)
        
        if (labs<llowbin(1)) cycle
        if (labs>=llowbin(npmax+1)) cycle
        
        np=1
        do while (labs>=llowbin(np))
           np=np+1
        enddo
        np=np-1
        
        call matcomp(lx,ly,dmat,immax,lside)
        
        if ((j1==0).or.(i1==j1)) then
           fac=1.
        else
           fac=2.
        endif
        
        do np2=1,npmax
           do im=1,immax
              mat_part(np,np2,im)=mat_part(np,np2,im)+dmat(np2,im)*fac
           enddo
        enddo
        num_part(np)=num_part(np)+fac
        
     enddo
  enddo
  
  call MPI_REDUCE(num_part,num,npmax,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  
  mat=0.d0
  do im=1,immax
     do np2=1,npmax
        call MPI_REDUCE(mat_part(:,np2,im),mat(:,np2,im),npmax,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     enddo
  enddo
  
  if (rank==0) then
     
     do np=1,npmax
        do np2=1,npmax
           do im=1,immax
              if (num(np)>0) then
                 mat(np,np2,im)=mat(np,np2,im)/num(np)
              else
                 if (np==np2.and.im/=2) mat(np,np2,im)=1.d0
              endif
           enddo
        enddo
     enddo
     
     do np=1,npmax
        do np2=1,npmax
           mat3(np,np2)=mat(np,np2,1)
           mat3(np,np2+npmax)=mat(np,np2,2)
           mat3(np+npmax,np2)=mat(np,np2,2)
           mat3(np+npmax,np2+npmax)=mat(np,np2,1)
           mat3(np+npmax2,np2+npmax2)=mat(np,np2,1)-mat(np,np2,2)
        enddo
     enddo
     
  endif
  
  if (rank==0) then
     
     cform='(2(i4,1x),(1pe15.8,1x))'

     write(*,*) 'output file: ',outf(1:len_trim(outf))
     open(2,file=outf,status='unknown')
     do np=1,npmax3
        do np2=1,npmax3
           write(2,cform) np,np2,mat3(np,np2)
        enddo
        write(2,*)
     enddo
     close(2)
     
  endif

  write(*,*) 'finish calculation: rank',rank,jstart,jend
  call MPI_FINALIZE(ierror)

contains

  subroutine matcomp(qx,qy,dmat,immax,lside)

    implicit none
    
    integer ::  ik,is,ismax,np
    integer, intent(in) :: immax
    real ::  qx,qy,q,cosi,sini,ds,theta,kpara,kperp,labs,pbl,qbl,s2,c2
    real ::  kx,ky,winx,winy,win,rx,ry,gx,gy
    integer :: ikmax1,ikmax
    integer, parameter :: iklim=ni*100
    real ::  kf,knyq,kmax2,dk1,kmax1,dk2,k(iklim),dk(iklim)
    double precision ::  dmat(npmax,immax)
    double precision, intent(in) :: lside(2)
    
    q=sqrt(qx*qx+qy*qy)
    cosi=qx/q
    sini=qy/q
    
    pbl=q*(q+1)/2./pi

!!!
!  integral range of k is set from 0 to knyq/5.
!  From 0 to kf*4, binning width is kf/100. 
!  From kf*4 to knyq/5, binning width is kf/4. 

    kf=2.*pi/lside(1)
    knyq=kf*(ni/2)

    kmax2=knyq/5.
    dk1=kf/100.
    ikmax1=400
    kmax1=dk1*ikmax1
    dk2=kf/4.
    ikmax=ikmax1+nint((kmax2-kmax1)/dk2)

    do ik=1,ikmax1
       dk(ik)=dk1
       k(ik)=dk1*(ik-0.5)
    enddo

    do ik=ikmax1+1,ikmax
       dk(ik)=dk2
       k(ik)=kmax1+dk2*(ik-ikmax1-0.5)
    enddo

    dmat=0.d0
    
    do ik=1,ikmax
       ismax=2*max(ik,200)
       ds=2.*pi/ismax
       do is=1,ismax
          theta=ds*(is-1.)
          kpara=k(ik)*cos(theta)
          kperp=k(ik)*sin(theta)
          labs=sqrt(q*q+k(ik)*k(ik)+2.*kpara*q)
          

          qbl=2.*pi/labs/(labs+1)

          s2=(kperp/labs)**2
          c2=1.-s2

          kx=kpara*cosi-kperp*sini
          ky=kpara*sini+kperp*cosi
          
          rx=nint((qx+kx)/kf)*kf
          ry=nint((qy+ky)/kf)*kf

          labs=sqrt(rx*rx+ry*ry)
          if (labs<llowbin(1)) cycle
          if (labs>=llowbin(npmax+1)) cycle

          np=1
          do while (labs>=llowbin(np))
             np=np+1
          enddo
          np=np-1
          
          gx=kx/kf
          gy=ky/kf
          call recwin(gx,ni,winx)
          call recwin(gy,nj,winy)
          winx=winx/kf
          winy=winy/kf
          win=winx*winy
          dmat(np,1)=dmat(np,1)+k(ik)*dk(ik)*ds*(c2-s2)*(c2-s2)*win*qbl
          dmat(np,2)=dmat(np,2)+k(ik)*dk(ik)*ds*4*s2*c2*win*qbl
          dmat(np,3)=dmat(np,3)+k(ik)*dk(ik)*ds*(c2-s2)*win*qbl
          dmat(np,4)=dmat(np,4)+k(ik)*dk(ik)*ds*win*qbl
       enddo
    enddo

    dmat=dmat*pbl

  end subroutine matcomp

  subroutine jstartend(size,rank,jstart,jend,lside)

    implicit none

    integer :: i,j,i1,j1
    integer, intent(in) :: size,rank
    integer, intent(out) :: jstart,jend
    integer :: nmode,neach
    double precision :: lx,ly,labs
    double precision, intent(in) :: lside(2)

    nmode=0
    do j=1,nj/2+1
       j1=j-1
       do i=j,ni/2+1
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
    do j=1,nj/2+1
       j1=j-1
       do i=j,ni/2+1
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
    if (rank==size-1) jend=nj/2+1

  end subroutine jstartend
     
end program Modecoupmat_square

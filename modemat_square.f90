program modecoupmat_square

  use fieldinfo_
  use sub_
  implicit none

  integer :: np,np2,i,j,i1,j1
  real ::  lx,ly,labs,num(npmax),fac
  double precision :: dmat(npmax2),mat(npmax,npmax2)
  character ::  outf*100

  call getarg(1,outf)
  call binninginit

  mat=0.d0
  num=0.

  write(*,*) 'Start computing the convolution matrix due to a square patch of sky'

  do j=1,nj/2+1
     j1 = j-1
     write(*,*) j,'/',nj/2+1
     do i=j,ni/2+1
        i1 = i-1

        lx=float(i1)/rad1*(2.*pi)
        ly=float(j1)/rad2*(2.*pi)
        labs=sqrt(lx*lx+ly*ly)

        if (labs<llowbin(1)) cycle
        if (labs>=llowbin(npmax+1)) cycle

        np=1
        do while (labs>=llowbin(np))
           np=np+1
        enddo
        np=np-1

        call matcomp(lx,ly,dmat)

        if ((j1==0).or.(i1==j1)) then
           fac=1.
        else
           fac=2.
        endif

        do np2=1,npmax2
           mat(np,np2)=mat(np,np2)+dmat(np2)*fac
        enddo
        num(np)=num(np)+fac

     enddo
  enddo

  do np=1,npmax
     do np2=1,npmax2
        if (num(np)>0) then
           mat(np,np2)=mat(np,np2)/num(np)
        else
           if (np==np2) mat(np,np2)=1.d0
        endif
     enddo
  enddo

  call outmat(outf,mat)

contains

  subroutine matcomp(qx,qy,dmat)

    use fieldinfo_
    implicit none
    
    integer ::  ik,is,ismax,np
    real ::  qx,qy,q,cosi,sini,ds,theta,kpara,kperp,labs,pbl,qbl,s2,c2
    real ::  kx,ky,winx,winy,win,rx,ry,gx,gy
    integer :: ikmax1,ikmax
    integer, parameter :: iklim=ni*100
    real ::  kf,knyq,kmax,dk1,kmax1,dk2,k(iklim),dk(iklim)
    double precision ::  dmat(npmax2)
    
    q=sqrt(qx*qx+qy*qy)
    cosi=qx/q
    sini=qy/q

    pbl=q*q/2./pi

!  integral range of k is set from 0 to knyq/5.
!  From 0 to kf*4, binning width is kf/100. 
!  From kf*4 to knyq/5, binning width is kf/4. 

    kf=2.*pi/rad1
    knyq=kf*(ni/2)

    kmax=knyq/5.
    dk1=kf/100.
    ikmax1=400
    kmax1=dk1*ikmax1
    dk2=kf/4.
    ikmax=ikmax1+nint((kmax-kmax1)/dk2)

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
          
          qbl=2.*pi/labs/labs

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
          dmat(np)=dmat(np)+k(ik)*dk(ik)*ds*(c2-s2)*(c2-s2)*win*qbl
          dmat(np+npmax)=dmat(np+npmax)+k(ik)*dk(ik)*ds*4*s2*c2*win*qbl
       enddo
    enddo

    dmat=dmat*pbl

  end subroutine matcomp

  subroutine recwin(s,nd,win)
  
    implicit none
    integer :: nd
    real :: s,win,spi
    real, parameter :: pi=3.14159265358979
  
    if (abs(s)<1.e-10) then
       win=1.
    else
       spi=s*pi
       win=(sin(spi)/(spi))/(sin(spi/nd)/(spi/nd))
       win=win*win
    endif
    
  end subroutine recwin

  subroutine recwin2(s,nd,win)
    
    implicit none
    integer :: nd
    real :: s,win,spi
    real, parameter :: pi=3.14159265358979
    
    if (abs(s)<1.e-10) then
       win=1./4.
    else
       spi=s*pi
       win=(sin(spi/2.)/(spi))/(sin(spi/nd)/(spi/nd))
       win=win*win
    endif
    
  end subroutine recwin2

  subroutine welchwin(s,nd,win)
    
    implicit none
    integer :: nd
    real :: s,win,spi,spind
    real, parameter :: pi=3.14159265358979
    
    if (abs(s)<1.e-10) then
       win=4./9.
    else
       spi=s*pi
       spind=spi/nd
       win=(cos(spi)-(sin(spi)/sin(spind)/nd))/(sin(spind)*nd)**2
       win=4.*win*win
    endif
    
  end subroutine welchwin

  subroutine outmat(outf,mat)

    use fieldinfo_
    implicit none
    integer :: np,np2
    double precision :: mat(npmax,npmax2)    
    character :: outf*100

    write(*,*) 'output file: ',outf(1:len_trim(outf))
    
    open(2,file=outf,status='unknown')
    do np=1,npmax
       do np2=1,npmax2
          write(2,*) np,np2,mat(np,np2)
       enddo
       do np2=npmax2+1,npmax3
          write(2,*) np,np2,0.d0
       enddo
       write(2,*)
    enddo
    
    do np=npmax+1,npmax2
       do np2=1,npmax
          write(2,*) np,np2,mat(np-npmax,np2+npmax)
       enddo
       do np2=npmax+1,npmax2
          write(2,*) np,np2,mat(np-npmax,np2-npmax)
       enddo
       do np2=npmax2+1,npmax3
          write(2,*) np,np2,0.d0
       enddo
       write(2,*)
    enddo
    
    do np=npmax2+1,npmax3
       do np2=1,npmax2
          write(2,*) np,np2,0.d0
       enddo
       do np2=npmax2+1,npmax3
          write(2,*) np,np2,mat(np-npmax2,np2-npmax2)-mat(np-npmax,np2-npmax2)
       enddo
       write(2,*)
    enddo
    close(2)

  end subroutine outmat
  
end program Modecoupmat_square

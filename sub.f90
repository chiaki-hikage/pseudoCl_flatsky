module sub_

contains

  subroutine binninginit
    
    use fieldinfo_
    implicit none
    
    integer ::  np,i,j,i1,j1
    real ::  dlogl,lx,ly,labs
    real, dimension(npmax) :: lbinm,num
    
    if (bkind=='lin') then
       do np=1,npmax
          dlbin(np)=(lmax-lmin)/npmax
          lbinm(np)=lmin+dlbin(np)*(np-0.5)
          llowbin(np)=lmin+dlbin(np)*(np-1)
       enddo
    else
       do np=1,npmax
          dlogl=log(lmax/lmin)/npmax
          dlbin(np)=lmin*(exp(dlogl*np)-exp(dlogl*(np-1)))
          lbinm(np)=lmin*exp(dlogl*(np-0.5))
          llowbin(np)=lmin*exp(dlogl*(np-1))
       enddo
    endif
    llowbin(npmax+1)=lmax

    lbin=0.
    num=0.
    
    do j=1,nj
       j1 = j-1
       if (j1 > nj/2)  j1 = j1 - nj
       do i=1,ni
          i1 = i-1
          lx=float(i1)*(2.*pi)/rad1
          ly=float(j1)*(2.*pi)/rad2
          labs=sqrt(lx*lx+ly*ly)
          
          if (labs<llowbin(1)) cycle
          if (labs>=llowbin(npmax+1)) cycle
          
          np=1
          do while (labs>=llowbin(np))
             np=np+1
          enddo
          np=np-1
          if (i==1) then
             lbin(np)=lbin(np)+labs
             num(np)=num(np)+1.
          else
             lbin(np)=lbin(np)+2.*labs
             num(np)=num(np)+2.
          endif
       enddo
    enddo
    lbin=lbin/num
    
  end subroutine binninginit

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

  subroutine readcl1(inf,cl)
    
    use fieldinfo_
    implicit none
    integer :: np
    real :: x
    real, pointer :: cl(:)
    character, intent(in) :: inf*100

    if (inf=="none") then
       cl=0.
    else
       open(1,file=inf,status='old')
       do np=1,npmax
          read(1,*) x,cl(np)
       enddo
       close(1)
    endif

  end subroutine readcl1

  subroutine readcl3(inf,cl)
    
    use fieldinfo_
    implicit none
    integer :: np,np2,np3
    real :: x
    real, pointer :: cl(:)
    character, intent(in) :: inf*100

    if (inf=="none") then
       cl=0.
    else
       open(1,file=inf,status='old')
       do np=1,npmax
          np2=np+npmax
          np3=np+npmax2
          read(1,*) x,cl(np),cl(np2),cl(np3)
       enddo
       close(1)
    endif

  end subroutine readcl3

  subroutine outcl3(outf,cl)

    use fieldinfo_
    implicit none

    integer :: np,np2,np3
    real, pointer :: cl(:)
    character, intent(in) :: outf*100

    write(*,*) 'output cl:',outf(1:len_trim(outf))
    open(2,file=outf,status='unknown')
    do np=1,npmax
       np2=np+npmax
       np3=np+npmax2
       write(2,'(4(1pe12.5,1x))') lbin(np),cl(np),cl(np2),cl(np3)
    enddo
    close(2)
    
  end subroutine outcl3

  subroutine outcl1(outf,cl)

    use fieldinfo_
    implicit none

    integer :: np
    real, pointer :: cl(:)
    character, intent(in) :: outf*100

    write(*,*) 'output cl:',outf(1:len_trim(outf))
    open(2,file=outf,status='unknown')
    do np=1,npmax
       write(2,*) lbin(np),cl(np)
    enddo
    close(2)
    
  end subroutine outcl1

  subroutine outinitpow(outf,clkappa,clshear)

    use fieldinfo_
    implicit none

    integer :: np,np2,np3
    real, pointer :: clkappa(:),clshear(:)
    character, intent(in) :: outf*100

    open(2,file=outf,status='unknown')
    do np=1,npmax
       np2=np+npmax
       np3=np+npmax2
       write(2,"(7(1pe12.5,1x))") lbin(np),clkappa(np),0.,0., &
            clshear(np),clshear(np2),clshear(np3)
    enddo
    close(2)
    deallocate(clkappa,clshear)
    
  end subroutine outinitpow

  subroutine randist(pos,idum)

!  create random distribution in range of [0,ni] and [0,nj]
    
    use fieldinfo_
    implicit none

    integer :: n,idum
    real, pointer :: pos(:,:)

    allocate(pos(ntot,2))
    do n=1,ntot
       pos(n,1)=ran2(idum)*ni
       pos(n,2)=ran2(idum)*nj
    enddo

  end subroutine randist

  subroutine gridpos(pos)

    use fieldinfo_
    implicit none
    integer :: n,i,j
    real, pointer :: pos(:,:)

    allocate(pos(ntot,2))
    n=0
    do j=1,nj
       do i=1,ni
          n=n+1
          pos(n,1)=i-0.5
          pos(n,2)=j-0.5
       enddo
    enddo

  end subroutine gridpos

  subroutine outsheardat(outf,pos,sdata)

!  output shear data
!    1st line: position(1:ntot,1:2)
!    2nd line: shear(1:ntot)
    
    use fieldinfo_
    implicit none
  
    integer :: n
    real :: rpart,ipart
    real, pointer :: pos(:,:)
    complex, pointer :: sdata(:)
    character(100) :: outf
  
    write(*,*) 'output shear data:',outf
    open(2,file=outf,status='unknown',form='unformatted')
    write(2) pos
    write(2) sdata
    close(2)
    
  end subroutine outsheardat

  subroutine readmat(inf,mat)
    
    use fieldinfo_
    implicit none
    
    integer :: np1,np2,i,j
    double precision ::   mat(npmax3,npmax3)
    character :: inf*100
    
    mat=0.d0
    if (inf(1:4)=='none') then
       do np2=1,npmax3
          do np1=1,npmax3
             if (np1==np2) mat(np1,np2)=1.d0
          enddo
       enddo
    else
       write(*,*) 'read mode-mode coupling matrix from'
       write(*,*) inf
       open(1,file=inf,status='old')
       do np2=1,npmax3
          do np1=1,npmax3
             read(1,*) i,j,mat(i,j)
          enddo
          read(1,*)
       enddo
       close(1)
    endif
    
  end subroutine readmat

  subroutine calcpow(shear,pow)

!   shear(ni,nj): shear field (ni x nj)
!   pow(npmax3): power spectrum for EE, BB, and EB

    use fieldinfo_
    
    implicit none
    
    integer ::  i,j,i1,j1,np,np1,np2,nn(2)
    complex, pointer :: shear(:,:)
    real, allocatable ::   g1(:,:),g2(:,:)
    real ::  lx,ly,labs,area,cos1,sin1,cos2,sin2,cellarea,pbl
    real ::  num(npmax)
    real, pointer :: pow(:)
    real ::  g1r,g1i,g2r,g2i,elr,eli,blr,bli,powee,powbb,poweb
    character ::  outf*100
    
    write(*,*)  'Calculate Power Spectrum'
    nn(1)=ni
    nn(2)=nj
    area=rad1*rad2
    cellarea=area/float(ni)/float(nj)
    
    allocate(g1(ni*2,nj),g2(ni*2,nj))
    do j=1,nj
       do i=1,ni
          g1(2*i-1,j)=real(shear(i,j))*cellarea
          g1(2*i,j)=0.
          g2(2*i-1,j)=imag(shear(i,j))*cellarea
          g2(2*i,j)=0.
       enddo
    enddo
    
    num=0.
    allocate(pow(npmax3))
    pow=0.
    
    call fourn(g1,nn,2,-1)
    call fourn(g2,nn,2,-1)
    
    do j=1,nj
       j1 = j-1
       if (j1 > nj/2)  j1 = j1 - nj
       do i=1,ni/2+1
          i1 = i-1
          lx=float(i1)/rad1*(2.*pi)
          ly=float(j1)/rad2*(2.*pi)
          labs=sqrt(lx*lx+ly*ly)
          pbl=labs*labs/2./pi
          
          if (labs<llowbin(1)) cycle
          if (labs>=llowbin(npmax+1)) cycle
          
          cos1=lx/labs
          sin1=ly/labs
          cos2=cos1*cos1-sin1*sin1
          sin2=2.*sin1*cos1

          g1r = g1(2*i-1,j)
          g1i = g1(2*i,j)
          g2r = g2(2*i-1,j)
          g2i = g2(2*i,j)

          elr =  g1r*cos2 + g2r*sin2
          eli =  g1i*cos2 + g2i*sin2
          blr = -g1r*sin2 + g2r*cos2
          bli = -g1i*sin2 + g2i*cos2
          
          powee=elr*elr+eli*eli
          powbb=blr*blr+bli*bli
          poweb=elr*blr+eli*bli

          np=1
          do while (labs>=llowbin(np))
             np=np+1
          enddo
          np=np-1
          
          if (i==1) then
             num(np)=num(np)+1.
             pow(np)=pow(np)+powee*pbl
             pow(np+npmax)=pow(np+npmax)+powbb*pbl
             pow(np+npmax2)=pow(np+npmax2)+poweb*pbl
          else
             num(np)=num(np)+2.
             pow(np)=pow(np)+2.*powee*pbl
             pow(np+npmax)=pow(np+npmax)+2.*powbb*pbl
             pow(np+npmax2)=pow(np+npmax2)+2.*poweb*pbl
          endif
          
       enddo
    enddo
    
    do np=1,npmax
       if (num(np)>0.) then
          pow(np)=pow(np)/num(np)/area
          pow(np+npmax)=pow(np+npmax)/num(np)/area
          pow(np+npmax2)=pow(np+npmax2)/num(np)/area
       endif
    enddo
    
    deallocate(g1,g2)
    
  end subroutine calcpow
  
  subroutine calcpow_scalar(rhoin,pow)
    
    use fieldinfo_
    
    implicit none
    
    integer ::  i,j,i1,j1,np,np1,np2,nn(2)
    real, pointer ::  rhoin(:,:)
    real, allocatable :: rho(:,:)
    real ::  num(npmax),rhor,rhoi,lx,ly,labs,area,cellarea,pbl
    real ::  powrho,lnyq
    real, pointer ::  pow(:)
    character :: outf*100
    
    allocate(rho(2*ni,nj))
    nn(1)=ni
    nn(2)=nj
    area=rad1*rad2
    cellarea=area/float(ni)/float(nj)
    
    do j=1,nj
       do i=1,ni
          rho(2*i-1,j)=rhoin(i,j)*cellarea
          rho(2*i,j)=0.
       enddo
    enddo
    
    num=0.
    allocate(pow(npmax))    
    pow=0.
    
    call fourn(rho,nn,2,-1)
    
    lnyq=2.*pi/rad1*(float(ni)/2.)
    
    do j=1,nj
       j1 = j-1
       if (j1 > nj/2)  j1 = j1 - nj
       do i=1,ni/2+1
          i1 = i-1
          lx=float(i1)/rad1*(2.*pi)
          ly=float(j1)/rad2*(2.*pi)
          labs=sqrt(lx*lx+ly*ly)
          pbl=labs*labs/2./pi
          
          if (labs<llowbin(1)) cycle
          if (labs>=llowbin(npmax+1)) cycle
          
          rhor=rho(2*i-1,j)
          rhoi=rho(2*i,j)
          powrho=rhor*rhor+rhoi*rhoi
          
          np=1
          do while (labs>=llowbin(np))
             np=np+1
          enddo
          np=np-1
          
          if (i==1) then
             num(np)=num(np)+1.
             pow(np)=pow(np)+powrho*pbl
          else
             num(np)=num(np)+2.
             pow(np)=pow(np)+2.*powrho*pbl
          endif
          
       enddo
    enddo
    
    do np=1,npmax
       if (num(np)>0.) then
          pow(np)=pow(np)/num(np)/area
       endif
    enddo
    deallocate(rho)
    
  end subroutine calcpow_scalar

  subroutine readinitspec(inf,cl,ilmax)
    
    use fieldinfo_
    implicit none
    
    integer, intent(in) :: ilmax
    integer :: il,l
    real, intent(out) :: cl(ilmax)
    character, intent(in) :: inf*100
    
    write(*,*) 'Read Initial Spectrum from',inf
    
    cl=0.
    open(1,file=inf,status='old')
    do il=1,ilmax
       read(1,*) l,cl(l)
    enddo
    close(1)     
    
  end subroutine readinitspec
  
  subroutine deconv_kappa(cl,mat)

    use fieldinfo_
    implicit none
    
    integer :: np,np2
    real, pointer ::  cl(:)
    real, allocatable :: cl2(:)
    double precision,dimension(npmax3,npmax3),intent(in) :: mat
    double precision,dimension(npmax,npmax) :: arr,invarr

    allocate(cl2(npmax))
    do np=1,npmax
       cl2(np)=cl(np)
       cl(np)=0.
       do np2=1,npmax
          arr(np,np2)=mat(np,np2)+mat(np+npmax,np2)
       enddo
    enddo
    call calcinvarray(arr,invarr,npmax)
    cl=matmul(invarr,cl2)
    deallocate(cl2)

  end subroutine deconv_kappa

  subroutine conv_kappa(cl,mat)

    use fieldinfo_
    implicit none
    
    integer :: np,np2
    real, pointer ::  cl(:)
    real, allocatable :: cl2(:)
    double precision,dimension(npmax3,npmax3),intent(in) :: mat

    allocate(cl2(npmax))
    do np=1,npmax
       cl2(np)=cl(np)
       cl(np)=0.
    enddo
    do np=1,npmax
       do np2=1,npmax
          cl(np)=cl(np)+(mat(np,np2)+mat(np+npmax,np2))*cl2(np2)
       enddo
    enddo
    deallocate(cl2)

  end subroutine conv_kappa
 
  subroutine deconv_shear(cl,mat,mat2)
    
    use fieldinfo_
    implicit none
    
    real, pointer ::  cl(:)
    real, allocatable :: cl2(:)
    double precision,dimension(npmax3,npmax3),intent(in) :: mat,mat2
    double precision :: imat(npmax3,npmax3)
    
    allocate(cl2(npmax3))
    cl2=cl
    call calcinvarray(matmul(mat2,mat),imat,npmax3)
    cl=matmul(imat,cl2)
    call conv_kappa(cl,mat)
    deallocate(cl2)
    
  end subroutine deconv_shear

  subroutine conv_shear(cl,mat,mat2)
    
    use fieldinfo_
    implicit none
    
    real, pointer ::  cl(:)
    real, allocatable :: cl2(:)
    double precision,dimension(npmax3,npmax3),intent(in) :: mat,mat2
    
    allocate(cl2(npmax3))
    cl2=cl
    cl=matmul(mat2,matmul(mat,cl2))
    deallocate(cl2)
    
  end subroutine conv_shear

  subroutine addnoise(sdata,idum)

    use fieldinfo_
    implicit none

    integer :: idum
    integer :: n
    real :: rpart,ipart
    complex, pointer :: sdata(:)

    write(*,*) 'add noise'
    do n=1,ntot
       rpart=normsinv(dble(ran2(idum)))*intshear
       ipart=normsinv(dble(ran2(idum)))*intshear
       sdata(n)=sdata(n)+cmplx(rpart,ipart)
    enddo

  end subroutine addnoise

  subroutine rotateshear(sdata,idum)
    
    use fieldinfo_
    implicit none
    
    integer :: n,idum
    complex, pointer :: sdata(:)
    real :: phase
    
    write(*,*) 'rotate shape'
    do n=1,ntot
       phase=2*pi*ran2(idum)
       sdata(n)=sdata(n)*cmplx(cos(phase),sin(phase))
    enddo
    
  end subroutine rotateshear

  real function dist(a,b,n)
    
    integer :: i,n
    real :: a(n),b(n),d
    
    d=0.
    do i=1,n
       d=d+(a(i)-b(i))*(a(i)-b(i))
    enddo
    
    dist=sqrt(d)
    
  end function dist

  SUBROUTINE spline(x,y,n,yp1,ypn,y2)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n
    REAL, INTENT(in) :: yp1,ypn
    REAL, DIMENSION(n), INTENT(in) :: x,y
    REAL, DIMENSION(n), INTENT(out) :: y2
    INTEGER, parameter :: NMAX=1500
    INTEGER :: i,k
    REAL :: p,qn,sig,un
    REAL, DIMENSION(NMAX) :: u
    if (yp1>.99e30) then
       y2(1)=0.
       u(1)=0.
    else
       y2(1)=-0.5
       u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    do i=2,n-1
       sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
       p=sig*y2(i-1)+2.
       y2(i)=(sig-1.)/p
       u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1))) &
            /(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    if (ypn>.99e30) then
       qn=0.
       un=0.
    else
       qn=0.5
       un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif
    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
    do k=n-1,1,-1
       y2(k)=y2(k)*y2(k+1)+u(k)
    enddo
    
  END SUBROUTINE spline

  SUBROUTINE splint(xa,ya,y2a,n,x,y)

    IMPLICIT NONE
    INTEGER, INTENT(in) :: n
    REAL, INTENT(in) :: x
    REAL, INTENT(in), DIMENSION(n) :: xa, y2a, ya
    REAL, INTENT(out) :: y
    INTEGER :: k,khi,klo
    REAL :: a,b,h
    
    klo=1
    khi=n
    do
       k=(khi+klo)/2
       if(xa(k)>x)then
          khi=k
       else
          klo=k
       endif
       if (khi-klo<=1) exit
    end do
    h=xa(khi)-xa(klo)
    if (h==0.) then
       write(*,*) 'bad xa input in splint'
       stop
    endif
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
    
  END SUBROUTINE splint


  SUBROUTINE fourn(data,nn,ndim,isign)

!
!   Replaces DATA by its NDIM-dimensional Fourier transform,
!  if ISIGN is input as 1. NN is an integer array of length NDIM,
!  containing the lengths of each dimension (number of complex
!  values), which MUST all be powers of 2. DATA is a real array
!  of length twice the product of these lengths, in which the data
!  are stored as in a multidimensional complex FORTRAN array. If
!  ISIGN is input as -1, DATA is replaced by its inverse transform
!  times the product of the lengths of all dimensions.
!
!     adopted from Numerical Recipes chapter 12, p.451.
! ====================================================================
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ndim, isign, nn(ndim)
    REAL, INTENT(INOUT) :: data(2*product(nn))
    DOUBLE PRECISION, PARAMETER :: pi=3.141592653589793238d0
    DOUBLE PRECISION :: wr, wi, wpr, wpi, wtemp, theta
    INTEGER :: ntot, idim, nprev, n, nrem, ibit, ifp1, ifp2
    INTEGER :: i1, i2, i3, ip1, ip2, ip3, i2rev, i3rev, k1, k2
    REAL :: tempr, tempi
!
!  Compute total number of complex values.
!  ---------------------------------------
    ntot = product(nn)
!
!  Main loop over the dimensions.
!  ------------------------------
    nprev = 1
    do idim=1,ndim
       n = nn(idim)
       nrem  = ntot/(n*nprev)
       ip1   = 2*nprev
       ip2   = ip1*n
       ip3   = ip2*nrem
       i2rev = 1
!     
!  This is the bit reversal section of the routine.
!  ------------------------------------------------
       do i2 = 1,ip2,ip1
          if (i2<i2rev) then
             do i1 = i2, i2+ip1-2,2
                do i3 = i1, ip3, ip2
                   i3rev = i2rev + i3 - i2
                   tempr = data(i3)
                   tempi = data(i3+1)
                   data(i3) = data(i3rev)
                   data(i3+1) = data(i3rev+1)
                   data(i3rev) = tempr
                   data(i3rev+1) = tempi
                enddo
             enddo
          endif
          ibit = ip2/2
          do while (ibit>ip1.and.i2rev>ibit)
             i2rev = i2rev - ibit
             ibit  = ibit/2
          end do
          i2rev = i2rev + ibit
       end do
!
!  Here begins the Danielson-Lanczos section of the routine.
!  ---------------------------------------------------------
       ifp1 = ip1
       do while (ifp1<ip2)
          ifp2  = 2*ifp1
          theta = isign*2.0*pi/(ifp2/ip1)
          wpr   = -2.d0*dsin(0.5d0*theta)**2
          wpi   = dsin(theta)
          wr    = 1.d0
          wi    = 0.d0
          do i3=1,ifp1,ip1
             do i1=i3,i3+ip1-2,2
                do i2=i1,ip3,ifp2
                   k1    = i2
                   k2    = k1 + ifp1
                   tempr = sngl(wr)*data(k2)   - sngl(wi)*data(k2+1)
                   tempi = sngl(wr)*data(k2+1) + sngl(wi)*data(k2)
                   data(k2)   = data(k1) - tempr
                   data(k2+1) = data(k1+1) - tempi
                   data(k1)   = data(k1) + tempr
                   data(k1+1) = data(k1+1) + tempi
                enddo
             enddo
             wtemp = wr
             wr    = wr*wpr - wi*wpi + wr
             wi    = wi*wpr + wtemp*wpi + wi
          enddo
          ifp1 = ifp2
       end do
       nprev = n*nprev
    enddo
    
  END SUBROUTINE fourn


  subroutine ave_sum(y,nmax,mean,var,nave)
    
    implicit none
    integer, intent(in) :: nmax,nave
    real, dimension(nmax), intent(in) :: y
    real, dimension(nmax), intent(inout) :: mean, var
    
    if (nave==1) then
       mean=y
       var=y*y
    else
       mean=mean+y
       var=var+y*y
    endif
    
  end subroutine ave_sum


  subroutine ave_div(mean,err,nmax,nave)
    
    implicit none
    
    integer, intent(in) :: nmax, nave
    real, dimension(nmax), intent(inout) ::  mean, err
    
    if (nave==1) then
       mean=mean
       err=0.
    else
       mean=mean/float(nave)
       err=err-float(nave)*mean**2
       err=err/float(nave-1)
       err=sqrt(err)
    endif
    
  end subroutine ave_div


  subroutine calcinvarray(covorg,incov,nmax)
    
    implicit none
    
    integer, intent(in) ::  nmax
    double precision, intent(in) :: covorg(nmax,nmax)
    double precision, intent(out) :: incov(nmax,nmax)
    double precision :: cov(nmax,nmax)
    integer ::  k, nx, ny, indx(nmax)
    double precision :: d, col(nmax)
    
    cov=covorg
    call ludcmp(cov,nmax,nmax,indx,d)
    do ny=1,nmax
       col=0.d0
       col(ny)=1.d0
       call lubksb(cov,nmax,nmax,indx,col)
       do nx=1,nmax
          incov(nx,ny)=col(nx)
       enddo
    enddo

    call invarraycheck(covorg,incov,nmax,nmax)
    
  end subroutine calcinvarray
  

  SUBROUTINE ludcmp(a,n,np,indx,d)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, np
    INTEGER, INTENT(OUT) :: indx(n)
    INTEGER, PARAMETER :: NMAX=1300
    DOUBLE PRECISION, INTENT(OUT) :: d
    DOUBLE PRECISION, INTENT(INOUT) :: a(np,np)
    DOUBLE PRECISION, PARAMETER :: TINY=1.0d-20
    INTEGER :: i,imax,j,k
    DOUBLE PRECISION :: aamax,dum,sum,vv(NMAX)
    
    d=1.
    do i=1,n
       aamax=0.d0
       do j=1,n
          if (abs(a(i,j))>aamax) aamax=abs(a(i,j))
       enddo
       if (aamax==0.) then
          write(*,*) 'singular matrix in ludcmp'
          stop
       endif
       vv(i)=1./aamax
    enddo
    do j=1,n
       do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
             sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
       enddo
       aamax=0.
       do i=j,n
          sum=a(i,j)
          do k=1,j-1
             sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
             imax=i
             aamax=dum
          endif
       enddo
       if (j/=imax)then
          do k=1,n
             dum=a(imax,k)
             a(imax,k)=a(j,k)
             a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
       endif
       indx(j)=imax
       if (a(j,j)==0.) a(j,j)=TINY
       if (j/=n) then
          dum=1./a(j,j)
          do i=j+1,n
             a(i,j)=a(i,j)*dum
          enddo
       endif
    enddo
    
  END SUBROUTINE ludcmp


  SUBROUTINE lubksb(a,n,np,indx,b)
    INTEGER, INTENT(IN) :: n,np
    INTEGER, INTENT(IN) :: indx(n)
    DOUBLE PRECISION, INTENT(IN) :: a(np,np)
    DOUBLE PRECISION, INTENT(INOUT) :: b(n)
    INTEGER :: i,ii,j,ll
    DOUBLE PRECISION :: sum
    
    ii=0
    do i=1,n
       ll=indx(i)
       sum=b(ll)
       b(ll)=b(i)
       if (ii/=0) then
          do j=ii,i-1
             sum=sum-a(i,j)*b(j)
          enddo
       else if (sum/=0.) then
          ii=i
       endif
       b(i)=sum
    enddo
    do i=n,1,-1
       sum=b(i)
       do j=i+1,n
          sum=sum-a(i,j)*b(j)
       enddo
       b(i)=sum/a(i,i)
    enddo
    
  END SUBROUTINE lubksb


  subroutine invarraycheck(a,b,nmax,ncheck)
    
    implicit none
    
    integer, intent(in) :: nmax,ncheck
    double precision, dimension(nmax,nmax), intent(in) :: a, b
    integer ::  nx,ny
    double precision :: e
    double precision, parameter :: tolerance=1.d-8
    
    do ny=1,ncheck
       do nx=1,ncheck
          e=dot_product(a(nx,:),b(:,ny))
          if (nx==ny) then
             if (abs(e-1.d0)>tolerance) then
                write(*,*) nx,ny,e
                write(*,*) 'accuracy of invarray is bad'
             endif
          else
             if (abs(e)>tolerance) then
                write(*,*) nx,ny,e
                write(*,*) 'accuracy of invarray is bad'
             endif
          endif
       enddo
    enddo
    
  end subroutine invarraycheck
  
  double precision function normsinv(p)

    implicit none
    
    double precision, parameter ::  &
         A1=-3.969683028665376e+01, &
         A2=2.209460984245205e+02 , &
         A3=-2.759285104469687e+02, &
         A4=1.383577518672690e+02, &
         A5=-3.066479806614716e+01, &
         A6=2.506628277459239e+00, &
         B1=-5.447609879822406e+01, &
         B2= 1.615858368580409e+02, &
         B3= -1.556989798598866e+02, &
         B4=   6.680131188771972e+01, &
         B5= -1.328068155288572e+01, &
         C1= -7.784894002430293e-03, &
         C2= -3.223964580411365e-01, &
         C3= -2.400758277161838e+00, &
         C4= -2.549732539343734e+00, &
         C5=  4.374664141464968e+00, &
         C6=  2.938163982698783e+00, &
         D1=  7.784695709041462e-03, &
         D2=  3.224671290700398e-01, &
         D3=  2.445134137142996e+00, &
         D4=  3.754408661907416e+00, &
         P_LOW= 0.02425, &
         P_HIGH=  0.97575
    
    double precision, intent(in) :: p
    double precision :: x, q, r, u, e
    if ((0<p).and.(p<P_LOW)) then
       q = sqrt(-2*log(p))
       x = (((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6)/((((D1*q+D2)*q+D3)*q+D4)*q+1)
    elseif ((P_LOW<=p).and.(p<=P_HIGH)) then
       q = p - 0.5
       r = q*q
       x = (((((A1*r+A2)*r+A3)*r+A4)*r+A5)*r+A6)*q &
            /(((((B1*r+B2)*r+B3)*r+B4)*r+B5)*r+1)
    elseif ((P_HIGH<p).and.(p<1)) then
       q = sqrt(-2*log(1-p))
       x = -(((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) &
            / ((((D1*q+D2)*q+D3)*q+D4)*q+1)
    endif
    
    normsinv=x

  end function normsinv

  REAL FUNCTION ran2(idum)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: idum
    INTEGER, PARAMETER :: IM1=2147483563,IM2=2147483399,IMM1=IM1-1
    INTEGER, PARAMETER :: IA1=40014,IA2=40692,IQ1=53668,IQ2=52774
    INTEGER, PARAMETER :: IR1=12211,IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB
    REAL, PARAMETER :: AM=1./IM1,EPS=1.2e-7,RNMX=1.-EPS
    INTEGER :: idum2=123456789,j,k,iy=0
    INTEGER, DIMENSION(NTAB) :: iv=0
    
    if (idum.le.0) then
       idum=max(-idum,1)
       idum2=idum
       do j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
       enddo
       iy=iv(1)
    endif
    k=idum/IQ1
    idum=IA1*(idum-k*IQ1)-k*IR1
    if (idum.lt.0) idum=idum+IM1
    k=idum2/IQ2
    idum2=IA2*(idum2-k*IQ2)-k*IR2
    if (idum2.lt.0) idum2=idum2+IM2
    j=1+iy/NDIV
    iy=iv(j)-idum2
    iv(j)=idum
    if(iy.lt.1)iy=iy+IMM1
    ran2=min(AM*iy,RNMX)
    return
  END FUNCTION ran2
  
end module sub_

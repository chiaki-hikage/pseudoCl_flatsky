program modecoupmat_inside

!  compute the convolution matrix due to a mask inside a square patch
!   inf: input file of mask field [1:ni,1:nj] (real)
!   outf: output field of mode coupling matrix of EE, BB, EB power
!          [1:npmax*3,1:npmax*3] (real)

  use fieldinfo_
  use sub_
  implicit none

  integer :: i,j,i1,j1,id,jd,id1,jd1,im,jm,nn(2),np1,np2
  real, pointer :: mask(:,:)
  real ::  num(npmax),fmask(ni*2,nj)
  real ::  l,area,cellarea,pbl,qbl,powm,lnyq
  real ::  l1,l2,l3,g1,g2,g3,cosj,cosj2,sinj,sinj2,l2x,l2y
  double precision ::  mmat(npmax2,npmax2),dmmat_auto,dmmat_cros
  character :: inf*100,outf*100

  call getarg(1,inf)
  call getarg(2,outf)

  call binninginit
  call readmask(inf,mask)

  nn(1)=ni
  nn(2)=nj
  area=rad1*rad2
  cellarea=area/float(ni)/float(nj)

  fmask=0.
  do j=1,nj
     do i=1,ni
        fmask(2*i-1,j)=mask(i,j)*cellarea
     enddo
  enddo

  deallocate(mask)

  call fourn(fmask,nn,2,-1)

  num=0.
  mmat=0.d0

  lnyq=2.*pi/rad1*(float(ni)/2.)

  write(*,*) 'start calculation of mode coupling matrix due to the simulated mask'
  do j=1,nj
     write(*,*) j,'/',nj
     j1 = j-1
     if (j1>nj/2)  j1 = j1 - nj
     do i=1,ni/2+1
        i1 = i-1

        g1=sqrt(float(i1*i1+j1*j1))
        l1=g1/rad1*(2.*pi)
        if (l1<llowbin(1)) cycle
        if (l1>=llowbin(npmax+1)) cycle
        pbl=l1*l1/(2.*pi)

        np1=1
        do while (l1>=llowbin(np1))
           np1=np1+1
        enddo
        np1=np1-1

        if ((i1==0).or.(i1==ni/2)) then
           num(np1)=num(np1)+1.
        else
           num(np1)=num(np1)+2.
        endif
        
        do jd=1,nj
           jd1 = jd-1
           if (jd1>nj/2)  jd1 = jd1 - nj
           do id=1,ni
              id1 = id-1
              if (id1>ni/2)  id1 = id1 - ni

              l2x=float(id1)/rad1*2.*pi
              l2y=float(jd1)/rad2*2.*pi
              g2=sqrt(float(id1*id1+jd1*jd1))
              l2=g2/rad1*(2.*pi)
              if (l2<llowbin(1)) cycle
              if (l2>=llowbin(npmax+1)) cycle
              qbl=(2.*pi)/l2/l2

              np2=1
              do while (l2>=llowbin(np2))
                 np2=np2+1
              enddo
              np2=np2-1

              im=i1-id1
              jm=j1-jd1
!
!!! comment out if aliasing effect is included
!
!              if (abs(im)>ni/2) cycle
!              if (abs(jm)>nj/2) cycle

              im=mod(im+ni,ni)+1
              jm=mod(jm+nj,nj)+1

              powm=(fmask(2*im-1,jm)**2+fmask(2*im,jm)**2)/area

              cosj=(i1*id1+j1*jd1)/g1/g2
              cosj2=2.*cosj*cosj-1.
              sinj=(i1*jd1-j1*id1)/g1/g2
              sinj2=2.*sinj*cosj

              dmmat_auto=pbl*qbl*powm*cosj2*cosj2
              dmmat_cros=pbl*qbl*powm*sinj2*sinj2

              if ((i1>0).and.(i1<ni/2)) then
                 dmmat_auto=dmmat_auto*2
                 dmmat_cros=dmmat_cros*2
              endif

              mmat(np1,np2)=mmat(np1,np2)+dmmat_auto
              mmat(np1+npmax,np2+npmax)=mmat(np1+npmax,np2+npmax)+dmmat_auto
              mmat(np1,np2+npmax)=mmat(np1,np2+npmax)+dmmat_cros
              mmat(np1+npmax,np2)=mmat(np1+npmax,np2)+dmmat_cros
           enddo
        enddo
     enddo
  enddo

  do np2=1,npmax2
     do np1=1,npmax
        if (num(np1)>0.) then
           mmat(np1,np2)=mmat(np1,np2)/num(np1)/area
           mmat(np1+npmax,np2)=mmat(np1+npmax,np2)/num(np1)/area
        else
           if (np2==np1) mmat(np1,np2)=1.d0
           if (np2/=np1) mmat(np1,np2)=0.d0
        endif
     enddo
  enddo

  open(2,file=outf,status='unknown')
  do np2=1,npmax3
     do np1=1,npmax3
        if ((np1.le.npmax2).and.(np2.le.npmax2)) then
           write(2,*) np1,np2,mmat(np1,np2)
        elseif ((np1>npmax2).and.(np2>npmax2)) then
           i=mod(np1-1,npmax)+1
           j=mod(np2-1,npmax)+1
           write(2,*) np1,np2,mmat(i,j)-mmat(i+npmax,j)
        else
           write(2,*) np1,np2,0.d0
        endif
     enddo
     write(2,*)
  enddo
  close(2)

contains

  subroutine readmask(inf,mask)
    
    use fieldinfo_
    use sub_
    implicit none

    real, pointer :: mask(:,:)
    character :: inf*100

    allocate(mask(ni,nj))
    open(1,file=inf,status='old',form='unformatted')
    read(1) mask
    close(1)

  end subroutine readmask

end program modecoupmat_inside

program Gauss_shear_distribution

  use fieldinfo_
  use sub_
  implicit none

  integer :: idum,n,irun
  real, pointer ::  clkappa(:),clshear(:)
  real, pointer :: pos(:,:),kappa(:,:)
  complex, pointer :: sdata(:),shear(:,:)
  character :: outfdat*100,outfpow*100,crun*10,pcond*1

  call getarg(1,outfdat)
  call getarg(2,outfpow)
  call getarg(3,crun)
  read(crun,*) irun
  call getarg(4,pcond)

  call binninginit
  call makegashear(irun,shear,kappa,pcond)
  call calcpow_scalar(kappa,clkappa)
  deallocate(kappa)
  call calcpow(shear,clshear)
  call outinitpow(outfpow,clkappa,clshear)
  idum=-1
  call randist(pos,idum)
  call shearassign(shear,pos,sdata)
  deallocate(shear)
  call outsheardat(outfdat,pos,sdata)

contains

  subroutine makegashear(nreal,shear,kappa,pcond)
  
    use fieldinfo_
    use sub_
    
    implicit none

    integer, intent(in) :: nreal
    integer ::  i,j,ic,jc,i1,j1,np,idum,nn(2),fac
    real, pointer ::  g1(:,:),g2(:,:),ka(:,:),kappa(:,:)
    complex, pointer ::  shear(:,:)
    real :: clk,elr,eli,blr,bli,labs,lx,ly,area,c1,s1,c2,s2,rad1f,rad2f
    character, intent(in) :: pcond*1

    integer, parameter :: nin=69
    real :: x(nin),y(nin),y2(nin)
    real, parameter :: d2nl=1.e30
    character(100) :: infcl='clin.dat'

    idum=-nreal-1 
    if (pcond=='y') then
       fac=1
       write(*,*) 'Make a Gaussian Shear Field'
    else
       fac=2
       write(*,*) 'To break periodic boundary condition, &
            make a larger size Gaussian Shear Field and &
            use a part of it'
    endif
    write(*,*) 'initial seed: ',idum

    nn(1)=ni*fac
    nn(2)=nj*fac
    rad1f=rad1*fac
    rad2f=rad2*fac
    area=rad1f*rad2f
    
    allocate(g1(nn(1)*2,nn(2)),g2(nn(1)*2,nn(2)),ka(nn(1)*2,nn(2)))
    
    ka=0.
    g1=0.
    g2=0.
    
    open(1,file=infcl,status='old')
    do n=1,nin
       read(1,*) x(n),y(n)
       y(n)=y(n)/x(n)/x(n)*(2.*pi)
    enddo
    close(1)
    call spline(x,y,nin,d2nl,d2nl,y2)
    
    do j=1,nn(2)
       j1 = j-1
       if (j1>nn(2)/2)  j1 = j1 - nn(2)
       do i=1,nn(1)/2+1
          i1 = i-1
          if ((i==1).and.(j==1)) cycle
          if ((i==1).and.(j>nn(2)/2)) cycle
          lx=float(i1)/(rad1f)*(2.*pi)
          ly=float(j1)/(rad2f)*(2.*pi)
          labs=sqrt(lx*lx+ly*ly)
          c1=lx/labs
          s1=ly/labs
          c2=c1*c1-s1*s1
          s2=2*c1*s1
          
          call splint(x,y,y2,nin,labs,clk)
          
          elr=sqrt(clk/2.)*normsinv(dble(ran2(idum)))*sqrt(area)
          eli=sqrt(clk/2.)*normsinv(dble(ran2(idum)))*sqrt(area)
          
          blr=0.
          bli=0.
          
          ka(2*i-1,j) = elr
          ka(2*i,j) = eli
          g1(2*i-1,j) = c2*elr - s2*blr
          g1(2*i,j)   = c2*eli - s2*bli
          g2(2*i-1,j) = s2*elr + c2*blr
          g2(2*i,j)   = s2*eli + c2*bli        
          
          if (((i1==0).or.(i1==nn(1)/2)).and.((j1==0).or.(j1==nn(2)/2))) then
             ka(2*i,j) = 0.
             g1(2*i,j) = 0.
             g2(2*i,j) = 0.
          endif
          
          ic = mod(nn(1)-i1,nn(1))+1
          jc = mod(nn(2)-j1,nn(2))+1
          ka(2*ic-1,jc) = ka(2*i-1,j)
          ka(2*ic,jc) = -ka(2*i,j)
          g1(2*ic-1,jc) =  g1(2*i-1,j)
          g1(2*ic,jc)   = -g1(2*i,j)
          g2(2*ic-1,jc) =  g2(2*i-1,j)
          g2(2*ic,jc)   = -g2(2*i,j)      
       enddo
    enddo
    
    call fourn(ka,nn,2,1)
    call fourn(g1,nn,2,1)
    call fourn(g2,nn,2,1)
    
    allocate(kappa(ni,nj))
    
    do j=1,nj
       do i=1,ni
          kappa(i,j)=ka(2*i-1,j)/area
       enddo
    enddo
    deallocate(ka)
    
    allocate(shear(ni,nj))
    do j=1,nj
       do i=1,ni
          shear(i,j)=cmplx(g1(2*i-1,j),g2(2*i-1,j))/area
       enddo
    enddo
    deallocate(g1,g2)
    
  end subroutine makegashear

  subroutine shearassign(shear,pos,sdata)

    use fieldinfo_

    implicit none

    integer :: i,j,n
    real, pointer :: pos(:,:)
    complex, pointer :: shear(:,:),sdata(:)

    allocate(sdata(ntot))
    do n=1,ntot
       i=int(pos(n,1))+1
       j=int(pos(n,2))+1
       sdata(n)=shear(i,j)
    enddo

  end subroutine shearassign
  
end program Gauss_shear_distribution

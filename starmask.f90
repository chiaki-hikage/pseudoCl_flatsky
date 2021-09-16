program make_starmask

  use fieldinfo_
  use sub_

  implicit none

!     idum: set arbitrary value for initial seed for ran2

  integer :: idum=-152103

!     inf: input file of shear distribution
!    outf: mask field [1:ni,1:nj](real)
!    outf2: flag data describing which data is masked or not [1:ntot] (logical*1)
!           if (flag(n)=.true.) n-th data is masked

  character(100) :: inf,outf,outf2

  integer :: i,j
  real, pointer :: mask(:,:),pos(:,:)
  logical*1, pointer :: flag(:)

  call getarg(1,inf)
  call getarg(2,outf)
  call getarg(3,outf2)

  call readdata(inf,pos)
  call maskdata(pos,flag,idum)
  call densassign(pos,flag,mask)
  call outmask(mask,outf)
  call outflag(outf2,flag)

contains
  
  subroutine readdata(inf,pos)

    use fieldinfo_
    implicit none

    integer :: n
    real, pointer :: pos(:,:)
    character :: inf*100

    write(*,*) 'read data from ',inf(1:len_trim(inf))
    allocate(pos(ntot,2))
    open(1,file=inf,status='old',form='unformatted')
    read(1) pos
    close(1)

  end subroutine readdata

  subroutine outflag(outf,flag)

    logical*1, pointer :: flag(:)
    character :: outf*100

    write(*,*) 'output flag:',outf(1:len_trim(outf))
    open(2,file=outf,status='unknown',form='unformatted')
    write(2) flag
    close(2)
    deallocate(flag)

  end subroutine outflag


  subroutine densassign(pos,flag,mask)

    use fieldinfo_
    use sub_
    implicit none

    integer :: i,j,i1,j1,idum,n
    real :: x,y,fx,fy
    real, pointer :: mask(:,:),pos(:,:)
    logical*1, pointer :: flag(:)

    allocate(mask(ni,nj))
    mask=0.
    do n=1,ntot
       if (flag(n)) cycle
       i=int(pos(n,1))+1
       j=int(pos(n,2))+1
       mask(i,j)=mask(i,j)+1.
    enddo
    close(1)
    deallocate(pos)

  end subroutine densassign

  subroutine maskdata(pos,flag,idum)

    use fieldinfo_
    use maskinfo_
    use sub_
    implicit none

    integer :: i,j,n,nmask,nm,next
    real :: hpos(2),hsize,ssize(2)  ! in unit of grid
    integer, intent(inout) :: idum
    real, pointer :: pos(:,:)
    integer, pointer :: num(:,:),pnum(:,:,:)
    logical*1, pointer :: flag(:)

    allocate(flag(ntot))
    flag=.false.
    call fieldassign(pos,num,pnum)
    next=0  
    nmask=nint(ntot*fmask)
    write(*,*) nmask
    do while (next<nmask)
       hpos(1)=ran2(idum)*ni
       hpos(2)=ran2(idum)*nj
       hsize=(hmin+(hmax-hmin)*ran2(idum)**4)*amin2grid
       call starmask(hsize,hpos,pos,num,pnum,flag,next)
       if (hsize>rsmin*amin2grid) then
          ssize(1)=hsize*xrec
          ssize(2)=hsize*yrec
          call recmask(ssize,hpos,pos,num,pnum,flag,next)
       endif
    enddo
    call zpad(pos,num,pnum,flag,next,idum)
    deallocate(num,pnum)

    write(*,*) 'unmasked fraction:',1-float(next)/float(ntot)
    
  end subroutine maskdata

  subroutine outmask(mask,outf)

    implicit none
    real, pointer :: mask(:,:)
    character(100) :: outf

    write(*,*) 'output mask:',outf(1:len_trim(outf))
    open(2,file=outf,status='unknown',form='unformatted')
    write(2) mask
    close(2)
    deallocate(mask)
    
  end subroutine outmask

  subroutine fieldassign(pos,num,pnum)

    implicit none
    integer :: i,j,n,nmax
    real, pointer :: pos(:,:)
    integer, pointer :: num(:,:),pnum(:,:,:)
    
    allocate(num(ni,nj))
    num=0
    do n=1,ntot
       i=int(pos(n,1))+1
       j=int(pos(n,2))+1
       num(i,j)=num(i,j)+1
    enddo
    nmax=0
    do j=1,nj
       do i=1,ni
          if (num(i,j)>nmax) nmax=num(i,j)
       enddo
    enddo
    write(*,*) 'maximum number of galaxies per pixel',nmax
    allocate(pnum(ni,nj,nmax))
    num=0
    do n=1,ntot
       i=int(pos(n,1))+1
       j=int(pos(n,2))+1
       num(i,j)=num(i,j)+1
       pnum(i,j,num(i,j))=n
    enddo

  end subroutine fieldassign

  subroutine starmask(hsize,hpos,pos,num,pnum,flag,next)
    
    use fieldinfo_
    use sub_
    implicit none

    integer :: i,j,imin,imax,jmin,jmax,n
    integer, intent(inout) :: next
    real :: r
    real, intent(in) :: hsize,hpos(2)
    real, pointer :: pos(:,:)
    integer, pointer :: num(:,:),pnum(:,:,:)
    logical*1, pointer :: flag(:)

    imin=int(hpos(1)-hsize)
    imax=int(hpos(1)+hsize)+1
    jmin=int(hpos(2)-hsize)
    jmax=int(hpos(2)+hsize)+1
    do i=imin,imax
       do j=jmin,jmax
          if ((i<1).or.(i>ni)) cycle
          if ((j<1).or.(j>nj)) cycle
          do n=1,num(i,j)
             if (flag(pnum(i,j,n))) cycle
             r=dist(hpos,pos(pnum(i,j,n),:),2)
             if (r<hsize) then
                next=next+1
                flag(pnum(i,j,n))=.true.
             endif
          enddo
       enddo
    enddo

  end subroutine starmask

  subroutine recmask(ssize,hpos,pos,num,pnum,flag,next)

    use fieldinfo_
    implicit none

    integer :: i,j,imin,imax,jmin,jmax,n
    integer, intent(inout) :: next
    real, intent(in) :: ssize(2),hpos(2)
    real :: dx,dy
    real, pointer :: pos(:,:)
    integer, pointer :: num(:,:),pnum(:,:,:)
    logical*1, pointer :: flag(:)

    imin=int(hpos(1)-ssize(1))
    imax=int(hpos(1)+ssize(1))+1
    jmin=int(hpos(2)-ssize(2))
    jmax=int(hpos(2)+ssize(2))+1
    do i=imin,imax
       do j=jmin,jmax
          if ((i<1).or.(i>ni)) cycle
          if ((j<1).or.(j>nj)) cycle
          do n=1,num(i,j)
             if (flag(pnum(i,j,n))) cycle
             dx=abs(pos(pnum(i,j,n),1)-hpos(1))
             dy=abs(pos(pnum(i,j,n),2)-hpos(2))
             if ((dx<ssize(1)).and.(dy<ssize(2))) then
                next=next+1
                flag(pnum(i,j,n))=.true.
             endif
          enddo
       enddo
    enddo
    
  end subroutine recmask

  subroutine zpad(pos,num,pnum,flag,next,idum)

    use fieldinfo_
    use maskinfo_
    use sub_
    implicit none

    integer :: i,j,n,i2,i2max
    integer, intent(inout) :: idum,next
    integer, parameter :: nisim=2048,njsim=nisim
    real :: x
    real, pointer :: pos(:,:)
    integer, pointer :: num(:,:),pnum(:,:,:)
    logical*1, pointer :: flag(:)

    i2max=max(nint(1./(zwidth*ni)),1)
    do i=1,ni
       do i2=1,i2max
          if (ran2(idum)<=zfrac) then
             x=(i-1)+(i2-0.5)/float(i2max)
             do j=1,nj
                do n=1,num(i,j)
                   if (flag(pnum(i,j,n))) cycle
                   if (abs(pos(pnum(i,j,n),1)-x)<zwidth*ni/2.) then
                      flag(pnum(i,j,n))=.true.
                      next=next+1
                   endif
                enddo
             enddo
          endif
       enddo
    enddo

  end subroutine zpad

end program make_starmask

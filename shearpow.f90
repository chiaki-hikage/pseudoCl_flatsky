program calc_shearspec

! calculate power spectrum from input shear field
!   infd: input file of shear data
!   inff: input file of flag data desribing the data is masked or not
!   outf: output file of shear power spectrum
!   addnoi: if addnoi="y", add noise to the shear data
!   rotate: if rotate="y", rotate shear to estimate noise power spectrum

  use fieldinfo_
  use sub_

  implicit none

  integer :: np,np2,np3,irun,idum
  real, pointer ::  pos(:,:),cl(:)
  complex, pointer ::  shear(:,:),sdata(:)
  character(100) :: infd,outf,inff
  character(1) :: addnoi,rotate
  character(10) :: arg
  logical*1, pointer :: flag(:)
  
  call getarg(1,infd)
  call getarg(2,outf)
  call getarg(3,addnoi)
  call getarg(4,inff)
  call getarg(5,rotate)
  call getarg(6,arg)
  read(arg,*) irun
  idum=-irun-100

  call binninginit

  call readsheardata(infd,pos,sdata)
  if (addnoi=="y") call addnoise(sdata,idum)
  if (rotate=="y") call rotateshear(sdata,idum)
  call readflag(inff,flag)
  call gridassign(pos,sdata,shear,flag)
  deallocate(pos,sdata,flag)
  allocate(cl(npmax3))
  call calcpow(shear,cl)
  deallocate(shear)

  open(2,file=outf,status='unknown')
  do np=1,npmax
     np2=np+npmax
     np3=np+npmax2
     write(2,"(4(1pe12.5,1x))") lbin(np),cl(np),cl(np2),cl(np3)
  enddo
  close(2)
  deallocate(cl)

contains

  subroutine readflag(inf,flag)

    logical*1, pointer :: flag(:)
    character, intent(in) :: inf*100

    allocate(flag(ntot))
    if (inf=="none") then
       flag=.false.
    else
       open(1,file=inf,status='old',form='unformatted')
       read(1) flag
       close(1)
    endif
    
  end subroutine readflag

  subroutine readsheardata(inf,pos,sdata)

    use fieldinfo_
    implicit none

    integer :: n
    complex, pointer :: sdata(:)
    real, pointer :: pos(:,:)
    character, intent(in) :: inf*100
    
    allocate(pos(ntot,2),sdata(ntot))
    open(1,file=inf,status='old',form='unformatted')
    read(1) pos
    read(1) sdata
    close(1)

  end subroutine readsheardata

  subroutine gridassign(pos,sdata,shear,flag)

    integer :: i,j,i1,j1,n
    real :: fx,fy
    real, pointer :: pos(:,:)
    complex, pointer :: shear(:,:),sdata(:)
    logical*1, pointer :: flag(:)

    allocate(shear(ni,nj))
    shear=(0.,0.)
    do n=1,ntot
       if (flag(n)) cycle
       i=int(pos(n,1))+1
       j=int(pos(n,2))+1
       shear(i,j)=shear(i,j)+sdata(n)
    enddo
    close(1)

  end subroutine gridassign
  
end program calc_shearspec

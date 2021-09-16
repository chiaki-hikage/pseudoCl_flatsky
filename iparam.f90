program initparam

  integer :: iarg
  character :: arg*10,outf*100


  outf='fieldinfo.f90'
  open(2,file=outf,status='unknown')

  write(2,'(a17)') 'MODULE fieldinfo_'
  write(*,*) 'Grid number of each side (e.g., 512)'
  read(*,*) arg

  write(2,*) 'integer, parameter :: ni='//arg(1:len_trim(arg))//',nj=ni'
  write(*,*) 'Side length [arcmin] (e.g., 300)'
  read(*,*) arg
  write(2,*) 'real, parameter ::  leni='//arg(1:len_trim(arg))//',lenj=leni' 
  write(2,*) 'real, parameter :: pi=3.14159265358979' 
  write(2,*) 'real, parameter :: rad1=leni/10800.*pi, rad2=rad1' 
  write(2,*) 'real, parameter :: amin2grid=float(ni)/leni'
  write(*,*) 'binning of power spectrum'
  write(*,*) '1. linear  2. logarithmic'
  read(*,*) iarg
  if (iarg==1) then
     arg='lin'
  else
     arg='log'
  endif
  write(2,*) 'character(3), parameter :: bkind="'//arg(1:len_trim(arg))//'"'
  write(*,*) 'lowest wavenumber [unit:fundamental mode 2pi/l] (e.g, 1.)'
  read(*,*)  arg
  write(2,*) 'real, parameter :: lmin=2.*pi/rad1*'//arg(1:len_trim(arg))
  write(*,*) 'highest wavenumber [unit:Nyquist freq. (2pi/l)*(ngrid/2)] &
       (e.g., 0.8)'
  read(*,*)  arg
  write(2,*)  'real, parameter :: lmax=2.*pi/rad1*(ni/2)*'//arg(1:len_trim(arg))
  write(*,*) 'number of bins (e.g., 30)' 
  read(*,*)  arg
  write(2,*) 'integer, parameter :: npmax='//arg(1:len_trim(arg))//',&
       npmax2=npmax*2,npmax3=npmax*3' 
  write(*,*) 'number density of source galaxies per square arcmin (e.g., 30)'
  read(*,*) arg
  write(2,*) 'real, parameter :: rhoasec='//arg(1:len_trim(arg))
  write(2,*) 'integer, parameter :: ntot=nint(rhoasec*leni*lenj)'
  write(*,*) 'intrinsic shear (e.g., 0.22)'
  read(*,*) arg
  write(2,*) 'real, parameter :: intshear='//arg(1:len_trim(arg))
  write(2,*) 'real :: lbin(npmax),dlbin(npmax),llowbin(npmax+1)' 
  write(2,'(a21)') 'end module fieldinfo_'

  close(2)

end program initparam

program convolve_or_deconvolve_field

! (de)convolve power spectrum using mode-coupling matrices
!   infp: input file of power spectrum
!   infn: input file of noise power spectrum (if "none", noise power set as 0)
!   outf: output file of (de)convolved power spectrum
!   infmsquare: input file of mode coupling matrix due to a square patch of sky
!   infminside: input file of mode coupling matrix due to a mask inside the square
!   ctype: if "d", deconvolve power, elseif "c", convolve power

  use fieldinfo_
  use sub_

  implicit none

  real, pointer :: cl(:),cln(:)
  double precision,dimension(npmax3,npmax3) :: mat,mat2
  character(100) :: infp,infn,infmsquare,infminside,outf
  character(1) :: ctype

  call getarg(1,infp)
  call getarg(2,infn)
  call getarg(3,outf)
  call getarg(4,infmsquare)
  call getarg(5,infminside)
  call getarg(6,ctype)

  call binninginit

  call readmat(infmsquare,mat)
  call readmat(infminside,mat2)

  allocate(cl(npmax3),cln(npmax3))
  call readcl3(infp,cl)
  call readcl3(infn,cln)
  cl=cl-cln
  if (ctype=="d") then
     call deconv_shear(cl,mat,mat2)
  elseif (ctype=="c") then
     call conv_shear(cl,mat,mat2)
  endif

  call outcl3(outf,cl)
  deallocate(cl,cln)
  
end program convolve_or_deconvolve_field

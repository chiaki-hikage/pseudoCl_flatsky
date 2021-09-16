MODULE maskinfo_

!   Simulate a mask for lensing data including
!
!    bright stars
!       circles with the radius ranging from "hmin" to "hmax"

  real, parameter :: hmin=0.2, hmax=2. ! in unit of arcmin

!    saturation spike
!       rectangles with xrec*radius by yrec*radius
!       for bright stars with radius > rsmin

  real, parameter :: xrec=0.2, yrec=5.
  real, parameter :: rsmin=0.3  ! in unit of arcmin

!    fmask: mask fraction by bright stars and saturation spike

  real, parameter :: fmask=0.2

!    other mask
!       zero padding in y direction with a probability of zfrac
!       width of zero paddign is given by zwidth

  real, parameter :: zfrac=0.05, zwidth=1./2048.

END MODULE maskinfo_

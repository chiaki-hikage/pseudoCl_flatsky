MODULE fieldinfo_
 integer, parameter :: ni=256,nj=ni
 real, parameter ::  leni=300,lenj=leni
 real, parameter :: pi=3.14159265358979
 real, parameter :: rad1=leni/10800.*pi, rad2=rad1
 real, parameter :: amin2grid=float(ni)/leni
 character(3), parameter :: bkind="lin"
 real, parameter :: lmin=2.*pi/rad1*0.5
 real, parameter :: lmax=2.*pi/rad1*(ni/2)*1.
 integer, parameter :: npmax=50,npmax2=npmax*2,npmax3=npmax*3
 real, parameter :: rhoasec=30
 integer, parameter :: ntot=nint(rhoasec*leni*lenj)
 real, parameter :: intshear=0.22
 real :: lbin(npmax),dlbin(npmax),llowbin(npmax+1)
end module fieldinfo_

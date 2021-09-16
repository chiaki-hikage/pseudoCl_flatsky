program collect_power

! collect power spectrum data
! 1st col: l
! 2nd col: input power (E-mode)
! 3rd col: deconvolved power (E-mode)
! 4th col: convolved power (E-mode)
! 5th col: masked power (E-mode)
! 6-9th cols: same as 2-5th cols but for B-mode power

  use fieldinfo_
  use sub_
  implicit none

  integer :: np,np2
  real, pointer :: cl1(:),cl2(:),cl3(:),cl4(:)
  character(100) :: inf1,inf2,inf3,inf4,outf

  call getarg(1,inf1)
  call getarg(2,inf2)
  call getarg(3,inf3)
  call getarg(4,inf4)
  call getarg(5,outf)

  call binninginit

  allocate(cl1(npmax3),cl2(npmax3),cl3(npmax3),cl4(npmax3))
  call readcl3(inf1,cl1)
  call readcl3(inf2,cl2)
  call readcl3(inf3,cl3)
  call readcl3(inf4,cl4)

  open(2,file=outf,status='unknown')
  do np=1,npmax
     np2=np+npmax
     write(2,'(9(1pe12.5,1x))') lbin(np),cl1(np),cl2(np),cl3(np),cl4(np), &
          cl1(np2),cl2(np2),cl3(np2),cl4(np2)
  enddo
  close(2)
  deallocate(cl1,cl2,cl3,cl4)

end program collect_power

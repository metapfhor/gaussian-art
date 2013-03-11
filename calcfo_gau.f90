subroutine GAUcalcforce(nat,posa,boxl,tmp_force, pot_energy)
  integer, intent(in)                          :: nat
  real(kind=8), intent(in), dimension(3*nat) :: posa
  real(kind=8), intent(inout), dimension(3*nat), target :: tmp_force

  call mkinput(nat,posa)

  call runGaussian()

  call readforce(nat,tmp_force)

end subroutine GAUcalcforce

subroutine mkinput(nat,posa)

  use defs
  implicit none

  integer, intent(in)                          :: nat
  real(kind=8), intent(in), dimension(3*nat) :: posa
  integer :: i
  character(len=80) :: temp


  call getenv('temporary_gaussian',temp)
  open(unit=FGAU,file='temp.com',status='replace')


  write(FGAU,*) '%chk=temp.chk'
  write(FGAU,*) '# opt rhf/3-21g nosymm force '

 do i=0,nat-1
  write(FGAU,'(8X,A2,4X,F16.10,F16.10,F16.10)')type_name(typat(i)),posa(3*i+1),posa(3*i+2),posa(3*i+3)
 enddo
 close(FGAU)

end subroutine mkinput

subroutine runGaussian()

  call system("perl getGaussianForces.pl")

end subroutine runGaussian

subroutine readforce(nat,force)

  use defs , only: FFORCES
  implicit none

  integer, intent(in)                          :: nat
  real(kind=8), intent(out), dimension(3*nat) :: force
  integer :: i

  open(unit=FFORCES,file='temp.forces',status='replace')
  do i=0,nat-1
    read(FFORCES,'(F16.10,F16.10,F16.10)')force(3*i+1),force(3*i+2),force(3*i+3)
  enddo
  close(FFORCES)

end subroutine readforce

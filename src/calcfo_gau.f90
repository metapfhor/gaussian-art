subroutine GAUcalcforce(nat,posa,boxl,tmp_force, pot_energy)
  integer, intent(in)                          :: nat
  real(kind=8), intent(in), dimension(3*nat) :: posa
  real(kind=8), intent(inout), dimension(3*nat), target :: tmp_force

  call mkinput(nat,posa)

  call runGaussian()

  call readforce(nat,tmp_force,pot_energy)

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
  open(unit=FXYZ,file='temp.xyz',status='unknown',access='append')

  write(FGAU,*) '%chk=temp.chk'
  write(FGAU,*) '#rhf/3-21g nosymm force '
  write(FGAU,*) ' '
  write(FGAU,*) ' '
  write(FGAU,*) ' '
  write(FGAU,*) '0 1'

  write(FXYZ,*) nat
  write(FXYZ,*) 'MOLECULAR TITLE'


 do i=1,nat
  write(FGAU,'(8X,A2,4X,F16.10,F16.10,F16.10)')type_name(typat(i)),posa(i),posa(nat+i),pos(2*nat+i)
  write(FXYZ,'(8X,A2,4X,F16.10,F16.10,F16.10)')type_name(typat(i)),posa(i),posa(nat+i),pos(2*nat+i)
 enddo
 write(FGAU,*) ' '
 close(FGAU)

end subroutine mkinput

subroutine runGaussian()

  call system("perl getGaussianForces.pl")

end subroutine runGaussian

subroutine readforce(nat,force,energy)

  use defs , only: FFORCES
  implicit none

  integer, intent(in)                          :: nat
  real(kind=8), intent(out), dimension(3*nat) :: force
  real(kind=8), intent(out) :: energy
  integer :: i
  character(len=80) :: temp


  open(unit=FFORCES,file='temp.forces',status='old')
    read(FFORCES,'(ES16.10)')energy
  do i=1,nat
    read(FFORCES,'(ES16.10,A3,ES16.10,A3,ES16.10)')force(i),temp,force(nat+i),temp,force(2*nat+i)
  enddo
  close(FFORCES)

end subroutine readforce

!> @file
!!   Subroutine which calls the right for type inside art
!! @author
!!   Written by Laurent Karim Beland, UdeM 2011!!
!!   Copyright (C) 2010-2011 BigDFT group, Normand Mousseau
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Subroutine which calls the right for type inside art
subroutine calcforce(nat, posa, boxl, forca, energy, evalf_number, conv )

   use defs, only : energy_type
   implicit none

   !Arguments
   integer,      intent(in)                            :: nat
   real(kind=8), intent(in),  dimension(3*nat)         :: posa
   real(kind=8), dimension(3), intent(inout)           :: boxl
   real(kind=8), intent(out), dimension(3*nat)         :: forca
   real(kind=8), intent(out)                           :: energy
   integer,      intent(inout)                         :: evalf_number
   logical,      intent(in)                            :: conv


   if(energy_type == "SWP")  then
    call SWcalcforce(nat,posa,boxl,forca, energy)
    evalf_number = evalf_number +1
   elseif(energy_type== "GAU") then
   	call GAUcalcforce(nat,posa,boxl,forca, energy)
		evalf_number = evalf_number +1
   endif

END SUBROUTINE calcforce

!> @file
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! Modified by:
!! -EM 2010, see ~/AUTHORS
!! -Laurent Karim Beland, UdeM, 2011. For working with QM/MM !!
!! -EM 2011, see ~/AUTHORS
!!
!> ART initialize
!!   Initialization of art method
!!   It relaxes it into a local minimum without volume optimization
!!   The initial configuration is a "reference configuration", it will be the reference 
!!   configuration until a new event is accepted. 
!! 
subroutine initialize()

  use defs
  use lanczos_defs,  only : LANCZOS_MIN 
  use saddles,       only : g_pos, GUESSFILE
  implicit none

  !Local variables
  integer :: i, ierror
  character(len=20) :: dummy, fname
  character(len=4)  :: scounter
  logical           :: flag, success

  integer           :: nat_test
  integer, allocatable               :: typ_a(:)    ! Atomic type
  integer, allocatable               :: const_a(:)  ! Constraints
  real(kind=8), allocatable, target  :: pos_a(:)    ! Working positions of the atoms
  real(kind=8),dimension(3)          :: boxref_     ! Reference box from posinp file

  character(len=1)                    :: boundary_b
  integer, allocatable                :: typ_b(:)    ! Atomic type
  integer, allocatable                :: const_b(:)  ! Constraints
  real(kind=8), allocatable, target   :: pos_b(:)    ! Working positions of the atoms
  real(kind=8), dimension(:), pointer :: xb, yb, zb
  real(kind=8), dimension(3) :: boxl
  !_______________________
  
  ! Read the counter in order to continue the run where it stopped or for refine
  ! Format:
  ! Counter:     1000

  ! Laurent Modification: Open up temp.xyz to clean it out
  open(unit=FXYZ,file='temp.xyz',status='replace',access='sequential',position='rewind')
  !write(FXYZ,*)''
  !flush(FXYZ)
  !close(FXYZ)

  if (.not. restart) then 
     inquire( file = COUNTER, exist = flag )
     if ( flag .and. iproc == 0 ) then 
        open(unit=FCOUNTER,file=COUNTER,status='old',action='read',iostat=ierror)
        read(FCOUNTER,'(A12,I6)') dummy, mincounter 
        close(FCOUNTER)
     else
        mincounter = 1000
     end if
  end if

  ! we read the initial/reference configuration
  if ( new_event .or. restart ) then   
     ! Read initial atomic file
     ! is it neccesary for restart ?
     
     nat_test = NATOMS
     allocate(typ_a(NATOMS))
     allocate(pos_a(3*NATOMS))
     allocate(const_a(NATOMS))
     call init_all_atoms( nat_test, typ_a, pos_a, const_a, boxref_, boundary, nproc, iproc, trim(REFCONFIG) )
  else if ( eventtype == 'REFINE_AND_RELAX' .or. eventtype == 'REFINE_SADDLE' ) then
     ! Read reference atomic file

     nat_test = NATOMS
     allocate(typ_a(NATOMS))
     allocate(pos_a(3*NATOMS))
     allocate(const_a(NATOMS))
     call init_all_atoms( nat_test, typ_a, pos_a, const_a, boxref_, boundary, nproc, iproc, trim(REFCONFIG) )
  end if
    
  !assign the data from the atomic file
  if ( .not. restart ) then
     typat(:)   = typ_a(:)
     pos(:)     = pos_a(:)
     constr(:)  = const_a(:)
     boxref(:)  = boxref_(:)
     refcounter = mincounter
     box = boxref
     
     deallocate(typ_a)
     deallocate(pos_a)
     deallocate(const_a)
     
     ! We read the position for the presumed saddle point.
     if ( eventtype == "GUESS_DIRECTION" ) then 

        nat_test = NATOMS
        allocate(typ_b(NATOMS))
        allocate(pos_b(3*NATOMS))
        allocate(const_b(NATOMS))
        call init_all_atoms( nat_test, typ_b, pos_b, const_b, boxref_, boundary_b, nproc, iproc, trim(GUESSFILE) )
        
        ! Let's check if it is in the same conditions as posinp.
        write(*,*) "eventtype: ", eventtype
        if ( boundary /= boundary_b ) then
           if ( iproc == 0 ) write(*,*) "GUESS: Different type of boundary conditions"
           call end_art()
        end if
        do i = 1, NATOMS
           if ( typ_b(i) /= typat(i) ) then
              if ( iproc == 0 ) write(*,*) "GUESS: Different type of atoms"
              call end_art()
           end if
           if ( const_b(i) /= constr(i) ) then
              if ( iproc == 0 ) write(*,*) "GUESS: Different type of constraints"
              call end_art()
           end if
        end do
        ! g_pos in module saddles
        allocate(g_pos(vecsize))
        g_pos(:) = 0.0d0
        g_pos(:) = pos_b(:)
        
        deallocate(typ_b)
        deallocate(pos_b)
        deallocate(const_b)
     end if
  endif

  ! write
  if ( iproc == 0 ) then 
   open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
   write(FLOG,'(1X,A34,I17)') ' - Mincounter                   : ', mincounter
   close(FLOG)
  end if

  ! We rescale the coordinates. For what ?? 
  scalaref = 1.0d0
  scala = scalaref

  call initialize_potential()         ! Initialize Potential (CORE) 

                                      ! for output files
  if ( iproc == 0 ) call convert_to_chain( refcounter, 4, scounter )
  fname = FINAL // scounter
  conf_initial = fname
                                      ! If this is a new event we relax 
  If_ne: if ( new_event .and. (.not. restart) ) then 
     write(*,*) "go into min_converge"
     call min_converge( success )     ! Converge the configuration to a local minimum
     write(*,*) "success: ", success

     posref = pos                     ! New reference configuration.
     ref_energy = total_energy
 
     if ( iproc == 0 ) then
                                      ! THIS iS NOT NECESSARY IN BIGDFT
        call write_refconfig( )       ! Write reference in REFCONFIG. 
        call store( fname )           ! Store the configuration into fname.

        open( unit = FLOG, file = LOGFILE, status = 'unknown',&
        & action = 'write', position = 'append', iostat = ierror )
        write(*,*) 'BART: Configuration stored in file (initialize.f90)',fname ,' success: ',success
        write(FLOG,'(1X,A34,A17)') ' - Configuration stored in file : ', trim(fname)
        if ( .not. success ) then
          write(FLOG,'(1X,A)') "ERROR: Initial configurations is not a minimum"
          call end_art()  
        end if
        close(FLOG)
     end if
     mincounter = mincounter + 1
                                      ! if dual_search we dont do this check at the
                                      ! beginning. It is not well defined.
     if ( success .and. .not. dual_search ) then
        if ( LANCZOS_MIN .or. setup_initial ) call check_min( 'M' ) 
     end if 

  else if ( (.not. new_event) .and. (.not. restart) ) then

                                      ! once we have the total energy we copy as reference values
     posref = pos
     ref_energy = total_energy

     if ( iproc == 0 ) then
        call store( fname )           ! Store the reference configuration into fname.

        open( unit = FLOG, file = LOGFILE, status = 'unknown',&
        & action = 'write', position = 'append', iostat = ierror )
        write(*,*) 'BART: Ref. Configuration stored in file ',fname
        write(FLOG,'(1X,A34,A17)') ' - Configuration stored in file : ', trim(fname)
        close(FLOG)
     end if
     mincounter = mincounter + 1

     ! now we read the positions of the configuration we want to refine it
     nat_test = NATOMS
     allocate(typ_b(NATOMS))
     allocate(pos_b(3*NATOMS))
     allocate(const_b(NATOMS))
     call init_all_atoms( nat_test, typ_b, pos_b, const_b, boxref_, boundary_b, nproc, iproc, 'posinp' )

     ! Let's check if it is in the same conditions as the reference.
     if ( boundary /= boundary_b ) then
        if ( iproc == 0 ) write(*,*) "posinp: Different boundary condition"
        call end_art()
     end if

     do i = 1, NATOMS
        if ( typ_b(i) /= typat(i) ) then
           if ( iproc == 0 ) write(*,*) "posinp: Different type of atoms"
           call end_art()
        end if

        if ( const_b(i) /= constr(i) ) then
           if ( iproc == 0 ) write(*,*) "posinp: Different type of constraints"
           call end_art()
        end if

     end do
                                      ! our real starting point
     pos(:)  = pos_b(:)
     box(:)  = boxref_(:)

     write(*,*) 'box: ', box
     write(*,*) 'positions: ', pos(1), pos(1+NATOMS), pos(1+2*NATOMS)

     deallocate(typ_b)
     deallocate(pos_b)
                                      ! we need the total energy for this configuration
     boxl = box * scala               ! we compute at constant volume.

     call calcforce( NATOMS, pos, boxl, force, total_energy, evalf_number, .false. )

  end if If_ne

END SUBROUTINE initialize

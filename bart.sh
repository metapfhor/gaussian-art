#!/bin/bash
#____________________________ ATOMS 
export gnrm=1.0E-05
export Inflection=4 
export N_atoms=8                        # Number of atoms in the problem
export Type_1=C                           # Element type for the xyz output. Add type2, etc. if needed
export Type_2=H
#_____________________________ ART 
export Event_type=NEW                  # Either 'NEW', 'REFINE_SADDLE' when further converging a saddle point
                                        # Or "REFINE_AND_RELAX", to refine at the saddle
			                            # and check the final minimum
export Temperature=0.5   # Temperature in eV, if negative always reject the event
export Max_Number_Events=100   # Maximum number of events
export Type_of_Events=local   # Initial move for events - global or local
export Radius_Initial_Deformation=3.5   # Cutoff for local-move (in angstroems)
# export Central_Atom=241    # Number of the atom around which the initial move takes place
export sym_break_dist=0.001   # Breaks the symmetry of the crystal by randomly displacing
                                               # all atoms by this distance
export Forcefield=GAU    # Choice of forcefield - currently allows SWP 


export Activation_MaxIter=100   # Maximum number of iteraction for reaching the saddle point
export Increment_Size=0.08   # Overall scale for the increment moves
export Force_Threshold_Perp_Rel=0.5   # Threshold for perpendicular relaxation

#_____________________________ HARMONIC WELL
export Initial_Step_Size=0.05   # Size of initial displacement, in A
export Basin_Factor=2.4   # Factor multiplying Increment_Size for leaving the basin
export Max_Perp_Moves_Basin=3   # Maximum number of perpendicular steps leaving basin
export Min_Number_KSteps=2   # Min. number of ksteps before calling lanczos 
export Eigenvalue_Threshold=-1.5   # Eigenvalue threshold for leaving basin
export Max_Iter_Basin=20   # Maximum number of iteraction for leaving the basin (kter)
#_____________________________ LANCZOS
export Lanczos_of_minimum=.False.    # Calculation of the Hessian for each minimum
export Number_Lanczos_Vectors=15   # Number of vectors included in lanczos procedure
export delta_disp_Lanczos=0.00025   # Step of the numerical derivative of forces in lanczos (Ang)
#_____________________________ CONVERGENCE
export Exit_Force_Threshold=0.10   # Threshold for convergence at saddle point
export Prefactor_Push_Over_Saddle=0.15   # Fraction of displacement over the saddle
export Save_Conf_Int=.False.   # Save the configuration at every step?
#_____________________________ DIIS
export Iterative=.False.   # Iterative use of Lanczos & DIIS
export Use_DIIS=.False.   # Use DIIS for the final convergence to saddle
export DIIS_Force_Threshold=1.5   # Force threshold for call DIIS
export DIIS_Memory=5   # Number of vectors kepts in memory for algorithm
export DIIS_Check_Eigenvector=.True.   # Check that the final state is indeed a saddle
export DIIS_Step_size=0.03   # prefactor multiplying forces
export FACTOR_DIIS=5.0   # times Increment_Size, max allowed diis step size
export MAX_DIIS=150   # max diis iterations per call
#_____________________________ INPUT 
export FILECOUNTER=filecounter   # File tracking  the file (event) number - facultative
export REFCONFIG=ethylene  # Reference configuration (actual local minimum)
#_____________________________ OUTPUT 
export LOGFILE=log.file   # General output for message
export EVENTSLIST=events.list   # list of events with success or failure
export Write_restart_file=.False.   # It is useful only for ab-initio
export RESTART=restart.dat   # current data for restarting event
export Write_JMOL=.False.   # Writes the configuration in jmol format.
###### RUN THE SIMULATION #######################################################################
gaussian-art/arttest.exe

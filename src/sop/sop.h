#undef finnis_sinclair

!impose Dirichlet boundary conditions
!#define dirichlet_particles 1

!Account for the variations in D when computing the gradient of Psi
#undef dPsi_dD

!# define pressure_stuff
!debug_verbosity
! 0 -> production runs
! 1 -> write little, on std output
! 2 -> write more, on std output (incl. gradient descent files)
! 3 -> write more, on std output AND dump some files
!#define debug_verbosity 3
#define debug_verbosity 2

!writeout verbosity
! 0 -> production runs (writeout final state of Phi only)
! 1 -> writeout derivatives too
! 2 -> writeout D, rcp, etc...
! 3 -> writeout nvlist and files during gradient descent
#define writeout_verbosity 1

!to writeout .xyz files with high precision
#define highprec_writeout 1

!to always compute first derivatives, even when not strictly needed
#define compute_all 1

!mainly for debugging purposes: makes sure that the (random) initial
!particle positions does not depend on the number of processors
!(not a memory-efficient implementation at all)
#define same_random_sequence_nproc 1


!Use random numbers to create/insert new particles
!These random numbers are generated locally (on each proc)
!This means that using a different number of processors lead to different results
#define __USE_RANDOMNUMBERS 1
!Rubbish implementation, use only for debugging 
!#undef __USE_RANDOMNUMBERS

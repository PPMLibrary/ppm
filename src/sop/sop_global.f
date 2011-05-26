!======================================================================!
! contains global variables (for all procedures using this module)
!======================================================================!
MODULE sop_global

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

#define __KIND __DOUBLE_PRECISION

#include "sop.h"

USE ppm_module_user
USE ppm_module_typedef
USE ppm_module_ctrl

#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER,PRIVATE :: prec = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER,PRIVATE :: prec = ppm_kind_double
#endif



TYPE sop_t_diff_eq
    !!! derived type to define a differential operator
    !!! The latter is assumed to be a linear combination
    !!! of multivariate derivatives

    INTEGER                           :: nterms
    !!! number of terms
    REAL(prec),DIMENSION(:),  POINTER :: coeffs
    !!! coefficients in front of each term
    REAL(prec),DIMENSION(:,:),POINTER :: order_d
    !!! degree of the derivative for each term

END TYPE sop_t_diff_eq

TYPE sop_t_opts
    !!! derived type for optional arguments to the adaptivity routine

    LOGICAL           :: level_set
    REAL(prec)        :: scale_D
    REAL(prec)        :: param_nb
    REAL(prec)        :: nb_width
    REAL(prec)        :: nb_width2
    REAL(prec)        :: nb_width_kill
    REAL(prec)        :: adaptivity_criterion
    REAL(prec)        :: fuse_radius
    INTEGER           :: order_approx
    !!! order of approximation for the interpolation kernels
    TYPE(sop_t_diff_eq),POINTER :: diff_eq
    !!! RHS of the differential equation that needs to be computed
    !!! at the interpolation step

END TYPE sop_t_opts


LOGICAL                  :: write_pdb    ! writeout pdb files
LOGICAL                  :: write_xyz    ! writeout xyz files

!====================================================================!
! Variable describing the 'state' of the system
!====================================================================!
LOGICAL            :: Psi_saturates ! true when the forces have
!          to be truncated to prevent particles from flying across the domain
!                (happens only when 2 particles are very close to each other)

REAL(prec)                         :: Psi_global ! mean interaction 
!                   potential of the particles
REAL(prec)                         :: Psi_global_old 
REAL(prec)                         :: Psi_max  
REAL(prec)                         :: rcp_over_D

!cutoff that was given as argument to mktopo
! the size of the ghostlayers should never be bigger than this
! (without re-constructing the topology)
REAL(prec)                         :: cutoff_topology


!====================================================================!
! random numbers
!====================================================================!
INTEGER                               :: seedsize
INTEGER,  DIMENSION(:), ALLOCATABLE   :: seed
REAL(prec), DIMENSION(:), ALLOCATABLE :: randnb
INTEGER                               :: randnb_i

!====================================================================!
! verlet lists
!====================================================================!
INTEGER                                             :: nneigh_critical
INTEGER                                             :: nneigh_theoretical
! number of neighbours above which something must be wrong
INTEGER                                             :: nneigh_toobig

!====================================================================!
! maths
!====================================================================!
REAL(prec),PARAMETER :: PI = 2._prec*ACOS(0._prec)
REAL(prec)           :: PI2

!====================================================================!
! tolerances in numerical schemes
!====================================================================!
REAL(prec)           :: tolerance_moment_conditions
REAL(prec)           :: tolerance_LSE
REAL(prec)           :: tolerance_singular_values
REAL(prec)           :: tolerance_livecheck

!====================================================================!
! numeric parameters
!====================================================================!
REAL(prec)                            :: scale_D        ! resolution scale
! maximum scaling distance allowed
REAL(prec)                            :: maximum_D      
! minimum resolution
REAL(prec)                            :: minimum_D      
! stopping criterion for particle adaptation
REAL(prec)                            :: adaptivity_criterion    
!counter for the number of gradient descent steps performed
INTEGER                               :: nb_grad_desc_steps 
!separation distance below which 2 particles fuse
REAL(prec)                            :: fuse_radius 

!some parameters used in various functions
REAL(prec)                            :: param_a !used in module_funcs
REAL(prec)                            :: param_d0, param_d1,param_p0 

REAL(prec)                            :: attractive_radius0 
!distance below which particles attract each other

!private:: define_sop_args

CONTAINS

SUBROUTINE define_sop_args

    CALL arg(fuse_radius, 'fuse_radius', &
        ctrl_name = 'fuse_radius',       &
        long_flag      =  '--sop-fuse-r',                 &
        default   =  0._prec )
        
    CALL arg(attractive_radius0, 'attractive_radius0', &
        ctrl_name = 'attractive_radius',       &
        long_flag      =  '--sop-attract-rad',                 &
        default   =  0._prec )
        
    CALL arg(adaptivity_criterion, 'adaptivity_criterion', &
        ctrl_name = 'adaptivity_criterion',       &
        long_flag      =  '--sop-adapt-crit',                 &
        default   =  8._prec )
        
    CALL arg(nneigh_theoretical, 'nneigh_theoretical', &
        ctrl_name = 'nneigh_theoretical',       &
        long_flag      =  '--nneigh-theoretical',                 &
        default   =  24 )
        
    CALL arg(rcp_over_D, 'rcp_over_D', &
        ctrl_name = 'rcp_over_D',     &
        long_flag      =  '--rcp-over-D',                 &
        default   =  2._prec ,         &
        min       =  0._prec           )
        
    CALL arg(maximum_D, 'maximum_D', &
        ctrl_name = 'maximum_D',     &
        long_flag      =  '--sop-max-D',                 &
        default   =  1._prec ,         &
        min       =  0._prec           )
        
    CALL arg(scale_D, 'scale_D', &
        ctrl_name = 'scale_D',     &
        long_flag      =  '--scale-D',                 &
        default   =  1._prec ,         &
        min       =  0._prec           )
        
    CALL arg(minimum_D, 'minimum_D', &
        ctrl_name = 'minimum_D',     &
        long_flag      =  '--sop-min-D',                 &
        default   =  0._prec ,         &
        min       =  0._prec            )
        
    CALL arg(write_pdb, 'write_pdb', &
        ctrl_name = 'write_pdb',     &
        long_flag      =  '--sop-write-pdb',                 &
        default   =  .FALSE.)
        
    CALL arg(write_xyz, 'write_xyz', &
        ctrl_name = 'write_xyz',     &
        long_flag =  '--sop-write-xyz',                 &
        default   =  .FALSE.)
        
END SUBROUTINE define_sop_args


END MODULE sop_global

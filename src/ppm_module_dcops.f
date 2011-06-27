!--------------------------------------------------------------------------
!  Module       :                    ppm_module_dcops
!--------------------------------------------------------------------------
!
!  Purpose      : Compute generalized DC operators to arbitrary order
!
!  Remarks      :
!
!  References   :
!
!  Revisions    :
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
!  Sylvain Reboux
!  ETH Zurich
!  CH-8092 Zurich, Switzerland
!--------------------------------------------------------------------------

!-------------------------------------------------------------------------
!  Define types
!-------------------------------------------------------------------------
#define __SINGLE_PRECISION         1
#define __DOUBLE_PRECISION         2
#define __INTEGER                  3
#define __LOGICAL                  4
#define __SINGLE_PRECISION_COMPLEX 5
#define __DOUBLE_PRECISION_COMPLEX 6
#define __UNIFORM 1
#define __ADAPTIVE 2

!temporary hack.
!Comment out (or remove...) once ppm_module_particles
! (which contains the derived type ppm_t_particles)
!can handle single_precision
!#define __PARTICLES__FINISHED_SINGLE_PREC 1


MODULE ppm_module_dcops

!!! This module provides routines that compute the weights of the DC operators
!!! that approximate partial derivatives of any degree with any order of
!!! consistency.
!!! The routines are not overloaded for 2D and 3D and must be called
!!! directly.
!!! Four different flavours exist: uniform vs adaptive and with vs without
!!! the ppm_t_particles data structure
!!! The difference between uniform and adaptive lies solely in the cutoff
!!! parameter (a scalar in one case and a vector in the other). It is used
!!! within the routine as a scaling factor for the inter-particle distances.
!!! 
!!! The routines return the coefficients of the DC operators
!!! (stored in the array eta). By default, it is assumed that the particles
!!! are irregularly distributed and that one wants to compute a given
!!! partial derivative (the differentiation order for each independant variable
!!! is given by the array order_deriv). If the optional flag islaplacian is 
!!! set to true, order_deriv is not used and an approximation of the Laplacian
!!! operator is returned instead. If iscartesian is set to true, particles
!!! are assumed to be located on the nodes of a Cartesian grid (there are then
!!! fewer coefficients to compute and the routine is faster).
!!! When isinterp is true, it is expected that two sets of particles are
!!! given on input. One of them contains the field values. The routine computes
!!! the operators that approximate the derivative of this field when evaluated
!!! on the second set of particles. If order_deriv is set to zero, one obtains
!!! the operators that interpolate the field from one set of particles to the
!!! other. (NOT available without ppm_t_particles, for now).
!!! Some usage examples can be found in the test client

!this module is compiled only if either the BLAS libraries or the MKL libraries
!can be found. It is left empty otherwise
#ifdef __DCOPS

    !-------------------------------------------------------------------------
    !  Modules
    !-------------------------------------------------------------------------
    USE ppm_module_data, ONLY: ppm_rank,ppm_dim
    USE ppm_module_typedef
    USE ppm_module_write
    USE ppm_module_error
    USE ppm_module_substart
    USE ppm_module_substop

    INTERFACE ppm_dcops_2d
        MODULE PROCEDURE ppm_dcops_adapt_2d_s
        MODULE PROCEDURE ppm_dcops_unif_2d_s
        MODULE PROCEDURE ppm_dcops_adapt_2d_d
        MODULE PROCEDURE ppm_dcops_unif_2d_d
    END INTERFACE

    INTERFACE ppm_dcops_3d
        MODULE PROCEDURE ppm_dcops_adapt_3d_s
        MODULE PROCEDURE ppm_dcops_unif_3d_s
        MODULE PROCEDURE ppm_dcops_adapt_3d_d
        MODULE PROCEDURE ppm_dcops_unif_3d_d
    END INTERFACE

    INTERFACE solveLSE
        MODULE PROCEDURE solveLSEs
        MODULE PROCEDURE solveLSEd
    END INTERFACE

    INTERFACE solveLSE_2
        MODULE PROCEDURE solveLSE_2s
        MODULE PROCEDURE solveLSE_2d
    END INTERFACE

    INTERFACE solveLSE_n
        MODULE PROCEDURE solveLSE_ns
        MODULE PROCEDURE solveLSE_nd
    END INTERFACE

    INTERFACE primitive
        MODULE PROCEDURE primitive_s
        MODULE PROCEDURE primitive_d
    END INTERFACE

    INTERFACE ppm_matrix_svd
        MODULE PROCEDURE ppm_matrix_svd_s
        MODULE PROCEDURE ppm_matrix_svd_d
    END INTERFACE

    !INTERFACE ppm_dcop_check_vandermonde
        !MODULE PROCEDURE ppm_dcop_check_vandermonde_s
        !MODULE PROCEDURE ppm_dcop_check_vandermonde_d
    !END INTERFACE

!-------------------------------------------------------------------------
! Public subroutine
!-------------------------------------------------------------------------
PUBLIC :: solveLSE,solveLSE_2,solveLSE_n, &
        & ppm_dcops_2d,ppm_dcops_3d,      &
        & ppm_part_dcops_2d,ppm_part_dcops_3d

CONTAINS

#define __KIND __DOUBLE_PRECISION
#include "dcop/ppm_dcops_helpers.f"

#define __KIND __DOUBLE_PRECISION
#define __DIM 2
#define __MODE __ADAPTIVE
#include "dcop/ppm_dcops.f"
#define __KIND __DOUBLE_PRECISION
#define __DIM 2
#define __MODE __UNIFORM
#include "dcop/ppm_dcops.f"

#define __KIND __DOUBLE_PRECISION
#define __DIM 3
#define __MODE __ADAPTIVE
#include "dcop/ppm_dcops.f"
#define __KIND __DOUBLE_PRECISION
#define __DIM 3
#define __MODE __UNIFORM
#include "dcop/ppm_dcops.f"

#define __KIND __DOUBLE_PRECISION
#define __DIM 2
#include "dcop/ppm_part_dcops_build.f"
#define __KIND __DOUBLE_PRECISION
#define __DIM 3
#include "dcop/ppm_part_dcops_build.f"

#define __KIND __DOUBLE_PRECISION
#define __DIM 2
#include "dcop/ppm_part_dcops_build2.f"
#define __KIND __DOUBLE_PRECISION
#define __DIM 3
#include "dcop/ppm_part_dcops_build2.f"

#define __KIND __DOUBLE_PRECISION
#include "dcop/ppm_part_dcops.f"


#define __KIND __SINGLE_PRECISION
#include "dcop/ppm_dcops_helpers.f"

#define __KIND __SINGLE_PRECISION
#define __DIM 2
#define __MODE __ADAPTIVE
#include "dcop/ppm_dcops.f"
#define __KIND __SINGLE_PRECISION
#define __DIM 2
#define __MODE __UNIFORM
#include "dcop/ppm_dcops.f"

#define __KIND __SINGLE_PRECISION
#define __DIM 3
#define __MODE __ADAPTIVE
#include "dcop/ppm_dcops.f"
#define __KIND __SINGLE_PRECISION
#define __DIM 3
#define __MODE __UNIFORM
#include "dcop/ppm_dcops.f"

#ifdef __PARTICLES__FINISHED_SINGLE_PREC
#define __KIND __SINGLE_PRECISION
#define __DIM 2
#include "dcop/ppm_part_dcops_build.f"
#define __KIND __SINGLE_PRECISION
#define __DIM 3
#include "dcop/ppm_part_dcops_build.f"

#define __KIND __SINGLE_PRECISION
#define __DIM 2
#include "dcop/ppm_part_dcops_build2.f"
#define __KIND __SINGLE_PRECISION
#define __DIM 3
#include "dcop/ppm_part_dcops_build2.f"
#endif

#undef __KIND

SUBROUTINE particles_dcop_compute(Particles,eta_id,info,c,min_sv)

    USE ppm_module_data, ONLY: ppm_dim,ppm_rank
    USE ppm_module_particles_typedef
    USE ppm_module_write
    IMPLICIT NONE

    !---------------------------------------------------------
    ! arguments
    !---------------------------------------------------------
    TYPE(ppm_t_particles),   POINTER,    INTENT(INOUT)   :: Particles
    !!! particles
    INTEGER,                             INTENT(IN   )   :: eta_id
    !!! id of the operator kernel
    INTEGER,                             INTENT(  OUT)   :: info
    !!! non-zero on output if some error occurred
    !---------------------------------------------------------
    ! Optional arguments
    !---------------------------------------------------------
    REAL(ppm_kind_double),OPTIONAL                       :: c
    !!! ratio h/epsilon (default is 1.0)
    REAL(ppm_kind_double),OPTIONAL   ,  INTENT(  OUT)    :: min_sv
    !!! smallest singular value
    !---------------------------------------------------------
    ! local variables
    !---------------------------------------------------------
    CHARACTER(LEN = ppm_char)               :: caller = 'particles_dcop_compute'
    CHARACTER(LEN = ppm_char)               :: cbuf
    REAL(KIND(1.D0))                        :: t0
    LOGICAL                                 :: interp

    !-------------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------------
    info = 0 ! change if error occurs
    CALL substart(caller,t0,info)

    !-------------------------------------------------------------------------
    ! Check arguments
    !-------------------------------------------------------------------------
    IF (.NOT. ASSOCIATED(Particles)) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            & 'Particles not defined',&
            __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT. ASSOCIATED(Particles%ops)) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            & 'No operator data structure found, use particles_dcop_define first',&
            __LINE__,info)
        GOTO 9999
    ENDIF
    IF (eta_id.LE.0 .OR. eta_id.GT.Particles%ops%max_opsid) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            & 'Invalid eta_id, use particles_dcop_define first',&
            __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT. ASSOCIATED(Particles%ops%desc(eta_id)%degree)) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            & 'Operator not found, use particles_dcop_define first',&
            __LINE__,info)
        GOTO 9999
    ENDIF
    IF (Particles%ops%desc(eta_id)%is_computed) THEN
        WRITE(cbuf,*) 'WARNING: The operator with id ',eta_id,&
            & ' and name *',TRIM(ADJUSTL(Particles%ops%desc(eta_id)%name)),&
            &'* seems to have already been computed. Unnecessary call to',&
            &' particles_dcop_compute()'
        CALL ppm_write(ppm_rank,caller,cbuf,info)
    ENDIF

    interp = Particles%ops%desc(eta_id)%interp
    IF (interp) THEN
        IF (.NOT. ASSOCIATED(Particles%Particles_cross)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,&
                'Need to specify which set of particles &
                &   (particles_cross) should be used for interpolation',&
                & __LINE__,info)
            GOTO 9999
        ENDIF
        IF (.NOT. (Particles%neighlists_cross)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,&
                'Please compute xset neighbor lists first',&
                __LINE__,info)
            GOTO 9999
        ENDIF
    ELSE
        IF (.NOT. Particles%neighlists) THEN
            info = ppm_error_error
            CALL ppm_error(999,caller,   &
                & 'Compute neighbor lists first',&
                __LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF

    !-------------------------------------------------------------------------
    ! Compute the DC operator
    !-------------------------------------------------------------------------
    IF (ppm_dim .EQ. 2) THEN
        CALL ppm_dcop_compute2d(Particles,eta_id,info,interp,c,min_sv)
    ELSE
        CALL ppm_dcop_compute3d(Particles,eta_id,info,interp,c,min_sv)
    ENDIF
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            & 'ppm_part_dcops failed',&
            __LINE__,info)
        GOTO 9999
    ENDIF

    !-------------------------------------------------------------------------
    ! Update states
    !-------------------------------------------------------------------------
    Particles%ops%desc(eta_id)%is_computed = .TRUE.

    !-------------------------------------------------------------------------
    ! Finalize
    !-------------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE particles_dcop_compute



FUNCTION factorial_m(multi_ind,ndim)
    INTEGER                 :: factorial_m,i,ndim
    INTEGER,DIMENSION(ndim) :: multi_ind

    factorial_m = 1
    DO i=1,ndim
        factorial_m = factorial_m * factorial(multi_ind(i))
    ENDDO
END FUNCTION factorial_m

FUNCTION factorial(n)
    INTEGER    :: n,i
    INTEGER    :: factorial
    factorial = 1
    DO i=1,n
        factorial = i*factorial
    ENDDO
    RETURN
END FUNCTION factorial

FUNCTION binomial(n,k)
    INTEGER    :: n,k
    INTEGER    :: binomial
    binomial = factorial(n)/(factorial(k)*factorial(n-k))
END FUNCTION binomial

#endif

END MODULE ppm_module_dcops

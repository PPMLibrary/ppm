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

#define interface_s_d(a) \
    INTERFACE a ;\
        MODULE PROCEDURE a/**/_s; \
        MODULE PROCEDURE a/**/_d; \
    END INTERFACE

    interface_s_d(solveLSE)
    interface_s_d(solveLSE_2)
    interface_s_d(solveLSE_n)
    interface_s_d(primitive)
    interface_s_d(ppm_matrix_svd)
    interface_s_d(particles_dcop_compute)

    !interface_s_d(ppm_dcop_check_vandermonde)

!-------------------------------------------------------------------------
! Public subroutine
!-------------------------------------------------------------------------
PUBLIC :: solveLSE,solveLSE_2,solveLSE_n, &
        & ppm_dcops_2d,ppm_dcops_3d,      &
        & ppm_part_dcops_2d,ppm_part_dcops_3d, &
        & particles_dcop_compute

CONTAINS

#define __KIND __SINGLE_PRECISION
#define  DTYPE(a) a/**/_s
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single

#include "dcop/ppm_dcops_helpers.f"

#define __DIM 2
#define __MODE __ADAPTIVE
#include "dcop/ppm_dcops.f"

#define __DIM 2
#define __MODE __UNIFORM
#include "dcop/ppm_dcops.f"

#define __DIM 3
#define __MODE __ADAPTIVE
#include "dcop/ppm_dcops.f"

#define __DIM 3
#define __MODE __UNIFORM
#include "dcop/ppm_dcops.f"

#define __DIM 2
#include "dcop/ppm_dcop_compute.f"

#define __DIM 3
#include "dcop/ppm_dcop_compute.f"

#undef  __KIND
#undef  DTYPE
#undef  DEFINE_MK


#define __KIND __DOUBLE_PRECISION
#define  DTYPE(a) a/**/_d
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double

#include "dcop/ppm_dcops_helpers.f"

#define __DIM 2
#define __MODE __ADAPTIVE
#include "dcop/ppm_dcops.f"

#define __DIM 2
#define __MODE __UNIFORM
#include "dcop/ppm_dcops.f"

#define __DIM 3
#define __MODE __ADAPTIVE
#include "dcop/ppm_dcops.f"

#define __DIM 3
#define __MODE __UNIFORM
#include "dcop/ppm_dcops.f"

#define __DIM 2
#include "dcop/ppm_dcop_compute.f"

#define __DIM 3
#include "dcop/ppm_dcop_compute.f"

#undef __KIND
#undef  DTYPE
#undef  DEFINE_MK




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

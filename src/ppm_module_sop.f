MODULE ppm_module_sop
!!! Self-Organizing Particles
!!! This module extends the Variable-Blob Particles (ppm_module_vbp_typdef)

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2


USE ppm_module_interfaces
USE ppm_module_topo_typedef
USE ppm_module_particles_typedef
USE ppm_module_vbp_typedef
USE ppm_module_operator_typedef

IMPLICIT NONE

PRIVATE

!----------------------------------------------------------------------
! Global variables 
!----------------------------------------------------------------------

#define  DTYPE(a) a/**/_s
#define  CTYPE(a) a/**/_sc
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#include "sop/sop_typedef.f"

#define  DTYPE(a) a/**/_d
#define  CTYPE(a) a/**/_dc
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "sop/sop_typedef.f"

INTEGER, PRIVATE, DIMENSION(3)    :: ldc
!!! Number of elements in all dimensions for allocation

PUBLIC :: ppm_t_sop_s, ppm_t_sop_d

INTERFACE sop_init_opts
    MODULE PROCEDURE sop_init_opts_s
    MODULE PROCEDURE sop_init_opts_d
END INTERFACE

CONTAINS

#define DTYPE(a) a/**/_s
#define __KIND __SINGLE_PRECISION
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#include "sop/sop_typeproc.f"
#include "sop/ppm_sop_helpers.f"
#include "sop/sop_adapt_particles.f"
#undef  DEFINE_MK

#undef  DTYPE
#undef  __KIND


#define DTYPE(a) a/**/_d
#define __KIND __DOUBLE_PRECISION
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#include "sop/sop_typeproc.f"
#include "sop/ppm_sop_helpers.f"
#include "sop/sop_adapt_particles.f"
#undef  DEFINE_MK

#undef  DTYPE
#undef  __KIND

END MODULE ppm_module_sop


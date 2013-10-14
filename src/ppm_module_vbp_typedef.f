      MODULE ppm_module_vbp_typedef
      !!! Variable-blob particles
      !!! (basically particles for which the cutoff radius can be different
      !!! for each particle)
      !!! The type extends the vanilla particle data type

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      USE ppm_module_interfaces
      USE ppm_module_topo_typedef
      USE ppm_module_particles_typedef
      
      IMPLICIT NONE
      
      PRIVATE
      
      !----------------------------------------------------------------------
      ! Global variables
      !----------------------------------------------------------------------
#define  DTYPE(a) a/**/_s
#define  CTYPE(a) a/**/_sc
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#include "vbp/vbp_typedef.f"

#define  DTYPE(a) a/**/_d
#define  CTYPE(a) a/**/_dc
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "vbp/vbp_typedef.f"

      INTEGER, PRIVATE, DIMENSION(3) :: ldc
      !!! Number of elements in all dimensions for allocation
      
      PUBLIC :: ppm_t_vbp_s, ppm_t_vbp_d
      
      CONTAINS

#define DTYPE(a) a/**/_s
#define __KIND __SINGLE_PRECISION
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#include "vbp/vbp_typeproc.f"
#undef  DEFINE_MK

#undef  DTYPE
#undef  __KIND


#define DTYPE(a) a/**/_d
#define __KIND __DOUBLE_PRECISION
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#include "vbp/vbp_typeproc.f"
#undef  DEFINE_MK

#undef  DTYPE
#undef  __KIND

      END MODULE ppm_module_vbp_typedef


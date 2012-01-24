MODULE ppm_module_particles_typedef

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2


USE ppm_module_alloc
USE ppm_module_typedef
USE ppm_module_data, ONLY: ppm_rank,ppm_dim
USE ppm_module_error
USE ppm_module_write
USE ppm_module_substart
USE ppm_module_substop

IMPLICIT NONE

PRIVATE

!PPM internal parameters used only to access entries in the
!particle data structure.
INTEGER,PARAMETER   :: ppm_part_ghosts = 1
INTEGER,PARAMETER   :: ppm_part_partial = 2
INTEGER,PARAMETER   :: ppm_part_reqput = 3
INTEGER,PARAMETER   :: ppm_part_areinside = 4
INTEGER,PARAMETER   :: ppm_part_cartesian = 5
INTEGER,PARAMETER   :: ppm_param_length_partflags = 5

!PPM internal parameters used only to access entries in the
!particle's property data structure.
INTEGER,PARAMETER   :: ppm_ppt_ghosts = 1
INTEGER,PARAMETER   :: ppm_ppt_partial = 2
INTEGER,PARAMETER   :: ppm_ppt_reqput = 3
INTEGER,PARAMETER   :: ppm_ppt_map_parts = 4
INTEGER,PARAMETER   :: ppm_ppt_map_ghosts = 5
INTEGER,PARAMETER   :: ppm_param_length_pptflags = 5


#define  DTYPE(a) a/**/_s
#define  CTYPE(a) a/**/_sc
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#include "part/particles_typedef.inc"

#define  DTYPE(a) a/**/_d
#define  CTYPE(a) a/**/_dc
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "part/particles_typedef.inc"

CHARACTER(LEN=ppm_char)         :: cbuf
CHARACTER(LEN=ppm_char)         :: line_of_stars='**************************'
INTEGER, PRIVATE, DIMENSION(3)    :: ldc
!!! Number of elements in all dimensions for allocation

PUBLIC :: ppm_t_particles_s, ppm_t_particles_d

CONTAINS

#include "part/ppm_particles_helpers.f"

#define DTYPE(a) a/**/_s
#define __KIND __SINGLE_PRECISION
#define __DIM 1
#define __TYPE INTEGER
#define DATANAME data_1d_i
#include "part/ppm_particles_get.f"
#define __TYPE INTEGER(ppm_kind_int64)
#define DATANAME data_1d_li
#include "part/ppm_particles_get.f"
#define __TYPE REAL(ppm_kind_single)
#define DATANAME data_1d_r
#include "part/ppm_particles_get.f"
#define __TYPE COMPLEX(ppm_kind_single)
#define DATANAME data_1d_c
#include "part/ppm_particles_get.f"
#define __TYPE LOGICAL
#define DATANAME data_1d_l
#include "part/ppm_particles_get.f"
#undef  __DIM

#define __DIM 2
#define __TYPE INTEGER
#define DATANAME data_2d_i
#include "part/ppm_particles_get.f"
#define __TYPE INTEGER(ppm_kind_int64)
#define DATANAME data_2d_li
#include "part/ppm_particles_get.f"
#define __TYPE REAL(ppm_kind_single)
#define DATANAME data_2d_r
#include "part/ppm_particles_get.f"
#define __TYPE COMPLEX(ppm_kind_single)
#define DATANAME data_2d_c
#include "part/ppm_particles_get.f"
#define __TYPE LOGICAL
#define DATANAME data_2d_l
#include "part/ppm_particles_get.f"
#undef  __DIM

#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#include "part/ppm_part_OOprocs.f"

#undef  DTYPE
#undef  __KIND


#define DTYPE(a) a/**/_d
#define __KIND __DOUBLE_PRECISION
#define __DIM 1
#define __TYPE INTEGER
#define DATANAME data_1d_i
#include "part/ppm_particles_get.f"
#define __TYPE INTEGER(ppm_kind_int64)
#define DATANAME data_1d_li
#include "part/ppm_particles_get.f"
#define __TYPE REAL(ppm_kind_double)
#define DATANAME data_1d_r
#include "part/ppm_particles_get.f"
#define __TYPE COMPLEX(ppm_kind_double)
#define DATANAME data_1d_c
#include "part/ppm_particles_get.f"
#define __TYPE LOGICAL
#define DATANAME data_1d_l
#include "part/ppm_particles_get.f"
#undef  __DIM

#define __DIM 2
#define __TYPE INTEGER
#define DATANAME data_2d_i
#include "part/ppm_particles_get.f"
#define __TYPE INTEGER(ppm_kind_int64)
#define DATANAME data_2d_li
#include "part/ppm_particles_get.f"
#define __TYPE REAL(ppm_kind_double)
#define DATANAME data_2d_r
#include "part/ppm_particles_get.f"
#define __TYPE COMPLEX(ppm_kind_double)
#define DATANAME data_2d_c
#include "part/ppm_particles_get.f"
#define __TYPE LOGICAL
#define DATANAME data_2d_l
#include "part/ppm_particles_get.f"
#undef  __DIM

#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#include "part/ppm_part_OOprocs.f"

#undef  DTYPE
#undef  __KIND

END MODULE ppm_module_particles_typedef


MODULE ppm_module_particles_typedef

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

#define __REAL_SINGLE 3 
#define __REAL_DOUBLE 4 
#define __COMPLEX_SINGLE 5 
#define __COMPLEX_DOUBLE 6 
#define __INTEGER 7 
#define __LONGINT 8 
#define __LOGICAL 9 
#define __CHAR 10 

#define __crash_on_null_pointers  1

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
INTEGER,PARAMETER   :: ppm_part_neighlists = 6
INTEGER,PARAMETER   :: ppm_param_length_partflags = 6

!PPM internal parameters used only to access entries in the
!particle's property data structure.
INTEGER,PARAMETER   :: ppm_ppt_ghosts = 1
INTEGER,PARAMETER   :: ppm_ppt_partial = 2
INTEGER,PARAMETER   :: ppm_ppt_reqput = 3
INTEGER,PARAMETER   :: ppm_ppt_map_parts = 4
INTEGER,PARAMETER   :: ppm_ppt_map_ghosts = 5
INTEGER,PARAMETER   :: ppm_param_length_pptflags = 5

!PPM internal parameters used only to access entries in the
!particle's property data structure.
INTEGER,PARAMETER   :: ppm_ops_inc_ghosts = 1
INTEGER,PARAMETER   :: ppm_ops_interp = 2
INTEGER,PARAMETER   :: ppm_ops_iscomputed = 3
INTEGER,PARAMETER   :: ppm_ops_isdefined = 4
INTEGER,PARAMETER   :: ppm_ops_vector = 5
INTEGER,PARAMETER   :: ppm_param_length_opsflags = 5


!User parameters 
INTEGER, PARAMETER :: ppm_param_part_init_cartesian = 1
INTEGER, PARAMETER :: ppm_param_part_init_random = 2

!----------------------------------------------------------------------
! Global variables 
!----------------------------------------------------------------------

INTEGER                               :: ppm_particles_seedsize
INTEGER,  DIMENSION(:  ), POINTER     :: ppm_particles_seed => NULL()

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
PUBLIC :: ppm_t_sop_s, ppm_t_sop_d


CONTAINS


#include "part/ppm_particles_helpers.f"

#define DTYPE(a) a/**/_s
#define __KIND __SINGLE_PRECISION
#define __DIM 1
#define __TYPE INTEGER
#define __MYTYPE __INTEGER
#define DATANAME data_1d_i
#include "part/ppm_particles_get.f"
#define __TYPE INTEGER(ppm_kind_int64)
#define __MYTYPE __LONGINT
#define DATANAME data_1d_li
#include "part/ppm_particles_get.f"
#define __TYPE REAL(ppm_kind_single)
#define __MYTYPE __REAL_SINGLE
#define DATANAME data_1d_r
#include "part/ppm_particles_get.f"
#define __TYPE COMPLEX(ppm_kind_single)
#define __MYTYPE __COMPLEX_SINGLE
#define DATANAME data_1d_c
#include "part/ppm_particles_get.f"
#define __TYPE LOGICAL
#define __MYTYPE __LOGICAL
#define DATANAME data_1d_l
#include "part/ppm_particles_get.f"
#undef  __DIM

#define __DIM 2
#define __TYPE INTEGER
#define __MYTYPE __INTEGER
#define DATANAME data_2d_i
#include "part/ppm_particles_get.f"
#define __TYPE INTEGER(ppm_kind_int64)
#define __MYTYPE __LONGINT
#define DATANAME data_2d_li
#include "part/ppm_particles_get.f"
#define __TYPE REAL(ppm_kind_single)
#define __MYTYPE __REAL_SINGLE
#define DATANAME data_2d_r
#include "part/ppm_particles_get.f"
#define __TYPE COMPLEX(ppm_kind_single)
#define __MYTYPE __COMPLEX_SINGLE
#define DATANAME data_2d_c
#include "part/ppm_particles_get.f"
#define __TYPE LOGICAL
#define __MYTYPE __LOGICAL
#define DATANAME data_2d_l
#include "part/ppm_particles_get.f"
#undef  __DIM

#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#include "part/ppm_part_OOprocs.f"

#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#include "part/ppm_dcop_helpers.f"
#define __DIM 2
#include "part/ppm_dcop_compute.f"
#define __DIM 3
#include "part/ppm_dcop_compute.f"
#undef DEFINE_MK

#include "part/container_procedures.inc"

#undef  DTYPE
#undef  __KIND


#define DTYPE(a) a/**/_d
#define __KIND __DOUBLE_PRECISION
#define __DIM 1
#define __TYPE INTEGER
#define __MYTYPE __INTEGER
#define DATANAME data_1d_i
#include "part/ppm_particles_get.f"
#define __TYPE INTEGER(ppm_kind_int64)
#define __MYTYPE __LONGINT
#define DATANAME data_1d_li
#include "part/ppm_particles_get.f"
#define __TYPE REAL(ppm_kind_double)
#define __MYTYPE __REAL_DOUBLE
#define DATANAME data_1d_r
#include "part/ppm_particles_get.f"
#define __TYPE COMPLEX(ppm_kind_double)
#define __MYTYPE __COMPLEX_DOUBLE
#define DATANAME data_1d_c
#include "part/ppm_particles_get.f"
#define __TYPE LOGICAL
#define __MYTYPE __LOGICAL
#define DATANAME data_1d_l
#include "part/ppm_particles_get.f"
#undef  __DIM

#define __DIM 2
#define __TYPE INTEGER
#define __MYTYPE __INTEGER
#define DATANAME data_2d_i
#include "part/ppm_particles_get.f"
#define __TYPE INTEGER(ppm_kind_int64)
#define __MYTYPE __LONGINT
#define DATANAME data_2d_li
#include "part/ppm_particles_get.f"
#define __TYPE REAL(ppm_kind_double)
#define __MYTYPE __REAL_DOUBLE
#define DATANAME data_2d_r
#include "part/ppm_particles_get.f"
#define __TYPE COMPLEX(ppm_kind_double)
#define __MYTYPE __COMPLEX_DOUBLE
#define DATANAME data_2d_c
#include "part/ppm_particles_get.f"
#define __TYPE LOGICAL
#define __MYTYPE __LOGICAL
#define DATANAME data_2d_l
#include "part/ppm_particles_get.f"
#undef  __DIM

#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#include "part/ppm_part_OOprocs.f"

#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#include "part/ppm_dcop_helpers.f"
#define __DIM 2
#include "part/ppm_dcop_compute.f"
#define __DIM 3
#include "part/ppm_dcop_compute.f"
#undef DEFINE_MK

#include "part/container_procedures.inc"

#undef  DTYPE
#undef  __KIND

END MODULE ppm_module_particles_typedef


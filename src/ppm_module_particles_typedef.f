MODULE ppm_module_particles_typedef

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

#define __REAL 3 
#define __COMPLEX 4 
#define __INTEGER 5 
#define __LONGINT 6 
#define __LOGICAL 7 
#define __CHAR 8 

USE ppm_module_alloc
USE ppm_module_interfaces
USE ppm_module_topo_typedef
USE ppm_module_mapping_typedef
USE ppm_module_data
USE ppm_module_error
USE ppm_module_write
USE ppm_module_substart
USE ppm_module_substop

IMPLICIT NONE

!----------------------------------------------------------------------
! Global variables 
!----------------------------------------------------------------------

INTEGER                               :: ppm_particles_seedsize
INTEGER,  DIMENSION(:  ), POINTER     :: ppm_particles_seed => NULL()

#define  DTYPE(a) a/**/_s
#define  CTYPE(a) a/**/_sc
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#include "part/particles_typedef.f"

#define  DTYPE(a) a/**/_d
#define  CTYPE(a) a/**/_dc
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "part/particles_typedef.f"

CHARACTER(LEN=ppm_char)         :: cbuf
INTEGER, PRIVATE, DIMENSION(3)  :: ldc
!!! Number of elements in all dimensions for allocation

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
#define __MYTYPE __REAL
#define DATANAME data_1d_r
#include "part/ppm_particles_get.f"
#define __TYPE COMPLEX(ppm_kind_single)
#define __MYTYPE __COMPLEX
#define DATANAME data_1d_c
#include "part/ppm_particles_get.f"
#define __TYPE LOGICAL
#define __MYTYPE __LOGICAL
#define DATANAME data_1d_l
#include "part/ppm_particles_get.f"
#undef  __DIM

#define __DIM 2
#define __TYPE INTEGER
#define __MYTYPE __LONGINT
#define DATANAME data_2d_i
#include "part/ppm_particles_get.f"
#define __TYPE INTEGER(ppm_kind_int64)
#define __MYTYPE __LONGINT
#define DATANAME data_2d_li
#include "part/ppm_particles_get.f"
#define __TYPE REAL(ppm_kind_single)
#define __MYTYPE __REAL
#define DATANAME data_2d_r
#include "part/ppm_particles_get.f"
#define __TYPE COMPLEX(ppm_kind_single)
#define __MYTYPE __COMPLEX
#define DATANAME data_2d_c
#include "part/ppm_particles_get.f"
#define __TYPE LOGICAL
#define __MYTYPE __LOGICAL
#define DATANAME data_2d_l
#include "part/ppm_particles_get.f"
#undef  __DIM

#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#include "part/particles_typeproc.f"
#include "part/ppm_part_neighlists_get.f"
#include "part/part_interp_to_mesh.f"
#undef  DEFINE_MK
# define DEFINE_MK() parameter(mk,<#integer#>,ppm_kind_single)
#include "part/part_remesh.f"
#undef  DEFINE_MK

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
#define __MYTYPE __REAL
#define DATANAME data_1d_r
#include "part/ppm_particles_get.f"
#define __TYPE COMPLEX(ppm_kind_double)
#define __MYTYPE __COMPLEX
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
#define __MYTYPE __REAL
#define DATANAME data_2d_r
#include "part/ppm_particles_get.f"
#define __TYPE COMPLEX(ppm_kind_double)
#define __MYTYPE __COMPLEX
#define DATANAME data_2d_c
#include "part/ppm_particles_get.f"
#define __TYPE LOGICAL
#define __MYTYPE __LOGICAL
#define DATANAME data_2d_l
#include "part/ppm_particles_get.f"
#undef  __DIM

#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#include "part/particles_typeproc.f"
#include "part/ppm_part_neighlists_get.f"
#include "part/part_interp_to_mesh.f"
#undef  DEFINE_MK
# define DEFINE_MK() parameter(mk,<#integer#>,ppm_kind_double)
#include "part/part_remesh.f"
#undef  DEFINE_MK

#undef  DTYPE
#undef  __KIND

END MODULE ppm_module_particles_typedef


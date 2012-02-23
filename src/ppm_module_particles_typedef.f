      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_particles_typdef
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License 
      ! as published by the Free Software Foundation, either 
      ! version 3 of the License, or (at your option) any later 
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

#define __REAL 3 
#define __COMPLEX 4 
#define __INTEGER 5 
#define __LONGINT 6 
#define __LOGICAL 7 
#define __CHAR 8 

#define __crash_on_null_pointers  1
#undef __WITH_KDTREE

      MODULE ppm_module_particles_typedef
      !!! Declares particle set data types
      !!!
      !!! [NOTE]
      !!! Most of the declared variables in this module should not be accessed
      !!! directly by the PPM client developer, they are used internally in the
      !!! library.

         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_alloc
         USE ppm_module_data
         USE ppm_module_container_typedef
         USE ppm_module_mapping_typedef
         USE ppm_module_data, ONLY: ppm_rank,ppm_dim,ppm_comm
         USE ppm_module_error
         USE ppm_module_write
         USE ppm_module_substart
         USE ppm_module_substop

         IMPLICIT NONE
         PRIVATE
         !----------------------------------------------------------------------
         ! Global parameters
         !----------------------------------------------------------------------

         !User parameters 
         INTEGER,PARAMETER,PUBLIC   :: ppm_param_part_init_cartesian = 1
         INTEGER,PARAMETER,PUBLIC   :: ppm_param_part_init_random = 2


         !PPM internal parameters used only to access entries in the
         !particle data structure.
         INTEGER,PARAMETER,PUBLIC   :: ppm_part_ghosts = 1
         INTEGER,PARAMETER,PUBLIC   :: ppm_part_partial = 2
         INTEGER,PARAMETER,PUBLIC   :: ppm_part_reqput = 3
         INTEGER,PARAMETER,PUBLIC   :: ppm_part_areinside = 4
         INTEGER,PARAMETER,PUBLIC   :: ppm_part_cartesian = 5
         INTEGER,PARAMETER,PUBLIC   :: ppm_part_neighlists = 6
         INTEGER,PARAMETER,PUBLIC   :: ppm_part_global_index = 7
         INTEGER,PARAMETER,PUBLIC   :: ppm_param_length_partflags = 7

         !PPM internal parameters used only to access entries in the
         !particle's property data structure.
         INTEGER,PARAMETER,PUBLIC   :: ppm_ppt_ghosts = 1
         INTEGER,PARAMETER,PUBLIC   :: ppm_ppt_partial = 2
         INTEGER,PARAMETER,PUBLIC   :: ppm_ppt_reqput = 3
         INTEGER,PARAMETER,PUBLIC   :: ppm_ppt_map_parts = 4
         INTEGER,PARAMETER,PUBLIC   :: ppm_ppt_map_ghosts = 5
         INTEGER,PARAMETER,PUBLIC   :: ppm_param_length_pptflags = 5

         !PPM internal parameters used only to access entries in the
         !particle's property data structure.
         INTEGER,PARAMETER,PUBLIC   :: ppm_ops_inc_ghosts = 1
         INTEGER,PARAMETER,PUBLIC   :: ppm_ops_interp = 2
         INTEGER,PARAMETER,PUBLIC   :: ppm_ops_iscomputed = 3
         INTEGER,PARAMETER,PUBLIC   :: ppm_ops_isdefined = 4
         INTEGER,PARAMETER,PUBLIC   :: ppm_ops_vector = 5
         INTEGER,PARAMETER,PUBLIC   :: ppm_param_length_opsflags = 5

         !PPM internal parameters for default storage IDs of some DS.
         INTEGER, PARAMETER,PUBLIC :: ppm_param_default_nlID = 1


         !----------------------------------------------------------------------
         ! Global variables 
         !----------------------------------------------------------------------

         INTEGER                               :: ppm_particles_seedsize
         INTEGER,  DIMENSION(:  ), POINTER     :: ppm_particles_seed => NULL()

         !----------------------------------------------------------------------
         ! Type declaration
         !----------------------------------------------------------------------

#define  DTYPE(a) a/**/_s
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#include "part/particles_typedef.inc"

#define  DTYPE(a) a/**/_d
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "part/particles_typedef.inc"

         CHARACTER(LEN=ppm_char)         :: cbuf
         CHARACTER(LEN=ppm_char)         :: line_of_stars='********************'
         INTEGER, PRIVATE, DIMENSION(3)  :: ldc
         !!! Number of elements in all dimensions for allocation

         PUBLIC :: ppm_t_particles_s, ppm_t_particles_d
         PUBLIC :: ppm_t_sop_s, ppm_t_sop_d
         PUBLIC :: ppm_t_part_prop_s, ppm_t_part_prop_d

         !----------------------------------------------------------------------
         ! Type-bound procedures
         !----------------------------------------------------------------------


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
#define __MYTYPE __INTEGER
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

#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#include "part/ppm_dcop_helpers.f"
#define __DIM 2
#include "part/ppm_dcop_compute.f"
#define __DIM 3
#include "part/ppm_dcop_compute.f"
#undef DEFINE_MK

#include "part/ppm_part_neighlists_get.f"

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

#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#include "part/ppm_dcop_helpers.f"
#define __DIM 2
#include "part/ppm_dcop_compute.f"
#define __DIM 3
#include "part/ppm_dcop_compute.f"
#undef DEFINE_MK

#include "part/ppm_part_neighlists_get.f"

#undef  DTYPE
#undef  __KIND


#define DTYPE(a) a/**/_s
#define DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single

#define __DIM 1
#define __KIND __SINGLE_PRECISION
#include "part/map_part_pop.f"
#define __KIND __INTEGER
#include "part/map_part_pop.f"
#define __KIND __LOGICAL
#include "part/map_part_pop.f"
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "part/map_part_pop.f"
#undef __DIM

#define __DIM 2
#define __KIND __SINGLE_PRECISION
#include "part/map_part_pop.f"
#define __KIND __INTEGER
#include "part/map_part_pop.f"
#define __KIND __LOGICAL
#include "part/map_part_pop.f"
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "part/map_part_pop.f"
#undef __DIM

#define __DIM 1
#define __KIND __SINGLE_PRECISION
#include "part/map_part_push.f"
#define __KIND __INTEGER
#include "part/map_part_push.f"
#define __KIND __LOGICAL
#include "part/map_part_push.f"
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "part/map_part_push.f"
#undef __DIM

#define __DIM 2
#define __KIND __SINGLE_PRECISION
#include "part/map_part_push.f"
#define __KIND __INTEGER
#include "part/map_part_push.f"
#define __KIND __LOGICAL
#include "part/map_part_push.f"
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "part/map_part_push.f"
#undef __DIM

#define __KIND __SINGLE_PRECISION
#include "part/map_part_send.f"

#undef  DTYPE
#undef  DEFINE_MK

#define DTYPE(a) a/**/_d
#define DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double

#define __DIM 1
#define __KIND __DOUBLE_PRECISION
#include "part/map_part_pop.f"
#define __KIND __INTEGER
#include "part/map_part_pop.f"
#define __KIND __LOGICAL
#include "part/map_part_pop.f"
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "part/map_part_pop.f"
#undef __DIM

#define __DIM 2
#define __KIND __DOUBLE_PRECISION
#include "part/map_part_pop.f"
#define __KIND __INTEGER
#include "part/map_part_pop.f"
#define __KIND __LOGICAL
#include "part/map_part_pop.f"
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "part/map_part_pop.f"
#undef __DIM

#define __DIM 1
#define __KIND __DOUBLE_PRECISION
#include "part/map_part_push.f"
#define __KIND __INTEGER
#include "part/map_part_push.f"
#define __KIND __LOGICAL
#include "part/map_part_push.f"
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "part/map_part_push.f"
#undef __DIM

#define __DIM 2
#define __KIND __DOUBLE_PRECISION
#include "part/map_part_push.f"
#define __KIND __INTEGER
#include "part/map_part_push.f"
#define __KIND __LOGICAL
#include "part/map_part_push.f"
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "part/map_part_push.f"
#undef __DIM

#define __KIND __DOUBLE_PRECISION
#include "part/map_part_send.f"

#undef  DTYPE
#undef  DEFINE_MK

         END MODULE ppm_module_particles_typedef


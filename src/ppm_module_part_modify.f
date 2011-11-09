      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_part_modify
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
     
      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION         1
#define __DOUBLE_PRECISION         2
#define __INTEGER                  3
#define __LOGICAL                  4
#define __SINGLE_PRECISION_COMPLEX 5
#define __DOUBLE_PRECISION_COMPLEX 6

      MODULE ppm_module_part_modify
      !!! This module provides the basic mapping routines for particles, namely
      !!! the push, sned and pop routine.
      !!!
      !!! The module also holds the needed work arrays for these routines.
         
         USE ppm_module_data, ONLY: ppm_kind_single,ppm_kind_double
         PRIVATE :: ppm_kind_single,ppm_kind_double

         !----------------------------------------------------------------------
         !  types
         !----------------------------------------------------------------------

         TYPE ppm_t_part_modify
             INTEGER                            :: Nrnew
             !number of new real particles 
             INTEGER                            :: Ngnew
             !number of new ghost particles 
             INTEGER                            :: Ngsendnew
             !number of new real particles that are ghosts for other procs
             INTEGER, DIMENSION(:), POINTER     :: idx_real_new => NULL()
             !indeces of new real particles 
             INTEGER, DIMENSION(:), POINTER     :: idx_ghost_new => NULL()
             !indeces of new ghost particles 
             INTEGER, DIMENSION(:), POINTER     :: idx_ghost_send_new => NULL()
             !indeces of new real particles that are ghosts for other procs
             INTEGER, DIMENSION(:), POINTER     :: psend_new_offset => NULL()
         END TYPE

         !----------------------------------------------------------------------
         !  data
         !----------------------------------------------------------------------

         TYPE(ppm_t_part_modify)                :: modify

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __DIM 1
#define __KIND __SINGLE_PRECISION
#include "part/ppm_part_modify_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "part/ppm_part_modify_pop.f"
#undef __KIND
#define __KIND __INTEGER
#include "part/ppm_part_modify_pop.f"
#undef __KIND
#define __KIND __LOGICAL
#include "part/ppm_part_modify_pop.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "part/ppm_part_modify_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "part/ppm_part_modify_pop.f"
#undef __KIND
#undef __DIM

#define __DIM 2
#define __KIND __SINGLE_PRECISION
#include "part/ppm_part_modify_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "part/ppm_part_modify_pop.f"
#undef __KIND
#define __KIND __INTEGER
#include "part/ppm_part_modify_pop.f"
#undef __KIND
#define __KIND __LOGICAL
#include "part/ppm_part_modify_pop.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "part/ppm_part_modify_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "part/ppm_part_modify_pop.f"
#undef __KIND
#undef __DIM

#define __DIM 1
#define __KIND __SINGLE_PRECISION
#include "part/ppm_part_modify_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "part/ppm_part_modify_push.f"
#undef __KIND
#define __KIND __INTEGER
#include "part/ppm_part_modify_push.f"
#undef __KIND
#define __KIND __LOGICAL
#include "part/ppm_part_modify_push.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "part/ppm_part_modify_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "part/ppm_part_modify_push.f"
#undef __KIND
#undef __DIM

#define __DIM 2
#define __KIND __SINGLE_PRECISION
#include "part/ppm_part_modify_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "part/ppm_part_modify_push.f"
#undef __KIND
#define __KIND __INTEGER
#include "part/ppm_part_modify_push.f"
#undef __KIND
#define __KIND __LOGICAL
#include "part/ppm_part_modify_push.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "part/ppm_part_modify_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "part/ppm_part_modify_push.f"
#undef __KIND
#undef __DIM

#include "part/ppm_part_modify_send.f"

#define __KIND __SINGLE_PRECISION
#include "part/ppm_part_split_compute.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "part/ppm_part_split_compute.f"

#define __DIM 1
#define __KIND __SINGLE_PRECISION
#include "part/ppm_part_split_apply.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "part/ppm_part_split_apply.f"
#undef __KIND
#define __KIND __INTEGER
#include "part/ppm_part_split_apply.f"
#undef __KIND
#define __KIND __LOGICAL
#include "part/ppm_part_split_apply.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "part/ppm_part_split_apply.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "part/ppm_part_split_apply.f"
#undef __KIND
#undef __DIM

#define __DIM 2
#define __KIND __SINGLE_PRECISION
#include "part/ppm_part_split_apply.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "part/ppm_part_split_apply.f"
#undef __KIND
#define __KIND __INTEGER
#include "part/ppm_part_split_apply.f"
#undef __KIND
#define __KIND __LOGICAL
#include "part/ppm_part_split_apply.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "part/ppm_part_split_apply.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "part/ppm_part_split_apply.f"
#undef __KIND
#undef __DIM

      END MODULE ppm_module_part_modify

      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_map_part_pop
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
#define __NORMAL                   7
#define __ADD                      8

      MODULE ppm_module_map_part
      !!! This module provides the basic mapping routines for particles, namely
      !!! the push, sned and pop routine.
      !!!
      !!! The module also holds the needed work arrays for these routines.
         
         USE ppm_module_map_part_util
         USE ppm_module_map_part_ghost
         USE ppm_module_map_part_global
         USE ppm_module_map_part_partial
         USE ppm_module_data, ONLY: ppm_kind_single,ppm_kind_double
         PRIVATE :: ppm_kind_single,ppm_kind_double

         !----------------------------------------------------------------------
         !  Work lists
         !----------------------------------------------------------------------
         REAL(ppm_kind_single), DIMENSION(:), POINTER :: sends => NULL()
         REAL(ppm_kind_single), DIMENSION(:), POINTER :: recvs => NULL()
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: sendd => NULL()
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: recvd => NULL()

         TYPE ppm_t_plists
             INTEGER, DIMENSION(:), POINTER   :: nsend => NULL()
             INTEGER, DIMENSION(:), POINTER   :: nrecv => NULL()
             INTEGER, DIMENSION(:), POINTER   :: psend => NULL()
             INTEGER, DIMENSION(:), POINTER   :: precv => NULL()
         END TYPE

         TYPE(ppm_t_plists), TARGET, SAVE    :: plists_normal
         TYPE(ppm_t_plists), TARGET, SAVE    :: plists_add

         INTEGER, DIMENSION(:,:), POINTER :: pp    => NULL()
         INTEGER, DIMENSION(:,:), POINTER :: qq    => NULL()

         !PRIVATE :: sends,recvs,sendd,recvd,nsend,nrecv,psend,precv,pp,qq
         PRIVATE :: sends,recvs,sendd,recvd,pp,qq,plists_normal,plists_add
         !PRIVATE :: old_nsendlist,old_buffer_set

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_pop
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_pop
            ! scalar (1d) particle data
            MODULE PROCEDURE ppm_map_part_pop_1dd
            MODULE PROCEDURE ppm_map_part_pop_1ds
            MODULE PROCEDURE ppm_map_part_pop_1di
            MODULE PROCEDURE ppm_map_part_pop_1dl
            MODULE PROCEDURE ppm_map_part_pop_1ddc
            MODULE PROCEDURE ppm_map_part_pop_1dsc

            ! vector (2d) particle data
            MODULE PROCEDURE ppm_map_part_pop_2dd
            MODULE PROCEDURE ppm_map_part_pop_2ds
            MODULE PROCEDURE ppm_map_part_pop_2di
            MODULE PROCEDURE ppm_map_part_pop_2dl
            MODULE PROCEDURE ppm_map_part_pop_2ddc
            MODULE PROCEDURE ppm_map_part_pop_2dsc
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_push
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_push
            ! scalar (1d) particle data
            MODULE PROCEDURE ppm_map_part_push_1dd
            MODULE PROCEDURE ppm_map_part_push_1ds
            MODULE PROCEDURE ppm_map_part_push_1di
            MODULE PROCEDURE ppm_map_part_push_1dl
            MODULE PROCEDURE ppm_map_part_push_1ddc
            MODULE PROCEDURE ppm_map_part_push_1dsc

            ! vector (2d) particle data
            MODULE PROCEDURE ppm_map_part_push_2dd
            MODULE PROCEDURE ppm_map_part_push_2ds
            MODULE PROCEDURE ppm_map_part_push_2di
            MODULE PROCEDURE ppm_map_part_push_2dl
            MODULE PROCEDURE ppm_map_part_push_2ddc
            MODULE PROCEDURE ppm_map_part_push_2dsc
         END INTERFACE
         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_send
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_send
            MODULE PROCEDURE ppm_map_part_send
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_pop_add
         !----------------------------------------------------------------------
         INTERFACE ppm_part_modify_pop
            ! scalar (1d) particle data
            MODULE PROCEDURE ppm_map_part_pop_add1dd
            MODULE PROCEDURE ppm_map_part_pop_add1ds
            MODULE PROCEDURE ppm_map_part_pop_add1di
            MODULE PROCEDURE ppm_map_part_pop_add1dl
            MODULE PROCEDURE ppm_map_part_pop_add1ddc
            MODULE PROCEDURE ppm_map_part_pop_add1dsc

            ! vector (2d) particle data
            MODULE PROCEDURE ppm_map_part_pop_add2dd
            MODULE PROCEDURE ppm_map_part_pop_add2ds
            MODULE PROCEDURE ppm_map_part_pop_add2di
            MODULE PROCEDURE ppm_map_part_pop_add2dl
            MODULE PROCEDURE ppm_map_part_pop_add2ddc
            MODULE PROCEDURE ppm_map_part_pop_add2dsc
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_push
         !----------------------------------------------------------------------
         INTERFACE ppm_part_modify_push
            ! scalar (1d) particle data
            MODULE PROCEDURE ppm_map_part_push_add1dd
            MODULE PROCEDURE ppm_map_part_push_add1ds
            MODULE PROCEDURE ppm_map_part_push_add1di
            MODULE PROCEDURE ppm_map_part_push_add1dl
            MODULE PROCEDURE ppm_map_part_push_add1ddc
            MODULE PROCEDURE ppm_map_part_push_add1dsc

            ! vector (2d) particle data
            MODULE PROCEDURE ppm_map_part_push_add2dd
            MODULE PROCEDURE ppm_map_part_push_add2ds
            MODULE PROCEDURE ppm_map_part_push_add2di
            MODULE PROCEDURE ppm_map_part_push_add2dl
            MODULE PROCEDURE ppm_map_part_push_add2ddc
            MODULE PROCEDURE ppm_map_part_push_add2dsc
         END INTERFACE
         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_send_add
         !----------------------------------------------------------------------
         INTERFACE ppm_part_modify_send
            MODULE PROCEDURE ppm_map_part_send_add
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __VARIANT __NORMAL
#define __DIM 1
#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __INTEGER
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __LOGICAL
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "map/ppm_map_part_pop.f"
#undef __KIND
#undef __DIM

#define __DIM 2
#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __INTEGER
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __LOGICAL
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "map/ppm_map_part_pop.f"
#undef __KIND
#undef __DIM

#define __DIM 1
#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __INTEGER
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __LOGICAL
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "map/ppm_map_part_push.f"
#undef __KIND
#undef __DIM

#define __DIM 2
#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __INTEGER
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __LOGICAL
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "map/ppm_map_part_push.f"
#undef __KIND
#undef __DIM

#include "map/ppm_map_part_send.f"

#undef __VARIANT


#define __VARIANT __ADD
#define __DIM 1
#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __INTEGER
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __LOGICAL
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "map/ppm_map_part_pop.f"
#undef __KIND
#undef __DIM

#define __DIM 2
#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __INTEGER
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __LOGICAL
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "map/ppm_map_part_pop.f"
#undef __KIND
#undef __DIM


#define __DIM 1
#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __INTEGER
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __LOGICAL
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "map/ppm_map_part_push.f"
#undef __KIND
#undef __DIM

#define __DIM 2
#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __INTEGER
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __LOGICAL
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "map/ppm_map_part_push.f"
#undef __KIND
#undef __DIM


#include "map/ppm_map_part_send.f"

#undef __VARIANT

      END MODULE ppm_module_map_part

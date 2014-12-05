      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_map_part_util
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
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

      MODULE ppm_module_map_part_util
      !!! This module provides various particle mapping utility and helper
      !!! routines and their work-arrays.

         IMPLICIT NONE

         PRIVATE

         !----------------------------------------------------------------------
         !  Work lists
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:), POINTER :: plist_des => NULL()
         INTEGER, DIMENSION(:), POINTER :: plist_act => NULL()
         INTEGER, DIMENSION(:), POINTER :: plist_exc => NULL()
         INTEGER, DIMENSION(:), POINTER :: srlist1   => NULL()
         INTEGER, DIMENSION(:), POINTER :: srlist2   => NULL()
         INTEGER                        :: nlist1,nlist2,nlist3,nlist4
         INTEGER ,DIMENSION(:), POINTER :: ilist1    => NULL()
         INTEGER ,DIMENSION(:), POINTER :: ilist2    => NULL()
         INTEGER ,DIMENSION(:), POINTER :: ilist3    => NULL()
         INTEGER ,DIMENSION(:), POINTER :: ilist4    => NULL()

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_eqdistrib
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_eqdistrib
            MODULE PROCEDURE ppm_map_part_eqdistrib_d
            MODULE PROCEDURE ppm_map_part_eqdistrib_s
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_send
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_cancel
            MODULE PROCEDURE ppm_map_part_cancel
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to the ppm_map_part_get_sub
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_get_sub
            MODULE PROCEDURE ppm_map_part_get_sub_s
            MODULE PROCEDURE ppm_map_part_get_sub_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_load
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_load
            MODULE PROCEDURE ppm_map_part_load
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_ring_shift
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_ring_shift
            MODULE PROCEDURE ppm_map_part_ring_shift_s
            MODULE PROCEDURE ppm_map_part_ring_shift_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_store
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_store
            MODULE PROCEDURE ppm_map_part_store
         END INTERFACE

         PUBLIC :: ppm_map_part_eqdistrib
         PUBLIC :: ppm_map_part_cancel
         PUBLIC :: ppm_map_part_get_sub
         PUBLIC :: ppm_map_part_load
         PUBLIC :: ppm_map_part_ring_shift
         PUBLIC :: ppm_map_part_store

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_eqdistrib.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_eqdistrib.f"
#undef __KIND

#include "map/ppm_map_part_cancel.f"

#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_get_sub.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_get_sub.f"
#undef __KIND

#include "map/ppm_map_part_load.f"

#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_ring_shift.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_ring_shift.f"
#undef __KIND

#include "map/ppm_map_part_store.f"

      END MODULE ppm_module_map_part_util

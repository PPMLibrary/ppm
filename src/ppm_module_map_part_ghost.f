      !--*- f90 -*--------------------------------------------------------------
      !  Module       :            ppm_module_map_part_ghost
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

      MODULE ppm_module_map_part_ghost
      !!! This module provides the particle ghost mapping routines and holds the
      !!! temporary work arrays.

         IMPLICIT NONE

         PRIVATE

         !----------------------------------------------------------------------
         !  Work lists
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:), POINTER :: ilist1 => NULL()
         INTEGER, DIMENSION(:), POINTER :: ilist2 => NULL()
         INTEGER, DIMENSION(:), POINTER :: ighost => NULL()
         INTEGER                        :: prev_allocsize

         LOGICAL, DIMENSION(:), POINTER :: lghost => NULL()

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_ghost_get
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_ghost_get
            MODULE PROCEDURE ppm_map_part_ghost_get_d
            MODULE PROCEDURE ppm_map_part_ghost_get_s
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_ghost_pop
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_ghost_pop
            ! scalar (1d) particle data
            MODULE PROCEDURE ppm_map_part_ghost_pop_1dd
            MODULE PROCEDURE ppm_map_part_ghost_pop_1ds
            MODULE PROCEDURE ppm_map_part_ghost_pop_1di
            MODULE PROCEDURE ppm_map_part_ghost_pop_1dl
            MODULE PROCEDURE ppm_map_part_ghost_pop_1ddc
            MODULE PROCEDURE ppm_map_part_ghost_pop_1dsc

            ! vector (2d) particle data
            MODULE PROCEDURE ppm_map_part_ghost_pop_2dd
            MODULE PROCEDURE ppm_map_part_ghost_pop_2ds
            MODULE PROCEDURE ppm_map_part_ghost_pop_2di
            MODULE PROCEDURE ppm_map_part_ghost_pop_2dl
            MODULE PROCEDURE ppm_map_part_ghost_pop_2ddc
            MODULE PROCEDURE ppm_map_part_ghost_pop_2dsc
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_ghost_put
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_ghost_put
            MODULE PROCEDURE ppm_map_part_ghost_put
         END INTERFACE

         PUBLIC :: ppm_map_part_ghost_get
         PUBLIC :: ppm_map_part_ghost_pop
         PUBLIC :: ppm_map_part_ghost_put

      !----------------------------------------------------------------------
      !  include the source
      !----------------------------------------------------------------------
      CONTAINS

#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_ghost_get.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_ghost_get.f"
#undef __KIND

#define __DIM 1
#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __INTEGER
#include "map/ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __LOGICAL
#include "map/ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "map/ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "map/ppm_map_part_ghost_pop.f"
#undef __KIND
#undef __DIM

#define __DIM 2
#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __INTEGER
#include "map/ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __LOGICAL
#include "map/ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "map/ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "map/ppm_map_part_ghost_pop.f"
#undef __KIND
#undef __DIM

#include "map/ppm_map_part_ghost_put.f"

      END MODULE ppm_module_map_part_ghost

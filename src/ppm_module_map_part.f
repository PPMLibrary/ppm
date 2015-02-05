      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_map_part
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

      MODULE ppm_module_map_part
      !!! This module provides the basic mapping routines for particles, namely
      !!! the push, send and pop routine.
      !!!
      !!! The module also holds the needed work arrays for these routines.

      USE ppm_module_data, ONLY: ppm_kind_single,ppm_kind_double
      USE ppm_module_map_part_util, ONLY : ppm_map_part_eqdistrib,    &
      &   ppm_map_part_cancel,ppm_map_part_get_sub,ppm_map_part_load, &
      &   ppm_map_part_ring_shift,ppm_map_part_store
      USE ppm_module_map_part_ghost, ONLY : ppm_map_part_ghost_get, &
      &   ppm_map_part_ghost_pop,ppm_map_part_ghost_put
      USE ppm_module_map_part_global,ONLY : ppm_map_part_global, &
      &   ppm_map_part_remap
      USE ppm_module_map_part_partial, ONLY : ppm_map_part_partial
      IMPLICIT NONE

      PRIVATE

      !----------------------------------------------------------------------
      !  Work lists
      !----------------------------------------------------------------------
      REAL(ppm_kind_single), DIMENSION(:),   POINTER :: sends => NULL()
      REAL(ppm_kind_single), DIMENSION(:),   POINTER :: recvs => NULL()
      REAL(ppm_kind_double), DIMENSION(:),   POINTER :: sendd => NULL()
      REAL(ppm_kind_double), DIMENSION(:),   POINTER :: recvd => NULL()

      INTEGER,               DIMENSION(:),   POINTER :: nsend => NULL()
      INTEGER,               DIMENSION(:),   POINTER :: nrecv => NULL()
      INTEGER,               DIMENSION(:),   POINTER :: psend => NULL()
      INTEGER,               DIMENSION(:),   POINTER :: precv => NULL()
      INTEGER,               DIMENSION(:,:), POINTER :: pp    => NULL()
      INTEGER,               DIMENSION(:,:), POINTER :: qq    => NULL()
      INTEGER                                        :: old_nsendlist = 0
      INTEGER                                        :: old_buffer_set = 0

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

      INTERFACE ppm_map_part_isend
         MODULE PROCEDURE ppm_map_part_isend
      END INTERFACE

      PUBLIC :: ppm_map_part_pop
      PUBLIC :: ppm_map_part_push
      PUBLIC :: ppm_map_part_send
      PUBLIC :: ppm_map_part_isend
      !----------------------------------------------------------------------
      ! PUBLIC var from ppm_module_map_part_util
      !----------------------------------------------------------------------
      PUBLIC :: ppm_map_part_eqdistrib
      PUBLIC :: ppm_map_part_cancel
      PUBLIC :: ppm_map_part_get_sub
      PUBLIC :: ppm_map_part_load
      PUBLIC :: ppm_map_part_ring_shift
      PUBLIC :: ppm_map_part_store
      !----------------------------------------------------------------------
      ! PUBLIC var from ppm_module_map_part_ghost
      !----------------------------------------------------------------------
      PUBLIC :: ppm_map_part_ghost_get
      PUBLIC :: ppm_map_part_ghost_pop
      PUBLIC :: ppm_map_part_ghost_put
      !----------------------------------------------------------------------
      ! PUBLIC var from ppm_module_map_part_global
      !----------------------------------------------------------------------
      PUBLIC :: ppm_map_part_global
      PUBLIC :: ppm_map_part_remap
      !----------------------------------------------------------------------
      ! PUBLIC var from ppm_module_map_part_partial
      !----------------------------------------------------------------------
      PUBLIC :: ppm_map_part_partial
      !----------------------------------------------------------------------
      !  include the source
      !----------------------------------------------------------------------
      CONTAINS

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
#include "map/ppm_map_part_isend.f"

      END MODULE ppm_module_map_part

      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                  ppm_module_topo_alloc
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
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_topo_alloc
      !!! This module contains the topology allocation and deallocation
      !!! subroutines.

         IMPLICIT NONE

         !----------------------------------------------------------------------
         !  Define interfaces to topology allocation routine
         !----------------------------------------------------------------------
         INTERFACE ppm_topo_alloc
            MODULE PROCEDURE ppm_topo_alloc
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to topology deallocation routine
         !----------------------------------------------------------------------
         INTERFACE ppm_topo_dealloc
            MODULE PROCEDURE ppm_topo_dealloc
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "topo/ppm_topo_alloc.f"

#include "topo/ppm_topo_dealloc.f"

      END MODULE ppm_module_topo_alloc

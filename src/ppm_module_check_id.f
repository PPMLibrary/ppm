      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_check_topoid
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

      MODULE ppm_module_check_id
      !!! This module provides the utility routines,
      !!! which check topology and mesh IDs for their validity.
         IMPLICIT NONE
         !----------------------------------------------------------------------
         !  Define interface to the topoid check routine
         !----------------------------------------------------------------------
         INTERFACE ppm_check_topoid
         !!! checks topology ID
            MODULE PROCEDURE ppm_check_topoid
         END INTERFACE

         !!----------------------------------------------------------------------
         !!  Define interface to the meshid check routine
         !!----------------------------------------------------------------------
         !INTERFACE ppm_check_meshid
         !!!! checks mesh ID
            !MODULE PROCEDURE ppm_check_meshid
         !END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "topo/ppm_check_topoid.f"

#include "topo/ppm_check_meshid.f"

      END MODULE ppm_module_check_id

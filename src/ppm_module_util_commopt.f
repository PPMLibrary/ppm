      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_util_commopt
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

      MODULE ppm_module_util_commopt
      !!! This module provides the utility routines.

         IMPLICIT NONE

         !----------------------------------------------------------------------
         !  Define interface to the communication optimization routine
         !----------------------------------------------------------------------
         INTERFACE ppm_util_commopt
             MODULE PROCEDURE ppm_util_commopt
         END INTERFACE

         INTERFACE ppm_util_commopt_cart
             MODULE PROCEDURE ppm_util_commopt_cart
         END INTERFACE
         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "util/ppm_util_commopt.f"

#include "util/ppm_util_commopt_cart.f"
      END MODULE ppm_module_util_commopt

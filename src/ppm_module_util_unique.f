      !-------------------------------------------------------------------------
      !  Module       :                 ppm_module_util_unique
      !-------------------------------------------------------------------------
      ! Copyright (c) 2016 CSE Lab (ETH Zurich), MOSAIC Group (MPI-CBG Dresden),
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
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION           1
#define __DOUBLE_PRECISION           2
#define __INTEGER                    3
#define __LONGINT                    4

      MODULE ppm_module_util_unique
      !!! This module provides the utility unique routines.

      IMPLICIT NONE
      PRIVATE

      !----------------------------------------------------------------------
      !  Define interface to the list sorting routine
      !----------------------------------------------------------------------
      INTERFACE ppm_util_unique
        MODULE PROCEDURE ppm_util_unique_s
        MODULE PROCEDURE ppm_util_unique_d
        MODULE PROCEDURE ppm_util_unique_i
        MODULE PROCEDURE ppm_util_unique_li
        ! From a list of values generates a sorted unique list of values
      END INTERFACE

      PUBLIC :: ppm_util_unique

      !----------------------------------------------------------------------
      !  include the source
      !----------------------------------------------------------------------
      CONTAINS
#define __KIND __SINGLE_PRECISION
#include "util/ppm_util_unique.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "util/ppm_util_unique.f"
#undef __KIND

#define __KIND __INTEGER
#include "util/ppm_util_unique.f"
#undef __KIND

#define __KIND __LONGINT
#include "util/ppm_util_unique.f"
#undef __KIND

      END MODULE ppm_module_util_unique

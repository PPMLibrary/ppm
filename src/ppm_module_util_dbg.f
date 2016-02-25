      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_util_dbg
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
#define __SCALAR 3
#define __VECTOR 4

      MODULE ppm_module_util_dbg
      !!! This Module provides a simple means to visualize particles and
      !!! domain decompositions for debugging and monitoring purposes.
      !!!
      !!! The module is used in conjunction with the ppmdbg.py script.
      !!! It creates two files with names ppm_dbg_####.sub and ppm_dbg_####.dat
      !!! That contain the domain decomposition and particle data respectively

         IMPLICIT NONE

         !----------------------------------------------------------------------
         !  Define interfaces to the main routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_dbg_print
            MODULE PROCEDURE dbg_print_sca_s
            MODULE PROCEDURE dbg_print_sca_d
            MODULE PROCEDURE dbg_print_vec_s
            MODULE PROCEDURE dbg_print_vec_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS
#define __CTAG __SCALAR
#define __KIND __SINGLE_PRECISION
#include "util/ppm_dbg_print.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "util/ppm_dbg_print.f"
#undef __KIND
#undef __CTAG

#define __CTAG __VECTOR
#define __KIND __SINGLE_PRECISION
#include "util/ppm_dbg_print.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "util/ppm_dbg_print.f"
#undef __KIND
#undef __CTAG


      END MODULE ppm_module_util_dbg

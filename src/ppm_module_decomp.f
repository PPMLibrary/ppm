      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_decomp
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

      MODULE ppm_module_decomp
      !!! This module provides the PPM decomposition routines.
      !!!
      !!! PPM core provides 4 different decomposition routines that are called
      !!! within the topology routines.

         IMPLICIT NONE

         PRIVATE

         !----------------------------------------------------------------------
         !  Local work arrays used in ppm_decomp_pruned_cell
         !----------------------------------------------------------------------
         INTEGER , DIMENSION(:), POINTER :: npbx  => NULL()
         INTEGER , DIMENSION(:), POINTER :: npbxg => NULL()

         !----------------------------------------------------------------------
         !  Define interface to box split routine
         !----------------------------------------------------------------------
         INTERFACE ppm_decomp_boxsplit
            MODULE PROCEDURE decomp_bsplit_s
            MODULE PROCEDURE decomp_bsplit_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to cartesian decomposition routine
         !----------------------------------------------------------------------
         INTERFACE ppm_decomp_cartesian
            MODULE PROCEDURE decomp_cart_s
            MODULE PROCEDURE decomp_cart_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to pruned cell decomposition routine
         !----------------------------------------------------------------------
         INTERFACE ppm_decomp_pruned_cell
            MODULE PROCEDURE decomp_pcell_s
            MODULE PROCEDURE decomp_pcell_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to tree decomposition routine
         !----------------------------------------------------------------------
         INTERFACE ppm_decomp_tree
            MODULE PROCEDURE decomp_tree_s
            MODULE PROCEDURE decomp_tree_d
         END INTERFACE


         PUBLIC :: ppm_decomp_boxsplit
         PUBLIC :: ppm_decomp_cartesian
         PUBLIC :: ppm_decomp_pruned_cell
         PUBLIC :: ppm_decomp_tree

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "decomp/ppm_decomp_boxsplit.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "decomp/ppm_decomp_boxsplit.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION
#include "decomp/ppm_decomp_cartesian.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "decomp/ppm_decomp_cartesian.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION
#include "decomp/ppm_decomp_pruned_cell.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "decomp/ppm_decomp_pruned_cell.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION
#include "decomp/ppm_decomp_tree.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "decomp/ppm_decomp_tree.f"
#undef __KIND

END MODULE ppm_module_decomp

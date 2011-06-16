      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                   ppm_module_tree
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
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __TREE             3
#define __DECOMP           4

      MODULE ppm_module_tree
      !!! This module provides the tree routine.
         !----------------------------------------------------------------------
         !  Define interfaces to the routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_tree
            MODULE PROCEDURE ppm_tree_ts
            MODULE PROCEDURE ppm_tree_td
            MODULE PROCEDURE ppm_tree_ds
            MODULE PROCEDURE ppm_tree_dd
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __TYPE __TREE
#define __KIND __SINGLE_PRECISION
#include "tree/ppm_tree.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "tree/ppm_tree.f"
#undef __KIND
#undef __TYPE

#define __TYPE __DECOMP
#define __KIND __SINGLE_PRECISION
#include "tree/ppm_tree.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "tree/ppm_tree.f"
#undef __KIND
#undef __TYPE

      END MODULE ppm_module_tree

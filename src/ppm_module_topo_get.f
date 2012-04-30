      !--*- f90 -*--------------------------------------------------------------
      !  Module       :               ppm_module_topo_get
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

      MODULE ppm_module_topo_get
      !!! This module provides the routines to read an internally stored topology.
         !----------------------------------------------------------------------
         !  Define interfaces to topology store routine
         !----------------------------------------------------------------------
         INTERFACE ppm_topo_get
            MODULE PROCEDURE ppm_topo_get
         END INTERFACE
         
         INTERFACE ppm_topo_get_decomp
            MODULE PROCEDURE ppm_topo_get_decomp_s
            MODULE PROCEDURE ppm_topo_get_decomp_d
         END INTERFACE
         
         INTERFACE ppm_topo_get_meshinfo
            MODULE PROCEDURE ppm_topo_get_meshinfo
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "topo/ppm_topo_get.f"

#define __KIND __SINGLE_PRECISION
#include "topo/ppm_topo_get_decomp.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "topo/ppm_topo_get_decomp.f"
#undef __KIND

#include "topo/ppm_topo_get_meshinfo.f"

      END MODULE ppm_module_topo_get

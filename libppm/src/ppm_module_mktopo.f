      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_mktopo
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Lab (ETH Zurich), 
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

      MODULE ppm_module_mktopo
      !!! This module provides the interface to the routines that create
      !!! topologies.
         !----------------------------------------------------------------------
         !  Define interfaces to the main topology routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_mktopo
            MODULE PROCEDURE ppm_topo_mkpart_s
            MODULE PROCEDURE ppm_topo_mkpart_d
            MODULE PROCEDURE ppm_topo_mkfield_s
            MODULE PROCEDURE ppm_topo_mkfield_d
            MODULE PROCEDURE ppm_topo_mkgeom_s
            MODULE PROCEDURE ppm_topo_mkgeom_d
            MODULE PROCEDURE ppm_topo_mktree_s
            MODULE PROCEDURE ppm_topo_mktree_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "topo/ppm_topo_mkpart.f"
#include "topo/ppm_topo_mkfield.f"
#include "topo/ppm_topo_mkgeom.f"
#include "topo/ppm_topo_mktree.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "topo/ppm_topo_mkpart.f"
#include "topo/ppm_topo_mkfield.f"
#include "topo/ppm_topo_mkgeom.f"
#include "topo/ppm_topo_mktree.f"
#undef __KIND

      END MODULE ppm_module_mktopo

      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                   ppm_module_map
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

      MODULE ppm_module_map
      !!! This module contains all user-callable routines for mapping
      !!! particles and fields.
      !!!
      !!! You may either `USE` this module or directly the modules for the
      !!! subroutines you will be using.
         ! TODO: map_field_init added from Petros codebase (this might
         ! need revising, maybe transparent handling as done for ghost_init)
         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------

         USE ppm_module_map_connect
      !   USE ppm_module_map_field
         USE ppm_module_map_part
         USE ppm_module_impose_part_bc
      CONTAINS

#include "map/ppm_map_type_isactive.f"

      END MODULE ppm_module_map

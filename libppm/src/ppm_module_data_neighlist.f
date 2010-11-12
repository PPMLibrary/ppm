      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_data_neighlist
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
     
      MODULE ppm_module_data_neighlist
      !!! This module provides data used by the neighbor search routines.
      !!!
      !!! [NOTE]
      !!! The variables declared in this module should not be accessed by the
      !!! PPM client developer. They are managed internally by the library.
         !----------------------------------------------------------------------
         !  Define data TYPEs
         !----------------------------------------------------------------------
         ! Pointer to cell list (needed to make lists of cell lists)
         TYPE ppm_type_ptr_to_clist
             ! particle index list
             INTEGER, DIMENSION(:), POINTER    :: lpdx
             ! first particle in each cell
             INTEGER, DIMENSION(:), POINTER    :: lhbx
         END TYPE

      END MODULE ppm_module_data_neighlist

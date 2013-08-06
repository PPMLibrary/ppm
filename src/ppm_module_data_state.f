      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                ppm_module_data_state
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

      MODULE ppm_module_data_state
      !!! This module contains  data structures and definitions that
      !!! are PRIVATE to the ppm_map_part_store and ppm_map_part_load routines.
      !!!
      !!! [NOTE]
      !!! The variables declared in this modules should not be directly accessed
      !!! by the user. They are managed interally by the library.

         IMPLICIT NONE

         !----------------------------------------------------------------------
         !
         !----------------------------------------------------------------------
         INTEGER :: ppm_map_type_state
         INTEGER :: ppm_nrecvlist_state
         INTEGER :: ppm_nsendlist_state
         INTEGER :: ppm_nsendbuffer_state
         INTEGER :: ppm_buffer_set_state

         INTEGER, DIMENSION(:), POINTER :: ppm_psendbuffer_state => NULL()
         INTEGER, DIMENSION(:), POINTER :: ppm_buffer2part_state => NULL()

         INTEGER, DIMENSION(:), POINTER :: ppm_irecvlist_state => NULL()
         INTEGER, DIMENSION(:), POINTER :: ppm_isendlist_state => NULL()

      END MODULE ppm_module_data_state

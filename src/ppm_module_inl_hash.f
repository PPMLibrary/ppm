     !--*- f90 -*--------------------------------------------------------------
     !  Module   :                  ppm_module_inl_hash
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
      MODULE ppm_module_inl_hash
      !!! This module provides the utility to insert index of a cell and its position
      !!! on 'borders' array, as the range of indices of cells can be very too large
      !!! to allocate memory space for whole range. Hence, hash table is the
      !!! workaround for the redundancy in terms of memory consumption.

         USE ppm_module_data
         USE ppm_module_alloc
         USE ppm_module_error
        !-------------------------------------------------------------------------
        !  Declaration of parameters
        !-------------------------------------------------------------------------
        INTEGER,                 PARAMETER :: htable_null = -1
        !!! NULL value for hash table
        INTEGER(ppm_kind_int64), PARAMETER :: seed1 = 738235926
        !!! Hardcoded seed value taken from MurmurHash
        INTEGER(ppm_kind_int64), PARAMETER :: seed2 = 1243832038
        !!! Hardcoded seed value taken from MurmurHash

        TYPE ppm_htable
          !---------------------------------------------------------------------
          !  Declaration of arrays
          !---------------------------------------------------------------------
          INTEGER(ppm_kind_int64), DIMENSION(:), POINTER :: keys => NULL()
          !!! Array for keeping hash table keys.
          INTEGER,                 DIMENSION(:), POINTER :: borders_pos => NULL()
          !!! Array for keeping positions of cells on "borders" array.

          !---------------------------------------------------------------------
          !  Declaration of variables
          !--------------------------------------------------------------------
          INTEGER                                        :: nrow = 0
          !!! number of rows in hash table
!          CONTAINS
!              PROCEDURE :: create => create_htable
!              PROCEDURE :: destroy => destroy_htable
!              PROCEDURE :: insert => hash_insert
!              FUNCTION :: search => hash_search
        END TYPE

        PRIVATE :: seed1, seed2
!        PRIVATE :: create_htable, destroy_htable, hash_insert, hash_search
!        PRIVATE :: h_func, h_key

        CONTAINS

#include "neighlist/ppm_inl_hash.f"

      END MODULE ppm_module_inl_hash

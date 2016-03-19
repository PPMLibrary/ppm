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
         IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Declaration of parameters
        !-------------------------------------------------------------------------
        PRIVATE

        INTEGER,                 PARAMETER :: htable_null    = -HUGE(1)
        INTEGER(ppm_kind_int64), PARAMETER :: htable_null_li = -HUGE(1_ppm_kind_int64)
        !!! NULL value for hash table
        INTEGER,                 PARAMETER :: seed1 = 738235926
        !!! Hardcoded seed value taken from MurmurHash
        INTEGER,                 PARAMETER :: seed2 = 1243832038
        !!! Hardcoded seed value taken from MurmurHash

        ! TODO
        ! Extend the hashtable to account for longer keys (or arbitrary keys)
        ! This is useful for ppm applications

        TYPE ppm_htable
          !---------------------------------------------------------------------
          !  Declaration of arrays
          !---------------------------------------------------------------------
          INTEGER(ppm_kind_int64), DIMENSION(:), POINTER :: keys => NULL()
          !!! Array for keeping hash table keys.
          !!! Any key should not be bigger than 4 Bytes (32bits) value
          INTEGER,                 DIMENSION(:), POINTER :: borders_pos => NULL()
          !!! Array for keeping positions of cells on "borders" array.
          !---------------------------------------------------------------------
          !  Declaration of variables
          !--------------------------------------------------------------------
          INTEGER                                        :: nrow = 0
          !!! number of rows in hash table
         CONTAINS
          PROCEDURE :: create => create_htable
          PROCEDURE :: destroy => destroy_htable

          PROCEDURE :: h_func
          PROCEDURE :: h_key1
          PROCEDURE :: h_key2
          GENERIC   :: h_key => h_key1,h_key2

          PROCEDURE :: hash_insert
          PROCEDURE :: hash_insert_
          GENERIC   :: insert => hash_insert,hash_insert_

          PROCEDURE :: hash_search
          PROCEDURE :: hash_search_
          GENERIC   :: search => hash_search,hash_search_

          PROCEDURE :: hash_remove
          PROCEDURE :: hash_remove_
          GENERIC   :: remove => hash_remove,hash_remove_

          PROCEDURE :: grow => grow_htable
          PROCEDURE :: shrink => shrink_htable

          PROCEDURE :: size => hash_size
        END TYPE

        PUBLIC :: ppm_htable,htable_null

        CONTAINS

#include "neighlist/ppm_inl_hash.f"

      END MODULE ppm_module_inl_hash

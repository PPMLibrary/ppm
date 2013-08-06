      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_container_typedef
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
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_container_typedef
      !!! Declares container data types
      !!!
      !!! [NOTE]
      !!! Most of the declared variables in this module should not be accessed
      !!! directly by the PPM client developer, they are used internally in the
      !!! library.



         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_alloc
         USE ppm_module_data, ONLY: ppm_rank,ppm_dim,ppm_comm,ppm_char
         USE ppm_module_error
         USE ppm_module_write
         USE ppm_module_substart
         USE ppm_module_substop

         IMPLICIT NONE

         !----------------------------------------------------------------------
         ! Type declaration
         !----------------------------------------------------------------------
         !List of integers, stored as an array

#include "cont/container_typedef.inc"

         TYPE, EXTENDS(ppm_t_container) :: idList
             INTEGER, DIMENSION(:), POINTER :: vec => NULL()

         END TYPE idList
         !----------------------------------------------------------------------
         ! Global variables
         !----------------------------------------------------------------------
         CHARACTER(LEN=ppm_char),PRIVATE :: cbuf
         INTEGER, PRIVATE, DIMENSION(3)  :: ldc
         !!! Number of elements in all dimensions for allocation


!         !----------------------------------------------------------------------
!         ! Type-bound procedures
!         !----------------------------------------------------------------------
!         CONTAINS
!
!
!#define DTYPE(a) a/**/_s
!#include "cont/container_typeproc.f"
!
!#define DTYPE(a) a/**/_d
!#include "cont/container_typeproc.f"



         END MODULE ppm_module_container_typedef


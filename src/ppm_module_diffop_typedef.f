      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_diffop_typedef
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

#define __REAL 3 
#define __COMPLEX 4 
#define __INTEGER 5 
#define __LONGINT 6 
#define __LOGICAL 7 
#define __CHAR 8 

#define __crash_on_null_pointers  1
#undef __WITH_KDTREE

      MODULE ppm_module_diffop_typedef
      !!! Declares differential operator data types
      !!!
      !!! [NOTE]
      !!! Most of the declared variables in this module should not be accessed
      !!! directly by the PPM client developer, they are used internally in the
      !!! library.

         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_alloc
         USE ppm_module_data
         !USE ppm_module_typedef
         USE ppm_module_container_typedef
         USE ppm_module_error
         USE ppm_module_write
         USE ppm_module_substart
         USE ppm_module_substop

         IMPLICIT NONE
         !----------------------------------------------------------------------
         ! Global parameters
         !----------------------------------------------------------------------

         !----------------------------------------------------------------------
         ! Type declaration
         !----------------------------------------------------------------------

#define  DTYPE(a) a/**/_s
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#include "diffop/diffop_typedef.inc"

#define  DTYPE(a) a/**/_d
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "diffop/diffop_typedef.inc"

         !----------------------------------------------------------------------
         ! Type-bound procedures
         !----------------------------------------------------------------------


         CONTAINS


#define DTYPE(a) a/**/_s
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#include "diffop/diffop_typeproc.f"

#define DTYPE(a) a/**/_d
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#include "diffop/diffop_typeproc.f"


         END MODULE ppm_module_diffop_typedef


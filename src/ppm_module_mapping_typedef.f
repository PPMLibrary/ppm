      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_mapping_typedef
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
     
#define __SINGLE_PRECISION         1
#define __DOUBLE_PRECISION         2
#define __INTEGER                  3
#define __LOGICAL                  4
#define __SINGLE_PRECISION_COMPLEX 5
#define __DOUBLE_PRECISION_COMPLEX 6

      MODULE ppm_module_mapping_typedef
      !!! This module defines the mapping data type and provides 
      !!!  the basic mapping routines for particles and meshes; namely
      !!!  the push, send and pop routine.
         
      
         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_map_part_util
         USE ppm_module_map_part_ghost
         USE ppm_module_map_part_global
         USE ppm_module_map_part_partial
         USE ppm_module_data, ONLY: ppm_kind_single,ppm_kind_double
         USE ppm_module_substart
         USE ppm_module_substop

         IMPLICIT NONE
         PRIVATE

         PUBLIC :: ppm_ptr_t_part_mapping_s,ppm_ptr_t_part_mapping_d
         PUBLIC :: ppm_t_part_mapping_s,ppm_t_part_mapping_d


         !----------------------------------------------------------------------
         !  Type declaration
         !----------------------------------------------------------------------

#define  DTYPE(a) a/**/_s
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#include "map/mapping_typedef.inc"

#define  DTYPE(a) a/**/_d
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "map/mapping_typedef.inc"

         !----------------------------------------------------------------------
         !  Type-bound procedures
         !----------------------------------------------------------------------
         CONTAINS

#define  DTYPE(a) a/**/_s
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#include "map/mapping_typeproc.f"

#define  DTYPE(a) a/**/_d
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#include "map/mapping_typeproc.f"

      END MODULE ppm_module_mapping_typedef

!--*- f90 -*--------------------------------------------------------------
!  Module   :                   ppm_module_field_typedef
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

MODULE ppm_module_field_typedef
!!! Declares field data type

!----------------------------------------------------------------------
!  Modules
!----------------------------------------------------------------------
USE ppm_module_substart
USE ppm_module_substop
USE ppm_module_topo_typedef
USE ppm_module_interfaces
USE ppm_module_data
USE ppm_module_alloc
USE ppm_module_error
USE ppm_module_util_functions

IMPLICIT NONE

!----------------------------------------------------------------------
! Internal parameters
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Global variables 
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Module variables 
!----------------------------------------------------------------------
INTEGER, PRIVATE, DIMENSION(3)  :: ldc

!----------------------------------------------------------------------
! Type declaration
!----------------------------------------------------------------------
TYPE, EXTENDS(ppm_t_mesh_data_):: ppm_t_mesh_data
    CONTAINS
    PROCEDURE :: create => mesh_data_create
    PROCEDURE :: destroy => mesh_data_destroy

END TYPE ppm_t_mesh_data
!----------------------------------------------------------------------
! Container for mesh_data
!----------------------------------------------------------------------
#define CONTAINER ppm_c_mesh_data
#define __CONTAINER(a) ppm_c_mesh_data_/**/a
#define VEC_TYPE ppm_t_mesh_data
#include "cont/collection_template.inc"

TYPE,EXTENDS(ppm_t_field_) :: ppm_t_field
    CONTAINS
    PROCEDURE :: create => field_create
    PROCEDURE :: destroy => field_destroy
    PROCEDURE :: discretize_on => field_discretize_on
END TYPE ppm_t_field

!----------------------------------------------------------------------
!  Type-bound procedures
!----------------------------------------------------------------------
CONTAINS

#define CONTAINER ppm_c_mesh_data
#define __CONTAINER(a) ppm_c_mesh_data_/**/a
#define VEC_TYPE ppm_t_mesh_data
#include "cont/collection_typeproc.f"

#include "field/field_typeproc.f"


END MODULE ppm_module_field_typedef

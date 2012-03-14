!--*- f90 -*--------------------------------------------------------------
!  Module   :                   ppm_module_mesh_typedef
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

MODULE ppm_module_mesh_typedef
!!! Declares mesh data types
!!!
!!! [NOTE]
!!! Most of the declared variables in this module should not be accessed
!!! directly by the PPM client developer, they are used internally in the
!!! library.

!----------------------------------------------------------------------
!  Modules
!----------------------------------------------------------------------
USE ppm_module_substart
USE ppm_module_substop
USE ppm_module_data
USE ppm_module_alloc
USE ppm_module_error
USE ppm_module_util_functions
USE ppm_module_interfaces

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

!!----------------------------------------------------------------------
!! Patches (contains the actual data arrays for this field)
!!----------------------------------------------------------------------
TYPE,EXTENDS(ppm_t_subpatch_data_) :: ppm_t_subpatch_data
    CONTAINS
    PROCEDURE :: create    => subpatch_data_create
    PROCEDURE :: destroy   => subpatch_data_destroy
END TYPE

!----------------------------------------------------------------------
! Container for lists of (pointers to) subpatch_data
!----------------------------------------------------------------------
#define CONTAINER ppm_c_subpatch_data
#define __CONTAINER(a) ppm_c_subpatch_data_/**/a
#define VEC_TYPE ppm_t_subpatch_data
#include "cont/collection_template.inc"

TYPE,EXTENDS(ppm_t_subpatch_) :: ppm_t_subpatch
    CONTAINS
    PROCEDURE  :: create    => subpatch_create
    PROCEDURE  :: destroy   => subpatch_destroy

    PROCEDURE  :: subpatch_get_field_2d_rd
    PROCEDURE  :: subpatch_get_field_3d_rd
    GENERIC    :: get_field => &
        &                   subpatch_get_field_2d_rd,&
        &                   subpatch_get_field_3d_rd
    !PROCEDURE  :: get => subpatch_get
END TYPE
!----------------------------------------------------------------------
! Container for lists of (pointers to) subpatch
!----------------------------------------------------------------------
#define CONTAINER ppm_c_subpatch
#define __CONTAINER(a) ppm_c_subpatch_/**/a
#define VEC_TYPE ppm_t_subpatch
#include "cont/collection_template.inc"



TYPE,EXTENDS(ppm_t_equi_mesh_) :: ppm_t_equi_mesh
    CONTAINS
    PROCEDURE  :: create    => equi_mesh_create
    PROCEDURE  :: destroy   => equi_mesh_destroy
    PROCEDURE  :: add_patch => equi_mesh_add_patch
END TYPE

!----------------------------------------------------------------------
! Container for lists of (pointers to) Meshes
!----------------------------------------------------------------------
#define CONTAINER ppm_c_meshes
#define __CONTAINER(a) ppm_c_meshes_/**/a
#define VEC_TYPE ppm_t_equi_mesh
#include "cont/collection_template.inc"

!----------------------------------------------------------------------
! DATA STORAGE for the meshes
!----------------------------------------------------------------------
TYPE(ppm_c_meshes)                 :: ppm_mesh
!!! container for PPM meshes data structures
!----------------------------------------------------------------------
!  Type-bound procedures
!----------------------------------------------------------------------
CONTAINS

#define CONTAINER ppm_c_meshes
#define __CONTAINER(a) ppm_c_meshes_/**/a
#define VEC_TYPE ppm_t_equi_mesh
#include "cont/collection_typeproc.f"

#define CONTAINER ppm_c_subpatch_data
#define __CONTAINER(a) ppm_c_subpatch_data_/**/a
#define VEC_TYPE ppm_t_subpatch_data
#include "cont/collection_typeproc.f"

#define CONTAINER ppm_c_subpatch
#define __CONTAINER(a) ppm_c_subpatch_/**/a
#define VEC_TYPE ppm_t_subpatch
#include "cont/collection_typeproc.f"


#include "mesh/mesh_typeproc.f"

!#define  DTYPE(a) a/**/_s
!#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
!#include "mesh/mesh_data_typeproc.f"
!
!#define  DTYPE(a) a/**/_d
!#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
!#include "mesh/mesh_data_typeproc.f"
!

END MODULE ppm_module_mesh_typedef

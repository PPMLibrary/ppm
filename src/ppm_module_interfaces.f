!--*- f90 -*--------------------------------------------------------------
!  Module   :                   ppm_module_typedef
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

MODULE ppm_module_interfaces
!!! Declares all data types 
!!! The derived types are declared as abstract. They contain all the
!!! different fields as well as the interfaces to the type-bound
!!! procedures.
!!! Each type is then extended in another module, which also
!!! contains the implementations of the type-bound procedures.
!!! This allows for type-bound procedures to take arguments of any
!!! other derived type (because they only need to use the module
!!! with the abstract types and there is no risk of circular
!!! dependencies between modules).

!----------------------------------------------------------------------
!  Modules
!----------------------------------------------------------------------
USE ppm_module_container_typedef
USE ppm_module_data
USE ppm_module_util_functions

IMPLICIT NONE

!----------------------------------------------------------------------
! Internal parameters
!----------------------------------------------------------------------
!PPM internal parameters used only to access entries in the
!mesh data structures.
INTEGER,PARAMETER,PUBLIC   :: ppm_mdata_ghosts = 1
INTEGER,PARAMETER,PUBLIC   :: ppm_mdata_reqput = 2
INTEGER,PARAMETER,PUBLIC   :: ppm_mdata_cartesian = 3
INTEGER,PARAMETER,PUBLIC   :: ppm_mdata_lflags = 3

!PPM internal parameters used only to access entries in the
!mesh data structures.
INTEGER,PARAMETER,PUBLIC   :: ppm_mesh_ghosts = 1
INTEGER,PARAMETER,PUBLIC   :: ppm_mesh_reqput = 2
INTEGER,PARAMETER,PUBLIC   :: ppm_mesh_cartesian = 3
INTEGER,PARAMETER,PUBLIC   :: ppm_mesh_length_partflags = 3

!PPM internal parameters used only to access entries in the
!particle's property data structure.
INTEGER,PARAMETER,PUBLIC   :: ppm_mdat_ghosts = 1
INTEGER,PARAMETER,PUBLIC   :: ppm_mdat_reqput = 2
INTEGER,PARAMETER,PUBLIC   :: ppm_mdat_map_ghosts = 3
INTEGER,PARAMETER,PUBLIC   :: ppm_param_length_mdatflags = 3

!User parameters 
INTEGER,PARAMETER,PUBLIC   :: ppm_param_part_init_cartesian = 1
INTEGER,PARAMETER,PUBLIC   :: ppm_param_part_init_random = 2


!PPM internal parameters used only to access entries in the
!particle data structure.
INTEGER,PARAMETER,PUBLIC   :: ppm_part_ghosts = 1
INTEGER,PARAMETER,PUBLIC   :: ppm_part_partial = 2
INTEGER,PARAMETER,PUBLIC   :: ppm_part_reqput = 3
INTEGER,PARAMETER,PUBLIC   :: ppm_part_areinside = 4
INTEGER,PARAMETER,PUBLIC   :: ppm_part_cartesian = 5
INTEGER,PARAMETER,PUBLIC   :: ppm_part_neighlists = 6
INTEGER,PARAMETER,PUBLIC   :: ppm_part_global_index = 7
INTEGER,PARAMETER,PUBLIC   :: ppm_param_length_partflags = 7

!PPM internal parameters used only to access entries in the
!particle's property data structure.
INTEGER,PARAMETER,PUBLIC   :: ppm_ppt_ghosts = 1
INTEGER,PARAMETER,PUBLIC   :: ppm_ppt_partial = 2
INTEGER,PARAMETER,PUBLIC   :: ppm_ppt_reqput = 3
INTEGER,PARAMETER,PUBLIC   :: ppm_ppt_map_parts = 4
INTEGER,PARAMETER,PUBLIC   :: ppm_ppt_map_ghosts = 5
INTEGER,PARAMETER,PUBLIC   :: ppm_param_length_pptflags = 5

!PPM internal parameters for default storage IDs of some DS.
INTEGER, PARAMETER,PUBLIC :: ppm_param_default_nlID = 1

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
#define  DTYPE(a) a/**/_s
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#include "map/mapping_abstract_typedef.f"

#define  DTYPE(a) a/**/_d
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "map/mapping_abstract_typedef.f"

#define  DTYPE(a) a/**/_s
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#include "part/particles_abstract_typedef.f"

#define  DTYPE(a) a/**/_d
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "part/particles_abstract_typedef.f"

TYPE,ABSTRACT ::  ppm_t_mesh_data_
    !!! Data structure for mesh data. 
    !!! (Contained inside a ppm_t_field, so relates to one specific
    !!!  field, denoted by fieldID)
    !!! 
    !!! Contains pointers to the data and bookkeeping information
    !!! for each mesh on which fieldID has been discretized.

    INTEGER                                          :: meshID
    !!! ID of the mesh on which fieldID is discretized
    INTEGER                                          :: lda
    !!! number of components (1 for scalar fields)
    INTEGER                                          :: p_idx
    !!! Storage index for of the data array on the subpatches of this mesh
    !!! (A mesh stores data for several fields. Each subpatch thus 
    !!!  contains a list of data arrays, corresponding to the different
    !!!  fields. The p_idx index gives the location of the data array
    !!!  corresponding to fieldID within this list)
    LOGICAL,DIMENSION(ppm_mdata_lflags)              :: flags = .FALSE.
    !!! Booleans used to track the state of this discretization.

    CONTAINS
    PROCEDURE(mesh_data_create_),DEFERRED :: create
    PROCEDURE(mesh_data_destroy_),DEFERRED :: destroy
END TYPE
!----------------------------------------------------------------------
! Container for mesh_data
!----------------------------------------------------------------------
#define CONTAINER ppm_c_mesh_data_
#define __CONTAINER(a) ppm_c_mesh_data__/**/a
#define VEC_TYPE ppm_t_mesh_data_
#include "cont/collection_abstract_template.inc"

TYPE,ABSTRACT :: ppm_t_field_
    !!! Data structure for fields 
    !!! A field represents a mathematical concept (e.g. velocity, or
    !!! vorticity) and links to its discretized representation on meshes 
    !!! and/or on particles.

    INTEGER                                          :: fieldID
    !!! global identifier 
    CHARACTER(LEN=ppm_char)                          :: name
    !!! string description
    INTEGER                                          :: lda
    !!! number of components (1 for scalar fields)
    !!!
    !!! pointers to arrays where the scalar-value properties are stored
    CLASS(ppm_c_mesh_data_),ALLOCATABLE              :: M
    !!! Collection of pointers to the data and bookkeeping information
    !!! for each mesh on which this field has been discretized.
    ! CLASS(ppm_c_part_data_),ALLOCATABLE            :: P
    !    !!! Collection of pointers to the data and bookkeeping information
    !    !!! for each particle set on which this field has been discretized.

    CONTAINS
    PROCEDURE(field_create_),DEFERRED :: create
    PROCEDURE(field_destroy_),DEFERRED :: destroy
    PROCEDURE(field_discretize_on_),DEFERRED :: discretize_on
END TYPE ppm_t_field_

!!----------------------------------------------------------------------
!! Patches (contains the actual data arrays for this field)
!!----------------------------------------------------------------------
TYPE,ABSTRACT :: ppm_t_subpatch_data_
    !!! pointers to arrays where the data are stored
    INTEGER, DIMENSION(:,:), POINTER          :: data_2d_i => NULL()
    !!! if the data is 2d int
    INTEGER, DIMENSION(:,:,:), POINTER        :: data_3d_i => NULL()
    !!! if the data is 3d int
    INTEGER, DIMENSION(:,:,:,:), POINTER      :: data_4d_i => NULL()
    !!! if the data is 4d int
    REAL(ppm_kind_single), DIMENSION(:,:), POINTER         :: data_2d_rs => NULL()
    !!! if the data is 2d real
    REAL(ppm_kind_single), DIMENSION(:,:,:), POINTER       :: data_3d_rs => NULL()
    !!! if the data is 3d real
    REAL(ppm_kind_single), DIMENSION(:,:,:,:), POINTER     :: data_4d_rs => NULL()
    !!! if the data is 4d real
    COMPLEX(ppm_kind_single), DIMENSION(:,:), POINTER      :: data_2d_cs => NULL()
    !!! if the data is 2d complex
    COMPLEX(ppm_kind_single), DIMENSION(:,:,:), POINTER    :: data_3d_cs => NULL()
    !!! if the data is 3d complex
    COMPLEX(ppm_kind_single), DIMENSION(:,:,:,:), POINTER  :: data_4d_cs => NULL()
    !!! if the data is 4d complex
    REAL(ppm_kind_double), DIMENSION(:,:), POINTER         :: data_2d_rd => NULL()
    !!! if the data is 2d real
    REAL(ppm_kind_double), DIMENSION(:,:,:), POINTER       :: data_3d_rd => NULL()
    !!! if the data is 3d real
    REAL(ppm_kind_double), DIMENSION(:,:,:,:), POINTER     :: data_4d_rd => NULL()
    !!! if the data is 4d real
    COMPLEX(ppm_kind_double), DIMENSION(:,:), POINTER      :: data_2d_cd => NULL()
    !!! if the data is 2d complex
    COMPLEX(ppm_kind_double), DIMENSION(:,:,:), POINTER    :: data_3d_cd => NULL()
    !!! if the data is 3d complex
    COMPLEX(ppm_kind_double), DIMENSION(:,:,:,:), POINTER  :: data_4d_cd => NULL()
    !!! if the data is 4d complex
    LOGICAL, DIMENSION(:,:), POINTER          :: data_2d_l => NULL()
    !!! if the data is 2d logical
    LOGICAL, DIMENSION(:,:,:), POINTER        :: data_3d_l => NULL()
    !!! if the data is 3d logical
    LOGICAL, DIMENSION(:,:,:,:), POINTER      :: data_4d_l => NULL()
    !!! if the data is 4d logical

    CONTAINS
    PROCEDURE(subpatch_data_create_),DEFERRED  :: create
    PROCEDURE(subpatch_data_destroy_),DEFERRED :: destroy 
END TYPE

!----------------------------------------------------------------------
! Container for lists of (pointers to) subpatch_data
!----------------------------------------------------------------------
#define CONTAINER ppm_c_subpatch_data_
#define __CONTAINER(a) ppm_c_subpatch_data__/**/a
#define VEC_TYPE ppm_t_subpatch_data_
#include "cont/collection_abstract_template.inc"

TYPE,ABSTRACT :: ppm_t_subpatch_
    INTEGER               :: meshID
    !!! MeshID to which this subpatch belongs
    INTEGER, DIMENSION(:),POINTER :: istart => NULL()
    !!! Lower-left coordinates
    INTEGER, DIMENSION(:),POINTER :: iend => NULL()
    !!! Upper-right coordinates
    CLASS(ppm_c_subpatch_data_),ALLOCATABLE :: subpatch_data
    !!! container for the data arrays for each property discretized
    !!! on this mesh
    CONTAINS
    PROCEDURE(subpatch_create_),DEFERRED   :: create
    PROCEDURE(subpatch_destroy_),DEFERRED  :: destroy
    PROCEDURE(subpatch_get_field_2d_rd_), DEFERRED :: subpatch_get_field_2d_rd
    PROCEDURE(subpatch_get_field_3d_rd_), DEFERRED :: subpatch_get_field_3d_rd
    GENERIC :: get_field => subpatch_get_field_2d_rd,subpatch_get_field_3d_rd
    !PROCEDURE  :: get => subpatch_get
END TYPE
!----------------------------------------------------------------------
! Container for lists of (pointers to) subpatch
!----------------------------------------------------------------------
#define CONTAINER ppm_c_subpatch_
#define __CONTAINER(a) ppm_c_subpatch__/**/a
#define VEC_TYPE ppm_t_subpatch_
#include "cont/collection_abstract_template.inc"

TYPE ppm_t_ptr_subpatch
    CLASS(ppm_t_subpatch_), POINTER :: t => NULL()
END TYPE

TYPE ppm_t_A_subpatch
    TYPE(ppm_t_ptr_subpatch),DIMENSION(:), POINTER :: subpatch => NULL()
END TYPE

TYPE ppm_t_mesh_maplist
    !!! TODO: check what this is used for (imported from Petros code
    INTEGER, POINTER  :: target_topoid => NULL()
    !!! target topology ID
    INTEGER, POINTER  :: target_meshid => NULL()
    !!! target mesh ID
    INTEGER, POINTER  :: nsendlist => NULL()
    !!! send rank list size
    INTEGER, POINTER  :: nrecvlist => NULL()
    !!! recv rank size
    INTEGER, DIMENSION(:), POINTER  :: isendlist => NULL()
    !!! send rank lists
    INTEGER, DIMENSION(:), POINTER  :: irecvlist => NULL()
    !!! recv rank lists
    INTEGER, DIMENSION(:), POINTER  :: isendfromsub => NULL()
    !!! source send sublist
    INTEGER, DIMENSION(:,:), POINTER  :: isendblkstart => NULL()
    !!! send block start list
    INTEGER, DIMENSION(:,:), POINTER  :: isendblksize  => NULL()
    !!! send block size list
    INTEGER, DIMENSION(:), POINTER  :: irecvtosub   => NULL()
    !!! recv send sublist
    INTEGER, DIMENSION(:,:), POINTER  :: irecvblkstart => NULL()
    !!! recv block start list
    INTEGER, DIMENSION(:,:), POINTER  :: irecvblksize => NULL()
    !!! recv block size list
    INTEGER, DIMENSION(:), POINTER  :: psendbuffer => NULL()
    !!! send buffer pointer
    INTEGER, DIMENSION(:), POINTER  :: precvbuffer => NULL()
    !!! recv buffer pointer
END TYPE

TYPE,ABSTRACT :: ppm_t_equi_mesh_
    !!! Type for equispaced cartesian meshes on subs
    INTEGER                           :: ID
    !!! ID of the mesh in the belonging topology
    !!! It is the same as its index in the ppm_t_topo%mesh array
    INTEGER                           :: topoid
    !!! ID of the topology for which this mesh is defined

    INTEGER, DIMENSION(:), POINTER    :: Nm    => NULL()
    !!! global number of mesh points in computational domain

    REAL(ppm_kind_double),DIMENSION(:),POINTER :: Offset => NULL()
    !!! Offset

    REAL(ppm_kind_double),DIMENSION(:),POINTER :: h => NULL()
    !!! mesh spacing

    !TODO : delete those 2 guys
    INTEGER, DIMENSION(:,:), POINTER    :: nnodes    => NULL()
    INTEGER, DIMENSION(:,:), POINTER    :: istart    => NULL()

    CLASS(ppm_c_subpatch_),ALLOCATABLE            :: subpatch
    !!! container for subdomains patches 

    TYPE(ppm_t_A_subpatch), DIMENSION(:), POINTER :: patch => NULL()
    !!! array of arrays of pointers to the subpatches for each patch

    TYPE(ppm_t_A_subpatch), DIMENSION(:), POINTER :: sub => NULL()
    !!! array of arrays of pointers to the subpatches for each sub.


    !------------------------------------------------------------------
    !  Mesh ghosts mappings
    !------------------------------------------------------------------
    LOGICAL                          :: ghost_initialized = .FALSE.
    !!! is .TRUE. if the ghost mappings have been initialized
    !!! else, .FALSE.
    INTEGER, DIMENSION(:),   POINTER :: ghost_fromsub => NULL()
    !!! list of source subs of ghost mesh blocks (globel sub number).
    !!! These are the owner subs of the actual real mesh points
    !!! 1st index: meshblock ID
    INTEGER, DIMENSION(:),   POINTER :: ghost_tosub   => NULL()
    !!! list of target subs of ghost mesh blocks (globel sub number).
    !!! These are the subs a block will serve as a ghost on.
    !!! 1st index: meshblock ID
    INTEGER, DIMENSION(:,:), POINTER :: ghost_blkstart => NULL()
    !!! start (lower-left corner) of ghost mesh block in GLOBAL
    !!! mesh coordinates. First index: x,y[,z], 2nd: meshblock ID
    INTEGER, DIMENSION(:,:), POINTER :: ghost_blksize  => NULL()
    !!! size (in grid points) of ghost blocks. 1st index: x,y[,z], 2nd:
    !!! meshblock ID
    INTEGER, DIMENSION(:)  , POINTER :: ghost_blk      => NULL()
    !!! mesh ghost block list. 1st index: target processor
    INTEGER                          :: ghost_nsend
    !!! number of mesh blocks to be sent as ghosts
    INTEGER                          :: ghost_nrecv
    !!! number of mesh blocks to be recvd as ghosts
    INTEGER, DIMENSION(:), POINTER   :: ghost_recvtosub => NULL()
    !!! list of target subs for ghost mesh blocks to be received,
    !!! i.e. being ghost on the local processor (globel sub number).
    !!! These are the subs where the blocks will serve as ghosts
    !!! 1st index: meshblock ID
    INTEGER, DIMENSION(:,:), POINTER  :: ghost_recvblkstart => NULL()
    !!! start (lower-left corner) of received ghost mesh block in
    !!! GLOBAL  mesh coordinates. 1st index: x,y[,z], 2nd: meshblock ID
    INTEGER, DIMENSION(:,:), POINTER  :: ghost_recvblksize => NULL()
    !!! size (in grid points) of recvd ghost blocks.
    !!! 1st index: x,y[,z], 2nd: meshblock ID
    INTEGER, DIMENSION(:)  , POINTER  :: ghost_recvblk => NULL()
    !!! mesh ghost block receive list. 1st index: target processor


    TYPE(ppm_t_mesh_maplist), POINTER :: mapping => NULL()

    CONTAINS
    PROCEDURE(equi_mesh_create_),   DEFERRED :: create
    PROCEDURE(equi_mesh_destroy_),  DEFERRED :: destroy
    PROCEDURE(equi_mesh_add_patch_),DEFERRED :: add_patch
END TYPE


!----------------------------------------------------------------------
! Container for Meshes
!----------------------------------------------------------------------
#define CONTAINER ppm_c_meshes_
#define __CONTAINER(a) ppm_c_meshes__/**/a
#define VEC_TYPE ppm_t_equi_mesh_
#include "cont/collection_abstract_template.inc"

!----------------------------------------------------------------------
!  INTERFACES
!----------------------------------------------------------------------

INTERFACE
#define  DTYPE(a) a/**/_s
#include "map/mapping_interfaces.f"

#define  DTYPE(a) a/**/_d
#include "map/mapping_interfaces.f"

#define  DTYPE(a) a/**/_s
#define  MK ppm_kind_single
#include "part/particles_interfaces.f"

#define  DTYPE(a) a/**/_d
#define  MK ppm_kind_double
#include "part/particles_interfaces.f"

#define DTYPE(a) a/**/_s
#define __KIND __SINGLE_PRECISION
#define __DIM 1
#define __TYPE INTEGER
#define __MYTYPE __INTEGER
#define DATANAME data_1d_i
#include "part/particles_get_interfaces.f"
#define __TYPE INTEGER(ppm_kind_int64)
#define __MYTYPE __LONGINT
#define DATANAME data_1d_li
#include "part/particles_get_interfaces.f"
#define __TYPE REAL(ppm_kind_single)
#define __MYTYPE __REAL
#define DATANAME data_1d_r
#include "part/particles_get_interfaces.f"
#define __TYPE COMPLEX(ppm_kind_single)
#define __MYTYPE __COMPLEX
#define DATANAME data_1d_c
#include "part/particles_get_interfaces.f"
#define __TYPE LOGICAL
#define __MYTYPE __LOGICAL
#define DATANAME data_1d_l
#include "part/particles_get_interfaces.f"
#undef  __DIM

#define __DIM 2
#define __TYPE INTEGER
#define __MYTYPE __INTEGER
#define DATANAME data_2d_i
#include "part/particles_get_interfaces.f"
#define __TYPE INTEGER(ppm_kind_int64)
#define __MYTYPE __LONGINT
#define DATANAME data_2d_li
#include "part/particles_get_interfaces.f"
#define __TYPE REAL(ppm_kind_single)
#define __MYTYPE __REAL
#define DATANAME data_2d_r
#include "part/particles_get_interfaces.f"
#define __TYPE COMPLEX(ppm_kind_single)
#define __MYTYPE __COMPLEX
#define DATANAME data_2d_c
#include "part/particles_get_interfaces.f"
#define __TYPE LOGICAL
#define __MYTYPE __LOGICAL
#define DATANAME data_2d_l
#include "part/particles_get_interfaces.f"
#undef  __DIM
#undef DTYPE
#undef __KIND

#define DTYPE(a) a/**/_d
#define __KIND __DOUBLE_PRECISION
#define __DIM 1
#define __TYPE INTEGER
#define __MYTYPE __INTEGER
#define DATANAME data_1d_i
#include "part/particles_get_interfaces.f"
#define __TYPE INTEGER(ppm_kind_int64)
#define __MYTYPE __LONGINT
#define DATANAME data_1d_li
#include "part/particles_get_interfaces.f"
#define __TYPE REAL(ppm_kind_double)
#define __MYTYPE __REAL
#define DATANAME data_1d_r
#include "part/particles_get_interfaces.f"
#define __TYPE COMPLEX(ppm_kind_double)
#define __MYTYPE __COMPLEX
#define DATANAME data_1d_c
#include "part/particles_get_interfaces.f"
#define __TYPE LOGICAL
#define __MYTYPE __LOGICAL
#define DATANAME data_1d_l
#include "part/particles_get_interfaces.f"
#undef  __DIM

#define __DIM 2
#define __TYPE INTEGER
#define __MYTYPE __INTEGER
#define DATANAME data_2d_i
#include "part/particles_get_interfaces.f"
#define __TYPE INTEGER(ppm_kind_int64)
#define __MYTYPE __LONGINT
#define DATANAME data_2d_li
#include "part/particles_get_interfaces.f"
#define __TYPE REAL(ppm_kind_double)
#define __MYTYPE __REAL
#define DATANAME data_2d_r
#include "part/particles_get_interfaces.f"
#define __TYPE COMPLEX(ppm_kind_double)
#define __MYTYPE __COMPLEX
#define DATANAME data_2d_c
#include "part/particles_get_interfaces.f"
#define __TYPE LOGICAL
#define __MYTYPE __LOGICAL
#define DATANAME data_2d_l
#include "part/particles_get_interfaces.f"
#undef  __DIM
#undef DTYPE
#undef __KIND

#include "field/field_interfaces.f"

#include "mesh/mesh_interfaces.f"


#define CONTAINER ppm_c_mesh_data_
#define __CONTAINER(a) ppm_c_mesh_data__/**/a
#define VEC_TYPE ppm_t_mesh_data_
#include "cont/collection_interfaces.f"

#define CONTAINER ppm_c_meshes_
#define __CONTAINER(a) ppm_c_meshes__/**/a
#define VEC_TYPE ppm_t_equi_mesh_
#include "cont/collection_interfaces.f"

#define CONTAINER ppm_c_subpatch_data_
#define __CONTAINER(a) ppm_c_subpatch_data__/**/a
#define VEC_TYPE ppm_t_subpatch_data_
#include "cont/collection_interfaces.f"

#define CONTAINER ppm_c_subpatch_
#define __CONTAINER(a) ppm_c_subpatch__/**/a
#define VEC_TYPE ppm_t_subpatch_
#include "cont/collection_interfaces.f"
         

END INTERFACE

END MODULE ppm_module_interfaces

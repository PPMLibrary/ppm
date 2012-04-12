!minclude ppm_header(ppm_module_interfaces)

#define __REAL 3 
#define __COMPLEX 4 
#define __INTEGER 5 
#define __LONGINT 6 
#define __LOGICAL 7 
#define __CHAR 8 

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
! Global parameters
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

!PPM internal parameters used only to access entries in the
!particle's property data structure.
INTEGER,PARAMETER,PUBLIC   :: ppm_ops_inc_ghosts = 1
INTEGER,PARAMETER,PUBLIC   :: ppm_ops_interp = 2
INTEGER,PARAMETER,PUBLIC   :: ppm_ops_iscomputed = 3
INTEGER,PARAMETER,PUBLIC   :: ppm_ops_isdefined = 4
INTEGER,PARAMETER,PUBLIC   :: ppm_ops_vector = 5
INTEGER,PARAMETER,PUBLIC   :: ppm_param_length_opsflags = 5


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

TYPE,ABSTRACT ::  ppm_t_mesh_discr_info_
    !!! (Contained inside a ppm_t_field, so relates to one specific
    !!!  field, denoted by fieldID)
    !!! Data structure containing info about the current discretization
    !!! of fieldID on a given Mesh
    !!! 
    !!! Contains pointers to the data and bookkeeping information
    !!! for a mesh on which fieldID has been discretized.

    INTEGER                                          :: meshID = 0
    !!! id of the mesh on which fieldID is discretized
    INTEGER                                          :: lda = 0
    !!! number of components (1 for scalar fields)
    INTEGER                                          :: p_idx = 0
    !!! Storage index for the subpatch_data object which contains the data where
    !!! fieldID has been discretized on this mesh.
    !!! (A mesh stores data for several fields. Each subpatch thus 
    !!!  contains a list of data arrays, corresponding to the different
    !!!  fields. The data_ptr points to the location of the data array
    !!!  corresponding to fieldID within this list)
    LOGICAL,DIMENSION(ppm_mdata_lflags)              :: flags = .FALSE.
    !!! Booleans used to track the state of this discretization.

    CONTAINS
    PROCEDURE(mesh_discr_info_create_), DEFERRED :: create
    PROCEDURE(mesh_discr_info_destroy_),DEFERRED :: destroy
END TYPE
! Container for mesh_discr_info
minclude define_abstract_collection_type(ppm_t_mesh_discr_info_)

TYPE,ABSTRACT :: ppm_t_field_
    !!! Data structure for fields 
    !!! A field represents a mathematical concept (e.g. velocity, or
    !!! vorticity) and links to its discretized representation on meshes 
    !!! and/or on particles.

    INTEGER                                         :: ID = 0
    !!! global identifier 
    CHARACTER(LEN=ppm_char)                         :: name
    !!! string description
    INTEGER                                         :: lda = 0
    !!! number of components (1 for scalar fields)
    !!!
    !!! pointers to arrays where the scalar-value properties are stored
    CLASS(ppm_c_mesh_discr_info_),POINTER           :: M => NULL()
    !!! Collection of pointers to the data and bookkeeping information
    !!! for each mesh on which this field has been discretized.
    ! CLASS(ppm_c_part_discr_info_),POINTER           :: P => NULL()
    !    !!! Collection of pointers to the data and bookkeeping information
    !    !!! for each particle set on which this field has been discretized.

    CONTAINS
    PROCEDURE(field_create_),       DEFERRED :: create
    PROCEDURE(field_destroy_),      DEFERRED :: destroy
    PROCEDURE(field_discretize_on_),DEFERRED :: discretize_on
    PROCEDURE(field_set_rel_),      DEFERRED :: set_rel
    PROCEDURE(field_map_ghost_push_),DEFERRED:: map_ghost_push
    PROCEDURE(field_map_ghost_pop_),DEFERRED:: map_ghost_pop
END TYPE ppm_t_field_
! Container for fields
minclude define_abstract_collection_type(ppm_t_field_)

!!----------------------------------------------------------------------
!! Patches (contains the actual data arrays for this field)
!!----------------------------------------------------------------------
TYPE,ABSTRACT :: ppm_t_subpatch_data_
    !!! pointers to arrays where the data are stored
    INTEGER                                                :: fieldID = 0
    !!! ID of the field that is discretized here
    INTEGER                                                :: datatype = 0
    !!! Data type of the data being discretized
    INTEGER, DIMENSION(:,:), POINTER                       :: data_2d_i => NULL()
    !!! if the data is 2d int
    INTEGER, DIMENSION(:,:,:), POINTER                     :: data_3d_i => NULL()
    !!! if the data is 3d int
    INTEGER, DIMENSION(:,:,:,:), POINTER                   :: data_4d_i => NULL()
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
    LOGICAL, DIMENSION(:,:), POINTER                       :: data_2d_l => NULL()
    !!! if the data is 2d logical
    LOGICAL, DIMENSION(:,:,:), POINTER                     :: data_3d_l => NULL()
    !!! if the data is 3d logical
    LOGICAL, DIMENSION(:,:,:,:), POINTER                   :: data_4d_l => NULL()
    !!! if the data is 4d logical

    CONTAINS
    PROCEDURE(subpatch_data_create_), DEFERRED :: create
    PROCEDURE(subpatch_data_destroy_),DEFERRED :: destroy 
END TYPE
! Container for lists of (pointers to) subpatch_data
minclude define_abstract_collection_type(ppm_t_subpatch_data_)

TYPE,ABSTRACT :: ppm_t_subpatch_
    !!! intersection of a user-defined patch and a subdomain
    INTEGER                       :: meshID = 0
    !!! ID of the mesh to which this subpatch belongs
    INTEGER, DIMENSION(:),POINTER :: istart   => NULL()
    !!! Lower-left coordinates
    INTEGER, DIMENSION(:),POINTER :: iend     => NULL()
    !!! Upper-right coordinates
    INTEGER, DIMENSION(:),POINTER :: nnodes   => NULL()
    !!! number of nodes in each direction
    INTEGER, DIMENSION(:),POINTER :: istart_g => NULL()
    !!! Lower-left coordinates of the patch associated with this subpatch
    !!! (used for mapping routines)
    INTEGER, DIMENSION(:),POINTER :: iend_g   => NULL()
    !!! Upper-right coordinates of the patch associated with this subpatch
    !!! (used for mapping routines)
    CLASS(ppm_c_subpatch_data_),POINTER :: subpatch_data => NULL()
    !!! container for the data arrays for each property discretized
    !!! on this mesh
    CONTAINS
    PROCEDURE(subpatch_create_), DEFERRED  :: create
    PROCEDURE(subpatch_destroy_),DEFERRED  :: destroy
    PROCEDURE(subpatch_get_field_2d_rd_), DEFERRED :: subpatch_get_field_2d_rd
    PROCEDURE(subpatch_get_field_3d_rd_), DEFERRED :: subpatch_get_field_3d_rd
    GENERIC :: get_field => subpatch_get_field_2d_rd,subpatch_get_field_3d_rd
    !PROCEDURE  :: get => subpatch_get
END TYPE
! Container for lists of (pointers to) subpatch
minclude define_abstract_collection_type(ppm_t_subpatch_)


TYPE,ABSTRACT :: ppm_t_A_subpatch_
    !!! Bookkeeping derived-type. Contains an array of pointers to subpatches
    INTEGER                                        :: patchid = 0
    !!! Id of the patch
    INTEGER                                        :: nsubpatch = 0
    !!! Number of subpatches for this patch
    TYPE(ppm_t_ptr_subpatch),DIMENSION(:),POINTER  :: subpatch => NULL()
    !!! Pointers to each subpatches for this patch
    CONTAINS
    PROCEDURE(subpatch_A_create_) ,DEFERRED  :: create
    PROCEDURE(subpatch_A_destroy_),DEFERRED  :: destroy
END TYPE
! Container for subpatch pointer arrays
minclude define_abstract_collection_type(ppm_t_A_subpatch_)


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

TYPE,ABSTRACT ::  ppm_t_field_info_
    !!! (Contained inside a ppm_t_equi_mesh, so relates to one specific
    !!!  field, denoted by meshID)
    !!! Data structure containing info about a given field currently 
    !!! discretized this Mesh
    !!! 
    !!! Contains pointers to the field itself as well as some
    !!!  bookkeeping information 

    INTEGER                                 :: fieldID = 0
    !!! pointer to a field that is discretized on this mesh

    CONTAINS
    PROCEDURE(field_info_create_), DEFERRED :: create
    PROCEDURE(field_info_destroy_),DEFERRED :: destroy
END TYPE
! Container for field_info
minclude define_abstract_collection_type(ppm_t_field_info_)

TYPE ppm_t_subpatch_ptr_array
    INTEGER                                       :: size = 0
    INTEGER                                       :: nsubpatch = 0
    TYPE(ppm_t_ptr_subpatch),DIMENSION(:),POINTER :: vec => NULL()
END TYPE


TYPE,ABSTRACT :: ppm_t_equi_mesh_
    !!! Type for equispaced cartesian meshes on subs
   
    INTEGER                           :: ID = 0
    !!! ID of the mesh in the belonging topology
    !!! It is the same as its index in the ppm_t_topo%mesh array
    INTEGER                           :: topoid = 0
    !!! ID of the topology for which this mesh is defined

    INTEGER, DIMENSION(:), POINTER    :: Nm    => NULL()
    !!! global number of mesh points in computational domain

    REAL(ppm_kind_double),DIMENSION(:),POINTER :: Offset => NULL()
    !!! Offset

    REAL(ppm_kind_double),DIMENSION(:),POINTER :: h => NULL()
    !!! mesh spacing

    INTEGER, DIMENSION(:), POINTER    :: ghostsize   => NULL()
    !!! size of the ghost layer, in number of mesh points

    !TODO : delete this guy once the transition to the new DS is done
    INTEGER, DIMENSION(:,:), POINTER    :: nnodes    => NULL()

    INTEGER, DIMENSION(:,:), POINTER    :: istart    => NULL()
    INTEGER, DIMENSION(:,:), POINTER    :: iend    => NULL()

    CLASS(ppm_c_subpatch_),POINTER            :: subpatch => NULL()
    !!! container for subdomains patches 

    INTEGER                                   :: npatch = 0
    !!! Number of patches (NOT subpatches)
    CLASS(ppm_c_A_subpatch_),POINTER          :: patch => NULL()
    !!! array of arrays of pointers to the subpatches for each patch

    TYPE(ppm_t_subpatch_ptr_array),DIMENSION(:),POINTER :: &
                                                 subpatch_by_sub => NULL()
    !!! pointers to the subpatches contained in each sub.

    CLASS(ppm_c_field_info_),POINTER          :: field_ptr => NULL()
    !!! Pointers to the fields that are currently discretized on this mesh

    CLASS(ppm_t_mesh_mapping_s_),POINTER      :: mapping_s => NULL()
    CLASS(ppm_t_mesh_mapping_d_),POINTER      :: mapping_d => NULL()


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
    INTEGER, DIMENSION(:,:), POINTER :: ghost_patchid   => NULL()
    !!! list of patches of ghost mesh blocks (globel sub number).
    !!! 1st index: patch ID
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
    INTEGER, DIMENSION(:,:),POINTER  :: ghost_recvpatchid => NULL()
    !!! list of patches (global indices) for ghost mesh blocks to be received,
    !!! i.e. being ghost on the local processor (globel sub number).
    !!! 1st index: patch ID
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
    PROCEDURE(equi_mesh_create_),         DEFERRED :: create
    PROCEDURE(equi_mesh_destroy_),        DEFERRED :: destroy
    PROCEDURE(equi_mesh_def_patch_),      DEFERRED :: def_patch
    PROCEDURE(equi_mesh_set_rel_),        DEFERRED :: set_rel
    PROCEDURE(equi_mesh_def_uniform_),    DEFERRED :: def_uniform
    PROCEDURE(equi_mesh_new_subpatch_data_ptr_),&
      &                                   DEFERRED :: new_subpatch_data_ptr 
    PROCEDURE(equi_mesh_list_of_fields_), DEFERRED :: list_of_fields
    PROCEDURE(equi_mesh_block_intersect_),DEFERRED :: block_intersect
    PROCEDURE(equi_mesh_map_ghost_init_), DEFERRED :: map_ghost_init
    PROCEDURE(equi_mesh_map_ghost_get_),  DEFERRED :: map_ghost_get
    PROCEDURE(equi_mesh_map_ghost_push_), DEFERRED :: map_ghost_push
    PROCEDURE(equi_mesh_map_ghost_pop_),  DEFERRED :: map_ghost_pop
    PROCEDURE(equi_mesh_map_send_),       DEFERRED :: map_send
END TYPE
! Container for meshes
minclude define_abstract_collection_type(ppm_t_equi_mesh_)


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


!----------------------------------------------------------------------
! Interfaces for collections type-bound procedures
!----------------------------------------------------------------------
minclude define_abstract_collection_interfaces(ppm_t_equi_mesh_)
minclude define_abstract_collection_interfaces(ppm_t_A_subpatch_)
minclude define_abstract_collection_interfaces(ppm_t_mesh_discr_info_)
minclude define_abstract_collection_interfaces(ppm_t_field_info_)
minclude define_abstract_collection_interfaces(ppm_t_field_)
minclude define_abstract_collection_interfaces(ppm_t_subpatch_data_)
minclude define_abstract_collection_interfaces(ppm_t_subpatch_)

END INTERFACE

END MODULE ppm_module_interfaces

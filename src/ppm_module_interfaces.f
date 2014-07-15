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
      USE ppm_module_options

      IMPLICIT NONE

      !----------------------------------------------------------------------
      ! Global parameters
      !----------------------------------------------------------------------
      !PPM internal parameters used only to access entries in the
      !mesh data structures.
      INTEGER, PARAMETER, PUBLIC :: ppm_mdata_ghosts            = 1
      INTEGER, PARAMETER, PUBLIC :: ppm_mdata_reqput            = 2
      INTEGER, PARAMETER, PUBLIC :: ppm_mdata_cartesian         = 3
      INTEGER, PARAMETER, PUBLIC :: ppm_param_length_mdataflags = 3

      !PPM internal parameters used only to access entries in the
      !mesh data structures.
      INTEGER, PARAMETER, PUBLIC :: ppm_mesh_ghosts            = 1
      INTEGER, PARAMETER, PUBLIC :: ppm_mesh_reqput            = 2
      INTEGER, PARAMETER, PUBLIC :: ppm_mesh_cartesian         = 3
      INTEGER, PARAMETER, PUBLIC :: ppm_param_length_meshflags = 3

      !PPM internal parameters used only to access entries in the
      !particle's property data structure.
      INTEGER, PARAMETER, PUBLIC :: ppm_mdat_ghosts            = 1
      INTEGER, PARAMETER, PUBLIC :: ppm_mdat_reqput            = 2
      INTEGER, PARAMETER, PUBLIC :: ppm_mdat_map_ghosts        = 3
      INTEGER, PARAMETER, PUBLIC :: ppm_param_length_mdatflags = 3

      !User parameters
      INTEGER, PARAMETER, PUBLIC :: ppm_param_part_init_cartesian = 1
      INTEGER, PARAMETER, PUBLIC :: ppm_param_part_init_random    = 2


      !PPM internal parameters used only to access entries in the
      !particle data structure.
      INTEGER, PARAMETER, PUBLIC :: ppm_part_ghosts            = 1
      ! consider adding (?)
      !INTEGER, PARAMETER, PUBLIC :: ppm_part_ghosts_exist = 2
      INTEGER, PARAMETER, PUBLIC :: ppm_part_partial           = 2
      INTEGER, PARAMETER, PUBLIC :: ppm_part_reqput            = 3
      INTEGER, PARAMETER, PUBLIC :: ppm_part_areinside         = 4
      INTEGER, PARAMETER, PUBLIC :: ppm_part_cartesian         = 5
      INTEGER, PARAMETER, PUBLIC :: ppm_part_neighlists        = 6
      INTEGER, PARAMETER, PUBLIC :: ppm_part_global_index      = 7
      INTEGER, PARAMETER, PUBLIC :: ppm_param_length_partflags = 7

      INTEGER, PARAMETER, PUBLIC :: ppm_pdata_lflags = 3

      !PPM internal parameters used only to access entries in the
      !particle's property data structure.
      INTEGER, PARAMETER, PUBLIC :: ppm_ppt_ghosts            = 1
      INTEGER, PARAMETER, PUBLIC :: ppm_ppt_partial           = 2
      INTEGER, PARAMETER, PUBLIC :: ppm_ppt_reqput            = 3
      INTEGER, PARAMETER, PUBLIC :: ppm_ppt_map_parts         = 4
      INTEGER, PARAMETER, PUBLIC :: ppm_ppt_map_ghosts        = 5
      INTEGER, PARAMETER, PUBLIC :: ppm_param_length_pptflags = 5

      !PPM internal parameters used only to access entries in the
      !particle's property data structure.
      INTEGER, PARAMETER, PUBLIC :: ppm_ops_inc_ghosts        = 1
      INTEGER, PARAMETER, PUBLIC :: ppm_ops_interp            = 2
      INTEGER, PARAMETER, PUBLIC :: ppm_ops_iscomputed        = 3
      INTEGER, PARAMETER, PUBLIC :: ppm_ops_vector            = 4
      INTEGER, PARAMETER, PUBLIC :: ppm_param_length_opsflags = 4

      !----------------------------------------------------------------------
      ! Global variables
      !----------------------------------------------------------------------
      INTEGER :: ppm_nb_meshes    = 0
      INTEGER :: ppm_nb_part_sets = 0
      !----------------------------------------------------------------------
      ! Module variables
      !----------------------------------------------------------------------
      INTEGER, PRIVATE, DIMENSION(3) :: ldc
      !!! Generic type for all main PPM types

      !----------------------------------------------------------------------
      ! Type declaration
      !----------------------------------------------------------------------


      TYPE, ABSTRACT :: ppm_t_main_abstr
          !!! Generic type for all main PPM types
      END TYPE ppm_t_main_abstr
minclude ppm_create_collection(main_abstr,main_abstr,generate="concrete",vec=true)

      TYPE, ABSTRACT, EXTENDS(ppm_t_main_abstr) :: ppm_t_discr_kind
          !!! Discretization kinds (Particles and Meshes)
          INTEGER                 :: ID = 0
          !!! ID of the mesh or particle set in the belonging topology
          CHARACTER(LEN=ppm_char) :: name
          !!! name (optional...)
      END TYPE ppm_t_discr_kind
minclude ppm_create_collection(discr_kind,discr_kind,generate="concrete",vec=true,def_ptr=true)


      TYPE, ABSTRACT, EXTENDS(ppm_t_main_abstr) :: ppm_t_discr_data
          !!! Data (discretized on either Particles or Meshes)
          INTEGER                                       :: data_type
          !!! Data type for this property
          !!! One of: ppm_type_int, ppm_type_real, ppm_type_comp, ppm_type_logical
          !!!
          CLASS(ppm_t_main_abstr), POINTER              :: field_ptr => NULL()
          !!! Pointer to the field for which this is a discretization
          CLASS(ppm_t_discr_kind), POINTER              :: discr     => NULL()
          !!! Pointer to the discretization to which this data belongs
          CHARACTER(LEN=ppm_char)                       :: name
          !!! Name for this property
          LOGICAL, DIMENSION(ppm_param_length_pptflags) :: flags     = .FALSE.
          !!! Logical flags (applicable to either particle data or mesh data or both)
          !!!    ppm_ppt_ghosts
          !!!          true if ghost values are up-to-date
          !!!    ppm_ppt_partial
          !!!          true if there is a one-to-one mapping with the particles
          !!!    ppm_ppt_reqput
          !!!    ppm_ppt_map_parts
          !!!          true if partial mappings are desired for this property (default)
          !!!          (if false, the array for this property is not reallocated when
          !!!           particles move to a different processor or when they are
          !!!           interpolated from one distribution to another)
          !!!    ppm_ppt_map_ghosts
          !!!          true if ghost mappings are desired for this property (default)
          INTEGER                                       :: lda       = 0
          !!! Leading dimension of the data array
          !!!
      CONTAINS
          PROCEDURE :: has_ghosts => discr_data_has_ghosts
      END TYPE ppm_t_discr_data
minclude ppm_create_collection(discr_data,discr_data,generate="concrete",vec=true,def_ptr=true)

      TYPE, ABSTRACT, EXTENDS(ppm_t_main_abstr) :: ppm_t_operator_
          !!! Generic differential operator
          !!! (It only contains semantic information on the operator)
          INTEGER, DIMENSION(:),               POINTER :: degree => NULL()
          !!! degree of each term in the linear combination of differential ops
          REAL(ppm_kind_double), DIMENSION(:), POINTER :: coeffs => NULL()
          !!! array where the coefficients in linear combinations of
          !!! differential ops are stored
          INTEGER                                      :: nterms = 0
          !!! number of terms
          CHARACTER(LEN=ppm_char)                      :: name
          !!! name of the vector-valued property

      CONTAINS
          PROCEDURE(operator_create_),        DEFERRED :: create
          PROCEDURE(operator_destroy_),       DEFERRED :: destroy
          PROCEDURE(operator_discretize_on_), DEFERRED :: discretize_on
      END TYPE ppm_t_operator_
minclude ppm_create_collection(operator_,operator_,generate="abstract")


      TYPE, ABSTRACT :: ppm_t_operator_discr_
          !!! discretized operator
          CLASS(ppm_t_discr_kind), POINTER              :: discr_src => NULL()
          !!! Pointer to the discretization (mesh or particles) that this operator
          !!! takes data from
          CLASS(ppm_t_discr_kind), POINTER              :: discr_to  => NULL()
          !!! Pointer to the discretization (mesh or particles) that this operator
          !!! returns values on
          INTEGER, DIMENSION(:),   POINTER              :: order     => NULL()
          !!! Order of approximation for each term of the differential operator
          CLASS(ppm_t_operator_),  POINTER              :: op_ptr    => NULL()

          LOGICAL, DIMENSION(ppm_param_length_opsflags) :: flags     = .FALSE.
          !!! logical flags
          !!!    ppm_ops_inc_ghosts
          !!!           true if the operator should be computed for ghost
          !!!           particles too.  Note that the resulting values
          !!!           will be wrong for the ghost particles
          !!!           that have some neighbours outside the ghost layers.
          !!!           Default is false.
          !!!    ppm_ops_interp
          !!!          true if the op interpolates data from one set of particles
          !!!    ppm_ops_iscomputed
          !!!          true if the operator has been computed and is uptodate
          !!!    ppm_ops_vector
          !!!          true if each term represents a component (ie the result
          !!!          of the operator should be a vector field, like for gradients)
          !!!          false if the components are added up (like for the divergence)

      CONTAINS
          PROCEDURE(operator_discr_destroy_), DEFERRED :: destroy
          !PROCEDURE(operator_discr_compute_), DEFERRED :: compute
      END TYPE ppm_t_operator_discr_
minclude ppm_create_collection(operator_discr_,operator_discr_,generate="abstract")

      TYPE, EXTENDS(ppm_t_operator_discr_) :: ppm_t_operator_discr
      CONTAINS
          PROCEDURE :: create  => operator_discr_create
          PROCEDURE :: destroy => operator_discr_destroy
          PROCEDURE :: compute => operator_discr_compute
      END TYPE ppm_t_operator_discr
minclude ppm_create_collection(operator_discr,operator_discr,generate="extend")


#define  DTYPE(a) a/**/_s
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#include "map/mapping_abstract_typedef.f"
#include "operator/operator_discr_abstract_typedef.f"
#include "part/particles_abstract_typedef.f"
#undef  DTYPE
#undef  MK
#undef  _MK


#define  DTYPE(a) a/**/_d
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "map/mapping_abstract_typedef.f"
#include "operator/operator_discr_abstract_typedef.f"
#include "part/particles_abstract_typedef.f"
#undef  DTYPE
#undef  MK
#undef  _MK

      TYPE, ABSTRACT ::  ppm_t_discr_info_
          !!! (Contained inside a ppm_t_field, so relates to one specific
          !!!  field, denoted by fieldID)
          !!! Data structure containing info about the current discretization
          !!! of fieldID on a given Mesh or particle set
          !!!
          !!! Contains pointers to the data and bookkeeping information
          !!! for a mesh on which fieldID has been discretized.

          INTEGER                                         :: discrID    = 0
          !!! id of the discretization kind on which fieldID is discretized
          CLASS(ppm_t_discr_kind), POINTER                :: discr_ptr  => NULL()
          !!! pointer to the mesh
          CLASS(ppm_t_discr_data), POINTER                :: discr_data => NULL()
          !!! pointer to the data
          INTEGER                                         :: lda        = 0
          !!! number of components (1 for scalar fields)
          INTEGER                                         :: p_idx      = 0
          !!! For meshes:
          !!!     Storage index for the subpatch_data object which contains the data where
          !!! fieldID has been discretized on this mesh.
          !!! (A mesh stores data for several fields. Each subpatch thus
          !!!  contains a list of data arrays, corresponding to the different
          !!!  fields. The data_ptr points to the location of the data array
          !!!  corresponding to fieldID within this list)
          !!! For particle sets:
          !!!     not used
          LOGICAL, DIMENSION(ppm_param_length_mdataflags) :: flags = .FALSE.
          !!! Booleans used to track the state of this discretization.

      CONTAINS
          PROCEDURE(discr_info_create_),  DEFERRED :: create
          PROCEDURE(discr_info_destroy_), DEFERRED :: destroy
      END TYPE ppm_t_discr_info_
minclude ppm_create_collection(discr_info_,discr_info_,generate="abstract")

      TYPE, ABSTRACT, EXTENDS(ppm_t_main_abstr) :: ppm_t_field_
          !!! Data structure for fields
          !!! A field represents a mathematical concept (e.g. velocity, or
          !!! vorticity) and links to its discretized representation on meshes
          !!! and/or on particles.

          INTEGER                           :: ID         = 0
          !!! global identifier
          CHARACTER(LEN=ppm_char)           :: name
          !!! string description
          INTEGER                           :: data_type  = 0
          !!! data type
          !!! One of:
          !!!     ppm_type_int
          !!!     ppm_type_longint
          !!!     ppm_type_real
          !!!     ppm_type_comp
          !!!     ppm_type_logical
          !!!
          INTEGER                           :: lda        = 0
          !!! number of components (1 for scalar fields)
          !!!
          !!! pointers to arrays where the scalar-value properties are stored
          CLASS(ppm_c_discr_info_), POINTER :: discr_info => NULL()
          !!! Collection of pointers to the data and bookkeeping information
          !!! for each mesh or particle set on which this field has been discretized.

      CONTAINS
          PROCEDURE(field_create_),            DEFERRED :: create
          PROCEDURE(field_destroy_),           DEFERRED :: destroy
          PROCEDURE(field_set_rel_discr_),     DEFERRED :: set_rel_discr

          PROCEDURE(field_map_ghost_push_),    DEFERRED :: map_ghost_push
          PROCEDURE(field_map_ghost_pop_),     DEFERRED :: map_ghost_pop
          PROCEDURE(field_map_push_),          DEFERRED :: map_push
          PROCEDURE(field_map_pop_),           DEFERRED :: map_pop

          PROCEDURE(field_is_discretized_on_), DEFERRED :: is_discretized_on
          PROCEDURE(field_discretize_on_),     DEFERRED :: discretize_on
          PROCEDURE(field_get_pid_),           DEFERRED :: get_pid
          PROCEDURE(field_get_discr_),         DEFERRED :: get_discr
      END TYPE ppm_t_field_
minclude ppm_create_collection(field_,field_,generate="abstract")
minclude ppm_create_collection(field_,field_,generate="abstract",vec=true,def_ptr=false)

      !!----------------------------------------------------------------------
      !! Patches (contains the actual data arrays for this field)
      !!----------------------------------------------------------------------
      TYPE, ABSTRACT :: ppm_t_subpatch_data_
          !!! pointers to arrays where the data are stored
          INTEGER                                               :: fieldID = 0
          !!! ID of the field that is discretized here
          CLASS(ppm_t_discr_data),                      POINTER :: discr_data => NULL()
          !!! Pointer to a data structure that contains all information about
          !!! the discretization of a given field on this mesh (it has pointers to all
          !!! subpatches, for example)
          INTEGER,                  DIMENSION(:,:),     POINTER :: data_2d_i => NULL()
          !!! if the data is 2d int
          INTEGER,                  DIMENSION(:,:,:),   POINTER :: data_3d_i => NULL()
          !!! if the data is 3d int
          INTEGER,                  DIMENSION(:,:,:,:), POINTER :: data_4d_i => NULL()
          !!! if the data is 4d int
          REAL(ppm_kind_single),    DIMENSION(:,:),     POINTER :: data_2d_rs => NULL()
          !!! if the data is 2d real
          REAL(ppm_kind_single),    DIMENSION(:,:,:),   POINTER :: data_3d_rs => NULL()
          !!! if the data is 3d real
          REAL(ppm_kind_single),    DIMENSION(:,:,:,:), POINTER :: data_4d_rs => NULL()
          !!! if the data is 4d real
          COMPLEX(ppm_kind_single), DIMENSION(:,:),     POINTER :: data_2d_cs => NULL()
          !!! if the data is 2d complex
          COMPLEX(ppm_kind_single), DIMENSION(:,:,:),   POINTER :: data_3d_cs => NULL()
          !!! if the data is 3d complex
          COMPLEX(ppm_kind_single), DIMENSION(:,:,:,:), POINTER :: data_4d_cs => NULL()
          !!! if the data is 4d complex
          REAL(ppm_kind_double),    DIMENSION(:,:),     POINTER :: data_2d_rd => NULL()
          !!! if the data is 2d real
          REAL(ppm_kind_double),    DIMENSION(:,:,:),   POINTER :: data_3d_rd => NULL()
          !!! if the data is 3d real
          REAL(ppm_kind_double),    DIMENSION(:,:,:,:), POINTER :: data_4d_rd => NULL()
          !!! if the data is 4d real
          COMPLEX(ppm_kind_double), DIMENSION(:,:),     POINTER :: data_2d_cd => NULL()
          !!! if the data is 2d complex
          COMPLEX(ppm_kind_double), DIMENSION(:,:,:),   POINTER :: data_3d_cd => NULL()
          !!! if the data is 3d complex
          COMPLEX(ppm_kind_double), DIMENSION(:,:,:,:), POINTER :: data_4d_cd => NULL()
          !!! if the data is 4d complex
          LOGICAL,                  DIMENSION(:,:),     POINTER :: data_2d_l => NULL()
          !!! if the data is 2d logical
          LOGICAL,                  DIMENSION(:,:,:),   POINTER :: data_3d_l => NULL()
          !!! if the data is 3d logical
          LOGICAL,                  DIMENSION(:,:,:,:), POINTER :: data_4d_l => NULL()
          !!! if the data is 4d logical

      CONTAINS
          PROCEDURE(subpatch_data_create_),  DEFERRED :: create
          PROCEDURE(subpatch_data_destroy_), DEFERRED :: destroy
      END TYPE ppm_t_subpatch_data_
minclude ppm_create_collection(subpatch_data_,subpatch_data_,generate="abstract")
minclude ppm_create_collection(subpatch_data_,subpatch_data_,generate="abstract",vec=true,def_ptr=false,generate="abstract")

      TYPE, ABSTRACT, EXTENDS(ppm_t_discr_data) :: ppm_t_mesh_discr_data_
          CLASS(ppm_v_subpatch_data_), POINTER  :: subpatch => NULL()
      CONTAINS
          PROCEDURE(mesh_discr_data_create_),  DEFERRED :: create
          PROCEDURE(mesh_discr_data_destroy_), DEFERRED :: destroy
      END TYPE ppm_t_mesh_discr_data_
minclude ppm_create_collection(mesh_discr_data_,mesh_discr_data_,generate="abstract")
minclude ppm_create_collection(mesh_discr_data_,mesh_discr_data_,generate="abstract",vec=true,def_ptr=false)

      TYPE, ABSTRACT :: ppm_t_subpatch_
          !!! intersection of a user-defined patch and a subdomain
          INTEGER                                      :: meshID        =  0
          !!! ID of the mesh to which this subpatch belongs
          INTEGER                                      :: isub          =  0
          !!! subdomain that contains this subpatch (local id)
          CLASS(ppm_t_discr_kind),             POINTER :: mesh          => NULL()
          !!! Pointer to the mesh to which this subpatch belongs
          INTEGER,               DIMENSION(:), POINTER :: istart        => NULL()
          !!! Lower-left coordinates on the global mesh
          INTEGER,               DIMENSION(:), POINTER :: iend          => NULL()
          !!! Upper-right coordinates on the global mesh
          REAL(ppm_kind_double), DIMENSION(:), POINTER :: start         => NULL()
          !!! Lower-left absolute coordinates of the subpatch
          REAL(ppm_kind_double), DIMENSION(:), POINTER :: end           => NULL()
          !!! Upper-right absolute coordinates of the subpatch
          REAL(ppm_kind_double), DIMENSION(:), POINTER :: start_ext     => NULL()
          !!! Lower-left absolute coordinates of the subpatch extended with the
          !!! mesh-wide ghostsize. This corresponds to the "zone of influence" of
          !!! the subpatch (i.e., particles within that volume will interact with
          !!! the nodes of the subpatch during interpolation)
          REAL(ppm_kind_double), DIMENSION(:), POINTER :: end_ext       => NULL()
          !!! Upper-right absolute coordinates of the subpatch extended with the
          !!! mesh-wide ghostsize. This corresponds to the "zone of influence" of
          !!! the subpatch (i.e., particles within that volume will interact with
          !!! the nodes of the subpatch during interpolation)
          REAL(ppm_kind_double), DIMENSION(:), POINTER :: start_red     => NULL()
          !!! Lower-left absolute coordinates of the subpatch reduced by the
          !!! mesh-wide ghostsize. This corresponds to the area where
          !!! the field is well described by centered interpolation kernels
          !!! (i.e., we can used centred m2p interpolation kernels for particles
          !!! located inside that region)
          REAL(ppm_kind_double), DIMENSION(:), POINTER :: end_red       => NULL()
          !!! Upper-right absolute coordinates of the subpatch reduced by the
          !!! mesh-wide ghostsize. This corresponds to the area where
          !!! the field is well described by centered interpolation kernels
          !!! (i.e., we can used centred m2p interpolation kernels for particles
          !!! located inside that region)
          INTEGER,               DIMENSION(:), POINTER :: bc            => NULL()
          !!! boundary conditions on a subpatch:
          !!!
          !!! - west  : 1
          !!! - east  : 2
          !!! - south : 3
          !!! - north : 4
          !!! - bottom: 5
          !!! - top   : 6
          !!!
          !!! index 1: the index of the 4 or 6 faces in 2 and 3 D
          !!!
          !!! states:
          !!!
          !!! - value: 0 the face is internal
          !!! - value: 1 otherwise
          INTEGER,               DIMENSION(:), POINTER :: nnodes        => NULL()
          !!! number of (real) nodes in each direction
          INTEGER,               DIMENSION(:), POINTER :: lo_a          => NULL()
          !!! lbound of the allocated data array (subpatch extended with its ghostlayer)
          INTEGER,               DIMENSION(:), POINTER :: hi_a          => NULL()
          !!! ubound of the allocated data array (subpatch extended with its ghostlayer)
          INTEGER,               DIMENSION(:), POINTER :: istart_p      => NULL()
          !!! Lower-left coordinates of the patch associated with this subpatch
          !!! (used for mapping routines)
          INTEGER,               DIMENSION(:), POINTER :: iend_p        => NULL()
          !!! Upper-right coordinates of the patch associated with this subpatch
          !!! (used for mapping routines)
          INTEGER,               DIMENSION(:), POINTER :: ghostsize     => NULL()
          !!! Extent of the ghost layer in all 4 (resp. 6) directions in
          !!! 2D (resp. 3D). It is measured in number of mesh nodes.
          !!! - west  : 1
          !!! - east  : 2
          !!! - south : 3
          !!! - north : 4
          !!! - bottom: 5
          !!! - top   : 6
          CLASS(ppm_c_subpatch_data_),         POINTER :: subpatch_data => NULL()
          !!! container for the data arrays for each property discretized
          !!! on this mesh
      CONTAINS
          PROCEDURE(subpatch_create_),    DEFERRED  :: create
          PROCEDURE(subpatch_destroy_),   DEFERRED  :: destroy
          PROCEDURE(subpatch_get_pos2d_), DEFERRED  :: get_pos2d
          PROCEDURE(subpatch_get_pos3d_), DEFERRED  :: get_pos3d
          GENERIC :: get_pos => get_pos2d, get_pos3d

          PROCEDURE(subpatch_get_field_2d_i_),  DEFERRED :: subpatch_get_field_2d_i
          PROCEDURE(subpatch_get_field_3d_i_),  DEFERRED :: subpatch_get_field_3d_i
          PROCEDURE(subpatch_get_field_4d_i_),  DEFERRED :: subpatch_get_field_4d_i
          PROCEDURE(subpatch_get_field_2d_rs_), DEFERRED :: subpatch_get_field_2d_rs
          PROCEDURE(subpatch_get_field_3d_rs_), DEFERRED :: subpatch_get_field_3d_rs
          PROCEDURE(subpatch_get_field_4d_rs_), DEFERRED :: subpatch_get_field_4d_rs
          PROCEDURE(subpatch_get_field_2d_rd_), DEFERRED :: subpatch_get_field_2d_rd
          PROCEDURE(subpatch_get_field_3d_rd_), DEFERRED :: subpatch_get_field_3d_rd
          PROCEDURE(subpatch_get_field_4d_rd_), DEFERRED :: subpatch_get_field_4d_rd
          PROCEDURE(subpatch_get_field_2d_cs_), DEFERRED :: subpatch_get_field_2d_cs
          PROCEDURE(subpatch_get_field_3d_cs_), DEFERRED :: subpatch_get_field_3d_cs
          PROCEDURE(subpatch_get_field_4d_cs_), DEFERRED :: subpatch_get_field_4d_cs
          PROCEDURE(subpatch_get_field_2d_cd_), DEFERRED :: subpatch_get_field_2d_cd
          PROCEDURE(subpatch_get_field_3d_cd_), DEFERRED :: subpatch_get_field_3d_cd
          PROCEDURE(subpatch_get_field_4d_cd_), DEFERRED :: subpatch_get_field_4d_cd
          PROCEDURE(subpatch_get_field_2d_l_),  DEFERRED :: subpatch_get_field_2d_l
          PROCEDURE(subpatch_get_field_3d_l_),  DEFERRED :: subpatch_get_field_3d_l
          PROCEDURE(subpatch_get_field_4d_l_),  DEFERRED :: subpatch_get_field_4d_l

          GENERIC :: get_field =>              &
          &          subpatch_get_field_2d_i,  &
          &          subpatch_get_field_3d_i,  &
          &          subpatch_get_field_4d_i,  &
          &          subpatch_get_field_2d_rs, &
          &          subpatch_get_field_3d_rs, &
          &          subpatch_get_field_4d_rs, &
          &          subpatch_get_field_2d_rd, &
          &          subpatch_get_field_3d_rd, &
          &          subpatch_get_field_4d_rd, &
          &          subpatch_get_field_2d_cs, &
          &          subpatch_get_field_3d_cs, &
          &          subpatch_get_field_4d_cs, &
          &          subpatch_get_field_2d_cd, &
          &          subpatch_get_field_3d_cd, &
          &          subpatch_get_field_4d_cd, &
          &          subpatch_get_field_2d_l,  &
          &          subpatch_get_field_3d_l,  &
          &          subpatch_get_field_4d_l

          !PROCEDURE  :: get => subpatch_get
      END TYPE ppm_t_subpatch_
minclude ppm_create_collection(subpatch_,subpatch_,generate="abstract")

      TYPE, ABSTRACT :: ppm_t_A_subpatch_
          !!! Bookkeeping derived-type. Contains an array of pointers to subpatches
          INTEGER                                        :: patchid   = 0
          !!! Id of the patch
          INTEGER                                        :: nsubpatch = 0
          !!! Number of subpatches for this patch
          TYPE(ppm_t_ptr_subpatch),DIMENSION(:), POINTER :: subpatch  => NULL()
          !!! Pointers to each subpatches for this patch
      CONTAINS
          PROCEDURE(subpatch_A_create_) , DEFERRED  :: create
          PROCEDURE(subpatch_A_destroy_), DEFERRED  :: destroy
      END TYPE ppm_t_A_subpatch_
minclude ppm_create_collection(A_subpatch_,A_subpatch_,generate="abstract")

      TYPE ppm_t_mesh_maplist
          !!! TODO: check what this is used for (imported from Petros code
          INTEGER,                 POINTER :: target_topoid => NULL()
          !!! target topology ID
          INTEGER,                 POINTER :: target_meshid => NULL()
          !!! target mesh ID
          INTEGER,                 POINTER :: nsendlist     => NULL()
          !!! send rank list size
          INTEGER,                 POINTER :: nrecvlist     => NULL()
          !!! recv rank size
          INTEGER, DIMENSION(:),   POINTER :: isendlist     => NULL()
          !!! send rank lists
          INTEGER, DIMENSION(:),   POINTER :: irecvlist     => NULL()
          !!! recv rank lists
          INTEGER, DIMENSION(:),   POINTER :: isendfromsub  => NULL()
          !!! source send sublist
          INTEGER, DIMENSION(:,:), POINTER :: isendblkstart => NULL()
          !!! send block start list
          INTEGER, DIMENSION(:,:), POINTER :: isendblksize  => NULL()
          !!! send block size list
          INTEGER, DIMENSION(:),   POINTER :: irecvtosub    => NULL()
          !!! recv send sublist
          INTEGER, DIMENSION(:,:), POINTER :: irecvblkstart => NULL()
          !!! recv block start list
          INTEGER, DIMENSION(:,:), POINTER :: irecvblksize  => NULL()
          !!! recv block size list
          INTEGER, DIMENSION(:),   POINTER :: psendbuffer   => NULL()
          !!! send buffer pointer
          INTEGER, DIMENSION(:),   POINTER :: precvbuffer   => NULL()
          !!! recv buffer pointer
      END TYPE ppm_t_mesh_maplist

      TYPE ppm_t_subpatch_ptr_array
          INTEGER                                         :: size      = 0
          INTEGER                                         :: nsubpatch = 0
          TYPE(ppm_t_ptr_subpatch), DIMENSION(:), POINTER :: vec       => NULL()
      END TYPE ppm_t_subpatch_ptr_array


      TYPE, ABSTRACT, EXTENDS(ppm_t_discr_kind) :: ppm_t_equi_mesh_
          !!! Type for equispaced cartesian meshes on subs

          INTEGER                                        :: topoid            = 0
          !!! ID of the topology for which this mesh is defined

          INTEGER,               DIMENSION(:),   POINTER :: Nm                => NULL()
          !!! global number of mesh points in computational domain

          REAL(ppm_kind_double), DIMENSION(:),   POINTER :: Offset            => NULL()
          !!! Offset

          REAL(ppm_kind_double), DIMENSION(:),   POINTER :: h                 => NULL()
          !!! mesh spacing

          INTEGER,               DIMENSION(:),   POINTER :: ghostsize         => NULL()
          !!! size of the ghost layer, in number of mesh points

          !TODO : delete this guy once the transition to the new DS is done
          INTEGER,               DIMENSION(:,:), POINTER :: nnodes            => NULL()

          INTEGER,               DIMENSION(:,:), POINTER :: istart            => NULL()
          INTEGER,               DIMENSION(:,:), POINTER :: iend              => NULL()

          CLASS(ppm_c_subpatch_),                POINTER :: subpatch          => NULL()
          !!! container for subdomains patches (these contain all the data
          !!! discretized on this mesh)

          CLASS(ppm_v_mesh_discr_data_),         POINTER :: mdata             => NULL()
          !!! pointers to the discretized data on this mesh, stored by field.
          !!! (for each field corresponds a ppm_t_mesh_discr_data object, which
          !!! contains pointers to the data stored on each subpatch)

          INTEGER                                        :: npatch            = 0
          !!! Number of patches (NOT subpatches)
          CLASS(ppm_c_A_subpatch_),              POINTER :: patch             => NULL()
          !!! array of arrays of pointers to the subpatches for each patch

          TYPE(ppm_t_subpatch_ptr_array), &
          &                      DIMENSION(:),   POINTER :: subpatch_by_sub   => NULL()
          !!! pointers to the subpatches contained in each sub.

          CLASS(ppm_v_main_abstr),               POINTER :: field_ptr         => NULL()
          !!! Pointers to the fields that are currently discretized on this mesh

          CLASS(ppm_t_mesh_mapping_s_),          POINTER :: mapping_s         => NULL()
          CLASS(ppm_t_mesh_mapping_d_),          POINTER :: mapping_d         => NULL()


          !------------------------------------------------------------------
          !  Mesh ghosts mappings
          !------------------------------------------------------------------
          LOGICAL                                        :: ghost_initialized = .FALSE.
          !!! is .TRUE. if the ghost mappings have been initialized
          !!! else, .FALSE.
          INTEGER,               DIMENSION(:),   POINTER :: ghost_fromsub     => NULL()
          !!! list of source subs of ghost mesh blocks (globel sub number).
          !!! These are the owner subs of the actual real mesh points
          !!! 1st index: meshblock ID
          INTEGER,               DIMENSION(:),   POINTER :: ghost_tosub       => NULL()
          !!! list of target subs of ghost mesh blocks (globel sub number).
          !!! These are the subs a block will serve as a ghost on.
          !!! 1st index: meshblock ID
          INTEGER,               DIMENSION(:,:), POINTER :: ghost_patchid     => NULL()
          !!! list of patches of ghost mesh blocks (globel sub number).
          !!! 1st index: patch ID
          INTEGER,               DIMENSION(:,:), POINTER :: ghost_blkstart    => NULL()
          !!! start (lower-left corner) of ghost mesh block in GLOBAL
          !!! mesh coordinates. First index: x,y[,z], 2nd: meshblock ID
          INTEGER,               DIMENSION(:,:), POINTER :: ghost_blksize     => NULL()
          !!! size (in grid points) of ghost blocks. 1st index: x,y[,z], 2nd:
          !!! meshblock ID
          INTEGER,               DIMENSION(:),   POINTER :: ghost_blk         => NULL()
          !!! mesh ghost block list. 1st index: target processor
          INTEGER                                        :: ghost_nsend
          !!! number of mesh blocks to be sent as ghosts
          INTEGER                                        :: ghost_nrecv
          !!! number of mesh blocks to be recvd as ghosts
          INTEGER,               DIMENSION(:),   POINTER :: ghost_recvtosub   => NULL()
          !!! list of target subs for ghost mesh blocks to be received,
          !!! i.e. being ghost on the local processor (globel sub number).
          !!! These are the subs where the blocks will serve as ghosts
          !!! 1st index: meshblock ID
          INTEGER,               DIMENSION(:,:), POINTER :: ghost_recvpatchid => NULL()
          !!! list of patches (global indices) for ghost mesh blocks to be received,
          !!! i.e. being ghost on the local processor (globel sub number).
          !!! 1st index: patch ID
          INTEGER,               DIMENSION(:,:), POINTER :: ghost_recvblkstart=> NULL()
          !!! start (lower-left corner) of received ghost mesh block in
          !!! GLOBAL  mesh coordinates. 1st index: x,y[,z], 2nd: meshblock ID
          INTEGER,               DIMENSION(:,:), POINTER :: ghost_recvblksize => NULL()
          !!! size (in grid points) of recvd ghost blocks.
          !!! 1st index: x,y[,z], 2nd: meshblock ID
          INTEGER,               DIMENSION(:),   POINTER :: ghost_recvblk     => NULL()
          !!! mesh ghost block receive list. 1st index: target processor

          TYPE(ppm_t_mesh_maplist),              POINTER :: mapping           => NULL()

      CONTAINS
          PROCEDURE(equi_mesh_create_),         DEFERRED :: create
          PROCEDURE(equi_mesh_destroy_),        DEFERRED :: destroy
          PROCEDURE(equi_mesh_create_prop_),    DEFERRED :: create_prop
          PROCEDURE(equi_mesh_prop_zero_),      DEFERRED :: zero
          PROCEDURE(equi_mesh_def_patch_),      DEFERRED :: def_patch
          PROCEDURE(equi_mesh_def_uniform_),    DEFERRED :: def_uniform
          PROCEDURE(equi_mesh_new_subpatch_data_ptr_),&
          &                                     DEFERRED :: new_subpatch_data_ptr
          PROCEDURE(equi_mesh_list_of_fields_), DEFERRED :: list_of_fields
          PROCEDURE(equi_mesh_block_intersect_),DEFERRED :: block_intersect
          PROCEDURE(equi_mesh_map_ghost_init_), DEFERRED :: map_ghost_init
          PROCEDURE(equi_mesh_map_ghost_get_),  DEFERRED :: map_ghost_get
          PROCEDURE(equi_mesh_map_ghost_put_),  DEFERRED :: map_ghost_put
          PROCEDURE(equi_mesh_map_push_),       DEFERRED :: map_ghost_push
          PROCEDURE(equi_mesh_map_pop_),        DEFERRED :: map_ghost_pop
          PROCEDURE(equi_mesh_map_send_),       DEFERRED :: map_send
          PROCEDURE(equi_mesh_map_isend_),      DEFERRED :: map_isend

          PROCEDURE(equi_mesh_map_global_),     DEFERRED :: map
          PROCEDURE(equi_mesh_map_push_),       DEFERRED :: map_push
          PROCEDURE(equi_mesh_map_pop_),        DEFERRED :: map_pop

          PROCEDURE(equi_mesh_print_vtk_),      DEFERRED :: print_vtk
          PROCEDURE(equi_mesh_m2p_),            DEFERRED :: interp_to_part
      END TYPE ppm_t_equi_mesh_
minclude ppm_create_collection(equi_mesh_,equi_mesh_,generate="abstract")

      TYPE, EXTENDS(ppm_t_main_abstr) :: ppm_t_var_discr_pair
        CLASS(ppm_t_main_abstr), POINTER :: var   => NULL()
        CLASS(ppm_t_discr_kind), POINTER :: discr => NULL()
      END TYPE ppm_t_var_discr_pair
minclude ppm_create_collection(var_discr_pair,var_discr_pair,vec=true,generate="concrete")

      TYPE, EXTENDS(ppm_t_main_abstr) :: ppm_t_field_discr_pair
        CLASS(ppm_t_field_),     POINTER :: field => NULL()
        CLASS(ppm_t_discr_kind), POINTER :: discr => NULL()
      END TYPE ppm_t_field_discr_pair
minclude ppm_create_collection(field_discr_pair,field_discr_pair,vec=true,generate="concrete")

      !----------------------------------------------------------------------
      !  INTERFACES
      !----------------------------------------------------------------------
      INTERFACE
      !----------------------------------------------------------------------
      ! Interfaces for collections type-bound procedures
      !----------------------------------------------------------------------
minclude ppm_create_collection_interfaces(equi_mesh_,equi_mesh_)
minclude ppm_create_collection_interfaces(A_subpatch_,A_subpatch_)
minclude ppm_create_collection_interfaces(discr_info_,discr_info_)
minclude ppm_create_collection_interfaces(field_,field_)
minclude ppm_create_collection_interfaces(field_,field_,vec=true)
minclude ppm_create_collection_interfaces(operator_,operator_)
minclude ppm_create_collection_interfaces(operator_discr_,operator_discr_)
minclude ppm_create_collection_interfaces(subpatch_data_,subpatch_data_)
minclude ppm_create_collection_interfaces(subpatch_data_,subpatch_data_,vec=true)
minclude ppm_create_collection_interfaces(subpatch_,subpatch_)
minclude ppm_create_collection_interfaces(mesh_discr_data_,mesh_discr_data_)
minclude ppm_create_collection_interfaces(mesh_discr_data_,mesh_discr_data_,vec=true)
      !minclude ppm_create_collection_interfaces(ppm_t_discr_kind_,vec=true)

#define  DTYPE(a) a/**/_s
#include "map/mapping_interfaces.f"

#define  DTYPE(a) a/**/_d
#include "map/mapping_interfaces.f"

#define  DTYPE(a) a/**/_s
#define  MK ppm_kind_single
#include "operator/dcop_interfaces.f"
#include "part/particles_interfaces.f"
#undef  DTYPE
#undef  MK

#define  DTYPE(a) a/**/_d
#define  MK ppm_kind_double
#include "operator/dcop_interfaces.f"
#include "part/particles_interfaces.f"
#undef  DTYPE
#undef  MK

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

#include "operator/operator_interfaces.f"

      END INTERFACE

      CONTAINS

minclude ppm_create_collection_procedures(operator_discr,operator_discr_)
minclude ppm_create_collection_procedures(discr_kind,discr_kind,vec=true)
minclude ppm_create_collection_procedures(discr_data,discr_data,vec=true)
minclude ppm_create_collection_procedures(main_abstr,main_abstr,vec=true)
minclude ppm_create_collection_procedures(var_discr_pair,var_discr_pair,vec=true)
minclude ppm_create_collection_procedures(field_discr_pair,field_discr_pair,vec=true)
minclude ppm_create_collection_procedures(part_prop_s_,part_prop_s_,vec=true)
minclude ppm_create_collection_procedures(part_prop_d_,part_prop_d_,vec=true)

      !CREATE (DUMMY ROUTINE)
      SUBROUTINE operator_discr_create(this,Op,Part_src,Part_to,info,&
      &          with_ghosts,vector,interp,order,order_v,prop)
          CLASS(ppm_t_operator_discr)                               :: this
          CLASS(ppm_t_operator_),            TARGET,  INTENT(IN   ) :: Op
          CLASS(ppm_t_discr_kind),           TARGET,  INTENT(IN   ) :: Part_src
          CLASS(ppm_t_discr_kind),           TARGET,  INTENT(IN   ) :: Part_to
          INTEGER,                                    INTENT(  OUT) :: info
          LOGICAL,                 OPTIONAL,          INTENT(IN   ) :: with_ghosts
          LOGICAL,                 OPTIONAL,          INTENT(IN   ) :: vector
          LOGICAL,                 OPTIONAL,          INTENT(IN   ) :: interp
          INTEGER,                 OPTIONAL,          INTENT(IN   ) :: order
          INTEGER, DIMENSION(:),   OPTIONAL, POINTER, INTENT(IN   ) :: order_v
          CLASS(ppm_t_discr_data), OPTIONAL, TARGET                 :: prop
          start_subroutine("operator_discr_create")
          fail("this dummy routine should not be called")
          end_subroutine()
      END SUBROUTINE operator_discr_create
      !DESTROY (DUMMY ROUTINE)
      SUBROUTINE operator_discr_destroy(this,info)
          CLASS(ppm_t_operator_discr)               :: this
          INTEGER,                    INTENT(  OUT) :: info
          start_subroutine("operator_discr_destroy")
          fail("this dummy routine should not be called")
          end_subroutine()
      END SUBROUTINE operator_discr_destroy
      !COMPUTE (DUMMY ROUTINE)
      SUBROUTINE operator_discr_compute(this,Field_src,Field_to,info)
          CLASS(ppm_t_operator_discr)                    :: this
          CLASS(ppm_t_main_abstr), TARGET, INTENT(IN   ) :: Field_src
          CLASS(ppm_t_main_abstr), TARGET, INTENT(INOUT) :: Field_to
          INTEGER,                         INTENT(  OUT) :: info
          start_subroutine("operator_discr_compute")
          fail("this dummy routine should not be called")
          end_subroutine()
      END SUBROUTINE operator_discr_compute

      FUNCTION discr_data_has_ghosts(this) RESULT(res)
          CLASS (ppm_t_discr_data)  :: this
          LOGICAL                   :: res
          res = this%flags(ppm_ppt_ghosts)
      END FUNCTION discr_data_has_ghosts

      END MODULE ppm_module_interfaces

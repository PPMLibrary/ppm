!minclude ppm_header(ppm_module_field_typedef)

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
INTEGER                            :: ppm_nb_fields = 0

!----------------------------------------------------------------------
! Module variables 
!----------------------------------------------------------------------
INTEGER, PRIVATE, DIMENSION(3)  :: ldc

!----------------------------------------------------------------------
! Type declaration
!----------------------------------------------------------------------
TYPE, EXTENDS(ppm_t_mesh_discr_info_):: ppm_t_mesh_discr_info
    CONTAINS
    PROCEDURE :: create => mesh_discr_info_create
    PROCEDURE :: destroy => mesh_discr_info_destroy

END TYPE ppm_t_mesh_discr_info
minclude define_collection_type(ppm_t_mesh_discr_info)

TYPE,EXTENDS(ppm_t_field_) :: ppm_t_field
    CONTAINS
    PROCEDURE :: create => field_create
    PROCEDURE :: destroy => field_destroy
    PROCEDURE :: discretize_on => field_discretize_on
    PROCEDURE :: set_rel => field_set_rel
    PROCEDURE :: map_ghost_push => field_map_ghost_push
    PROCEDURE :: map_ghost_pop  => field_map_ghost_pop
END TYPE ppm_t_field
minclude define_collection_type(ppm_t_field)



!----------------------------------------------------------------------
!  Type-bound procedures
!----------------------------------------------------------------------
CONTAINS

minclude define_collection_procedures(ppm_t_mesh_discr_info)
minclude define_collection_procedures(ppm_t_field)

!CREATE
SUBROUTINE field_create(this,lda,name,info,init_func)
    !!! Constructor for subdomain data data structure
    CLASS(ppm_t_field)                      :: this
    INTEGER,                     INTENT(IN) :: lda
    !!! number of components
    CHARACTER(LEN=*),            INTENT(IN) :: name
    INTEGER,                    INTENT(OUT) :: info
    REAL(ppm_kind_double),EXTERNAL,POINTER,OPTIONAL,INTENT(IN) :: init_func
    !!! support for initialisation function not finished (need to think
    !!! about data types...)

    start_subroutine("field_create")

    ppm_nb_fields = ppm_nb_fields + 1
    this%ID = ppm_nb_fields
    this%name = TRIM(ADJUSTL(name))
    this%lda = lda
    IF (ASSOCIATED(this%M)) THEN
        fail("Seems like this field was alrady allocated - Call destroy() first?")
    ENDIF
    ALLOCATE(ppm_c_mesh_discr_info::this%M,STAT=info)
    or_fail_alloc("could not allocate this%M with ppm_c_mesh_discr_info type")

    end_subroutine()
END SUBROUTINE field_create
!DESTROY
SUBROUTINE field_destroy(this,info)
    !!! Destructor for subdomain data data structure
    CLASS(ppm_t_field)                 :: this
    INTEGER,               INTENT(OUT) :: info

    start_subroutine("field_destroy")

    this%ID = 0
    this%name = ''
    this%lda = 0

    !Destroy the bookkeeping entries in the fields that are
    !discretized on this mesh
    !TODO !!!

    destroy_collection_ptr(this%M)

    end_subroutine()
END SUBROUTINE field_destroy
!CREATE
SUBROUTINE mesh_discr_info_create(this,meshID,lda,p_idx,flags,info)
    !!! Constructor for subdomain data data structure
    CLASS(ppm_t_mesh_discr_info)             :: this
    INTEGER,                      INTENT(IN) :: meshID
    INTEGER,                      INTENT(IN) :: lda
    !!! number of components
    INTEGER,                      INTENT(IN) :: p_idx
    !!! index in the subpatch_data arrays where the data is stored
    LOGICAL,DIMENSION(ppm_mdata_lflags)     :: flags
    INTEGER,                    INTENT(OUT) :: info

    start_subroutine("mesh_discr_info_create")

    this%meshID = meshID
    this%lda    = lda
    this%p_idx  = p_idx
    this%flags  = flags

    end_subroutine()
END SUBROUTINE mesh_discr_info_create
!DESTROY
SUBROUTINE mesh_discr_info_destroy(this,info)
    !!! Destructor for subdomain data data structure
    CLASS(ppm_t_mesh_discr_info)             :: this
    INTEGER,               INTENT(OUT) :: info

    start_subroutine("mesh_discr_info_destroy")
    
    this%meshID = 0
    this%lda = 0
    this%p_idx = 0
    this%flags = .FALSE.

    end_subroutine()
END SUBROUTINE mesh_discr_info_destroy


SUBROUTINE field_discretize_on(this,mesh,info,datatype)
    !!! Allocate field on a mesh
    !!! If the field has a procedure for initialization (e.g. an
    !!! initial condition), then the field is also initialized.
    !!! If the mesh has patches, the field is allocated only on these
    !!! patches. If no patches have been defined, it is assumed that
    !!! that the user expects the field to be allocated on the whole domain.
    !!! A single patch is then defined, covering all the subdomains.
    CLASS(ppm_t_field)                 :: this
    CLASS(ppm_t_equi_mesh_)            :: mesh
    !!! mesh onto which this field is to be discretized
    INTEGER,               INTENT(OUT)  :: info
    INTEGER, OPTIONAL                   :: datatype
    !!! By default, the type is assumed to be real, double-precision.

    CLASS(ppm_t_mesh_discr_info_),    POINTER :: mdinfo => NULL()
    CLASS(ppm_t_subpatch_),     POINTER :: p => NULL()
    CLASS(ppm_t_subpatch_data_),POINTER :: subpdat => NULL()
    INTEGER                             :: dtype,p_idx

    start_subroutine("field_discretize_on")

    IF (PRESENT(datatype)) THEN
        dtype = datatype
    ELSE
        dtype = ppm_type_real_double
    END IF

    !Check whether this field has already been initialized
    IF (this%ID.LE.0 .OR. this%lda.LE.0) THEN
        fail("Field needs to be initialized before calling discretized. Call ThisField%create() first")
    ENDIF

    !Check that the mesh contains patches onto which the data
    ! should be allocated. If not, create a single patch that
    ! covers the whole domain.
    IF (.NOT.ASSOCIATED(mesh%subpatch)) THEN
        fail("Mesh not allocated. Use mesh%create() first.")
    ELSE
        IF (mesh%subpatch%nb.LE.0) THEN
            CALL mesh%def_uniform(info)
                or_fail("failed to create a uniform patch data structure")
        ENDIF
    ENDIF

    !Create a new data array on the mesh to store this field
    p => mesh%subpatch%begin()
    DO WHILE (ASSOCIATED(p))
        ! create a new subpatch_data object
        !ALLOCATE(ppm_t_subpatch_data::subpdat,STAT=info)
        subpdat => mesh%new_subpatch_data_ptr(info)
            or_fail_alloc("could not get a new ppm_t_subpatch_data pointer")

        CALL subpdat%create(this%ID,dtype,this%lda,p%nnodes,info)
            or_fail("could not create new subpatch_data")

        IF (.NOT.ASSOCIATED(p%subpatch_data)) THEN
           fail("p%subpatch_data not allocated")
        ENDIF

        CALL p%subpatch_data%push(subpdat,info)
            or_fail("could not add new subpatch_data to subpatch collection")

        p => mesh%subpatch%next()
    ENDDO

    !on each subpatch, the new data is the last element of 
    !the subpatch_data collection
    p => mesh%subpatch%begin()
    IF (ASSOCIATED(p)) THEN
        p_idx = p%subpatch_data%max_id
    ELSE
        !this proc has no subpatches for this mesh
        p_idx = -1
    ENDIF


    !Update the bookkeeping table to store the relationship between
    ! the mesh and the field.
    CALL this%set_rel(mesh,p_idx,info)
        or_fail("failed to log the relationship between this field and that mesh")

    CALL mesh%set_rel(this,info)
        or_fail("failed to log the relationship between this mesh and that field")

    end_subroutine()
END SUBROUTINE field_discretize_on

SUBROUTINE field_set_rel(this,mesh,p_idx,info)
    !!! Create bookkeeping data structure to log the relationship between
    !!! the field and a mesh
    CLASS(ppm_t_field)                  :: this
    CLASS(ppm_t_equi_mesh_)             :: mesh
    !!! mesh that this field is discretized on
    INTEGER,                INTENT(IN ) :: p_idx
    !!! index in the mesh data structure where the data for this field is stored
    INTEGER,                INTENT(OUT) :: info

    TYPE(ppm_t_mesh_discr_info),POINTER :: mdinfo => NULL()
    LOGICAL,DIMENSION(ppm_mdata_lflags) :: flags

    start_subroutine("field_set_rel")

    flags = .FALSE.

    ALLOCATE(mdinfo,STAT=info)
        or_fail_alloc("could not allocate new mesh_discr_info pointer")

    CALL mdinfo%create(mesh%ID,this%lda,p_idx,flags,info)
        or_fail("could not create new mesh_discr_info object")

    IF (.NOT.ASSOCIATED(this%M)) THEN
        ALLOCATE(ppm_c_mesh_discr_info::this%M,STAT=info)
            or_fail_alloc("could not allocate mesh_discr_info collection")
    ENDIF

    IF (this%M%vec_size.LT.mesh%ID) THEN
        CALL this%M%grow_size(info)
            or_fail("could not grow_size this%M")
    ENDIF
    !if the collection is still too small, it means we should consider
    ! using a hash table for meshIDs...
    IF (this%M%vec_size.LT.mesh%ID) THEN
        fail("The id of the mesh that we are trying to store in a Field collection seems to large. Why not implement a hash function?")
    ENDIF

    IF (this%M%exists(mesh%ID)) THEN
        fail("It seems like this Mesh had already been discretized on this field. Are you sure you want to do that a second time?")
    ELSE
        !add mesh_discr_info to the collection, at index meshID
        !Update the counters nb,max_id and min_id accordingly
        !(this is not very clean without hash table)
        this%M%vec(mesh%ID)%t => mdinfo
        this%M%nb = this%M%nb + 1
        IF (this%M%nb.EQ.1) THEN
            this%M%max_id = mesh%ID
            this%M%min_id = mesh%ID
        ELSE
            this%M%max_id = MAX(this%M%max_id,mesh%ID)
            this%M%min_id = MIN(this%M%min_id,mesh%ID)
        ENDIF
    ENDIF

    end_subroutine()
END SUBROUTINE field_set_rel

SUBROUTINE field_map_ghost_push(this,mesh,info)
    !!! Push field data into the buffers of a mesh for ghost mappings
    !!! The field must of course be stored on this mesh
    !!! (for now) we assume that mesh%map_ghost_get() has already been called
    CLASS(ppm_t_field)                  :: this
    CLASS(ppm_t_equi_mesh_)             :: mesh
    !!! mesh that this field is discretized on
    INTEGER,                INTENT(OUT) :: info

    start_subroutine("field_map_ghost_push")


    CALL mesh%map_ghost_push(this,info)
        or_fail("mesh%map_ghost_push")

    end_subroutine()
END SUBROUTINE field_map_ghost_push

SUBROUTINE field_map_ghost_pop(this,mesh,info)
    !!! Pop field data from the buffers of a mesh for ghost mappings
    !!! The field must of course be stored on this mesh
    CLASS(ppm_t_field)                  :: this
    CLASS(ppm_t_equi_mesh_)             :: mesh
    !!! mesh that this field is discretized on
    INTEGER,                INTENT(OUT) :: info

    start_subroutine("field_map_ghost_pop")


    CALL mesh%map_ghost_pop(this,info)
        or_fail("mesh%map_ghost_pop")

    end_subroutine()
END SUBROUTINE field_map_ghost_pop

END MODULE ppm_module_field_typedef

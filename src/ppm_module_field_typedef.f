ppm_header(ppm_module_field_typedef)

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
define_collection_type(ppm_t_mesh_data)

TYPE,EXTENDS(ppm_t_field_) :: ppm_t_field
    CONTAINS
    PROCEDURE :: create => field_create
    PROCEDURE :: destroy => field_destroy
    PROCEDURE :: discretize_on => field_discretize_on
END TYPE ppm_t_field
define_collection_type(ppm_t_field)



!----------------------------------------------------------------------
!  Type-bound procedures
!----------------------------------------------------------------------
CONTAINS

define_collection_procedures(ppm_t_mesh_data)
define_collection_procedures(ppm_t_field)

!CREATE
SUBROUTINE field_create(this,fieldID,lda,name,info)
    !!! Constructor for subdomain data data structure
    CLASS(ppm_t_field)                      :: this
    INTEGER,                     INTENT(IN) :: fieldID
    INTEGER,                     INTENT(IN) :: lda
    !!! number of components
    CHARACTER(LEN=*),            INTENT(IN) :: name
    INTEGER,                    INTENT(OUT) :: info


    REAL(KIND(1.D0))                   :: t0
    CHARACTER(LEN=ppm_char)            :: caller = 'field_create'


    CALL substart(caller,t0,info)


    this%fieldID = fieldID
    this%name = TRIM(ADJUSTL(name))
    this%lda = lda
    ALLOCATE(ppm_c_mesh_data::this%M,STAT=info)
    or_fail_alloc("could not allocate this%M with ppm_c_mesh_data type")

    CALL substop(caller,t0,info)

    9999  CONTINUE

END SUBROUTINE field_create
!DESTROY
SUBROUTINE field_destroy(this,info)
    !!! Destructor for subdomain data data structure
    CLASS(ppm_t_field)                 :: this
    INTEGER,               INTENT(OUT) :: info

    REAL(KIND(1.D0))                   :: t0
    CHARACTER(LEN=ppm_char)            :: caller = 'field_destroy'

    CALL substart(caller,t0,info)

    this%fieldID = 0
    this%name = ''
    this%lda = 0
    CALL this%M%destroy(info)

    CALL substop(caller,t0,info)
    9999  CONTINUE
END SUBROUTINE field_destroy
!CREATE
SUBROUTINE mesh_data_create(this,meshID,lda,p_idx,flags,info)
    !!! Constructor for subdomain data data structure
    CLASS(ppm_t_mesh_data)                      :: this
    INTEGER,                     INTENT(IN) :: meshID
    INTEGER,                     INTENT(IN) :: lda
    !!! number of components
    INTEGER,                     INTENT(IN) :: p_idx
    LOGICAL,DIMENSION(ppm_mdata_lflags)     :: flags
    INTEGER,                    INTENT(OUT) :: info


    REAL(KIND(1.D0))                   :: t0
    CHARACTER(LEN=ppm_char)            :: caller = 'mesh_data_create'


    CALL substart(caller,t0,info)


    this%meshID = meshID
    this%lda = lda
    this%p_idx = p_idx
    this%flags = flags


    CALL substop(caller,t0,info)

    9999  CONTINUE

END SUBROUTINE mesh_data_create
!DESTROY
SUBROUTINE mesh_data_destroy(this,info)
    !!! Destructor for subdomain data data structure
    CLASS(ppm_t_mesh_data)             :: this
    INTEGER,               INTENT(OUT) :: info

    REAL(KIND(1.D0))                   :: t0
    CHARACTER(LEN=ppm_char)            :: caller = 'mesh_data_destroy'

    CALL substart(caller,t0,info)

    this%meshID = 0
    this%lda = 0
    this%p_idx = 0
    this%flags = .FALSE.

    CALL substop(caller,t0,info)
    9999  CONTINUE
END SUBROUTINE mesh_data_destroy


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

    REAL(KIND(1.D0))                    :: t0
    CHARACTER(LEN=ppm_char)             :: caller = 'field_discretize_on'
    CLASS(ppm_t_mesh_data_),    POINTER :: mdat => NULL()
    CLASS(ppm_t_subpatch_),     POINTER :: p => NULL()
    CLASS(ppm_t_subpatch_data_),POINTER :: subpdat => NULL()
    INTEGER                             :: dtype,p_idx
    INTEGER,DIMENSION(:),       POINTER :: Nmp => NULL()
    LOGICAL,DIMENSION(ppm_mdata_lflags) :: flags

    CALL substart(caller,t0,info)

    IF (PRESENT(datatype)) THEN
        dtype = datatype
    ELSE
        dtype = ppm_type_real_double
    END IF

    ALLOCATE(Nmp(ppm_dim),STAT=info)
    or_fail_alloc("could not allocate Nmp")

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
        or_fail_alloc("could not allocate ppm_t_subpatch_data pointer")

        !Nmp(1:ppm_dim) = p%iend(1:ppm_dim) - p%istart(1:ppm_dim)
        CALL subpdat%create(dtype,this%lda,p%nnodes,info)
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
    flags = .FALSE.

    ALLOCATE(ppm_t_mesh_data::mdat,STAT=info)
    or_fail_alloc("could not allocate new mesh_data pointer")

    CALL mdat%create(mesh%ID,this%lda,p_idx,flags,info)
    or_fail("could not create new mesh_data object")

    IF (.NOT.ASSOCIATED(this%M)) THEN
        ALLOCATE(ppm_c_mesh_data::this%M,STAT=info)
        or_fail_alloc("could not allocate mesh_data collection")
    ENDIF
    IF (.NOT.this%M%exists(mesh%ID)) THEN
        CALL this%M%push(mdat,info)
        or_fail("could not add new mesh_data to M collection")
    ELSE
        SELECT TYPE(md => mdat)
        TYPE IS (ppm_t_mesh_data)
            SELECT TYPE(v => this%M%vec(mesh%ID)%t)
            TYPE IS (ppm_t_mesh_data)
                v = md
            CLASS DEFAULT
                fail("wrong type - implementation error")
            END SELECT
        CLASS DEFAULT
            fail("wrong type - implementation error")
        END SELECt
    ENDIF

    DEALLOCATE(Nmp,STAT=info)
    or_fail_dealloc("could not deallocate Nmp")

    CALL substop(caller,t0,info)
    9999  CONTINUE
END SUBROUTINE field_discretize_on


END MODULE ppm_module_field_typedef

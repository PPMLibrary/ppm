ppm_header(ppm_module_mesh_typedef)

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
INTEGER                            :: ppm_nb_meshes = 0
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
define_collection_type(ppm_t_subpatch_data)

TYPE,EXTENDS(ppm_t_subpatch_) :: ppm_t_subpatch
    CONTAINS
    PROCEDURE  :: create    => subpatch_create
    PROCEDURE  :: destroy   => subpatch_destroy

    PROCEDURE  :: subpatch_get_field_2d_rd
    PROCEDURE  :: subpatch_get_field_3d_rd
!    GENERIC    :: get_field => subpatch_get_field_3d_rd,&
!        &                   subpatch_get_field_2d_rd
    !PROCEDURE  :: get => subpatch_get
END TYPE
define_collection_type(ppm_t_subpatch)


TYPE,EXTENDS(ppm_t_A_subpatch_) :: ppm_t_A_subpatch
    CONTAINS
    PROCEDURE :: create  => subpatch_A_create
    PROCEDURE :: destroy => subpatch_A_destroy
END TYPE
define_collection_type(ppm_t_A_subpatch)


TYPE,EXTENDS(ppm_t_field_info_) ::  ppm_t_field_info
    CONTAINS
    PROCEDURE :: create  => field_info_create
    PROCEDURE :: destroy => field_info_destroy
END TYPE
define_collection_type(ppm_t_field_info)

TYPE,EXTENDS(ppm_t_equi_mesh_) :: ppm_t_equi_mesh
    CONTAINS
    PROCEDURE  :: create    => equi_mesh_create
    PROCEDURE  :: destroy   => equi_mesh_destroy
    PROCEDURE  :: def_patch => equi_mesh_def_patch
    PROCEDURE  :: def_uniform           => equi_mesh_def_uniform
    PROCEDURE  :: new_subpatch_data_ptr => equi_mesh_new_subpatch_data_ptr
    PROCEDURE  :: list_of_fields        => equi_mesh_list_of_fields
    PROCEDURE  :: set_rel   => equi_mesh_set_rel
END TYPE
define_collection_type(ppm_t_equi_mesh)

!----------------------------------------------------------------------
! DATA STORAGE for the meshes
!----------------------------------------------------------------------
TYPE(ppm_c_equi_mesh)              :: ppm_mesh



!----------------------------------------------------------------------
!  Type-bound procedures
!----------------------------------------------------------------------
CONTAINS

!Procedures for collections of derived types
define_collection_procedures(ppm_t_equi_mesh)
define_collection_procedures(ppm_t_subpatch_data)
define_collection_procedures(ppm_t_field_info)
define_collection_procedures(ppm_t_subpatch)
define_collection_procedures(ppm_t_A_subpatch)

SUBROUTINE subpatch_get_field_3d_rd(this,wp,Field,info)
    !!! Constructor for subdomain data data structure
    CLASS(ppm_t_subpatch)                :: this
    CLASS(ppm_t_field_)                  :: Field
    REAL(ppm_kind_double),DIMENSION(:,:,:),POINTER :: wp
    INTEGER,                 INTENT(OUT) :: info

    start_subroutine("subpatch_data_create")

    !Direct access to the data arrays 
    wp => this%subpatch_data%vec(Field%M%vec(this%meshID)%t%p_idx)%t%data_3d_rd
    IF (.NOT.ASSOCIATED(wp)) THEN
        IF (Field%lda+ppm_dim.NE.3) THEN
            fail("wrong dimensions for pointer arg wp")
        ENDIF
    ENDIF

    end_subroutine()
END SUBROUTINE

SUBROUTINE subpatch_get_field_2d_rd(this,wp,Field,info)
    !!! Constructor for subdomain data data structure
    CLASS(ppm_t_subpatch)                :: this
    CLASS(ppm_t_field_)                  :: Field
    REAL(ppm_kind_double),DIMENSION(:,:),POINTER :: wp
    INTEGER,                 INTENT(OUT) :: info

    INTEGER                            :: iopt, ndim

    start_subroutine("subpatch_data_create")

    !Direct access to the data arrays 
    wp => this%subpatch_data%vec(Field%M%vec(this%meshID)%t%p_idx)%t%data_2d_rd
    IF (.NOT.ASSOCIATED(wp)) THEN
        IF (Field%lda+ppm_dim.NE.2) THEN
            fail("wrong dimensions for pointer arg wp")
        ENDIF
    ENDIF

    end_subroutine()
END SUBROUTINE

SUBROUTINE subpatch_data_create(this,fieldID,datatype,lda,Nmp,info)
    !!! Constructor for subdomain data data structure
    CLASS(ppm_t_subpatch_data)              :: this
    INTEGER,                     INTENT(IN) :: fieldID
    !!! ID of the field that is discretized on this mesh patch
    INTEGER,                     INTENT(IN) :: datatype
    !!! data type of the data
    INTEGER,                     INTENT(IN) :: lda
    !!! number of data components per mesh node
    INTEGER,DIMENSION(:),POINTER,INTENT(IN) :: Nmp
    !!! number of mesh nodes in each dimension on this patch
    INTEGER,                    INTENT(OUT) :: info

    INTEGER                            :: iopt, ndim

    start_subroutine("subpatch_data_create")

    this%fieldID = fieldID
    this%datatype = datatype

    iopt   = ppm_param_alloc_grow

    ! Determine allocation size of the data array
    IF (MINVAL(Nmp(1:ppm_dim)) .LE. 0) THEN
        or_fail("invalid size for patch data. This patch should be deleted")
    ENDIF

    ldc(1:ppm_dim) = Nmp(1:ppm_dim)

    ndim = ppm_dim
    IF (lda.GE.2) THEN
        ndim = ndim +1
        ldc(ndim) = lda
    ENDIF

    SELECT CASE (ndim)
    CASE (2)
        SELECT CASE (datatype)
        CASE (ppm_type_int)
            CALL ppm_alloc(this%data_2d_i,ldc,iopt,info)
        CASE (ppm_type_real_single)
            CALL ppm_alloc(this%data_2d_rs,ldc,iopt,info)
        CASE (ppm_type_real_double)
            CALL ppm_alloc(this%data_2d_rd,ldc,iopt,info)
        CASE (ppm_type_comp_single)
            CALL ppm_alloc(this%data_2d_cs,ldc,iopt,info)
        CASE (ppm_type_comp_double)
            CALL ppm_alloc(this%data_2d_cd,ldc,iopt,info)
        CASE (ppm_type_logical)
            CALL ppm_alloc(this%data_2d_l,ldc,iopt,info)
        CASE DEFAULT
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,caller,   &
                &        'invalid type for mesh patch data',__LINE__,info)
        END SELECT
    CASE (3)
        SELECT CASE (datatype)
        CASE (ppm_type_int)
            CALL ppm_alloc(this%data_3d_i,ldc,iopt,info)
        CASE (ppm_type_real_single)
            CALL ppm_alloc(this%data_3d_rs,ldc,iopt,info)
        CASE (ppm_type_real_double)
            CALL ppm_alloc(this%data_3d_rd,ldc,iopt,info)
        CASE (ppm_type_comp_single)
            CALL ppm_alloc(this%data_3d_cs,ldc,iopt,info)
        CASE (ppm_type_comp_double)
            CALL ppm_alloc(this%data_3d_cd,ldc,iopt,info)
        CASE (ppm_type_logical)
            CALL ppm_alloc(this%data_3d_l,ldc,iopt,info)
        CASE DEFAULT
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,caller,   &
                &        'invalid type for mesh patch data',__LINE__,info)
        END SELECT
    CASE (4)
        SELECT CASE (datatype)
        CASE (ppm_type_int)
            CALL ppm_alloc(this%data_4d_i,ldc,iopt,info)
        CASE (ppm_type_real_single)
            CALL ppm_alloc(this%data_4d_rs,ldc,iopt,info)
        CASE (ppm_type_real_double)
            CALL ppm_alloc(this%data_4d_rd,ldc,iopt,info)
        CASE (ppm_type_comp_single)
            CALL ppm_alloc(this%data_4d_cs,ldc,iopt,info)
        CASE (ppm_type_comp_double)
            CALL ppm_alloc(this%data_4d_cd,ldc,iopt,info)
        CASE (ppm_type_logical)
            CALL ppm_alloc(this%data_4d_l,ldc,iopt,info)
        CASE DEFAULT
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,caller,   &
                &        'invalid type for mesh patch data',__LINE__,info)
        END SELECT

    END SELECT

    or_fail_alloc('allocating mesh patch data failed')

    end_subroutine()
END SUBROUTINE subpatch_data_create
!DESTROY
SUBROUTINE subpatch_data_destroy(this,info)
    !!! Destructor for subdomain data data structure
    CLASS(ppm_t_subpatch_data)         :: this
    INTEGER,               INTENT(OUT) :: info

    INTEGER                            :: iopt

    start_subroutine("patch_data_destroy")

    this%fieldID = 0
    this%datatype = 0

    iopt = ppm_param_dealloc
    CALL ppm_alloc(this%data_2d_i,ldc,iopt,info)
    CALL ppm_alloc(this%data_3d_i,ldc,iopt,info)
    CALL ppm_alloc(this%data_4d_i,ldc,iopt,info)

    CALL ppm_alloc(this%data_2d_l,ldc,iopt,info)
    CALL ppm_alloc(this%data_3d_l,ldc,iopt,info)
    CALL ppm_alloc(this%data_4d_l,ldc,iopt,info)

    CALL ppm_alloc(this%data_2d_rs,ldc,iopt,info)
    CALL ppm_alloc(this%data_3d_rs,ldc,iopt,info)
    CALL ppm_alloc(this%data_4d_rs,ldc,iopt,info)
    CALL ppm_alloc(this%data_2d_rd,ldc,iopt,info)
    CALL ppm_alloc(this%data_3d_rd,ldc,iopt,info)
    CALL ppm_alloc(this%data_4d_rd,ldc,iopt,info)

    CALL ppm_alloc(this%data_2d_cs,ldc,iopt,info)
    CALL ppm_alloc(this%data_3d_cs,ldc,iopt,info)
    CALL ppm_alloc(this%data_4d_cs,ldc,iopt,info)
    CALL ppm_alloc(this%data_2d_cd,ldc,iopt,info)
    CALL ppm_alloc(this%data_3d_cd,ldc,iopt,info)
    CALL ppm_alloc(this%data_4d_cd,ldc,iopt,info)

    end_subroutine()
END SUBROUTINE subpatch_data_destroy

!CREATE
SUBROUTINE subpatch_create(p,meshID,istart,iend,info)
    !!! Destructor for subdomain data data structure
    CLASS(ppm_t_subpatch)              :: p
    INTEGER                            :: meshID
    INTEGER,DIMENSION(:)               :: istart
    INTEGER,DIMENSION(:)               :: iend
    INTEGER,               INTENT(OUT) :: info

    INTEGER                            :: iopt

    start_subroutine("subpatch_create")

    IF (.NOT.ASSOCIATED(p%istart)) THEN
        ALLOCATE(p%istart(ppm_dim),STAT=info)
            or_fail_alloc("could not allocate p%istart")
    ENDIF
    IF (.NOT.ASSOCIATED(p%iend)) THEN
        ALLOCATE(p%iend(ppm_dim),STAT=info)
            or_fail_alloc("could not allocate p%iend")
    ENDIF
    IF (.NOT.ASSOCIATED(p%nnodes)) THEN
        ALLOCATE(p%nnodes(ppm_dim),STAT=info)
            or_fail_alloc("could not allocate p%nnodes")
    ENDIF

    p%meshID = meshID
    p%istart = istart
    p%iend   = iend
    p%nnodes(1:ppm_dim) = 1 + iend(1:ppm_dim) - istart(1:ppm_dim)
    IF (.NOT.ASSOCIATED(p%subpatch_data)) THEN
        ALLOCATE(ppm_c_subpatch_data::p%subpatch_data,STAT=info)
            or_fail_alloc("could not allocate p%subpatch_data")
    ENDIF

    end_subroutine()
END SUBROUTINE subpatch_create

!DESTROY
SUBROUTINE subpatch_destroy(p,info)
    !!! Destructor for subdomain data data structure
    CLASS(ppm_t_subpatch)              :: p
    INTEGER,               INTENT(OUT) :: info

    INTEGER                            :: iopt

    start_subroutine("subpatch_destroy")

    iopt = ppm_param_dealloc
    CALL ppm_alloc(p%nnodes,ldc,iopt,info)
        or_fail_dealloc("p%nnodes")
    CALL ppm_alloc(p%istart,ldc,iopt,info)
        or_fail_dealloc("p%istart")
    CALL ppm_alloc(p%iend,ldc,iopt,info)
        or_fail_dealloc("p%iend")

    destroy_collection_ptr(p%subpatch_data)

    end_subroutine()
END SUBROUTINE subpatch_destroy

!CREATE
SUBROUTINE subpatch_A_create(this,vecsize,info,patchid)
    !!! Destructor for subdomain data data structure
    CLASS(ppm_t_A_subpatch)            :: this
    INTEGER                            :: vecsize
    INTEGER,               INTENT(OUT) :: info
    INTEGER,OPTIONAL,      INTENT(IN)  :: patchid

    start_subroutine("subpatch_A_create")

    IF (ASSOCIATED(this%subpatch)) THEN
        CALL this%destroy(info)
        or_fail_dealloc("could not destroy this ppm_t_A_subpatch object")
    ENDIF
    IF (.NOT.ASSOCIATED(this%subpatch)) THEN
        ALLOCATE(ppm_t_ptr_subpatch::this%subpatch(vecsize),STAT=info)
        or_fail_alloc("could not allocate this%subpatch")
    ENDIF

    this%nsubpatch = 0
    IF (PRESENT(patchid)) THEN
        this%patchid = patchid
    ELSE
        this%patchid = 0
    ENDIF

    end_subroutine()
END SUBROUTINE subpatch_A_create

!DESTROY
SUBROUTINE subpatch_A_destroy(this,info)
    !!! Destructor for subdomain data data structure
    CLASS(ppm_t_A_subpatch)            :: this
    INTEGER,               INTENT(OUT) :: info

    start_subroutine("subpatch_A_destroy")

    this%patchid = 0
    this%nsubpatch = 0
    IF (ASSOCIATED(this%subpatch)) THEN
        DEALLOCATE(this%subpatch,STAT=info)
        or_fail_dealloc("could not deallocate this%subpatch")
    ENDIF
    NULLIFY(this%subpatch)

    end_subroutine()
END SUBROUTINE subpatch_A_destroy


SUBROUTINE equi_mesh_def_patch(this,patch,info,patchid)
    !!! Add a patch to a mesh
    USE ppm_module_topo_typedef

    CLASS(ppm_t_equi_mesh)                  :: this
    REAL(ppm_kind_double),DIMENSION(:)      :: patch
    !!! Positions of the corners of the patch
    !!! (x1,y1,z1,x2,y2,z2), where 1 is the lower-left-bottom corner
    !!! and 2 is the upper-right-top corner.
    INTEGER                 , INTENT(  OUT) :: info
    !!! Returns status, 0 upon success
    INTEGER, OPTIONAL                       :: patchid
    !!! id of the patch, if we want one.
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    TYPE(ppm_t_topo), POINTER :: topo => NULL()
    INTEGER                   :: i,isub,nsubpatch,id,pid,meshID
    CLASS(ppm_t_subpatch_),  POINTER :: p => NULL()
    CLASS(ppm_t_A_subpatch_),POINTER :: A_p => NULL()
    INTEGER,              DIMENSION(ppm_dim) :: istart,iend
    REAL(ppm_kind_double),DIMENSION(ppm_dim) :: h,Offset

    start_subroutine("mesh_def_patch")

    !-------------------------------------------------------------------------
    !  Check arguments
    !-------------------------------------------------------------------------
    IF (.NOT.ASSOCIATED(this%subpatch)) THEN
        fail("Mesh not allocated. Call Mesh%create() first?")
    ENDIF
    !-------------------------------------------------------------------------
    !  get mesh parameters
    !-------------------------------------------------------------------------
    h = this%h
    Offset = this%Offset
    meshID = this%ID
    topo => ppm_topo(this%topoid)%t
    !-------------------------------------------------------------------------
    !  Allocate bookkeeping arrays (pointers between patches and subpatches)
    !-------------------------------------------------------------------------

    nsubpatch = 0
    ALLOCATE(ppm_t_A_subpatch::A_p,STAT=info)
    or_fail_alloc("could not allocate ppm_t_A_subpatch pointer")

    SELECT TYPE(t => A_p)
    TYPE IS (ppm_t_A_subpatch)
        CALL t%create(topo%nsublist,info,patchid)
        or_fail("could not initialize ppm_t_A_subpatch pointer")
    END SELECT

    IF (PRESENT(patchid)) THEN
        pid = patchid
    ELSE
        pid = 0
    ENDIF
    !-------------------------------------------------------------------------
    !  intersect each subdomain with the patch and store the corresponding
    !  subpatch in the mesh data structure
    !-------------------------------------------------------------------------
    ! loop through all subdomains on this processor
    DO i = 1,topo%nsublist
        isub = topo%isublist(i)
        ! check if the subdomain overlaps with the patch
        IF (ALL(patch(1:ppm_dim).LT.topo%max_subd(1:ppm_dim,isub)) .AND. &
            ALL(patch(ppm_dim+1:2*ppm_dim).GT.topo%min_subd(1:ppm_dim,isub)))&
               THEN 

            ! finds the mesh nodes which are contained in the overlap region
            istart(1:ppm_dim) = 1 + CEILING((MAX(patch(1:ppm_dim),&
                topo%min_subd(1:ppm_dim,isub))-Offset(1:ppm_dim))/h(1:ppm_dim))
            iend(1:ppm_dim)= 1 + FLOOR((MIN(patch(ppm_dim+1:2*ppm_dim),&
                topo%max_subd(1:ppm_dim,isub))-Offset(1:ppm_dim))/h(1:ppm_dim))

            ! create a new subpatch object
            ALLOCATE(ppm_t_subpatch::p,STAT=info)
                or_fail_alloc("could not allocate ppm_t_subpatch pointer")

            CALL p%create(meshID,istart,iend,info)
                or_fail("could not create new subpatch")

            nsubpatch = nsubpatch+1
            A_p%subpatch(nsubpatch)%t => p

            ! add it to the list of subpatches on this mesh
            CALL this%subpatch%push(p,info,id)
                or_fail("could not add new subpatch to mesh")

            !TODO keeping track of the list of subpatches for each subdomain
            !CALL this%sub%push(info,id)
                !or_fail("could not add new subpatch_ptr_array to mesh%sub")
            !this%sub%vec(id)%t%nsubpatch = nsubpatch
            !this%sub%vec(id)%t%patchid = pid

        ENDIF
    ENDDO
    CALL this%patch%push(A_p,info,id)
        or_fail("could not add new subpatch_ptr_array to mesh%patch")
    this%patch%vec(id)%t%nsubpatch = nsubpatch
    this%patch%vec(id)%t%patchid = pid

    end_subroutine()
END SUBROUTINE equi_mesh_def_patch

SUBROUTINE equi_mesh_def_uniform(this,info,patchid)
    !!! Add a uniform patch to a mesh
    CLASS(ppm_t_equi_mesh)                  :: this
    INTEGER                 , INTENT(  OUT) :: info
    !!! Returns status, 0 upon success
    INTEGER, OPTIONAL                       :: patchid
    !!! id of the (uniform) patch, if we want one.
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    REAL(ppm_kind_double),DIMENSION(2*ppm_dim) :: patch

    start_subroutine("mesh_def_uniform")

    !create a huge patch
    patch(1:ppm_dim)           = -HUGE(1._ppm_kind_double)
    patch(ppm_dim+1:2*ppm_dim) =  HUGE(1._ppm_kind_double)

    !and add it to the mesh (it will compute the intersection
    ! between this infinite patch and the (hopefully) finite
    ! computational domain)
    CALL this%def_patch(patch,info,patchid)
        or_fail("failed to add patch")

    !TODO add some checks for the finiteness of the computational domain

    end_subroutine()
END SUBROUTINE equi_mesh_def_uniform

SUBROUTINE equi_mesh_create(this,topoid,Offset,info,Nm,h)
    !!! Constructor for the cartesian mesh object
    !-------------------------------------------------------------------------
    !  Modules
    !-------------------------------------------------------------------------
    USE ppm_module_topo_typedef
    USE ppm_module_check_id
    IMPLICIT NONE
    !-------------------------------------------------------------------------
    !  Includes
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    CLASS(ppm_t_equi_mesh)                  :: this
    !!! cartesian mesh object
    INTEGER                 , INTENT(IN   ) :: topoid
    !!! Topology ID for which mesh has been created 
    REAL(ppm_kind_double), DIMENSION(:), INTENT(IN   ) :: Offset
    !!! Offset in each dimension
    INTEGER                 , INTENT(  OUT) :: info
    !!! Returns status, 0 upon success
    INTEGER,DIMENSION(:),              OPTIONAL,INTENT(IN   ) :: Nm
    !!! Global number of mesh points in the whole comput. domain
    !!! Makes sense only if the computational domain is bounded.
    !!! Note: Exactly one of Nm and h should be specified
    REAL(ppm_kind_double),DIMENSION(:),OPTIONAL,INTENT(IN   ) :: h
    !!! Mesh spacing
    !!! Note: Exactly one of Nm and h should be specified
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER , DIMENSION(3)    :: ldc
    INTEGER                   :: iopt,ld,ud,kk,i,j,isub
    LOGICAL                   :: valid
    TYPE(ppm_t_topo), POINTER :: topo => NULL()

    start_subroutine("equi_mesh_create")

    !-------------------------------------------------------------------------
    !  Check arguments
    !-------------------------------------------------------------------------
    IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
    ENDIF

    !This mesh is defined for a given topology
    this%topoid = topoid

    !dumb way of creating a global ID for this mesh
    !TODO find something better? (needed if one creates and destroy
    ! many meshes)
    ppm_nb_meshes = ppm_nb_meshes + 1
    this%ID = ppm_nb_meshes 

    topo => ppm_topo(topoid)%t

    !macro_test( 1 )
    !macro_test( "1" )
    !macro_test( "(1)" )
    !macro_test( 1+2 )
    !macro_test( (1) )
    !macro_test( 1/(1+2) )

    !check_equal(a,b,"this is an error message")

    !check_equal(SIZE(Nm,1),ppm_dim,"invalid size for Nm")

    !-------------------------------------------------------------------------
    !  (Re)allocate memory for the internal mesh list and Arrays at meshid
    !-------------------------------------------------------------------------
    iopt   = ppm_param_alloc_fit

    ldc(1) = ppm_dim
    CALL ppm_alloc(this%Nm,ldc,iopt,info)
    or_fail_alloc('Nm')

    CALL ppm_alloc(this%Offset,ldc,iopt,info)
    or_fail_alloc('Offset')

    CALL ppm_alloc(this%h,ldc,iopt,info)
    or_fail_alloc('h')

    IF (.NOT.ASSOCIATED(this%subpatch)) THEN
        ALLOCATE(ppm_c_subpatch::this%subpatch,STAT=info)
        or_fail_alloc("could not allocate this%subpatch")
    ELSE
        fail("subpatch collection is already allocated. Call destroy() first?")
    ENDIF

    IF (.NOT.ASSOCIATED(this%patch)) THEN
        ALLOCATE(ppm_c_A_subpatch::this%patch,STAT=info)
        or_fail_alloc("could not allocate this%patch")
    ELSE
        fail("patch collection is already allocated. Call destroy() first?")
    ENDIF

    IF (.NOT.ASSOCIATED(this%sub)) THEN
        ALLOCATE(ppm_c_A_subpatch::this%sub,STAT=info)
        or_fail_alloc("could not allocate this%sub")
    ELSE
        fail("sub collection is already allocated. Call destroy() first?")
    ENDIF

    !-------------------------------------------------------------------------
    !  Store the mesh information
    !-------------------------------------------------------------------------
    this%Offset(1:ppm_dim) = Offset(1:ppm_dim)

    IF (PRESENT(h)) THEN
        this%h(1:ppm_dim) = h(1:ppm_dim)
    ELSE
        this%Nm(1:ppm_dim) = Nm(1:ppm_dim)
        IF (ASSOCIATED(topo%max_physs)) THEN
            this%h(1:ppm_dim) = (topo%max_physs(1:ppm_dim) - &
                topo%min_physs(1:ppm_dim))/(Nm(1:ppm_dim)-1)
        ELSE
            this%h(1:ppm_dim) = (topo%max_physd(1:ppm_dim) - &
                topo%min_physd(1:ppm_dim))/(Nm(1:ppm_dim)-1)
        ENDIF
    ENDIF
    


    !-------------------------------------------------------------------------
    !  Return
    !-------------------------------------------------------------------------
    end_subroutine()
    RETURN

    CONTAINS
    SUBROUTINE check
        CALL ppm_check_topoid(topoid,valid,info)
        IF (.NOT. valid) THEN
            fail("topoid not valid",ppm_err_argument,info,8888)
        ENDIF
        IF (ppm_topo(topoid)%t%nsubs .LE. 0) THEN
            fail("nsubs mush be >0",ppm_err_argument,info,8888)
        ENDIF
        IF (SIZE(Offset,1) .NE. ppm_dim) THEN
            fail("invalid size for Offset. Should be ppm_dim",ppm_err_argument,info,8888)
        ENDIF

        IF (PRESENT(Nm)) THEN
            IF (PRESENT(h)) THEN
                fail("cannot specify both Nm and h. Choose only one.",ppm_err_argument,info,8888)
            ENDIF
            !TODO: check that the domain is finite
            IF (SIZE(Nm,1) .NE. ppm_dim) THEN
                fail("invalid size for Nm. Should be ppm_dim",ppm_err_argument,info,8888)
            ENDIF
            DO i=1,ppm_dim
                IF (Nm(i) .LT. 2) THEN
                    fail("Nm must be >1 in all space dimensions",ppm_err_argument,info,8888)
                ENDIF
            ENDDO

        ENDIF
        IF (PRESENT(h)) THEN
            IF (SIZE(h,1) .NE. ppm_dim) THEN
                fail("invalid size for h. Should be ppm_dim",ppm_err_argument,info,8888)
            ENDIF
            IF (ANY (h .LE. ppm_myepsd)) THEN
                fail("h must be >0 in all space dimensions",ppm_err_argument,info,8888)
            ENDIF
        ENDIF
        8888     CONTINUE
    END SUBROUTINE check
END SUBROUTINE equi_mesh_create

SUBROUTINE equi_mesh_destroy(this,info)
    !!! Destructor for the cartesian mesh object
    !-------------------------------------------------------------------------
    !  Modules
    !-------------------------------------------------------------------------
    USE ppm_module_topo_typedef
    USE ppm_module_check_id
    IMPLICIT NONE
    !-------------------------------------------------------------------------
    !  Includes
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    CLASS(ppm_t_equi_mesh)                  :: this
    !!! cartesian mesh object
    INTEGER                 , INTENT(  OUT) :: info
    !!! Returns status, 0 upon success
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER , DIMENSION(3)    :: ldc
    INTEGER                   :: iopt,ld,ud,kk,i,j,isub
    LOGICAL                   :: valid
    TYPE(ppm_t_topo), POINTER :: topo => NULL()

    start_subroutine("equi_mesh_destroy")

    !-------------------------------------------------------------------------
    !  (Re)allocate memory for the internal mesh list and Arrays at meshid
    !-------------------------------------------------------------------------
    ldc = 1
    iopt   = ppm_param_dealloc
    CALL ppm_alloc(this%Offset,ldc,iopt,info)
        or_fail_dealloc("Offset")
    CALL ppm_alloc(this%Nm,ldc,iopt,info)
        or_fail_dealloc("Nm")
    CALL ppm_alloc(this%h,ldc,iopt,info)
        or_fail_dealloc("h")

    destroy_collection_ptr(this%subpatch)
    destroy_collection_ptr(this%patch)
    destroy_collection_ptr(this%sub)


    !Destroy the bookkeeping entries in the fields that are
    !discretized on this mesh
    !TODO !!!


    destroy_collection_ptr(this%field_ptr)

    this%ID = 0

    end_subroutine()
END SUBROUTINE equi_mesh_destroy

FUNCTION equi_mesh_new_subpatch_data_ptr(this,info) RESULT(sp)
    !!! returns a pointer to a new subpatch_data object
    IMPLICIT NONE
    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    CLASS(ppm_t_equi_mesh)                  :: this
    !!! cartesian mesh object
    CLASS(ppm_t_subpatch_data_),POINTER     :: sp
    INTEGER                 , INTENT(  OUT) :: info
    !!! Returns status, 0 upon success
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER , DIMENSION(3)    :: ldc
    start_subroutine("equi_mesh_new_subpatch_data_ptr")

    ALLOCATE(ppm_t_subpatch_data::sp,STAT=info)
    or_fail_alloc("could not allocate ppm_t_subpatch_data pointer")

    end_subroutine()
END FUNCTION equi_mesh_new_subpatch_data_ptr

!CREATE
SUBROUTINE field_info_create(this,fieldID,info)
    !!! Constructor for subdomain data data structure
    CLASS(ppm_t_field_info)               :: this
    INTEGER,                  INTENT(IN)  :: fieldID
    INTEGER,                  INTENT(OUT) :: info

    start_subroutine("field_info_create")

    this%fieldID = fieldID

    end_subroutine()
END SUBROUTINE field_info_create
!DESTROY
SUBROUTINE field_info_destroy(this,info)
    !!! Destructor for subdomain data data structure
    CLASS(ppm_t_field_info)             :: this
    INTEGER,               INTENT(OUT) :: info

    start_subroutine("field_info_destroy")

    this%fieldID = 0

    end_subroutine()
END SUBROUTINE field_info_destroy

FUNCTION equi_mesh_list_of_fields(this,info) RESULT(fids)
    !!! Returns a pointer to an array containing the IDs of all the
    !!! fields that are currently discretized on this mesh
    CLASS(ppm_t_equi_mesh)          :: this
    INTEGER,DIMENSION(:),POINTER    :: fids
    INTEGER,            INTENT(OUT) :: info

    INTEGER                         :: i,j

    start_function("equi_mesh_list_of_fields")

    IF (.NOT.ASSOCIATED(this%field_ptr)) THEN
        fids => NULL()
        RETURN
    ENDIF
    IF (this%field_ptr%nb.LE.0) THEN
        fids => NULL()
        RETURN
    ENDIF

    ALLOCATE(fids(this%field_ptr%nb),STAT=info)
        or_fail_alloc("fids")

    j=1
    DO i=this%field_ptr%min_id,this%field_ptr%max_id
        IF (ASSOCIATED(this%field_ptr%vec(i)%t)) THEN
            fids(j) = this%field_ptr%vec(i)%t%fieldID
            j=j+1
        ENDIF
    ENDDO

    RETURN

    end_function()
END FUNCTION

SUBROUTINE equi_mesh_set_rel(this,field,info)
    !!! Create bookkeeping data structure to log the relationship between
    !!! the mesh and a field
    CLASS(ppm_t_equi_mesh)             :: this
    CLASS(ppm_t_field_)                :: field
    !!! this mesh is discretized on that field
    INTEGER,               INTENT(OUT) :: info

    TYPE(ppm_t_field_info),POINTER     :: fdinfo => NULL()

    start_subroutine("equi_mesh_set_rel")


    ALLOCATE(fdinfo,STAT=info)
        or_fail_alloc("could not allocate new field_info pointer")


    CALL fdinfo%create(field%ID,info)
        or_fail("could not create new field_info object")

    IF (.NOT.ASSOCIATED(this%field_ptr)) THEN
        ALLOCATE(ppm_c_field_info::this%field_ptr,STAT=info)
            or_fail_alloc("could not allocate field_info collection")
    ENDIF

    IF (this%field_ptr%vec_size.LT.field%ID) THEN
        CALL this%field_ptr%grow_size(info)
            or_fail("could not grow_size this%field_ptr")
    ENDIF
    !if the collection is still too small, it means we should consider
    ! using a hash table for meshIDs...
    IF (this%field_ptr%vec_size.LT.field%ID) THEN
        fail("The id of the mesh that we are trying to store in a Field collection seems to large. Why not implement a hash function?")
    ENDIF

    IF (this%field_ptr%exists(field%ID)) THEN
        fail("It seems like this Mesh already has a pointer to that field. Are you sure you want to do that a second time?")
    ELSE
        !add field_info to the collection, at index fieldID
        !Update the counters nb,max_id and min_id accordingly
        !(this is not very clean without hash table)
        this%field_ptr%vec(field%ID)%t => fdinfo
        this%field_ptr%nb = this%field_ptr%nb + 1
        IF (this%field_ptr%nb.EQ.1) THEN
            this%field_ptr%max_id = field%ID
            this%field_ptr%min_id = field%ID
        ELSE
            this%field_ptr%max_id = MAX(this%field_ptr%max_id,field%ID)
            this%field_ptr%min_id = MIN(this%field_ptr%min_id,field%ID)
        ENDIF
    ENDIF
    

    end_subroutine()
END SUBROUTINE


END MODULE ppm_module_mesh_typedef

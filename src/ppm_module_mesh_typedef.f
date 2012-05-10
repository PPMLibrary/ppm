!minclude ppm_header(ppm_module_mesh_typedef)

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
! Module variables 
!----------------------------------------------------------------------
INTEGER, PRIVATE, DIMENSION(4)  :: ldc

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
minclude ppm_create_collection(subpatch_data,subpatch_data,generate="extend")
minclude ppm_create_collection(subpatch_data,subpatch_data,generate="extend",vec=true)

TYPE,EXTENDS(ppm_t_mesh_discr_data_) :: ppm_t_mesh_discr_data
    CONTAINS
    PROCEDURE :: create    => mesh_discr_data_create
    PROCEDURE :: destroy   => mesh_discr_data_destroy
END TYPE
minclude ppm_create_collection(mesh_discr_data,mesh_discr_data,generate="extend")
minclude ppm_create_collection(mesh_discr_data,mesh_discr_data,vec=true,generate="extend")

TYPE,EXTENDS(ppm_t_subpatch_) :: ppm_t_subpatch
    CONTAINS
    PROCEDURE  :: create    => subpatch_create
    PROCEDURE  :: destroy   => subpatch_destroy
    PROCEDURE  :: get_pos   => subpatch_get_pos

    PROCEDURE  :: subpatch_get_field_2d_rd
    PROCEDURE  :: subpatch_get_field_3d_rd
    PROCEDURE  :: subpatch_get_field_4d_rd
END TYPE
minclude ppm_create_collection(subpatch,subpatch,generate="extend")

TYPE,EXTENDS(ppm_t_A_subpatch_) :: ppm_t_A_subpatch
    CONTAINS
    PROCEDURE :: create  => subpatch_A_create
    PROCEDURE :: destroy => subpatch_A_destroy
END TYPE
minclude ppm_create_collection(A_subpatch,A_subpatch,generate="extend")

TYPE,EXTENDS(ppm_t_equi_mesh_) :: ppm_t_equi_mesh
    CONTAINS
    PROCEDURE  :: create                => equi_mesh_create
    PROCEDURE  :: destroy               => equi_mesh_destroy
    PROCEDURE  :: create_prop           => equi_mesh_create_prop
    PROCEDURE  :: def_patch             => equi_mesh_def_patch
    PROCEDURE  :: set_rel               => equi_mesh_set_rel
    PROCEDURE  :: def_uniform           => equi_mesh_def_uniform
    PROCEDURE  :: new_subpatch_data_ptr => equi_mesh_new_subpatch_data_ptr
    PROCEDURE  :: list_of_fields        => equi_mesh_list_of_fields
    PROCEDURE  :: block_intersect       => equi_mesh_block_intersect
    PROCEDURE  :: map_ghost_init        => equi_mesh_map_ghost_init
    PROCEDURE  :: map_ghost_get         => equi_mesh_map_ghost_get
    PROCEDURE  :: map_ghost_push        => equi_mesh_map_ghost_push
    PROCEDURE  :: map_ghost_pop         => equi_mesh_map_ghost_pop
    PROCEDURE  :: map_send              => equi_mesh_map_send
    PROCEDURE  :: print_vtk             => equi_mesh_print_vtk
END TYPE
minclude ppm_create_collection(equi_mesh,equi_mesh,generate="extend")

!----------------------------------------------------------------------
! DATA STORAGE for the meshes
!----------------------------------------------------------------------
TYPE(ppm_c_equi_mesh)              :: ppm_mesh

!------------------------------------------------
! TODO: stuff that should be moved somewhere else:
!------------------------------------------------
!used to be in ppm_module_map_field_ghost.f
         INTEGER, DIMENSION(:  ), POINTER :: isendfromsub  => NULL()
         INTEGER, DIMENSION(:  ), POINTER :: isendtosub    => NULL()
         INTEGER, DIMENSION(:,:), POINTER :: isendpatchid    => NULL()
         INTEGER, DIMENSION(:  ), POINTER :: sendbuf       => NULL()
         INTEGER, DIMENSION(:  ), POINTER :: recvbuf       => NULL()
         INTEGER, DIMENSION(:,:), POINTER :: isendblkstart => NULL()
         INTEGER, DIMENSION(:,:), POINTER :: isendblksize  => NULL()
         INTEGER, DIMENSION(:,:), POINTER :: ioffset       => NULL()
         ! sorted (according to proc-proc interaction order) offset list)
         INTEGER, DIMENSION(:,:), POINTER :: mesh_ghost_offset => NULL()

         PRIVATE :: isendfromsub,isendtosub,sendbuf,recvbuf,isendblkstart
         PRIVATE :: isendblksize,ioffset,mesh_ghost_offset,isendpatchid

         INTEGER, DIMENSION(:), POINTER :: invsublist => NULL()
         INTEGER, DIMENSION(:), POINTER :: sublist    => NULL()

         PRIVATE :: invsublist,sublist

!used to be in ppm_module_map_field.f
         REAL(ppm_kind_single), DIMENSION(:), POINTER :: sends => NULL()
         REAL(ppm_kind_single), DIMENSION(:), POINTER :: recvs => NULL()
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: sendd => NULL()
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: recvd => NULL()
         INTEGER, DIMENSION(:), POINTER   :: nsend => NULL()
         INTEGER, DIMENSION(:), POINTER   :: nrecv => NULL()
         INTEGER, DIMENSION(:), POINTER   :: psend => NULL()
         INTEGER, DIMENSION(:), POINTER   :: precv => NULL()
         INTEGER, DIMENSION(:,:), POINTER :: pp => NULL()
         INTEGER, DIMENSION(:,:), POINTER :: qq => NULL()

         PRIVATE :: sends,recvs,sendd,recvd,nsend,nrecv,psend,precv,qq,pp
         
!----------------------------------------------------------------------
!  Type-bound procedures
!----------------------------------------------------------------------
CONTAINS

!Procedures for collections of derived types
minclude ppm_create_collection_procedures(subpatch_data,subpatch_data_)
minclude ppm_create_collection_procedures(subpatch_data,subpatch_data_,vec=true)
minclude ppm_create_collection_procedures(mesh_discr_data,mesh_discr_data_)
minclude ppm_create_collection_procedures(mesh_discr_data,mesh_discr_data_,vec=true)
minclude ppm_create_collection_procedures(subpatch,subpatch_)
minclude ppm_create_collection_procedures(A_subpatch,A_subpatch_)
minclude ppm_create_collection_procedures(equi_mesh,equi_mesh_)

SUBROUTINE subpatch_get_field_4d_rd(this,wp,Field,info)
    !!! Returns a pointer to the data array for a given field on this subpatch
    CLASS(ppm_t_subpatch)                :: this
    CLASS(ppm_t_field_)                  :: Field
    REAL(ppm_kind_double),DIMENSION(:,:,:,:),POINTER :: wp
    INTEGER,                 INTENT(OUT) :: info
    INTEGER                              :: p_idx

    start_subroutine("subpatch_get_field_4d")

    check_true("Field%lda+ppm_dim.GT.4",&
        "wrong dimensions for pointer arg wp")

    !Direct access to the data arrays 

    check_associated(this%mesh)
    
    p_idx = Field%get_pid(this%mesh)

    check_associated("this%subpatch_data")
    stdout("p_idx = ",p_idx)
    check_associated("this%subpatch_data%vec(p_idx)%t")

    wp => this%subpatch_data%vec(p_idx)%t%data_4d_rd

    check_associated(wp)

    end_subroutine()
END SUBROUTINE

SUBROUTINE subpatch_get_field_3d_rd(this,wp,Field,info)
    !!! Returns a pointer to the data array for a given field on this subpatch
    CLASS(ppm_t_subpatch)                :: this
    CLASS(ppm_t_field_)                  :: Field
    REAL(ppm_kind_double),DIMENSION(:,:,:),POINTER :: wp
    INTEGER,                 INTENT(OUT) :: info
    INTEGER                              :: p_idx

    start_subroutine("subpatch_get_field_3d")

    check_true("Field%lda+ppm_dim.GT.3",&
        "wrong dimensions for pointer arg wp")

    !Direct access to the data arrays 

    check_associated(this%mesh)
    
    p_idx = Field%get_pid(this%mesh)

    check_associated("this%subpatch_data")
    stdout("p_idx = ",p_idx)
    check_associated("this%subpatch_data%vec(p_idx)%t")

    wp => this%subpatch_data%vec(p_idx)%t%data_3d_rd

    check_associated(wp)

    end_subroutine()
END SUBROUTINE

SUBROUTINE subpatch_get_field_2d_rd(this,wp,Field,info)
    !!! Returns a pointer to the data array for a given field on this subpatch
    CLASS(ppm_t_subpatch)                :: this
    CLASS(ppm_t_field_)                  :: Field
    REAL(ppm_kind_double),DIMENSION(:,:),POINTER :: wp
    INTEGER,                 INTENT(OUT) :: info
    INTEGER                              :: p_idx

    INTEGER                            :: iopt, ndim

    start_subroutine("subpget_field_2d")

    check_true("Field%lda+ppm_dim.GT.2",&
        "wrong dimensions for pointer arg wp")

    !Direct access to the data arrays 
    p_idx = Field%get_pid(this%mesh)

    check_associated("this%subpatch_data")
    stdout("p_idx = ",p_idx)
    check_associated("this%subpatch_data%vec(p_idx)%t")

    wp => this%subpatch_data%vec(p_idx)%t%data_2d_rd

    check_associated(wp)

    end_subroutine()
END SUBROUTINE

SUBROUTINE subpatch_data_create(this,discr_data,sp,info)
    !!! Constructor for subdomain_data data structure
    CLASS(ppm_t_subpatch_data)              :: this
    CLASS(ppm_t_mesh_discr_data_),TARGET,  INTENT(IN) :: discr_data
    !!! field that is discretized on this mesh patch
    CLASS(ppm_t_subpatch_),      INTENT(IN) :: sp
    !!! subpatch to which this subpatch_data belongs
    INTEGER,                    INTENT(OUT) :: info

    INTEGER                                 :: ndim, datatype
    INTEGER,DIMENSION(ppm_dim+1)            :: hi,lo

    start_subroutine("subpatch_data_create")

    datatype =  discr_data%data_type
    this%discr_data => discr_data

    ! Determine allocation size of the data array
    IF (MINVAL(sp%nnodes(1:ppm_dim)) .LE. 0) THEN
        or_fail("invalid size for patch data. This patch should be deleted")
    ENDIF


    ndim = ppm_dim
    IF (discr_data%lda.LT.2) THEN
        !lo(1:ndim) = sp%istart_g(1:ppm_dim)
        !hi(1:ndim) = sp%iend_g(1:ppm_dim)
        lo(1) = 1 - sp%ghostsize(1)
        hi(1) = sp%nnodes(1) + sp%ghostsize(2)
        lo(2) = 1 - sp%ghostsize(3)
        hi(2) = sp%nnodes(2) + sp%ghostsize(4)
        IF (ppm_dim.EQ.3) THEN
            lo(3) = 1 - sp%ghostsize(5)
            hi(3) = sp%nnodes(3) + sp%ghostsize(6)
        ENDIF
    ELSE
        ndim = ndim +1
        lo(1) = 1
        hi(1) = discr_data%lda
        !lo(2:ndim) = sp%istart_g(1:ppm_dim)
        !hi(2:ndim) = sp%iend_g(1:ppm_dim)
        lo(2) = 1 - sp%ghostsize(1)
        hi(2) = sp%nnodes(1) + sp%ghostsize(2)
        lo(3) = 1 - sp%ghostsize(3)
        hi(3) = sp%nnodes(2) + sp%ghostsize(4)
        IF (ppm_dim.EQ.3) THEN
            lo(4) = 1 - sp%ghostsize(5)
            hi(4) = sp%nnodes(3) + sp%ghostsize(6)
        ENDIF
    ENDIF

    SELECT CASE (ndim)
    CASE (2)
        SELECT CASE (datatype)
        CASE (ppm_type_int)
            alloc_pointer_with_bounds2("this%data_2d_i",lo,hi)
        CASE (ppm_type_real_single)
            alloc_pointer_with_bounds2("this%data_2d_rs",lo,hi)
        CASE (ppm_type_real)
            alloc_pointer_with_bounds2("this%data_2d_rd",lo,hi)
        CASE (ppm_type_comp_single)
            alloc_pointer_with_bounds2("this%data_2d_cs",lo,hi)
        CASE (ppm_type_comp)
            alloc_pointer_with_bounds2("this%data_2d_cd",lo,hi)
        CASE (ppm_type_logical)
            alloc_pointer_with_bounds2("this%data_2d_l",lo,hi)
        CASE DEFAULT
            stdout("datatype = ",datatype)
            stdout("ndim = ",ndim)
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,caller,   &
                &        'invalid type for mesh patch data',__LINE__,info)
        END SELECT
    CASE (3)
        SELECT CASE (datatype)
        CASE (ppm_type_int)
            alloc_pointer_with_bounds3("this%data_3d_i",lo,hi)
        CASE (ppm_type_real_single)
            alloc_pointer_with_bounds3("this%data_3d_rs",lo,hi)
        CASE (ppm_type_real)
            alloc_pointer_with_bounds3("this%data_3d_rd",lo,hi)
        CASE (ppm_type_comp_single)
            alloc_pointer_with_bounds3("this%data_3d_cs",lo,hi)
        CASE (ppm_type_comp)
            alloc_pointer_with_bounds3("this%data_3d_cd",lo,hi)
        CASE (ppm_type_logical)
            alloc_pointer_with_bounds3("this%data_3d_l",lo,hi)
        CASE DEFAULT
            stdout("datatype = ")
            stdout(datatype)
            stdout("ndim = ")
            stdout(ndim)
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,caller,   &
                &        'invalid type for mesh patch data',__LINE__,info)
        END SELECT
    CASE (4)
        SELECT CASE (datatype)
        CASE (ppm_type_int)
            alloc_pointer_with_bounds4("this%data_4d_i",lo,hi)
        CASE (ppm_type_real_single)
            alloc_pointer_with_bounds4("this%data_4d_rs",lo,hi)
        CASE (ppm_type_real)
            !CALL ppm_alloc(this%data_4d_rd,ldc,iopt,info)
            alloc_pointer_with_bounds4("this%data_4d_rd",lo,hi)
        CASE (ppm_type_comp_single)
            alloc_pointer_with_bounds4("this%data_4d_cs",lo,hi)
        CASE (ppm_type_comp)
            alloc_pointer_with_bounds4("this%data_4d_cd",lo,hi)
        CASE (ppm_type_logical)
            alloc_pointer_with_bounds4("this%data_4d_l",lo,hi)
        CASE DEFAULT
            stdout("datatype = ")
            stdout(datatype)
            stdout("ndim = ")
            stdout(ndim)
            fail("invalid type for mesh patch data")
        END SELECT

    END SELECT

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
    this%discr_data => NULL()

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
SUBROUTINE subpatch_create(p,mesh,istart,iend,istart_p,iend_p,ghostsize,info)
    !!! Constructor for subpatch data structure
    CLASS(ppm_t_subpatch)              :: p
    CLASS(ppm_t_equi_mesh_),TARGET     :: mesh
    INTEGER,DIMENSION(:)               :: istart
    !!! Lower left corner of the subpatch
    INTEGER,DIMENSION(:)               :: iend
    !!! Upper right corner of the subpatch
    INTEGER,DIMENSION(:)               :: istart_p
    !!! Lower left corner of the patch
    INTEGER,DIMENSION(:)               :: iend_p
    !!! Upper right corner of the patch
    INTEGER,DIMENSION(:)               :: ghostsize
    !!! ghostlayer size in each direction
    INTEGER,               INTENT(OUT) :: info
    !!! return status. On success 0

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
    IF (.NOT.ASSOCIATED(p%istart_p)) THEN
        ALLOCATE(p%istart_p(ppm_dim),STAT=info)
            or_fail_alloc("could not allocate p%istart_p")
    ENDIF
    IF (.NOT.ASSOCIATED(p%iend_p)) THEN
        ALLOCATE(p%iend_p(ppm_dim),STAT=info)
            or_fail_alloc("could not allocate p%iend_p")
    ENDIF
    IF (.NOT.ASSOCIATED(p%istart_g)) THEN
        ALLOCATE(p%istart_g(ppm_dim),STAT=info)
            or_fail_alloc("could not allocate p%istart_g")
    ENDIF
    IF (.NOT.ASSOCIATED(p%iend_g)) THEN
        ALLOCATE(p%iend_g(ppm_dim),STAT=info)
            or_fail_alloc("could not allocate p%iend_g")
    ENDIF
    IF (.NOT.ASSOCIATED(p%nnodes)) THEN
        ALLOCATE(p%nnodes(ppm_dim),STAT=info)
            or_fail_alloc("could not allocate p%nnodes")
    ENDIF
    IF (.NOT.ASSOCIATED(p%ghostsize)) THEN
        ALLOCATE(p%ghostsize(2*ppm_dim),STAT=info)
            or_fail_alloc("could not allocate p%ghostsize")
    ENDIF

    p%meshID = mesh%ID
    p%mesh   => mesh
    p%istart = istart
    p%iend   = iend
    p%istart_p = istart_p
    p%iend_p   = iend_p
    p%nnodes(1:ppm_dim) = 1 + iend(1:ppm_dim) - istart(1:ppm_dim)
    p%ghostsize(1:2*ppm_dim) = ghostsize(1:2*ppm_dim)
    !p%istart_g(1:ppm_dim) = p%istart(1:ppm_dim)-ghostsize(1:ppm_dim)
    !p%iend_g(1:ppm_dim)   = p%iend(1:ppm_dim)+ghostsize(ppm_dim+1:2*ppm_dim)

    p%istart_g(1) = p%istart(1) - ghostsize(1)
    p%iend_g(1)   = p%iend(1)   + ghostsize(2)
    p%istart_g(2) = p%istart(2) - ghostsize(3)
    p%iend_g(2)   = p%iend(2)   + ghostsize(4)

    IF (ppm_dim.EQ.3) THEN
        p%istart_g(3) = p%istart(3) - ghostsize(5)
        p%iend_g(3)   = p%iend(3)   + ghostsize(6)
    ENDIF

    IF (.NOT.ASSOCIATED(p%subpatch_data)) THEN
        ALLOCATE(ppm_c_subpatch_data::p%subpatch_data,STAT=info)
            or_fail_alloc("could not allocate p%subpatch_data")
    ENDIF

    end_subroutine()
END SUBROUTINE subpatch_create

!DESTROY
SUBROUTINE subpatch_destroy(p,info)
    !!! Destructor for subpatch data structure
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
    CALL ppm_alloc(p%istart_g,ldc,iopt,info)
        or_fail_dealloc("p%istart_g")
    CALL ppm_alloc(p%iend_g,ldc,iopt,info)
        or_fail_dealloc("p%iend_g")
    CALL ppm_alloc(p%istart_p,ldc,iopt,info)
        or_fail_dealloc("p%istart_p")
    CALL ppm_alloc(p%iend_p,ldc,iopt,info)
        or_fail_dealloc("p%iend_p")
    CALL ppm_alloc(p%ghostsize,ldc,iopt,info)
        or_fail_dealloc("p%ghostsize")

    p%meshID = 0
    p%mesh => NULL()

    destroy_collection_ptr(p%subpatch_data)

    end_subroutine()
END SUBROUTINE subpatch_destroy

!GET POSITION
PURE FUNCTION subpatch_get_pos(p,i,j,k) RESULT (pos)
    !!! Return position of mesh node i,j,k
    CLASS(ppm_t_subpatch),  INTENT(IN) :: p
    INTEGER,                INTENT(IN) :: i
    INTEGER,                INTENT(IN) :: j
    INTEGER,OPTIONAL,       INTENT(IN) :: k
    REAL(ppm_kind_double),DIMENSION(ppm_dim) :: pos

    SELECT TYPE(mesh => p%mesh)
    TYPE IS (ppm_t_equi_mesh)
        pos(1) = (p%istart(1)+i)*mesh%h(1)
        pos(2) = (p%istart(2)+j)*mesh%h(2)
        IF (PRESENT(k)) pos(3) = (p%istart(3)+k)*mesh%h(3)
        pos(1:ppm_dim) = pos(1:ppm_dim) + mesh%offset(1:ppm_dim)
    END SELECT

END FUNCTION subpatch_get_pos

!CREATE
SUBROUTINE subpatch_A_create(this,vecsize,info,patchid)
    !!! Constructor for subpatch Arrays data structure
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
    !!! Destructor for subpatch Arrays data structure
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


SUBROUTINE equi_mesh_def_patch(this,patch,info,patchid,infinite)
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
    LOGICAL, OPTIONAL                       :: infinite
    !!! true if the patch should cover the whole domain
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    TYPE(ppm_t_topo), POINTER :: topo => NULL()
    INTEGER                   :: i,j,isub,id,pid
    INTEGER                   :: size2,size_tmp,nsubpatch,nsubpatchi
    CLASS(ppm_t_subpatch_),  POINTER :: p => NULL()
    CLASS(ppm_t_A_subpatch_),POINTER :: A_p => NULL()
    INTEGER,              DIMENSION(ppm_dim) :: istart,iend
    INTEGER,              DIMENSION(ppm_dim) :: istart_p,iend_p
    INTEGER,              DIMENSION(2*ppm_dim) :: ghostsize
    REAL(ppm_kind_double),DIMENSION(ppm_dim) :: h,Offset
    TYPE(ppm_t_ptr_subpatch),DIMENSION(:),POINTER :: tmp_array => NULL()
    LOGICAL                                  :: linfinite

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
    topo => ppm_topo(this%topoid)%t

    check_false("associated(topo%min_subs)",&
        "Mesh_def_patch is not yet implemented for single precision")

    !-------------------------------------------------------------------------
    ! Extent of the patch on the global mesh (before it is divided into
    ! subpatches)
    !-------------------------------------------------------------------------
    IF (PRESENT(infinite)) THEN
        linfinite = infinite
    ELSE
        linfinite = .FALSE.
    ENDIF

    IF (linfinite) THEN
        istart_p(1:ppm_dim) = -HUGE(1)/2
        iend_p(1:ppm_dim)   =  HUGE(1)/2
    ELSE
        istart_p(1:ppm_dim) = 1 + &
            CEILING((   patch(1:ppm_dim)     - Offset(1:ppm_dim))/h(1:ppm_dim))
        iend_p(1:ppm_dim)   = 1 + &
            FLOOR((patch(ppm_dim+1:2*ppm_dim)- Offset(1:ppm_dim))/h(1:ppm_dim))
    ENDIF

    !-------------------------------------------------------------------------
    !  Allocate bookkeeping arrays (pointers between patches and subpatches)
    !-------------------------------------------------------------------------

    ALLOCATE(ppm_t_A_subpatch::A_p,STAT=info)
        or_fail_alloc("could not allocate ppm_t_A_subpatch pointer")

    !SELECT TYPE(t => A_p)
    !TYPE IS (ppm_t_A_subpatch)
        !CALL t%create(topo%nsublist,info,patchid)
    CALL A_p%create(topo%nsublist,info,patchid)
        or_fail("could not initialize ppm_t_A_subpatch pointer")
    !END SELECT

    IF (PRESENT(patchid)) THEN
        pid = patchid
    ELSE
        pid = 0
    ENDIF
    !-------------------------------------------------------------------------
    !  intersect each subdomain with the patch and store the corresponding
    !  subpatch in the mesh data structure
    !-------------------------------------------------------------------------
    ! loop through all subdomains on this processor to allocate some book-
    ! keeping arrays.
    ! TODO: clean this up....
    size_tmp=0
    DO i = 1,topo%nsublist
        isub = topo%isublist(i)
        ! check if the subdomain overlaps with the patch
        IF (ALL(patch(1:ppm_dim).LT.topo%max_subd(1:ppm_dim,isub)) .AND. &
            ALL(patch(ppm_dim+1:2*ppm_dim).GT.topo%min_subd(1:ppm_dim,isub)))&
               THEN 
           ASSOCIATE (sarray => this%subpatch_by_sub(isub))
           IF (sarray%nsubpatch.GE.sarray%size) THEN
               size2 = MAX(2*sarray%size,5)
               IF (size_tmp.LT.size2) THEN
                   dealloc_pointer(tmp_array)
                   ALLOCATE(tmp_array(size2),STAT=info)
                       or_fail_alloc("tmp_array")
               ENDIF
               DO j=1,sarray%nsubpatch
                   tmp_array(j)%t => sarray%vec(j)%t
               ENDDO
               dealloc_pointer(sarray%vec)

               ALLOCATE(sarray%vec(size2),STAT=info)
                   or_fail_alloc("sarray%vec")
               DO j=1,sarray%nsubpatch
                   sarray%vec(j)%t => tmp_array(j)%t
               ENDDO
               sarray%size = size2
           ENDIF
           END ASSOCIATE
        ENDIF
    ENDDO
    dealloc_pointer(tmp_array)

    nsubpatch = 0
    ! loop through all subdomains on this processor
    DO i = 1,topo%nsublist
        isub = topo%isublist(i)
        nsubpatchi = this%subpatch_by_sub(isub)%nsubpatch
        ! check if the subdomain overlaps with the patch
        IF (ALL(patch(1:ppm_dim).LT.topo%max_subd(1:ppm_dim,isub)) .AND. &
            ALL(patch(ppm_dim+1:2*ppm_dim).GT.topo%min_subd(1:ppm_dim,isub)))&
               THEN 

            ! finds the mesh nodes which are contained in the overlap region
            istart(1:ppm_dim) = 1 + CEILING((MAX(patch(1:ppm_dim),&
                topo%min_subd(1:ppm_dim,isub))-Offset(1:ppm_dim))/h(1:ppm_dim))
            iend(1:ppm_dim)   = 1 + FLOOR((  MIN(patch(ppm_dim+1:2*ppm_dim),&
                topo%max_subd(1:ppm_dim,isub))-Offset(1:ppm_dim))/h(1:ppm_dim))

            ! create a new subpatch object
            ALLOCATE(ppm_t_subpatch::p,STAT=info)
                or_fail_alloc("could not allocate ppm_t_subpatch pointer")

            ! determine ghostlayer size for this subpatch
            ! (there is a ghostlayer if and only if the border of the subpatch
            ! does not coincide with a border of the patch itself - in that
            ! case, the width of the ghostlayer is truncated by the "mesh-wide"
            ! ghostsize parameter) 
            
            !TODO: deal with computational domain boundaries (e.g.
            ! when the patch is infinite)
            ghostsize(1) = MIN(istart(1)-istart_p(1),this%ghostsize(1))
            ghostsize(2) = MIN(iend_p(1)-iend(1)    ,this%ghostsize(1))
            ghostsize(3) = MIN(istart(2)-istart_p(2),this%ghostsize(2))
            ghostsize(4) = MIN(iend_p(2)-iend(2)    ,this%ghostsize(2))
            IF (ppm_dim.EQ.3) THEN
                ghostsize(5) = MIN(istart(3)-istart_p(3),this%ghostsize(3))
                ghostsize(6) = MIN(iend_p(3)-iend(3)    ,this%ghostsize(3))
            ENDIF

            CALL p%create(this,istart,iend,istart_p,iend_p,ghostsize,info)
                or_fail("could not create new subpatch")

            nsubpatch = nsubpatch+1
            A_p%subpatch(nsubpatch)%t => p

            !add a pointer to this subpatch
            nsubpatchi = nsubpatchi + 1
            this%subpatch_by_sub(isub)%vec(nsubpatchi)%t => p
            this%subpatch_by_sub(isub)%nsubpatch = nsubpatchi

            ! put the subpatch object in the collection of subpatches on this mesh
            CALL this%subpatch%push(p,info,id)
                or_fail("could not add new subpatch to mesh")
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
    CALL this%def_patch(patch,info,patchid,infinite=.TRUE.)
        or_fail("failed to add patch")

    !TODO add some checks for the finiteness of the computational domain

    end_subroutine()
END SUBROUTINE equi_mesh_def_uniform

SUBROUTINE equi_mesh_create(this,topoid,Offset,info,Nm,h,ghostsize)
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
    INTEGER,DIMENSION(:),              OPTIONAL,INTENT(IN   ) :: ghostsize
    !!! size of the ghost layer, in number of mesh nodes for each dimension
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER , DIMENSION(3)    :: ldc
    INTEGER                   :: iopt,ld,ud,kk,i,j,isub,nsubs
    LOGICAL                   :: valid
    INTEGER, DIMENSION(ppm_dim)   :: Nc
    REAL(ppm_kind_double), DIMENSION(ppm_dim)  :: len_phys,rat
    REAL(ppm_kind_double)                      :: lmyeps
    CHARACTER(LEN=ppm_char)                    :: msg
    TYPE(ppm_t_topo), POINTER                  :: topo => NULL()

    start_subroutine("equi_mesh_create")

    !-------------------------------------------------------------------------
    !  Check arguments
    !-------------------------------------------------------------------------
    IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
    ENDIF

    lmyeps = ppm_myepsd


    !This mesh is defined for a given topology
    this%topoid = topoid

    !dumb way of creating a global ID for this mesh
    !TODO find something better? (needed if one creates and destroy
    ! many meshes)
    ppm_nb_meshes = ppm_nb_meshes + 1
    this%ID = ppm_nb_meshes 

    topo => ppm_topo(topoid)%t

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

    CALL ppm_alloc(this%ghostsize,ldc,iopt,info)
        or_fail_alloc('ghostsize')

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

    nsubs = topo%nsubs
    IF (.NOT.ASSOCIATED(this%subpatch_by_sub)) THEN
        ALLOCATE(this%subpatch_by_sub(nsubs),STAT=info)
            or_fail_alloc("could not allocate this%subpatch_by_sub")
    ELSE
        fail("subpatch_by_sub is already allocated. Call destroy() first?")
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

    IF (PRESENT(ghostsize)) THEN
        this%ghostsize(1:ppm_dim) = ghostsize(1:ppm_dim)
    ELSE
        this%ghostsize(1:ppm_dim) = 0
    ENDIF

    !-------------------------------------------------------------------------
    !  Allocate memory for the subdomain indices
    !-------------------------------------------------------------------------
    iopt = ppm_param_alloc_fit
    ldc(1) = ppm_dim
    ldc(2) = nsubs
    CALL ppm_alloc(this%istart,ldc,iopt,info)
    or_fail_alloc("could not allocate this%istart")
    CALL ppm_alloc(this%iend,ldc,iopt,info)
    or_fail_alloc("could not allocate this%iend")

    !-------------------------------------------------------------------------
    !  Determine the start indices in the global mesh for each
    ! subdomain
    !-------------------------------------------------------------------------
    DO i=1,nsubs
            this%istart(1:ppm_dim,i) = 1 + CEILING((     &
                topo%min_subd(1:ppm_dim,i)-Offset(1:ppm_dim))/this%h(1:ppm_dim))
            this%iend(1:ppm_dim,i)   = 1 + FLOOR((       &
                topo%max_subd(1:ppm_dim,i)-Offset(1:ppm_dim))/this%h(1:ppm_dim))
    ENDDO

    !DO i=1,nsubs
        !this%istart(1:ppm_dim,i) = &
            !NINT((topo%min_subd(1:ppm_dim,i)-&
            !topo%min_physd(1:ppm_dim))/this%h(1:ppm_dim)) + 1
    !ENDDO
    !!-------------------------------------------------------------------------
    !!  Determine number of mesh points for each sub. 
    !!  2D and 3D case have separate loops for vectorization.
    !!-------------------------------------------------------------------------
    !IF (ppm_dim .EQ. 3) THEN
        !DO i=1,nsubs
            !len_phys(1) = topo%max_subd(1,i)-topo%min_subd(1,i)
            !len_phys(2) = topo%max_subd(2,i)-topo%min_subd(2,i)
            !len_phys(3) = topo%max_subd(3,i)-topo%min_subd(3,i)
            !rat(1) = len_phys(1)/this%h(1)
            !rat(2) = len_phys(2)/this%h(2)
            !rat(3) = len_phys(3)/this%h(3)
            !Nc(1)  = NINT(rat(1))
            !Nc(2)  = NINT(rat(2))
            !Nc(3)  = NINT(rat(3))
            !IF (ABS(rat(1)-REAL(Nc(1),ppm_kind_double)) .GT. lmyeps*rat(1)) THEN
                !WRITE(msg,'(2(A,F12.6))') 'in dimension 1: sub_length=',  &
                    !&                len_phys(1),' mesh_spacing=',this%h(1)
                !fail(msg,ppm_err_subs_incomp)
            !ENDIF
            !IF (ABS(rat(2)-REAL(Nc(2),ppm_kind_double)) .GT. lmyeps*rat(2)) THEN
                !WRITE(msg,'(2(A,F12.6))') 'in dimension 2: sub_length=',  &
                    !&                len_phys(2),' mesh_spacing=',this%h(2)
                !fail(msg,ppm_err_subs_incomp)
            !ENDIF
            !IF (ABS(rat(3)-REAL(Nc(3),ppm_kind_double)) .GT. lmyeps*rat(3)) THEN
                !WRITE(msg,'(2(A,F12.6))') 'in dimension 3: sub_length=',  &
                    !&                len_phys(3),' mesh_spacing=',this%h(3)
                !fail(msg,ppm_err_subs_incomp)
            !ENDIF
            !this%iend(1,i) = this%istart(1,i) + Nc(1) + 1 
            !this%iend(2,i) = this%istart(2,i) + Nc(2) + 1 
            !this%iend(3,i) = this%istart(3,i) + Nc(3) + 1 
        !ENDDO
    !ELSE
        !DO i=1,nsubs
            !len_phys(1) = topo%max_subd(1,i)-topo%min_subd(1,i)
            !len_phys(2) = topo%max_subd(2,i)-topo%min_subd(2,i)
            !rat(1) = len_phys(1)/this%h(1)
            !rat(2) = len_phys(2)/this%h(2)
            !Nc(1)  = NINT(rat(1))
            !Nc(2)  = NINT(rat(2))
            !IF (ABS(rat(1)-REAL(Nc(1),ppm_kind_double)) .GT. lmyeps*rat(1)) THEN
                !WRITE(msg,'(2(A,F12.6))') 'in dimension 1: sub_length=',  &
                    !&                len_phys(1),' mesh_spacing=',this%h(1)
                !fail(msg,ppm_err_subs_incomp)
            !ENDIF
            !IF (ABS(rat(2)-REAL(Nc(2),ppm_kind_double)) .GT. lmyeps*rat(2)) THEN
                !WRITE(msg,'(2(A,F12.6))') 'in dimension 2: sub_length=',  &
                    !&                len_phys(2),' mesh_spacing=',this%h(2)
                !fail(msg,ppm_err_subs_incomp)
            !ENDIF
            !this%iend(1,i) = this%istart(1,i) + Nc(1) + 1 
            !this%iend(2,i) = this%istart(2,i) + Nc(2) + 1 
        !ENDDO
    !ENDIF

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
            fail("invalid size for Offset. Should be ppm_dim",&
                ppm_err_argument,info,8888)
        ENDIF

        IF (PRESENT(Nm)) THEN
            IF (PRESENT(h)) THEN
                fail("cannot specify both Nm and h. Choose only one.",&
                    ppm_err_argument,info,8888)
            ENDIF
            !TODO: check that the domain is finite
            IF (SIZE(Nm,1) .NE. ppm_dim) THEN
                fail("invalid size for Nm. Should be ppm_dim",&
                    ppm_err_argument,info,8888)
            ENDIF
            DO i=1,ppm_dim
                IF (Nm(i) .LT. 2) THEN
                    fail("Nm must be >1 in all space dimensions",&
                        ppm_err_argument,info,8888)
                ENDIF
            ENDDO

        ENDIF
        IF (PRESENT(h)) THEN
            IF (SIZE(h,1) .NE. ppm_dim) THEN
                fail("invalid size for h. Should be ppm_dim",&
                    ppm_err_argument,info,8888)
            ENDIF
            IF (ANY (h .LE. ppm_myepsd)) THEN
                fail("h must be >0 in all space dimensions",&
                    ppm_err_argument,info,8888)
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

    CALL ppm_alloc(this%istart,ldc,iopt,info)
        or_fail_dealloc("istart")
    CALL ppm_alloc(this%iend,ldc,iopt,info)
        or_fail_dealloc("iend")

    destroy_collection_ptr(this%subpatch)
    destroy_collection_ptr(this%patch)
    IF (ASSOCIATED(this%subpatch_by_sub)) THEN
        DEALLOCATE(this%subpatch_by_sub,STAT=info)
            or_fail_dealloc("subpatch_by_sub")
    ENDIF


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


FUNCTION equi_mesh_list_of_fields(this,info) RESULT(fids)
    !!! Returns a pointer to an array containing the IDs of all the
    !!! fields that are currently discretized on this mesh
    CLASS(ppm_t_equi_mesh)          :: this
    INTEGER,DIMENSION(:),POINTER    :: fids

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


    CALL fdinfo%create(field,info)
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


SUBROUTINE equi_mesh_map_ghost_push(this,field,info)
    !!! Push field data onto the mesh mappings buffers
    CLASS(ppm_t_equi_mesh)             :: this
    CLASS(ppm_t_field_)                :: field
    !!! this mesh is discretized on that field
    INTEGER,               INTENT(OUT) :: info

    INTEGER                            :: p_idx
    REAL(ppm_kind_double),DIMENSION(:,:),POINTER    :: wp2_dummy => NULL()
    REAL(ppm_kind_double),DIMENSION(:,:,:),POINTER  :: wp3_dummy => NULL()
    REAL(ppm_kind_double),DIMENSION(:,:,:,:),POINTER:: wp4_dummy => NULL()

    start_subroutine("equi_mesh_map_ghost_push")


    !p_idx = field%M%vec(this%ID)%t%p_idx
    p_idx = field%get_pid(this)

    IF (ppm_dim.EQ.2) THEN
        IF (field%lda.EQ.1) THEN
            CALL ppm_map_field_push_2d_sca_d(this,wp2_dummy,p_idx,info)
            or_fail("map_field_push_2d")
        ELSE
            CALL ppm_map_field_push_2d_vec_d(this,wp3_dummy,field%lda,p_idx,info)
            or_fail("map_field_push_2d")
        ENDIF
    ELSE
        IF (field%lda.EQ.1) THEN
            CALL ppm_map_field_push_3d_sca_d(this,wp3_dummy,p_idx,info)
            or_fail("map_field_push_3d")
        ELSE
            CALL ppm_map_field_push_3d_vec_d(this,wp4_dummy,field%lda,p_idx,info)
            or_fail("map_field_push_3d")
        ENDIF
    ENDIF

    end_subroutine()
END SUBROUTINE equi_mesh_map_ghost_push


SUBROUTINE equi_mesh_map_ghost_pop(this,field,info)
    !!! Push field data onto the mesh mappings buffers
    CLASS(ppm_t_equi_mesh)             :: this
    CLASS(ppm_t_field_)                :: field
    !!! this mesh is discretized on that field
    INTEGER,               INTENT(OUT) :: info

    INTEGER                            :: p_idx
    REAL(ppm_kind_double),DIMENSION(:,:),POINTER    :: wp2_dummy => NULL()
    REAL(ppm_kind_double),DIMENSION(:,:,:),POINTER  :: wp3_dummy => NULL()
    REAL(ppm_kind_double),DIMENSION(:,:,:,:),POINTER:: wp4_dummy => NULL()

    start_subroutine("equi_mesh_map_ghost_pop")


    !p_idx = field%M%vec(this%ID)%t%p_idx
    p_idx = field%get_pid(this)

    IF (ppm_dim.EQ.2) THEN
        IF (field%lda.EQ.1) THEN
            CALL ppm_map_field_pop_2d_sca_d(this,wp2_dummy,p_idx,info)
            or_fail("map_field_pop_2d")
        ELSE
            CALL ppm_map_field_pop_2d_vec_d(this,wp3_dummy,field%lda,p_idx,info)
            or_fail("map_field_pop_2d")
        ENDIF
    ELSE
        IF (field%lda.EQ.1) THEN
            CALL ppm_map_field_pop_3d_sca_d(this,wp3_dummy,p_idx,info)
            or_fail("map_field_pop_3d")
        ELSE
            CALL ppm_map_field_pop_3d_vec_d(this,wp4_dummy,field%lda,p_idx,info)
            or_fail("map_field_pop_3d")
        ENDIF
    ENDIF

    end_subroutine()
END SUBROUTINE equi_mesh_map_ghost_pop

SUBROUTINE mesh_discr_data_create(this,field,info)
    CLASS(ppm_t_mesh_discr_data)            :: this
    CLASS(ppm_t_field_),TARGET, INTENT(IN)  :: field
    INTEGER,                    INTENT(OUT) :: info
    start_subroutine("mesh_discr_data_create")


    this%data_type = field%data_type
    this%name = field%name
    this%lda = field%lda
    this%field_ptr => field
    ALLOCATE(ppm_v_subpatch_data::this%subpatch,STAT=info)
        or_fail_alloc("this%subpatch")

    end_subroutine()
END SUBROUTINE
SUBROUTINE mesh_discr_data_destroy(this,info)
    CLASS(ppm_t_mesh_discr_data)        :: this
    INTEGER,                INTENT(OUT) :: info
    start_subroutine("mesh_discr_data_destroy")

    this%data_type = 0
    this%field_ptr => NULL()
    this%lda = 0
    this%name = ""

    end_subroutine()
END SUBROUTINE

SUBROUTINE equi_mesh_create_prop(this,field,discr_data,info,p_idx)
    CLASS(ppm_t_equi_mesh)                            :: this
    CLASS(ppm_t_field_),                  INTENT(IN)  :: field
    CLASS(ppm_t_mesh_discr_data_),POINTER,INTENT(OUT) :: discr_data
    INTEGER,                              INTENT(OUT) :: info
    INTEGER,OPTIONAL,                     INTENT(OUT) :: p_idx

    CLASS(ppm_t_subpatch_),     POINTER :: p => NULL()
    CLASS(ppm_t_subpatch_data_),POINTER :: subpdat => NULL()

    start_subroutine("equi_mesh_create_prop")


    IF (.NOT.ASSOCIATED(this%mdata)) THEN
        ALLOCATE(ppm_v_mesh_discr_data::this%mdata,STAT=info)
        or_fail_alloc("mdata")
    ENDIF


    ALLOCATE(ppm_t_mesh_discr_data::discr_data,STAT=info)
            or_fail_alloc("discr_data")
    CALL discr_data%create(field,info)
            or_fail("discr_data%create()")

    !Create a new data array on the mesh to store this field
    p => this%subpatch%begin()
    DO WHILE (ASSOCIATED(p))
        ! create a new subpatch_data object
        subpdat => this%new_subpatch_data_ptr(info)
            or_fail_alloc("could not get a new ppm_t_subpatch_data pointer")

        CALL subpdat%create(discr_data,p,info)
            or_fail("could not create new subpatch_data")

        CALL discr_data%subpatch%push(subpdat,info)
            or_fail("could not add new subpatch_data to discr_data")

        CALL p%subpatch_data%push(subpdat,info,p_idx)
            or_fail("could not add new subpatch_data to subpatch collection")

        p => this%subpatch%next()
    ENDDO

    CALL this%mdata%push(discr_data,info)
        or_fail("could not add new discr_data to discr%mdata")


    !check that the returned pointer makes sense
    check_associated("discr_data")


    end_subroutine()
END SUBROUTINE equi_mesh_create_prop


SUBROUTINE equi_mesh_print_vtk(this,filename,info)
    USE ppm_module_io_vtk
    CLASS(ppm_t_equi_mesh)                            :: this
    CHARACTER(LEN=*)                                  :: filename
    INTEGER,                              INTENT(OUT) :: info

    start_subroutine("equi_mesh_print_vtk")

    IF (ppm_dim.EQ.2) THEN
        CALL ppm_vtk_fields_2d(filename,this,info)
        or_fail("ppm_vtk_fields_2d")
    ELSE
        CALL ppm_vtk_fields_3d(filename,this,info)
        or_fail("ppm_vtk_fields_3d")
    ENDIF

    end_subroutine()
END SUBROUTINE equi_mesh_print_vtk

#include "mesh/mesh_block_intersect.f"
#include "mesh/mesh_map_ghost_init.f"
#include "mesh/mesh_map_ghost_get.f"

#include "mesh/mesh_map_send.f"

#define __SFIELD 1
#define __VFIELD 2
#define __DOUBLE_PRECISION 17
#define __DIM __VFIELD
#define __KIND __DOUBLE_PRECISION
#include "mesh/mesh_map_push_2d.f"
#include "mesh/mesh_map_pop_2d.f"
#include "mesh/mesh_map_push_3d.f"
#include "mesh/mesh_map_pop_3d.f"
#undef __DIM
#define __DIM __SFIELD
#include "mesh/mesh_map_push_2d.f"
#include "mesh/mesh_map_pop_2d.f"
#include "mesh/mesh_map_push_3d.f"
#include "mesh/mesh_map_pop_3d.f"
#undef __DIM
#undef __SFIELD
#undef __VFIELD
#undef __DOUBLE_PRECISION

END MODULE ppm_module_mesh_typedef

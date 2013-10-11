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
          PROCEDURE :: create    => subpatch_create
          PROCEDURE :: destroy   => subpatch_destroy
          PROCEDURE :: get_pos2d => subpatch_get_pos2d
          PROCEDURE :: get_pos3d => subpatch_get_pos3d

          PROCEDURE :: subpatch_get_field_2d_rs
          PROCEDURE :: subpatch_get_field_3d_rs
          PROCEDURE :: subpatch_get_field_4d_rs
          PROCEDURE :: subpatch_get_field_2d_rd
          PROCEDURE :: subpatch_get_field_3d_rd
          PROCEDURE :: subpatch_get_field_4d_rd
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
          PROCEDURE :: create                => equi_mesh_create
          PROCEDURE :: destroy               => equi_mesh_destroy
          PROCEDURE :: create_prop           => equi_mesh_create_prop
          PROCEDURE :: zero                  => equi_mesh_prop_zero
          PROCEDURE :: def_patch             => equi_mesh_def_patch
          PROCEDURE :: def_uniform           => equi_mesh_def_uniform
          PROCEDURE :: new_subpatch_data_ptr => equi_mesh_new_subpatch_data_ptr
          PROCEDURE :: list_of_fields        => equi_mesh_list_of_fields
          PROCEDURE :: block_intersect       => equi_mesh_block_intersect
          PROCEDURE :: map_ghost_init        => equi_mesh_map_ghost_init
          PROCEDURE :: map_ghost_get         => equi_mesh_map_ghost_get
          PROCEDURE :: map_ghost_put         => equi_mesh_map_ghost_put
          PROCEDURE :: map_ghost_push        => equi_mesh_map_push
          PROCEDURE :: map_ghost_pop         => equi_mesh_map_pop
          PROCEDURE :: map_send              => equi_mesh_map_send

          PROCEDURE :: map                   => equi_mesh_map_global
          PROCEDURE :: map_push              => equi_mesh_map_push
          PROCEDURE :: map_pop               => equi_mesh_map_pop

          PROCEDURE :: print_vtk             => equi_mesh_print_vtk
          PROCEDURE :: interp_to_part        => equi_mesh_m2p
      END TYPE
minclude ppm_create_collection(equi_mesh,equi_mesh,generate="extend")

      !----------------------------------------------------------------------
      ! DATA STORAGE for the meshes
      !----------------------------------------------------------------------
      TYPE(ppm_c_equi_mesh) :: ppm_mesh

      !------------------------------------------------
      ! TODO: stuff that should be moved somewhere else:
      !------------------------------------------------
      !used to be in ppm_module_map_field.f
      REAL(ppm_kind_single), DIMENSION(:),   PRIVATE, POINTER :: sends => NULL()
      REAL(ppm_kind_single), DIMENSION(:),   PRIVATE, POINTER :: recvs => NULL()
      REAL(ppm_kind_double), DIMENSION(:),   PRIVATE, POINTER :: sendd => NULL()
      REAL(ppm_kind_double), DIMENSION(:),   PRIVATE, POINTER :: recvd => NULL()

      INTEGER,               DIMENSION(:),   PRIVATE, POINTER :: nsend => NULL()
      INTEGER,               DIMENSION(:),   PRIVATE, POINTER :: nrecv => NULL()
      INTEGER,               DIMENSION(:),   PRIVATE, POINTER :: psend => NULL()
      INTEGER,               DIMENSION(:),   PRIVATE, POINTER :: precv => NULL()

      INTEGER,               DIMENSION(:,:), PRIVATE, POINTER :: pp => NULL()
      INTEGER,               DIMENSION(:,:), PRIVATE, POINTER :: qq => NULL()

      INTEGER,               DIMENSION(:  ), PRIVATE, POINTER :: irecvfromsub  => NULL()
      INTEGER,               DIMENSION(:  ), PRIVATE, POINTER :: irecvtosub    => NULL()
      INTEGER,               DIMENSION(:,:), PRIVATE, POINTER :: irecvpatchid  => NULL()
      INTEGER,               DIMENSION(:,:), PRIVATE, POINTER :: irecvblkstart => NULL()
      INTEGER,               DIMENSION(:,:), PRIVATE, POINTER :: irecvblksize  => NULL()

      !used to be in ppm_module_map_field_ghost.f
      INTEGER,               DIMENSION(:  ), PRIVATE, POINTER :: isendfromsub  => NULL()
      INTEGER,               DIMENSION(:  ), PRIVATE, POINTER :: isendtosub    => NULL()
      INTEGER,               DIMENSION(:,:), PRIVATE, POINTER :: isendpatchid  => NULL()
      INTEGER,               DIMENSION(:,:), PRIVATE, POINTER :: isendblkstart => NULL()
      INTEGER,               DIMENSION(:,:), PRIVATE, POINTER :: isendblksize  => NULL()
      INTEGER,               DIMENSION(:,:), PRIVATE, POINTER :: ioffset       => NULL()

      INTEGER,               DIMENSION(:  ), PRIVATE, POINTER :: sendbuf       => NULL()
      INTEGER,               DIMENSION(:  ), PRIVATE, POINTER :: recvbuf       => NULL()

      ! sorted (according to proc-proc interaction order) offset list)
      INTEGER,               DIMENSION(:,:), PRIVATE, POINTER :: mesh_ghost_offset => NULL()

      INTEGER,               DIMENSION(:),   PRIVATE, POINTER :: invsublist => NULL()
      INTEGER,               DIMENSION(:),   PRIVATE, POINTER :: sublist    => NULL()


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

      ! home-made templating system for the get_field routines
minclude ppm_get_field_template(2,s)
minclude ppm_get_field_template(3,s)
minclude ppm_get_field_template(4,s)
minclude ppm_get_field_template(2,d)
minclude ppm_get_field_template(3,d)
minclude ppm_get_field_template(4,d)

      SUBROUTINE subpatch_data_create(this,discr_data,sp,info)
          !!! Constructor for subdomain_data data structure

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_subpatch_data)                           :: this
          CLASS(ppm_t_mesh_discr_data_),TARGET,  INTENT(IN   ) :: discr_data
          !!! field that is discretized on this mesh patch
          CLASS(ppm_t_subpatch_),                INTENT(IN   ) :: sp
          !!! subpatch to which this subpatch_data belongs
          INTEGER,                               INTENT(  OUT) :: info

          INTEGER                                              :: ndim, datatype
          INTEGER, DIMENSION(ppm_dim+1)                        :: hi,lo

          start_subroutine("subpatch_data_create")

          datatype =  discr_data%data_type
          this%discr_data => discr_data

          ! Determine allocation size of the data array
          IF (MINVAL(sp%nnodes(1:ppm_dim)) .LE. 0) THEN
              or_fail("invalid size for patch data. This patch should be deleted")
          ENDIF


          ndim = ppm_dim
          IF (discr_data%lda.LT.2) THEN
              lo(1:ppm_dim) = sp%lo_a(1:ppm_dim)
              hi(1:ppm_dim) = sp%hi_a(1:ppm_dim)
          ELSE
              ndim = ndim +1
              lo(1) = 1
              hi(1) = discr_data%lda
              lo(2:ndim) = sp%lo_a(1:ppm_dim)
              hi(2:ndim) = sp%hi_a(1:ppm_dim)
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
                  &    'invalid type for mesh patch data',__LINE__,info)
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
                  &    'invalid type for mesh patch data',__LINE__,info)
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

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
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
      SUBROUTINE subpatch_create(p,mesh,isub,istart,iend,pstart,pend,&
              istart_p,iend_p,ghostsize,bcdef,info)
          !!! Constructor for subpatch data structure

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_subpatch)              :: p
          CLASS(ppm_t_equi_mesh_),TARGET     :: mesh
          !! mesh
          INTEGER                            :: isub
          !! subdomain
          INTEGER,DIMENSION(:)               :: istart
          !!! Lower left corner of the subpatch (integer coord on the mesh)
          INTEGER,DIMENSION(:)               :: iend
          !!! Upper right corner of the subpatch (integer coord on the mesh)
          REAL(ppm_kind_double),DIMENSION(:) :: pstart
          !!! Lower left corner of the subpatch (absolute coord)
          REAL(ppm_kind_double),DIMENSION(:) :: pend
          !!! Upper right corner of the subpatch (absolute coord)
          INTEGER,DIMENSION(:)               :: istart_p
          !!! Lower left corner of the patch
          INTEGER,DIMENSION(:)               :: iend_p
          !!! Upper right corner of the patch
          INTEGER,DIMENSION(:)               :: ghostsize
          !!! ghostlayer size in each direction
          INTEGER,DIMENSION(:)               :: bcdef
          !!! boundary conditions
          INTEGER,               INTENT(OUT) :: info
          !!! return status. On success 0

          INTEGER                            :: iopt,i

          start_subroutine("subpatch_create")

          !Allocate arrays
          alloc_pointer("p%istart",ppm_dim)
          alloc_pointer("p%iend",ppm_dim)
          alloc_pointer("p%start",ppm_dim)
          alloc_pointer("p%end",ppm_dim)
          alloc_pointer("p%start_ext",ppm_dim)
          alloc_pointer("p%end_ext",ppm_dim)
          alloc_pointer("p%start_red",ppm_dim)
          alloc_pointer("p%end_red",ppm_dim)
          alloc_pointer("p%istart_p",ppm_dim)
          alloc_pointer("p%iend_p",ppm_dim)
          alloc_pointer("p%lo_a",ppm_dim)
          alloc_pointer("p%hi_a",ppm_dim)
          alloc_pointer("p%nnodes",ppm_dim)
          alloc_pointer("p%ghostsize",'2*ppm_dim')
          alloc_pointer("p%bc",'2*ppm_dim')

          p%meshID   = mesh%ID
          p%mesh     => mesh
          p%isub     = isub
          p%istart   = istart
          p%iend     = iend
          p%start    = pstart
          p%end      = pend
          p%istart_p = istart_p
          p%iend_p   = iend_p
          p%nnodes(1:ppm_dim) = 1 + iend(1:ppm_dim) - istart(1:ppm_dim)
          p%ghostsize(1:2*ppm_dim) = ghostsize(1:2*ppm_dim)
          p%bc(1:2*ppm_dim) = bcdef(1:2*ppm_dim)

          !Define allocated size of the subpatch
          p%lo_a(1)     = 1 - ghostsize(1)
          p%hi_a(1)     = p%nnodes(1)   + ghostsize(2)
          p%lo_a(2)     = 1 - ghostsize(3)
          p%hi_a(2)     = p%nnodes(2)   + ghostsize(4)
          IF (ppm_dim.EQ.3) THEN
             p%lo_a(3) = 1 - ghostsize(5)
             p%hi_a(3) = p%nnodes(3)   + ghostsize(6)
          ENDIF

          DO i=1,ppm_dim
              !Extended subpatch (abs. coord)
              p%start_ext(i) = pstart(i)- mesh%ghostsize(i) * mesh%h(i)
              p%end_ext(i)   = pend(i)  + mesh%ghostsize(i) * mesh%h(i)

              !Reduced subpatch (abs. coord)
              p%start_red(i) = pstart(i)-(p%ghostsize(2*i-1)-mesh%ghostsize(i))*mesh%h(i)
              p%end_red(i)   = pend(i)  -(mesh%ghostsize(i)-p%ghostsize(2*i)) * mesh%h(i)
          ENDDO


          IF (.NOT.ASSOCIATED(p%subpatch_data)) THEN
             ALLOCATE(ppm_c_subpatch_data::p%subpatch_data,STAT=info)
             or_fail_alloc("could not allocate p%subpatch_data")
          ENDIF

          end_subroutine()

      END SUBROUTINE subpatch_create

      !DESTROY
      SUBROUTINE subpatch_destroy(p,info)
          !!! Destructor for subpatch data structure

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
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
          CALL ppm_alloc(p%start,ldc,iopt,info)
          or_fail_dealloc("p%start")
          CALL ppm_alloc(p%end,ldc,iopt,info)
          or_fail_dealloc("p%end")
          CALL ppm_alloc(p%start_ext,ldc,iopt,info)
          or_fail_dealloc("p%start_ext")
          CALL ppm_alloc(p%end_ext,ldc,iopt,info)
          or_fail_dealloc("p%end_ext")
          CALL ppm_alloc(p%start_red,ldc,iopt,info)
          or_fail_dealloc("p%start_red")
          CALL ppm_alloc(p%end_red,ldc,iopt,info)
          or_fail_dealloc("p%end_red")
          CALL ppm_alloc(p%lo_a,ldc,iopt,info)
          or_fail_dealloc("p%lo_a")
          CALL ppm_alloc(p%hi_a,ldc,iopt,info)
          or_fail_dealloc("p%hi_a")
          CALL ppm_alloc(p%istart_p,ldc,iopt,info)
          or_fail_dealloc("p%istart_p")
          CALL ppm_alloc(p%iend_p,ldc,iopt,info)
          or_fail_dealloc("p%iend_p")
          CALL ppm_alloc(p%ghostsize,ldc,iopt,info)
          or_fail_dealloc("p%ghostsize")

          p%meshID = 0
          p%isub = 0
          p%mesh => NULL()

          destroy_collection_ptr(p%subpatch_data)

          end_subroutine()

      END SUBROUTINE subpatch_destroy

      !GET POSITION
      PURE FUNCTION subpatch_get_pos2d(p,i,j) RESULT (pos)
          !!! Return position of mesh node i,j,k

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_subpatch),  INTENT(IN) :: p
          INTEGER,                INTENT(IN) :: i
          INTEGER,                INTENT(IN) :: j
          REAL(ppm_kind_double),DIMENSION(ppm_dim) :: pos

          SELECT TYPE(mesh => p%mesh)
          TYPE IS (ppm_t_equi_mesh)
              !numbering of mesh nodes goes from 1 to N
              pos(1) = (p%istart(1)+i-2)*mesh%h(1) + mesh%offset(1)
              pos(2) = (p%istart(2)+j-2)*mesh%h(2) + mesh%offset(2)
          END SELECT

      END FUNCTION subpatch_get_pos2d
      PURE FUNCTION subpatch_get_pos3d(p,i,j,k) RESULT (pos)
          !!! Return position of mesh node i,j,k

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_subpatch),  INTENT(IN) :: p
          INTEGER,                INTENT(IN) :: i
          INTEGER,                INTENT(IN) :: j
          INTEGER,                INTENT(IN) :: k
          REAL(ppm_kind_double),DIMENSION(ppm_dim) :: pos

          SELECT TYPE(mesh => p%mesh)
          TYPE IS (ppm_t_equi_mesh)
              !numbering of mesh nodes goes from 1 to N
              pos(1) = (p%istart(1)+i-2)*mesh%h(1) + mesh%offset(1)
              pos(2) = (p%istart(2)+j-2)*mesh%h(2) + mesh%offset(2)
              pos(3) = (p%istart(3)+k-2)*mesh%h(3) + mesh%offset(3)
          END SELECT

      END FUNCTION subpatch_get_pos3d

      !CREATE
      SUBROUTINE subpatch_A_create(this,vecsize,info,patchid)
          !!! Constructor for subpatch Arrays data structure

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_A_subpatch)            :: this
          INTEGER                            :: vecsize
          INTEGER,               INTENT(OUT) :: info
          INTEGER,OPTIONAL,      INTENT(IN)  :: patchid

          start_subroutine("subpatch_A_create")

          IF (ASSOCIATED(this%subpatch)) THEN
             CALL this%destroy(info)
             or_fail_dealloc("could not destroy this ppm_t_A_subpatch object")
          ENDIF

          ALLOCATE(ppm_t_ptr_subpatch::this%subpatch(vecsize),STAT=info)
          or_fail_alloc("could not allocate this%subpatch")

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

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
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

      SUBROUTINE equi_mesh_create(this,topoid,Offset,info,Nm,h,ghostsize,name)
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
          CLASS(ppm_t_equi_mesh),                        INTENT(INOUT) :: this
          !!! cartesian mesh object
          INTEGER,                                       INTENT(IN   ) :: topoid
          !!! Topology ID for which mesh has been created
          REAL(ppm_kind_double), DIMENSION(:),           INTENT(IN   ) :: Offset
          !!! Offset in each dimension
          INTEGER,                                       INTENT(  OUT) :: info
          !!! Returns status, 0 upon success
          INTEGER,               DIMENSION(:), OPTIONAL, INTENT(IN   ) :: Nm
          !!! Global number of mesh points in the whole comput. domain
          !!! Makes sense only if the computational domain is bounded.
          !!! Note: Exactly one of Nm and h should be specified
          REAL(ppm_kind_double), DIMENSION(:), OPTIONAL, INTENT(IN   ) :: h
          !!! Mesh spacing
          !!! Note: Exactly one of Nm and h should be specified
          INTEGER,               DIMENSION(:), OPTIONAL, INTENT(IN   ) :: ghostsize
          !!! size of the ghost layer, in number of mesh nodes for each dimension
          CHARACTER(LEN=*),                    OPTIONAL, INTENT(IN   ) :: name
          !!! name of this mesh
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          TYPE(ppm_t_topo), POINTER :: topo => NULL()

          REAL(ppm_kind_double), DIMENSION(ppm_dim) :: len_phys,rat

          INTEGER                     :: iopt,ld,ud,kk
          INTEGER                     :: i,j,k,isub,nsubs
          INTEGER, DIMENSION(ppm_dim) :: Nc
          INTEGER, DIMENSION(3)       :: ldc

          CHARACTER(LEN=ppm_char) :: msg

          LOGICAL :: valid

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

          !By default, there are no patches defined on this mesh yet.
          this%npatch = 0

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

          IF (.NOT.ASSOCIATED(this%field_ptr)) THEN
              ALLOCATE(this%field_ptr,STAT=info)
              or_fail_alloc("could not allocate this%field_ptr")
          ELSE
              fail("field_ptr is already allocated. Call destroy() first?")
          ENDIF

          !-------------------------------------------------------------------------
          !  Store the mesh information
          !-------------------------------------------------------------------------
          this%Offset(1:ppm_dim) = Offset(1:ppm_dim)

          IF (PRESENT(h)) THEN
              this%h(1:ppm_dim) = h(1:ppm_dim)
              IF (ASSOCIATED(topo%max_physs)) THEN
                 this%Nm(1:ppm_dim) = FLOOR((topo%max_physs(1:ppm_dim) - &
                 &                    topo%min_physs(1:ppm_dim))/(h(1:ppm_dim))) + 1
              ELSE
                 this%Nm(1:ppm_dim) = FLOOR((topo%max_physd(1:ppm_dim) - &
                 &                    topo%min_physd(1:ppm_dim))/(h(1:ppm_dim))) + 1
              ENDIF
          ELSE
              this%Nm(1:ppm_dim) = Nm(1:ppm_dim)
              IF (ASSOCIATED(topo%max_physs)) THEN
                 this%h(1:ppm_dim) = (topo%max_physs(1:ppm_dim) - &
                 &                   topo%min_physs(1:ppm_dim))/  &
                 &                   REAL(Nm(1:ppm_dim)-1,ppm_kind_single)
                 !check for round-off problems and fix them if necessary
                 DO k=1,ppm_dim
                    DO WHILE (topo%min_physs(k)+(this%Nm(k)-1)*this%h(k).LT.topo%max_physs(k))
                       this%h(k)=this%h(k)+EPSILON(this%h(k))
                    ENDDO
                 ENDDO
                 check_true(<#ALL(topo%min_physs(1:ppm_dim)+(Nm(1:ppm_dim)-1)*this%h(1:ppm_dim).GE.topo%max_physs(1:ppm_dim))#>,"round-off problem in mesh creation")
              ELSE
                 this%h(1:ppm_dim) = (topo%max_physd(1:ppm_dim) - &
                 &                   topo%min_physd(1:ppm_dim))/  &
                 &                   REAL(Nm(1:ppm_dim)-1,ppm_kind_double)
                 !check for round-off problems and fix them if necessary
                 DO k=1,ppm_dim
                    DO WHILE (topo%min_physd(k)+(this%Nm(k)-1)*this%h(k).LT.topo%max_physd(k))
                       this%h(k)=this%h(k)+EPSILON(this%h(k))
                    ENDDO
                 ENDDO
                 check_true(<#ALL(topo%min_physd(1:ppm_dim)+(Nm(1:ppm_dim)-1)*this%h(1:ppm_dim).GE.topo%max_physd(1:ppm_dim))#>,"round-off problem in mesh creation")
              ENDIF
          ENDIF

          IF (PRESENT(ghostsize)) THEN
              this%ghostsize(1:ppm_dim) = ghostsize(1:ppm_dim)
          ELSE
              this%ghostsize(1:ppm_dim) = 0
          ENDIF

          IF (PRESENT(name)) THEN
              this%name = TRIM(ADJUSTL(name))
          ELSE
              this%name = "default_mesh_name"
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
          IF (ASSOCIATED(topo%min_subd)) Then
             DO i=1,nsubs
                this%istart(1:ppm_dim,i) = 1 + CEILING((topo%min_subd(1:ppm_dim,i) &
                &                        - Offset(1:ppm_dim))/this%h(1:ppm_dim))
                this%iend(1:ppm_dim,i)   = 1 + FLOOR((topo%max_subd(1:ppm_dim,i)   &
                &                        - Offset(1:ppm_dim))/(this%h(1:ppm_dim)   &
                &                        - EPSILON(this%h(1:ppm_dim))))
                !WARNING: this is a hack to resolve a round-off issue when h is such
                !that a node falls epsilon away from the subdomain boundary

                !stdout("#isub = ",i," BEFORE CHOP")
                !stdout("sub%istart = ",'this%istart(1:ppm_dim,i)')
                !stdout("sub%iend = ",'this%iend(1:ppm_dim,i)')
                !stdout("sub%bc = ",'topo%subs_bc(:,i)')

                !Decide what to do if a mesh node falls on a subdomain boundary
                DO k=1,ppm_dim
                   ! For internal boundaries, mesh nodes are not duplicated. If a
                   ! node is on the boundary, only the West/South/Bottom one
                   ! is kept.
                   ! For periodic boundary conditions, subdomains that have a mesh
                   ! node right on the East, North, or Top domain boundary are
                   ! reduced by 1 mesh node so that real mesh nodes are not
                   ! duplicated.
                   IF ((this%iend(k,i)-1).GE.(topo%max_subd(k,i)-Offset(k))/this%h(k)) THEN
                      IF (topo%subs_bc(2*k,i).NE.1  .OR.  &
                      &   (topo%subs_bc(2*k,i).EQ.1 .AND. &
                      &   topo%bcdef(2*k).EQ.ppm_param_bcdef_periodic)) THEN
                         this%iend(k,i) = this%iend(k,i) -1
                      ENDIF
                   ENDIF
                ENDDO

                !Check that this subdomain is large enough and thus
                !compatible with this mesh and its resolution h.
                check_true(<#ALL((topo%max_subd(1:ppm_dim,i)-topo%min_subd(1:ppm_dim,i)).GT.(this%h(1:ppm_dim)*this%ghostsize(1:ppm_dim)+ppm_myepsd))#>,&
                "Grid spacing h (times ghostsize) has to be stricly smaller than any subdomain.")
                !stdout("#isub = ",i," AFTER CHOP")
                !stdout("sub%istart = ",'this%istart(1:ppm_dim,i)')
                !stdout("sub%iend = ",'this%iend(1:ppm_dim,i)')
             ENDDO
          ELSE
             DO i=1,nsubs
                this%istart(1:ppm_dim,i) = 1 + CEILING((topo%min_subs(1:ppm_dim,i) &
                &                        - Offset(1:ppm_dim))/this%h(1:ppm_dim))
                this%iend(1:ppm_dim,i)   = 1 + FLOOR((topo%max_subs(1:ppm_dim,i)   &
                &                        - Offset(1:ppm_dim))/(this%h(1:ppm_dim)   &
                &                        - EPSILON(this%h(1:ppm_dim))))
                !WARNING: this is a hack to resolve a round-off issue when h is such
                !that a node falls epsilon away from the subdomain boundary

                !stdout("#isub = ",i," BEFORE CHOP")
                !stdout("sub%istart = ",'this%istart(1:ppm_dim,i)')
                !stdout("sub%iend = ",'this%iend(1:ppm_dim,i)')
                !stdout("sub%bc = ",'topo%subs_bc(:,i)')

                !Decide what to do if a mesh node falls on a subdomain boundary
                DO k=1,ppm_dim
                   ! For internal boundaries, mesh nodes are not duplicated. If a
                   ! node is on the boundary, only the West/South/Bottom one
                   ! is kept.
                   ! For periodic boundary conditions, subdomains that have a mesh
                   ! node right on the East, North, or Top domain boundary are
                   ! reduced by 1 mesh node so that real mesh nodes are not
                   ! duplicated.
                   IF ((this%iend(k,i)-1).GE.(topo%max_subs(k,i)-Offset(k))/this%h(k)) THEN
                      IF (topo%subs_bc(2*k,i).NE.1  .OR.  &
                      &   (topo%subs_bc(2*k,i).EQ.1 .AND. &
                      &   topo%bcdef(2*k).EQ.ppm_param_bcdef_periodic)) THEN
                         this%iend(k,i) = this%iend(k,i) -1
                      ENDIF
                   ENDIF
                ENDDO

                !Check that this subdomain is large enough and thus
                !compatible with this mesh and its resolution h.
                check_true(<#ALL((topo%max_subs(1:ppm_dim,i)-topo%min_subs(1:ppm_dim,i)).GT.(this%h(1:ppm_dim)*this%ghostsize(1:ppm_dim)+ppm_myepss))#>,&
                "Grid spacing h (times ghostsize) has to be stricly smaller than any subdomain.")
                !stdout("#isub = ",i," AFTER CHOP")
                !stdout("sub%istart = ",'this%istart(1:ppm_dim,i)')
                !stdout("sub%iend = ",'this%iend(1:ppm_dim,i)')
             ENDDO
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
          CLASS(ppm_t_equi_mesh), INTENT(INOUT) :: this
          !!! cartesian mesh object
          INTEGER,                INTENT(  OUT) :: info
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
          CALL ppm_alloc(this%Nm,ldc,iopt,info)
          or_fail_dealloc("Nm")
          CALL ppm_alloc(this%Offset,ldc,iopt,info)
          or_fail_dealloc("Offset")
          CALL ppm_alloc(this%h,ldc,iopt,info)
          or_fail_dealloc("h")
          CALL ppm_alloc(this%ghostsize,ldc,iopt,info)
          or_fail_dealloc("ghostsize")
          CALL ppm_alloc(this%nnodes,ldc,iopt,info)
          or_fail_dealloc("nnodes")
          CALL ppm_alloc(this%istart,ldc,iopt,info)
          or_fail_dealloc("istart")
          CALL ppm_alloc(this%iend,ldc,iopt,info)
          or_fail_dealloc("iend")

          destroy_collection_ptr(this%subpatch)
          destroy_collection_ptr(this%patch)
          destroy_collection_ptr(this%mdata)

          dealloc_pointer("this%subpatch_by_sub")

          !Destroy the bookkeeping entries in the fields that are
          !discretized on this mesh
          !TODO !!!

          destroy_collection_ptr(this%field_ptr)

          dealloc_pointers(this%ghost_fromsub,this%ghost_tosub,this%ghost_patchid,this%ghost_blkstart,this%ghost_blksize,this%ghost_blk,this%ghost_recvtosub,this%ghost_recvpatchid,this%ghost_recvblkstart,this%ghost_recvblksize,this%ghost_recvblk)

          this%ID = 0
          this%topoid = 0
          this%npatch = 0
          this%ghost_initialized = .FALSE.

          end_subroutine()
      END SUBROUTINE equi_mesh_destroy

      SUBROUTINE equi_mesh_create_prop(this,field,discr_data,info,p_idx)

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_equi_mesh),                 INTENT(INOUT) :: this
          CLASS(ppm_t_field_),                    INTENT(IN   ) :: field
          CLASS(ppm_t_mesh_discr_data_), POINTER, INTENT(  OUT) :: discr_data
          INTEGER,                                INTENT(  OUT) :: info
          INTEGER,OPTIONAL,                       INTENT(  OUT) :: p_idx

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
          check_associated(discr_data)

          end_subroutine()
      END SUBROUTINE equi_mesh_create_prop

      SUBROUTINE equi_mesh_prop_zero(this,Field,info)
          !!! Zero a property allocated on this mesh

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_equi_mesh), INTENT(INOUT) :: this
          CLASS(ppm_t_field_),    INTENT(IN   ) :: Field
          INTEGER,                INTENT(  OUT) :: info

          INTEGER :: lda

          start_subroutine("equi_mesh_prop_zero")

          lda = Field%lda

          SELECT CASE (ppm_dim)
          CASE (2)
              SELECT CASE (lda)
              CASE (1)
              foreach n in equi_mesh(this) with sca_fields(Field) prec(ppm_kind_double)
                  for all
                      Field_n = REAL(0,ppm_kind_double)
              end foreach
              CASE (2)
              foreach n in equi_mesh(this) with vec_fields(Field) prec(ppm_kind_double)
                  for all
                      Field_n(1) = REAL(0,ppm_kind_double)
                      Field_n(2) = REAL(0,ppm_kind_double)
              end foreach
              CASE (3)
              foreach n in equi_mesh(this) with vec_fields(Field) prec(ppm_kind_double)
                  for all
                      Field_n(1) = REAL(0,ppm_kind_double)
                      Field_n(2) = REAL(0,ppm_kind_double)
                      Field_n(3) = REAL(0,ppm_kind_double)
              end foreach
              CASE (4)
              foreach n in equi_mesh(this) with vec_fields(Field) prec(ppm_kind_double)
                  for all
                      Field_n(1) = REAL(0,ppm_kind_double)
                      Field_n(2) = REAL(0,ppm_kind_double)
                      Field_n(3) = REAL(0,ppm_kind_double)
                      Field_n(4) = REAL(0,ppm_kind_double)
              end foreach
              CASE DEFAULT
              foreach n in equi_mesh(this) with vec_fields(Field) prec(ppm_kind_double)
                  for all
                      Field_n(1:lda) = REAL(0,ppm_kind_double)
              end foreach
              END SELECT

          CASE DEFAULT
              SELECT CASE (lda)
              CASE (1)
              foreach n in equi_mesh(this) with sca_fields(Field) indices(i,j,k) prec(ppm_kind_double)
                  for all
                      Field_n = REAL(0,ppm_kind_double)
              end foreach
              CASE (2)
              foreach n in equi_mesh(this) with vec_fields(Field) indices(i,j,k) prec(ppm_kind_double)
                  for all
                      Field_n(1) = REAL(0,ppm_kind_double)
                      Field_n(2) = REAL(0,ppm_kind_double)
              end foreach
              CASE (3)
              foreach n in equi_mesh(this) with vec_fields(Field) indices(i,j,k) prec(ppm_kind_double)
                  for all
                      Field_n(1) = REAL(0,ppm_kind_double)
                      Field_n(2) = REAL(0,ppm_kind_double)
                      Field_n(3) = REAL(0,ppm_kind_double)
              end foreach
              CASE (4)
              foreach n in equi_mesh(this) with vec_fields(Field) indices(i,j,k) prec(ppm_kind_double)
                  for all
                      Field_n(1) = REAL(0,ppm_kind_double)
                      Field_n(2) = REAL(0,ppm_kind_double)
                      Field_n(3) = REAL(0,ppm_kind_double)
                      Field_n(4) = REAL(0,ppm_kind_double)
              end foreach
              CASE DEFAULT
              foreach n in equi_mesh(this) with vec_fields(Field) indices(i,j,k) prec(ppm_kind_double)
                  for all
                      Field_n(1:lda) = REAL(0,ppm_kind_double)
              end foreach
              END SELECT
          END SELECT

          end_subroutine()
      END SUBROUTINE equi_mesh_prop_zero

      SUBROUTINE equi_mesh_def_patch(this,patch,info,patchid,infinite,bcdef)
          !!! Add a patch to a mesh
          !!! The patch corners are given in terms of their absolute coordinates,
          !!! but the patch is then shrunk to the nearest mesh nodes.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_topo_typedef

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------

          CLASS(ppm_t_equi_mesh),                     INTENT(INOUT) :: this
          REAL(ppm_kind_double), DIMENSION(:),        INTENT(INOUT) :: patch
          !!! Positions of the corners of the patch
          !!! (x1,y1,z1,x2,y2,z2), where 1 is the lower-left-bottom corner
          !!! and 2 is the upper-right-top corner.
          INTEGER,                                    INTENT(  OUT) :: info
          !!! Returns status, 0 upon success
          INTEGER, OPTIONAL,                          INTENT(IN   ) :: patchid
          !!! id of the patch, if we want one.
          LOGICAL, OPTIONAL,                          INTENT(IN   ) :: infinite
          !!! true if the patch should cover the whole domain
          INTEGER, OPTIONAL,    DIMENSION(2*ppm_dim), INTENT(IN   ) :: bcdef
          !!! Boundary conditions. Default is free-space boundary conditions if
          !!! strictly inside the computational domain, or the bcdef specified
          !!! for the comput. domain if not.
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          TYPE(ppm_t_ptr_subpatch), DIMENSION(:), POINTER :: tmp_array => NULL()
          TYPE(ppm_t_topo),                       POINTER :: topo => NULL()

          CLASS(ppm_t_subpatch_),                 POINTER :: p => NULL()
          CLASS(ppm_t_A_subpatch_),               POINTER :: A_p => NULL()

          REAL(ppm_kind_double), DIMENSION(ppm_dim) :: h,Offset
          REAL(ppm_kind_double), DIMENSION(ppm_dim) :: pstart,pend

          INTEGER                        :: i,j,k,isub,jsub,id,pid
          INTEGER                        :: size2,size_tmp,iopt
          INTEGER                        :: nsubpatch,nsubpatchi
          INTEGER, DIMENSION(ppm_dim)    :: istart,iend
          INTEGER, DIMENSION(ppm_dim)    :: istart_p,iend_p
          INTEGER, DIMENSION(ppm_dim)    :: istart_d,iend_d
          INTEGER, DIMENSION(2*ppm_dim)  :: ghostsize
          INTEGER, DIMENSION(2*ppm_dim)  :: bc
          INTEGER, DIMENSION(1)          :: ldu
          INTEGER, DIMENSION(:), POINTER :: indsub => NULL()

          LOGICAL :: linfinite

          start_subroutine("mesh_def_patch")

          !-------------------------------------------------------------------------
          !  Check arguments
          !-------------------------------------------------------------------------
          check_associated(<#this%subpatch#>,&
          & "Mesh not allocated. Call Mesh%create() first?")

          !-------------------------------------------------------------------------
          !  get mesh parameters
          !-------------------------------------------------------------------------
          h = this%h
          Offset = this%Offset
          topo => ppm_topo(this%topoid)%t

          check_false(<#ASSOCIATED(topo%min_subs)#>,&
          & "Mesh_def_patch is not yet implemented for single precision")

          !-------------------------------------------------------------------------
          ! Extent of the patch on the global mesh (before it is divided into
          ! subpatches)
          !-------------------------------------------------------------------------
          IF (PRESENT(infinite)) THEN
             linfinite = infinite
          ELSE
             linfinite = .FALSE.
          ENDIF

          !Bounds for the mesh nodes that are inside the patch
          IF (linfinite) THEN
             istart_p(1:ppm_dim) = -HUGE(1)/2
             iend_p(1:ppm_dim)   =  HUGE(1)/2
          ELSE
             istart_p(1:ppm_dim) = 1 + &
             & CEILING((   patch(1:ppm_dim)     - Offset(1:ppm_dim))/h(1:ppm_dim))
             iend_p(1:ppm_dim)   = 1 + &
             & FLOOR((patch(ppm_dim+1:2*ppm_dim)- Offset(1:ppm_dim))/h(1:ppm_dim))
          ENDIF

          !Re-define the patch boundaries so that its corners fall on mesh nodes
          patch(1:ppm_dim)           = (istart_p(1:ppm_dim) - 1) * &
          &                             h(1:ppm_dim) + Offset(1:ppm_dim)
          patch(ppm_dim+1:2*ppm_dim) = (iend_p(1:ppm_dim) - 1)   * &
          &                             h(1:ppm_dim) + Offset(1:ppm_dim)

          !Bounds for the mesh nodes that are inside the computational
          !domain
          istart_d(1:ppm_dim) = 1 + CEILING( (topo%min_physd(1:ppm_dim)     &
          &                       -          Offset(1:ppm_dim))/h(1:ppm_dim) )
          iend_d(1:ppm_dim)   = 1 + FLOOR( (topo%max_physd(1:ppm_dim)       &
          &                       -        Offset(1:ppm_dim))/(h(1:ppm_dim) &
          &                       -        EPSILON(h(1:ppm_dim))) )

          !stdout("istart_d = ",istart_d)
          !stdout("iend_d = ",iend_d)
          !stdout("topo%bcdef = ",'topo%bcdef')
          !stdout("before roundoff: ",'(topo%max_physd(1:ppm_dim)- Offset(1:ppm_dim))/h(1:ppm_dim)')

          !-------------------------------------------------------------------------
          !  Allocate bookkeeping arrays (pointers between patches and subpatches)
          !-------------------------------------------------------------------------
          ALLOCATE(ppm_t_A_subpatch::A_p,STAT=info)
          or_fail_alloc("could not allocate ppm_t_A_subpatch pointer")

          IF (PRESENT(patchid)) THEN
             CALL A_p%create(topo%nsublist,info,patchid)
             or_fail("could not initialize ppm_t_A_subpatch pointer")

             pid = patchid
          ELSE
             CALL A_p%create(topo%nsublist,info)
             or_fail("could not initialize ppm_t_A_subpatch pointer")

             pid = 0
          ENDIF

          !-------------------------------------------------------------------------
          !  Build the indsub list to diffrentiate local sub indeices
          !  from global ones
          !-------------------------------------------------------------------------
          iopt   = ppm_param_alloc_fit
          ldu(1) = topo%nsubs
          CALL ppm_alloc(indsub,ldu,iopt,info)
          or_fail_alloc("indsub")

          indsub=-1
          DO i=1,topo%nsublist
             isub=topo%isublist(i)
             indsub(isub)=1
          ENDDO

          !-------------------------------------------------------------------------
          !  intersect each subdomain with the patch and store the corresponding
          !  subpatch in the mesh data structure
          !-------------------------------------------------------------------------
          ! loop through all subdomains on this processor to allocate some book-
          ! keeping arrays.
          size_tmp=0
          !DO i = 1,topo%nsublist
          DO isub=1,topo%nsubs
             !isub = topo%isublist(i)
             !----------------------------------------------------------------
             ! check if the subdomain overlaps with the patch
             ! if so, then there might be a subpatch created for that subdomain.
             ! Grow the size of the array of pointers this%subpach_by_sub.
             ! (simpler to do it this way if we dont know how many patches
             ! are going to be added)
             !----------------------------------------------------------------
             IF (ALL(patch(1:ppm_dim).LT.topo%max_subd(1:ppm_dim,isub)) .AND. &
             &   ALL(patch(ppm_dim+1:2*ppm_dim).GE.topo%min_subd(1:ppm_dim,isub))) THEN
                ! yaser: subpach_by_sub should contain global index
                ! so I would use isub instead of i
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
!           ! loop through all subdomains on this processor
!           sub: DO i = 1,topo%nsublist
          ! loop through all subdomains on this topology
          sub: DO isub=1,topo%nsubs
              !isub = topo%isublist(i)
              !how many subpatches do we already have here
              ! yaser: subpach_by_sub should contain global index
              ! so I would use isub instead of i
              nsubpatchi = this%subpatch_by_sub(isub)%nsubpatch
              !----------------------------------------------------------------
              ! check if the subdomain overlaps with the patch
              !----------------------------------------------------------------
              IF (ALL(patch(1:ppm_dim).LT.topo%max_subd(1:ppm_dim,isub)) .AND. &
              &   ALL(patch(ppm_dim+1:2*ppm_dim).GE.topo%min_subd(1:ppm_dim,isub))) THEN
                  !----------------------------------------------------------------
                  !Finds the mesh nodes which are contained in the overlap region
                  !----------------------------------------------------------------

                  !----------------------------------------------------------------
                  !Absolute positions of the corners (stored for use
                  !during m2p interpolation, where we need to check whether
                  !a particle is within a given subpatch)
                  !----------------------------------------------------------------
                  pstart(1:ppm_dim) = MAX(patch(1:ppm_dim),&
                  &   topo%min_subd(1:ppm_dim,isub))-Offset(1:ppm_dim)
                  pend(1:ppm_dim)   = MIN(patch(ppm_dim+1:2*ppm_dim),&
                  &   topo%max_subd(1:ppm_dim,isub))-Offset(1:ppm_dim)

                  !----------------------------------------------------------------
                  !Coordinates on the grid
                  !----------------------------------------------------------------
                  istart(1:ppm_dim) = 1 + CEILING(pstart(1:ppm_dim)/h(1:ppm_dim))
                  iend(1:ppm_dim)   = 1 + FLOOR(  pend  (1:ppm_dim)/h(1:ppm_dim))

                  !----------------------------------------------------------------
                  !Intersect these coordinates with those of the subdomain
                  !(some nodes may have been removed from the subdomain, e.g. to
                  !implement some boundary conditions and/or avoid node duplication.
                  !----------------------------------------------------------------
                  istart(1:ppm_dim) = MAX(istart(1:ppm_dim),this%istart(1:ppm_dim,isub))
                  iend(1:ppm_dim)   = MIN(iend(1:ppm_dim),this%iend(1:ppm_dim,isub))

                  !----------------------------------------------------------------
                  ! Specify boundary conditions for the subpatch
                  !----------------------------------------------------------------
                  DO k=1,ppm_dim
                     IF (istart(k).EQ.istart_d(k) ) THEN
                        bc(2*k-1) = topo%bcdef(2*k-1)
                     ELSE
                        IF (istart(k).EQ.istart_p(k) ) THEN
                           IF (PRESENT(bcdef)) THEN
                              bc(2*k-1) = bcdef(2*k-1)
                           ELSE
                              bc(2*k-1) = ppm_param_bcdef_freespace
                           ENDIF
                        ELSE
                           bc(2*k-1) = -1
                        ENDIF
                     ENDIF
                     IF (iend(k).EQ.iend_d(k)) THEN
                        bc(2*k) = topo%bcdef(2*k)
                     ELSE
                        IF (iend(k).EQ.iend_p(k) ) THEN
                           IF (PRESENT(bcdef)) THEN
                              bc(2*k) = bcdef(2*k)
                           ELSE
                              bc(2*k) = ppm_param_bcdef_freespace
                           ENDIF
                        ELSE
                           bc(2*k) = -1
                        ENDIF
                     ENDIF
                     ! For periodic boundary conditions, subpatches that touch
                     ! the East, North, or Top domain boundary are
                     ! reduced by 1 mesh node so that real mesh nodes are not
                     ! duplicated.
                     IF (bc(2*k).EQ.ppm_param_bcdef_periodic) THEN
                        iend(k) = iend(k) -1
                     ENDIF
                     IF (iend(k)+1.EQ.iend_d(k)) THEN
                        bc(2*k) = topo%bcdef(2*k)
                     ENDIF
                  ENDDO
                  !----------------------------------------------------------------
                  !Check that the subpatch contains at least one mesh nodes
                  !Otherwise, exit loop
                  ! (a subpatch comprises all the mesh nodes that are WITHIN
                  !  the patch. If a patch overlap with a subdomain by less than
                  !  h, it could be that this overlap regions has actually
                  !  no mesh nodes)
                  !----------------------------------------------------------------
                  IF (.NOT.ALL(istart(1:ppm_dim).LE.iend(1:ppm_dim))) THEN
                     CYCLE sub
                  ENDIF

                  !----------------------------------------------------------------
                  ! create a new subpatch object
                  !----------------------------------------------------------------
                  ALLOCATE(ppm_t_subpatch::p,STAT=info)
                  or_fail_alloc("could not allocate ppm_t_subpatch pointer")

                  !----------------------------------------------------------------
                  ! determine ghostlayer size for this subpatch
                  ! (there is a ghostlayer if and only if the border of the subpatch
                  ! does not coincide with a border of the patch itself - in that
                  ! case, the width of the ghostlayer is truncated by the "mesh-wide"
                  ! ghostsize parameter)
                  !----------------------------------------------------------------
                  ghostsize(1) = MIN(istart(1)-istart_p(1),this%ghostsize(1))
                  ghostsize(2) = MIN(iend_p(1)-iend(1)    ,this%ghostsize(1))
                  ghostsize(3) = MIN(istart(2)-istart_p(2),this%ghostsize(2))
                  ghostsize(4) = MIN(iend_p(2)-iend(2)    ,this%ghostsize(2))
                  IF (ppm_dim.EQ.3) THEN
                     ghostsize(5) = MIN(istart(3)-istart_p(3),this%ghostsize(3))
                     ghostsize(6) = MIN(iend_p(3)-iend(3)    ,this%ghostsize(3))
                  ENDIF

                  CALL p%create(this,isub,istart,iend,pstart,pend,&
                  &    istart_p,iend_p,ghostsize,bc,info)
                  or_fail("could not create new subpatch")

                  !----------------------------------------------------------------
                  !add a pointer to this subpatch
                  !----------------------------------------------------------------
                  nsubpatchi = nsubpatchi + 1
                  ! yaser: subpach_by_sub should contain global index
                  ! so I would use isub instead of i
                  this%subpatch_by_sub(isub)%vec(nsubpatchi)%t => p
                  this%subpatch_by_sub(isub)%nsubpatch = nsubpatchi

                  SELECT CASE (indsub(isub))
                  CASE (1)
                     ! subdomains on this processor
                     nsubpatch = nsubpatch+1
                     A_p%subpatch(nsubpatch)%t => p

                     !----------------------------------------------------------------
                     ! put the subpatch object in the collection of subpatches on this mesh
                     !----------------------------------------------------------------
                     CALL this%subpatch%push(p,info)
                     or_fail("could not add new subpatch to mesh")
                  CASE(-1)
                     p=>NULL()
                  END SELECT
              ENDIF
          ENDDO sub

          CALL this%patch%push(A_p,info,id)
          or_fail("could not add new subpatch_ptr_array to mesh%patch")

          this%patch%vec(id)%t%nsubpatch = nsubpatch
          this%patch%vec(id)%t%patchid = pid

          !Increment the number of patches defined on this mesh
          this%npatch = this%npatch + 1

          !The ghost mesh nodes have not been computed
          this%ghost_initialized = .FALSE.

          iopt = ppm_param_dealloc
          CALL ppm_alloc(indsub,ldu,iopt,info)
          or_fail_dealloc("indsub")

          end_subroutine()
      END SUBROUTINE equi_mesh_def_patch


      SUBROUTINE equi_mesh_def_uniform(this,info,patchid)
          !!! Add a uniform patch to a mesh (the patch covers the whole computational
          !!! domain, effectively mimicking a normal usual mesh without patches).

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_equi_mesh), INTENT(INOUT) :: this
          INTEGER,                INTENT(  OUT) :: info
          !!! Returns status, 0 upon success
          INTEGER, OPTIONAL,      INTENT(IN   ) :: patchid
          !!! id of the (uniform) patch, if we want one.
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double),DIMENSION(2*ppm_dim) :: patch

          start_subroutine("mesh_def_uniform")

          !----------------------------------------------------------------
          ! Create a huge patch
          !----------------------------------------------------------------
          patch(1:ppm_dim)           = -HUGE(1._ppm_kind_double)
          patch(ppm_dim+1:2*ppm_dim) =  HUGE(1._ppm_kind_double)

          !----------------------------------------------------------------
          !and add it to the mesh (it will compute the intersection
          ! between this infinite patch and the (hopefully) finite
          ! computational domain)
          !----------------------------------------------------------------
          CALL this%def_patch(patch,info,patchid,infinite=.TRUE.)
          or_fail("failed to add patch")

          !TODO add some checks for the finiteness of the computational domain

          end_subroutine()
      END SUBROUTINE equi_mesh_def_uniform


      FUNCTION equi_mesh_new_subpatch_data_ptr(this,info) RESULT(sp)
          !!! returns a pointer to a new subpatch_data object

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_equi_mesh),      INTENT(IN   ) :: this
          !!! cartesian mesh object
          CLASS(ppm_t_subpatch_data_), POINTER       :: sp
          INTEGER,                     INTENT(  OUT) :: info
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

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_equi_mesh), INTENT(INOUT) :: this
          INTEGER,                INTENT(  OUT) :: info
          INTEGER, DIMENSION(:),  POINTER       :: fids

          INTEGER                         :: i,j
          CLASS(ppm_t_main_abstr),POINTER :: f => NULL()

          CHARACTER(LEN=ppm_char) :: caller = "equi_mesh_list_of_fields"

          info = 0

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

          !j=1
          !DO i=this%field_ptr%min_id,this%field_ptr%max_id
              !IF (ASSOCIATED(this%field_ptr%vec(i)%t)) THEN
                  !fids(j) = this%field_ptr%vec(i)%t%fieldID
                  !j=j+1
              !ENDIF
          !ENDDO
          j=1
          f => this%field_ptr%begin()
          DO WHILE(ASSOCIATED(f))
              SELECT TYPE(f)
              CLASS IS (ppm_t_field_)
                  fids(j) = f%ID
                  j=j+1
              END SELECT
              f => this%field_ptr%next()
          ENDDO

          RETURN

          end_function()
      END FUNCTION


      SUBROUTINE equi_mesh_map_push(this,field,info)
          !!! Push field data onto the mesh mappings buffers

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_equi_mesh), INTENT(INOUT) :: this
          CLASS(ppm_t_field_),    INTENT(IN   ) :: field
          !!! this field is discretized on that mesh
          INTEGER,                INTENT(OUT) :: info

          REAL(ppm_kind_double), DIMENSION(:,:),     POINTER :: wp2_dummy => NULL()
          REAL(ppm_kind_double), DIMENSION(:,:,:),   POINTER :: wp3_dummy => NULL()
          REAL(ppm_kind_double), DIMENSION(:,:,:,:), POINTER :: wp4_dummy => NULL()

          INTEGER :: p_idx

          start_subroutine("equi_mesh_map_push")

          !p_idx = field%M%vec(this%ID)%t%p_idx
          p_idx = field%get_pid(this)

          SELECT CASE (ppm_dim)
          CASE (2)
              IF (field%lda.EQ.1) THEN
                  CALL mesh_map_push_2d_sca_d(this,wp2_dummy,p_idx,info)
                  or_fail("map_field_push_2d")
              ELSE
                  CALL mesh_map_push_2d_vec_d(this,wp3_dummy,field%lda,p_idx,info)
                  or_fail("map_field_push_2d")
              ENDIF
          CASE DEFAULT
              IF (field%lda.EQ.1) THEN
                  CALL mesh_map_push_3d_sca_d(this,wp3_dummy,p_idx,info)
                  or_fail("map_field_push_3d")
              ELSE
                  CALL mesh_map_push_3d_vec_d(this,wp4_dummy,field%lda,p_idx,info)
                  or_fail("map_field_push_3d")
              ENDIF
          END SELECT

          end_subroutine()
      END SUBROUTINE equi_mesh_map_push


      SUBROUTINE equi_mesh_map_pop(this,field,info)
          !!! Push field data onto the mesh mappings buffers

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_equi_mesh), INTENT(INOUT) :: this
          CLASS(ppm_t_field_),    INTENT(INOUT) :: field
          !!! this mesh is discretized on that field
          INTEGER,                INTENT(  OUT) :: info

          REAL(ppm_kind_double), DIMENSION(:,:),     POINTER :: wp2_dummy => NULL()
          REAL(ppm_kind_double), DIMENSION(:,:,:),   POINTER :: wp3_dummy => NULL()
          REAL(ppm_kind_double), DIMENSION(:,:,:,:), POINTER :: wp4_dummy => NULL()

          INTEGER :: p_idx

          start_subroutine("equi_mesh_map_pop")

          !p_idx = field%M%vec(this%ID)%t%p_idx
          p_idx = field%get_pid(this)

          IF (ppm_dim.EQ.2) THEN
              IF (field%lda.EQ.1) THEN
                  CALL mesh_map_pop_2d_sca_d(this,wp2_dummy,p_idx,info)
                  or_fail("map_field_pop_2d")
              ELSE
                  CALL mesh_map_pop_2d_vec_d(this,wp3_dummy,field%lda,p_idx,info)
                  or_fail("map_field_pop_2d")
              ENDIF
          ELSE
              IF (field%lda.EQ.1) THEN
                  CALL mesh_map_pop_3d_sca_d(this,wp3_dummy,p_idx,info)
                  or_fail("map_field_pop_3d")
              ELSE
                  CALL mesh_map_pop_3d_vec_d(this,wp4_dummy,field%lda,p_idx,info)
                  or_fail("map_field_pop_3d")
              ENDIF
          ENDIF

          end_subroutine()
      END SUBROUTINE equi_mesh_map_pop

      SUBROUTINE mesh_discr_data_create(this,field,info)

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_mesh_discr_data), INTENT(INOUT) :: this
          CLASS(ppm_t_field_),TARGET,   INTENT(IN   ) :: field
          INTEGER,                      INTENT(  OUT) :: info

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

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_mesh_discr_data), INTENT(INOUT) :: this
          INTEGER,                      INTENT(  OUT) :: info
          start_subroutine("mesh_discr_data_destroy")

          this%data_type = 0
          this%field_ptr => NULL()
          this%lda = 0
          this%name = ""

          end_subroutine()
      END SUBROUTINE

      SUBROUTINE equi_mesh_print_vtk(this,filename,info)

          USE ppm_module_io_vtk

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------

          CLASS(ppm_t_equi_mesh), INTENT(INOUT) :: this
          CHARACTER(LEN=*),       INTENT(IN   ) :: filename
          INTEGER,                INTENT(  OUT) :: info

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
#include "mesh/mesh_map_ghost_put.f"
#include "mesh/mesh_map_send.f"
#include "mesh/mesh_map_global.f"

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

#include "mesh/mesh_interp_to_part.f"

      END MODULE ppm_module_mesh_typedef

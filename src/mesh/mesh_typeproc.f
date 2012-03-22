!-------------------------------------------------------------------------
!  Subroutine   :                  mesh_typeproc
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
SUBROUTINE subpatch_get_field_3d_rd(this,wp,Field,info)
    !!! Constructor for subdomain data data structure
    CLASS(ppm_t_subpatch)                :: this
    CLASS(ppm_t_field_)                  :: Field
    REAL(ppm_kind_double),DIMENSION(:,:,:),POINTER :: wp
    INTEGER,                 INTENT(OUT) :: info

    !Direct access to the data arrays 
    ! TODO? different API?
    wp => this%subpatch_data%vec(Field%M%vec(this%meshID)%t%p_idx)%t%data_3d_rd


END SUBROUTINE

SUBROUTINE subpatch_get_field_2d_rd(this,wp,Field,info)
    !!! Constructor for subdomain data data structure
    CLASS(ppm_t_subpatch)                :: this
    CLASS(ppm_t_field_)                  :: Field
    REAL(ppm_kind_double),DIMENSION(:,:),POINTER :: wp
    INTEGER,                 INTENT(OUT) :: info

    !Direct access to the data arrays 
    ! TODO? different API?
    wp => this%subpatch_data%vec(Field%M%vec(this%meshID)%t%p_idx)%t%data_2d_rd


END SUBROUTINE

SUBROUTINE subpatch_data_create(pdata,datatype,lda,Nmp,info)
    !!! Constructor for subdomain data data structure
    CLASS(ppm_t_subpatch_data)              :: pdata
    INTEGER,                     INTENT(IN) :: datatype
    INTEGER,                     INTENT(IN) :: lda
    !!! number of data components per mesh node
    INTEGER,DIMENSION(:),POINTER,INTENT(IN) :: Nmp
    !!! number of mesh nodes in each dimension on this patch
    INTEGER,                    INTENT(OUT) :: info

    REAL(KIND(1.D0))                   :: t0
    INTEGER                            :: iopt, ndim
    CHARACTER(LEN=ppm_char)            :: caller = 'subpatch_data_create'

    CALL substart(caller,t0,info)


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
            CALL ppm_alloc(pdata%data_2d_i,ldc,iopt,info)
        CASE (ppm_type_real_single)
            CALL ppm_alloc(pdata%data_2d_rs,ldc,iopt,info)
        CASE (ppm_type_real_double)
            CALL ppm_alloc(pdata%data_2d_rd,ldc,iopt,info)
        CASE (ppm_type_comp_single)
            CALL ppm_alloc(pdata%data_2d_cs,ldc,iopt,info)
        CASE (ppm_type_comp_double)
            CALL ppm_alloc(pdata%data_2d_cd,ldc,iopt,info)
        CASE (ppm_type_logical)
            CALL ppm_alloc(pdata%data_2d_l,ldc,iopt,info)
        CASE DEFAULT
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,caller,   &
                &        'invalid type for mesh patch data',__LINE__,info)
        END SELECT
    CASE (3)
        SELECT CASE (datatype)
        CASE (ppm_type_int)
            CALL ppm_alloc(pdata%data_3d_i,ldc,iopt,info)
        CASE (ppm_type_real_single)
            CALL ppm_alloc(pdata%data_3d_rs,ldc,iopt,info)
        CASE (ppm_type_real_double)
            CALL ppm_alloc(pdata%data_3d_rd,ldc,iopt,info)
        CASE (ppm_type_comp_single)
            CALL ppm_alloc(pdata%data_3d_cs,ldc,iopt,info)
        CASE (ppm_type_comp_double)
            CALL ppm_alloc(pdata%data_3d_cd,ldc,iopt,info)
        CASE (ppm_type_logical)
            CALL ppm_alloc(pdata%data_3d_l,ldc,iopt,info)
        CASE DEFAULT
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,caller,   &
                &        'invalid type for mesh patch data',__LINE__,info)
        END SELECT
    CASE (4)
        SELECT CASE (datatype)
        CASE (ppm_type_int)
            CALL ppm_alloc(pdata%data_4d_i,ldc,iopt,info)
        CASE (ppm_type_real_single)
            CALL ppm_alloc(pdata%data_4d_rs,ldc,iopt,info)
        CASE (ppm_type_real_double)
            CALL ppm_alloc(pdata%data_4d_rd,ldc,iopt,info)
        CASE (ppm_type_comp_single)
            CALL ppm_alloc(pdata%data_4d_cs,ldc,iopt,info)
        CASE (ppm_type_comp_double)
            CALL ppm_alloc(pdata%data_4d_cd,ldc,iopt,info)
        CASE (ppm_type_logical)
            CALL ppm_alloc(pdata%data_4d_l,ldc,iopt,info)
        CASE DEFAULT
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,caller,   &
                &        'invalid type for mesh patch data',__LINE__,info)
        END SELECT

    END SELECT

    or_fail_alloc('allocating mesh patch data failed')

    CALL substop(caller,t0,info)

    9999  CONTINUE

END SUBROUTINE subpatch_data_create
!DESTROY
SUBROUTINE subpatch_data_destroy(pdata,info)
    !!! Destructor for subdomain data data structure
    CLASS(ppm_t_subpatch_data)     :: pdata
    INTEGER,               INTENT(OUT) :: info

    REAL(KIND(1.D0))                   :: t0
    INTEGER                            :: iopt
    CHARACTER(LEN=ppm_char)            :: caller = 'patch_data_destroy'

    CALL substart(caller,t0,info)

    iopt = ppm_param_dealloc
    CALL ppm_alloc(pdata%data_2d_i,ldc,iopt,info)
    CALL ppm_alloc(pdata%data_3d_i,ldc,iopt,info)
    CALL ppm_alloc(pdata%data_4d_i,ldc,iopt,info)

    CALL ppm_alloc(pdata%data_2d_l,ldc,iopt,info)
    CALL ppm_alloc(pdata%data_3d_l,ldc,iopt,info)
    CALL ppm_alloc(pdata%data_4d_l,ldc,iopt,info)

    CALL ppm_alloc(pdata%data_2d_rs,ldc,iopt,info)
    CALL ppm_alloc(pdata%data_3d_rs,ldc,iopt,info)
    CALL ppm_alloc(pdata%data_4d_rs,ldc,iopt,info)
    CALL ppm_alloc(pdata%data_2d_rd,ldc,iopt,info)
    CALL ppm_alloc(pdata%data_3d_rd,ldc,iopt,info)
    CALL ppm_alloc(pdata%data_4d_rd,ldc,iopt,info)

    CALL ppm_alloc(pdata%data_2d_cs,ldc,iopt,info)
    CALL ppm_alloc(pdata%data_3d_cs,ldc,iopt,info)
    CALL ppm_alloc(pdata%data_4d_cs,ldc,iopt,info)
    CALL ppm_alloc(pdata%data_2d_cd,ldc,iopt,info)
    CALL ppm_alloc(pdata%data_3d_cd,ldc,iopt,info)
    CALL ppm_alloc(pdata%data_4d_cd,ldc,iopt,info)

    CALL substop(caller,t0,info)
    9999  CONTINUE
END SUBROUTINE subpatch_data_destroy

!CREATE
SUBROUTINE subpatch_create(p,meshid,istart,iend,info)
    !!! Destructor for subdomain data data structure
    CLASS(ppm_t_subpatch)              :: p
    INTEGER                            :: meshid
    INTEGER,DIMENSION(:)               :: istart
    INTEGER,DIMENSION(:)               :: iend
    INTEGER,               INTENT(OUT) :: info

    REAL(KIND(1.D0))                   :: t0
    INTEGER                            :: iopt
    CHARACTER(LEN=ppm_char)            :: caller = 'subpatch_create'

    CALL substart(caller,t0,info)

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

    p%meshID=meshid
    p%istart = istart
    p%iend   = iend
    p%nnodes(1:ppm_dim) = 1 + iend(1:ppm_dim) - istart(1:ppm_dim)
    IF (.NOT.ASSOCIATED(p%subpatch_data)) THEN
        ALLOCATE(ppm_c_subpatch_data::p%subpatch_data,STAT=info)
        or_fail_alloc("could not allocate p%subpatch_data")
    ENDIF

    CALL substop(caller,t0,info)
    9999  CONTINUE
END SUBROUTINE subpatch_create

!DESTROY
SUBROUTINE subpatch_destroy(p,info)
    !!! Destructor for subdomain data data structure
    CLASS(ppm_t_subpatch)              :: p
    INTEGER,               INTENT(OUT) :: info

    REAL(KIND(1.D0))                   :: t0
    INTEGER                            :: iopt
    CHARACTER(LEN=ppm_char)            :: caller = 'subpatch_destroy'

    CALL substart(caller,t0,info)

    iopt = ppm_param_dealloc
    CALL ppm_alloc(p%nnodes,ldc,iopt,info)
    or_fail_dealloc("p%nnodes")
    CALL ppm_alloc(p%istart,ldc,iopt,info)
    or_fail_dealloc("p%istart")
    CALL ppm_alloc(p%iend,ldc,iopt,info)
    or_fail_dealloc("p%iend")
    IF (ASSOCIATED(p%subpatch_data)) THEN
        CALL p%subpatch_data%destroy(info)
        or_fail_dealloc("p%subpatch_data")
    ENDIF

    p%meshID=0

    CALL substop(caller,t0,info)
    9999  CONTINUE
END SUBROUTINE subpatch_destroy

!CREATE
SUBROUTINE subpatch_A_create(this,vecsize,info,patchid)
    !!! Destructor for subdomain data data structure
    CLASS(ppm_t_A_subpatch)            :: this
    INTEGER                            :: vecsize
    INTEGER,               INTENT(OUT) :: info
    INTEGER,OPTIONAL,      INTENT(IN)  :: patchid

    REAL(KIND(1.D0))                   :: t0
    CHARACTER(LEN=ppm_char)            :: caller = 'subpatch_A_create'

    CALL substart(caller,t0,info)

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

    CALL substop(caller,t0,info)
    9999  CONTINUE
END SUBROUTINE subpatch_A_create

!DESTROY
SUBROUTINE subpatch_A_destroy(this,info)
    !!! Destructor for subdomain data data structure
    CLASS(ppm_t_A_subpatch)            :: this
    INTEGER,               INTENT(OUT) :: info

    REAL(KIND(1.D0))                   :: t0
    CHARACTER(LEN=ppm_char)            :: caller = 'subpatch_A_destroy'

    CALL substart(caller,t0,info)

    this%patchid = 0
    this%nsubpatch = 0
    IF (ASSOCIATED(this%subpatch)) THEN
        DEALLOCATE(this%subpatch,STAT=info)
        or_fail_dealloc("could not deallocate this%subpatch")
    ENDIF
    NULLIFY(this%subpatch)

    CALL substop(caller,t0,info)
    9999  CONTINUE
END SUBROUTINE subpatch_A_destroy


SUBROUTINE equi_mesh_add_patch(this,patch,info,patchid)
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
    REAL(KIND(1.D0))                   :: t0
    TYPE(ppm_t_topo), POINTER :: topo => NULL()
    CHARACTER(LEN=ppm_char)   :: caller = 'mesh_add_patch'
    INTEGER                   :: i,meshid,isub,nsubpatch,id
    CLASS(ppm_t_subpatch_),  POINTER :: p => NULL()
    CLASS(ppm_t_A_subpatch_),POINTER :: A_p => NULL()
    INTEGER,              DIMENSION(ppm_dim) :: istart,iend
    REAL(ppm_kind_double),DIMENSION(ppm_dim) :: h,Offset
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)

    !-------------------------------------------------------------------------
    !  Check arguments
    !-------------------------------------------------------------------------
    !-------------------------------------------------------------------------
    !  get mesh parameters
    !-------------------------------------------------------------------------
    h = this%h
    Offset = this%Offset
    meshid = this%ID
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

            !SELECT TYPE(pp => p)
            !TYPE IS (ppm_t_subpatch)
                CALL p%create(meshid,istart,iend,info)
                or_fail("could not create new subpatch")
            !END SELECT

            nsubpatch = nsubpatch+1
            A_p%subpatch(nsubpatch)%t => p

            ! add it to the list of subpatches on this mesh
            CALL this%subpatch%push(p,info,id)
            or_fail("could not add new subpatch to mesh")

        ENDIF
    ENDDO
    CALL this%patch%push(A_p,info,id)
    this%patch%vec(id)%t%nsubpatch = nsubpatch

    9999 CONTINUE
    CALL substop(caller,t0,info)
END SUBROUTINE equi_mesh_add_patch

SUBROUTINE equi_mesh_create(this,topoid,Offset,info,Nm,h)
    !!! Creator for the cartesian mesh object
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
    REAL(ppm_kind_double)     :: t0
    LOGICAL                   :: valid
    TYPE(ppm_t_topo), POINTER :: topo => NULL()
    CHARACTER(LEN=ppm_char)   :: caller = 'equi_mesh_create'
    !-------------------------------------------------------------------------
    !  Externals
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)

    !-------------------------------------------------------------------------
    !  Check arguments
    !-------------------------------------------------------------------------
    IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
    ENDIF

    !This mesh is defined for a given topology
    this%topoid = topoid

    !Global id, internal to ppm
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
    9999 CONTINUE
    CALL substop(caller,t0,info)
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
    REAL(ppm_kind_double)     :: t0
    LOGICAL                   :: valid
    TYPE(ppm_t_topo), POINTER :: topo => NULL()
    CHARACTER(LEN=ppm_char)   :: caller = 'equi_mesh_destroy'
    !-------------------------------------------------------------------------
    !  Externals
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)

    !-------------------------------------------------------------------------
    !  (Re)allocate memory for the internal mesh list and Arrays at meshid
    !-------------------------------------------------------------------------
    ldc = 1
    iopt   = ppm_param_dealloc
    CALL ppm_alloc(this%Offset,ldc,iopt,info)
    or_fail_dealloc('Offset')
    CALL ppm_alloc(this%Nm,ldc,iopt,info)
    or_fail_dealloc('Nm')
    CALL ppm_alloc(this%h,ldc,iopt,info)
    or_fail_dealloc('h')

    IF (ASSOCIATED(this%subpatch)) THEN
        CALL this%subpatch%destroy(info)
        or_fail_dealloc('subpatch object')
        DEALLOCATE(this%subpatch,STAT=info)
        or_fail_dealloc('subpatch')
    ENDIF
    IF (ASSOCIATED(this%patch)) THEN
        CALL this%patch%destroy(info)
        or_fail_dealloc('patch object')
        DEALLOCATE(this%patch,STAT=info)
        or_fail_dealloc('patch')
    ENDIF
    IF (ASSOCIATED(this%sub)) THEN
        CALL this%sub%destroy(info)
        or_fail_dealloc('sub object')
        DEALLOCATE(this%sub,STAT=info)
        or_fail_dealloc('sub')
    ENDIF

    !-------------------------------------------------------------------------
    !  Return
    !-------------------------------------------------------------------------
    9999 CONTINUE
    CALL substop(caller,t0,info)
    RETURN
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
    REAL(ppm_kind_double)     :: t0
    CHARACTER(LEN=ppm_char)   :: caller = 'equi_mesh_new_subpatch_data_ptr'
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)

    ALLOCATE(ppm_t_subpatch_data::sp,STAT=info)
    or_fail_alloc("could not allocate ppm_t_subpatch_data pointer")

    !-------------------------------------------------------------------------
    !  Return
    !-------------------------------------------------------------------------
    9999 CONTINUE
    CALL substop(caller,t0,info)
    RETURN
END FUNCTION equi_mesh_new_subpatch_data_ptr

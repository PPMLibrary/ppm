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
    INTEGER :: info

    !Direct access to the data arrays 
    ! TODO? different API?
    wp => this%subpatch_data%vec(Field%M%vec(this%meshID)%p_idx)%data_3d_rd


END SUBROUTINE

SUBROUTINE subpatch_get_field_2d_rd(this,wp,Field,info)
    !!! Constructor for subdomain data data structure
    CLASS(ppm_t_subpatch)                :: this
    CLASS(ppm_t_field_)                  :: Field
    REAL(ppm_kind_double),DIMENSION(:,:),POINTER :: wp
    INTEGER :: info

    !Direct access to the data arrays 
    ! TODO? different API?
    wp => this%subpatch_data%vec(Field%M%vec(this%meshID)%p_idx)%data_2d_rd


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
    CHARACTER(LEN=ppm_char)            :: caller = 'sub_create'
    INTEGER,DIMENSION(:),POINTER :: istart,iend



    CALL substart(caller,t0,info)


    iopt   = ppm_param_alloc_grow

    ! Determine allocation size of the data array
    IF (MINVAL(Nmp(1:ppm_dim)) .LE. 0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_argument,caller,   &
            & 'invalid size for patch data. This patch should be deleted',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    ldc(1:ppm_dim) = Nmp(1:ppm_dim)

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

    p%meshID=meshid
    p%istart = istart
    p%iend = iend
    ALLOCATE(ppm_c_subpatch_data::p%subpatch_data)

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
    CALL ppm_alloc(p%istart,ldc,iopt,info)
    or_fail_dealloc("p%istart")
    CALL ppm_alloc(p%iend,ldc,iopt,info)
    or_fail_dealloc("p%iend")
    CALL p%subpatch_data%destroy(info)
    or_fail_dealloc("p%subpatch_data")

    p%meshID=0

    CALL substop(caller,t0,info)
    9999  CONTINUE
END SUBROUTINE subpatch_destroy

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
    INTEGER                   :: i,meshid,isub
    CLASS(ppm_t_subpatch_),POINTER :: p => NULL()
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

    !-------------------------------------------------------------------------
    !  intersect each subdomain with the patch and store the corresponding
    !  subpatch in the mesh data structure
    !-------------------------------------------------------------------------
    topo => ppm_topo(this%topoid)%t
    ! loop through all subdomains on this processor
    DO i = 1,topo%nsublist
        isub = topo%isublist(i)
        ! check if the subdomain overlaps with the patch
        IF (ALL(patch(1:ppm_dim).LT.topo%max_subd(1:ppm_dim,isub)) .AND. &
            ALL(patch(ppm_dim+1:2*ppm_dim).GT.topo%min_subd(1:ppm_dim,isub)))&
               THEN 

            ! finds the mesh nodes which are contained in the overlap region
            istart(1:ppm_dim) = CEILING((MAX(patch(1:ppm_dim),&
                topo%min_subd(1:ppm_dim,isub))-Offset(1:ppm_dim))/h(1:ppm_dim))
            iend(1:ppm_dim)= FLOOR((MIN(patch(ppm_dim+1:2*ppm_dim),&
                topo%max_subd(ppm_dim+1:2*ppm_dim,isub))-Offset(1:ppm_dim))&
                /h(1:ppm_dim))

            ! create a new subpatch object
            CALL p%create(meshid,istart,iend,info)
            or_fail("could not create new subpatch")

            ! add it to the list of subpatches on this mesh
            CALL this%subpatch%push(p,info)
            or_fail("could not add new subpatch to mesh")

        ENDIF
    ENDDO

    9999 CONTINUE
    CALL substop(caller,t0,info)
END SUBROUTINE equi_mesh_add_patch

SUBROUTINE equi_mesh_create(this,topoid,Nm,Offset,info)
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
    INTEGER , DIMENSION(:),   INTENT(IN   ) :: Nm
    !!! Global number of mesh points in the whole comput. domain
    REAL(ppm_kind_double), DIMENSION(:), INTENT(IN   ) :: Offset
    !!! Offset in each dimension
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

    this%topoid = topoid

    topo => ppm_topo(topoid)%t

    !macro_test( 1 )
    !macro_test( "1" )
    !macro_test( "(1)" )
    !macro_test( 1+2 )
    !macro_test( (1) )
    !macro_test( 1/(1+2) )

    !check_equal(a,b,"this is an error message")

    !check_equal(SIZE(Nm,1),ppm_dim,"invalid size for Nm")

    IF (SIZE(Nm,1) .NE. ppm_dim) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,caller,  &
            &       'invalid size for Nm. Should be ppm_dim',__LINE__,info)
        GOTO 9999
    ENDIF

    IF (SIZE(Offset,1) .NE. ppm_dim) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,caller,  &
            &       'invalid size for Offset. Should be ppm_dim',__LINE__,info)
        GOTO 9999
    ENDIF

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

    !-------------------------------------------------------------------------
    !  Store the mesh information
    !-------------------------------------------------------------------------

    this%Nm(1:ppm_dim) = Nm(1:ppm_dim)
    this%Offset(1:ppm_dim) = Offset(1:ppm_dim)


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
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,  &
                &            'topoid not valid',__LINE__,info)
            GOTO 8888
        ENDIF
        IF (ppm_topo(topoid)%t%nsubs .LE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,  &
                &            'nsubs must be >0',__LINE__,info)
            GOTO 8888
        ENDIF
        DO i=1,ppm_dim
            IF (Nm(i) .LT. 2) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,caller,  &
                    &                'Nm must be >1 in all space dimensions',__LINE__,info)
                GOTO 8888
            ENDIF
        ENDDO
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

    CALL this%subpatch%destroy(info)
    or_fail_dealloc('subpatch')

    !TODO !!!!!
!    CALL this%patch%destroy(info)
!    or_fail_dealloc('patch')
!    ENDIF
!    CALL this%sub%destroy(info)
!    or_fail_dealloc('sub')

    !-------------------------------------------------------------------------
    !  Return
    !-------------------------------------------------------------------------
    9999 CONTINUE
    CALL substop(caller,t0,info)
    RETURN
END SUBROUTINE equi_mesh_destroy

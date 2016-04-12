test_suite ppm_module_field_typedef

  USE ppm_module_mesh_typedef
  USE ppm_module_topo_typedef
  USE ppm_module_data
  USE ppm_module_mktopo
  USE ppm_module_finalize
  USE ppm_module_interfaces

  INTEGER, PARAMETER              :: debug = 0
  INTEGER, PARAMETER              :: MK = KIND(1.0d0) !KIND(1.0e0)
#ifdef __MPI
  INTEGER, PARAMETER              :: comm = MPI_COMM_WORLD
#endif
  INTEGER                         :: ndim
  INTEGER                         :: nspec
  INTEGER                         :: rank
  INTEGER                         :: nproc
  INTEGER                         :: decomp
  INTEGER                         :: assig
  INTEGER                         :: tolexp
  REAL(MK)                        :: tol
  INTEGER                         :: info
  INTEGER                         :: topoid=-1
  REAL(MK),DIMENSION(:,:),POINTER :: xp => NULL()
  REAL(MK),DIMENSION(:,:),POINTER :: wp => NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: min_phys => NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: max_phys => NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: h => NULL()
  INTEGER, DIMENSION(:),  POINTER :: ighostsize => NULL()
  REAL(MK), DIMENSION(:), POINTER :: ghostsize => NULL()
  INTEGER                         :: i,j,k,p_i,ai,aj,it,isub
  INTEGER, DIMENSION(:  ),POINTER :: bcdef => NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: cost => NULL()
  INTEGER, DIMENSION(:,:),POINTER :: istart => NULL()
  INTEGER, DIMENSION(:,:),POINTER :: ndata => NULL()
  INTEGER, DIMENSION(:  ),POINTER :: nm => NULL()
  INTEGER                         :: np,mp
  INTEGER                         :: kernel
  INTEGER                         :: seedsize
  INTEGER,  DIMENSION(:),ALLOCATABLE :: seed
  REAL(MK),DIMENSION(:,:),  POINTER:: field2d_1=>NULL(),field2d_2=>NULL()
  REAL(MK),DIMENSION(:,:,:),POINTER:: field3d_1=>NULL(),field3d_2=>NULL()
  REAL(MK),DIMENSION(:,:,:),POINTER:: field4d_1=>NULL(),field4d_2=>NULL()
  INTEGER, DIMENSION(2)            :: maxndata
  INTEGER, DIMENSION(:  ), POINTER :: isublist => NULL()
  INTEGER                          :: nsublist
  INTEGER                          :: ipatch,mypatchid
  TYPE(ppm_t_topo),      POINTER   :: topo => NULL()
  REAL(MK)                         :: sca_ghostsize


  !-------------------------- init testsuit -------------------------------------
  init
    USE ppm_module_data
    USE ppm_module_init

    tol = 100.0_MK*EPSILON(1.0_MK)
    tolexp = INT(LOG10(EPSILON(1.0_MK)))
    ndim = 2
    nspec = 1
    ALLOCATE(bcdef(2*ndim))
    ALLOCATE(min_phys(ndim),max_phys(ndim),ghostsize(ndim),&
    &    ighostsize(ndim),nm(ndim),h(ndim),STAT=info)


#ifdef __MPI
    CALL MPI_Comm_rank(comm,rank,info)
    CALL MPI_Comm_size(comm,nproc,info)
#else
    nproc = 1
    rank = 0
#endif

    CALL ppm_init(ndim,mk,tolexp,0,debug,info,99)

  end init
  !------------------------------------------------------------------------------

  !------------------------- finalize testsuit ----------------------------------
  finalize
    USE ppm_module_finalize

    CALL ppm_finalize(info)

    DEALLOCATE(min_phys,max_phys,ghostsize,nm)
  end finalize
  !------------------------------------------------------------------------------

  !------------------------------ test setup ------------------------------------
  setup
    NULLIFY(xp)
    NULLIFY(wp)
    NULLIFY(cost)

    np = 400*nproc
    mp = 0

    CALL RANDOM_SEED(SIZE=seedsize)
    ALLOCATE(seed(seedsize))
    DO i=1,seedsize
        seed(i)=10+i*i*(rank+1)
    ENDDO
    CALL RANDOM_SEED(PUT=seed)

    bcdef(1:2*ndim) = ppm_param_bcdef_freespace
    kernel = ppm_param_rmsh_kernel_mp4
    DO i=1,ndim
       min_phys(i) = 0.0_MK
       max_phys(i) = 1.0_MK
       ighostsize(i) = 2
       ghostsize(i) = 0.05_MK
    ENDDO
  end setup
  !------------------------------------------------------------------------------


  !--------------------------- test teardown ------------------------------------
  teardown
    IF (ASSOCIATED(xp)) DEALLOCATE(xp)
    IF (ASSOCIATED(wp)) DEALLOCATE(wp)
    IF (ASSOCIATED(cost)) DEALLOCATE(cost)
    IF (ALLOCATED(seed)) DEALLOCATE(seed)
  end teardown
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  test field_uniform_basics
    REAL(MK),DIMENSION(ndim)        :: offset
    TYPE(ppm_t_field)               :: Vort,Veloc
    TYPE(ppm_t_equi_mesh),  TARGET  :: Mesh1,Mesh2
    CLASS(ppm_t_subpatch_), POINTER :: p => NULL()
    INTEGER,DIMENSION(:),   POINTER :: ilist => NULL()

    start_subroutine("field_uniform_basics")
    !----------------
    ! make topology
    !----------------
    decomp = ppm_param_decomp_cuboid
    assig  = ppm_param_assign_internal
    topoid = 0

    ALLOCATE(nm(ndim),STAT=info)
    nm(1:ndim) = 16*nproc

    sca_ghostsize = 0.05_MK

    CALL ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,sca_ghostsize,cost,info)
    Assert_Equal(info,0)

    !--------------------------
    !Define Fields
    !--------------------------
    CALL Vort%create(1,info,name="Vorticity") !scalar field
    Assert_Equal(info,0)
    CALL Veloc%create(ndim,info,name="Velocity") !vector field
    Assert_Equal(info,0)

    !--------------------------
    !Create Mesh
    !--------------------------
    offset = 0.0_MK
    CALL Mesh1%create(topoid,offset,info,Nm=Nm)
    Assert_Equal(info,0)

    CALL Mesh1%def_uniform(info)
    Assert_Equal(info,0)
    Assert_True(ASSOCIATED(Mesh1%subpatch))

    p => Mesh1%subpatch%begin()
    DO WHILE(ASSOCIATED(p))
       Assert_True(ASSOCIATED(p%subpatch_data))
       p => Mesh1%subpatch%next()
    ENDDO

    !--------------------------
    !Create data arrays on the mesh for the vorticity and velocity fields
    !--------------------------
    CALL Vort%discretize_on(Mesh1,info)
    Assert_Equal(info,0)

    CALL Veloc%discretize_on(Mesh1,info)
    Assert_Equal(info,0)

    !--------------------------
    ! Iterate through patches and initialize the data arrays
    !--------------------------
    p => Mesh1%subpatch%begin()

    DO WHILE (ASSOCIATED(p))
       CALL p%get_field(Vort,field2d_1,info)
       CALL p%get_field(Veloc,field3d_1,info)

       DO i = 1,p%nnodes(1)
          DO j = 1,p%nnodes(2)
             field2d_1(i,j) = COS(i*h(1)+j)
             field3d_1(1,i,j) = SIN(field2d_1(i,j))
             field3d_1(2,i,j) = COS(field2d_1(i,j))
          ENDDO
       ENDDO
       p => Mesh1%subpatch%next()
    ENDDO

    !Second version
    DO ipatch = 1,Mesh1%subpatch%nb
       p => Mesh1%subpatch%vec(ipatch)%t
       CALL p%get_field(Vort,field2d_1,info)
       CALL p%get_field(Veloc,field3d_1,info)

       DO i = 1,p%nnodes(1)
          DO j = 1,p%nnodes(2)
             field2d_1(i,j) = COS(i*h(1)+j)
             field3d_1(1,i,j) = SIN(field2d_1(i,j))
             field3d_1(2,i,j) = COS(field2d_1(i,j))
          ENDDO
       ENDDO
    ENDDO

    field2d_1 => NULL()
    field3d_1 => NULL()


    ilist => Mesh1%list_of_fields(info)
    IF (ASSOCIATED(ilist)) THEN
       stdout("mesh1%list: ",ilist)
    ENDIF
    ilist => Mesh2%list_of_fields(info)
    IF (ASSOCIATED(ilist)) THEN
       stdout("mesh2%list: ",ilist)
    ENDIF

    !--------------------
    ! Remove a patch
    !--------------------
    !Mesh1%remove_patch(patch_ID = 2)
    CALL Mesh1%subpatch%destroy(info)
    Assert_Equal(info,0)

    CALL Vort%destroy(info)
    Assert_Equal(info,0)
    CALL Veloc%destroy(info)
    Assert_Equal(info,0)
    CALL Mesh1%destroy(info)
    Assert_Equal(info,0)
    CALL Mesh2%destroy(info)
    Assert_Equal(info,0)

    p=>NULL()

    end_subroutine()
  end test

  !------------------------------------------------------------------------------
  test field_basics

    IMPLICIT NONE
    CLASS(ppm_t_subpatch_), POINTER :: p => NULL()

    TYPE(ppm_t_field)             :: Vort,Veloc
    TYPE(ppm_t_equi_mesh), TARGET :: Mesh1,Mesh2

    REAL(MK),DIMENSION(4)    :: my_patch
    REAL(MK),DIMENSION(ndim) :: offset

    !----------------
    ! make topology
    !----------------
    decomp = ppm_param_decomp_cuboid
    assig  = ppm_param_assign_internal
    topoid = 0

    ALLOCATE(nm(ndim),STAT=info)
    nm(1:ndim) = 16*nproc

    sca_ghostsize = 0.05_MK

    CALL ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,sca_ghostsize,cost,info)
    Assert_Equal(info,0)

    !--------------------------
    !Define Fields
    !--------------------------
    CALL Vort%create(1,info,name="Vorticity")
    Assert_Equal(info,0)
    CALL Veloc%create(ndim,info,name="Velocity")
    Assert_Equal(info,0)

    !--------------------------
    !Create Mesh
    !--------------------------
    offset = 0.0_MK
    CALL Mesh1%create(topoid,offset,info,Nm=Nm)
    Assert_Equal(info,0)

    mypatchid = 1
    my_patch(1:4) = (/0.5,0.1,5.1,10.0/)

    CALL Mesh1%def_patch(my_patch,info,mypatchid)
    Assert_Equal(info,0)
    Assert_True(ASSOCIATED(Mesh1%subpatch))

    p => Mesh1%subpatch%begin()
    DO WHILE(ASSOCIATED(p))
       Assert_True(ASSOCIATED(p%subpatch_data))
       p => Mesh1%subpatch%next()
    ENDDO

    !--------------------------
    !Create data arrays on the mesh for the vorticity and velocity fields
    !--------------------------
    CALL Vort%discretize_on(Mesh1,info)
    Assert_Equal(info,0)

    CALL Veloc%discretize_on(Mesh1,info)
    Assert_Equal(info,0)

    !--------------------------
    ! Iterate through patches and initialize the data arrays
    !--------------------------
    p => Mesh1%subpatch%begin()
    DO WHILE (ASSOCIATED(p))
       CALL p%get_field(Vort,field2d_1,info)
       CALL p%get_field(Veloc,field3d_1,info)

       DO i = 1,p%nnodes(1)
          DO j = 1,p%nnodes(2)
             field2d_1(i,j) = COS(i*h(1)+j)
             field3d_1(1:ndim,i,j) = 17.4
          ENDDO
       ENDDO
       p => Mesh1%subpatch%next()
    ENDDO

    !Second version
    DO ipatch = 1,Mesh1%subpatch%nb
       p => Mesh1%subpatch%vec(ipatch)%t
       CALL p%get_field(Vort,field2d_1,info)
       CALL p%get_field(Veloc,field3d_1,info)

       DO i = 1,p%nnodes(1)
          DO j = 1,p%nnodes(2)
             field2d_1(i,j) = COS(i*h(1)+j)
             field3d_1(1:ndim,i,j) = 17.4
          ENDDO
       ENDDO
    ENDDO

    !--------------------
    ! Remove a patch
    !--------------------
    CALL Mesh1%subpatch%destroy(info)
    Assert_Equal(info,0)

    CALL Mesh1%destroy(info)
    Assert_Equal(info,0)

  end test
  !------------------------------------------------------------------------------

end test_suite

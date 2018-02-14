test_suite ppm_module_particles_big

  USE ppm_module_particles_typedef
  USE ppm_module_topo_typedef
  USE ppm_module_field_typedef
  USE ppm_module_operator_typedef
  USE ppm_module_sop
  USE ppm_module_interfaces
  USE ppm_module_data

  INTEGER, PARAMETER              :: debug = 0
  INTEGER, PARAMETER              :: MK = KIND(1.0d0) !KIND(1.0e0)
  REAL(MK),PARAMETER              :: tol=EPSILON(1.0_MK)*100
  REAL(MK),PARAMETER              :: pi = ACOS(-1.0_MK)
  REAL(MK),PARAMETER              :: skin = 0.0_MK
  INTEGER,PARAMETER               :: ndim=3
  INTEGER                         :: decomp,assig,tolexp
  INTEGER                         :: info,comm,rank,nproc,topoid
  INTEGER                         :: np_global = 3000
  REAL(MK),PARAMETER              :: cutoff = 0.15_MK
  REAL(MK),DIMENSION(:,:),POINTER :: xp=>NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: min_phys=>NULL(),max_phys=>NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: len_phys=>NULL()
  INTEGER                         :: i,j,k,ip,wp_id
  INTEGER                         :: nstep
  INTEGER,DIMENSION(3)            :: ldc
  INTEGER, DIMENSION(6)           :: bcdef
  REAL(MK),DIMENSION(:  ),POINTER :: cost=>NULL()
  INTEGER                         :: isymm = 0
  REAL(MK)                        :: t0,t1,t2,t3
  INTEGER                         :: seedsize
  INTEGER,  DIMENSION(:),ALLOCATABLE :: seed
  INTEGER, DIMENSION(:),POINTER   :: nvlist=>NULL()
  INTEGER, DIMENSION(:,:),POINTER :: vlist=>NULL()
  REAL(MK)                        :: err

  INTEGER, DIMENSION(:), POINTER                 :: wp_1i => NULL()
  INTEGER, DIMENSION(:,:), POINTER               :: wp_2i => NULL()
  INTEGER(ppm_kind_int64),DIMENSION(:),  POINTER :: wp_1li => NULL()
  INTEGER(ppm_kind_int64),DIMENSION(:,:),POINTER :: wp_2li => NULL()
  REAL(MK), DIMENSION(:),   POINTER              :: wp_1r => NULL()
  REAL(MK), DIMENSION(:,:), POINTER              :: wp_2r => NULL()
  COMPLEX(MK), DIMENSION(:),   POINTER           :: wp_1c => NULL()
  COMPLEX(MK), DIMENSION(:,:), POINTER           :: wp_2c => NULL()
  LOGICAL, DIMENSION(:),   POINTER               :: wp_1l => NULL()

  INTEGER, DIMENSION(:),ALLOCATABLE              :: degree,order
  REAL(ppm_kind_double),DIMENSION(:),ALLOCATABLE :: coeffs
  INTEGER                                        :: nterms

  init
    USE ppm_module_init
    USE ppm_module_mktopo

    ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),STAT=info)

    min_phys(1:ndim) = 0.0_MK
    max_phys(1:ndim) = 1.0_MK
    len_phys(1:ndim) = max_phys-min_phys
    bcdef(1:6) = ppm_param_bcdef_periodic

#ifdef __MPI
    comm = MPI_COMM_WORLD
    CALL MPI_Comm_rank(comm,rank,info)
    CALL MPI_Comm_size(comm,nproc,info)
#else
    rank = 0
    nproc = 1
#endif
    tolexp = INT(LOG10(EPSILON(1.0_MK)))+10
    CALL ppm_init(ndim,mk,tolexp,comm,debug,info,99)

    CALL RANDOM_SEED(SIZE=seedsize)
    ALLOCATE(seed(seedsize))
    DO i=1,seedsize
        seed(i)=10+i*i*(rank+1)
    ENDDO
    CALL RANDOM_SEED(PUT=seed)

    !----------------
    ! make topology
    !----------------
    decomp = ppm_param_decomp_cuboid
    !decomp = ppm_param_decomp_xpencil
    assig  = ppm_param_assign_internal

    topoid = 0

    CALL ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,cutoff,cost,info)
  end init

  finalize
    USE ppm_module_finalize

    CALL ppm_finalize(info)

    DEALLOCATE(min_phys,max_phys,len_phys)
  end finalize

  setup
  end setup

  teardown
  end teardown

  test initialize_random
    ! test initialization of particles sampled unif. random
    TYPE(ppm_t_particles_d),TARGET :: Part1

    CALL Part1%initialize(np_global,info,topoid=topoid,distrib=ppm_param_part_init_random)
    Assert_Equal(info,0)
    Assert_True(ASSOCIATED(Part1%xp))

    CALL Part1%get_xp(xp,info)
    Assert_Equal(info,0)
    CALL Part1%set_xp(xp,info,read_only=.TRUE.)
    Assert_Equal(info,0)

    CALL Part1%map(info,global=.TRUE.,topoid=topoid)
    Assert_Equal(info,0)
    CALL Part1%map_ghosts(info)
    Assert_Equal(info,0)

    CALL Part1%destroy(info)
    Assert_Equal(info,0)

  end test

  test initialize_cart
    ! test initialization of particles on a grid
    TYPE(ppm_t_particles_d),TARGET   :: Part1
    CLASS(ppm_t_discr_data),POINTER  :: Prop1=>NULL(),Prop2=>NULL()
    TYPE(ppm_t_field)                :: Field1


    CALL Part1%initialize(np_global,info,topoid=topoid)
    Assert_Equal(info,0)
    Assert_True(ASSOCIATED(Part1%xp))

    CALL Field1%create(3,info,name="Field1")
    Assert_Equal(info,0)
    CALL Field1%discretize_on(Part1,info)
    Assert_Equal(info,0)

    CALL Part1%get_xp(xp,info)
    Assert_Equal(info,0)
    CALL Part1%set_xp(xp,info,read_only=.TRUE.)
    Assert_Equal(info,0)

    CALL Part1%map(info,global=.TRUE.,topoid=topoid)
    Assert_Equal(info,0)

    CALL Part1%map_ghosts(info)
    Assert_Equal(info,0)

    !creating/destroying properties of different types
    CALL Part1%create_prop(info,discr_data=Prop1,dtype=ppm_type_longint,name='test_li',with_ghosts=.TRUE.)
    Assert_Equal(info,0)
    CALL Part1%get(Prop1,wp_1li,info,with_ghosts=.TRUE.)
    Assert_Equal(info,0)
    CALL Part1%destroy_prop(Prop1,info)
    Assert_Equal(info,0)

    CALL Part1%create_prop(info,discr_data=Prop1,dtype=ppm_type_int,name='test_i',with_ghosts=.TRUE.)
    Assert_Equal(info,0)
    CALL Part1%get(Prop1,wp_1i,info)
    Assert_Equal(info,0)

    CALL Part1%create_prop(info,discr_data=Prop1,dtype=ppm_type_int,name='test_2i',lda=10,zero=.TRUE.)
    Assert_Equal(info,0)
    CALL Part1%get(Prop1,wp_2i,info)
    Assert_Equal(MINVAL(wp_2i),0)
    Assert_Equal(info,0)
    Assert_Equal(MAXVAL(wp_2i),0)
    Assert_Equal(info,0)


    CALL Part1%create_prop(info,discr_data=Prop1,dtype=ppm_type_real,name='test_r',zero=.TRUE.)
    Assert_Equal(info,0)
    CALL Part1%get(Prop1,wp_1r,info)
    Assert_Equal(info,0)

    CALL Part1%create_prop(info,discr_data=Prop2,dtype=ppm_type_LOGICAL,name='test_l',zero=.TRUE.)
    Assert_Equal(info,0)
    CALL Part1%get(Prop2,wp_1l,info)
    Assert_Equal(info,0)

    DO i=1,25
        CALL Part1%create_prop(info,discr_data=Prop1,dtype=ppm_type_comp,lda=3,zero=.TRUE.,name="prop1")
        Assert_Equal(info,0)
        CALL Part1%create_prop(info,discr_data=Prop2,dtype=ppm_type_real,lda=1,name="prop2")
        Assert_Equal(info,0)
        CALL Part1%get(Prop1,wp_2c,info)
        Assert_Equal(info,0)

        CALL Part1%destroy_prop(Prop1,info)
        Assert_Equal(info,0)
        SELECT TYPE(Prop1)
        CLASS IS(ppm_t_part_prop_d)
        Assert_Equal(Part1%props%get_id(Prop1),-1)
        END SELECT
        SELECT TYPE(Prop2)
        CLASS IS(ppm_t_part_prop_d)
        Assert_True(Part1%props%get_id(Prop2).GT.0)
        wp_id = Part1%props%get_id(Prop2)
        Assert_Equal(Part1%props%vec(wp_id)%t%lda,1)
        END SELECT

        CALL Part1%destroy_prop(Prop2,info)
        Assert_Equal(info,0)

    ENDDO

    !Set up a velocity field and a scalar test function on the particles
    CALL Part1%create_prop(info,discr_data=Prop1,dtype=ppm_type_real,lda=3,name="velocity")
    Assert_Equal(info,0)
    CALL Part1%create_prop(info,discr_data=Prop2,dtype=ppm_type_real,lda=1,name="testf")
    Assert_Equal(info,0)
    CALL Part1%get(Prop1,wp_2r,info)
    Assert_Equal(info,0)
    CALL Part1%get(Prop2,wp_1r,info)
    Assert_Equal(info,0)
    CALL Part1%get_xp(xp,info)
    Assert_Equal(info,0)
    DO ip=1,Part1%Npart
        wp_2r(1:ndim,ip) = COS((10._MK*xp(1:ndim,ip))**2)
        wp_1r(ip) = f0_test(xp(1:ndim,ip),ndim)
    ENDDO
    wp_2r = COS(wp_2r) * Part1%ghostlayer

    !Move the particles with this displacement field
    CALL Part1%move(wp_2r,info)
    Assert_Equal(info,0)


    !Apply boundary conditions and remap the particles
    CALL Part1%apply_bc(info)
    Assert_Equal(info,0)
    CALL Part1%map(info)
    Assert_Equal(info,0)

    !Get the new ghosts
    CALL Part1%map_ghosts(info)
    Assert_Equal(info,0)

    !Move the particles back to their initial positions
    CALL Part1%get(Prop1,wp_2r,info)
    wp_2r = -wp_2r
    CALL Part1%move(wp_2r,info)
    Assert_Equal(info,0)

    !Re-apply boundary conditions and remap the particles
    CALL Part1%apply_bc(info)
    Assert_Equal(info,0)
    CALL Part1%map(info)
    Assert_Equal(info,0)

    CALL Part1%map_ghosts(info)
    Assert_Equal(info,0)

    !Compare values and check that they are still the same
    CALL Part1%get_xp(xp,info)
    Assert_Equal(info,0)
    CALL Part1%get(Prop2,wp_1r,info)
    Assert_Equal(info,0)
    err = 0.0_MK
    DO ip=1,Part1%Npart
        err = MAX(err,ABS(wp_1r(ip) - f0_test(xp(1:ndim,ip),ndim)))
    ENDDO
    Assert_Equal_Within(err,0,tol)

    !CALL Part1%print_info(info)
    !Assert_Equal(info,0)
    CALL Part1%destroy(info)
    Assert_Equal(info,0)
    CALL Field1%destroy(info)
    Assert_Equal(info,0)
  end test

  test neighlists
    ! test neighbour lists
    TYPE(ppm_t_particles_d),TARGET :: Part1

    CALL Part1%destroy(info)
    Assert_Equal(info,0)
    CALL Part1%initialize(np_global,info,topoid=topoid)
    Assert_Equal(info,0)

    CALL Part1%map(info,global=.TRUE.,topoid=topoid)
    Assert_Equal(info,0)
    CALL Part1%map_ghosts(info)
    Assert_Equal(info,0)

    CALL Part1%comp_neighlist(info)
    Assert_Equal(info,0)

    Assert_True(Part1%has_neighlist())
  end test

  test sop_type
    ! test procedures for sop data structures
    TYPE(ppm_t_sop_d),        TARGET  :: Part1_A
    CLASS(ppm_t_discr_data),  POINTER :: Prop1=>NULL(),Prop2=>NULL()

    CALL Part1_A%initialize(np_global,info,topoid=topoid)
    Assert_Equal(info,0)


    CALL Part1_A%map(info,global=.TRUE.,topoid=topoid)
    Assert_Equal(info,0)

    CALL Part1_A%map_ghosts(info)
    Assert_Equal(info,0)

    CALL Part1_A%comp_neighlist(info)
    Assert_Equal(info,0)

    Assert_True(Part1_A%has_neighlist())

    !creating/destroying properties of different types
    CALL Part1_A%create_prop(info,discr_data=Prop1,dtype=ppm_type_longint,name="test_li",with_ghosts=.TRUE.)
    Assert_Equal(info,0)
    CALL Part1_A%get(Prop1,wp_1li,info,with_ghosts=.TRUE.)
    Assert_Equal(info,0)


    CALL Part1_A%destroy_prop(Prop1,info)
    Assert_Equal(info,0)

    CALL Part1_A%create_prop(info,discr_data=Prop1,dtype=ppm_type_int,name="test_i",zero=.FALSE.)
    Assert_Equal(info,0)
    CALL Part1_A%get(Prop1,wp_1i,info)
    Assert_Equal(info,0)

    !Set up a velocity field and a scalar test function on the particles
    CALL Part1_A%create_prop(info,discr_data=Prop1,dtype=ppm_type_real,lda=3,name="velocity")
    Assert_Equal(info,0)
    CALL Part1_A%create_prop(info,discr_data=Prop2,dtype=ppm_type_real,lda=1,name="testf")
    Assert_Equal(info,0)
    CALL Part1_A%get(Prop1,wp_2r,info)
    Assert_Equal(info,0)
    CALL Part1_A%get(Prop2,wp_1r,info)
    Assert_Equal(info,0)
    CALL Part1_A%get_xp(xp,info)
    Assert_Equal(info,0)
    DO ip=1,Part1_A%Npart
        wp_2r(1:ndim,ip) = COS((10._MK*xp(1:ndim,ip))**2)
        wp_1r(ip) = f0_test(xp(1:ndim,ip),ndim)
    ENDDO
    wp_2r = COS(wp_2r) * Part1_A%ghostlayer

    !Move the particles with this displacement field
    CALL Part1_A%move(wp_2r,info)
    Assert_Equal(info,0)

    !Apply boundary conditions and remap the particles
    CALL Part1_A%apply_bc(info)
    Assert_Equal(info,0)
    CALL Part1_A%map(info)
    Assert_Equal(info,0)

    !Get the new ghosts
    CALL Part1_A%map_ghosts(info)
    Assert_Equal(info,0)


    !Move the particles back to their initial positions
    CALL Part1_A%get(Prop1,wp_2r,info)
    wp_2r = -wp_2r
    CALL Part1_A%move(wp_2r,info)
    Assert_Equal(info,0)

    !Re-apply boundary conditions and remap the particles
    CALL Part1_A%apply_bc(info)
    Assert_Equal(info,0)
    CALL Part1_A%map(info)
    Assert_Equal(info,0)

    CALL Part1_A%map_ghosts(info)
    Assert_Equal(info,0)

    !Compare values and check that they are still the same
    CALL Part1_A%get_xp(xp,info)
    Assert_Equal(info,0)
    CALL Part1_A%get(Prop2,wp_1r,info)
    Assert_Equal(info,0)
    err = 0.0_MK
    DO ip=1,Part1_A%Npart
        err = MAX(err,ABS(wp_1r(ip) - f0_test(xp(1:ndim,ip),ndim)))
    ENDDO
    Assert_Equal_Within(err,0,tol)

    !CALL Part1_A%print_info(info)
    !Assert_Equal(info,0)
    CALL Part1_A%destroy(info)
    Assert_Equal(info,0)
  end test

  !-------------------------------------------------------------
  ! test function
  !-------------------------------------------------------------
  PURE FUNCTION f0_test(pos,ndim)
    IMPLICIT NONE
    REAL(MK)                              :: f0_test
    INTEGER                 ,  INTENT(IN) :: ndim
    REAL(MK), DIMENSION(ndim), INTENT(IN) :: pos

    f0_test=SIN(2._MK*pi*pos(1))*COS(2._MK*pi*pos(2))*SIN(2._MK*pi*pos(ndim))
    RETURN
  END FUNCTION f0_test



end test_suite

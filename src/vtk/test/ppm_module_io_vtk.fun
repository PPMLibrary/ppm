test_suite ppm_module_io_vtk

  INTEGER, PARAMETER              :: debug = 0
  INTEGER, PARAMETER              :: MK = KIND(1.0d0) !kind(1.0e0)
  REAL(MK),PARAMETER              :: pi = 3.1415926535897931_MK
  INTEGER,PARAMETER               :: ndim=3
  INTEGER                         :: decomp,assig,tolexp
  REAL(MK)                        :: tol,min_rcp,max_rcp
  INTEGER                         :: info,comm,rank,nproc
  REAL(MK),DIMENSION(:  ),POINTER :: min_phys => NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: max_phys => NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: len_phys => NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: ghostlayer => NULL()
  INTEGER                         :: i,j,k,sum1,sum2
  INTEGER, DIMENSION(ndim*2)      :: bcdef
  REAL(MK),DIMENSION(:  ),POINTER :: cost => NULL()
  INTEGER                         :: seedsize
  INTEGER,  DIMENSION(:),ALLOCATABLE :: seed
  REAL(MK), DIMENSION(:),ALLOCATABLE :: randnb
  INTEGER                          :: isymm = 0
  LOGICAL                          :: lsymm = .false.,ok
  REAL(MK)                         :: t0,t1,t2,t3
  REAL(MK)                         :: eps

  init

    USE ppm_module_data
    USE ppm_module_topo_typedef
    USE ppm_module_init

    ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),ghostlayer(2*ndim),STAT=info)

    min_phys(1:ndim) = 0.0_MK
    max_phys(1:ndim) = 1.0_MK
    len_phys(1:ndim) = max_phys-min_phys
    ghostlayer(1:2*ndim) = max_rcp
    bcdef(1:ndim) = ppm_param_bcdef_periodic

    eps = EPSILON(1.0_MK)
    tolexp = int(log10(EPSILON(1.0_MK)))

#ifdef __MPI
    comm = MPI_COMM_WORLD
    CALL MPI_Comm_rank(comm,rank,info)
    CALL MPI_Comm_size(comm,nproc,info)
#else
    rank = 0
    nproc = 1
#endif

    CALL ppm_init(ndim,mk,tolexp,0,debug,info,99)
  end init


  finalize
    USE ppm_module_finalize

    CALL ppm_finalize(info)

    DEALLOCATE(min_phys,max_phys,len_phys)

  end finalize


  setup

    CALL RANDOM_SEED(size=seedsize)
    ALLOCATE(seed(seedsize))
    DO i=1,seedsize
      seed(i)=10+i*i*(rank+1)
    ENDDO
    CALL RANDOM_SEED(put=seed)

  end setup

  teardown

    DEALLOCATE(seed)

  end teardown

  test vtkparticles
    USE ppm_module_topo_typedef
    USE ppm_module_interfaces
    USE ppm_module_particles_typedef
    USE ppm_module_field_typedef
    USE ppm_module_MKtopo
    USE ppm_module_map
    USE ppm_module_topo_check
    USE ppm_module_util_dbg
    USE ppm_module_test

    INTEGER                         :: topoid
    INTEGER                         :: npart
    INTEGER                         :: mpart
    INTEGER                         :: newnpart
    INTEGER                         :: oldip,ip = -1
    REAL(MK),DIMENSION(:,:),POINTER :: xp => NULL()
    REAL(MK),DIMENSION(:),POINTER   :: wp => NULL()
    INTEGER ,DIMENSION(:),POINTER   :: ra => NULL()
    REAL(MK)                        :: cutoff = 0.001_MK
    REAL(MK)                        :: h
    TYPE(ppm_t_particles_d),TARGET  :: particles
    TYPE(ppm_t_field)               :: Field1
    CLASS(ppm_t_discr_data),POINTER :: Prop1 => NULL()
    INTEGER                         :: i,j
    CHARACTER(LEN=32)               :: fname

    start_subroutine("vtk")

    npart=100
    bcdef(1:ndim) = ppm_param_bcdef_freespace
    decomp = ppm_param_decomp_cuboid
    assig  = ppm_param_assign_internal

    topoid = 0
    CALL ppm_MKtopo(topoid,decomp,assig,min_phys,max_phys,bcdef,cutoff,cost,info)
    Assert_Equal(info,0)

    CALL Field1%create(ndim,info,name="F_vec") !vector field
    Assert_Equal(info,0)

    !initialize particles on a grid
    CALL particles%initialize(npart,info,topoid=topoid,distrib=ppm_param_part_init_cartesian)
    Assert_Equal(info,0)

    !define a property on the particle set by discretizing a Field on it
    CALL Field1%discretize_on(particles,info)
    Assert_Equal(info,0)

    !define a property directly
    CALL particles%create_prop(info,discr_data=Prop1,dtype=ppm_type_real,name='test_r',zero=.TRUE.)
    Assert_Equal(info,0)

    foreach p in particles(particles) with positions(x) sca_fields(F2=Prop1) vec_fields(F1=Field1)
      F1_p(1:ndim) = 10._MK * x_p(1:ndim)
      F2_p         = 42._MK
    end foreach

    fname = 'test'
    !CALL ppm_vtk_particles(fname,particles,info)
    !Assert_Equal(info,0)

    end_subroutine()
  end test

end test_suite

test_suite ppm_module_timestats

  INTEGER, PARAMETER              :: debug = 0
  INTEGER, PARAMETER              :: MK = KIND(1.0d0) !kind(1.0e0)
  INTEGER,PARAMETER               :: ndim=2
  INTEGER                         :: decomp,assig,tolexp
  REAL(MK)                        :: tol
  INTEGER                         :: info,comm,rank,nproc
  INTEGER                         :: topoid
  INTEGER                         :: np = 1000
  REAL(MK),DIMENSION(:,:),POINTER :: xp => NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: min_phys => NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: max_phys => NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: len_phys => NULL()
  INTEGER                         :: i,j,k
  INTEGER, DIMENSION(6)           :: bcdef
  REAL(MK),DIMENSION(:  ),POINTER :: cost => NULL()
  REAL(MK)                        :: latency,bandwidth
  LOGICAL                         :: ok

  init
    USE ppm_module_data
    USE ppm_module_init

    ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),STAT=info)

    min_phys(1:ndim) = 0.0_MK
    max_phys(1:ndim) = 1.0_MK
    len_phys(1:ndim) = max_phys-min_phys
    bcdef(1:6) = ppm_param_bcdef_periodic
    tol = EPSILON(1.0_MK)
    tolexp = int(log10(EPSILON(1.0_MK)))

    NULLIFY(xp)

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
    ALLOCATE(xp(ndim,np),STAT=info)
  end setup

  teardown
    DEALLOCATE(xp,STAT=info)
  end teardown

  test netstat
    ! test netstat

    USE ppm_module_data
    USE ppm_module_MKtopo
    USE ppm_module_topo_check
    USE ppm_module_test
    IMPLICIT NONE

    INTEGER :: test1,test2
    REAL(MK) :: time

    !----------------
    ! create particles
    !----------------

    CALL part_init(xp,np,min_phys,max_phys,info)


    !----------------
    ! make topology
    !----------------
    decomp = ppm_param_decomp_cuboid
    !decomp = ppm_param_decomp_xpencil
    assig  = ppm_param_assign_internal

    topoid = 0

    CALL ppm_tstats_setup(2,info)
    Assert_Equal(info,0)
    CALL ppm_tstats_add('test1',test1,info)
    Assert_Equal(info,0)
    CALL ppm_tstats_add('test2',test2,info)
    Assert_Equal(info,0)

    CALL ppm_tstats_tic(test1,1,info)
    Assert_Equal(info,0)
    CALL ppm_tstats_tic(test2,1,info)
    Assert_Equal(info,0)
    CALL ppm_MKtopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef,0.1_MK,cost,info)
    Assert_Equal(info,0)
    CALL ppm_tstats_toc(test2,1,time,info)
    Assert_Equal(info,0)
    CALL ppm_tstats_toc(test1,1,time,info)
    Assert_Equal(info,0)

    CALL ppm_tstats_collect('time.dat',info)
    Assert_Equal(info,0)

    IF (ppm_rank.EQ.0) THEN
       Open(UNIT=42, file='time.dat')
       Close(UNIT=42, STATUS='Delete')
    ENDIF

  end test
end test_suite

test_suite ppm_module_alloc

  INTEGER, PARAMETER              :: debug = 0
  INTEGER, PARAMETER              :: MK = KIND(1.0d0) !KIND(1.0e0)
  INTEGER, PARAMETER              :: ndim=2
  INTEGER                         :: tolexp
  INTEGER                         :: nproc,rank,comm
  REAL(MK),DIMENSION(:,:),POINTER :: xp => NULL()
  REAL(MK),DIMENSION(:),POINTER   :: wp => NULL()
  REAL(MK)                        :: t0
  REAL(MK)                        :: tol
  INTEGER                         :: info

  init
    USE ppm_module_init

    tol = EPSILON(1.0_MK)
    tolexp = INT(LOG10(EPSILON(1.0_MK)))

    NULLIFY(xp)

#ifdef __MPI
    comm = MPI_COMM_WORLD
    CALL MPI_Comm_rank(comm,rank,info)
    CALL MPI_Comm_size(comm,nproc,info)
#else
    rank = 0
    nproc = 1
#endif
    CALL ppm_init(ndim,mk,tolexp,comm,debug,info,99)
  end init

  finalize
    USE ppm_module_finalize

    CALL ppm_finalize(info)
  end finalize

  test alloc
    ! test simple allocation
    USE ppm_module_data
    IMPLICIT NONE
    INTEGER, DIMENSION(1) :: lda
    INTEGER               :: iopt
    INTEGER               :: info

    iopt = ppm_param_alloc_fit
    lda(1) = 1024
    CALL ppm_alloc(wp,lda,iopt,info)
    Assert_Equal(info,0)

  end test

  test dealloc
    ! test simple allocation
    USE ppm_module_data
    IMPLICIT NONE
    INTEGER, DIMENSION(1) :: lda
    INTEGER               :: iopt
    INTEGER               :: info

    iopt = ppm_param_dealloc
    lda(1) = 1024
    CALL ppm_alloc(wp,lda,iopt,info)
    Assert_Equal(info,0)
  end test
end test_suite

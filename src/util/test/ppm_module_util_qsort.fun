test_suite ppm_module_util_qsort

  INTEGER, PARAMETER :: debug = 0
  INTEGER, PARAMETER :: MK = KIND(1.0D0) !KIND(1.0E0)
  INTEGER, PARAMETER :: ndim=2
  INTEGER            :: tolexp
  INTEGER            :: info,comm

  init
    USE ppm_module_init

#ifdef __MPI
    comm = MPI_COMM_WORLD
#endif

    tolexp = INT(LOG10(EPSILON(1.0_MK)))
    CALL ppm_init(ndim,MK,tolexp,comm,debug,info,99)
  end init

  finalize
    USE ppm_module_finalize

    CALL ppm_finalize(info)
  end finalize

  setup
  end setup

  teardown
  end teardown

  test qsort
    USE ppm_module_data
    USE ppm_module_alloc
    USE ppm_module_error
    USE ppm_module_write
    USE ppm_module_util_time

    REAL(ppm_kind_double) :: t1,t2,t3,t4
    REAL                  :: temp

    INTEGER, DIMENSION(:), ALLOCATABLE :: arrayA
    INTEGER, DIMENSION(:), POINTER     :: arrayranking
    INTEGER                            :: nsize
    INTEGER                            :: i

    CHARACTER(LEN=ppm_char) :: caller='ppm_module_util_qsort_test'

    nsize=26500

    ALLOCATE(arrayA(nsize),STAT=info)
    Assert_Equal(info,0)

    DO i=1,nsize
       CALL RANDOM_NUMBER(temp)
       arrayA(i)=INT(temp*1000000)
    ENDDO

    NULLIFY(arrayranking)
    CALL ppm_util_time(t1)
    CALL ppm_util_time(t1)
    CALL ppm_util_qsort(arrayA,arrayranking,info,nsize)
    CALL ppm_util_time(t2)
    Assert_Equal(info,0)

    IF (ppm_rank.EQ.0) THEN
       stdout("time spent on qsort1=",'t2-t1')
    ENDIF

    CALL ppm_util_time(t1)
    CALL ppm_util_qsort(arrayA,info,nsize)
    CALL ppm_util_time(t2)
    Assert_Equal(info,0)

    IF (ppm_rank.EQ.0) THEN
       stdout("time spent on qsort2 in-place=",'t2-t1')
    ENDIF

    DEALLOCATE(arrayA,STAT=info)
    Assert_Equal(info,0)

    CALL ppm_alloc(arrayranking,(/0/),ppm_param_dealloc,info)
    Assert_Equal(info,0)

  end test
end test_suite

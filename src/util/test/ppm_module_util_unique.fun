test_suite ppm_module_util_unique

  INTEGER, PARAMETER :: debug = 0
  INTEGER, PARAMETER :: MK = 8
  INTEGER, PARAMETER :: ndim=2
  INTEGER            :: tolexp
  INTEGER            :: info

  init
    USE ppm_module_init

    tolexp = INT(LOG10(EPSILON(1.0_MK)))
    CALL ppm_init(ndim,MK,tolexp,0,debug,info,99)
  end init

  finalize
    USE ppm_module_finalize

    CALL ppm_finalize(info)
  end finalize

  setup
  end setup

  teardown
  end teardown

  test unique
    USE ppm_module_data
    USE ppm_module_alloc
    USE ppm_module_error
    USE ppm_module_write
    USE ppm_module_util_time

    REAL(ppm_kind_double) :: t1,t2
    REAL(MK)              :: tmp

    REAL(MK), DIMENSION(:), ALLOCATABLE :: arrayAR
    REAL(MK), DIMENSION(:), POINTER     :: UarrayAR

    INTEGER, DIMENSION(:), ALLOCATABLE :: arrayA
    INTEGER, DIMENSION(:), POINTER     :: UarrayA
    INTEGER                            :: nsize,unsize
    INTEGER                            :: i

    CHARACTER(LEN=ppm_char) :: caller='ppm_module_util_unique_test'

    CALL RANDOM_NUMBER(tmp)
    nsize=INT(tmp*100000)

    ALLOCATE(arrayA(nsize),arrayAR(nsize*2),STAT=info)
    Assert_Equal(info,0)

    DO i=1,nsize
       CALL RANDOM_NUMBER(tmp)
       arrayA(i)=INT(tmp*1000000)
    ENDDO

    !Warm-up
    CALL ppm_util_time(t1)

    NULLIFY(UarrayA)
    CALL ppm_util_time(t1)
    CALL ppm_util_unique(arrayA,UarrayA,info,unsize)
    CALL ppm_util_time(t2)
    Assert_Equal(info,0)

    IF (ppm_rank.EQ.0) THEN
       stdout("The time spent on the unique routine for an INTEGER array of size ",nsize," is =",'t2-t1')
       stdout("The size of the output unique array is =",unsize)
    ENDIF

    ! Free memory
    DEALLOCATE(arrayA,STAT=info)
    Assert_Equal(info,0)

    CALL ppm_alloc(UarrayA,(/0/),ppm_param_dealloc,info)
    Assert_Equal(info,0)

    CALL RANDOM_NUMBER(arrayAR)

    NULLIFY(UarrayAR)
    CALL ppm_util_time(t1)
    CALL ppm_util_unique(arrayAR,UarrayAR,info,unsize)
    CALL ppm_util_time(t2)
    Assert_Equal(info,0)

    IF (ppm_rank.EQ.0) THEN
       stdout("The time spent on the unique routine for a REAL array of size ",'nsize*2'," is =",'t2-t1')
       stdout("The size of the output unique array is =",unsize)
    ENDIF

    ! Free memory
    DEALLOCATE(arrayAR,STAT=info)
    Assert_Equal(info,0)

    CALL ppm_alloc(UarrayAR,(/0/),ppm_param_dealloc,info)
    Assert_Equal(info,0)

  end test
end test_suite

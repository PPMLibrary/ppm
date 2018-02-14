test_suite ppm_module_inl_hash

  USE ppm_module_substart
  USE ppm_module_substop

  INTEGER, PARAMETER :: debug = 0
  INTEGER, PARAMETER :: MK = KIND(1.0D0) !KIND(1.0E0)
  INTEGER, PARAMETER :: ndim=2
  INTEGER            :: tolexp
  INTEGER            :: info,comm
  INTEGER, PARAMETER :: size = 200

  TYPE(ppm_htable), POINTER :: ht

  init
    USE ppm_module_init

    tolexp = INT(LOG10(EPSILON(1.0_MK)))

#ifdef __MPI
    comm = MPI_COMM_WORLD
#endif

    CALL ppm_init(ndim,MK,tolexp,comm,debug,info,99)
  end init

  setup
    NULLIFY(ht)
    ALLOCATE(ht)
  end setup

  teardown
    IF (ASSOCIATED(ht)) THEN
       DEALLOCATE(ht)
       NULLIFY(ht)
    ENDIF
  end teardown

  finalize
    USE ppm_module_finalize

    CALL ppm_finalize(info)
  end finalize

  test basic
    USE ppm_module_data
    USE ppm_module_substart
    USE ppm_module_substop
    IMPLICIT NONE
    start_subroutine("basic")
    ! create
    CALL ht%create(size, info)
    Assert_True(info.EQ.0)

    ! basic usage
    Assert_True(ht%search(7_8).EQ.htable_null)
    CALL ht%insert(12_8, 144, info)
    Assert_True(info.EQ.0)
    Assert_True(ht%search(12_8).EQ.144)

    Assert_True(ht%search(7_8).EQ.htable_null)
    CALL ht%insert(7_8, 49, info)
    Assert_True(info.EQ.0)
    Assert_True(ht%search(7_8).EQ.49)

    CALL ht%remove(7_8,info)
    Assert_True(info.EQ.0)
    Assert_True(ht%search(7_8).EQ.htable_null)

    CALL ht%insert(7_8, 49, info)
    CALL ht%insert(7_8, 56, info)
    Assert_True(ht%search(7_8).EQ.56)

    ! destroy
    CALL ht%destroy(info)
    Assert_True(info.EQ.0)

    end_subroutine()
  end test

  test basic1
    USE ppm_module_data
    USE ppm_module_substart
    USE ppm_module_substop
    IMPLICIT NONE
    INTEGER :: i
    start_subroutine("basic1")
    ! create
    CALL ht%create(10000, info)
    Assert_True(info.EQ.0)

    ! basic usage
    Assert_True(ht%search(7_8).EQ.htable_null)

    DO i=1,9990,10
       CALL ht%insert(INT(i,KIND=8), i*i+1, info)
       Assert_True(info.EQ.0)
       !check the interface for both INTEGER(kind=4) and INTEGER(kind=8)
       CALL ht%insert(i+1, i*i+1+1, info)
       Assert_True(info.EQ.0)
    ENDDO

    CALL ht%insert(12000,120,info)
    Assert_True(info.EQ.0)
    CALL ht%insert(10300,103,info)
    Assert_True(info.EQ.0)
    CALL ht%insert(13000,130,info)
    Assert_True(info.EQ.0)
    CALL ht%insert(14500,145,info)
    Assert_True(info.EQ.0)
    CALL ht%insert(16000,160,info)
    Assert_True(info.EQ.0)
    CALL ht%insert(220350,220,info)
    Assert_True(info.EQ.0)

    DO i=1,9990,10
       Assert_True(ht%search(INT(i,KIND=8)).EQ.i*i+1)
       Assert_True(ht%search(i+1).EQ.i*i+1+1)
    ENDDO

    Assert_True(ht%search(2_8).EQ.3)
    Assert_True(ht%search(3_8).EQ.htable_null)
    Assert_True(ht%search(10300_8).EQ.103)
    Assert_True(ht%search(220350_8).EQ.220)
    Assert_True(ht%search(13000_8).EQ.130)

    CALL ht%grow(info)
    Assert_True(info.EQ.0)
    Assert_True(ht%nrow.EQ.32768)

    DO i=1,9990,10
       Assert_True(ht%search(INT(i,KIND=8)).EQ.i*i+1)
       Assert_True(ht%search(i+1).EQ.i*i+1+1)
    ENDDO

    Assert_True(ht%search(2_8).EQ.3)
    Assert_True(ht%search(3_8).EQ.htable_null)
    Assert_True(ht%search(10300_8).EQ.103)
    Assert_True(ht%search(220350_8).EQ.220)
    Assert_True(ht%search(13000_8).EQ.130)

    ! destroy
    CALL ht%destroy(info)
    Assert_True(info.EQ.0)


    !check the table growth when the size is not big enough
    CALL ht%create(5, info)
    Assert_True(info.EQ.0)

    CALL ht%insert(385,1, info)
    Assert_True(info.EQ.0)
    CALL ht%insert(96,2, info)
    Assert_True(info.EQ.0)
    CALL ht%insert(420,3, info)
    Assert_True(info.EQ.0)
    CALL ht%insert(122,4, info)
    Assert_True(info.EQ.0)
    CALL ht%insert(432,5, info)
    Assert_True(info.EQ.0)
    CALL ht%insert(131,6, info)
    Assert_True(info.EQ.0)
    CALL ht%insert(470,7, info)
    Assert_True(info.EQ.0)
    CALL ht%insert(157,8, info)
    Assert_True(info.EQ.0)
    CALL ht%insert(738235926,9, info)
    Assert_True(info.EQ.0)
    CALL ht%insert(1243832038_8,11, info)
    Assert_True(info.EQ.0)
    CALL ht%insert(HUGE(1)+2_8,15, info)
    Assert_True(info.EQ.0)
    CALL ht%insert(INT(HUGE(1),8)*2_8-200_8,56, info)
    Assert_True(info.EQ.0)

    Assert_True(ht%search(385_8).EQ.1)
    Assert_True(ht%search(96_8) .EQ.2)
    Assert_True(ht%search(420_8).EQ.3)
    Assert_True(ht%search(122_8).EQ.4)
    Assert_True(ht%search(738235926_8).EQ.9)
    Assert_True(ht%search(1243832038).EQ.11)
    Assert_True(ht%search(432_8).EQ.5)
    Assert_True(ht%search(131_8).EQ.6)
    Assert_True(ht%search(INT(HUGE(1),8)*2_8-200_8).EQ.56)
    Assert_True(ht%search(470_8).EQ.7)
    Assert_True(ht%search(157_8).EQ.8)
    Assert_True(ht%search(INT(HUGE(1),8)+2_8).EQ.15)

    CALL ht%grow(info)
    Assert_True(info.EQ.0)

    CALL ht%grow(info)
    Assert_True(info.EQ.0)

    CALL ht%shrink(info,3)
    Assert_True(info.EQ.0)

    ! destroy
    CALL ht%destroy(info)
    Assert_True(info.EQ.0)
    end_subroutine()
  end test

end test_suite

test_suite ppm_module_inl_hash
#ifdef __MPI
  include "mpif.h"
#endif

  TYPE(ppm_htable), POINTER :: ht

  INTEGER            :: info
  INTEGER, PARAMETER :: size = 200

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

  test basic
    USE ppm_module_data
    USE ppm_module_substart
    USE ppm_module_substop
    IMPLICIT NONE
    start_subroutine("basic")
    ! create
    CALL ht%create(size, info)
    assert_true(info.EQ.0)

    ! basic usage
    assert_true(ht%search(7_8).EQ.htable_null)
    CALL ht%insert(12_8, 144, info)
    assert_true(info.EQ.0)
    assert_true(ht%search(12_8).EQ.144)

    assert_true(ht%search(7_8).EQ.htable_null)
    CALL ht%insert(7_8, 49, info)
    assert_true(info.EQ.0)
    assert_true(ht%search(7_8).EQ.49)

    CALL ht%remove(7_8,info)
    assert_true(info.EQ.0)
    assert_true(ht%search(7_8).EQ.htable_null)

    CALL ht%insert(7_8, 49, info)
    CALL ht%insert(7_8, 56, info)
    assert_true(ht%search(7_8).EQ.56)

    ! destroy
    CALL ht%destroy(info)
    assert_true(info.EQ.0)

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
    assert_true(info.EQ.0)

    ! basic usage
    assert_true(ht%search(7_8).EQ.htable_null)

    DO i=1,9990,10
       CALL ht%insert(INT(i,KIND=8), i*i+1, info)
       assert_true(info.EQ.0)
       !check the interface for both integer(kind=4) and integer(kind=8)
       CALL ht%insert(i+1, i*i+1+1, info)
       assert_true(info.EQ.0)
    ENDDO

    CALL ht%insert(12000,120,info)
    assert_true(info.EQ.0)
    CALL ht%insert(10300,103,info)
    assert_true(info.EQ.0)
    CALL ht%insert(13000,130,info)
    assert_true(info.EQ.0)
    CALL ht%insert(14500,145,info)
    assert_true(info.EQ.0)
    CALL ht%insert(16000,160,info)
    assert_true(info.EQ.0)
    CALL ht%insert(220350,220,info)
    assert_true(info.EQ.0)

    DO i=1,9990,10
       assert_true(ht%search(INT(i,KIND=8)).EQ.i*i+1)
       assert_true(ht%search(i+1).EQ.i*i+1+1)
    ENDDO

    assert_true(ht%search(2_8).EQ.3)
    assert_true(ht%search(3_8).EQ.htable_null)
    assert_true(ht%search(10300_8).eq.103)
    assert_true(ht%search(220350_8).eq.220)
    assert_true(ht%search(13000_8).eq.130)

    CALL ht%grow(info)
    assert_true(info.EQ.0)
    assert_true(ht%nrow.eq.32768)

    DO i=1,9990,10
       assert_true(ht%search(INT(i,KIND=8)).EQ.i*i+1)
       assert_true(ht%search(i+1).EQ.i*i+1+1)
    ENDDO

    assert_true(ht%search(2_8).EQ.3)
    assert_true(ht%search(3_8).EQ.htable_null)
    assert_true(ht%search(10300_8).eq.103)
    assert_true(ht%search(220350_8).eq.220)
    assert_true(ht%search(13000_8).eq.130)

    ! destroy
    CALL ht%destroy(info)
    assert_true(info.EQ.0)


    !check the table growth when the size is not big enough
    CALL ht%create(5, info)
    assert_true(info.EQ.0)

    CALL ht%insert(385,1, info)
    assert_true(info.EQ.0)
    CALL ht%insert(96,2, info)
    assert_true(info.EQ.0)
    CALL ht%insert(420,3, info)
    assert_true(info.EQ.0)
    CALL ht%insert(122,4, info)
    assert_true(info.EQ.0)
    CALL ht%insert(432,5, info)
    assert_true(info.EQ.0)
    CALL ht%insert(131,6, info)
    assert_true(info.EQ.0)
    CALL ht%insert(470,7, info)
    assert_true(info.EQ.0)
    CALL ht%insert(157,8, info)
    assert_true(info.EQ.0)

    assert_true(ht%search(385_8).EQ.1)
    assert_true(ht%search(96_8) .EQ.2)
    assert_true(ht%search(420_8).EQ.3)
    assert_true(ht%search(122_8).EQ.4)
    assert_true(ht%search(432_8).EQ.5)
    assert_true(ht%search(131_8).EQ.6)
    assert_true(ht%search(470_8).EQ.7)
    assert_true(ht%search(157_8).EQ.8)

    CALL ht%grow(info)
    assert_true(info.EQ.0)

    CALL ht%grow(info)
    assert_true(info.EQ.0)

    CALL ht%shrink(info,3)
    assert_true(info.EQ.0)

    ! destroy
    CALL ht%destroy(info)
    assert_true(info.EQ.0)
    end_subroutine()
  end test

end test_suite

test_suite ppm_module_inl_hash
#ifdef __MPI
  include "mpif.h"
#endif

  type(ppm_htable), pointer :: ht

  integer            :: info
  integer, parameter :: size = 200

  setup
  nullify(ht)
  allocate(ht)
  end setup

  teardown
  if (associated(ht)) then
     deallocate(ht)
     nullify(ht)
  endif
  end teardown

  test basic
    USE ppm_module_data
    use ppm_module_substart
    use ppm_module_substop
    implicit none
    start_subroutine("basic")
    ! create
    call ht%create(size, info)
    assert_true(info .eq. 0)

    ! basic usage
    assert_true(ht%search(7_8) .eq. htable_null)
    call ht%insert(12_8, 144, info)
    assert_true(info .eq. 0)
    assert_true(ht%search(12_8) .eq. 144)

    assert_true(ht%search(7_8) .eq. htable_null)
    call ht%insert(7_8, 49, info)
    assert_true(info .eq. 0)
    assert_true(ht%search(7_8) .eq. 49)

    call ht%remove(7_8,info)
    assert_true(info .eq. 0)
    assert_true(ht%search(7_8) .eq. htable_null)

    call ht%insert(7_8, 49, info)
    call ht%insert(7_8, 56, info)
    assert_true(ht%search(7_8) .eq. 56)

    ! destroy
    call ht%destroy(info)
    assert_true(info .eq. 0)

    end_subroutine()
  end test

  test basic1
    USE ppm_module_data
    use ppm_module_substart
    use ppm_module_substop
    implicit none
    integer :: i
    start_subroutine("basic1")
    ! create
    call ht%create(10000, info)
    assert_true(info .eq. 0)

    ! basic usage
    assert_true(ht%search(7_8) .eq. htable_null)

    DO i=1,9990,10
       call ht%insert(INT(i,KIND=8), i*i+1, info)
       assert_true(info .eq. 0)
       !check the interface for both integer(kind=4) and integer(kind=8)
       call ht%insert(i+1, i*i+1+1, info)
       assert_true(info .eq. 0)
    ENDDO

    call ht%insert(12000,120,info)
    assert_true(info .eq. 0)
    call ht%insert(10300,103,info)
    assert_true(info .eq. 0)
    call ht%insert(13000,130,info)
    assert_true(info .eq. 0)
    call ht%insert(14500,145,info)
    assert_true(info .eq. 0)
    call ht%insert(16000,160,info)
    assert_true(info .eq. 0)
    call ht%insert(220350,220,info)
    assert_true(info .eq. 0)

    DO i=1,9990,10
       assert_true(ht%search(INT(i,KIND=8)) .eq. i*i+1)
       assert_true(ht%search(i+1) .eq. i*i+1+1)
    ENDDO

    assert_true(ht%search(2_8) .eq. 3)
    assert_true(ht%search(3_8) .eq. htable_null)
    assert_true(ht%search(10300_8).eq.103)
    assert_true(ht%search(220350_8).eq.220)
    assert_true(ht%search(13000_8).eq.130)

    ! destroy
    call ht%destroy(info)
    assert_true(info .eq. 0)

    end_subroutine()
  end test

end test_suite

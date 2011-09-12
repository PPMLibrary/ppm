test_suite ppm_module_inl_hash

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

type(ppm_htable)   :: ht
integer            :: info
integer, parameter :: size = 200

  test basic
    ! create
    call create_htable(ht, size, info)
    assert_true(info .eq. 0)

    ! basic usage
    assert_true(hash_search(ht, 7_8) .eq. htable_null)
    call hash_insert(ht, 12_8, 144, info)
    assert_true(info .eq. 0)
    assert_true(hash_search(ht, 12_8) .eq. 144)

    assert_true(hash_search(ht, 7_8) .eq. htable_null)
    call hash_insert(ht, 7_8, 49, info)
    assert_true(info .eq. 0)
    assert_true(hash_search(ht, 7_8) .eq. 49)

    ! destroy
    call destroy_htable(ht, info)
    assert_true(info .eq. 0)
  end test

  test weird_bug
    integer, dimension(9), parameter :: seed = (/11,   5234,   451,  &
                                                 173,  423,    7473, & 
                                                 9257, 817259, 819257 /)
    integer(8) :: k, inserts1, inserts2
    integer    :: i
    integer    :: o
    integer(4) :: rsize
    integer    :: hsize1 = 1000
    integer    :: hsize2 = 200
    integer    :: range  = 180 ! < hsize!!!
    integer    :: repeat = 10000
    real       :: r

    call random_seed(size=rsize)
    assert_true(rsize .le. 9)

    call create_htable(ht, hsize1, info)
    assert_true(info .eq. 0)

    ! basic usage
    assert_true(hash_search(ht, 1_8) .eq. htable_null)
    call hash_insert(ht, 1_8, 1, info)
    assert_true(info .eq. 0)
    assert_true(hash_search(ht, 1_8) .eq. 1)

    call random_seed(put=seed(1:rsize))
    inserts1 = 1
    do i=1,repeat
       ! get random number
       call random_number(r)
       k = ceiling(r*range)

       o = hash_search(ht,k)
       if (o .eq. htable_null) then
          call hash_insert(ht, k, 1, info)
          assert_true(info .eq. 0)
          inserts1 = inserts1 + 1
       end if
    end do

    ! destroy
    call destroy_htable(ht, info)
    assert_true(info .eq. 0)    

    ! create again
    call create_htable(ht, hsize2, info)
    assert_true(info .eq. 0)

    ! basic usage
    assert_true(hash_search(ht, 1_8) .eq. htable_null)
    call hash_insert(ht, 1_8, 1, info)
    assert_true(info .eq. 0)
    assert_true(hash_search(ht, 1_8) .eq. 1)

    call random_seed(put=seed(1:rsize))
    inserts2 = 1
    do i=1,repeat
       ! get random number
       call random_number(r)
       k = ceiling(r*range)

       o = hash_search(ht,k)
       if (o .eq. htable_null) then
          call hash_insert(ht, k, i, info)
          if (info .ne. 0) then
              write(*,*) "k : ", k, " | i : ", i, " | o : ", o, " | info : ", info
          end if
          assert_true(info .eq. 0)
          inserts2 = inserts2 + 1
       end if
    end do

    assert_true(inserts1 .eq. inserts2)

    ! destroy
    call destroy_htable(ht, info)
    assert_true(info .eq. 0)    
  end test

end test_suite

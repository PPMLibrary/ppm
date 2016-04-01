      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_neighlist_hash
      !-------------------------------------------------------------------------
      ! Copyright (c) 2011 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modIFy
      ! it under the terms of the GNU Lesser General Public License
      ! as published by the Free Software Foundation, either
      ! version 3 of the License, or (at your option) any later
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE create_htable(table,nelement,info)
      !!! Given number of rows of the table, creates the hash table. Number of
      !!! rows will be greater than nelement, as we use the first value that is
      !!! power of 2 and greater than nelement.
      !!!
      !!! [WARNING]
      !!! If you allocate a hashtable with more than 2^31-1 elements the hash
      !!! function will most probably produce incorrect hash keys and fail.
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_htable)      :: table
      !!! The hashtable to create.

      INTEGER, INTENT(IN   ) :: nelement
      !!! Number of desired elements.
      INTEGER, INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER(ppm_kind_int64), DIMENSION(1) :: lda
      INTEGER(ppm_kind_int64)               :: i
      INTEGER                               :: iopt

      CHARACTER(LEN=ppm_char) :: caller='create_htable'

      CALL substart(caller,t0,info)

      IF (nelement.LE.0) THEN
         table%nrow = 1_ppm_kind_int64
      ELSE
         table%nrow = 2_ppm_kind_int64**(CEILING(LOG(REAL(nelement))/LOG(2.0)))
      ENDIF

      lda(1)=table%nrow
      iopt  =ppm_param_alloc_fit
      !---------------------------------------------------------------------
      !  Allocate array for hash table keys.
      !---------------------------------------------------------------------
      CALL ppm_alloc(table%keys,lda,iopt,info)
      or_fail_alloc('htable_keys',ppm_error=ppm_error_fatal)

      !---------------------------------------------------------------------
      !  Allocate array for positions on "borders" array.
      !---------------------------------------------------------------------
      CALL ppm_alloc(table%borders_pos,lda,iopt,info)
      or_fail_alloc('htable_borders_pos',ppm_error=ppm_error_fatal)

      !---------------------------------------------------------------------
      !  Set everything to NULL.
      !---------------------------------------------------------------------
      FORALL (i=1_ppm_kind_int64:table%nrow)
         table%keys(i)        = htable_null_li
         table%borders_pos(i) = htable_null
      END FORALL

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE create_htable

      SUBROUTINE destroy_htable(table,info)

      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_htable)      :: table
      !!! The hashtable to create.

      INTEGER, INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)   :: t0
      INTEGER                 :: iopt
      INTEGER, DIMENSION(1)   :: lda

      CHARACTER(LEN=ppm_char) :: caller='destroy_htable'

      CALL substart(caller,t0,info)

      iopt=ppm_param_dealloc
      lda =0
      !---------------------------------------------------------------------
      !  Deallocate array for hash table keys.
      !---------------------------------------------------------------------
      CALL ppm_alloc(table%keys,lda,iopt,info)
      or_fail_dealloc('htable_keys',ppm_error=ppm_error_fatal)

      !---------------------------------------------------------------------
      !  Deallocate array for positions on "borders" array.
      !---------------------------------------------------------------------
      CALL ppm_alloc(table%borders_pos,lda,iopt,info)
      or_fail_dealloc('htable_borders_pos',ppm_error=ppm_error_fatal)

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE destroy_htable

      ! Note - This code makes an assumption about how your machine behaves -
      ! 1. sizeof(INTEGER(ppm_kind_int64)) == 8
      ! Limitation
      ! It will not produce the same results on little-endian and big-endian machines.
      ELEMENTAL FUNCTION h_func(table,key,seed) RESULT(hash_val)
        IMPLICIT NONE

        !---------------------------------------------------------------------
        !  Arguments
        !---------------------------------------------------------------------
        CLASS(ppm_htable),       INTENT(IN   ) :: table
        !!! The hashtable to create.

        INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
        !!! Input key (64 BITS)
        INTEGER(ppm_kind_int64), INTENT(IN   ) :: seed
        !!! Seed to be used for mixing
        INTEGER(ppm_kind_int64)                :: hash_val
        !!! Result of the hash function
        !---------------------------------------------------------------------
        !  Local variables and parameters
        !---------------------------------------------------------------------
        INTEGER(ppm_kind_int64), PARAMETER :: two32=4294967296_ppm_kind_int64
        INTEGER(ppm_kind_int64), PARAMETER :: m = 1540483477_ppm_kind_int64 !0x5bd1e995
        INTEGER(ppm_kind_int64)            :: h,k

        k=IBITS(key,0,32)

        ! This multiplication is safe.
        ! For the biggest 32 bit integer times m, it does not overflow.
        k=MOD(k*m,two32)
        k=IEOR(k,ISHFT(k,-24))
        k=MOD(k*m,two32)

        ! Initialize the hash to a 'random' value
        ! 8 is for 8 Bytes key len
        h=IEOR(seed,8_ppm_kind_int64)
        h=MOD(h*m,two32)
        h=IEOR(h,k)

        k=IBITS(key,32,32)

        k=MOD(k*m,two32)
        k=IEOR(k,ISHFT(k,-24))
        k=MOD(k*m,two32)

        h=MOD(h*m,two32)
        h=IEOR(h,k)

        ! Do a few final mixes of the hash to ensure the last few
        ! bytes are well-incorporated.
        h=IEOR(h,ISHFT(h,-13))
        h=MOD(h*m,two32)
        h=IEOR(h,ISHFT(h,-15))

        hash_val=IAND(h,table%nrow-1_ppm_kind_int64)
        RETURN
      END FUNCTION h_func

      ELEMENTAL FUNCTION h_key1(table,key,seed) RESULT(address)
        !!! Given the key value, returns corresponding address on the "borders" array.
        IMPLICIT NONE

        !---------------------------------------------------------------------
        !  Arguments
        !---------------------------------------------------------------------
        CLASS(ppm_htable),       INTENT(IN   ) :: table
        !!! The hashtable to create.

        INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
        !!! Input key
        INTEGER(ppm_kind_int64), INTENT(IN   ) :: seed
        !!! Seed to be used for mixing
        INTEGER(ppm_kind_int64)                :: address
        !!! Address that corresponds to given key
        !---------------------------------------------------------------------
        !  Local variables
        !---------------------------------------------------------------------
        address=table%h_func(key,seed)
        RETURN
      END FUNCTION h_key1

      ELEMENTAL FUNCTION h_key2(table,spot1,spot2,jump) RESULT(address)
        !!! Given the key and jump value, returns corresponding address on
        !!! "borders" array.

        IMPLICIT NONE

        !---------------------------------------------------------------------
        !  Arguments
        !---------------------------------------------------------------------
        CLASS(ppm_htable),       INTENT(IN   ) :: table
        !!! The hashtable to create.

        INTEGER(ppm_kind_int64), INTENT(IN   ) :: spot1
        !!! First spot value to avoid double computing
        INTEGER(ppm_kind_int64), INTENT(IN   ) :: spot2
        !!! Second spot value to avoid double computing
        INTEGER(ppm_kind_int64), INTENT(IN   ) :: jump
        !!! Jump value for double hashing
        INTEGER(ppm_kind_int64)                :: address
        !!! Address that corresponds to given key
        !---------------------------------------------------------------------
        !  Local variables
        !---------------------------------------------------------------------
        ! we need to work with 64 bit integers because
        ! jump*h_func might overflow, and as there is no unsigned
        ! integer the resulting value might be negative and not
        ! a valid array index
        address=MOD(spot1+spot2*jump,table%nrow)+1_ppm_kind_int64
      RETURN
      END FUNCTION h_key2

      SUBROUTINE hash_insert(table,key,value,info)
        !!! Given the key and the value, stores both in the hash table.
        !!!
        !!! [NOTE]
        !!! This routine needs to be very fast, therefor we skip the usual
        !!! chit-chat and get right to it. (-> no substart,substop unless
        !!! compiled with __DEBUG flag)
#ifdef __DEBUG
        USE ppm_module_substart
        USE ppm_module_substop
#endif
        IMPLICIT NONE

        !---------------------------------------------------------------------
        !  Arguments
        !---------------------------------------------------------------------
        CLASS(ppm_htable)                      :: table
        !!! The hashtable

        INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
        !!! Key to be stored (The 64 BITS key)
        INTEGER,                 INTENT(IN   ) :: value
        !!! Value that corresponds to given key
        INTEGER,                 INTENT(  OUT) :: info
        !!! Info for whether insertion was successful or not.
        !!! 0 if SUCCESSFUL and -1 otherwise.
        !---------------------------------------------------------------------
        !  Local variables
        !---------------------------------------------------------------------
        INTEGER(ppm_kind_int64) :: jump
        INTEGER(ppm_kind_int64) :: spot
        INTEGER(ppm_kind_int64) :: fspot
        INTEGER(ppm_kind_int64) :: sspot

#ifdef __DEBUG
        REAL(ppm_kind_double) :: t0

        CHARACTER(LEN=ppm_char) :: caller='hash_insert'

        CALL substart(caller,t0,info)
#else
        info=0
#endif
        jump=0_ppm_kind_int64

        ! Get the address corresponding to given key
        fspot=table%h_key(key,seed1)
        spot=fspot+1_ppm_kind_int64

        ! Keep on searching withing bounds of hash table.
        DO WHILE (jump.LT.table%nrow)
           ! If an empty slot found ...
           IF (table%keys(spot).EQ.htable_null_li) THEN
              ! Store the key and the corresponding value and RETURN.
              table%keys(spot) = key
              table%borders_pos(spot) = value
              RETURN
           !If the key is the same the value should be updated
           ELSE IF (table%keys(spot).EQ.key) THEN
              table%borders_pos(spot) = value
              RETURN
           ENDIF
           ! If the current slot is occupied, jump to next key that results
           ! in same hash function.
           jump=jump+1_ppm_kind_int64
           IF (jump.EQ.1_ppm_kind_int64) THEN
              ! Get the address corresponding to given key
              sspot=table%h_key(key,seed2)
           ENDIF

           spot=table%h_key(fspot,sspot,jump)
        ENDDO

        ! If NOT returned within the while-loop, that means our hash table is
        ! not sufficiently large. Hence, we should grow the table!
        CALL table%grow(info,key,value)
#ifdef __DEBUG
        or_fail("table%grow")

      9999 CONTINUE
        CALL substop(caller,t0,info)
#endif
        RETURN
      END SUBROUTINE hash_insert

      SUBROUTINE hash_insert_(table,key_,value,info)
        IMPLICIT NONE

        !---------------------------------------------------------------------
        !  Arguments
        !---------------------------------------------------------------------
        CLASS(ppm_htable)      :: table
        !!! The hashtable

        INTEGER, INTENT(IN   ) :: key_
        !!! Key to be stored
        INTEGER, INTENT(IN   ) :: value
        !!! Value that corresponds to given key
        INTEGER, INTENT(  OUT) :: info
        !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
        !!! and -1 otherwise.
        !---------------------------------------------------------------------
        !  Local variables
        !---------------------------------------------------------------------
        INTEGER(ppm_kind_int64) :: key

        key=INT(key_,KIND=ppm_kind_int64)
        CALL table%hash_insert(key,value,info)
        RETURN
      END SUBROUTINE hash_insert_

#ifdef __DEBUG
      FUNCTION hash_search(table,key) RESULT(value)
#else
      ELEMENTAL FUNCTION hash_search(table,key) RESULT(value)
#endif
      !!! Given the key, searchs the key on the hash table and returns the
      !!! corresponding value.
      !!!
      !!! [NOTE]
      !!! This routine needs to be very fast, therefor we skip the usual
      !!! chit-chat and get right to it. (-> no substart,substop unless
      !!! compiled with __DEBUG flag)
#ifdef __DEBUG
      USE ppm_module_substart
      USE ppm_module_substop
#endif
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_htable),       INTENT(IN   ) :: table
      !!! The hashtable to create.

      INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
      !!! Input key, which its corresponding value is asked for
      INTEGER                                :: value
      !!! Value corresponding to the input key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot

      LOGICAL :: KeyExist

#ifdef __DEBUG
      REAL(ppm_kind_double) :: t0

      INTEGER :: info
      ! the info and t0 are here only needed in debug mode an
      ! will not be propagated
      CALL substart('hash_search',t0,info)
#endif

      jump=0_ppm_kind_int64
      KeyExist=.TRUE.

      ! Get the other key that results in same hash key as for the input key.
      fspot=table%h_key(key,seed1)
      spot =fspot+1_ppm_kind_int64

      ! Keep on searching while we don't come across a NULL value or we don't
      ! exceed bounds of hash table.
      loop: DO WHILE(jump.LT.table%nrow)
          ! If key matches ...
          IF (table%keys(spot).EQ.key) THEN
             ! Set the return value and return
             value = table%borders_pos(spot)
             RETURN
          ELSE IF (table%keys(spot).EQ.htable_null_li) THEN
             IF (KeyExist) THEN
                KeyExist=.FALSE.
                IF (.NOT.ANY(key.EQ.table%keys)) EXIT loop
             ENDIF
          ENDIF

          ! Otherwise, keep on incrementing jump distance
          jump=jump+1_ppm_kind_int64
          IF (jump.EQ.1_ppm_kind_int64) THEN
             ! Get the address corresponding to given key
             sspot=table%h_key(key,seed2)
          ENDIF

          spot=table%h_key(fspot,sspot,jump)
      ENDDO loop

      value = htable_null

#ifdef __DEBUG
      CALL substop('hash_search',t0,info)
#endif
      RETURN
      END FUNCTION hash_search

#ifdef __DEBUG
      FUNCTION hash_search_(table,key_) RESULT(value)
#else
      ELEMENTAL FUNCTION hash_search_(table,key_) RESULT(value)
#endif
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_htable), INTENT(IN   ) :: table
      !!! The hashtable
      INTEGER,           INTENT(IN   ) :: key_
      !!! Input key, which the corresponding value is asked for
      INTEGER                          :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=INT(key_,KIND=ppm_kind_int64)
      value=table%hash_search(key)
      RETURN
      END FUNCTION hash_search_

      !TOCHECK
      SUBROUTINE hash_remove(table,key,info,existed)
      !!! Given the key, removes the elements in the hash table. Info is
      !!! set to -1 if the key was NOT found.
      !!!
      !!! [NOTE]
      !!! This routine needs to be very fast, therefor we skip the usual
      !!! chit-chat and get right to it. (-> no substart,substop unless
      !!! compiled with __DEBUG flag)
#ifdef __DEBUG
      USE ppm_module_substart
      USE ppm_module_substop
#endif
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_htable)                      :: table
      !!! The hashtable

      INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
      !!! Key to be removed
      INTEGER,                 INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      LOGICAL, OPTIONAL,       INTENT(IN   ) :: existed
      !!! User responsibility to sure whether the key exists or not!
      !!! If existed is true, the algorithm does not check its validity
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot

#ifdef __DEBUG
      REAL(ppm_kind_double) :: t0

      CALL substart('hash_remove',t0,info)
#else
      info=0
#endif
      IF (PRESENT(existed)) THEN
         IF (existed) THEN
            jump=0_ppm_kind_int64

            ! Get the address corresponding to given key
            fspot=table%h_key(key,seed1)
            spot =fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=htable_null
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  ! Get the address corresponding to given key
                  sspot=table%h_key(key,seed2)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info=ppm_error_fatal
         ENDIF
         RETURN
      ELSE
         IF (ANY(key.EQ.table%keys)) THEN
            jump=0_ppm_kind_int64

            ! Get the address corresponding to given key
            fspot=table%h_key(key,seed1)
            spot =fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=htable_null
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  ! Get the address corresponding to given key
                  sspot=table%h_key(key,seed2)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info=ppm_error_fatal
         ENDIF
      ENDIF
#ifdef __DEBUG
      9999 CONTINUE
      CALL substop('hash_remove',t0,info)
#endif
      RETURN
      END SUBROUTINE hash_remove

      SUBROUTINE hash_remove_(table, key_, info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_htable)                :: table
      !!! The hashtable

      INTEGER,           INTENT(IN   ) :: key_
      !!! Key to be removed
      INTEGER,           INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      LOGICAL, OPTIONAL, INTENT(IN   ) :: existed
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=INT(key_,KIND=ppm_kind_int64)
      IF (PRESENT(existed)) THEN
         CALL table%hash_remove(key,info,existed)
      ELSE
         CALL table%hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE hash_remove_

      SUBROUTINE grow_htable(table,info,key_,value_)
      !!! Based on the number of rows of the table, creates the hash table with
      !!! double size.
      !!!
      !!! [WARNING]
      !!! If you allocate a hashtable with more than 2^31-1 elements the hash
      !!! function will most probably produce incorrect hash keys and fail.
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_util_qsort
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_htable)                                :: table
      !!! The hashtable to grow.

      INTEGER,                           INTENT(  OUT) :: info
      INTEGER(ppm_kind_int64), OPTIONAL, INTENT(IN   ) :: key_
      !!! Key to be stored
      INTEGER,                 OPTIONAL, INTENT(IN   ) :: value_
      !!! Value that corresponds to given key
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER(ppm_kind_int64), DIMENSION(:), ALLOCATABLE :: keys_tmp
      INTEGER,                 DIMENSION(:), ALLOCATABLE :: borders_pos_tmp
      INTEGER,                 DIMENSION(:), POINTER     :: ranklist
      INTEGER                                            :: nsize
      INTEGER                                            :: i,j

      CHARACTER(LEN=ppm_char) :: caller='grow_htable'

      CALL substart(caller,t0,info)

      IF (table%nrow*2_ppm_kind_int64.GE.INT(ppm_big_i-1,ppm_kind_int64)) THEN
         !TOCHCECK
         fail("hashtable with more than 2^31-1 elements will fail",ppm_error=ppm_error_fatal)
      ENDIF

      nsize=INT(table%nrow)

      ALLOCATE(keys_tmp(nsize),SOURCE=table%keys,STAT=info)
      or_fail_alloc("keys_tmp")

      ALLOCATE(borders_pos_tmp(nsize),SOURCE=table%borders_pos,STAT=info)
      or_fail_alloc("keys_tmp & borders_pos_tmp")

      NULLIFY(ranklist)
      CALL ppm_util_qsort(keys_tmp,ranklist,info,nsize)
      or_fail("ppm_util_qsort")

      CALL table%destroy(info)
      or_fail("table%destroy")

      CALL table%create(nsize*2,info)
      or_fail("table%create")

      DO i=1,nsize
         j=ranklist(i)
         IF (keys_tmp(j).EQ.htable_null_li) CYCLE
         CALL table%insert(keys_tmp(j),borders_pos_tmp(j),info)
      ENDDO

      DEALLOCATE(keys_tmp,borders_pos_tmp,STAT=info)
      or_fail_dealloc("keys_tmp & borders_pos_tmp")

      !---------------------------------------------------------------------
      !  Deallocate ranklist array.
      !---------------------------------------------------------------------
      CALL ppm_alloc(ranklist,(/0/),ppm_param_dealloc,info)
      or_fail_dealloc('htable_keys',ppm_error=ppm_error_fatal)

      IF (PRESENT(key_)) THEN
         IF (PRESENT(value_)) THEN
            CALL table%insert(key_,value_,info)
         ENDIF
      ENDIF

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE grow_htable

      SUBROUTINE shrink_htable(table,info,shrinkage_ratio)
      !!! Based on the number of rows of the table, shrinks the hash table with
      !!! half a size or more.
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_util_qsort
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_htable)                :: table
      !!! The hashtable to shrink.

      INTEGER,           INTENT(  OUT) :: info
      INTEGER, OPTIONAL, INTENT(IN   ) :: shrinkage_ratio
      !!! OPTIONAL shrinkage_ratio (positive value).
      !!! If the size of hash table is shrinkage_ratio times bigger than the
      !!! real elements inside table, we reduce the table size
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER(ppm_kind_int64), DIMENSION(:), ALLOCATABLE :: keys_tmp
      INTEGER,                 DIMENSION(:), ALLOCATABLE :: borders_pos_tmp
      INTEGER                                            :: shrinkage_ratio_
      INTEGER,                 DIMENSION(:), POINTER     :: ranklist
      INTEGER                                            :: nsize,ssize
      INTEGER                                            :: i,j

      CHARACTER(LEN=ppm_char) :: caller='shrink_htable'

      CALL substart(caller,t0,info)

      shrinkage_ratio_=MERGE(shrinkage_ratio,4,PRESENT(shrinkage_ratio))
      shrinkage_ratio_=MERGE(4,shrinkage_ratio_,shrinkage_ratio_.LE.0)

      nsize=INT(table%nrow)
      ssize=table%size()

      IF (nsize.GE.shrinkage_ratio_*ssize) THEN
         ALLOCATE(keys_tmp(nsize),SOURCE=table%keys,STAT=info)
         or_fail_alloc("keys_tmp")

         ALLOCATE(borders_pos_tmp(nsize),SOURCE=table%borders_pos,STAT=info)
         or_fail_alloc("keys_tmp & borders_pos_tmp")

         NULLIFY(ranklist)
         CALL ppm_util_qsort(keys_tmp,ranklist,info,nsize)
         or_fail("ppm_util_qsort")

         CALL table%destroy(info)
         or_fail("table%destroy")

         CALL table%create(ssize,info)
         or_fail("table%create")

         DO i=1,nsize
            j=ranklist(i)
            IF (keys_tmp(j).EQ.htable_null_li) CYCLE
            CALL table%insert(keys_tmp(j),borders_pos_tmp(j),info)
         ENDDO

         DEALLOCATE(keys_tmp,borders_pos_tmp,STAT=info)
         or_fail_dealloc("keys_tmp & borders_pos_tmp")

         !---------------------------------------------------------------------
         !  Deallocate ranklist array.
         !---------------------------------------------------------------------
         CALL ppm_alloc(ranklist,(/0/),ppm_param_dealloc,info)
         or_fail_dealloc('htable_keys',ppm_error=ppm_error_fatal)
      ENDIF

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE shrink_htable

      FUNCTION hash_size(table) RESULT(value)
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_htable) :: table
      !!! The hashtable

      INTEGER           :: value
      !!! size of the arguments inside the hash table
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      IF (table%nrow.GT.0) THEN
         value=COUNT(table%keys.NE.htable_null_li)
      ELSE
         value=0
      ENDIF
      RETURN
      END FUNCTION hash_size



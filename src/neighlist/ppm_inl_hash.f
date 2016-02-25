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

      INTEGER, DIMENSION(1) :: lda
      INTEGER               :: iopt,i

      CHARACTER(LEN=ppm_char) :: caller='create_htable'

      CALL substart(caller,t0,info)

      IF (nelement.LE.0) THEN
         table%nrow = 1
      ELSE
         table%nrow = 2**(CEILING(LOG(REAL(nelement))/LOG(2.0)))
      ENDIF

      lda(1) = table%nrow
      iopt   = ppm_param_alloc_fit

      !---------------------------------------------------------------------
      !  Allocate array for hash table keys.
      !---------------------------------------------------------------------
      CALL ppm_alloc(table%keys, lda, iopt, info)
      or_fail_alloc('htable_keys',ppm_error=ppm_error_fatal)

      !---------------------------------------------------------------------
      !  Allocate array for positions on "borders" array.
      !---------------------------------------------------------------------
      CALL ppm_alloc(table%borders_pos, lda, iopt, info)
      or_fail_alloc('htable_borders_pos',ppm_error=ppm_error_fatal)

      !---------------------------------------------------------------------
      !  Set everything to NULL.
      !---------------------------------------------------------------------
      FORALL (i=1:lda(1))
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
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER, INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)   :: t0
      INTEGER                 :: iopt
      INTEGER, DIMENSION(1)   :: lda

      CHARACTER(LEN=ppm_char) :: caller='destroy_htable'

      CALL substart(caller,t0,info)

      iopt = ppm_param_dealloc
      lda = 0
      !---------------------------------------------------------------------
      !  Deallocate array for hash table keys.
      !---------------------------------------------------------------------
      CALL ppm_alloc(table%keys, lda, iopt, info)
      or_fail_dealloc('htable_keys',ppm_error=ppm_error_fatal)

      !---------------------------------------------------------------------
      !  Deallocate array for positions on "borders" array.
      !---------------------------------------------------------------------
      CALL ppm_alloc(table%borders_pos, lda, iopt, info)
      or_fail_dealloc('htable_borders_pos',ppm_error=ppm_error_fatal)

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE destroy_htable

      ELEMENTAL FUNCTION h_func(table, key, seed) RESULT(hash_val)

          IMPLICIT NONE

          !---------------------------------------------------------------------
          !  Arguments
          !---------------------------------------------------------------------
          CLASS(ppm_htable),       INTENT(IN   ) :: table
          !!! The hashtable to create. The pointer must not be NULL
          INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
          !!! Input key
          INTEGER(ppm_kind_int64), INTENT(IN   ) :: seed
          !!! Seed to be used for mixing
          INTEGER(ppm_kind_int64)                :: hash_val
          !!! Result of the hash function
          !!!
          !!! [NOTE]
          !!! we need to work with 64 bit integers because
          !!! jump*h_func might overflow, and as there is no unsigned
          !!! integer the resulting value might be negative and not
          !!! a valid array index

          !---------------------------------------------------------------------
          !  Local variables and parameters
          !---------------------------------------------------------------------
          INTEGER(ppm_kind_int64), PARAMETER  :: m = 1431374979_ppm_kind_int64
          INTEGER                             :: h
          INTEGER(ppm_kind_int64)             :: data
          INTEGER                             :: k
          INTEGER                             :: len

          len = 4

          h = IEOR(seed, len)
          data = key

          DO WHILE (len .GE. 4)
              k = IBITS(data, 0, 8) !data, pos, len. len = 1 always!
              k = IOR(k, ISHFT(IBITS(data,  8, 8),  8))
              k = IOR(k, ISHFT(IBITS(data, 16, 8), 16))
              k = IOR(k, ISHFT(IBITS(data, 24, 8), 24))

              k = k*m
              k = IEOR(k, ISHFT(k, -24))
              k = k*m

              h = h*m
              h = IEOR(h, k)

              data = data + 4
              len  = len - 1
          ENDDO

          SELECT CASE (len)
          CASE (3)
             h = IEOR(h, ISHFT(IBITS(data, 16, 8), 16))
             h = IEOR(h, ISHFT(IBITS(data,  8, 8),  8))
             h = IEOR(h, IBITS(data, 0, 8))
             h = h*m

          CASE (2)
             h = IEOR(h, ISHFT(IBITS(data,  8, 8),  8))
             h = IEOR(h, IBITS(data, 0, 8))
             h = h*m

          CASE (1)
             h = IEOR(h, IBITS(data, 0, 8))
             h = h*m

          END SELECT

          h = IEOR(h, ISHFT(h, -13))
          h = h*m
          h = IEOR(h, ISHFT(h, -15))

          hash_val = IAND(h, table%nrow - 1)
          RETURN
      END FUNCTION

      ELEMENTAL FUNCTION h_key(table, key, jump) RESULT(address)
      !!! Given the key and jump value, returns corresponding address on
      !!! "borders" array.

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_htable),       INTENT(IN   ) :: table
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
      !!! Input key, which corresponds to address requested
      INTEGER,                 INTENT(IN   ) :: jump
      !!! Jump value for double hashing
      INTEGER                                :: address
      !!! Address that corresponds to given key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: int_addr
      ! we need to work with 64 bit integers because
      ! jump*h_func might overflow, and as there is no unsigned
      ! integer the resulting value might be negative and not
      ! a valid array index

      int_addr = 1_ppm_kind_int64 + &
      & MOD((table%h_func(key,seed1)+jump*table%h_func(key,seed2)),table%nrow)
      address  = INT(int_addr)
      RETURN
      END FUNCTION h_key

      SUBROUTINE hash_insert(table, key, value, info)
      !!! Given the key and the value, stores both in the hash table. Info is
      !!! set to -1 if size of the hash table is not sufficient.
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
      !!! Key to be stored
      INTEGER,                 INTENT(IN   ) :: value
      !!! Value that corresponds to given key
      INTEGER,                 INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: jump
      INTEGER :: spot

#ifdef __DEBUG
      REAL(ppm_kind_double) :: t0

      CHARACTER(LEN=ppm_char) :: caller='hash_insert'

      CALL substart(caller,t0,info)
#else
      info = 0
#endif
      jump = 0
      ! Keep on searching withing bounds of hash table.
      DO WHILE (jump .LT. table%nrow)
         ! Get the address corresponding to given key
         spot = table%h_key(key, jump)
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
         jump = jump + 1
      ENDDO

      ! If NOT returned within the while-loop, that means our hash table is
      ! not sufficiently large. Hence, regrowth will be needed!
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
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
      !!! Input key, which the corresponding value is asked for
      INTEGER                                :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: jump
      INTEGER :: spot

#ifdef __DEBUG
      REAL(ppm_kind_double) :: t0

      INTEGER :: info
      ! the info and t0 are here only needed in debug mode an
      ! will not be propagated
      CALL substart('hash_search',t0,info)
#endif

      value = htable_null
      jump = 0

      ! Keep on searching while we don't come across a NULL value or we don't
      ! exceed bounds of hash table.
      loop: DO WHILE(jump .LT. table%nrow)
          ! Get the other key that results in same hash key as for the input
          ! key.
          spot = table%h_key(key, jump)
          ! If key matches ...
          IF (table%keys(spot).EQ.key) THEN
             ! Set the return value and return
             value = table%borders_pos(spot)
             RETURN
          ELSE IF (table%keys(spot).EQ.htable_null_li) THEN
             IF (.NOT.ANY(key.EQ.table%keys)) EXIT loop
          ENDIF
          ! Otherwise, keep on incrementing jump distance
          jump = jump + 1
      ENDDO loop

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

      !TODO check this sub
      !Yaser
      !This function PROBABLY suffers from a bug!
      SUBROUTINE hash_remove(table, key, info,existed)
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
      INTEGER :: jump
      INTEGER :: spot !,spot0

#ifdef __DEBUG
      REAL(ppm_kind_double) :: t0

      CALL substart('hash_remove',t0,info)
#else
      info = 0
#endif
      IF (PRESENT(existed)) THEN
         IF (existed) THEN
            jump = 0
            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! Get the address corresponding to given key
               spot = table%h_key(key, jump)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=htable_null
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump = jump + 1
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
         RETURN
      ELSE
         IF (ANY(key.EQ.table%keys)) THEN
            jump = 0
            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! Get the address corresponding to given key
               spot = table%h_key(key, jump)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=htable_null
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump = jump + 1
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
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
      INTEGER                                            :: nsize,i

      CHARACTER(LEN=ppm_char) :: caller='grow_htable'

      CALL substart(caller,t0,info)

      nsize=table%nrow

      ALLOCATE(keys_tmp(nsize),SOURCE=table%keys,STAT=info)
      or_fail_alloc("keys_tmp")

      ALLOCATE(borders_pos_tmp(nsize),SOURCE=table%borders_pos,STAT=info)
      or_fail_alloc("keys_tmp & borders_pos_tmp")

      nsize=table%nrow*2

      IF (nsize.GE.ppm_big_i-1) THEN
         !TOCHCECK
         fail("hashtable with more than 2^31-1 elements will fail",ppm_error=ppm_error_fatal)
      ENDIF

      CALL table%destroy(info)
      or_fail("table%destroy")

      CALL table%create(nsize,info)
      or_fail("table%create")

      DO i=1,SIZE(keys_tmp)
         IF (keys_tmp(i).EQ.htable_null_li) CYCLE
         CALL table%insert(keys_tmp(i),borders_pos_tmp(i),info)
      ENDDO

      DEALLOCATE(keys_tmp,borders_pos_tmp,STAT=info)
      or_fail_dealloc("keys_tmp & borders_pos_tmp")

      IF (PRESENT(key_)) THEN
         IF (PRESENT(value_)) THEN
            CALL table%insert(key_,value_,info)
         ENDIF
      ENDIF

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE grow_htable

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


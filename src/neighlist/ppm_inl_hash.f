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
      TYPE(ppm_htable), INTENT(INOUT) :: table
      !!! The hashtable to create.
      INTEGER, INTENT(IN   )   :: nelement
      !!! Number of desired elements.
      INTEGER, INTENT(  OUT)  :: info

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(1) :: lda
      INTEGER               :: iopt
      REAL(ppm_kind_double) :: t0

      CALL substart('create_htable',t0,info)

      table%nrow = 2**(CEILING(LOG(REAL(nelement))/LOG(2.0)))
      
      lda(1) = table%nrow
      iopt   = ppm_param_alloc_fit

      !---------------------------------------------------------------------
      !  Allocate array for hash table keys.
      !---------------------------------------------------------------------
      CALL ppm_alloc(table%keys, lda, iopt, info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'create_htable',     &
 &                      'htable_keys',__LINE__,info)
          GOTO 9999
      END IF

      !---------------------------------------------------------------------
      !  Allocate array for positions on "borders" array.
      !---------------------------------------------------------------------
      CALL ppm_alloc(table%borders_pos, lda, iopt, info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'create_htable',     &
 &                      'htable_borders_pos',__LINE__,info)
          GOTO 9999
      END IF

      !---------------------------------------------------------------------
      !  Set everything to NULL.
      !---------------------------------------------------------------------
      table%keys        = htable_null
      table%borders_pos = htable_null
      
 9999 CONTINUE
      CALL substop('create_htable',t0,info)
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
      TYPE(ppm_htable), INTENT(INOUT) :: table
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER, INTENT(OUT)    :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)   :: t0
      INTEGER                 :: iopt
      INTEGER, DIMENSION(1)   :: lda

      CALL substart('destroy_htable',t0,info)

      lda = 0
      iopt = ppm_param_dealloc

      !---------------------------------------------------------------------
      !  Deallocate array for hash table keys.
      !---------------------------------------------------------------------
      CALL ppm_alloc(table%keys, lda, iopt, info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_dealloc,'destroy_htable',     &
 &                      'htable_keys',__LINE__,info)
          GOTO 9999
      END IF

      !---------------------------------------------------------------------
      !  Deallocate array for positions on "borders" array.
      !---------------------------------------------------------------------
      CALL ppm_alloc(table%borders_pos, lda, iopt, info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_dealloc,'destroy_htable',     &
 &                      'htable_borders_pos',__LINE__,info)
          GOTO 9999
      END IF
 
 9999 CONTINUE
      CALL substop('destroy_htable',t0,info)
      END SUBROUTINE destroy_htable

      ELEMENTAL FUNCTION h_func(table, key, seed) result(hash_val)
          IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
          TYPE(ppm_htable), INTENT(IN) :: table
          !!! The hashtable to create. The pointer must not be NULL
          INTEGER(ppm_kind_int64), INTENT(in) :: key
          !!! Input key
          INTEGER(ppm_kind_int64), INTENT(in) :: seed
          !!! Seed to be used for mixing
          INTEGER(ppm_kind_int64)             :: hash_val
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
          INTEGER(ppm_kind_int64), PARAMETER  :: m = 1431374979
          INTEGER,                 PARAMETER  :: r = 24
          INTEGER                             :: h
          INTEGER(ppm_kind_int64)             :: data
          INTEGER                             :: k
          INTEGER                             :: len

          len = 4

          h = IEOR(seed, len)
          data = key

          DO WHILE(len .GE. 4)
              k = IBITS(data, 0, 8) !data, pos, len. len = 1 always!
              k = IOR(k, ISHFT(IBITS(data, 8, 8),  8))
              k = IOR(k, ISHFT(IBITS(data, 16, 8), 16))
              k = IOR(k, ISHFT(IBITS(data, 24, 8), 24))

              k = k*m
              k = IEOR(k, ISHFT(k, -1*r))
              k = k*m

              h = h*m
              h = IEOR(h, k)

              data = data + 4
              len  = len - 1
          END DO

          IF    (len .EQ. 3)    THEN
              h = IEOR(h, ISHFT(IBITS(data, 16, 8), 16))
              len = len - 1
          END IF

          IF(len .EQ. 2)    THEN
              h = IEOR(h, ISHFT(IBITS(data,  8, 8),  8))
              len = len - 1
          END IF

          IF(len .EQ. 1)    THEN
              h = IEOR(h, IBITS(data, 0, 8))
              h = h*m
          END IF

          h = IEOR(h, ISHFT(h, -13))
          h = h*m
          h = IEOR(h, ISHFT(h, -15))
          h = IAND(h, table%nrow - 1)

          hash_val = h
      END FUNCTION

      ELEMENTAL FUNCTION h_key(table, key, jump) result(address)
      !!! Given the key and jump value, returns corresponding address on
      !!! "borders" array.
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      TYPE(ppm_htable), INTENT(IN)        :: table
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER(ppm_kind_int64), INTENT(in) :: key
      !!! Input key, which corresponds to address requested
      INTEGER,                 INTENT(in) :: jump
      !!! Jump value for double hashing
      INTEGER(ppm_kind_int64)             :: int_addr
      ! we need to work with 64 bit integers because
      ! jump*h_func might overflow, and as there is no unsigned
      ! integer the resulting value might be negative and not
      ! a valid array index
      INTEGER                             :: address
      !!! Address that corresponds to given key
      int_addr = 1 + MOD((h_func(table, key, seed1) + jump*h_func(table, key,  &
 &                  seed2)), table%nrow)
      address = int_addr
      END FUNCTION h_key

      SUBROUTINE hash_insert(table, key, value, info)
      !!! Given the key and the value, stores both in the hash table. Info is
      !!! set to -1 if size of the hash table is not sufficient.
      !!!
      !!! [NOTE]
      !!! This routine needs to be very fast, therefor we skip the usual
      !!! chit-chat and get right to it. (-> no substart,substop)
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      TYPE(ppm_htable), INTENT(INOUT)      :: table
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER(ppm_kind_int64), INTENT(in)  :: key
      !!! Key to be stored
      INTEGER,                 INTENT(in)  :: value
      !!! Value that corresponds to given key
      INTEGER,                 INTENT(out) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER                              :: jump
      INTEGER                              :: spot

      info = 0
      jump = 0
      ! Keep on searching withing bounds of hash table.
      DO WHILE(jump .LT. table%nrow)
          ! Get the address corresponding to given key
          spot = h_key(table, key, jump)
          ! If an empty slot found ...
          IF(table%keys(spot) .EQ. htable_null)  THEN
              ! Store the key and the corresponding value and RETURN.
              table%keys(spot) = key
              table%borders_pos(spot) = value
              GOTO 9999
          ! If the current slot is occupied, jump to next key that results
          ! in same hash function.
          ELSE
              jump = jump + 1
          END IF
      END DO

      ! If NOT returned within the while-loop, that means our hash table is
      ! not sufficiently large. Hence, regrowth will be needed!
      info = -1

 9999 CONTINUE

      END SUBROUTINE hash_insert

      ELEMENTAL FUNCTION hash_search(table,key) result(value)
      !!! Given the key, searchs the key on the hash table and returns the
      !!! corresponding value.
      !!!
      !!! [NOTE]
      !!! This routine needs to be very fast, therefor we skip the usual
      !!! chit-chat and get right to it. (-> no substart,substop)
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      TYPE(ppm_htable), INTENT(IN)          :: table
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER(ppm_kind_int64), INTENT(in)   :: key
      !!! Input key, which the corresponding value is asked for
      INTEGER                               :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER                               :: jump
      INTEGER                               :: spot

      value = htable_null
      jump = 0

      spot = h_key(table, key, jump)
      ! Keep on searching while we don't come across a NULL value or we don't
      ! exceed bounds of hash table.
      DO WHILE((jump .LT. table%nrow) .AND. (table%keys(spot)  &
               &.NE. htable_null))
          ! If key matches ...
          IF(table%keys(spot) .EQ. key)  THEN
              ! Set the return value and return
              value = table%borders_pos(spot)
              RETURN
          ! Otherwise, keep on incrementing jump distance
          ELSE
              jump = jump + 1
          END IF
          ! Get the other key that results in same hash key as for the input
          ! key.
          spot = h_key(table, key, jump)
      END DO

      END FUNCTION hash_search

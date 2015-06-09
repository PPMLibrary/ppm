      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_map_part_store
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
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

      SUBROUTINE ppm_map_part_store(info)
      !!! This routine stores the state of the current/active particle mapping.
      !!!
      !!! This is especially useful for alternating ghost-get/ghost-put
      !!! mappings, which are used for implementing the communication in
      !!! symmetric particle interactions. The `ppm_map_part_ghost_get` mapping
      !!! can then be stored after its first call and reloaded after each
      !!! ghost-put mapping.

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_state
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  This routine is ALL integer, except the timing variable t0
      !  We pick double as the precision - no real reason for that
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(  OUT)  :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK) :: t0

      INTEGER, DIMENSION(1) :: ldu
      INTEGER               :: i,j,k
      INTEGER               :: iopt

      CHARACTER(LEN=ppm_char) :: caller='ppm_map_part_store'
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         IF (.NOT. ppm_initialized) THEN
            fail('Please call ppm_init first!')
         ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  Store the map type (is used in _send() line 229)
      !-------------------------------------------------------------------------
      ppm_map_type_state = ppm_map_type

      !-------------------------------------------------------------------------
      !  store the ppm_nsendlist: the number of processers to send to
      !-------------------------------------------------------------------------
      ppm_nsendlist_state = ppm_nsendlist

      !-------------------------------------------------------------------------
      !  store the ppm_nsendbuffer: the size of the send buffer - here zero
      !-------------------------------------------------------------------------
      ppm_nsendbuffer_state = 0

      !-------------------------------------------------------------------------
      !  store the ppm_buffer_set: the current number of sets stored - here zero
      !-------------------------------------------------------------------------
      ppm_buffer_set_state  = 0

      !-------------------------------------------------------------------------
      !  store the ppm_psendbuffer:
      !  first allocate memory for it
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist + 1
      CALL ppm_alloc(ppm_psendbuffer_state,ldu,iopt,info)
      or_fail_alloc('allocation of psendbuffer_state failed', &
      & ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  then store it: pointer to the first particle in the send buffer
      !-------------------------------------------------------------------------
      DO k=1,ppm_nsendlist+1
         ppm_psendbuffer_state(k) = ppm_psendbuffer(k)
      ENDDO

      !-------------------------------------------------------------------------
      !  store the ppm_buffer2part
      !  first allocate memory for it
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = SIZE(ppm_buffer2part)
      CALL ppm_alloc(ppm_buffer2part_state,ldu,iopt,info)
      or_fail_alloc('allocation of buffer2part_state failed', &
      & ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  then store it: pointer to the particle id of the j-th entry in the
      !  send buffer
      !-------------------------------------------------------------------------
      DO k=1,ppm_nsendlist
         DO j=ppm_psendbuffer(k),ppm_psendbuffer(k+1)-1
            ppm_buffer2part_state(j) = ppm_buffer2part(j)
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Store the ppm_nrecvlist and ppm_sendlist
      !-------------------------------------------------------------------------
      ppm_nrecvlist_state = ppm_nrecvlist
      ppm_nsendlist_state = ppm_nsendlist

      !-------------------------------------------------------------------------
      !  Store the receive list
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nrecvlist
      CALL ppm_alloc(ppm_irecvlist_state,ldu,iopt,info)
      or_fail_alloc('allocation of irecvlist_state failed', &
      & ppm_error=ppm_error_fatal)

      DO k=1,ppm_nrecvlist
         ppm_irecvlist_state(k) = ppm_irecvlist(k)
      ENDDO

      !-------------------------------------------------------------------------
      !  Store the send list
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(ppm_isendlist_state,ldu,iopt,info)
      or_fail_alloc('allocation of isendlist_state failed', &
      & ppm_error=ppm_error_fatal)

      DO k=1,ppm_nsendlist
         ppm_isendlist_state(k) = ppm_isendlist(k)
      ENDDO

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE ppm_map_part_store

      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_map_part_ghost_put
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

      SUBROUTINE ppm_map_part_ghost_put(topoid,info)
      !!! This routine puts back ghost particle values/properties to the
      !!! corresponding real particles. This is very useful in the case of
      !!! symmetric interactions as there the ghost particles are also updated.
      !!!
      !!! [IMPORTANT]
      !!! This routine can only be called after ghost particles have been
      !!! created using the `ppm_map_part_ghost_get` (+push/send/pop sequence)
      !!! routine.
      !!!
      !!! [WARNING]
      !!! It is an error to call `ppm_map_part_pop` in the
      !!! push-send-pop sequence of a `map_part_ghost_put` call. Instead, you
      !!! *must* call `ppm_map_part_ghost_pop` for popping the
      !!! properties pushed onto the put buffers.
      !!!
      !!! [TIP]
      !!! If you need to do alternating ghost-get and ghost-put sequences
      !!! you may want to use `ppm_map_part_store` and `ppm_map_part_load` to
      !!! store and load the internal buffers for ghost_get and spare yourself
      !!! the (costly) call to `ppm_map_part_ghost_get`. Positions can be pushed
      !!! using the `ppm_map_part_ghost_push` routine.
      !!!
      !!! [NOTE]
      !!! .Implementation Notes
      !!! ======================================================================
      !!! This routine swaps the send and receive lists of the
      !!! mapping routine to allow sending back the value/data
      !!! computed on the ghost particles. The routine must be
      !!! called after `ppm_map_part_send` and should be followed
      !!! by a `ppm_map_part_push`, `ppm_map_part_send` and `ppm_map_part_pop`
      !!! to get.
      !!!
      !!! This implementation only allows *one* ghost put at a
      !!! time; To allow multiple ghosts we have to pass
      !!! the ppm_user_hack to the ghost_pop OR to copy the
      !!! userpassed array to the internally stored ppm_ghosthack
      !!! used (and *not* passed) by the ghost_pop()
      !!! ======================================================================
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_topo_typedef
      USE ppm_module_alloc
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_check_id
      USE ppm_module_util_commopt
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  This routine is ALL integer, except the timing variable t0
      !  We pick double as the precision - no real reason for that
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER      :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN   )  :: topoid
      !!! Topology ID
      INTEGER, INTENT(  OUT)  :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(2)     :: ldu
      TYPE(ppm_t_topo), POINTER :: topo
      INTEGER                   :: i,j,k
      INTEGER                   :: iopt
      CHARACTER(ppm_char)       :: mesg
      LOGICAL                   :: valid
      REAL(MK)                  :: t0
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part_ghost_put',t0,info)


      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      topo => ppm_topo(topoid)%t

      ! if there is still some data left in the buffer, warn the user
      IF (ppm_buffer_set .GT. 0) THEN
        info = ppm_error_warning
        CALL ppm_error(ppm_err_map_incomp,'ppm_map_part_ghost_put',  &
     &       'Buffer was not empty. Possible loss of data!',__LINE__,info)
      ENDIF

      !----------------------------------------------------------------------
      !  first check if the optimal communication protocol is known
      !----------------------------------------------------------------------
      IF (.NOT.topo%isoptimized) THEN
        !-------------------------------------------------------------------
        !  if not: determine it before calling map_part_ghost_put
        !-------------------------------------------------------------------
        CALL ppm_util_commopt(topoid,info)
        IF (ppm_debug .GT. 1) THEN
            DO i=1,topo%nneighproc
                WRITE(mesg,'(A,I4)') 'have neighbor: ', topo%ineighproc(i)
                CALL ppm_write(ppm_rank,'ppm_map_part_ghost_put',mesg,info)
            ENDDO
            DO i=1,topo%ncommseq
                WRITE(mesg,'(A,I4)') 'communicate: ', topo%icommseq(i)
                CALL ppm_write(ppm_rank,'ppm_map_part_ghost_put',mesg,info)
            ENDDO
        ENDIF
      ENDIF


      !-------------------------------------------------------------------------
      !  Save the map type for the subsequent calls
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_param_map_ghost_put

      !-------------------------------------------------------------------------
      !  Find the required size of the ppm_ghosthack array: the total number of
      !  particles send times the ppm_ncommseq
      !-------------------------------------------------------------------------
      j = 0
      DO k=1,topo%ncommseq
         j = MAX(j,(ppm_psendbuffer(k+1) - ppm_psendbuffer(k)))
      ENDDO

      !-------------------------------------------------------------------------
      !  The 2nd leading dimension is plus one to catch all ghosts -)
      !-------------------------------------------------------------------------
      ldu(1)  = j + 2
      ldu(2)  = topo%ncommseq + 1

      !-------------------------------------------------------------------------
      !  Allocate the array
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_grow
      CALL ppm_alloc(ppm_ghosthack,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_put', &
     &        'allocation of ppm_ghosthack failed',__LINE__,info)
          GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      !  now store the data in ppm_ghosthack()
      !-------------------------------------------------------------------------
      i = 0
      DO k=1,ppm_nsendlist
         !----------------------------------------------------------------------
         !  Store the number of particles send to the k-th proc
         !----------------------------------------------------------------------
         ppm_ghosthack(1,k) = ppm_psendbuffer(k+1) - ppm_psendbuffer(k)
         DO j=1,ppm_ghosthack(1,k)
            !-------------------------------------------------------------------
            !  Increment the counter
            !-------------------------------------------------------------------
            i = i + 1

            !-------------------------------------------------------------------
            !  store the id of the particles send to the k-th proc
            !-------------------------------------------------------------------
            ppm_ghosthack(j+2,k) = ppm_buffer2part(i)
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  now store the data in ppm_ghosthack()
      !-------------------------------------------------------------------------
      DO k=1,ppm_nsendlist + 1
         !----------------------------------------------------------------------
         !  ppm_ghosthack(2,k) stores the pointer to the first particle received
         !  from the k-th proc
         !----------------------------------------------------------------------
         ppm_ghosthack(2,k) = ppm_precvbuffer(k)
      ENDDO

      !-------------------------------------------------------------------------
      !  Prepare what to send:
      !  we have to send the ghost back :-) so consequtive piece of data
      !  starting at Npart + 1 and ending at Mpart; the pointer to the first
      !  particle was stored in ppm_precvbuffer(k) and now in ppm_ghosthack(2,k).
      !  The buffer2part is easy as the data is consequtive
      !-------------------------------------------------------------------------
      i = -1
      DO k=1,ppm_nsendlist
         i = MAX(i,ppm_ghosthack(2,k+1))
      ENDDO
      ldu(1) = i
      CALL ppm_alloc(ppm_buffer2part,ldu,iopt,info)
      DO k=1,ppm_nsendlist
         !----------------------------------------------------------------------
         !  ppm_psendbuffer(k) stores the pointer to the first particle received
         !  from the k-th proc
         !----------------------------------------------------------------------
         ppm_psendbuffer(k) = ppm_ghosthack(2,k)
         DO j=ppm_ghosthack(2,k),ppm_ghosthack(2,k+1)-1
            ppm_buffer2part(j) = j
         ENDDO
      ENDDO
      ppm_psendbuffer(ppm_nsendlist+1) = ppm_ghosthack(2,ppm_nsendlist+1)

      !-------------------------------------------------------------------------
      !  Prepare what to receive:
      !  We will recieve the function value carried by the ghosts - these ghosts
      !  are REAL particles on our processor so we should NOT pop them as usual
      !  Notice the push does not care (after all it is just packing things for
      !  sending), neither does the send (it is just sending the buffer of and
      !  receiving one) - the pop cares - so much that we need a special pop
      !  Now in the full version, the ppm_ghosthack should be passed back out
      !  to the user (OR be given a 3rd index with an internal 'pset' entry).
      !  so either we pass the ppm_user_hack to the ghost_pop OR we copy the
      !  userpassed array to the internally stored ppm_ghosthack used
      !  (and NOT passed) by the ghost_pop() ... let us see !
      !     DO i
      !        ppm_ghosthack() = userpassed()
      !     ENDDO
      !  For the time being we assume ONE set of ghost and store this array in
      !  ppm_module_data. This piece of data will be used directly by the
      !  ghost_pop routine
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  flip the send and receive lists
      !-------------------------------------------------------------------------
      DO k=1,ppm_nrecvlist
         j                = ppm_irecvlist(k)
         ppm_irecvlist(k) = ppm_isendlist(k)
         ppm_isendlist(k) = j
      ENDDO

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_part_ghost_put',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_map_part_ghost_put',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (topoid .EQ. ppm_param_topo_undefined) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_no_topo,'ppm_map_part_ghost_put',  &
     &            'This routine needs a topology defined topo',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (topoid .NE. ppm_param_topo_undefined) THEN
            CALL ppm_check_topoid(topoid,valid,info)
            IF (.NOT. valid) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_map_part_ghost_put',  &
     &               'topoid out of range',__LINE__,info)
                GOTO 8888
            ENDIF
        ENDIF
 8888   CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_map_part_ghost_put

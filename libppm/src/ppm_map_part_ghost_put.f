      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_map_part_ghost_put
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine swaps the send and receive lists of the 
      !                 mapping routine to allow sending back the value/data 
      !                 compute on the ghost particles. The routine must be
      !                 called after ppm_map_part_send() and should be followed
      !                 by a ppm_map_part_push/send/pop to get 
      !
      !  Input        : 
      !                                    
      !  Output       : info         (I) : return status, 0 on success
      !
      !  Remarks      : This implementation only allows ONE ghost put at the
      !                 time; To allow multiple ghosts we either have to pass 
      !                 the ppm_user_hack to the ghost_pop OR to copy the 
      !                 userpassed array to the internally stored ppm_ghosthack 
      !                 used (and NOT passed) by the ghost_pop()
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_part_ghost_put.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2006/09/04 18:34:50  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.5  2006/04/06 14:52:35  walther
      !  cleanup of the arguments to the routine and of the comments.
      !
      !  Revision 1.4  2006/02/03 09:41:26  ivos
      !  Added the PRELIMINARY ghost_put functionality. Still needs clean-up,
      !  but should work.
      !
      !  Revision 1.3  2004/10/01 16:09:07  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.2  2004/07/26 07:42:45  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.1  2004/02/19 15:56:36  walther
      !  Initial implementation (not complete!).
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_map_part_ghost_put(info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  This routine is ALL integer, except the timing variable t0
      !  We pick double as the precision - no real reason for that 
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER, INTENT(OUT)  :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(2) :: ldu
      INTEGER               :: i,j,k,topoid
      INTEGER               :: iopt
      CHARACTER(ppm_char)   :: mesg
      REAL(MK)              :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part_ghost_put',t0,info)
      topoid = ppm_topoid

      !-------------------------------------------------------------------------
      !  Save the map type for the subsequent calls 
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_param_map_ghost_put

      !-------------------------------------------------------------------------
      !  Find the required size of the ppm_ghosthack array: the total number of 
      !  particles send times the ppm_ncommseq
      !-------------------------------------------------------------------------
      j = 0
      DO k=1,ppm_ncommseq(topoid)
         j = MAX(j,(ppm_psendbuffer(k+1) - ppm_psendbuffer(k)))
      ENDDO

      !-------------------------------------------------------------------------
      !  The 2nd leading dimension is plus one to catch all ghosts -)
      !-------------------------------------------------------------------------
      ldu(1)  = j + 2
      ldu(2)  = ppm_ncommseq(topoid) + 1

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
      END SUBROUTINE ppm_map_part_ghost_put

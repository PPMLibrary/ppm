      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_map_part_store
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine stores the state of the current/active
      !                 particle mapping.
      !
      !  Input        : 
      !                                    
      !  Output       : info         (I) : return status, 0 on success
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_part_store.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2006/10/10 21:11:49  walther
      !  *** empty log message ***
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_map_part_store(info)

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
      INTEGER, PARAMETER    :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER, INTENT(OUT)  :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(1) :: ldu
      INTEGER               :: i,j,k
      INTEGER               :: iopt
      REAL(MK)              :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part_store',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_load',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  Store the map type (is used in _send() line 229)
      !-------------------------------------------------------------------------
      ppm_map_type_state = ppm_map_type

      !-------------------------------------------------------------------------
      !  store the ppm_nsendlist: the number of processers to send to
      !-------------------------------------------------------------------------
      ppm_nsendlist_state   = ppm_nsendlist

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
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_store', &
     &        'allocation of psendbuffer_state failed',__LINE__,info)
          GOTO 9999
      ENDIF
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
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_store', &
     &        'allocation of buffer2part_state failed',__LINE__,info)
          GOTO 9999
      ENDIF
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
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_store', &
     &        'allocation of irecvlist_state failed',__LINE__,info)
          GOTO 9999
      ENDIF
      DO k=1,ppm_nrecvlist
         ppm_irecvlist_state(k) = ppm_irecvlist(k)
      ENDDO

      !-------------------------------------------------------------------------
      !  Store the send list
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(ppm_isendlist_state,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_state', &
     &        'allocation of isendlist_state failed',__LINE__,info)
          GOTO 9999
      ENDIF
      DO k=1,ppm_nsendlist
         ppm_isendlist_state(k) = ppm_isendlist(k)
      ENDDO

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_part_store',t0,info)
      RETURN
      END SUBROUTINE ppm_map_part_store

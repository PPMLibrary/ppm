      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_map_part_ghost
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine is the main interface routine to the 
      !                 mapping routines for the ghost particles.
      !
      !  Input        : pdata(:[,:]) (O) : the particle data (positions for
      !                                    ghost_get and ghost_put). Can be
      !                                    either 1d or 2d array of type
      !                                    single,double,integer,logical,single
      !                                    complex or double complex.
      !                                    Positions need to be a 2d array
      !                                    of type single or double.
      !                 [lda]        (I) : leading dimension of pdata (only
      !                                    important for push and pop) in
      !                                    case pdata is a 2d array. For 1d
      !                                    arrays, omit this argument.
      !                 Npart        (I) : the number of particles 
      !                 isymm        (I) : indicator for the use of symmetry 
      !                                    isymm > 0 use symmetry
      !                                    isymm = 0 do not use symmetry
      !                 ghostsize    (F) : the size of the ghost layer
      !                 maptype      (I) : the mapping type:
      !                                       ppm_param_map_ghost_put
      !                                       ppm_param_map_ghost_get
      !                                       ppm_param_map_send
      !                                       ppm_param_map_push 
      !                                       ppm_param_map_pop
      !                                       ppm_param_map_cancel
      !                                    
      !                 pushpp       (L) : OPTIONAL argument passed to
      !                                    ppm_map_part_push when xp are pushed
      !                                    
      !  Input/output : Mpart        (I) : the total number of particles 
      !                                    including the ghosts
      !                 info         (I) : return status, 0 on success
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_part_ghost.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.14  2006/10/10 20:40:40  walther
      !  Added the optional argument: pushpp.
      !
      !  Revision 1.13  2006/09/04 18:34:50  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.12  2006/04/06 14:55:19  walther
      !  The arguments to ppm_map_part_ghost_put() changed; 
      !  so updated the call to the routine.
      !
      !  Revision 1.11  2006/02/03 09:41:26  ivos
      !  Added the PRELIMINARY ghost_put functionality. Still needs clean-up,
      !  but should work.
      !
      !  Revision 1.10  2005/03/10 01:40:33  ivos
      !  Empty buffer is now only reported for ppm_debug.GT.0 and is no
      !  longer logged to prevent huge log files.
      !
      !  Revision 1.9  2004/10/01 16:09:06  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.8  2004/07/26 07:42:44  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.7  2004/07/20 06:44:56  walther
      !  Bug fix: added AUX versions for the LOGICAL version of the routine.
      !
      !  Revision 1.6  2004/07/19 07:41:53  ivos
      !  Overloaded particle push and pop operations for 1d data arrays.
      !  Added new routines to module and needed checks to interface
      !  routine.
      !
      !  Revision 1.5  2004/07/16 14:46:28  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.4  2004/07/16 14:07:59  ivos
      !  Added sequence and argument checks. New checks now allow multiple
      !  push-send-pop cycles per mapping.
      !
      !  Revision 1.3  2004/05/28 10:28:06  walther
      !  First functional release - without symmetry.
      !
      !  Revision 1.2  2004/03/22 09:43:16  walther
      !  Adding a call to ppm_util_commopt.
      !
      !  Revision 1.1  2004/02/19 15:56:35  walther
      !  Initial implementation (not complete!).
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if    __DIM == 1
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_ghost_1ds(pdata,Npart,Mpart,isymm,ghostsize, & 
     &                                  maptype,info,pushpp)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_ghost_1dd(pdata,Npart,Mpart,isymm,ghostsize, & 
     &                                  maptype,info,pushpp)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_ghost_1dsc(pdata,Npart,Mpart,isymm,ghostsize, & 
     &                                   maptype,info,pushpp)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_ghost_1ddc(pdata,Npart,Mpart,isymm,ghostsize, & 
     &                                   maptype,info,pushpp)
#elif  __KIND == __INTEGER
#if    __KIND_AUX == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_ghost_1dis(pdata,Npart,Mpart,isymm,ghostsize, & 
     &                                   maptype,info,pushpp)
#else
      SUBROUTINE ppm_map_part_ghost_1did(pdata,Npart,Mpart,isymm,ghostsize, & 
     &                                   maptype,info,pushpp)
#endif
#elif  __KIND == __LOGICAL
#if    __KIND_AUX == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_ghost_1dls(pdata,Npart,Mpart,isymm,ghostsize, & 
     &                                  maptype,info,pushpp)
#else
      SUBROUTINE ppm_map_part_ghost_1dld(pdata,Npart,Mpart,isymm,ghostsize, & 
     &                                  maptype,info,pushpp)
#endif 
#endif 

#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_ghost_2ds(pdata,lda,Npart,Mpart,isymm,ghostsize, & 
     &                                  maptype,info,pushpp)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_ghost_2dd(pdata,lda,Npart,Mpart,isymm,ghostsize, & 
     &                                  maptype,info,pushpp)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_ghost_2dsc(pdata,lda,Npart,Mpart,isymm,ghostsize, & 
     &                                   maptype,info,pushpp)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_ghost_2ddc(pdata,lda,Npart,Mpart,isymm,ghostsize, & 
     &                                   maptype,info,pushpp)
#elif  __KIND == __INTEGER
#if    __KIND_AUX == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_ghost_2dis(pdata,lda,Npart,Mpart,isymm,ghostsize, & 
     &                                   maptype,info,pushpp)
#else
      SUBROUTINE ppm_map_part_ghost_2did(pdata,lda,Npart,Mpart,isymm,ghostsize, & 
     &                                   maptype,info,pushpp)
#endif 
#elif  __KIND == __LOGICAL
#if    __KIND_AUX == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_ghost_2dls(pdata,lda,Npart,Mpart,isymm,ghostsize, & 
     &                                  maptype,info,pushpp)
#else
      SUBROUTINE ppm_map_part_ghost_2dld(pdata,lda,Npart,Mpart,isymm,ghostsize, & 
     &                                  maptype,info,pushpp)
#endif
#endif 
#endif 

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_part, ONLY: ppm_target_topoid
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_util_commopt
      USE ppm_module_write
      USE ppm_module_map_part_ghost_get
      USE ppm_module_map_part_ghost_put
      USE ppm_module_map_part_push
      USE ppm_module_map_part_send
      USE ppm_module_map_part_pop
      USE ppm_module_map_part_ghost_pop
      IMPLICIT NONE
#if     __KIND == __SINGLE_PRECISION | \
        __KIND == __SINGLE_PRECISION_COMPLEX | \
    __KIND_AUX == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if    __DIM == 1
#if    __KIND == __INTEGER
      INTEGER , DIMENSION(:  )  , POINTER       :: pdata
#elif  __KIND == __LOGICAL
      LOGICAL , DIMENSION(:  )  , POINTER       :: pdata
#elif  __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:  )  , POINTER    :: pdata
#else
      REAL(MK), DIMENSION(:  )  , POINTER       :: pdata
#endif

#elif  __DIM == 2
#if    __KIND == __INTEGER
      INTEGER , DIMENSION(:,:)  , POINTER       :: pdata
#elif  __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:)  , POINTER       :: pdata
#elif  __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:)  , POINTER    :: pdata
#else
      REAL(MK), DIMENSION(:,:)  , POINTER       :: pdata
#endif
#endif
      REAL(MK)                  , INTENT(IN   ) :: ghostsize
      INTEGER                   , INTENT(IN   ) :: Npart,maptype
#if    __DIM == 2
      INTEGER                   , INTENT(IN   ) :: lda
#endif
      LOGICAL , OPTIONAL        , INTENT(IN   ) :: pushpp
      INTEGER                   , INTENT(IN   ) :: isymm
      INTEGER                   , INTENT(  OUT) :: Mpart,info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
#if    __DIM == 1
      INTEGER, PARAMETER    :: lda = 1
#endif
      INTEGER               :: i,j,k
      CHARACTER(ppm_char)   :: mesg
      REAL(MK)              :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part_ghost',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_map_part_ghost',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (ppm_topoid .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_no_topo,'ppm_map_part_ghost',  &
     &            'No topology has been defined so far',__LINE__,info)
              GOTO 9999
          ENDIF
#if    __DIM == 2
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_ghost',  &
     &            'lda must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
          IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_ghost',  &
     &            'Npart must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (ghostsize .LT. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_ghost',  &
     &            'ghostsize must be >=0.0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Select the proper mapping routine for the ghosts
      !-------------------------------------------------------------------------
      IF     (maptype.EQ.ppm_param_map_ghost_get) THEN
         !----------------------------------------------------------------------
         !  Create the new ghosts 
         !----------------------------------------------------------------------
#if __KIND == __INTEGER | \
    __KIND == __LOGICAL | \
    __KIND == __SINGLE_PRECISION_COMPLEX | \
    __KIND == __DOUBLE_PRECISION_COMPLEX
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_type,'ppm_map_part_ghost',  &
     &       'Wrong type passed for mapping of ghosts',__LINE__,info)
         GOTO 9999
#else
#if __DIM == 1
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_type,'ppm_map_part_ghost',  &
     &       'Particle positions must be a 2d array',__LINE__,info)
         GOTO 9999
#else

         ! Store destination topoid
         IF (ppm_topoid .LT. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_no_topo,'ppm_map_part_ghost',  &
     &           'No topology has been defined so far',__LINE__,info)
             GOTO 9999
         ENDIF
         ppm_target_topoid = ppm_topoid

         ! if there is still some data left in the buffer, warn the user
         IF (ppm_buffer_set .GT. 0) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_map_incomp,'ppm_map_part_ghost',  &
     &           'Buffer was not empty. Possible loss of data!',__LINE__,info)
         ENDIF

         !----------------------------------------------------------------------
         !  first check if the optimal communication protocol is known
         !----------------------------------------------------------------------
         IF (.NOT.ppm_isoptimized(ppm_topoid)) THEN
            !-------------------------------------------------------------------
            !  if not: determine it before calling map_part_ghost_get
            !-------------------------------------------------------------------
            CALL ppm_util_commopt(ppm_topoid,info)
            IF (ppm_debug .GT. 1) THEN
               DO i=1,ppm_nneighlist(ppm_topoid)
                   WRITE(mesg,'(A,I4)') 'have neighbor: ',   &
     &                 ppm_ineighlist(i,ppm_topoid)
                   CALL ppm_write(ppm_rank,'ppm_map_part_ghost',mesg,info)
               ENDDO
               DO i=1,ppm_ncommseq(ppm_topoid)
                   WRITE(mesg,'(A,I4)') 'communicate: ',   &
     &                 ppm_icommseq(i,ppm_topoid)
                   CALL ppm_write(ppm_rank,'ppm_map_part_ghost',mesg,info)
               ENDDO
            ENDIF
         ENDIF

         !----------------------------------------------------------------------
         !  Map the ghost onto the subs and initialize the send buffers with 
         !  ghost particles
         !----------------------------------------------------------------------
         CALL ppm_map_part_ghost_get(pdata,lda,Npart,isymm,ghostsize,info)
         IF (info.NE.0) GOTO 9999
#endif
#endif
      ELSEIF (maptype.EQ.ppm_param_map_ghost_put) THEN
         !----------------------------------------------------------------------
         !  send the contribution computed on the ghosts back to their host 
         !  processor
         !----------------------------------------------------------------------
#if __KIND == __INTEGER | \
    __KIND == __LOGICAL | \
    __KIND == __SINGLE_PRECISION_COMPLEX | \
    __KIND == __DOUBLE_PRECISION_COMPLEX
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_type,'ppm_map_part_ghost',  &
     &       'Wrong type passed for mapping of ghosts',__LINE__,info)
         GOTO 9999
#else
#if __DIM == 1
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_type,'ppm_map_part_ghost',  &
     &       'Particle positions must be a 2d array',__LINE__,info)
         GOTO 9999
#else

         ! Store destination topoid
         IF (ppm_topoid .LT. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_no_topo,'ppm_map_part_ghost',  &
     &           'No topology has been defined so far',__LINE__,info)
             GOTO 9999
         ENDIF
         ppm_target_topoid = ppm_topoid

         ! if there is still some data left in the buffer, warn the user
         IF (ppm_buffer_set .GT. 0) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_map_incomp,'ppm_map_part_ghost',  &
     &           'Buffer was not empty. Possible loss of data!',__LINE__,info)
         ENDIF

         !----------------------------------------------------------------------
         !  first check if the optimal communication protocol is known
         !----------------------------------------------------------------------
         IF (.NOT.ppm_isoptimized(ppm_topoid)) THEN
            !-------------------------------------------------------------------
            !  if not: determine it before calling map_part_ghost_put
            !-------------------------------------------------------------------
            CALL ppm_util_commopt(ppm_topoid,info)
            IF (ppm_debug .GT. 1) THEN
               DO i=1,ppm_nneighlist(ppm_topoid)
                   WRITE(mesg,'(A,I4)') 'have neighbor: ',   &
     &                 ppm_ineighlist(i,ppm_topoid)
                   CALL ppm_write(ppm_rank,'ppm_map_part_ghost',mesg,info)
               ENDDO
               DO i=1,ppm_ncommseq(ppm_topoid)
                   WRITE(mesg,'(A,I4)') 'communicate: ',   &
     &                 ppm_icommseq(i,ppm_topoid)
                   CALL ppm_write(ppm_rank,'ppm_map_part_ghost',mesg,info)
               ENDDO
            ENDIF
         ENDIF

         !----------------------------------------------------------------------
         !  now, call the ghost_put to shift the send/recv lists to enable 
         !  subseq. push/send/pop of the ghost data
         !----------------------------------------------------------------------
         CALL ppm_map_part_ghost_put(info)
         IF (info.NE.0) GOTO 9999
#endif
#endif
      ELSEIF (maptype.EQ.ppm_param_map_push) THEN
         !----------------------------------------------------------------------
         !  push data onto these ghosts
         !----------------------------------------------------------------------
         ! warn the user if no mapping has yet been defined
         IF (ppm_target_topoid .LT. 0) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_nomap,'ppm_map_part_ghost',  &
     &           'Skipping push!',__LINE__,info)
             GOTO 9999
         ENDIF

         !----------------------------------------------------------------------
         !  when we push we do it with the optional argument
         !----------------------------------------------------------------------
         IF (PRESENT(pushpp)) THEN
            CALL ppm_map_part_push(pdata,lda,Npart,info,pushpp)
         ELSE
            CALL ppm_map_part_push(pdata,lda,Npart,info)
         ENDIF  
         IF (info.NE.0) GOTO 9999

      ELSEIF (maptype.EQ.ppm_param_map_pop) THEN
         !----------------------------------------------------------------------
         !  pop the particle data from the buffer
         !----------------------------------------------------------------------
         ! warn the user if no mapping has yet been defined
         IF (ppm_target_topoid .LT. 0) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_nomap,'ppm_map_part_ghost',  &
     &           'Skipping pop!',__LINE__,info)
             GOTO 9999
         ENDIF

         ! skip if the buffer is empty
         IF (ppm_buffer_set .LT. 1) THEN
             IF (ppm_debug .GT. 1) THEN
                 CALL ppm_write(ppm_rank,'ppm_map_part_ghost',  &
     &               'Buffer is empty: skipping pop!',info)
             ENDIF
             GOTO 9999
         ENDIF

         IF (ppm_map_type .EQ. ppm_param_map_ghost_put) THEN
             CALL ppm_map_part_ghost_pop(pdata,lda,Npart,Mpart,info)
         ELSE
             CALL ppm_map_part_pop(pdata,lda,Npart,Mpart,info)
         ENDIF
         IF (info.NE.0) GOTO 9999

      ELSEIF (maptype.EQ.ppm_param_map_send) THEN
         !----------------------------------------------------------------------
         !  send and receive the ghosts - only after this step will we know
         !  the total number of particle - including the ghosts
         !----------------------------------------------------------------------
         ! warn the user if no mapping has yet been defined
         IF (ppm_target_topoid .LT. 0) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_nomap,'ppm_map_part_ghost',  &
     &           'Skipping send!',__LINE__,info)
             GOTO 9999
         ENDIF

         ! skip if the buffer is empty
         IF (ppm_buffer_set .LT. 1) THEN
             IF (ppm_debug .GT. 1) THEN
                 CALL ppm_write(ppm_rank,'ppm_map_part_ghost',  &
     &               'Buffer is empty: skipping send!',info)
             ENDIF
             GOTO 9999
         ENDIF
         CALL ppm_map_part_send(Npart,Mpart,info)
         IF (info.NE.0) GOTO 9999

      ELSEIF (maptype.EQ.ppm_param_map_cancel) THEN
         !----------------------------------------------------------------------
         !  Cancel the mapping in progress and reset everything. No need to
         !  reset arrays since thay are alloc_fit-ed the next time a
         !  mapping is called. Their deallocation only happens in
         !  ppm_finalize anyway.
         !----------------------------------------------------------------------
         ppm_target_topoid = -1
         !----------------------------------------------------------------------
         !  These might actually not be needed since the mapping routines
         !  would in principle be expected to reset them. They are here for
         !  safety reasons.
         !----------------------------------------------------------------------
         ppm_buffer_set    = 0
         ppm_nsendbuffer   = 0
         ppm_nrecvbuffer   = 0

      ELSE
         !----------------------------------------------------------------------
         !  Unknow mapping
         !----------------------------------------------------------------------
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_map_part',    &
     &       'Unknown mapping action specified',__LINE__,info)
         GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_part_ghost',t0,info)
      RETURN
#if    __DIM == 1
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_ghost_1ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_ghost_1dd
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_ghost_1dsc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_ghost_1ddc
#elif  __KIND == __LOGICAL
#if __KIND_AUX == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_ghost_1dls
#else
      END SUBROUTINE ppm_map_part_ghost_1dld
#endif
#elif  __KIND == __INTEGER
#if __KIND_AUX == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_ghost_1dis
#else
      END SUBROUTINE ppm_map_part_ghost_1did
#endif
#endif

#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_ghost_2ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_ghost_2dd
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_ghost_2dsc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_ghost_2ddc
#elif  __KIND == __LOGICAL
#if __KIND_AUX == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_ghost_2dls
#else
      END SUBROUTINE ppm_map_part_ghost_2dld
#endif
#elif  __KIND == __INTEGER
#if __KIND_AUX == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_ghost_2dis
#else
      END SUBROUTINE ppm_map_part_ghost_2did
#endif
#endif
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_map_part_pop
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine pops the contents of the receive buffer.
      !
      !  Input        : pdata(:[,:])(O) particle data to be popped.
      !                                 Overloaded types: single, double,
      !                                 integer, logical, single complex,
      !                                 double complex. Can be either 1d or
      !                                 2d array.
      !                 lda         (I) the leading dimension of pdata.
      !                                 In the case of pdata being a 1d
      !                                 array use lda=1.
      !                 Npart       (I) the old number of particles 
      !                                 (on the processor)
      !                 Mpart       (I) the new number of particles 
      !                                 (on the processor)
      !
      !  Input/output :
      !
      !  Output       : info        (I) return status. 0 upon success.
      !
      !  Remarks      : The first part of the buffer contains the on processor
      !                 data. The packing could be performed more efficiently.
      !
      !                 Maybe remove the lda from the argument list (since
      !                 it is only used to check against ppm_buffer_dim
      !                 anyway) and check size(pdata,1) against
      !                 ppm_buffer_dim. QUESTION: does the size(,) stuff
      !                 also work if the array was allocated in a C client?
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_part_pop.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.22  2006/09/04 18:34:51  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.21  2005/02/16 05:06:29  ivos
      !  Bugfix in the unrolled versions.
      !
      !  Revision 1.20  2005/02/16 02:17:09  ivos
      !  Unrolled loops for 1,2,3,4,5 dimensional vectors to allow
      !  vectorization. Only vector mappings with lda.GT.5 do not
      !  vectorize. Scalar mappings do. Tested on NEC SX-5.
      !
      !  Revision 1.19  2004/10/01 16:09:07  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.18  2004/09/22 11:25:49  ivos
      !  bugfix: fixed broken conversions from single to double and
      !  other way round.
      !
      !  Revision 1.17  2004/09/21 15:19:07  hiebers
      !  erased quotation mark in revision comment
      !  : ----------------------------------------------------------------------
      !
      !  Revision 1.16  2004/09/21 14:38:27  ivos
      !  Added deallocate statements for ppm_sendbuffer and ppm_recvbuffer
      !  in order to save memory (Simone suggestion).
      !
      !  Revision 1.15  2004/07/26 07:42:45  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.14  2004/07/19 07:41:54  ivos
      !  Overloaded particle push and pop operations for 1d data arrays.
      !  Added new routines to module and needed checks to interface
      !  routine.
      !
      !  Revision 1.13  2004/07/15 16:22:04  ivos
      !  nsendbuffer is now reset to allow repeated push-send-pop cycles
      !  with the same mapping lists.
      !
      !  Revision 1.12  2004/06/10 16:20:01  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.11  2004/05/28 10:30:22  walther
      !  Bug fix: ppm_param_alloc_grow_preserve should be used when 
      !  reallocating the pdata. Now also support for ghosts.
      !
      !  Revision 1.10  2004/05/17 15:49:00  oingo
      !  Changed the check for Mpart from .LE. 0 to .LT. 0
      !
      !  Revision 1.9  2004/04/05 11:59:57  ivos
      !  Changed conversion from float to logical to use ppm_myeps.
      !
      !  Revision 1.8  2004/03/05 13:43:23  ivos
      !  bugfix: index errors corrected for COMPLEX version. COMPLEX version
      !  is now tested.
      !
      !  Revision 1.7  2004/02/24 17:12:40  ivos
      !  Added overloaded versions for single complex and double complex 
      !  particle data.
      !
      !  Revision 1.6  2004/02/20 14:58:33  walther
      !  Now checking that the type of the buffer matches the user type.
      !
      !  Revision 1.5  2004/01/23 17:24:16  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.4  2004/01/23 11:31:22  ivos
      !  Cleanup: (1) updated headers, (2) inserted ppm_error and ppm_write,
      !  (3) added argument checking, (4) added checks after every alloc.
      !
      !  Revision 1.3  2003/12/16 15:23:50  ivos
      !  bug fix: parallel version now also works in the limit case of 
      !  only 1 CPU.
      !
      !  Revision 1.2  2003/12/12 17:21:38  hiebers
      !  Bugfix: removed loop over processors in serial version
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      ! ATTN: DIM is the dimension of the pdata array and not the space
      ! dimension ppm_dim!
#if    __DIM ==1
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_pop_1ds(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_pop_1dd(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_pop_1dsc(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_pop_1ddc(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_part_pop_1di(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_part_pop_1dl(pdata,lda,Npart,Mpart,info)
#endif 

#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_pop_2ds(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_pop_2dd(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_pop_2dsc(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_pop_2ddc(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_part_pop_2di(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_part_pop_2dl(pdata,lda,Npart,Mpart,info)
#endif 
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if   __DIM == 1
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:  ), POINTER    :: pdata
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:  ), POINTER    :: pdata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:  ), POINTER :: pdata
#else
      REAL(MK), DIMENSION(:  ), POINTER    :: pdata
#endif

#elif __DIM == 2
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:,:), POINTER    :: pdata
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:), POINTER    :: pdata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:), POINTER :: pdata
#else
      REAL(MK), DIMENSION(:,:), POINTER    :: pdata
#endif
#endif
      INTEGER                 , INTENT(IN   ) :: lda,Npart,Mpart
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3) :: ldu
      INTEGER               :: k,ipart,bdim,ibuffer,btype
      INTEGER               :: iopt,edim,istart
      CHARACTER(ppm_char)   :: mesg
      REAL(MK)              :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part_pop',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (ppm_buffer_set .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_pop',  &
     &            'buffer is empty. Cannot pop.',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_pop',  &
     &            'Npart must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (Mpart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_pop',  &
     &            'Mpart must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
#if   __DIM == 2
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_pop',  &
     &            'lda must be >0 for vector data',__LINE__,info)
              GOTO 9999
          ENDIF
#elif __DIM == 1
          IF (lda .NE. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_pop',  &
     &            'lda must be =1 for scalar data',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that the required dimension fits the dimension of the buffer
      !-------------------------------------------------------------------------
      bdim = ppm_buffer_dim(ppm_buffer_set)
#if   __KIND == __SINGLE_PRECISION_COMPLEX | \
      __KIND == __DOUBLE_PRECISION_COMPLEX
      ! for complex, the effective dimension is half the data dimension
      edim = bdim/2
#else
      edim = bdim
#endif
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(2(A,I3))') 'bdim=',edim,'    lda=',lda
          CALL ppm_write(ppm_rank,'ppm_map_part_pop',mesg,info)
      ENDIF
#if   __DIM == 2
      IF (edim.NE.lda) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_dim,'ppm_map_part_pop',    &
     &       'leading dimension LDA is in error',__LINE__,info)
         GOTO 9999
      ENDIF 
#elif __DIM == 1
      IF (edim.NE.1) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_dim,'ppm_map_part_pop',    &
     &       'buffer does not contain 1d data!',__LINE__,info)
         GOTO 9999
      ENDIF 
#endif

      !-------------------------------------------------------------------------
      !  Check that the required type is identical to the type of the buffer
      !-------------------------------------------------------------------------
      btype = ppm_buffer_type(ppm_buffer_set)
#if    __KIND == __SINGLE_PRECISION
      IF (btype.NE.ppm_kind_single) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_part_pop',    &
     &       'trying to pop a non-single into single ',__LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __DOUBLE_PRECISION
      IF (btype.NE.ppm_kind_double) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_part_pop',    &
     &       'trying to pop a non-double into double ',__LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      IF (btype.NE.ppm_kind_single) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_part_pop',    &
     &       'trying to pop a non-single-complex into single-complex',  &
     &       __LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      IF (btype.NE.ppm_kind_double) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_part_pop',    &
     &       'trying to pop a non-double-complex into double-complex',  &
     &       __LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __INTEGER
      IF (btype.NE.ppm_integer) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_part_pop',    &
     &       'trying to pop a non-integer into integer ',__LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __LOGICAL
      IF (btype.NE.ppm_logical) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_part_pop',    &
     &       'trying to pop a non-logical into logical ',__LINE__,info)
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  (Re)allocate the particle data (if necessary)
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
#if   __DIM == 2 
      ldu(1) = edim
      ldu(2) = Mpart
#elif __DIM == 1
      ldu(1) = Mpart
#endif
      CALL ppm_alloc(pdata,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_pop',     &
     &        'particle data PDATA',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Decrement the pointer into the receive buffer
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(2(A,I9))') 'ppm_nrecvbuffer = ',ppm_nrecvbuffer,   &
     &        'Mpart*bdim = ',Mpart*bdim
          CALL ppm_write(ppm_rank,'ppm_map_part_pop',mesg,info)
      ENDIF
      IF (ppm_map_type.EQ.ppm_param_map_ghost_get) THEN
         ppm_nrecvbuffer = ppm_nrecvbuffer - (Mpart - Npart)*bdim 
      ELSE
         ppm_nrecvbuffer = ppm_nrecvbuffer - Mpart*bdim 
      ENDIF 

      !-------------------------------------------------------------------------
      !  Decrement the pointer into the send buffer to allow reuse by
      !  multiple sequential push-send-pop cycles.
      !-------------------------------------------------------------------------
      ppm_nsendbuffer = ppm_nsendbuffer - ppm_buffer_dim(ppm_buffer_set)*  &
     &    (ppm_psendbuffer(ppm_nsendlist+1)-1)

      !-------------------------------------------------------------------------
      !  loop over the processors in the ppm_isendlist() 
      !-------------------------------------------------------------------------
#ifdef __MPI
      ibuffer = ppm_nrecvbuffer 

      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(A,I9)') 'ibuffer = ',ibuffer
          CALL ppm_write(ppm_rank,'ppm_map_part_pop',mesg,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  compute the start of the loop: when ghosts are popped, we need to add
      !  them at the end of the current particle list, otherwise we overwrite  
      !  the current particle list
      !-------------------------------------------------------------------------
      IF (ppm_map_type.EQ.ppm_param_map_ghost_get) THEN
         istart = Npart + 1
      ELSE
         istart = 1
      ENDIF 
      !-------------------------------------------------------------------------
      !  DOUBLE PRECISION BUFFER
      !-------------------------------------------------------------------------
      IF (ppm_kind.EQ.ppm_kind_double) THEN
#if    __DIM == 2
         !----------------------------------------------------------------------
         !  Unrolled version of edim=1
         !----------------------------------------------------------------------
         IF (edim .EQ. 1) THEN
            DO ipart=istart,Mpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbufferd(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
#endif
            ENDDO
         !----------------------------------------------------------------------
         !  Unrolled version of edim=2
         !----------------------------------------------------------------------
         ELSEIF (edim .EQ. 2) THEN
            DO ipart=istart,Mpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = ppm_recvbufferd(ibuffer)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbufferd(ibuffer))
               ibuffer = ibuffer + 1
               pdata(2,ipart) = INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbufferd(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbufferd(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(2,ipart) = .TRUE.
               ELSE
                  pdata(2,ipart) = .FALSE.
               ENDIF 
#endif
            ENDDO
         !----------------------------------------------------------------------
         !  Unrolled version of edim=3
         !----------------------------------------------------------------------
         ELSEIF (edim .EQ. 3) THEN
            DO ipart=istart,Mpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = ppm_recvbufferd(ibuffer)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = ppm_recvbufferd(ibuffer)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbufferd(ibuffer))
               ibuffer = ibuffer + 1
               pdata(2,ipart) = INT(ppm_recvbufferd(ibuffer))
               ibuffer = ibuffer + 1
               pdata(3,ipart) = INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbufferd(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbufferd(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(2,ipart) = .TRUE.
               ELSE
                  pdata(2,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbufferd(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(3,ipart) = .TRUE.
               ELSE
                  pdata(3,ipart) = .FALSE.
               ENDIF 
#endif
            ENDDO
         !----------------------------------------------------------------------
         !  Unrolled version of edim=4
         !----------------------------------------------------------------------
         ELSEIF (edim .EQ. 4) THEN
            DO ipart=istart,Mpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(4,ipart) = REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = ppm_recvbufferd(ibuffer)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = ppm_recvbufferd(ibuffer)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = ppm_recvbufferd(ibuffer)
               ibuffer = ibuffer + 1
               pdata(4,ipart) = ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(4,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(4,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbufferd(ibuffer))
               ibuffer = ibuffer + 1
               pdata(2,ipart) = INT(ppm_recvbufferd(ibuffer))
               ibuffer = ibuffer + 1
               pdata(3,ipart) = INT(ppm_recvbufferd(ibuffer))
               ibuffer = ibuffer + 1
               pdata(4,ipart) = INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbufferd(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbufferd(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(2,ipart) = .TRUE.
               ELSE
                  pdata(2,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbufferd(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(3,ipart) = .TRUE.
               ELSE
                  pdata(3,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbufferd(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(4,ipart) = .TRUE.
               ELSE
                  pdata(4,ipart) = .FALSE.
               ENDIF 
#endif
            ENDDO
         !----------------------------------------------------------------------
         !  Unrolled version of edim=5
         !----------------------------------------------------------------------
         ELSEIF (edim .EQ. 5) THEN
            DO ipart=istart,Mpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(4,ipart) = REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(5,ipart) = REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = ppm_recvbufferd(ibuffer)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = ppm_recvbufferd(ibuffer)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = ppm_recvbufferd(ibuffer)
               ibuffer = ibuffer + 1
               pdata(4,ipart) = ppm_recvbufferd(ibuffer)
               ibuffer = ibuffer + 1
               pdata(5,ipart) = ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(4,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(5,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(4,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(5,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &             ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbufferd(ibuffer))
               ibuffer = ibuffer + 1
               pdata(2,ipart) = INT(ppm_recvbufferd(ibuffer))
               ibuffer = ibuffer + 1
               pdata(3,ipart) = INT(ppm_recvbufferd(ibuffer))
               ibuffer = ibuffer + 1
               pdata(4,ipart) = INT(ppm_recvbufferd(ibuffer))
               ibuffer = ibuffer + 1
               pdata(5,ipart) = INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbufferd(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbufferd(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(2,ipart) = .TRUE.
               ELSE
                  pdata(2,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbufferd(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(3,ipart) = .TRUE.
               ELSE
                  pdata(3,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbufferd(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(4,ipart) = .TRUE.
               ELSE
                  pdata(4,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbufferd(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(5,ipart) = .TRUE.
               ELSE
                  pdata(5,ipart) = .FALSE.
               ENDIF 
#endif
            ENDDO
         !----------------------------------------------------------------------
         !  For edim.GT.5 the vector length will be edim !!
         !----------------------------------------------------------------------
         ELSE
            DO ipart=istart,Mpart
               DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
                  ibuffer = ibuffer + 2
#else
                  ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                  pdata(k,ipart) = REAL(ppm_recvbufferd(ibuffer),       &
     &                ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                  pdata(k,ipart) = ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  pdata(k,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  pdata(k,ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                  pdata(k,ipart) = INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                  IF (ppm_recvbufferd(ibuffer) .GT.     &
     &               (1.0_ppm_kind_double-ppm_myepsd)) THEN
                     pdata(k,ipart) = .TRUE.
                  ELSE
                     pdata(k,ipart) = .FALSE.
                  ENDIF 
#endif
               ENDDO
            ENDDO
         ENDIF
#elif  __DIM == 1
         !----------------------------------------------------------------------
         !  Scalar version
         !----------------------------------------------------------------------
         DO ipart=istart,Mpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
            ibuffer = ibuffer + 2
#else
            ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
            pdata(ipart) = REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
            pdata(ipart) = ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
            pdata(ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &          ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
            pdata(ipart) = CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &          ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
            pdata(ipart) = INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
            IF (ppm_recvbufferd(ibuffer) .GT.     &
     &         (1.0_ppm_kind_double-ppm_myepsd)) THEN
               pdata(ipart) = .TRUE.
            ELSE
               pdata(ipart) = .FALSE.
            ENDIF 
#endif
         ENDDO
#endif
      !-------------------------------------------------------------------------
      !  SINGLE PRECISION BUFFER
      !-------------------------------------------------------------------------
      ELSE
#if    __DIM == 2
         !----------------------------------------------------------------------
         !  Unrolled verion for edim=1
         !----------------------------------------------------------------------
         IF (edim .EQ. 1) THEN
            DO ipart=istart,Mpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = ppm_recvbuffers(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbuffers(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
#endif
            ENDDO
         !----------------------------------------------------------------------
         !  Unrolled verion for edim=2
         !----------------------------------------------------------------------
         ELSEIF (edim .EQ. 2) THEN
            DO ipart=istart,Mpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = ppm_recvbuffers(ibuffer)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = ppm_recvbuffers(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(1,ipart) = REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbuffers(ibuffer))
               ibuffer = ibuffer + 1
               pdata(2,ipart) = INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbuffers(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffers(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(2,ipart) = .TRUE.
               ELSE
                  pdata(2,ipart) = .FALSE.
               ENDIF 
#endif
            ENDDO
         !----------------------------------------------------------------------
         !  Unrolled verion for edim=3
         !----------------------------------------------------------------------
         ELSEIF (edim .EQ. 3) THEN
            DO ipart=istart,Mpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = ppm_recvbuffers(ibuffer)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = ppm_recvbuffers(ibuffer)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = ppm_recvbuffers(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbuffers(ibuffer))
               ibuffer = ibuffer + 1
               pdata(2,ipart) = INT(ppm_recvbuffers(ibuffer))
               ibuffer = ibuffer + 1
               pdata(3,ipart) = INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbuffers(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffers(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(2,ipart) = .TRUE.
               ELSE
                  pdata(2,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffers(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(3,ipart) = .TRUE.
               ELSE
                  pdata(3,ipart) = .FALSE.
               ENDIF 
#endif
            ENDDO
         !----------------------------------------------------------------------
         !  Unrolled verion for edim=4
         !----------------------------------------------------------------------
         ELSEIF (edim .EQ. 4) THEN
            DO ipart=istart,Mpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = ppm_recvbuffers(ibuffer)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = ppm_recvbuffers(ibuffer)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = ppm_recvbuffers(ibuffer)
               ibuffer = ibuffer + 1
               pdata(4,ipart) = ppm_recvbuffers(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(4,ipart) = REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(4,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(4,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbuffers(ibuffer))
               ibuffer = ibuffer + 1
               pdata(2,ipart) = INT(ppm_recvbuffers(ibuffer))
               ibuffer = ibuffer + 1
               pdata(3,ipart) = INT(ppm_recvbuffers(ibuffer))
               ibuffer = ibuffer + 1
               pdata(4,ipart) = INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbuffers(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffers(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(2,ipart) = .TRUE.
               ELSE
                  pdata(2,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffers(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(3,ipart) = .TRUE.
               ELSE
                  pdata(3,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffers(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(4,ipart) = .TRUE.
               ELSE
                  pdata(4,ipart) = .FALSE.
               ENDIF 
#endif
            ENDDO
         !----------------------------------------------------------------------
         !  Unrolled verion for edim=5
         !----------------------------------------------------------------------
         ELSEIF (edim .EQ. 5) THEN
            DO ipart=istart,Mpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = ppm_recvbuffers(ibuffer)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = ppm_recvbuffers(ibuffer)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = ppm_recvbuffers(ibuffer)
               ibuffer = ibuffer + 1
               pdata(4,ipart) = ppm_recvbuffers(ibuffer)
               ibuffer = ibuffer + 1
               pdata(5,ipart) = ppm_recvbuffers(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(4,ipart) = REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(5,ipart) = REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(4,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(5,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(4,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(5,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &             ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbuffers(ibuffer))
               ibuffer = ibuffer + 1
               pdata(2,ipart) = INT(ppm_recvbuffers(ibuffer))
               ibuffer = ibuffer + 1
               pdata(3,ipart) = INT(ppm_recvbuffers(ibuffer))
               ibuffer = ibuffer + 1
               pdata(4,ipart) = INT(ppm_recvbuffers(ibuffer))
               ibuffer = ibuffer + 1
               pdata(5,ipart) = INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbuffers(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffers(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(2,ipart) = .TRUE.
               ELSE
                  pdata(2,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffers(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(3,ipart) = .TRUE.
               ELSE
                  pdata(3,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffers(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(4,ipart) = .TRUE.
               ELSE
                  pdata(4,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffers(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(5,ipart) = .TRUE.
               ELSE
                  pdata(5,ipart) = .FALSE.
               ENDIF 
#endif
            ENDDO
         !----------------------------------------------------------------------
         !  For edim.GT.5 the vector length will be edim !!
         !----------------------------------------------------------------------
         ELSE
            DO ipart=istart,Mpart
               DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
                  ibuffer = ibuffer + 2
#else
                  ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                  pdata(k,ipart) = ppm_recvbuffers(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
                  pdata(k,ipart) = REAL(ppm_recvbuffers(ibuffer),       &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  pdata(k,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  pdata(k,ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
                  pdata(k,ipart) = INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                  IF (ppm_recvbuffers(ibuffer) .GT.      &
     &               (1.0_ppm_kind_single-ppm_myepss)) THEN
                     pdata(k,ipart) = .TRUE.
                  ELSE
                     pdata(k,ipart) = .FALSE.
                  ENDIF 
#endif
               ENDDO
            ENDDO
         ENDIF
#elif  __DIM == 1
         !----------------------------------------------------------------------
         !  Scalar version
         !----------------------------------------------------------------------
         DO ipart=istart,Mpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
            ibuffer = ibuffer + 2
#else
            ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
            pdata(ipart) = ppm_recvbuffers(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
            pdata(ipart) = REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
            pdata(ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &          ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
            pdata(ipart) = CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &          ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
            pdata(ipart) = INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
            IF (ppm_recvbuffers(ibuffer) .GT.      &
     &         (1.0_ppm_kind_single-ppm_myepss)) THEN
               pdata(ipart) = .TRUE.
            ELSE
               pdata(ipart) = .FALSE.
            ENDIF 
#endif
         ENDDO
#endif
      ENDIF       ! ppm_kind
#endif

      !-------------------------------------------------------------------------
      !  Decrement the set counter
      !-------------------------------------------------------------------------
      ppm_buffer_set = ppm_buffer_set - 1  

      !-------------------------------------------------------------------------
      !  Deallocate the receive buffer if all sets have been poped
      !-------------------------------------------------------------------------
      IF (ppm_buffer_set .LT. 1) THEN
          iopt = ppm_param_dealloc
          IF (ppm_kind .EQ. ppm_kind_single) THEN
              CALL ppm_alloc(ppm_recvbuffers,ldu,iopt,info)
          ELSE
              CALL ppm_alloc(ppm_recvbufferd,ldu,iopt,info)
          ENDIF
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_alloc,'ppm_map_part_pop',     &
     &            'receive buffer PPM_RECVBUFFER',__LINE__,info)
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_part_pop',t0,info)
      RETURN
#if    __DIM == 1
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_pop_1ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_pop_1dd
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_pop_1dsc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_pop_1ddc
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_map_part_pop_1di
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_map_part_pop_1dl
#endif

#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_pop_2ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_pop_2dd
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_pop_2dsc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_pop_2ddc
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_map_part_pop_2di
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_map_part_pop_2dl
#endif
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_map_part_ring_shift
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Prepares the send and receive lists and the buffers for
      !                 sending the particle data along the ring
      !
      !  Input        : xp(:,:)  (F) : particle co-ordinates
      !                 lda      (I) : the leading dimension of xp
      !                 Npart    (I) : the number of particles
      !                 itarget  (I) : the processor to send to
      !                 isource  (I) : the processor to receive from
      !
      !  Input/output : info     (I) : return status. 0 on success.
      !
      !  Remarks      : The first part of the buffer contains the on processor
      !                 data. In this case all particles will be send to
      !                 itarget.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_part_ring_shift.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.13  2006/09/04 18:34:52  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.12  2004/10/01 16:09:08  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.11  2004/07/29 12:26:08  oingo
      !  Added new parameter lda for the leading dimension of xp
      !
      !  Revision 1.10  2004/07/26 11:48:08  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.9  2004/07/26 07:42:46  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.8  2004/07/19 11:14:33  oingo
      !  Added check for ppm_initialized and added check if the buffer is empty.
      !  If the buffer is not empty a warning will be raised (possible loss of
      !  data)
      !
      !  Revision 1.7  2004/07/13 10:38:40  oingo
      !  Bugfix: it was assumed that the leading dimension of xp is always equal
      !          to ppm_dim. Now SIZE(xp,1) is used to get the leading dimension
      !
      !  Revision 1.6  2004/06/10 16:20:02  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.5  2004/05/11 10:08:54  oingo
      !  Bugfix: The current topology ID will not be changed anymore by this
      !  routine
      !
      !  Revision 1.4  2004/05/03 12:26:16  oingo
      !  The ppm_target_topoid will now be set to 0 in order to be able to
      !  use the ppm_map_part interface to push, send and pop any further
      !  particle data
      !
      !  Revision 1.3  2004/04/26 10:56:32  oingo
      !  Added more argument checking
      !
      !  Revision 1.2  2004/04/22 11:39:48  oingo
      !  Bugfix: fixed a Fortran keyword conflict with a variable name
      !
      !  Revision 1.1  2004/04/22 08:20:13  oingo
      !  Initial release
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_ring_shift_s(xp,lda,Npart,itarget,isource,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_ring_shift_d(xp,lda,Npart,itarget,isource,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_part, ONLY: ppm_target_topoid
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      INTEGER                 , INTENT(IN   ) :: lda
      INTEGER                 , INTENT(IN   ) :: Npart
      INTEGER                 , INTENT(IN   ) :: itarget
      INTEGER                 , INTENT(IN   ) :: isource
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(1)       :: ldu
      INTEGER                     :: i
      INTEGER                     :: iopt,ibuffer,ipart
      REAL(MK)                    :: t0
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part_ring_shift',t0,info)
      
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         IF (.NOT. ppm_initialized) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_ppm_noinit,'ppm_map_part_ring_shift',  &
     &           'Please call ppm_init first!',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (lda .LT. 1) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,'ppm_map_part_ring_shift',  &
     &           'LDA must be >0',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (Npart .LT. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,'ppm_map_part_ring_shift',  &
     &           'NPART must be >=0',__LINE__,info)
            GOTO 9999
         ENDIF
         IF ((itarget .LT. 0) .OR. (itarget .GE. ppm_nproc)) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,'ppm_map_part_ring_shift',       &
     &           'ITARGET must be >=0 and smaller than the number processors',&
     &           __LINE__,info)
            GOTO 9999
         ENDIF
         IF ((isource .LT. 0) .OR. (isource .GE. ppm_nproc)) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,'ppm_map_part_ring_shift',       &
     &           'ISOURCE must be >=0 and smaller than the number processors',&
     &           __LINE__,info)
            GOTO 9999
         ENDIF
         IF (ppm_target_topoid .GE. 0) THEN
            info = ppm_error_warning
            CALL ppm_error(ppm_err_map_incomp,'ppm_map_part_ring_shift',  &
     &           'Discarded',__LINE__,info)
         ENDIF
         IF (ppm_topoid .LT. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_map_incomp,'ppm_map_part_ring_shift',  &
     &           'A topology must be defined first',__LINE__,info)
            GOTO 9999
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  If there is still some data left in the buffer, warn the user
      !-------------------------------------------------------------------------
      IF (ppm_buffer_set .GT. 0) THEN
         info = ppm_error_warning
         CALL ppm_error(ppm_err_map_incomp,'ppm_map_part_ring_shift',  &
     &        'Buffer was not empty. Possible loss of data!',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the pointer to the buffer; we only need to talk
      !  to 2 processors (including this one)
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = 3
      CALL ppm_alloc(ppm_psendbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring_shift',     &
     &        'particle send buffer PPM_PSENDBUFFER',__LINE__,info)
         GOTO 9999
      ENDIF

      iopt   = ppm_param_alloc_fit
      ldu(1) = Npart
      CALL ppm_alloc(ppm_buffer2part,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring_shift',     &
     &        'buffer-to-particles map PPM_BUFFER2PART',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store the number of buffer entries (this is the first)
      !-------------------------------------------------------------------------
      ppm_buffer_set = 1

      !-------------------------------------------------------------------------
      !  Allocate memory for the field registers that holds the dimension and
      !  type of the data
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_buffer_set
      CALL ppm_alloc(ppm_buffer_dim,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring_shift',     &
     &        'buffer dimensions PPM_BUFFER_DIM',__LINE__,info)
         GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_buffer_type,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring_shift',     &
     &        'buffer types PPM_BUFFER_TYPE',__LINE__,info)
         GOTO 9999
      ENDIF

      ppm_buffer_dim(ppm_buffer_set)  = lda
      IF (ppm_kind .EQ. ppm_kind_single) THEN
         ppm_buffer_type(ppm_buffer_set) = ppm_kind_single
      ELSE
         ppm_buffer_type(ppm_buffer_set) = ppm_kind_double
      ENDIF

      !-------------------------------------------------------------------------
      !  (Re)allocate memory for the buffer
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = lda * Npart
#if    __KIND == __SINGLE_PRECISION
      CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
#else
      CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
#endif

      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring_shift',     &
     &        'global send buffer PPM_SENDBUFFER',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the sendlist
      !-------------------------------------------------------------------------
      ppm_nsendlist = 2
      ppm_nrecvlist = 2

      iopt = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring_shift',     &
     &        'send list PPM_ISENDLIST',__LINE__,info)
         GOTO 9999
      ENDIF

      ldu(1) = ppm_nrecvlist
      CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring_shift',     &
     &        'receive list PPM_IRECVLIST',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Fill the sendbuffer with the particles that will stay on this
      !  processor (no particles will stay because the whole copy has
      !  to be send to the neighbour)
      !-------------------------------------------------------------------------
      ppm_psendbuffer(1) = 1
      ppm_psendbuffer(2) = 1

      !-------------------------------------------------------------------------
      !  Set the target and source processor
      !-------------------------------------------------------------------------
      ppm_isendlist(1) = ppm_rank
      ppm_isendlist(2) = itarget

      ppm_irecvlist(1) = ppm_rank
      ppm_irecvlist(2) = isource
      
      !-------------------------------------------------------------------------
      !  Fill the sendbuffer with the data that will be send away
      !-------------------------------------------------------------------------
      ibuffer = 0
      DO ipart = 1,Npart
         !----------------------------------------------------------------------
         !  Store the id of the data entity
         !----------------------------------------------------------------------
         ppm_buffer2part(ipart) = ipart

         !----------------------------------------------------------------------
         !  Store the data
         !----------------------------------------------------------------------
         DO i = 1,lda
            ibuffer                  = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
            ppm_sendbuffers(ibuffer) = xp(i,ipart)
#else
            ppm_sendbufferd(ibuffer) = xp(i,ipart)
#endif
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Everything has to go!
      !-------------------------------------------------------------------------
      ppm_psendbuffer(3) = Npart + 1

      !-------------------------------------------------------------------------
      !  Store the current size of the buffer
      !-------------------------------------------------------------------------
      ppm_nsendbuffer = ibuffer

      !-------------------------------------------------------------------------
      !  Set the target topoid to the current topology in order to be able to
      !  use the ppm_map_part interface to push, send and pop any further
      !  particle informations
      !------------------------------------------------------------------------
      ppm_target_topoid = ppm_topoid

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_part_ring_shift',t0,info)
      RETURN
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_ring_shift_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_ring_shift_d
#endif

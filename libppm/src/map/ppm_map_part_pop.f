      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_map_part_pop
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------


#if    __DIM ==1
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_pop_1ds(pdata,Npart,Mpart,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_pop_1dd(pdata,Npart,Mpart,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_pop_1dsc(pdata,Npart,Mpart,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_pop_1ddc(pdata,Npart,Mpart,info)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_part_pop_1di(pdata,Npart,Mpart,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_part_pop_1dl(pdata,Npart,Mpart,info)
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
      !!! This routine pops the contents of the receive buffer.
      !!!
      !!! [NOTE]
      !!! The first part of the buffer contains the on processor
      !!! data. The packing could be performed more efficiently.
      !!!
      !!! [NOTE]
      !!! Maybe remove the `lda` from the argument list (since
      !!! it is only used to check against ppm_buffer_dim
      !!! anyway) and check size(pdata,1) against
      !!! ppm_buffer_dim. _(does the size(,) stuff
      !!! also work if the array was allocated in a C client?)_
      !!!
      !!! [WARNING]
      !!! DIM is the dimension of the pdata array and not the space
      !!! dimension ppm_dim!
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
      !!! Particle data to be popped.
      !!! Can be either 1D or 2D array.
#if   __DIM == 2
      INTEGER                 , INTENT(IN   ) :: lda
      !!! The leading dimension of pdata.
#endif
      INTEGER                 , INTENT(IN   ) :: Npart
      !!! The old number of particles (on the processor)
      INTEGER                 , INTENT(IN   ) :: Mpart
      !!! The new number of particles (on the processor)
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3) :: ldu
      INTEGER               :: k,ipart,bdim,ibuffer,btype
      INTEGER               :: iopt,edim,istart
#if   __DIM == 1
      INTEGER, PARAMETER    :: lda = 1
#endif
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
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      ! skip if the buffer is empty
      IF (ppm_buffer_set .LT. 1) THEN
        IF (ppm_debug .GT. 1) THEN
            CALL ppm_write(ppm_rank,'ppm_map_part_pop',  &
     &          'Buffer is empty: skipping pop!',info)
        ENDIF
        GOTO 9999
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
      CONTAINS
      SUBROUTINE check
          IF (ppm_buffer_set .LT. 1) THEN
              info = ppm_error_notice
              CALL ppm_error(ppm_err_argument,'ppm_map_part_pop',  &
     &            'buffer is empty. Cannot pop.',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_pop',  &
     &            'Npart must be >=0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (Mpart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_pop',  &
     &            'Mpart must be >=0',__LINE__,info)
              GOTO 8888
          ENDIF
#if   __DIM == 2
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_pop',  &
     &            'lda must be >0 for vector data',__LINE__,info)
              GOTO 8888
          ENDIF
#elif __DIM == 1
          IF (lda .NE. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_pop',  &
     &            'lda must be =1 for scalar data',__LINE__,info)
              GOTO 8888
          ENDIF
#endif
 8888     CONTINUE
      END SUBROUTINE check

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

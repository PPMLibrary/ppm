      !-------------------------------------------------------------------------
      !  Subroutine   :                 map_part_pop
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
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


#if    __DIM ==1
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE DTYPE(map_part_pop_1d)(Pc,mapID,pdata,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE DTYPE(map_part_pop_1d)(Pc,mapID,pdata,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE DTYPE(map_part_pop_1dc)(Pc,mapID,pdata,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE DTYPE(map_part_pop_1dc)(Pc,mapID,pdata,info)
#elif  __KIND == __INTEGER
      SUBROUTINE DTYPE(map_part_pop_1di)(Pc,mapID,pdata,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE DTYPE(map_part_pop_1dl)(Pc,mapID,pdata,info)
#endif 

#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE DTYPE(map_part_pop_2d)(Pc,mapID,pdata,lda,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE DTYPE(map_part_pop_2d)(Pc,mapID,pdata,lda,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE DTYPE(map_part_pop_2dc)(Pc,mapID,pdata,lda,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE DTYPE(map_part_pop_2dc)(Pc,mapID,pdata,lda,info)
#elif  __KIND == __INTEGER
      SUBROUTINE DTYPE(map_part_pop_2di)(Pc,mapID,pdata,lda,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE DTYPE(map_part_pop_2dl)(Pc,mapID,pdata,lda,info)
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

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      IMPLICIT NONE
      DEFINE_MK()
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      CLASS(DTYPE(ppm_t_particles))        :: Pc
      INTEGER                              :: mapID
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
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(2) :: ldu
      INTEGER               :: k,ipart,bdim,ibuffer,btype
      INTEGER               :: iopt,edim,istart
      INTEGER               :: newNpart,oldNpart
#if   __DIM == 1
      INTEGER, PARAMETER    :: lda = 1
#endif
      CHARACTER(ppm_char)   :: mesg
      CHARACTER(ppm_char)   :: caller = 'map_part_pop'
      REAL(MK)              :: t0
      TYPE(DTYPE(ppm_t_part_mapping)), POINTER :: map => NULL()
      REAL(MK),DIMENSION(:),POINTER :: ppm_recvbuffer => NULL()
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)


      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      map => Pc%maps%vec(mapID)%t
      ppm_recvbuffer => map%ppm_recvbuffer

      ! skip if the buffer is empty
      IF (map%ppm_buffer_set .LT. 1) THEN
        info = ppm_error_notice
        IF (ppm_debug .GT. 1) THEN
            CALL ppm_error(ppm_err_buffer_empt,caller,    &
     &          'Buffer is empty: skipping pop!',__LINE__,info)
        ENDIF
        GOTO 9999
      ENDIF


      !-------------------------------------------------------------------------
      !  Check that the required dimension fits the dimension of the buffer
      !-------------------------------------------------------------------------
      bdim = map%ppm_buffer_dim(map%ppm_buffer_set)
#if   __KIND == __SINGLE_PRECISION_COMPLEX | \
      __KIND == __DOUBLE_PRECISION_COMPLEX
      ! for complex, the effective dimension is half the data dimension
      edim = bdim/2
#else
      edim = bdim
#endif
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(2(A,I3))') 'bdim=',edim,'    lda=',lda
          CALL ppm_write(ppm_rank,caller,mesg,info)
      ENDIF
#if   __DIM == 2
      IF (edim.NE.lda) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_dim,caller,    &
     &       'leading dimension LDA is in error',__LINE__,info)
         GOTO 9999
      ENDIF 
#elif __DIM == 1
      IF (edim.NE.1) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_dim,caller,    &
     &       'buffer does not contain 1d data!',__LINE__,info)
         GOTO 9999
      ENDIF 
#endif

      !-------------------------------------------------------------------------
      !  Check that the required type is identical to the type of the buffer
      !-------------------------------------------------------------------------
      btype = map%ppm_buffer_type(map%ppm_buffer_set)
#if    __KIND == __SINGLE_PRECISION
      IF (btype.NE.ppm_kind_single) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,caller,    &
     &       'trying to pop a non-single into single ',__LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __DOUBLE_PRECISION
      IF (btype.NE.ppm_kind_double) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,caller,    &
     &       'trying to pop a non-double into double ',__LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      IF (btype.NE.ppm_kind_single) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,caller,    &
     &       'trying to pop a non-single-complex into single-complex',  &
     &       __LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      IF (btype.NE.ppm_kind_double) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,caller,    &
     &       'trying to pop a non-double-complex into double-complex',  &
     &       __LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __INTEGER
      IF (btype.NE.ppm_integer) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,caller,    &
     &       'trying to pop a non-integer into integer ',__LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __LOGICAL
      IF (btype.NE.ppm_logical) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,caller,    &
     &       'trying to pop a non-logical into logical ',__LINE__,info)
      ENDIF
#endif
      !-------------------------------------------------------------------------
      !  (Re)allocate the particle data (if necessary)
      !-------------------------------------------------------------------------
      newNpart = map%newNpart
      oldNpart = map%oldNpart
      iopt   = ppm_param_alloc_grow_preserve
#if   __DIM == 2 
      ldu(1) = edim
      ldu(2) = newNpart
#elif __DIM == 1
      ldu(1) = newNpart
#endif
      CALL ppm_alloc(pdata,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,caller,     &
     &        'particle data PDATA',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Decrement the pointer into the receive buffer
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(2(A,I9))') 'ppm_nrecvbuffer = ',map%ppm_nrecvbuffer,   &
     &        'newNpart*bdim = ',newNpart*bdim
          CALL ppm_write(ppm_rank,caller,mesg,info)
      ENDIF
      IF (map%ppm_map_type.EQ.ppm_param_map_ghost_get) THEN
         map%ppm_nrecvbuffer = map%ppm_nrecvbuffer - (newNpart - oldNpart)*bdim 
      ELSE
         map%ppm_nrecvbuffer = map%ppm_nrecvbuffer - newNpart*bdim 
      ENDIF 

      !-------------------------------------------------------------------------
      !  Decrement the pointer into the send buffer to allow reuse by
      !  multiple sequential push-send-pop cycles.
      !-------------------------------------------------------------------------
      map%ppm_nsendbuffer = map%ppm_nsendbuffer - &
          map%ppm_buffer_dim(map%ppm_buffer_set)*  &
     &    (map%ppm_psendbuffer(map%ppm_nsendlist+1)-1)

      ibuffer = map%ppm_nrecvbuffer 

      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(A,I9)') 'ibuffer = ',ibuffer
          CALL ppm_write(ppm_rank,caller,mesg,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  compute the start of the loop: when ghosts are popped, we need to add
      !  them at the end of the current particle list, otherwise we overwrite  
      !  the current particle list
      !-------------------------------------------------------------------------
      IF (map%ppm_map_type.EQ.ppm_param_map_ghost_get) THEN
         istart = oldNpart + 1
      ELSE
         istart = 1
      ENDIF 
      !-------------------------------------------------------------------------
      !  loop over the processors in the ppm_isendlist() 
      !-------------------------------------------------------------------------
#ifdef __MPI
      !-------------------------------------------------------------------------
      !  DOUBLE PRECISION BUFFER
      !-------------------------------------------------------------------------
      IF (ppm_kind.EQ.ppm_kind_double) THEN
#if    __DIM == 2
         !----------------------------------------------------------------------
         !  Unrolled version of edim=1
         !----------------------------------------------------------------------
         IF (edim .EQ. 1) THEN
            DO ipart=istart,newNpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbuffer(ibuffer) .GT.     &
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
            DO ipart=istart,newNpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(2,ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbuffer(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.     &
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
            DO ipart=istart,newNpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(2,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(3,ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbuffer(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(2,ipart) = .TRUE.
               ELSE
                  pdata(2,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.     &
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
            DO ipart=istart,newNpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(4,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(4,ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(4,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(4,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(2,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(3,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(4,ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbuffer(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(2,ipart) = .TRUE.
               ELSE
                  pdata(2,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(3,ipart) = .TRUE.
               ELSE
                  pdata(3,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.     &
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
            DO ipart=istart,newNpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(4,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 1
               pdata(5,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(4,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(5,ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(4,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(5,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(4,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(5,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(2,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(3,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(4,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(5,ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbuffer(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(2,ipart) = .TRUE.
               ELSE
                  pdata(2,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(3,ipart) = .TRUE.
               ELSE
                  pdata(3,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.     &
     &            (1.0_ppm_kind_double-ppm_myepsd)) THEN
                  pdata(4,ipart) = .TRUE.
               ELSE
                  pdata(4,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.     &
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
            DO ipart=istart,newNpart
               DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
                  ibuffer = ibuffer + 2
#else
                  ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                  pdata(k,ipart) = REAL(ppm_recvbuffer(ibuffer),       &
     &                ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                  pdata(k,ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  pdata(k,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &                ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  pdata(k,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &                ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                  pdata(k,ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
                  IF (ppm_recvbuffer(ibuffer) .GT.     &
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
         DO ipart=istart,newNpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
            ibuffer = ibuffer + 2
#else
            ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
            pdata(ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
            pdata(ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
            pdata(ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &          ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
            pdata(ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &          ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
            pdata(ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
            IF (ppm_recvbuffer(ibuffer) .GT.     &
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
            DO ipart=istart,newNpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbuffer(ibuffer) .GT.      &
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
            DO ipart=istart,newNpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(1,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(2,ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbuffer(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.      &
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
            DO ipart=istart,newNpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(2,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(3,ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbuffer(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(2,ipart) = .TRUE.
               ELSE
                  pdata(2,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.      &
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
            DO ipart=istart,newNpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(4,ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(4,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(4,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(4,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(2,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(3,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(4,ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbuffer(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(2,ipart) = .TRUE.
               ELSE
                  pdata(2,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(3,ipart) = .TRUE.
               ELSE
                  pdata(3,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.      &
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
            DO ipart=istart,newNpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(1,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(4,ipart) = ppm_recvbuffer(ibuffer)
               ibuffer = ibuffer + 1
               pdata(5,ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(1,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(2,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(3,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(4,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 1
               pdata(5,ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(4,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
               ibuffer = ibuffer + 2
               pdata(5,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(1,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(2,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(3,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(4,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
               ibuffer = ibuffer + 2
               pdata(5,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &             ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
               pdata(1,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(2,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(3,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(4,ipart) = INT(ppm_recvbuffer(ibuffer))
               ibuffer = ibuffer + 1
               pdata(5,ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbuffer(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(1,ipart) = .TRUE.
               ELSE
                  pdata(1,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(2,ipart) = .TRUE.
               ELSE
                  pdata(2,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(3,ipart) = .TRUE.
               ELSE
                  pdata(3,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.      &
     &            (1.0_ppm_kind_single-ppm_myepss)) THEN
                  pdata(4,ipart) = .TRUE.
               ELSE
                  pdata(4,ipart) = .FALSE.
               ENDIF 
               ibuffer = ibuffer + 1
               IF (ppm_recvbuffer(ibuffer) .GT.      &
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
            DO ipart=istart,newNpart
               DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
                  ibuffer = ibuffer + 2
#else
                  ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                  pdata(k,ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
                  pdata(k,ipart) = REAL(ppm_recvbuffer(ibuffer),       &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  pdata(k,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &                ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  pdata(k,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &                ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
                  pdata(k,ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
                  IF (ppm_recvbuffer(ibuffer) .GT.      &
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
         DO ipart=istart,newNpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
            ibuffer = ibuffer + 2
#else
            ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
            pdata(ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
            pdata(ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
            pdata(ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &          ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
            pdata(ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &          ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
            pdata(ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
            IF (ppm_recvbuffer(ibuffer) .GT.      &
     &         (1.0_ppm_kind_single-ppm_myepss)) THEN
               pdata(ipart) = .TRUE.
            ELSE
               pdata(ipart) = .FALSE.
            ENDIF 
#endif
         ENDDO
#endif
      ENDIF       ! ppm_kind
      ! finish MPI
#else
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  DOUBLE PRECISION BUFFER
      !-------------------------------------------------------------------------
      IF (ppm_kind.EQ.ppm_kind_double) THEN
#if    __DIM == 2
        DO ipart=istart,newNpart
            DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
                ibuffer = ibuffer + 2
#else
                ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                pdata(k,ipart) = REAL(ppm_recvbuffer(ibuffer),       &
     &                ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                pdata(k,ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                pdata(k,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &                ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                pdata(k,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &                ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                pdata(k,ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
                IF (ppm_recvbuffer(ibuffer) .GT.     &
     &               (1.0_ppm_kind_double-ppm_myepsd)) THEN
                    pdata(k,ipart) = .TRUE.
                ELSE
                    pdata(k,ipart) = .FALSE.
                ENDIF 
#endif
               ENDDO
            ENDDO
#elif  __DIM == 1
         !----------------------------------------------------------------------
         !  Scalar version
         !----------------------------------------------------------------------
         DO ipart=istart,newNpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
            ibuffer = ibuffer + 2
#else
            ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
            pdata(ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
            pdata(ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
            pdata(ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &          ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
            pdata(ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &          ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
            pdata(ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
            IF (ppm_recvbuffer(ibuffer) .GT.     &
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
        DO ipart=istart,newNpart
            DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
                ibuffer = ibuffer + 2
#else
                ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                pdata(k,ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
                pdata(k,ipart) = REAL(ppm_recvbuffer(ibuffer),       &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                pdata(k,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &                ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                pdata(k,ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &                ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
                pdata(k,ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
                IF (ppm_recvbuffer(ibuffer) .GT.      &
     &               (1.0_ppm_kind_single-ppm_myepss)) THEN
                     pdata(k,ipart) = .TRUE.
                ELSE
                     pdata(k,ipart) = .FALSE.
                ENDIF 
#endif
            ENDDO
        ENDDO
#elif  __DIM == 1
         !----------------------------------------------------------------------
         !  Scalar version
         !----------------------------------------------------------------------
         DO ipart=istart,newNpart
#if    __KIND == __SINGLE_PRECISION_COMPLEX | \
       __KIND == __DOUBLE_PRECISION_COMPLEX
            ibuffer = ibuffer + 2
#else
            ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
            pdata(ipart) = ppm_recvbuffer(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
            pdata(ipart) = REAL(ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
            pdata(ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &          ppm_recvbuffer(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
            pdata(ipart) = CMPLX(ppm_recvbuffer(ibuffer-1),    &
     &          ppm_recvbuffer(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
            pdata(ipart) = INT(ppm_recvbuffer(ibuffer))
#elif  __KIND == __LOGICAL
            IF (ppm_recvbuffer(ibuffer) .GT.      &
     &         (1.0_ppm_kind_single-ppm_myepss)) THEN
               pdata(ipart) = .TRUE.
            ELSE
               pdata(ipart) = .FALSE.
            ENDIF 
#endif
         ENDDO
#endif
      ENDIF       ! ppm_kind


      ! finish non-MPI
#endif 

      !-------------------------------------------------------------------------
      !  Decrement the set counter
      !-------------------------------------------------------------------------
      map%ppm_buffer_set = map%ppm_buffer_set - 1  

      !-------------------------------------------------------------------------
      !  Deallocate the receive buffer if all sets have been poped
      !-------------------------------------------------------------------------
      IF (map%ppm_buffer_set .LT. 1) THEN
          iopt = ppm_param_dealloc
          CALL ppm_alloc(map%ppm_recvbuffer,ldu,iopt,info)
          ppm_recvbuffer => map%ppm_recvbuffer
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_alloc,caller,     &
     &            'receive buffer PPM_RECVBUFFER',__LINE__,info)
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. Pc%maps%exists(mapID)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_wrong_dim,caller,    &
                  &   'Invalid mapID: mapping does not exist.',__LINE__,info)
              GOTO 8888
          ENDIF 
          IF (Pc%maps%vec(mapID)%t%oldNpart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
     &            'Npart must be >=0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (Pc%maps%vec(mapID)%t%newNpart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
     &            'newNpart must be >=0',__LINE__,info)
              GOTO 8888
          ENDIF
#if   __DIM == 2
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
     &            'lda must be >0 for vector data',__LINE__,info)
              GOTO 8888
          ENDIF
#elif __DIM == 1
          IF (lda .NE. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
     &            'lda must be =1 for scalar data',__LINE__,info)
              GOTO 8888
          ENDIF
#endif
 8888     CONTINUE
      END SUBROUTINE check

#if    __DIM == 1
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE DTYPE(map_part_pop_1d)
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE DTYPE(map_part_pop_1d)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE DTYPE(map_part_pop_1dc)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE DTYPE(map_part_pop_1dc)
#elif  __KIND == __INTEGER
      END SUBROUTINE DTYPE(map_part_pop_1di)
#elif  __KIND == __LOGICAL
      END SUBROUTINE DTYPE(map_part_pop_1dl)
#endif

#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE DTYPE(map_part_pop_2d)
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE DTYPE(map_part_pop_2d)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE DTYPE(map_part_pop_2dc)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE DTYPE(map_part_pop_2dc)
#elif  __KIND == __INTEGER
      END SUBROUTINE DTYPE(map_part_pop_2di)
#elif  __KIND == __LOGICAL
      END SUBROUTINE DTYPE(map_part_pop_2dl)
#endif
#endif

#undef __KIND

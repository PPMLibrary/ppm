      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_map_part_pop
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


#if    __VARIANT == __NORMAL
#if    __DIM ==1
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_pop_1ds(pdata,Npart,newNpart,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_pop_1dd(pdata,Npart,newNpart,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_pop_1dsc(pdata,Npart,newNpart,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_pop_1ddc(pdata,Npart,newNpart,info)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_part_pop_1di(pdata,Npart,newNpart,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_part_pop_1dl(pdata,Npart,newNpart,info)
#endif 

#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_pop_2ds(pdata,lda,Npart,newNpart,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_pop_2dd(pdata,lda,Npart,newNpart,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_pop_2dsc(pdata,lda,Npart,newNpart,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_pop_2ddc(pdata,lda,Npart,newNpart,info)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_part_pop_2di(pdata,lda,Npart,newNpart,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_part_pop_2dl(pdata,lda,Npart,newNpart,info)
#endif
#endif
#elif  __VARIANT == __ADD
#if    __DIM ==1
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_pop_add1ds(pdata_old,pdata_add,Npart,Mpart,&
              newNpart,newMpart,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_pop_add1dd(pdata_old,pdata_add,Npart,Mpart,&
              newNpart,newMpart,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_pop_add1dsc(pdata_old,pdata_add,Npart,Mpart,&
              newNpart,newMpart,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_pop_add1ddc(pdata_old,pdata_add,Npart,Mpart,&
              newNpart,newMpart,info)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_part_pop_add1di(pdata_old,pdata_add,Npart,Mpart,&
              newNpart,newMpart,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_part_pop_add1dl(pdata_old,pdata_add,Npart,Mpart,&
              newNpart,newMpart,info)
#endif 

#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_pop_add2ds(pdata_old,pdata_add,lda,Npart,Mpart,&
              newNpart,newMpart,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_pop_add2dd(pdata_old,pdata_add,lda,Npart,Mpart,&
              newNpart,newMpart,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_pop_add2dsc(pdata_old,pdata_add,lda,Npart,Mpart,&
              newNpart,newMpart,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_pop_add2ddc(pdata_old,pdata_add,lda,Npart,Mpart,&
              newNpart,newMpart,info)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_part_pop_add2di(pdata_old,pdata_add,lda,Npart,Mpart,&
              newNpart,newMpart,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_part_pop_add2dl(pdata_old,pdata_add,lda,Npart,Mpart,&
              newNpart,newMpart,info)
#endif
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
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
#if    __VARIANT == __NORMAL
      USE ppm_module_data_buffers
#elif  __VARIANT == __ADD
      USE ppm_module_data_buffers_add
#endif
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if   __VARIANT == __NORMAL
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
#elif __VARIANT == __ADD
#if   __DIM == 1
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:  ), POINTER    :: pdata_old => NULL()
      INTEGER , DIMENSION(:  ), POINTER,INTENT(INOUT)    :: pdata_add => NULL()
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:  ), POINTER    :: pdata_old => NULL()
      LOGICAL , DIMENSION(:  ), POINTER,INTENT(INOUT)    :: pdata_add => NULL()
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:  ), POINTER :: pdata_old => NULL()
      COMPLEX(MK), DIMENSION(:  ), POINTER,INTENT(INOUT) :: pdata_add => NULL()
#else
      REAL(MK), DIMENSION(:  ), POINTER    :: pdata_old => NULL()
      REAL(MK), DIMENSION(:  ), POINTER,INTENT(INOUT)    :: pdata_add => NULL()
#endif

#elif __DIM == 2
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:,:), POINTER    :: pdata_old => NULL()
      INTEGER , DIMENSION(:,:), POINTER,INTENT(INOUT)    :: pdata_add => NULL()
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:), POINTER    :: pdata_old => NULL()
      LOGICAL , DIMENSION(:,:), POINTER,INTENT(INOUT)    :: pdata_add => NULL()
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:), POINTER :: pdata_old => NULL()
      COMPLEX(MK), DIMENSION(:,:), POINTER,INTENT(INOUT) :: pdata_add => NULL()
#else
      REAL(MK), DIMENSION(:,:), POINTER    :: pdata_old => NULL()
      REAL(MK), DIMENSION(:,:), POINTER,INTENT(INOUT)    :: pdata_add => NULL()
#endif
#endif
#endif
      
      !!! Particle data to be popped and to be added, respectively. 
      !!! Can be either 1D or 2D array.
#if   __DIM == 2
      INTEGER                 , INTENT(IN   ) :: lda
      !!! The leading dimension of pdata.
#endif
      INTEGER                 , INTENT(IN   ) :: Npart
      !!! The old number of particles (on the processor)
#if   __VARIANT == __NORMAL
      INTEGER                 , INTENT(IN   ) :: newNpart
      !!! The new number of particles (on the processor,incl ghosts)
#elif __VARIANT == __ADD
      INTEGER                 , INTENT(IN   ) :: Mpart
      !!! The old number of particles (on the processor,incl ghosts)
      INTEGER                 , INTENT(IN   ) :: newNpart
      !!! The new number of real particles (on the processor,excl. ghosts)
      INTEGER                 , INTENT(IN   ) :: newMpart
      !!! The new number of particles (on the processor,incl ghosts)
#endif
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(2) :: ldu
      INTEGER               :: i,j,ipart_add
      INTEGER               :: k,ipart,bdim,ibuffer,btype
      INTEGER               :: iopt,edim,istart,iend
#if   __DIM == 1
      INTEGER, PARAMETER    :: lda = 1
#endif
      CHARACTER(ppm_char)   :: mesg
#if   __VARIANT == __NORMAL
      CHARACTER(ppm_char)   :: caller = 'ppm_map_part_pop'
#elif __VARIANT == __ADD
      CHARACTER(ppm_char)   :: caller = 'ppm_part_modify_pop'
#endif
      REAL(MK)              :: t0
#if   __VARIANT == __ADD
#if   __DIM == 1
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:  ), POINTER    :: pdata => NULL()
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:  ), POINTER    :: pdata => NULL()
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:  ), POINTER :: pdata => NULL()
#else
      REAL(MK), DIMENSION(:  ), POINTER    :: pdata => NULL()
#endif

#elif __DIM == 2
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:,:), POINTER    :: pdata => NULL()
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:), POINTER    :: pdata => NULL()
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:), POINTER :: pdata => NULL()
#else
      REAL(MK), DIMENSION(:,:), POINTER    :: pdata => NULL()
#endif
#endif
#endif
#if   __VARIANT == __ADD
      INTEGER                              :: npart_added
#endif
      
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

      ! skip if the buffer is empty
      IF (ppm_buffer_set .LT. 1) THEN
        info = ppm_error_notice
        IF (ppm_debug .GT. 1) THEN
            CALL ppm_error(ppm_err_buffer_empt,caller,    &
     &          'Buffer is empty: skipping pop!',__LINE__,info)
        ENDIF
        GOTO 9999
      ENDIF

#if   __VARIANT == __ADD
      !number of particles to receive from the buffer
      IF (add_mode .EQ. ppm_param_add_ghost_particles) THEN
          npart_added = newMpart - Mpart - (newNpart - Npart) 
      ELSE
          npart_added = newMpart - Mpart - (newNpart - Npart) 
      ENDIF
#endif

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
      btype = ppm_buffer_type(ppm_buffer_set)
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
#if   __VARIANT == __NORMAL
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
#elif __VARIANT == __ADD
      iopt   = ppm_param_alloc_grow_preserve
      IF (add_mode .EQ. ppm_param_add_ghost_particles) THEN
#if   __DIM == 2 
          ldu(1) = edim
          ldu(2) = newMpart - Mpart
#elif __DIM == 1
          ldu(1) = newMpart - Mpart
#endif
          CALL ppm_alloc(pdata_add,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,caller,     &
                  &        'particle data PDATA_ADD',__LINE__,info)
              GOTO 9999
          ENDIF
      ELSE
#if   __DIM == 2 
          ldu(1) = edim
          ldu(2) = newMpart
#elif __DIM == 1
          ldu(1) = newMpart
#endif
          CALL ppm_alloc(pdata_old,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,caller,     &
                  &        'particle data PDATA_OLD',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      iopt   = ppm_param_alloc_grow
#if   __DIM == 2 
      ldu(1) = edim
      ldu(2) = npart_added
#elif __DIM == 1
      ldu(1) = npart_added
#endif
      CALL ppm_alloc(pdata,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,caller,     &
     &        'particle data PDATA',__LINE__,info)
          GOTO 9999
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Decrement the pointer into the receive buffer
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(2(A,I9))') 'ppm_nrecvbuffer = ',ppm_nrecvbuffer,   &
     &        ' newNpart*bdim = ',newNpart*bdim
          CALL ppm_write(ppm_rank,caller,mesg,info)
      ENDIF
#if   __VARIANT == __NORMAL
      IF (ppm_map_type.EQ.ppm_param_map_ghost_get) THEN
         ppm_nrecvbuffer = ppm_nrecvbuffer - (newNpart - Npart)*bdim 
      ELSE
         ppm_nrecvbuffer = ppm_nrecvbuffer - newNpart*bdim 
      ENDIF 
#elif __VARIANT == __ADD
      ppm_nrecvbuffer = ppm_nrecvbuffer - npart_added*bdim 
#endif

      !-------------------------------------------------------------------------
      !  Decrement the pointer into the send buffer to allow reuse by
      !  multiple sequential push-send-pop cycles.
      !-------------------------------------------------------------------------
      ppm_nsendbuffer = ppm_nsendbuffer - ppm_buffer_dim(ppm_buffer_set)*  &
     &    (ppm_psendbuffer(ppm_nsendlist+1)-1)

      ibuffer = ppm_nrecvbuffer 

      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(A,I9)') 'ibuffer = ',ibuffer
          CALL ppm_write(ppm_rank,caller,mesg,info)
          !FIXME remove next 2 lines
          WRITE(mesg,*) 'ppm_psendbuffer = ',ppm_psendbuffer
          CALL ppm_write(ppm_rank,caller,mesg,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  compute the start of the loop: when ghosts are popped, we need to add
      !  them at the end of the current particle list, otherwise we overwrite  
      !  the current particle list
      !-------------------------------------------------------------------------
#if   __VARIANT == __NORMAL
      IF (ppm_map_type.EQ.ppm_param_map_ghost_get) THEN
         istart = Npart + 1
      ELSE
         istart = 1
      ENDIF 
      iend = newNpart
#elif __VARIANT == __ADD
     istart = 1
     iend = npart_added !nb of particles in buffer
#endif
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
            DO ipart=istart,iend
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
            DO ipart=istart,iend
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
            DO ipart=istart,iend
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
            DO ipart=istart,iend
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
            DO ipart=istart,iend
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
            DO ipart=istart,iend
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
         DO ipart=istart,iend
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
            DO ipart=istart,iend
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
            DO ipart=istart,iend
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
            DO ipart=istart,iend
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
            DO ipart=istart,iend
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
            DO ipart=istart,iend
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
            DO ipart=istart,iend
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
         DO ipart=istart,iend
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
      ! finish MPI
#else
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  DOUBLE PRECISION BUFFER
      !-------------------------------------------------------------------------
      IF (ppm_kind.EQ.ppm_kind_double) THEN
#if    __DIM == 2
        DO ipart=istart,iend
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
#elif  __DIM == 1
         !----------------------------------------------------------------------
         !  Scalar version
         !----------------------------------------------------------------------
         DO ipart=istart,iend
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
        DO ipart=istart,iend
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
#elif  __DIM == 1
         !----------------------------------------------------------------------
         !  Scalar version
         !----------------------------------------------------------------------
         DO ipart=istart,iend
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


      ! finish non-MPI
#endif 

#if    __VARIANT == __ADD
      !-------------------------------------------------------------------------
      !  Weave the received data inside the existing one (so that the indexing 
      !  is preserved)
      !-------------------------------------------------------------------------
      IF (add_mode .EQ. ppm_param_add_ghost_particles) THEN
          !-------------------------------------------------------------------------
          !  The new data are real values
          !  we just have to append the data from the buffer at the end
          !  of pdata_add 
          !-------------------------------------------------------------------------
#if    __DIM == 1
          pdata_add(newNpart-Npart+1:newNpart-Npart+npart_added) = &
              pdata(1:npart_added)
#elif  __DIM == 2
          pdata_add(1:lda,newNpart-Npart+1:newNpart-Npart+npart_added) = &
              pdata(1:lda,1:npart_added)
#endif
      ELSE
          !-------------------------------------------------------------------------
          !  The new data are ghost values
          !  The data from the buffer go between newNpart+1 and newMpart. 
          !  This memory area consists of separate chuncks, each of which
          !  corresponds to a given neighboring processor. The new data has to
          !  be inserted so as to preserve this ordering.
          !  The new real data is inserted between Npart+1 and newNpart.
          !-------------------------------------------------------------------------
          ipart = Mpart + (newNpart-Npart)
          ipart_add = npart_added
          DO k=ppm_nsendlist,1,-1
              i=plists_normal%precv(k)
              j=plists_add%precv(k)
#if    __DIM == 1
              pdata_old(ipart+ipart_add-i-j+1:ipart+ipart_add-j) = &
                  pdata_old(ipart-i+1:ipart)
              pdata_old(ipart+ipart_add-j+1:ipart+ipart_add) = &
                  pdata(ipart_add-j+1:ipart_add)
#elif  __DIM == 2
              pdata_old(1:lda,ipart+ipart_add-i-j+1:ipart+ipart_add-j) = &
                  pdata_old(1:lda,ipart-i+1:ipart)
              pdata_old(1:lda,ipart+ipart_add-j+1:ipart+ipart_add) = &
                  pdata(1:lda,ipart_add-j+1:ipart_add)
#endif
              ipart = ipart - i
              ipart_add = ipart_add - j
          ENDDO ! loop over all processors in commseq
          DO i=1,newNpart-Npart
#if    __DIM == 1
          pdata_old(Npart+i) = pdata_add(modify%idx_real_new(i))
#elif  __DIM == 2
          pdata_old(1:lda,Npart+i) = pdata_add(1:lda,modify%idx_real_new(i))
#endif
          ENDDO
      ENDIF
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
          IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
     &            'Npart must be >=0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (newNpart .LT. 0) THEN
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

#if    __VARIANT == __NORMAL
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
#elif  __VARIANT == __ADD
#if    __DIM == 1
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_pop_add1ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_pop_add1dd
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_pop_add1dsc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_pop_add1ddc
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_map_part_pop_add1di
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_map_part_pop_add1dl
#endif

#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_pop_add2ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_pop_add2dd
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_pop_add2dsc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_pop_add2ddc
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_map_part_pop_add2di
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_map_part_pop_add2dl
#endif
#endif
#endif

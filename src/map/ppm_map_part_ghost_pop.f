      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_map_part_ghost_pop
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

#if    __DIM ==1
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_ghost_pop_1ds(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_ghost_pop_1dd(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_ghost_pop_1dsc(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_ghost_pop_1ddc(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_part_ghost_pop_1di(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_part_ghost_pop_1dl(pdata,lda,Npart,Mpart,info)
#endif

#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_ghost_pop_2ds(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_ghost_pop_2dd(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_ghost_pop_2dsc(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_ghost_pop_2ddc(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_part_ghost_pop_2di(pdata,lda,Npart,Mpart,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_part_ghost_pop_2dl(pdata,lda,Npart,Mpart,info)
#endif
#endif
      !!! This routine pops the contents of the receive buffer
      !!! when sending the ghosts back. Notice the value popped
      !!! will be *added* to array passed to the routine.
      !!!
      !!! [WARNING]
      !!! `lda` is the dimension of the pdata array and not the space
      !!! dimension ppm_dim!
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
      !  Includes
      !-------------------------------------------------------------------------
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
      INTEGER                 , INTENT(IN   ) :: lda
      !!! The leading dimension of pdata.
      !!! In the case of pdata being a 1D array use lda=1.
      INTEGER                 , INTENT(IN   ) :: Npart
      !!! The number of real particles (on the processor)
      INTEGER                 , INTENT(IN   ) :: Mpart
      !!! The number of real + ghost particles (on the processor)
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status. 0 upon success.
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK) :: t0

      INTEGER, DIMENSION(3) :: ldu
      INTEGER               :: k,ipart,bdim,ibuffer,btype
      INTEGER               :: iopt,edim,i,j

      CHARACTER(ppm_char) :: caller='ppm_map_part_ghost_pop'
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
      IF (ppm_debug.GT.0) THEN
         CALL check
         IF (info.NE.0) GOTO 9999
      ENDIF

      ! skip if the buffer is empty
      IF (ppm_buffer_set.LT.1) THEN
         IF (ppm_debug.GT.1) THEN
            fail('Buffer is empty: skipping pop!', &
            & ppm_err_buffer_empt,ppm_error=ppm_error_notice)
         ELSE
            info=ppm_error_notice
            GOTO 9999
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that the required dimension fits the dimension of the buffer
      !-------------------------------------------------------------------------
      bdim = ppm_buffer_dim(ppm_buffer_set)

#if   __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      ! for complex, the effective dimension is half the data dimension
      edim = bdim/2
#else
      edim = bdim
#endif
      IF (ppm_debug .GT. 1) THEN
         stdout_f('(2(A,I3))',"bdim=",edim,"    lda=",lda)
      ENDIF
#if   __DIM == 2
      IF (edim.NE.lda) THEN
         fail('leading dimension LDA is in error',ppm_err_wrong_dim)
      ENDIF
#elif __DIM == 1
      IF (edim.NE.1) THEN
         fail('buffer does not contain 1d data!',ppm_err_wrong_dim)
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Check that the required type is identical to the type of the buffer
      !-------------------------------------------------------------------------
      btype = ppm_buffer_type(ppm_buffer_set)

#if    __KIND == __SINGLE_PRECISION
      IF (btype.NE.ppm_kind_single) THEN
         fail('trying to pop a non-single into single ',ppm_err_wrong_prec)
      ENDIF
#elif  __KIND == __DOUBLE_PRECISION
      IF (btype.NE.ppm_kind_double) THEN
         fail('trying to pop a non-double into double ',ppm_err_wrong_prec)
      ENDIF
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      IF (btype.NE.ppm_kind_single) THEN
         fail('trying to pop a non-single-complex into single-complex ',ppm_err_wrong_prec)
      ENDIF
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      IF (btype.NE.ppm_kind_double) THEN
         fail('trying to pop a non-double-complex into double-complex ',ppm_err_wrong_prec)
      ENDIF
#elif  __KIND == __INTEGER
      IF (btype.NE.ppm_integer) THEN
         fail('trying to pop a non-integer into integer ',ppm_err_wrong_prec)
      ENDIF
#elif  __KIND == __LOGICAL
      IF (btype.NE.ppm_logical) THEN
         fail('trying to pop a non-logical into logical ',ppm_err_wrong_prec)
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Decrement the pointer into the receive buffer - this should be the
      !  total amount of particle send off to any processor time the size
      !  (dimension) of the data; see ppm_map_part_ghost_put for the contents
      !  of the ppm_ghosthack array
      !-------------------------------------------------------------------------
      i = 0
      DO k=1,ppm_nsendlist ! or ppm_ncommseq(topoid) or ppm_nrecvlist
         i = i + ppm_ghosthack(1,k)
      ENDDO
      ppm_nrecvbuffer = ppm_nrecvbuffer - i*bdim

      !-------------------------------------------------------------------------
      !  Decrement the pointer into the send buffer to allow reuse by
      !  multiple sequential push-send-pop cycles.
      ! JHW i dont know if this works with ghost - no time to check
      ! but should not matter here
      !-------------------------------------------------------------------------
      ppm_nsendbuffer = ppm_nsendbuffer - &
      & ppm_buffer_dim(ppm_buffer_set)*(ppm_psendbuffer(ppm_nsendlist+1)-1)

      !-------------------------------------------------------------------------
      !  loop over the processors in the ppm_isendlist()
      !-------------------------------------------------------------------------
#ifdef __MPI
      ibuffer = ppm_nrecvbuffer

      IF (ppm_debug .GT. 1) THEN
         stdout_f('(A,I9)',"ibuffer = ",ibuffer)
      ENDIF

      !-------------------------------------------------------------------------
      !  DOUBLE PRECISION BUFFER
      !-------------------------------------------------------------------------
      SELECT CASE (ppm_kind)
      CASE (ppm_kind_double)
#if    __DIM == 2
         !----------------------------------------------------------------------
         !  Unrolled version of edim=1
         !----------------------------------------------------------------------
         SELECT CASE (edim)
         CASE (1)
            !-------------------------------------------------------------------
            !  now empty the receive buffer:
            !-------------------------------------------------------------------
            DO k=1,ppm_nsendlist ! or ppm_nrecvlist or ppm_ncommseq(topoid)
               DO j=1,ppm_ghosthack(1,k)
                  !----------------------------------------------------------------
                  !  Get the id of the particle
                  !----------------------------------------------------------------
! if we are going to run this efficiently on ivos machine,
!we should put a vector directive here !
                  ipart = ppm_ghosthack(j+2,k)

#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                  ibuffer = ibuffer + 2
#else
                  ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + &
                  & CMPLX(ppm_recvbufferd(ibuffer-1),ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + &
                  & CMPLX(ppm_recvbufferd(ibuffer-1),ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                  pdata(1,ipart) = pdata(1,ipart) + INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
! JHW does not make sence to get the logical contrib from ghost or ?
                  IF (ppm_recvbufferd(ibuffer).GT.(1._ppm_kind_double-ppm_myepsd)) THEN
                     pdata(1,ipart) = .TRUE.
                  ELSE
                     pdata(1,ipart) = .FALSE.
                  ENDIF
#endif
               ENDDO
            ENDDO
         !----------------------------------------------------------------------
         !  Unrolled version of edim=2
         !----------------------------------------------------------------------
         CASE (2)
            DO k=1,ppm_nsendlist ! or ppm_nrecvlist or ppm_ncommseq(topoid)
               DO j=1,ppm_ghosthack(1,k)
                  !-------------------------------------------------------------
                  !  Get the id of the particle
                  !-------------------------------------------------------------
! if we are going to run this efficiently on ivos machine, we should put a vector directive here !
                  ipart = ppm_ghosthack(j+2,k)
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                  ibuffer = ibuffer + 2
#else
                  ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + ppm_recvbufferd(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + &
                  & CMPLX(ppm_recvbufferd(ibuffer-1),ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(2,ipart) = pdata(2,ipart) + &
                  & CMPLX(ppm_recvbufferd(ibuffer-1),ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + &
                  & CMPLX(ppm_recvbufferd(ibuffer-1),ppm_recvbufferd(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(2,ipart) = pdata(2,ipart) + &
                  & CMPLX(ppm_recvbufferd(ibuffer-1),ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                  pdata(1,ipart) = pdata(1,ipart) + INT(ppm_recvbufferd(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                  IF (ppm_recvbufferd(ibuffer).GT.(1._ppm_kind_double-ppm_myepsd)) THEN
                     pdata(1,ipart) = .TRUE.
                  ELSE
                     pdata(1,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbufferd(ibuffer).GT.(1._ppm_kind_double-ppm_myepsd)) THEN
                     pdata(2,ipart) = .TRUE.
                  ELSE
                     pdata(2,ipart) = .FALSE.
                  ENDIF
#endif
               ENDDO
            ENDDO
         !----------------------------------------------------------------------
         !  Unrolled version of edim=3
         !----------------------------------------------------------------------
         CASE (3)
            DO k=1,ppm_nsendlist ! or ppm_nrecvlist or ppm_ncommseq(topoid)
               DO j=1,ppm_ghosthack(1,k)
                  !-------------------------------------------------------------
                  !  Get the id of the particle
                  !-------------------------------------------------------------
! if we are going to run this efficiently on ivos machine, we should put a vector directive here !
                  ipart = ppm_ghosthack(j+2,k)
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                  ibuffer = ibuffer + 2
#else
                  ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
!JHWWRITE(40+ppm_rank,'(I8,6E12.4)') ipart,pdata(:,ipart),ppm_recvbufferd(ibuffer),ppm_recvbufferd(ibuffer+1),ppm_recvbufferd(ibuffer+2)
                  pdata(1,ipart) = pdata(1,ipart) + ppm_recvbufferd(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + ppm_recvbufferd(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + ppm_recvbufferd(ibuffer)

#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(2,ipart) = pdata(2,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(3,ipart) = pdata(3,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(2,ipart) = pdata(2,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(3,ipart) = pdata(3,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                  pdata(1,ipart) = pdata(1,ipart) + INT(ppm_recvbufferd(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + INT(ppm_recvbufferd(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                  IF (ppm_recvbufferd(ibuffer) .GT.     &
     &               (1.0_ppm_kind_double-ppm_myepsd)) THEN
                     pdata(1,ipart) = .TRUE.
                  ELSE
                     pdata(1,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbufferd(ibuffer) .GT.     &
     &               (1.0_ppm_kind_double-ppm_myepsd)) THEN
                     pdata(2,ipart) = .TRUE.
                  ELSE
                     pdata(2,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbufferd(ibuffer) .GT.     &
     &               (1.0_ppm_kind_double-ppm_myepsd)) THEN
                     pdata(3,ipart) = .TRUE.
                  ELSE
                     pdata(3,ipart) = .FALSE.
                  ENDIF
#endif
               ENDDO
            ENDDO
         !----------------------------------------------------------------------
         !  Unrolled version of edim=4
         !----------------------------------------------------------------------
         CASE (4)
            DO k=1,ppm_nsendlist ! or ppm_nrecvlist or ppm_ncommseq(topoid)
               DO j=1,ppm_ghosthack(1,k)
                  !-------------------------------------------------------------
                  !  Get the id of the particle
                  !-------------------------------------------------------------
! if we are going to run this efficiently on ivos machine, we should put a vector directive here !
                  ipart = ppm_ghosthack(j+2,k)
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                  ibuffer = ibuffer + 2
#else
                  ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  pdata(4,ipart) = pdata(4,ipart) + REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + ppm_recvbufferd(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + ppm_recvbufferd(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + ppm_recvbufferd(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(4,ipart) = pdata(4,ipart) + ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(2,ipart) = pdata(2,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(3,ipart) = pdata(3,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(4,ipart) = pdata(4,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(2,ipart) = pdata(2,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(3,ipart) = pdata(3,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(4,ipart) = pdata(4,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                  pdata(1,ipart) = pdata(1,ipart) + INT(ppm_recvbufferd(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + INT(ppm_recvbufferd(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + INT(ppm_recvbufferd(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(4,ipart) = pdata(4,ipart) + INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                  IF (ppm_recvbufferd(ibuffer) .GT.     &
     &               (1.0_ppm_kind_double-ppm_myepsd)) THEN
                     pdata(1,ipart) = .TRUE.
                  ELSE
                     pdata(1,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbufferd(ibuffer) .GT.     &
     &               (1.0_ppm_kind_double-ppm_myepsd)) THEN
                     pdata(2,ipart) = .TRUE.
                  ELSE
                     pdata(2,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbufferd(ibuffer) .GT.     &
     &               (1.0_ppm_kind_double-ppm_myepsd)) THEN
                     pdata(3,ipart) = .TRUE.
                  ELSE
                     pdata(3,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbufferd(ibuffer) .GT.     &
     &               (1.0_ppm_kind_double-ppm_myepsd)) THEN
                     pdata(4,ipart) = .TRUE.
                  ELSE
                     pdata(4,ipart) = .FALSE.
                  ENDIF
#endif
               ENDDO
            ENDDO
         !----------------------------------------------------------------------
         !  Unrolled version of edim=5
         !----------------------------------------------------------------------
         CASE (5)
            DO k=1,ppm_nsendlist ! or ppm_nrecvlist or ppm_ncommseq(topoid)
               DO j=1,ppm_ghosthack(1,k)
               !-------------------------------------------------------------
               !  Get the id of the particle
               !-------------------------------------------------------------
! if we are going to run this efficiently on ivos machine, we should put a vector directive here !
                  ipart = ppm_ghosthack(j+2,k)
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                  ibuffer = ibuffer + 2
#else
                  ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  pdata(4,ipart) = pdata(4,ipart) + REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  pdata(5,ipart) = pdata(5,ipart) + REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + ppm_recvbufferd(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + ppm_recvbufferd(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + ppm_recvbufferd(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(4,ipart) = pdata(4,ipart) + ppm_recvbufferd(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(5,ipart) = pdata(5,ipart) + ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(2,ipart) = pdata(2,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(3,ipart) = pdata(3,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(4,ipart) = pdata(4,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(5,ipart) = pdata(5,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(2,ipart) = pdata(2,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(3,ipart) = pdata(3,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(4,ipart) = pdata(4,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(5,ipart) = pdata(5,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                  pdata(1,ipart) = pdata(1,ipart) + INT(ppm_recvbufferd(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + INT(ppm_recvbufferd(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + INT(ppm_recvbufferd(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(4,ipart) = pdata(4,ipart) + INT(ppm_recvbufferd(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(5,ipart) = pdata(5,ipart) + INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                  IF (ppm_recvbufferd(ibuffer) .GT.     &
     &               (1.0_ppm_kind_double-ppm_myepsd)) THEN
                     pdata(1,ipart) = .TRUE.
                  ELSE
                     pdata(1,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbufferd(ibuffer) .GT.     &
     &               (1.0_ppm_kind_double-ppm_myepsd)) THEN
                     pdata(2,ipart) = .TRUE.
                  ELSE
                     pdata(2,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbufferd(ibuffer) .GT.     &
     &               (1.0_ppm_kind_double-ppm_myepsd)) THEN
                     pdata(3,ipart) = .TRUE.
                  ELSE
                     pdata(3,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbufferd(ibuffer) .GT.     &
     &               (1.0_ppm_kind_double-ppm_myepsd)) THEN
                     pdata(4,ipart) = .TRUE.
                  ELSE
                     pdata(4,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbufferd(ibuffer) .GT.     &
     &               (1.0_ppm_kind_double-ppm_myepsd)) THEN
                     pdata(5,ipart) = .TRUE.
                  ELSE
                     pdata(5,ipart) = .FALSE.
                  ENDIF
#endif
               ENDDO
            ENDDO

         !----------------------------------------------------------------------
         !  For edim.GT.5 the vector length will be edim !!
         !----------------------------------------------------------------------
         CASE DEFAULT
            DO k=1,ppm_nsendlist ! or ppm_nrecvlist or ppm_ncommseq(topoid)
               DO j=1,ppm_ghosthack(1,k)
               !-------------------------------------------------------------
               !  Get the id of the particle
               !-------------------------------------------------------------
! if we are going to run this efficiently on ivos machine, we should put a vector directive here !
                  ipart = ppm_ghosthack(j+2,k)
                  DO i=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                     ibuffer = ibuffer + 2
#else
                     ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                     pdata(i,ipart) = pdata(i,ipart) + REAL(ppm_recvbufferd(ibuffer),       &
     &                   ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                     pdata(i,ipart) = pdata(i,ipart) + ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                     pdata(i,ipart) = pdata(i,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                   ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                     pdata(i,ipart) = pdata(i,ipart) + CMPLX(ppm_recvbufferd(ibuffer-1),    &
     &                   ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                     pdata(i,ipart) = pdata(i,ipart) + INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                     IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                  (1.0_ppm_kind_double-ppm_myepsd)) THEN
                        pdata(i,ipart) = .TRUE.
                     ELSE
                        pdata(i,ipart) = .FALSE.
                     ENDIF
#endif
                  ENDDO
               ENDDO
            ENDDO

         END SELECT
#elif  __DIM == 1
         !----------------------------------------------------------------------
         !  Scalar version
         !----------------------------------------------------------------------
         DO k=1,ppm_nsendlist ! or ppm_nrecvlist or ppm_ncommseq(topoid)
            DO j=1,ppm_ghosthack(1,k)
               !-------------------------------------------------------------
               !  Get the id of the particle
               !-------------------------------------------------------------
               ! if we are going to run this efficiently on ivos machine,
               ! we should put a vector directive here !
               ipart = ppm_ghosthack(j+2,k)
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(ipart) = pdata(ipart) + REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(ipart) = pdata(ipart) + ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(ipart) = pdata(ipart) + &
               & CMPLX(ppm_recvbufferd(ibuffer-1),ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(ipart) = pdata(ipart) + &
               & CMPLX(ppm_recvbufferd(ibuffer-1),ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
               pdata(ipart) = pdata(ipart) + INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbufferd(ibuffer).GT.(1._ppm_kind_double-ppm_myepsd)) THEN
                  pdata(ipart) = .TRUE.
               ELSE
                  pdata(ipart) = .FALSE.
               ENDIF
#endif
            ENDDO
         ENDDO
#endif

      !-------------------------------------------------------------------------
      !  SINGLE PRECISION BUFFER
      !-------------------------------------------------------------------------
      CASE DEFAULT
#if    __DIM == 2
         !----------------------------------------------------------------------
         !  Unrolled verion for edim=1
         !----------------------------------------------------------------------
         SELECT CASE (edim)
         CASE (1)
            DO k=1,ppm_nsendlist ! or ppm_nrecvlist or ppm_ncommseq(topoid)
               DO j=1,ppm_ghosthack(1,k)
                  !-------------------------------------------------------------
                  !  Get the id of the particle
                  !-------------------------------------------------------------
! if we are going to run this efficiently on ivos machine, we should put a vector directive here !
                  ipart = ppm_ghosthack(j+2,k)
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                  ibuffer = ibuffer + 2
#else
                  ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + ppm_recvbuffers(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
                  pdata(1,ipart) = pdata(1,ipart) + INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                  IF (ppm_recvbuffers(ibuffer) .GT.      &
     &               (1.0_ppm_kind_single-ppm_myepss)) THEN
                     pdata(1,ipart) = .TRUE.
                  ELSE
                     pdata(1,ipart) = .FALSE.
                  ENDIF
#endif
               ENDDO
            ENDDO
         !----------------------------------------------------------------------
         !  Unrolled verion for edim=2
         !----------------------------------------------------------------------
         CASE (2)
            DO k=1,ppm_nsendlist ! or ppm_nrecvlist or ppm_ncommseq(topoid)
               DO j=1,ppm_ghosthack(1,k)
                  !-------------------------------------------------------------
                  !  Get the id of the particle
                  !-------------------------------------------------------------
! if we are going to run this efficiently on ivos machine, we should put a vector directive here !
                  ipart = ppm_ghosthack(j+2,k)
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                  ibuffer = ibuffer + 2
#else
                  ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + ppm_recvbuffers(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + ppm_recvbuffers(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  pdata(1,ipart) = pdata(2,ipart) + REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(2,ipart) = pdata(2,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(2,ipart) = pdata(2,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
                  pdata(1,ipart) = pdata(2,ipart) + INT(ppm_recvbuffers(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                  IF (ppm_recvbuffers(ibuffer) .GT.      &
     &               (1.0_ppm_kind_single-ppm_myepss)) THEN
                     pdata(1,ipart) = .TRUE.
                  ELSE
                     pdata(1,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbuffers(ibuffer) .GT.      &
     &               (1.0_ppm_kind_single-ppm_myepss)) THEN
                     pdata(2,ipart) = .TRUE.
                  ELSE
                     pdata(2,ipart) = .FALSE.
                  ENDIF
#endif
               ENDDO
            ENDDO
         !----------------------------------------------------------------------
         !  Unrolled verion for edim=3
         !----------------------------------------------------------------------
         CASE (3)
            DO k=1,ppm_nsendlist ! or ppm_nrecvlist or ppm_ncommseq(topoid)
               DO j=1,ppm_ghosthack(1,k)
                  !-------------------------------------------------------------
                  !  Get the id of the particle
                  !-------------------------------------------------------------
! if we are going to run this efficiently on ivos machine, we should put a vector directive here !
                  ipart = ppm_ghosthack(j+2,k)
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                  ibuffer = ibuffer + 2
#else
                  ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + ppm_recvbuffers(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + ppm_recvbuffers(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + ppm_recvbuffers(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(2,ipart) = pdata(2,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(3,ipart) = pdata(3,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(2,ipart) = pdata(2,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(3,ipart) = pdata(3,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
                  pdata(1,ipart) = pdata(1,ipart) + INT(ppm_recvbuffers(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + INT(ppm_recvbuffers(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                  IF (ppm_recvbuffers(ibuffer) .GT.      &
     &               (1.0_ppm_kind_single-ppm_myepss)) THEN
                     pdata(1,ipart) = .TRUE.
                  ELSE
                     pdata(1,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
               IF (ppm_recvbuffers(ibuffer) .GT.      &
     &               (1.0_ppm_kind_single-ppm_myepss)) THEN
                     pdata(2,ipart) = .TRUE.
                  ELSE
                     pdata(2,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbuffers(ibuffer) .GT.      &
     &               (1.0_ppm_kind_single-ppm_myepss)) THEN
                     pdata(3,ipart) = .TRUE.
                  ELSE
                     pdata(3,ipart) = .FALSE.
                  ENDIF
#endif
               ENDDO
            ENDDO
         !----------------------------------------------------------------------
         !  Unrolled verion for edim=4
         !----------------------------------------------------------------------
         CASE (4)
            DO k=1,ppm_nsendlist ! or ppm_nrecvlist or ppm_ncommseq(topoid)
               DO j=1,ppm_ghosthack(1,k)
                  !-------------------------------------------------------------
                  !  Get the id of the particle
                  !-------------------------------------------------------------
! if we are going to run this efficiently on ivos machine, we should put a vector directive here !
                  ipart = ppm_ghosthack(j+2,k)
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                  ibuffer = ibuffer + 2
#else
                  ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + ppm_recvbuffers(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + ppm_recvbuffers(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + ppm_recvbuffers(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(4,ipart) = pdata(4,ipart) + ppm_recvbuffers(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  pdata(4,ipart) = pdata(4,ipart) + REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(2,ipart) = pdata(2,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(3,ipart) = pdata(3,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(4,ipart) = pdata(4,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(2,ipart) = pdata(2,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(3,ipart) = pdata(3,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(4,ipart) = pdata(4,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
                  pdata(1,ipart) = pdata(1,ipart) + INT(ppm_recvbuffers(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + INT(ppm_recvbuffers(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + INT(ppm_recvbuffers(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(4,ipart) = pdata(4,ipart) + INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                  IF (ppm_recvbuffers(ibuffer) .GT.      &
     &               (1.0_ppm_kind_single-ppm_myepss)) THEN
                     pdata(1,ipart) = .TRUE.
                  ELSE
                     pdata(1,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbuffers(ibuffer) .GT.      &
     &               (1.0_ppm_kind_single-ppm_myepss)) THEN
                     pdata(2,ipart) = .TRUE.
                  ELSE
                     pdata(2,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbuffers(ibuffer) .GT.      &
     &               (1.0_ppm_kind_single-ppm_myepss)) THEN
                     pdata(3,ipart) = .TRUE.
                  ELSE
                     pdata(3,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbuffers(ibuffer) .GT.      &
     &               (1.0_ppm_kind_single-ppm_myepss)) THEN
                     pdata(4,ipart) = .TRUE.
                  ELSE
                     pdata(4,ipart) = .FALSE.
                  ENDIF
#endif
               ENDDO
            ENDDO
         !----------------------------------------------------------------------
         !  Unrolled verion for edim=5
         !----------------------------------------------------------------------
         CASE (5)
            DO k=1,ppm_nsendlist ! or ppm_nrecvlist or ppm_ncommseq(topoid)
               DO j=1,ppm_ghosthack(1,k)
                  !-------------------------------------------------------------
                  !  Get the id of the particle
                  !-------------------------------------------------------------
! if we are going to run this efficiently on ivos machine, we should put a vector directive here !
                  ipart = ppm_ghosthack(j+2,k)
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                  ibuffer = ibuffer + 2
#else
                  ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + ppm_recvbuffers(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + ppm_recvbuffers(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + ppm_recvbuffers(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(4,ipart) = pdata(4,ipart) + ppm_recvbuffers(ibuffer)
                  ibuffer = ibuffer + 1
                  pdata(5,ipart) = pdata(5,ipart) + ppm_recvbuffers(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
                  pdata(1,ipart) = pdata(1,ipart) + REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  pdata(4,ipart) = pdata(4,ipart) + REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  pdata(5,ipart) = pdata(5,ipart) + REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(2,ipart) = pdata(2,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(3,ipart) = pdata(3,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(4,ipart) = pdata(4,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_double)
                  ibuffer = ibuffer + 2
                  pdata(5,ipart) = pdata(5,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  pdata(1,ipart) = pdata(1,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(2,ipart) = pdata(2,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(3,ipart) = pdata(3,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(4,ipart) = pdata(4,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_single)
                  ibuffer = ibuffer + 2
                  pdata(5,ipart) = pdata(5,ipart) + CMPLX(ppm_recvbuffers(ibuffer-1),    &
     &                ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
                  pdata(1,ipart) = pdata(1,ipart) + INT(ppm_recvbuffers(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(2,ipart) = pdata(2,ipart) + INT(ppm_recvbuffers(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(3,ipart) = pdata(3,ipart) + INT(ppm_recvbuffers(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(4,ipart) = pdata(4,ipart) + INT(ppm_recvbuffers(ibuffer))
                  ibuffer = ibuffer + 1
                  pdata(5,ipart) = pdata(5,ipart) + INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                  IF (ppm_recvbuffers(ibuffer) .GT.      &
     &               (1.0_ppm_kind_single-ppm_myepss)) THEN
                     pdata(1,ipart) = .TRUE.
                  ELSE
                     pdata(1,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbuffers(ibuffer) .GT.      &
     &               (1.0_ppm_kind_single-ppm_myepss)) THEN
                     pdata(2,ipart) = .TRUE.
                  ELSE
                     pdata(2,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbuffers(ibuffer) .GT.      &
     &               (1.0_ppm_kind_single-ppm_myepss)) THEN
                     pdata(3,ipart) = .TRUE.
                  ELSE
                     pdata(3,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbuffers(ibuffer) .GT.      &
     &               (1.0_ppm_kind_single-ppm_myepss)) THEN
                     pdata(4,ipart) = .TRUE.
                  ELSE
                     pdata(4,ipart) = .FALSE.
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (ppm_recvbuffers(ibuffer) .GT.      &
     &               (1.0_ppm_kind_single-ppm_myepss)) THEN
                     pdata(5,ipart) = .TRUE.
                  ELSE
                     pdata(5,ipart) = .FALSE.
                  ENDIF
#endif
               ENDDO
            ENDDO

         !----------------------------------------------------------------------
         !  For edim.GT.5 the vector length will be edim !!
         !----------------------------------------------------------------------
         CASE DEFAULT
            DO k=1,ppm_nsendlist ! or ppm_nrecvlist or ppm_ncommseq(topoid)
               DO j=1,ppm_ghosthack(1,k)
                  !-------------------------------------------------------------
                  !  Get the id of the particle
                  !-------------------------------------------------------------
                  ! if we are going to run this efficiently on ivos machine,
                  ! we should put a vector directive here !
                  ipart = ppm_ghosthack(j+2,k)
                  DO i=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                     ibuffer = ibuffer + 2
#else
                     ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                     pdata(i,ipart) = pdata(i,ipart) + ppm_recvbuffers(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
                     pdata(i,ipart) = pdata(i,ipart) + REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                     pdata(i,ipart) = pdata(i,ipart) + &
                     & CMPLX(ppm_recvbuffers(ibuffer-1),ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                     pdata(i,ipart) = pdata(i,ipart) + &
                     & CMPLX(ppm_recvbuffers(ibuffer-1),ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
                     pdata(i,ipart) = pdata(i,ipart) + INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                     IF (ppm_recvbuffers(ibuffer).GT.(1._ppm_kind_single-ppm_myepss)) THEN
                        pdata(i,ipart) = .TRUE.
                     ELSE
                        pdata(i,ipart) = .FALSE.
                     ENDIF
#endif
                  ENDDO
               ENDDO
            ENDDO

         END SELECT
#elif  __DIM == 1
         !----------------------------------------------------------------------
         !  Scalar version
         !----------------------------------------------------------------------
         DO k=1,ppm_nsendlist ! or ppm_nrecvlist or ppm_ncommseq(topoid)
            DO j=1,ppm_ghosthack(1,k)
               !-------------------------------------------------------------
               !  Get the id of the particle
               !-------------------------------------------------------------
               ! if we are going to run this efficiently on ivos machine,
               ! we should put a vector directive here !
               ipart = ppm_ghosthack(j+2,k)
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
               ibuffer = ibuffer + 2
#else
               ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
               pdata(ipart) = pdata(ipart) + ppm_recvbuffers(ibuffer)
#elif  __KIND == __DOUBLE_PRECISION
               pdata(ipart) = pdata(ipart) + REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               pdata(ipart) = pdata(ipart) + &
               & CMPLX(ppm_recvbuffers(ibuffer-1),ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               pdata(ipart) = pdata(ipart) + &
               & CMPLX(ppm_recvbuffers(ibuffer-1),ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __INTEGER
               pdata(ipart) = pdata(ipart) + INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
               IF (ppm_recvbuffers(ibuffer).GT.(1._ppm_kind_single-ppm_myepss)) THEN
                  pdata(ipart) = .TRUE.
               ELSE
                  pdata(ipart) = .FALSE.
               ENDIF
#endif
            ENDDO
         ENDDO
#endif
      END SELECT ! ppm_kind
#endif

      !-------------------------------------------------------------------------
      !  Decrement the set counter
      !-------------------------------------------------------------------------
      ppm_buffer_set = ppm_buffer_set - 1

      !-------------------------------------------------------------------------
      !  Deallocate the receive buffer if all sets have been poped
      !-------------------------------------------------------------------------
      IF (ppm_buffer_set.LT.1) THEN
         iopt = ppm_param_dealloc
         IF (ppm_kind .EQ. ppm_kind_single) THEN
            CALL ppm_alloc(ppm_recvbuffers,ldu,iopt,info)
         ELSE
            CALL ppm_alloc(ppm_recvbufferd,ldu,iopt,info)
         ENDIF
         or_fail_dealloc('receive buffer PPM_RECVBUFFER')
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (ppm_buffer_set .LT. 1) THEN
            fail('buffer is empty. Cannot pop.',exit_point=8888)
          ENDIF
          IF (Npart .LT. 0) THEN
             fail('Npart must be >=0',exit_point=8888)
          ENDIF
          IF (Mpart .LT. 0) THEN
             fail('Mpart must be >=0',exit_point=8888)
          ENDIF
#if   __DIM == 2
          IF (lda .LT. 1) THEN
             fail('lda must be >0 for vector data',exit_point=8888)
          ENDIF
#elif __DIM == 1
          IF (lda .NE. 1) THEN
             fail('lda must be =1 for scalar data',exit_point=8888)
          ENDIF
#endif
      8888 CONTINUE
      END SUBROUTINE check
#if    __DIM == 1
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_ghost_pop_1ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_ghost_pop_1dd
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_ghost_pop_1dsc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_ghost_pop_1ddc
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_map_part_ghost_pop_1di
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_map_part_ghost_pop_1dl
#endif

#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_ghost_pop_2ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_ghost_pop_2dd
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_ghost_pop_2dsc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_ghost_pop_2ddc
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_map_part_ghost_pop_2di
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_map_part_ghost_pop_2dl
#endif
#endif

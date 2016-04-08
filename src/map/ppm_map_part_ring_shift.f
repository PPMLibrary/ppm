      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_map_part_ring_shift
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

#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_ring_shift_s(xp,lda,Npart,itarget,isource,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_ring_shift_d(xp,lda,Npart,itarget,isource,info)
#endif
      !!! Prepares the send and receive lists and the buffers for sending the
      !!! particle data along the ring.
      !!!
      !!! [NOTE]
      !!! The first part of the buffer contains the on processor data.
      !!! In this case all particles will be send to itarget.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
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
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      !!! Particle coordinates
      INTEGER                 , INTENT(IN   ) :: lda
      !!! The leading dimension of xp
      INTEGER                 , INTENT(IN   ) :: Npart
      !!! Number of particles
      INTEGER                 , INTENT(IN   ) :: itarget
      !!! Processor to send to
      INTEGER                 , INTENT(IN   ) :: isource
      !!! Processor to receive from
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER, DIMENSION(1) :: ldu
      INTEGER               :: i
      INTEGER               :: iopt,ibuffer,ipart

      CHARACTER(LEN=ppm_char) :: caller="ppm_map_part_ring_shift"
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  If there is still some data left in the buffer, warn the user
      !-------------------------------------------------------------------------
      IF (ppm_buffer_set .GT. 0) THEN
         fail("Buffer was not empty. Possible loss of data!",ppm_err_map_incomp, &
         & ppm_error=ppm_error_warning,exit_point=no)
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the pointer to the buffer; we only need to talk
      !  to 2 processors (including this one)
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = 3
      CALL ppm_alloc(ppm_psendbuffer,ldu,iopt,info)
      or_fail_alloc("particle send buffer PPM_PSENDBUFFER",ppm_error=ppm_error_fatal)

      iopt   = ppm_param_alloc_fit
      ldu(1) = Npart
      CALL ppm_alloc(ppm_buffer2part,ldu,iopt,info)
      or_fail_alloc("buffer-to-particles map PPM_BUFFER2PART",ppm_error=ppm_error_fatal)

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
      or_fail_alloc("buffer dimensions PPM_BUFFER_DIM",ppm_error=ppm_error_fatal)

      CALL ppm_alloc(ppm_buffer_type,ldu,iopt,info)
      or_fail_alloc("buffer types PPM_BUFFER_TYPE",ppm_error=ppm_error_fatal)

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
      or_fail_alloc("global send buffer PPM_SENDBUFFER",ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Allocate memory for the sendlist
      !-------------------------------------------------------------------------
      ppm_nsendlist = 2
      ppm_nrecvlist = 2

      iopt = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
      or_fail_alloc("send list PPM_ISENDLIST",ppm_error=ppm_error_fatal)

      ldu(1) = ppm_nrecvlist
      CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
      or_fail_alloc("receive list PPM_IRECVLIST",ppm_error=ppm_error_fatal)

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
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (.NOT. ppm_initialized) THEN
            fail("Please call ppm_init first!",ppm_err_ppm_noinit,exit_point=8888)
         ENDIF
         IF (lda .LT. 1) THEN
            fail("LDA must be >0",exit_point=8888,ppm_error=ppm_error_fatal)
         ENDIF
         IF (Npart .LT. 0) THEN
            fail("NPART must be >=0",exit_point=8888,ppm_error=ppm_error_fatal)
         ENDIF
         IF ((itarget .LT. 0) .OR. (itarget .GE. ppm_nproc)) THEN
            fail("ITARGET must be >=0 and smaller than the number processors",exit_point=8888,ppm_error=ppm_error_fatal)
         ENDIF
      8888 CONTINUE
      END SUBROUTINE check
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_ring_shift_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_ring_shift_d
#endif

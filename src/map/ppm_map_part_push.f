      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_map_part_push
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

#if    __DIM == 1
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_push_1ds(pdata,Npart,info,pushpp)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_push_1dd(pdata,Npart,info,pushpp)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_push_1dsc(pdata,Npart,info,pushpp)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_push_1ddc(pdata,Npart,info,pushpp)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_part_push_1di(pdata,Npart,info,pushpp)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_part_push_1dl(pdata,Npart,info,pushpp)
#endif
#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_push_2ds(pdata,lda,Npart,info,pushpp)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_push_2dd(pdata,lda,Npart,info,pushpp)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_push_2dsc(pdata,lda,Npart,info,pushpp)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_push_2ddc(pdata,lda,Npart,info,pushpp)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_part_push_2di(pdata,lda,Npart,info,pushpp)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_part_push_2dl(pdata,lda,Npart,info,pushpp)
#endif
#endif
      !!! This routine pushes particle data onto the send buffer.
      !!!
      !!! [NOTE]
      !!! The first part of the buffer contains the on processor data.
      !!! The packing could be performed more efficiently.
      !!!
      !!! [WARNING]
      !!! DIM is the dimension of the pdata array and not the space
      !!! dimension ppm_dim!
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
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
      INTEGER , DIMENSION(:  ), INTENT(IN   )    :: pdata
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:  ), INTENT(IN   )    :: pdata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:  ), INTENT(IN   ) :: pdata
#else
      REAL(MK), DIMENSION(:  ), INTENT(IN   )    :: pdata
#endif

#elif __DIM == 2
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:,:), INTENT(IN   )    :: pdata
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:), INTENT(IN   )    :: pdata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:), INTENT(IN   ) :: pdata
#else
      REAL(MK), DIMENSION(:,:), INTENT(IN   )    :: pdata
#endif
#endif
      !!! Particle data.
      !!! Can be either 1D or 2D array.
#if   __DIM == 2
      INTEGER                 , INTENT(IN   )    :: lda
      !!! The leading dimension of the pdata.
#endif
      INTEGER                 , INTENT(IN   )    :: Npart
      !!! Number of particles
      INTEGER                 , INTENT(  OUT)    :: info
      !!! Returns 0 upon success
      LOGICAL, OPTIONAL       , INTENT(IN   )    :: pushpp
      !!! If `TRUE` then pdata is assumed to contain the particle positions
      !!! (xp)
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER, DIMENSION(3) :: ldu
      INTEGER               :: i,j,k,ipart,ibuffer,icount
      INTEGER               :: iopt,ldb,incr
#if   __DIM == 1
      INTEGER, PARAMETER    :: lda = 1
#endif

      CHARACTER(LEN=ppm_char) :: caller='ppm_map_part_push'

      LOGICAL :: lpushpp
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
      !  The default is that we do not push the particle positions
      !-------------------------------------------------------------------------
      lpushpp = MERGE(pushpp,.FALSE.,PRESENT(pushpp))

      !-------------------------------------------------------------------------
      !  Increment the buffer set
      !-------------------------------------------------------------------------
      ppm_buffer_set = ppm_buffer_set + 1

      !-------------------------------------------------------------------------
      !  Allocate memory for the buffer dimension and type
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_grow_preserve
      ldu(1)  = ppm_buffer_set
      CALL ppm_alloc(ppm_buffer_dim ,ldu,iopt,info)
      or_fail_alloc('buffer dimensions PPM_BUFFER_DIM',ppm_error=ppm_error_fatal)

      CALL ppm_alloc(ppm_buffer_type,ldu,iopt,info)
      or_fail_alloc('buffer types PPM_BUFFER_TYPE',ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  A complex number is treated as two reals. Cannot change lda
      !  because it is INTENT(IN)
      !-------------------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      ldb = 2*lda
#else
      ldb = lda
#endif

      !-------------------------------------------------------------------------
      !  Store the dimension and type
      !-------------------------------------------------------------------------
      ppm_buffer_dim(ppm_buffer_set)  = ldb
#if    __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      ppm_buffer_type(ppm_buffer_set) = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_PRECISION_COMPLEX
      ppm_buffer_type(ppm_buffer_set) = ppm_kind_double
#elif __KIND == __INTEGER
      ppm_buffer_type(ppm_buffer_set) = ppm_integer
#elif __KIND == __LOGICAL
      ppm_buffer_type(ppm_buffer_set) = ppm_logical
#endif

      !-------------------------------------------------------------------------
      !  loop over the processors in the ppm_isendlist()
      !-------------------------------------------------------------------------
      ibuffer = ppm_nsendbuffer

      !-------------------------------------------------------------------------
      !  DOUBLE PRECISION BUFFER
      !-------------------------------------------------------------------------
      SELECT CASE (ppm_kind)
      CASE (ppm_kind_double)
         !----------------------------------------------------------------------
         !  (Re)allocate memory for the buffer
         !----------------------------------------------------------------------
         incr = 0
         DO i=1,ppm_nsendlist
            incr = incr + (ppm_psendbuffer(i+1)-ppm_psendbuffer(i))*ldb
         ENDDO

         ldu(1) = ppm_nsendbuffer + incr
         iopt   = ppm_param_alloc_grow_preserve
         CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
         or_fail_alloc('global send buffer PPM_SENDBUFFERD',ppm_error=ppm_error_fatal)

         DO i=1,ppm_nsendlist
            !-------------------------------------------------------------------
            !  access the particles belonging to the i-th processor in the
            !  sendlist
            !-------------------------------------------------------------------
            !-------------------------------------------------------------------
            !  Store the particle data in the buffer
            !-------------------------------------------------------------------
#if   __DIM == 2
            !-------------------------------------------------------------------
            !  Unrolled for lda=1
            !-------------------------------------------------------------------
            SELECT CASE (lda)
            CASE (1)
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                  ppm_sendbufferd(ibuffer) = pdata(1,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(1,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),ppm_kind_double)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Unrolled for lda=2
            !-------------------------------------------------------------------
            CASE (2)
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                  ppm_sendbufferd(ibuffer) = pdata(1,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(2,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(2,ipart)),  &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(1,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(2,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(2,ipart),ppm_kind_double)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(2,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Unrolled for lda=3
            !-------------------------------------------------------------------
            CASE (3)
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                  ppm_sendbufferd(ibuffer) = pdata(1,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(2,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(3,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(2,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(3,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(3,ipart)),  &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(1,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(2,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(3,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(3,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(3,ipart),ppm_kind_double)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(2,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(3,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Unrolled for lda=4
            !-------------------------------------------------------------------
            CASE (4)
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                  ppm_sendbufferd(ibuffer) = pdata(1,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(2,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(3,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(4,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(2,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(3,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(3,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(4,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(4,ipart)),  &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(1,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(2,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(3,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(3,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(4,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(4,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(3,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(4,ipart),ppm_kind_double)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(2,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(3,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(4,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Unrolled for lda=5
            !-------------------------------------------------------------------
            CASE (5)
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(5,ipart),   &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                  ppm_sendbufferd(ibuffer) = pdata(1,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(2,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(3,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(4,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(5,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(2,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(3,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(3,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(4,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(4,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(5,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(5,ipart)),  &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(1,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(2,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(3,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(3,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(4,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(4,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(5,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(5,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(3,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(4,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(5,ipart),ppm_kind_double)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(2,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(3,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(4,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(5,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Not unrolled for the rest. Vector length will be lda!!
            !-------------------------------------------------------------------
            CASE DEFAULT
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  DO k=1,lda
                     ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                     ppm_sendbufferd(ibuffer) = REAL(pdata(k,ipart),   &
     &                   ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                     ppm_sendbufferd(ibuffer) = pdata(k,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                     ppm_sendbufferd(ibuffer) =         &
     &                          REAL(pdata(k,ipart),ppm_kind_double)
                     ibuffer = ibuffer + 1
                     ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(k,ipart)),  &
     &                   ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                     ppm_sendbufferd(ibuffer) =     &
     &                          REAL(pdata(k,ipart),ppm_kind_double)
                     ibuffer = ibuffer + 1
                     ppm_sendbufferd(ibuffer) = AIMAG(pdata(k,ipart))
#elif  __KIND == __INTEGER
                     ppm_sendbufferd(ibuffer) = REAL(pdata(k,ipart),ppm_kind_double)
#elif  __KIND == __LOGICAL
                     IF (pdata(k,ipart)) THEN
                        ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                     ELSE
                        ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                     ENDIF
#endif
                  ENDDO
               ENDDO

            END SELECT
#elif  __DIM == 1
            !-------------------------------------------------------------------
            !  Scalar version
            !-------------------------------------------------------------------
            DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
               !----------------------------------------------------------------
               !  Get the particle id
               !----------------------------------------------------------------
               ipart = ppm_buffer2part(j)
               ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
               ppm_sendbufferd(ibuffer) = REAL(pdata(ipart),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
               ppm_sendbufferd(ibuffer) = pdata(ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               ppm_sendbufferd(ibuffer) = REAL(pdata(ipart),ppm_kind_double)
               ibuffer = ibuffer + 1
               ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(ipart)),  &
     &             ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               ppm_sendbufferd(ibuffer) = REAL(pdata(ipart),ppm_kind_double)
               ibuffer = ibuffer + 1
               ppm_sendbufferd(ibuffer) = AIMAG(pdata(ipart))
#elif  __KIND == __INTEGER
               ppm_sendbufferd(ibuffer) = REAL(pdata(ipart),ppm_kind_double)
#elif  __KIND == __LOGICAL
               IF (pdata(ipart)) THEN
                  ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
               ELSE
                  ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
               ENDIF
#endif
            ENDDO
#endif
         ENDDO                ! i=1,ppm_nsendlist
      !-------------------------------------------------------------------------
      !  SINGLE PRECISION BUFFER
      !-------------------------------------------------------------------------
      CASE DEFAULT
         !----------------------------------------------------------------------
         !  (Re)allocate memory for the buffer
         !----------------------------------------------------------------------
         incr = 0
         DO i=1,ppm_nsendlist
            incr = incr + (ppm_psendbuffer(i+1)-ppm_psendbuffer(i))*ldb
         ENDDO
         ldu(1) = ppm_nsendbuffer + incr
         iopt   = ppm_param_alloc_grow_preserve
         CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
         or_fail_alloc('global send buffer PPM_SENDBUFFERS',ppm_error=ppm_error_fatal)

         DO i=1,ppm_nsendlist
            !-------------------------------------------------------------------
            !  access the particles belonging to the i-th processor in the
            !  sendlist
            !-------------------------------------------------------------------
            !-------------------------------------------------------------------
            !  Store the particle data in the buffer
            !-------------------------------------------------------------------
#if   __DIM == 2
            !-------------------------------------------------------------------
            !  Unrolled for lda=1
            !-------------------------------------------------------------------
            SELECT CASE (lda)
            CASE (1)
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                  ppm_sendbuffers(ibuffer) = pdata(1,ipart)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(1,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbuffers(ibuffer) = TRANSFER(pdata(1,ipart),1.0_ppm_kind_single)
!                   & REAL(pdata(1,ipart),ppm_kind_single)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Unrolled for lda=2
            !-------------------------------------------------------------------
            CASE (2)
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                  ppm_sendbuffers(ibuffer) = pdata(1,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(2,ipart)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(2,ipart)),  &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(1,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(2,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbuffers(ibuffer) = TRANSFER(pdata(1,ipart),1.0_ppm_kind_single)
!                   REAL(pdata(1,ipart),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = TRANSFER(pdata(2,ipart),1.0_ppm_kind_single)
!                   REAL(pdata(2,ipart),ppm_kind_single)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(2,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Unrolled for lda=3
            !-------------------------------------------------------------------
            CASE (3)
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                  ppm_sendbuffers(ibuffer) = pdata(1,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(2,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(3,ipart)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(2,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(3,ipart)),  &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(1,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(2,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(3,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbuffers(ibuffer) = TRANSFER(pdata(1,ipart),1.0_ppm_kind_single)
!                   REAL(pdata(1,ipart),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = TRANSFER(pdata(2,ipart),1.0_ppm_kind_single)
!                   REAL(pdata(2,ipart),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = TRANSFER(pdata(3,ipart),1.0_ppm_kind_single)
!                   REAL(pdata(3,ipart),ppm_kind_single)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(2,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(3,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Unrolled for lda=4
            !-------------------------------------------------------------------
            CASE (4)
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                  ppm_sendbuffers(ibuffer) = pdata(1,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(2,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(3,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(4,ipart)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(2,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(3,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(4,ipart)),  &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(1,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(2,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(3,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(4,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbuffers(ibuffer) = TRANSFER(pdata(1,ipart),1.0_ppm_kind_single)
!                   REAL(pdata(1,ipart),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = TRANSFER(pdata(2,ipart),1.0_ppm_kind_single)
!                   REAL(pdata(2,ipart),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = TRANSFER(pdata(3,ipart),1.0_ppm_kind_single)
!                   REAL(pdata(3,ipart),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = TRANSFER(pdata(4,ipart),1.0_ppm_kind_single)
!                   REAL(pdata(4,ipart),ppm_kind_single)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(2,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(3,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(4,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Unrolled for lda=5
            !-------------------------------------------------------------------
            CASE (5)
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(5,ipart),   &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                  ppm_sendbuffers(ibuffer) = pdata(1,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(2,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(3,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(4,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(5,ipart)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(2,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(3,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(4,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(5,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(5,ipart)),  &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(1,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(2,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(3,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(4,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(5,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(5,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbuffers(ibuffer) = TRANSFER(pdata(1,ipart),1.0_ppm_kind_single)
!                   REAL(pdata(1,ipart),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = TRANSFER(pdata(2,ipart),1.0_ppm_kind_single)
!                   REAL(pdata(2,ipart),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = TRANSFER(pdata(3,ipart),1.0_ppm_kind_single)
!                   REAL(pdata(3,ipart),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = TRANSFER(pdata(4,ipart),1.0_ppm_kind_single)
!                   REAL(pdata(4,ipart),ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = TRANSFER(pdata(5,ipart),1.0_ppm_kind_single)
!                   REAL(pdata(5,ipart),ppm_kind_single)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(2,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(3,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(4,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(5,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
#endif
               ENDDO

            !-------------------------------------------------------------------
            !  Not unrolled. Vector length will be lda !!
            !-------------------------------------------------------------------
            CASE DEFAULT
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  DO k=1,lda
                     ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                     ppm_sendbuffers(ibuffer) = REAL(pdata(k,ipart),   &
     &                   ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                     ppm_sendbuffers(ibuffer) = pdata(k,ipart)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                     ppm_sendbuffers(ibuffer) = REAL(pdata(k,ipart),   &
     &                   ppm_kind_single)
                     ibuffer = ibuffer + 1
                     ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(k,ipart)),  &
     &                   ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                     ppm_sendbuffers(ibuffer) = REAL(pdata(k,ipart),   &
     &                   ppm_kind_single)
                     ibuffer = ibuffer + 1
                     ppm_sendbuffers(ibuffer) = AIMAG(pdata(k,ipart))
#elif  __KIND == __INTEGER
                     ppm_sendbuffers(ibuffer) = TRANSFER(pdata(k,ipart),1.0_ppm_kind_single)
!                      REAL(pdata(k,ipart),ppm_kind_single)
#elif  __KIND == __LOGICAL
                     IF (pdata(k,ipart)) THEN
                        ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                     ELSE
                        ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                     ENDIF
#endif
                  ENDDO
               ENDDO

            END SELECT
#elif  __DIM == 1
            !-------------------------------------------------------------------
            !  Scalar version
            !-------------------------------------------------------------------
            DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
               !----------------------------------------------------------------
               !  Get the particle id
               !----------------------------------------------------------------
               ipart = ppm_buffer2part(j)
               ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
               ppm_sendbuffers(ibuffer) = REAL(pdata(ipart),ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
               ppm_sendbuffers(ibuffer) = pdata(ipart)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               ppm_sendbuffers(ibuffer) = REAL(pdata(ipart),ppm_kind_single)
               ibuffer = ibuffer + 1
               ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(ipart)),  &
     &             ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               ppm_sendbuffers(ibuffer) = REAL(pdata(ipart),ppm_kind_single)
               ibuffer = ibuffer + 1
               ppm_sendbuffers(ibuffer) = AIMAG(pdata(ipart))
#elif  __KIND == __INTEGER
               ppm_sendbuffers(ibuffer) = TRANSFER(pdata(ipart),1.0_ppm_kind_single)
!                REAL(pdata(ipart),ppm_kind_single)
#elif  __KIND == __LOGICAL
               IF (pdata(ipart)) THEN
                  ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
               ELSE
                  ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
               ENDIF
#endif
            ENDDO
#endif
         ENDDO ! i=1,ppm_nsendlist
      END SELECT ! ppm_kind

      !-------------------------------------------------------------------------
      !  If we are pushing particle positions (if lpushpp is true) we need to
      !  add the offset to the particles; pushing the particles can only occur
      !  for 2D arrays in single or double precions: xp is USUALLY stored as
      !  xp(1:ppm_dim,1:Npart)
      !-------------------------------------------------------------------------
      IF (lpushpp) THEN
         !----------------------------------------------------------------------
         !  the particle positions are per construction (by a call to ghost_get)
         !  ALWAYS stored from ibuffer = 1 to nsendbuffer(2) - 1 and the
         !  ppm_ghost_offset is therefore also stored from 1 - nsendbuffer(2)-1
         !  thus icount = 1 and ibuffer = ppm_nsendlist (which has not yet been
         !  updated (see below)
         !----------------------------------------------------------------------
         icount  = 0
         ibuffer = ppm_nsendbuffer
         SELECT CASE (ppm_kind)
         CASE (ppm_kind_double)
            DO i=1,ppm_nsendlist
               IF (lda.EQ.2) THEN
                  !-------------------------------------------------------------
                  !  Unrolled for lda=2
                  !-------------------------------------------------------------
                  DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                     !----------------------------------------------------------
                     !  Get the particle id
                     !----------------------------------------------------------
                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbufferd(ibuffer) = &
                     & ppm_sendbufferd(ibuffer)*ppm_ghost_offset_facd(icount)+ &
                     & ppm_ghost_offsetd(icount)

                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbufferd(ibuffer) = &
                     & ppm_sendbufferd(ibuffer)*ppm_ghost_offset_facd(icount)+ &
                     & ppm_ghost_offsetd(icount)
                  ENDDO
               ELSEIF (lda.EQ.3) THEN
                  !-------------------------------------------------------------
                  !  Unrolled for lda=3
                  !-------------------------------------------------------------
                  DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                     !----------------------------------------------------------
                     !  Get the particle id
                     !----------------------------------------------------------
                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbufferd(ibuffer) = &
                     & ppm_sendbufferd(ibuffer)*ppm_ghost_offset_facd(icount)+ &
                     & ppm_ghost_offsetd(icount)

                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbufferd(ibuffer) = &
                     & ppm_sendbufferd(ibuffer)*ppm_ghost_offset_facd(icount)+ &
                     & ppm_ghost_offsetd(icount)

                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbufferd(ibuffer) = &
                     & ppm_sendbufferd(ibuffer)*ppm_ghost_offset_facd(icount)+ &
                     & ppm_ghost_offsetd(icount)
                  ENDDO
               ENDIF ! end of lda = 3
            ENDDO ! enddo of nsendlist

         CASE DEFAULT ! buffers are single precision
            DO i=1,ppm_nsendlist
               IF (lda.EQ.2) THEN
                  !-------------------------------------------------------------
                  !  Unrolled for lda=2
                  !-------------------------------------------------------------
                  DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                     !----------------------------------------------------------
                     !  Get the particle id
                     !----------------------------------------------------------
                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbuffers(ibuffer) = &
                     & ppm_sendbuffers(ibuffer)*ppm_ghost_offset_facs(icount)+ &
                     & ppm_ghost_offsets(icount)

                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbuffers(ibuffer) = &
                     & ppm_sendbuffers(ibuffer)*ppm_ghost_offset_facs(icount)+ &
                     & ppm_ghost_offsets(icount)
                  ENDDO
               ELSEIF (lda.EQ.3) THEN
                  !-------------------------------------------------------------
                  !  Unrolled for lda=3
                  !-------------------------------------------------------------
                  DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                     !----------------------------------------------------------
                     !  Get the particle id
                     !----------------------------------------------------------
                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbuffers(ibuffer) = &
                     & ppm_sendbuffers(ibuffer)*ppm_ghost_offset_facs(icount)+ &
                     & ppm_ghost_offsets(icount)

                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbuffers(ibuffer) = &
                     & ppm_sendbuffers(ibuffer)*ppm_ghost_offset_facs(icount)+ &
                     & ppm_ghost_offsets(icount)

                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbuffers(ibuffer) = &
                     & ppm_sendbuffers(ibuffer)*ppm_ghost_offset_facs(icount)+ &
                     & ppm_ghost_offsets(icount)
                  ENDDO
               ENDIF ! end of lda = 3
            ENDDO ! enddo of nsendlist

         END SELECT ! end of single/double precision buffers
      ENDIF ! end of lpushpp

      !-------------------------------------------------------------------------
      !  Update the particle buffer count
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
#if   __DIM == 2
          IF (lda .LT. 1) THEN
             fail('lda must be >0 for vector data',exit_point=8888)
          ENDIF
#elif __DIM == 1
          IF (lda .NE. 1) THEN
             fail('lda must be =1 for scalar data',exit_point=8888)
          ENDIF
#endif
          IF (Npart .LT. 0) THEN
             fail('Npart must be >=0',exit_point=8888)
          ENDIF
      8888 CONTINUE
      END SUBROUTINE check
#if   __DIM == 1
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_push_1ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_push_1dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_push_1dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_push_1ddc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_map_part_push_1di
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_map_part_push_1dl
#endif

#elif __DIM == 2
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_push_2ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_push_2dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_push_2dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_push_2ddc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_map_part_push_2di
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_map_part_push_2dl
#endif
#endif

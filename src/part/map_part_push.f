      !-------------------------------------------------------------------------
      !  Subroutine   :                  map_part_push
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

#if    __DIM == 1
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE DTYPE(map_part_push_1d)(Pc,mapID,pdata,info,pushpp)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE DTYPE(map_part_push_1d)(Pc,mapID,pdata,info,pushpp)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE DTYPE(map_part_push_1dc)(Pc,mapID,pdata,info,pushpp)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE DTYPE(map_part_push_1dc)(Pc,mapID,pdata,info,pushpp)
#elif  __KIND == __INTEGER
      SUBROUTINE DTYPE(map_part_push_1di)(Pc,mapID,pdata,info,pushpp)
#elif  __KIND == __LOGICAL
      SUBROUTINE DTYPE(map_part_push_1dl)(Pc,mapID,pdata,info,pushpp)
#endif
#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE DTYPE(map_part_push_2d)(Pc,mapID,pdata,lda,info,pushpp)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE DTYPE(map_part_push_2d)(Pc,mapID,pdata,lda,info,pushpp)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE DTYPE(map_part_push_2dc)(Pc,mapID,pdata,lda,info,pushpp)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE DTYPE(map_part_push_2dc)(Pc,mapID,pdata,lda,info,pushpp)
#elif  __KIND == __INTEGER
      SUBROUTINE DTYPE(map_part_push_2di)(Pc,mapID,pdata,lda,info,pushpp)
#elif  __KIND == __LOGICAL
      SUBROUTINE DTYPE(map_part_push_2dl)(Pc,mapID,pdata,lda,info,pushpp)
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
      CLASS(DTYPE(ppm_t_particles))              :: Pc
      INTEGER,                  INTENT(IN   )    :: mapID
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
#elif __KIND == __SINGLE_PRECISION_COMPLEX | \
      __KIND == __DOUBLE_PRECISION_COMPLEX
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
      INTEGER                 , INTENT(  OUT)    :: info
      !!! Returns 0 upon success
      LOGICAL, OPTIONAL       , INTENT(IN   )    :: pushpp
      !!! If `TRUE` then pdata is assumed to contain the particle positions
      !!! (xp)
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3) :: ldu
      INTEGER               :: i,j,k,ipart,ibuffer,icount
      INTEGER               :: iopt,ldb,incr
      REAL(MK)              :: t0
      LOGICAL               :: lpushpp
      CHARACTER(LEN=ppm_char) :: caller ='map_part_push'
#if   __DIM == 1
      INTEGER, PARAMETER    :: lda = 1
#endif
      TYPE(DTYPE(ppm_t_part_mapping)), POINTER :: map => NULL()
      REAL(MK),DIMENSION(:),POINTER :: ppm_sendbuffer => NULL()
      INTEGER, DIMENSION(:),POINTER :: ppm_buffer2part => NULL()
      REAL(MK),DIMENSION(:),POINTER :: ppm_ghost_offset => NULL()
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
      ppm_buffer2part => map%ppm_buffer2part
      ppm_ghost_offset => map%ppm_ghost_offset
      !-------------------------------------------------------------------------
      !  The default is that we do not push the particle positions
      !-------------------------------------------------------------------------
      IF (PRESENT(pushpp)) THEN
         lpushpp = pushpp
      ELSE
         lpushpp = .FALSE.
      ENDIF 

      !-------------------------------------------------------------------------
      !  Increment the buffer set 
      !-------------------------------------------------------------------------
      map%ppm_buffer_set = map%ppm_buffer_set + 1

      !-------------------------------------------------------------------------
      !  Allocate memory for the buffer dimension and type
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_grow_preserve
      ldu(1)  = map%ppm_buffer_set
      CALL ppm_alloc(map%ppm_buffer_dim ,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,caller,     &
     &        'buffer dimensions ppm_buffer_dim',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(map%ppm_buffer_type,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,caller,     &
     &        'buffer types PPM_BUFFER_TYPE',__LINE__,info)
          GOTO 9999
      ENDIF

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
      map%ppm_buffer_dim(map%ppm_buffer_set)  = ldb
#if    __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      map%ppm_buffer_type(map%ppm_buffer_set) = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_PRECISION_COMPLEX
      map%ppm_buffer_type(map%ppm_buffer_set) = ppm_kind_double
#elif __KIND == __INTEGER
      map%ppm_buffer_type(map%ppm_buffer_set) = ppm_integer
#elif __KIND == __LOGICAL
      map%ppm_buffer_type(map%ppm_buffer_set) = ppm_logical
#endif

      !-------------------------------------------------------------------------
      !  loop over the processors in the ppm_isendlist() 
      !-------------------------------------------------------------------------
      ibuffer = map%ppm_nsendbuffer

     !----------------------------------------------------------------------
     !  (Re)allocate memory for the buffer 
     !----------------------------------------------------------------------
     incr = 0 
     DO i=1,map%ppm_nsendlist
        incr = incr + (map%ppm_psendbuffer(i+1)-map%ppm_psendbuffer(i))*ldb
     ENDDO
     ldu(1) = map%ppm_nsendbuffer + incr
     iopt   = ppm_param_alloc_grow_preserve
     CALL ppm_alloc(map%ppm_sendbuffer,ldu,iopt,info)
     ppm_sendbuffer => map%ppm_sendbuffer
     IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,caller,     &
 &           'global send buffer ppm_sendbuffer',__LINE__,info)
         GOTO 9999
     ENDIF

     DO i=1,map%ppm_nsendlist
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
        IF (lda .EQ. 1) THEN
           DO j=map%ppm_psendbuffer(i),map%ppm_psendbuffer(i+1)-1
              !-------------------------------------------------------------
              !  Get the particle id 
              !-------------------------------------------------------------
              ipart = ppm_buffer2part(j)
              ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
#elif  __KIND == __DOUBLE_PRECISION
              ppm_sendbuffer(ibuffer) = pdata(1,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(AIMAG(pdata(1,ipart)),MK)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = AIMAG(pdata(1,ipart))
#elif  __KIND == __INTEGER
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
#elif  __KIND == __LOGICAL
              IF (pdata(1,ipart)) THEN
                 ppm_sendbuffer(ibuffer) = 1.0_MK
              ELSE
                 ppm_sendbuffer(ibuffer) = 0.0_MK
              ENDIF
#endif
           ENDDO
        !-------------------------------------------------------------------
        !  Unrolled for lda=2
        !-------------------------------------------------------------------
        ELSEIF (lda .EQ. 2) THEN
           DO j=map%ppm_psendbuffer(i),map%ppm_psendbuffer(i+1)-1
              !-------------------------------------------------------------
              !  Get the particle id 
              !-------------------------------------------------------------
              ipart = ppm_buffer2part(j)
              ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(2,ipart),MK) 
#elif  __KIND == __DOUBLE_PRECISION
              ppm_sendbuffer(ibuffer) = pdata(1,ipart)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = pdata(2,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(AIMAG(pdata(1,ipart)),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(2,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(AIMAG(pdata(2,ipart)),MK)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = AIMAG(pdata(1,ipart))
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) =     &
 &                       REAL(pdata(2,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = AIMAG(pdata(2,ipart))
#elif  __KIND == __INTEGER
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(2,ipart),MK)
#elif  __KIND == __LOGICAL
              IF (pdata(1,ipart)) THEN
                 ppm_sendbuffer(ibuffer) = 1.0_MK
              ELSE
                 ppm_sendbuffer(ibuffer) = 0.0_MK
              ENDIF
              ibuffer = ibuffer + 1
              IF (pdata(2,ipart)) THEN
                 ppm_sendbuffer(ibuffer) = 1.0_MK
              ELSE
                 ppm_sendbuffer(ibuffer) = 0.0_MK
              ENDIF
#endif
           ENDDO
        !-------------------------------------------------------------------
        !  Unrolled for lda=3
        !-------------------------------------------------------------------
        ELSEIF (lda .EQ. 3) THEN
           DO j=map%ppm_psendbuffer(i),map%ppm_psendbuffer(i+1)-1
              !-------------------------------------------------------------
              !  Get the particle id 
              !-------------------------------------------------------------
              ipart = ppm_buffer2part(j)
              ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(2,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(3,ipart),MK)
#elif  __KIND == __DOUBLE_PRECISION
              ppm_sendbuffer(ibuffer) = pdata(1,ipart)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = pdata(2,ipart)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = pdata(3,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(AIMAG(pdata(1,ipart)),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(2,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(AIMAG(pdata(2,ipart)),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(3,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(AIMAG(pdata(3,ipart)),MK)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = AIMAG(pdata(1,ipart))
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(2,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = AIMAG(pdata(2,ipart))
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(3,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = AIMAG(pdata(3,ipart))
#elif  __KIND == __INTEGER
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK) 
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(2,ipart),MK) 
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(3,ipart),MK)
#elif  __KIND == __LOGICAL
              IF (pdata(1,ipart)) THEN
                 ppm_sendbuffer(ibuffer) = 1.0_MK
              ELSE
                 ppm_sendbuffer(ibuffer) = 0.0_MK
              ENDIF
              ibuffer = ibuffer + 1
              IF (pdata(2,ipart)) THEN
                 ppm_sendbuffer(ibuffer) = 1.0_MK
              ELSE
                 ppm_sendbuffer(ibuffer) = 0.0_MK
              ENDIF
              ibuffer = ibuffer + 1
              IF (pdata(3,ipart)) THEN
                 ppm_sendbuffer(ibuffer) = 1.0_MK
              ELSE
                 ppm_sendbuffer(ibuffer) = 0.0_MK
              ENDIF
#endif
           ENDDO
        !-------------------------------------------------------------------
        !  Unrolled for lda=4
        !-------------------------------------------------------------------
        ELSEIF (lda .EQ. 4) THEN
           DO j=map%ppm_psendbuffer(i),map%ppm_psendbuffer(i+1)-1
              !-------------------------------------------------------------
              !  Get the particle id 
              !-------------------------------------------------------------
              ipart = ppm_buffer2part(j)
              ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(2,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(3,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(4,ipart),MK)
#elif  __KIND == __DOUBLE_PRECISION
              ppm_sendbuffer(ibuffer) = pdata(1,ipart)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = pdata(2,ipart)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = pdata(3,ipart)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = pdata(4,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(AIMAG(pdata(1,ipart)),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(2,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(AIMAG(pdata(2,ipart)),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(3,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(AIMAG(pdata(3,ipart)),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(4,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(AIMAG(pdata(4,ipart)),MK)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = AIMAG(pdata(1,ipart))
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(2,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = AIMAG(pdata(2,ipart))
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(3,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = AIMAG(pdata(3,ipart))
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(4,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = AIMAG(pdata(4,ipart))
#elif  __KIND == __INTEGER
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(2,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(3,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(4,ipart),MK)
#elif  __KIND == __LOGICAL
              IF (pdata(1,ipart)) THEN
                 ppm_sendbuffer(ibuffer) = 1.0_MK
              ELSE
                 ppm_sendbuffer(ibuffer) = 0.0_MK
              ENDIF
              ibuffer = ibuffer + 1
              IF (pdata(2,ipart)) THEN
                 ppm_sendbuffer(ibuffer) = 1.0_MK
              ELSE
                 ppm_sendbuffer(ibuffer) = 0.0_MK
              ENDIF
              ibuffer = ibuffer + 1
              IF (pdata(3,ipart)) THEN
                 ppm_sendbuffer(ibuffer) = 1.0_MK
              ELSE
                 ppm_sendbuffer(ibuffer) = 0.0_MK
              ENDIF
              ibuffer = ibuffer + 1
              IF (pdata(4,ipart)) THEN
                 ppm_sendbuffer(ibuffer) = 1.0_MK
              ELSE
                 ppm_sendbuffer(ibuffer) = 0.0_MK
              ENDIF
#endif
           ENDDO
        !-------------------------------------------------------------------
        !  Unrolled for lda=5
        !-------------------------------------------------------------------
        ELSEIF (lda .EQ. 5) THEN
           DO j=map%ppm_psendbuffer(i),map%ppm_psendbuffer(i+1)-1
              !-------------------------------------------------------------
              !  Get the particle id 
              !-------------------------------------------------------------
              ipart = ppm_buffer2part(j)
              ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(2,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(3,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(4,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(5,ipart),MK)
#elif  __KIND == __DOUBLE_PRECISION
              ppm_sendbuffer(ibuffer) = pdata(1,ipart)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = pdata(2,ipart)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = pdata(3,ipart)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = pdata(4,ipart)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = pdata(5,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(AIMAG(pdata(1,ipart)),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(2,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(AIMAG(pdata(2,ipart)),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(3,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(AIMAG(pdata(3,ipart)),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(4,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(AIMAG(pdata(4,ipart)),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(5,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(AIMAG(pdata(5,ipart)),MK)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = AIMAG(pdata(1,ipart))
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(2,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = AIMAG(pdata(2,ipart))
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(3,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = AIMAG(pdata(3,ipart))
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(4,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = AIMAG(pdata(4,ipart))
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(5,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = AIMAG(pdata(5,ipart))
#elif  __KIND == __INTEGER
              ppm_sendbuffer(ibuffer) = REAL(pdata(1,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(2,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(3,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(4,ipart),MK)
              ibuffer = ibuffer + 1
              ppm_sendbuffer(ibuffer) = REAL(pdata(5,ipart),MK)
#elif  __KIND == __LOGICAL
              IF (pdata(1,ipart)) THEN
                 ppm_sendbuffer(ibuffer) = 1.0_MK
              ELSE
                 ppm_sendbuffer(ibuffer) = 0.0_MK
              ENDIF
              ibuffer = ibuffer + 1
              IF (pdata(2,ipart)) THEN
                 ppm_sendbuffer(ibuffer) = 1.0_MK
              ELSE
                 ppm_sendbuffer(ibuffer) = 0.0_MK
              ENDIF
              ibuffer = ibuffer + 1
              IF (pdata(3,ipart)) THEN
                 ppm_sendbuffer(ibuffer) = 1.0_MK
              ELSE
                 ppm_sendbuffer(ibuffer) = 0.0_MK
              ENDIF
              ibuffer = ibuffer + 1
              IF (pdata(4,ipart)) THEN
                 ppm_sendbuffer(ibuffer) = 1.0_MK
              ELSE
                 ppm_sendbuffer(ibuffer) = 0.0_MK
              ENDIF
              ibuffer = ibuffer + 1
              IF (pdata(5,ipart)) THEN
                 ppm_sendbuffer(ibuffer) = 1.0_MK
              ELSE
                 ppm_sendbuffer(ibuffer) = 0.0_MK
              ENDIF
#endif
           ENDDO
        !-------------------------------------------------------------------
        !  Not unrolled for the rest. Vector length will be lda!!
        !-------------------------------------------------------------------
        ELSE
           DO j=map%ppm_psendbuffer(i),map%ppm_psendbuffer(i+1)-1
              !-------------------------------------------------------------
              !  Get the particle id 
              !-------------------------------------------------------------
              ipart = ppm_buffer2part(j)
              DO k=1,lda
                 ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                 ppm_sendbuffer(ibuffer) = REAL(pdata(k,ipart),MK)
#elif  __KIND == __DOUBLE_PRECISION
                 ppm_sendbuffer(ibuffer) = pdata(k,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                 ppm_sendbuffer(ibuffer) = REAL(pdata(k,ipart),MK)
                 ibuffer = ibuffer + 1
                 ppm_sendbuffer(ibuffer) = REAL(AIMAG(pdata(k,ipart)),MK)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                 ppm_sendbuffer(ibuffer) = REAL(pdata(k,ipart),MK)
                 ibuffer = ibuffer + 1
                 ppm_sendbuffer(ibuffer) = AIMAG(pdata(k,ipart))
#elif  __KIND == __INTEGER
                 ppm_sendbuffer(ibuffer) = REAL(pdata(k,ipart),MK)
#elif  __KIND == __LOGICAL
                 IF (pdata(k,ipart)) THEN
                    ppm_sendbuffer(ibuffer) = 1.0_MK
                 ELSE
                    ppm_sendbuffer(ibuffer) = 0.0_MK
                 ENDIF
#endif
              ENDDO
           ENDDO
        ENDIF
#elif  __DIM == 1
        !-------------------------------------------------------------------
        !  Scalar version
        !-------------------------------------------------------------------
        DO j=map%ppm_psendbuffer(i),map%ppm_psendbuffer(i+1)-1
           !----------------------------------------------------------------
           !  Get the particle id 
           !----------------------------------------------------------------
           ipart = ppm_buffer2part(j)
           ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
           ppm_sendbuffer(ibuffer) = REAL(pdata(ipart),MK)
#elif  __KIND == __DOUBLE_PRECISION
           ppm_sendbuffer(ibuffer) = pdata(ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
           ppm_sendbuffer(ibuffer) = REAL(pdata(ipart),MK)
           ibuffer = ibuffer + 1
           ppm_sendbuffer(ibuffer) = REAL(AIMAG(pdata(ipart)),MK)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
           ppm_sendbuffer(ibuffer) = REAL(pdata(ipart),MK)
           ibuffer = ibuffer + 1
           ppm_sendbuffer(ibuffer) = AIMAG(pdata(ipart))
#elif  __KIND == __INTEGER
           ppm_sendbuffer(ibuffer) = REAL(pdata(ipart),MK)
#elif  __KIND == __LOGICAL
           IF (pdata(ipart)) THEN
              ppm_sendbuffer(ibuffer) = 1.0_MK
           ELSE
              ppm_sendbuffer(ibuffer) = 0.0_MK
           ENDIF
#endif
        ENDDO
#endif
         ENDDO                ! i=1,map%ppm_nsendlist

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
         !  thus icount = 1 and ibuffer = map%ppm_nsendlist (which has not yet been
         !  updated (see below)
         !----------------------------------------------------------------------
         icount  = 0 
         ibuffer = map%ppm_nsendbuffer
        DO i=1,map%ppm_nsendlist
           IF (lda.EQ.2) THEN
              !-------------------------------------------------------------
              !  Unrolled for lda=2
              !-------------------------------------------------------------
              DO j=map%ppm_psendbuffer(i),map%ppm_psendbuffer(i+1)-1
                 !----------------------------------------------------------
                 !  Get the particle id
                 !----------------------------------------------------------
                 ibuffer = ibuffer + 1
                 icount  = icount  + 1
                 ppm_sendbuffer(ibuffer) = ppm_sendbuffer(ibuffer) &
 &                                        + ppm_ghost_offset(icount) 

                 ibuffer = ibuffer + 1
                 icount  = icount  + 1
                 ppm_sendbuffer(ibuffer) = ppm_sendbuffer(ibuffer) &
 &                                        + ppm_ghost_offset(icount) 
              ENDDO
           ELSEIF (lda.EQ.3) THEN
              !-------------------------------------------------------------
              !  Unrolled for lda=3
              !-------------------------------------------------------------
              DO j=map%ppm_psendbuffer(i),map%ppm_psendbuffer(i+1)-1
                 !----------------------------------------------------------
                 !  Get the particle id 
                 !----------------------------------------------------------
                 ibuffer = ibuffer + 1
                 icount  = icount  + 1
                 ppm_sendbuffer(ibuffer) = ppm_sendbuffer(ibuffer) &
 &                                        + ppm_ghost_offset(icount) 

                 ibuffer = ibuffer + 1
                 icount  = icount  + 1
                 ppm_sendbuffer(ibuffer) = ppm_sendbuffer(ibuffer) & 
 &                                        + ppm_ghost_offset(icount) 

                 ibuffer = ibuffer + 1
                 icount  = icount  + 1
                 ppm_sendbuffer(ibuffer) = ppm_sendbuffer(ibuffer) & 
 &                                        + ppm_ghost_offset(icount) 
              ENDDO
           ENDIF ! end of lda = 3
        ENDDO ! enddo of nsendlist
      ENDIF ! end of lpushpp

      !-------------------------------------------------------------------------
      !  Update the particle buffer count
      !-------------------------------------------------------------------------
      map%ppm_nsendbuffer = ibuffer

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
          IF (Pc%maps%vec(mapID)%t%oldNpart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
     &            'Npart must be >=0',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check

#if    __DIM == 1
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE DTYPE(map_part_push_1d)
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE DTYPE(map_part_push_1d)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE DTYPE(map_part_push_1dc)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE DTYPE(map_part_push_1dc)
#elif  __KIND == __INTEGER
      END SUBROUTINE DTYPE(map_part_push_1di)
#elif  __KIND == __LOGICAL
      END SUBROUTINE DTYPE(map_part_push_1dl)
#endif

#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE DTYPE(map_part_push_2d)
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE DTYPE(map_part_push_2d)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE DTYPE(map_part_push_2dc)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE DTYPE(map_part_push_2dc)
#elif  __KIND == __INTEGER
      END SUBROUTINE DTYPE(map_part_push_2di)
#elif  __KIND == __LOGICAL
      END SUBROUTINE DTYPE(map_part_push_2dl)
#endif
#endif
#undef __KIND

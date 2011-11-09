      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_part_modify_push
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
      SUBROUTINE ppm_part_modify_push_1ds(pdata,Npart,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_part_modify_push_1dd(pdata,Npart,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_part_modify_push_1dsc(pdata,Npart,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_part_modify_push_1ddc(pdata,Npart,info)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_part_modify_push_1di(pdata,Npart,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_part_modify_push_1dl(pdata,Npart,info)
#endif
#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_part_modify_push_2ds(pdata,lda,Npart,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_part_modify_push_2dd(pdata,lda,Npart,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_part_modify_push_2dsc(pdata,lda,Npart,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_part_modify_push_2ddc(pdata,lda,Npart,info)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_part_modify_push_2di(pdata,lda,Npart,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_part_modify_push_2dl(pdata,lda,Npart,info)
#endif
#endif
      !!! This routine pushes particle data onto modification buffer.
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
      INTEGER                 , INTENT(IN   )    :: Npart
      !!! Number of particles
      INTEGER                 , INTENT(  OUT)    :: info
      !!! Returns 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3) :: ldu
      INTEGER               :: i,j,k,ipart,ibuffer,icount
      INTEGER               :: iopt,ldb,incr
      REAL(MK)              :: t0
#if   __DIM == 1
      INTEGER, PARAMETER    :: lda = 1
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_part_modify_push',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF
      
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
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_part_modify_push',     &
     &        'buffer dimensions PPM_BUFFER_DIM',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_buffer_type,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_part_modify_push',     &
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
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_part_modify_push',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
#if   __DIM == 2
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_part_modify_push',  &
     &            'lda must be >0 for vector data',__LINE__,info)
              GOTO 8888
          ENDIF
#elif __DIM == 1
          IF (lda .NE. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_part_modify_push',  &
     &            'lda must be =1 for scalar data',__LINE__,info)
              GOTO 8888
          ENDIF
#endif
          IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_part_modify_push',  &
     &            'Npart must be >=0',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
#if   __DIM == 1
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_part_modify_push_1ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_part_modify_push_1dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_part_modify_push_1dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_part_modify_push_1ddc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_part_modify_push_1di
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_part_modify_push_1dl
#endif

#elif __DIM == 2
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_part_modify_push_2ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_part_modify_push_2dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_part_modify_push_2dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_part_modify_push_2ddc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_part_modify_push_2di
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_part_modify_push_2dl
#endif
#endif

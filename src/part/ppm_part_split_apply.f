      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_part_split_apply
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
      SUBROUTINE ppm_part_split_apply_1ds(pdata,pdata_subset,mode,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_part_split_apply_1dd(pdata,pdata_subset,mode,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_part_split_apply_1dsc(pdata,pdata_subset,mode,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_part_split_apply_1ddc(pdata,pdata_subset,mode,info)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_part_split_apply_1di(pdata,pdata_subset,mode,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_part_split_apply_1dl(pdata,pdata_subset,mode,info)
#endif
#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_part_split_apply_2ds(pdata,lda,pdata_subset,mode,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_part_split_apply_2dd(pdata,lda,pdata_subset,mode,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_part_split_apply_2dsc(pdata,lda,pdata_subset,mode,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_part_split_apply_2ddc(pdata,lda,pdata_subset,mode,info)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_part_split_apply_2di(pdata,lda,pdata_subset,mode,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_part_split_apply_2dl(pdata,lda,pdata_subset,mode,info)
#endif
#endif
      !!! This routine inserts new data into a data array following the 
      !!! splitting scheme:
      !!!   new real data goes between data_old(:,Np) and data_old(:,Np+1),
      !!!   new ghost data is appended at the end of data_old.
      !!! This routine should be used after ppm_part_split_compute, which
      !!! determine the splitting scheme, stored in _modify_
      !!! ==============================================================


      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_typedef
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_check_id
      USE ppm_module_util_commopt
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
      INTEGER , DIMENSION(:  ), INTENT(INOUT),POINTER    :: pdata
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:  ), INTENT(INOUT),POINTER    :: pdata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:  ), INTENT(INOUT),POINTER :: pdata
#else
      REAL(MK), DIMENSION(:  ), INTENT(INOUT),POINTER    :: pdata
#endif

#elif __DIM == 2
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:,:), INTENT(INOUT),POINTER    :: pdata
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:), INTENT(INOUT),POINTER    :: pdata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | \
      __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:), INTENT(INOUT),POINTER :: pdata
#else
      REAL(MK), DIMENSION(:,:), INTENT(INOUT),POINTER    :: pdata
#endif
#endif
      !!! The old existing data
#if   __DIM == 2
      INTEGER                 , INTENT(IN   )    :: lda
      !!! The leading dimension of the pdata.
#endif

#if   __DIM == 1
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:  ), INTENT(IN   ),POINTER    :: pdata_subset
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:  ), INTENT(IN   ),POINTER    :: pdata_subset
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:  ), INTENT(IN   ),POINTER :: pdata_subset
#else
      REAL(MK), DIMENSION(:  ), INTENT(IN   ),POINTER    :: pdata_subset
#endif

#elif __DIM == 2
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:,:), INTENT(IN   ),POINTER    :: pdata_subset
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:), INTENT(IN   ),POINTER    :: pdata_subset
#elif __KIND == __SINGLE_PRECISION_COMPLEX | \
      __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:), INTENT(IN   ),POINTER :: pdata_subset
#else
      REAL(MK), DIMENSION(:,:), INTENT(IN   ),POINTER    :: pdata_subset
#endif
#endif
      INTEGER                 , INTENT(IN   ) :: mode
      !!! One of:
      !!! ppm_param_add_real_particles
      !!! ppm_param_add_ghost_particles
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3) :: ldu
      INTEGER               :: i
      INTEGER               :: ipart,iopt
      REAL(MK)              :: t0
      CHARACTER(LEN=ppm_char):: caller='ppm_part_split_apply'
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
      IF (ppm_debug .GT. 1) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Grow pdata array
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow
#if   __DIM == 1
      IF (mode .EQ. ppm_param_add_real_particles) THEN
          ldu(1) = modify%Nrnew
      ELSE IF (mode .EQ. ppm_param_add_ghost_particles) THEN
          ldu(1) = modify%Ngnew
      ENDIF
#elif __DIM ==2
      ldu(1) = lda
      IF (mode .EQ. ppm_param_add_real_particles) THEN
          ldu(2) = modify%Nrnew
      ELSE IF (mode .EQ. ppm_param_add_ghost_particles) THEN
          ldu(2) = modify%Ngnew
      ENDIF
#endif
      CALL ppm_alloc(pdata_subset,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,caller,     &
     &        'pdata_subset',__LINE__,info)
          GOTO 9999
      ENDIF

      IF (mode .EQ. ppm_param_add_real_particles) THEN
          !---------------------------------------------------------------------
          !  Extract real particles
          !---------------------------------------------------------------------
          DO i=1,modify%Nrnew
#if   __DIM == 1
              pdata_subset(i)=pdata(modify%idx_real_new(i))
#elif __DIM ==2
              pdata_subset(1:lda,i)=pdata(1:lda,modify%idx_real_new(i))
#endif
          ENDDO

      ELSE IF (mode .EQ. ppm_param_add_ghost_particles) THEN
          !---------------------------------------------------------------------
          !  Extract ghost particles
          !---------------------------------------------------------------------
          DO i=1,modify%Ngnew
#if   __DIM == 1
              pdata_subset(i)=pdata(modify%idx_ghost_new(i))
#elif __DIM ==2
              pdata_subset(1:lda,i)=pdata(1:lda,modify%idx_ghost_new(i))
#endif
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,caller,  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
#if __DIM ==2
        IF (lda .LT. 1) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,  &
     &          'lda must be >0',__LINE__,info)
            GOTO 8888
        ENDIF
#endif
        IF (.NOT. ASSOCIATED(modify%idx_real_new)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,  &
     &          'unassociated pointer. Call ppm_part_split_compute',__LINE__,info)
            GOTO 8888
        ENDIF
        IF (.NOT. ASSOCIATED(modify%idx_ghost_new)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,  &
     &          'unassociated pointer. Call ppm_part_split_compute',__LINE__,info)
            GOTO 8888
        ENDIF
        IF ((SIZE(modify%idx_real_new).LT.modify%Nrnew)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,  &
     &          'Wrong size for idx_real_new. Call ppm_part_split_compute',&
     &            __LINE__,info)
            GOTO 8888
        ENDIF
        IF ((SIZE(modify%idx_ghost_new).LT.modify%Ngnew)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,  &
     &          'Wrong size for idx_ghost_new. Call ppm_part_split_compute',&
     &            __LINE__,info)
            GOTO 8888
        ENDIF
        IF (modify%Nrnew .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,  &
     &          'modify%Nrnew must be >=0',__LINE__,info)
            GOTO 8888
        ENDIF
        IF (modify%Ngnew .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,  &
     &          'modify%Ngnew must be >=0',__LINE__,info)
            GOTO 8888
        ENDIF
 8888     CONTINUE
      END SUBROUTINE check
#if    __DIM == 1
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_part_split_apply_1ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_part_split_apply_1dd
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_part_split_apply_1dsc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_part_split_apply_1ddc
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_part_split_apply_1di
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_part_split_apply_1dl
#endif
#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_part_split_apply_2ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_part_split_apply_2dd
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_part_split_apply_2dsc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_part_split_apply_2ddc
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_part_split_apply_2di
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_part_split_apply_2dl
#endif
#endif

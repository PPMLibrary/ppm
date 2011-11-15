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
      SUBROUTINE ppm_part_split_apply_1ds(pdata,Npart_old,Mpart_old,&
              pnew,Npart_new,Mpart_new,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_part_split_apply_1dd(pdata,Npart_old,Mpart_old,&
              pnew,Npart_new,Mpart_new,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_part_split_apply_1dsc(pdata,Npart_old,Mpart_old,&
              pnew,Npart_new,Mpart_new,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_part_split_apply_1ddc(pdata,Npart_old,Mpart_old,&
              pnew,Npart_new,Mpart_new,info)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_part_split_apply_1di(pdata,Npart_old,Mpart_old,&
              pnew,Npart_new,Mpart_new,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_part_split_apply_1dl(pdata,Npart_old,Mpart_old,&
              pnew,Npart_new,Mpart_new,info)
#endif
#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_part_split_apply_2ds(pdata,lda,Npart_old,Mpart_old,&
              pnew,Npart_new,Mpart_new,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_part_split_apply_2dd(pdata,lda,Npart_old,Mpart_old,&
              pnew,Npart_new,Mpart_new,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_part_split_apply_2dsc(pdata,lda,Npart_old,Mpart_old,&
              pnew,Npart_new,Mpart_new,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_part_split_apply_2ddc(pdata,lda,Npart_old,Mpart_old,&
              pnew,Npart_new,Mpart_new,info)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_part_split_apply_2di(pdata,lda,Npart_old,Mpart_old,&
              pnew,Npart_new,Mpart_new,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_part_split_apply_2dl(pdata,lda,Npart_old,Mpart_old,&
              pnew,Npart_new,Mpart_new,info)
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

      INTEGER                 , INTENT(IN   ) :: Npart_old
      !!! The number of particles (on the local processor)
      INTEGER                 , INTENT(IN   ) :: Mpart_old
      !!! The number of particles (including ghosts)
#if   __DIM == 1
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:  ), INTENT(IN   ),POINTER    :: pnew
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:  ), INTENT(IN   ),POINTER    :: pnew
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:  ), INTENT(IN   ),POINTER :: pnew
#else
      REAL(MK), DIMENSION(:  ), INTENT(IN   ),POINTER    :: pnew
#endif

#elif __DIM == 2
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:,:), INTENT(IN   ),POINTER    :: pnew
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:), INTENT(IN   ),POINTER    :: pnew
#elif __KIND == __SINGLE_PRECISION_COMPLEX | \
      __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:), INTENT(IN   ),POINTER :: pnew
#else
      REAL(MK), DIMENSION(:,:), INTENT(IN   ),POINTER    :: pnew
#endif
#endif
      INTEGER                 , INTENT(  OUT) :: Npart_new
      !!! The number of particles (on the local processor)
      INTEGER                 , INTENT(  OUT) :: Mpart_new
      !!! The number of particles (including ghosts)
      !!! The new data to be added
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3) :: ldu
      INTEGER               :: i,Nrnew,Ngnew,Nnew
      INTEGER               :: ipart,iopt
      REAL(MK)              :: t0
      LOGICAL               :: valid
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_part_split_apply',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      Nrnew = modify%Nrnew
      Ngnew = modify%Ngnew
      Nnew  = Nrnew + Ngnew


      !-------------------------------------------------------------------------
      !  Grow pdata array
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
#if   __DIM == 1
      ldu(1) = Mpart_old+Nnew
#elif __DIM ==2
      ldu(1) = lda
      ldu(2) = Mpart_old+Nnew
#endif
      CALL ppm_alloc(pdata,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_part_split_apply',     &
     &        'pdata',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Insert pnew into pdata according to split
      !-------------------------------------------------------------------------
      ! shift old ghost particle to leave space for new real ones.
      DO i=Mpart_old,Npart_old+1,-1
#if   __DIM == 1
          pdata(i+Nrnew)=pdata(i)
#elif __DIM ==2
          pdata(1:lda,i+Nrnew)=pdata(1:lda,i)
#endif
      ENDDO
      ! insert new real data
      DO i=1,Nrnew
#if   __DIM == 1
          pdata(Npart_old+i)=pnew(i)
#elif __DIM ==2
          pdata(1:lda,Npart_old+i)=pnew(1:lda,modify%idx_real_new(i))
#endif
      ENDDO
      ! append new ghost data
      DO i=1,Ngnew
#if   __DIM == 1
          pdata(Mpart_old+Nrnew+i)=pnew(modify%idx_ghost_new(i))
#elif __DIM ==2
          pdata(1:lda,Mpart_old+Nrnew+i)=pnew(1:lda,modify%idx_ghost_new(i))
#endif
      ENDDO

      Npart_new = Npart_old + Nrnew
      Mpart_new = Mpart_old + Nrnew + Ngnew


      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_part_split_apply',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_part_split_apply',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
#if __DIM ==2
        IF (lda .LT. 1) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_part_split_apply',  &
     &          'lda must be >0',__LINE__,info)
            GOTO 8888
        ENDIF
#endif
        IF (Npart_old .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_part_split_apply',  &
     &          'Npart_old must be >=0',__LINE__,info)
            GOTO 8888
        ENDIF
        IF (Mpart_old .LT. Npart_old) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_part_split_apply',  &
     &          'Mpart must be >= Npart',__LINE__,info)
            GOTO 8888
        ENDIF
        IF (.NOT. ASSOCIATED(modify%idx_real_new)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_part_split_apply',  &
     &          'unassociated pointer. Call ppm_part_split_compute',__LINE__,info)
            GOTO 8888
        ENDIF
        IF (.NOT. ASSOCIATED(modify%idx_ghost_new)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_part_split_apply',  &
     &          'unassociated pointer. Call ppm_part_split_compute',__LINE__,info)
            GOTO 8888
        ENDIF
        IF ((SIZE(modify%idx_real_new).LT.modify%Nrnew)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_part_split_apply',  &
     &          'Wrong size for idx_real_new. Call ppm_part_split_compute',&
     &            __LINE__,info)
            GOTO 8888
        ENDIF
        IF ((SIZE(modify%idx_ghost_new).LT.modify%Ngnew)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_part_split_apply',  &
     &          'Wrong size for idx_ghost_new. Call ppm_part_split_compute',&
     &            __LINE__,info)
            GOTO 8888
        ENDIF
        IF (modify%Nrnew .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_part_split_apply',  &
     &          'modify%Nrnew must be >=0',__LINE__,info)
            GOTO 8888
        ENDIF
        IF (modify%Ngnew .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_part_split_apply',  &
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

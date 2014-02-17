      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_alloc_3d
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

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE alloc_3d_s(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 3D real single arrays
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE alloc_3d_d(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 3D real double arrays
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE alloc_3d_sc(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 3D complex single arrays
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE alloc_3d_dc(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 3D complex double arrays
#elif __KIND == __INTEGER
      SUBROUTINE alloc_3d_i(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 3D integer arrays
#elif __KIND == __LONGINT
      SUBROUTINE alloc_3d_li(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 3D integer arrays
#elif __KIND == __LOGICAL
      SUBROUTINE alloc_3d_l(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 3D logical arrays
#endif
      !!! (pointers) based on the number of elements.
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
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      REAL(ppm_kind_single)   , DIMENSION(:,:,:), POINTER :: adata
#elif __KIND == __DOUBLE_PRECISION
      REAL(ppm_kind_double)   , DIMENSION(:,:,:), POINTER :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_single), DIMENSION(:,:,:), POINTER :: adata
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_double), DIMENSION(:,:,:), POINTER :: adata
#elif __KIND == __INTEGER
      INTEGER                 , DIMENSION(:,:,:), POINTER :: adata
#elif __KIND == __LONGINT
      INTEGER(ppm_kind_int64) , DIMENSION(:,:,:), POINTER :: adata
#elif __KIND == __LOGICAL
      LOGICAL                 , DIMENSION(:,:,:), POINTER :: adata
#endif
      !!! Pointer to array which is to be (re)allocated
      INTEGER, DIMENSION(:)   , INTENT(IN)    :: lda
      !!! Number of desired elements in all dimensions of array. (>0)
      INTEGER                 , INTENT(IN)    :: iopt
      !!! Allocation mode. One of:
      !!!
      !!! * ppm_param_alloc_fit
      !!! * ppm_param_alloc_fit_preserve
      !!! * ppm_param_alloc_grow
      !!! * ppm_param_alloc_grow_preserve
      !!! * ppm_param_dealloc
      INTEGER                 , INTENT(OUT)   :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      REAL(ppm_kind_single)   , DIMENSION(:,:,:), POINTER :: work
#elif __KIND == __DOUBLE_PRECISION
      REAL(ppm_kind_double)   , DIMENSION(:,:,:), POINTER :: work
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_single), DIMENSION(:,:,:), POINTER :: work
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_double), DIMENSION(:,:,:), POINTER :: work
#elif __KIND == __INTEGER
      INTEGER                 , DIMENSION(:,:,:), POINTER :: work
#elif __KIND == __LONGINT
      INTEGER(ppm_kind_int64) , DIMENSION(:,:,:), POINTER :: work
#elif __KIND == __LOGICAL
      LOGICAL                 , DIMENSION(:,:,:), POINTER :: work
#endif
      INTEGER, DIMENSION(3) :: ldb,ldc,lda_new
      INTEGER               :: i,j,k
      LOGICAL               :: lcopy,lalloc,lrealloc
#ifdef __DEBUG
      REAL(ppm_kind_double) :: t0
#endif
      CHARACTER(LEN=*), PARAMETER :: caller='ppm_alloc_3d'

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
#ifdef __DEBUG
      CALL substart(caller,t0,info)
#else
      info = 0
#endif

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         CALL ppm_alloc_argcheck(caller,iopt,lda,3,info)
         IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Point to proper work array
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      work => work_3ds
#elif __KIND == __DOUBLE_PRECISION
      work => work_3dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      work => work_3dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      work => work_3ddc
#elif __KIND == __INTEGER
      work => work_3di
#elif __KIND == __LONGINT
      work => work_3dli
#elif __KIND == __LOGICAL
      work => work_3dl
#endif

      !-------------------------------------------------------------------------
      !  Check the allocation type
      !-------------------------------------------------------------------------
      lcopy    = .FALSE.
      lalloc   = .FALSE.
      lrealloc = .FALSE.
      SELECT CASE (iopt)
      CASE (ppm_param_alloc_fit_preserve)
         !----------------------------------------------------------------------
         !  fit memory and preserve the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
         ! loop peeling
            ldb(1) = SIZE(adata,1)
            IF (ldb(1).NE.lda(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
               lda_new(1) = lda(1)
            ELSE
               lda_new(1) = ldb(1)
            ENDIF

            ldb(2) = SIZE(adata,2)
            IF (ldb(2).NE.lda(2)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
               lda_new(2) = lda(2)
            ELSE
               lda_new(2) = ldb(2)
            ENDIF

            ldb(3) = SIZE(adata,3)
            IF (ldb(3).NE.lda(3)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
               lda_new(3) = lda(3)
            ELSE
               lda_new(3) = ldb(3)
            ENDIF
         ELSE
            lalloc = .TRUE.
            lda_new(1) = lda(1)
            lda_new(2) = lda(2)
            lda_new(3) = lda(3)
         ENDIF

      CASE (ppm_param_alloc_fit)
         !----------------------------------------------------------------------
         !  fit memory but skip the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
         ! loop peeling
            ldb(1) = SIZE(adata,1)
            IF (ldb(1).NE.lda(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lda_new(1) = lda(1)
            ELSE
               lda_new(1) = ldb(1)
            ENDIF

            ldb(2) = SIZE(adata,2)
            IF (ldb(2).NE.lda(2)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lda_new(2) = lda(2)
            ELSE
               lda_new(2) = ldb(2)
            ENDIF

            ldb(3) = SIZE(adata,3)
            IF (ldb(3).NE.lda(3)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lda_new(3) = lda(3)
            ELSE
               lda_new(3) = ldb(3)
            ENDIF
         ELSE
            lalloc     = .TRUE.
            lda_new(1) = lda(1)
            lda_new(2) = lda(2)
            lda_new(3) = lda(3)
         ENDIF

      CASE (ppm_param_alloc_grow_preserve)
         !----------------------------------------------------------------------
         !  grow memory and preserve the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
         ! loop peeling
            ldb(1) = SIZE(adata,1)
            IF (ldb(1).LT.lda(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
            ENDIF

            ldb(2) = SIZE(adata,2)
            IF (ldb(2).LT.lda(2)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
            ENDIF

            ldb(3) = SIZE(adata,3)
            IF (ldb(3).LT.lda(3)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
            ENDIF
            IF (lrealloc) THEN
               ! no dimension must shrink
               lda_new(1) = MAX(ldb(1),lda(1))
               lda_new(2) = MAX(ldb(2),lda(2))
               lda_new(3) = MAX(ldb(3),lda(3))
            ELSE
               lda_new(1) = ldb(1)
               lda_new(2) = ldb(2)
               lda_new(3) = ldb(3)
            ENDIF
         ELSE
            lalloc     = .TRUE.
            lda_new(1) = lda(1)
            lda_new(2) = lda(2)
            lda_new(3) = lda(3)
         ENDIF

      CASE (ppm_param_alloc_grow)
         !----------------------------------------------------------------------
         !  grow memory but skip the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
         ! loop peeling
            ldb(1) = SIZE(adata,1)
            IF (ldb(1).LT.lda(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
            ENDIF

            ldb(2) = SIZE(adata,2)
            IF (ldb(2).LT.lda(2)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
            ENDIF

            ldb(3) = SIZE(adata,3)
            IF (ldb(3).LT.lda(3)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
            ENDIF
            IF (lrealloc) THEN
               ! no dimension must shrink
               lda_new(1) = MAX(ldb(1),lda(1))
               lda_new(2) = MAX(ldb(2),lda(2))
               lda_new(3) = MAX(ldb(3),lda(3))
            ELSE
               lda_new(1) = ldb(1)
               lda_new(2) = ldb(2)
               lda_new(3) = ldb(3)
            ENDIF
         ELSE
            lalloc     = .TRUE.
            lda_new(1) = lda(1)
            lda_new(2) = lda(2)
            lda_new(3) = lda(3)
         ENDIF

      CASE (ppm_param_dealloc)
         !----------------------------------------------------------------------
         !  deallocate
         !----------------------------------------------------------------------
         IF(ASSOCIATED(adata)) THEN
            DEALLOCATE(adata,STAT=info)
            NULLIFY(adata)
            or_fail_dealloc('DATA')
         ENDIF

      CASE DEFAULT
         !----------------------------------------------------------------------
         !  Unknown iopt
         !----------------------------------------------------------------------
         fail('unknown iopt')

      END SELECT

      !-------------------------------------------------------------------------
      !  Allocate new memory
      !-------------------------------------------------------------------------
      IF (lalloc) THEN
         ALLOCATE(work(lda_new(1),lda_new(2),lda_new(3)),STAT=info)
         or_fail_alloc('WORK',ppm_error=ppm_error_fatal)
      ENDIF

      !-------------------------------------------------------------------------
      !  Make a copy of the old contents
      !-------------------------------------------------------------------------
      IF (lcopy) THEN
         ldc(1) = MIN(lda_new(1),ldb(1))
         ldc(2) = MIN(lda_new(2),ldb(2))
         ldc(3) = MIN(lda_new(3),ldb(3))
         DO k=1,ldc(3)
            DO j=1,ldc(2)
               DO i=1,ldc(1)
                  work(i,j,k) = adata(i,j,k)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  If reallocating, deallocate the old data
      !-------------------------------------------------------------------------
      IF (lrealloc) THEN
         DEALLOCATE(adata,STAT=info)
         or_fail_dealloc('DATA',exit_point=no)
      ENDIF

      !-------------------------------------------------------------------------
      !  Set the pointer to the new array
      !-------------------------------------------------------------------------
      IF (lalloc) THEN
         adata => work
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
#ifdef __DEBUG
      CALL substop(caller,t0,info)
#endif
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE alloc_3d_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE alloc_3d_d
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE alloc_3d_sc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE alloc_3d_dc
#elif __KIND == __INTEGER
      END SUBROUTINE alloc_3d_i
#elif __KIND == __LONGINT
      END SUBROUTINE alloc_3d_li
#elif __KIND == __LOGICAL
      END SUBROUTINE alloc_3d_l
#endif


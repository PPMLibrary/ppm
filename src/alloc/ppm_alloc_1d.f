      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_alloc_1d
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
#if   __LKIND == __LDA64
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE alloc_1d_s_(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 1D real single arrays
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE alloc_1d_d_(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 1D real double arrays
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE alloc_1d_sc_(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 1D complex single arrays
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE alloc_1d_dc_(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 1D complex double arrays
#elif __KIND == __INTEGER
      SUBROUTINE alloc_1d_i_(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 1D integer arrays
#elif __KIND == __LONGINT
      SUBROUTINE alloc_1d_li_(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 1D 64bit integer arrays
#elif __KIND == __LOGICAL
      SUBROUTINE alloc_1d_l_(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 1D logical arrays
#endif
#else
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE alloc_1d_s(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 1D real single arrays
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE alloc_1d_d(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 1D real double arrays
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE alloc_1d_sc(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 1D complex single arrays
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE alloc_1d_dc(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 1D complex double arrays
#elif __KIND == __INTEGER
      SUBROUTINE alloc_1d_i(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 1D integer arrays
#elif __KIND == __LONGINT
      SUBROUTINE alloc_1d_li(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 1D 64bit integer arrays
#elif __KIND == __LOGICAL
      SUBROUTINE alloc_1d_l(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 1D logical arrays
#endif
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
      REAL(ppm_kind_single)     , DIMENSION(:), POINTER :: adata
#elif __KIND == __DOUBLE_PRECISION
      REAL(ppm_kind_double)     , DIMENSION(:), POINTER :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_single)  , DIMENSION(:), POINTER :: adata
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_double)  , DIMENSION(:), POINTER :: adata
#elif __KIND == __INTEGER
      INTEGER                   , DIMENSION(:), POINTER :: adata
#elif __KIND == __LONGINT
      INTEGER(ppm_kind_int64)   , DIMENSION(:), POINTER :: adata
#elif __KIND == __LOGICAL
      LOGICAL                   , DIMENSION(:), POINTER :: adata
#endif
      !!! Pointer to array which is to be (re)allocated.
#if   __LKIND == __LDA64
      INTEGER(ppm_kind_int64), DIMENSION(:), INTENT(IN   ) :: lda
#else
      INTEGER,                 DIMENSION(:), INTENT(IN   ) :: lda
#endif
      !!! Number of desired elements in leading dimension of array. (>0)
      INTEGER                   , INTENT(IN)    :: iopt
      !!! Allocation mode. One of:
      !!!
      !!! * ppm_param_alloc_fit
      !!! * ppm_param_alloc_fit_preserve
      !!! * ppm_param_alloc_grow
      !!! * ppm_param_alloc_grow_preserve
      !!! * ppm_param_dealloc
      INTEGER                   , INTENT(OUT)   :: info
      !!! Returns status, 0 upon success.
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
#ifdef __DEBUG
      REAL(ppm_kind_double) :: t0
#endif

#if   __KIND == __SINGLE_PRECISION
      REAL(ppm_kind_single)   , DIMENSION(:), POINTER :: work
#elif __KIND == __DOUBLE_PRECISION
      REAL(ppm_kind_double)   , DIMENSION(:), POINTER :: work
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_single), DIMENSION(:), POINTER :: work
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_double), DIMENSION(:), POINTER :: work
#elif __KIND == __INTEGER
      INTEGER                 , DIMENSION(:), POINTER :: work
#elif __KIND == __LONGINT
      INTEGER(ppm_kind_int64) , DIMENSION(:), POINTER :: work
#elif __KIND == __LOGICAL
      LOGICAL                 , DIMENSION(:), POINTER :: work
#endif
#if   __LKIND == __LDA64
      INTEGER(ppm_kind_int64) :: i,ldb,ldc
#else
      INTEGER :: i,ldb,ldc
#endif

      LOGICAL :: lcopy,lalloc,lrealloc

      CHARACTER(LEN=*), PARAMETER :: caller='ppm_alloc_1d'
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
#ifdef __DEBUG
      CALL substart(caller,t0,info)
#else
      info = 0
#endif

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL ppm_alloc_argcheck(caller,iopt,lda,1,info)
         IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Point to proper work array
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      work => work_1ds
#elif __KIND == __DOUBLE_PRECISION
      work => work_1dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      work => work_1dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      work => work_1ddc
#elif __KIND == __INTEGER
      work => work_1di
#elif __KIND == __LONGINT
      work => work_1dli
#elif __KIND == __LOGICAL
      work => work_1dl
#endif

      !-------------------------------------------------------------------------
      !  Check the allocation type
      !-------------------------------------------------------------------------
      lcopy    = .FALSE.
      lalloc   = .FALSE.
      lrealloc = .FALSE.
      SELECT CASE(iopt)
      CASE (ppm_param_alloc_fit_preserve)
         !----------------------------------------------------------------------
         !  fit memory and preserve the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
#if   __LKIND == __LDA64
            ldb = SIZE(adata,KIND=ppm_kind_int64)
#else
            ldb = SIZE(adata)
#endif
            IF (ldb.NE.lda(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
            ENDIF
         ELSE
            lalloc = .TRUE.
         ENDIF

      CASE (ppm_param_alloc_fit)
         !----------------------------------------------------------------------
         !  fit memory but skip the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
#if   __LKIND == __LDA64
            ldb = SIZE(adata,KIND=ppm_kind_int64)
#else
            ldb = SIZE(adata)
#endif
            IF (ldb.NE.lda(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
            ENDIF
         ELSE
            lalloc   = .TRUE.
         ENDIF

      CASE (ppm_param_alloc_grow_preserve)
         !----------------------------------------------------------------------
         !  grow memory and preserve the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
#if   __LKIND == __LDA64
            ldb = SIZE(adata,KIND=ppm_kind_int64)
#else
            ldb = SIZE(adata)
#endif
            IF (ldb.LT.lda(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
            ENDIF
         ELSE
            lalloc = .TRUE.
         ENDIF

      CASE (ppm_param_alloc_grow)
         !----------------------------------------------------------------------
         !  grow memory but skip the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
#if   __LKIND == __LDA64
            ldb = SIZE(adata,KIND=ppm_kind_int64)
#else
            ldb = SIZE(adata)
#endif
            IF (ldb.LT.lda(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
            ENDIF
         ELSE
            lalloc = .TRUE.
         ENDIF

      CASE (ppm_param_dealloc)
         !----------------------------------------------------------------------
         !  deallocate
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
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
         ALLOCATE(work(lda(1)),STAT=info)
         or_fail_alloc('WORK',ppm_error=ppm_error_fatal)
      ENDIF

      !-------------------------------------------------------------------------
      !  Make a copy of the present contents
      !-------------------------------------------------------------------------
      IF (lcopy) THEN
         ldc = MIN(ldb,lda(1))
#if   __LKIND == __LDA64
         DO i=1_ppm_kind_int64,ldc
#else
         DO i=1,ldc
#endif
            work(i) = adata(i)
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  If reallocating, deallocate old data
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
#if   __LKIND == __LDA64
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE alloc_1d_s_
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE alloc_1d_d_
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE alloc_1d_sc_
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE alloc_1d_dc_
#elif __KIND == __INTEGER
      END SUBROUTINE alloc_1d_i_
#elif __KIND == __LONGINT
      END SUBROUTINE alloc_1d_li_
#elif __KIND == __LOGICAL
      END SUBROUTINE alloc_1d_l_
#endif
#else
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE alloc_1d_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE alloc_1d_d
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE alloc_1d_sc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE alloc_1d_dc
#elif __KIND == __INTEGER
      END SUBROUTINE alloc_1d_i
#elif __KIND == __LONGINT
      END SUBROUTINE alloc_1d_li
#elif __KIND == __LOGICAL
      END SUBROUTINE alloc_1d_l
#endif
#endif

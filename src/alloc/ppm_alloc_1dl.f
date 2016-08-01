      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_alloc_1dl
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
      SUBROUTINE alloc_1dl_s_(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D real single arrays
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE alloc_1dl_d_(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D real double arrays
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE alloc_1dl_sc_(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D complex single arrays
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE alloc_1dl_dc_(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D complex double arrays
#elif __KIND == __INTEGER
      SUBROUTINE alloc_1dl_i_(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D integer arrays
#elif __KIND == __LONGINT
      SUBROUTINE alloc_1dl_li_(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D 64bit integer arrays
#elif __KIND == __LOGICAL
      SUBROUTINE alloc_1dl_l_(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D logical arrays
#endif
#else
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE alloc_1dl_s(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D real single arrays
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE alloc_1dl_d(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D real double arrays
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE alloc_1dl_sc(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D complex single arrays
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE alloc_1dl_dc(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D complex double arrays
#elif __KIND == __INTEGER
      SUBROUTINE alloc_1dl_i(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D integer arrays
#elif __KIND == __LONGINT
      SUBROUTINE alloc_1dl_li(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D 64bit integer arrays
#elif __KIND == __LOGICAL
      SUBROUTINE alloc_1dl_l(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D logical arrays
#endif
#endif
      !!! (pointers) based on absolute lower and upper index bounds.
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
      REAL(ppm_kind_single)   , DIMENSION(:), POINTER :: adata
#elif __KIND == __DOUBLE_PRECISION
      REAL(ppm_kind_double)   , DIMENSION(:), POINTER :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_single), DIMENSION(:), POINTER :: adata
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_double), DIMENSION(:), POINTER :: adata
#elif __KIND == __INTEGER
      INTEGER                 , DIMENSION(:), POINTER :: adata
#elif __KIND == __LONGINT
      INTEGER(ppm_kind_int64) , DIMENSION(:), POINTER :: adata
#elif __KIND == __LOGICAL
      LOGICAL                 , DIMENSION(:), POINTER :: adata
#endif
      !!! Pointer to array which is to be (re)allocated.
#if   __LKIND == __LDA64
      INTEGER(ppm_kind_int64), DIMENSION(:)   , INTENT(IN)    :: ldl
      !!! Lower index limit in leading dim.
      INTEGER(ppm_kind_int64), DIMENSION(:)   , INTENT(IN)    :: ldu
      !!! Upper index limit in leading dim. (>ldl(1)).
#else
      INTEGER, DIMENSION(:)   , INTENT(IN)    :: ldl
      !!! Lower index limit in leading dim.
      INTEGER, DIMENSION(:)   , INTENT(IN)    :: ldu
      !!! Upper index limit in leading dim. (>ldl(1)).
#endif
      INTEGER                 , INTENT(IN)    :: iopt
      !!! Allocation mode. One of:
      !!! * ppm_param_alloc_fit
      !!! * ppm_param_alloc_fit_preserve
      !!! * ppm_param_alloc_grow
      !!! * ppm_param_alloc_grow_preserve
      !!! * ppm_param_dealloc
      INTEGER                 , INTENT(OUT)   :: info
      !!! Returns status, 0 upon success.
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
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
      INTEGER(ppm_kind_int64) :: i,lda,ldb,ldc,ldd,ldl_new,ldu_new
#else
      INTEGER :: i,lda,ldb,ldc,ldd,ldl_new,ldu_new
#endif
      LOGICAL :: lcopy,lalloc,lrealloc
#ifdef __DEBUG
      REAL(ppm_kind_double) :: t0
#endif
      CHARACTER(LEN=*), PARAMETER :: caller='ppm_alloc_1dl'
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
      IF (ppm_debug.GT.0) THEN
         CALL ppm_alloc_argcheck(caller,iopt,ldl,1,info,ldu)
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
      SELECT CASE (iopt)
      CASE (ppm_param_alloc_fit_preserve)
         !----------------------------------------------------------------------
         !  fit memory and preserve the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
#if   __LKIND == __LDA64
            lda = LBOUND(adata,1,KIND=ppm_kind_int64)
            ldb = UBOUND(adata,1,KIND=ppm_kind_int64)
#else
            lda = LBOUND(adata,1)
            ldb = UBOUND(adata,1)
#endif
            IF (lda.NE.ldl(1).OR.ldb.NE.ldu(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
               ldl_new  = ldl(1)
               ldu_new  = ldu(1)
            ELSE
               ldl_new  = lda
               ldu_new  = ldb
            ENDIF
         ELSE
            lalloc  = .TRUE.
            ldl_new = ldl(1)
            ldu_new = ldu(1)
         ENDIF

      CASE (ppm_param_alloc_fit)
         !----------------------------------------------------------------------
         !  fit memory but skip the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
#if   __LKIND == __LDA64
            lda = LBOUND(adata,1,KIND=ppm_kind_int64)
            ldb = UBOUND(adata,1,KIND=ppm_kind_int64)
#else
            lda = LBOUND(adata,1)
            ldb = UBOUND(adata,1)
#endif
            IF (lda.NE.ldl(1).OR.ldb.NE.ldu(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               ldl_new  = ldl(1)
               ldu_new  = ldu(1)
            ELSE
               ldl_new  = lda
               ldu_new  = ldb
            ENDIF
         ELSE
            lalloc  = .TRUE.
            ldl_new = ldl(1)
            ldu_new = ldu(1)
         ENDIF

      CASE (ppm_param_alloc_grow_preserve)
         !----------------------------------------------------------------------
         !  grow memory and preserve the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
#if   __LKIND == __LDA64
            lda = LBOUND(adata,1,KIND=ppm_kind_int64)
            ldb = UBOUND(adata,1,KIND=ppm_kind_int64)
#else
            lda = LBOUND(adata,1)
            ldb = UBOUND(adata,1)
#endif
            IF (lda.GT.ldl(1).OR.ldb.LT.ldu(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
               ldl_new  = MIN(lda,ldl(1))
               ldu_new  = MAX(ldb,ldu(1))
            ELSE
               ldl_new  = lda
               ldu_new  = ldb
            ENDIF
         ELSE
            lalloc  = .TRUE.
            ldl_new = ldl(1)
            ldu_new = ldu(1)
         ENDIF

      CASE (ppm_param_alloc_grow)
         !----------------------------------------------------------------------
         !  grow memory but skip the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
#if   __LKIND == __LDA64
            lda = LBOUND(adata,1,KIND=ppm_kind_int64)
            ldb = UBOUND(adata,1,KIND=ppm_kind_int64)
#else
            lda = LBOUND(adata,1)
            ldb = UBOUND(adata,1)
#endif
            IF (lda.GT.ldl(1).OR.ldb.LT.ldu(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               ldl_new  = MIN(lda,ldl(1))
               ldu_new  = MAX(ldb,ldu(1))
            ELSE
               ldl_new  = lda
               ldu_new  = ldb
            ENDIF
         ELSE
            lalloc  = .TRUE.
            ldl_new = ldl(1)
            ldu_new = ldu(1)
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
         ALLOCATE(work(ldl_new:ldu_new),STAT=info)
         or_fail_alloc('WORK',ppm_error=ppm_error_fatal)
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate new memory
      !-------------------------------------------------------------------------
      IF (lcopy) THEN
         ldc = MAX(lda,ldl_new)
         ldd = MIN(ldb,ldu_new)
         DO i=ldc,ldd
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
      END SUBROUTINE alloc_1dl_s_
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE alloc_1dl_d_
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE alloc_1dl_sc_
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE alloc_1dl_dc_
#elif __KIND == __INTEGER
      END SUBROUTINE alloc_1dl_i_
#elif __KIND == __LONGINT
      END SUBROUTINE alloc_1dl_li_
#elif __KIND == __LOGICAL
      END SUBROUTINE alloc_1dl_l_
#endif
#else
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE alloc_1dl_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE alloc_1dl_d
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE alloc_1dl_sc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE alloc_1dl_dc
#elif __KIND == __INTEGER
      END SUBROUTINE alloc_1dl_i
#elif __KIND == __LONGINT
      END SUBROUTINE alloc_1dl_li
#elif __KIND == __LOGICAL
      END SUBROUTINE alloc_1dl_l
#endif
#endif

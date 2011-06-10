      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_alloc_2dl
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

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE alloc_2dl_s(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 2D real single arrays
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE alloc_2dl_d(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 2D real double arrays
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE alloc_2dl_sc(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 2D complex single arrays
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE alloc_2dl_dc(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 2D complex double arrays
#elif __KIND == __INTEGER
      SUBROUTINE alloc_2dl_i(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 2D integer arrays
#elif __KIND == __LONGINT
      SUBROUTINE alloc_2dl_li(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 2D 64bit integer arrays
#elif __KIND == __LOGICAL
      SUBROUTINE alloc_2dl_l(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 2D logical arrays
#endif
      !!! (pointers) based on absolute lower and upper index bounds.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

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
      REAL(ppm_kind_single)   , DIMENSION(:,:), POINTER :: adata
#elif __KIND == __DOUBLE_PRECISION
      REAL(ppm_kind_double)   , DIMENSION(:,:), POINTER :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_single), DIMENSION(:,:), POINTER :: adata
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_double), DIMENSION(:,:), POINTER :: adata
#elif __KIND == __INTEGER
      INTEGER                 , DIMENSION(:,:), POINTER :: adata
#elif __KIND == __LONGINT
      INTEGER(ppm_kind_int64) , DIMENSION(:,:), POINTER :: adata
#elif __KIND == __LOGICAL
      LOGICAL                 , DIMENSION(:,:), POINTER :: adata
#endif
      !!! Pointer to array which is to be (re)allocated.
      INTEGER, DIMENSION(:)   , INTENT(IN)    :: ldl
      !!! Lower index limits
      INTEGER, DIMENSION(:)   , INTENT(IN)    :: ldu
      !!! Upper index limits (>ldl(1:2)).
      INTEGER                 , INTENT(IN)    :: iopt
      !!! allocation mode. One of:
      !!!
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
      REAL(ppm_kind_single)   , DIMENSION(:,:), POINTER :: work
#elif __KIND == __DOUBLE_PRECISION
      REAL(ppm_kind_double)   , DIMENSION(:,:), POINTER :: work
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_single), DIMENSION(:,:), POINTER :: work
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_double), DIMENSION(:,:), POINTER :: work
#elif __KIND == __INTEGER
      INTEGER                 , DIMENSION(:,:), POINTER :: work
#elif __KIND == __LONGINT
      INTEGER(ppm_kind_int64) , DIMENSION(:,:), POINTER :: work
#elif __KIND == __LOGICAL
      LOGICAL                 , DIMENSION(:,:), POINTER :: work
#endif
      INTEGER, DIMENSION(2) :: lda,ldb,ldc,ldd,ldl_new,ldu_new
      INTEGER               :: i,j
      LOGICAL               :: lcopy,lalloc,lrealloc
      REAL(ppm_kind_double) :: t0
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      IF (ppm_debug.GE.3) THEN
          CALL substart('ppm_alloc_2dl',t0,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         CALL ppm_alloc_argcheck('ppm_alloc_2dl',iopt,ldl,2,info,ldu)
         IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Point to proper work array
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      work => work_2ds
#elif __KIND == __DOUBLE_PRECISION
      work => work_2dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      work => work_2dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      work => work_2ddc
#elif __KIND == __INTEGER
      work => work_2di
#elif __KIND == __LONGINT
      work => work_2dli
#elif __KIND == __LOGICAL
      work => work_2dl
#endif

      !-------------------------------------------------------------------------
      !  Check the allocation type
      !-------------------------------------------------------------------------
      lcopy    = .FALSE.
      lalloc   = .FALSE.
      lrealloc = .FALSE.
      IF     (iopt.EQ.ppm_param_alloc_fit_preserve) THEN
         !----------------------------------------------------------------------
         !  fit memory and preserve the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
         ! loop peeling
            lda(1) = LBOUND(adata,1)
            ldb(1) = UBOUND(adata,1)
            IF (lda(1).NE.ldl(1).OR.ldb(1).NE.ldu(1)) THEN
              lrealloc = .TRUE.
              lalloc   = .TRUE.
              lcopy    = .TRUE.
              ldl_new(1) = ldl(1)
              ldu_new(1) = ldu(1)
            ELSE
              ldl_new(1) = lda(1)
              ldu_new(1) = ldb(1)
            ENDIF

            lda(2) = LBOUND(adata,2)
            ldb(2) = UBOUND(adata,2)
            IF (lda(2).NE.ldl(2).OR.ldb(2).NE.ldu(2)) THEN
              lrealloc = .TRUE.
              lalloc   = .TRUE.
              lcopy    = .TRUE.
              ldl_new(2) = ldl(2)
              ldu_new(2) = ldu(2)
            ELSE
              ldl_new(2) = lda(2)
              ldu_new(2) = ldb(2)
            ENDIF
         ELSE
            lalloc = .TRUE.
            ldl_new(1) = ldl(1)
            ldl_new(2) = ldl(2)

            ldu_new(1) = ldu(1)
            ldu_new(2) = ldu(2)
         ENDIF
      ELSEIF (iopt.EQ.ppm_param_alloc_fit) THEN
         !----------------------------------------------------------------------
         !  fit memory but skip the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
         ! loop peeling
            lda(1) = LBOUND(adata,1)
            ldb(1) = UBOUND(adata,1)
            IF (lda(1).NE.ldl(1).OR.ldb(1).NE.ldu(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               ldl_new(1) = ldl(1)
               ldu_new(1) = ldu(1)
            ELSE
               ldl_new(1) = lda(1)
               ldu_new(1) = ldb(1)
            ENDIF

            lda(2) = LBOUND(adata,2)
            ldb(2) = UBOUND(adata,2)
            IF (lda(2).NE.ldl(2).OR.ldb(2).NE.ldu(2)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               ldl_new(2) = ldl(2)
               ldu_new(2) = ldu(2)
            ELSE
               ldl_new(2) = lda(2)
               ldu_new(2) = ldb(2)
            ENDIF
         ELSE
            lalloc     = .TRUE.
            ldl_new(1) = ldl(1)
            ldl_new(2) = ldl(2)

            ldu_new(1) = ldu(1)
            ldu_new(2) = ldu(2)
         ENDIF
      ELSEIF (iopt.EQ.ppm_param_alloc_grow_preserve) THEN
         !----------------------------------------------------------------------
         !  grow memory and preserve the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
         ! loop peeling
            lda(1) = LBOUND(adata,1)
            ldb(1) = UBOUND(adata,1)
            IF (lda(1).GT.ldl(1).OR.ldb(1).LT.ldu(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
            ENDIF

            lda(2) = LBOUND(adata,2)
            ldb(2) = UBOUND(adata,2)
            IF (lda(2).GT.ldl(2).OR.ldb(2).LT.ldu(2)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
            ENDIF
            IF (lrealloc) THEN
               ! no dimension must shrink
               ldu_new(1) = MAX(ldb(1),ldu(1))
               ldu_new(2) = MAX(ldb(2),ldu(2))

               ldl_new(1) = MIN(lda(1),ldl(1))
               ldl_new(2) = MIN(lda(2),ldl(2))
            ELSE
               ldl_new(1) = lda(1)
               ldl_new(2) = lda(2)

               ldu_new(1) = ldb(1)
               ldu_new(2) = ldb(2)
            ENDIF
         ELSE
            lalloc     = .TRUE.
            ldl_new(1) = ldl(1)
            ldl_new(2) = ldl(2)

            ldu_new(1) = ldu(1)
            ldu_new(2) = ldu(2)
         ENDIF
      ELSEIF (iopt.EQ.ppm_param_alloc_grow) THEN
         !----------------------------------------------------------------------
         !  grow memory but skip the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
         ! loop peeling
            lda(1) = LBOUND(adata,1)
            ldb(1) = UBOUND(adata,1)
            IF (lda(1).GT.ldl(1).OR.ldb(1).LT.ldu(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
            ENDIF

            lda(2) = LBOUND(adata,2)
            ldb(2) = UBOUND(adata,2)
            IF (lda(2).GT.ldl(2).OR.ldb(2).LT.ldu(2)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
            ENDIF
            IF (lrealloc) THEN
               ! no dimension must shrink
               ldu_new(1) = MAX(ldb(1),ldu(1))
               ldu_new(2) = MAX(ldb(2),ldu(2))

               ldl_new(1) = MIN(lda(1),ldl(1))
               ldl_new(2) = MIN(lda(2),ldl(2))
            ELSE
               ldl_new(1) = lda(1)
               ldl_new(2) = lda(2)

               ldu_new(1) = ldb(1)
               ldu_new(2) = ldb(2)
            ENDIF
         ELSE
            lalloc     = .TRUE.
            ldl_new(1) = ldl(1)
            ldl_new(2) = ldl(2)

            ldu_new(1) = ldu(1)
            ldu_new(2) = ldu(2)
         ENDIF
      ELSEIF (iopt.EQ.ppm_param_dealloc) THEN
         !----------------------------------------------------------------------
         !  deallocate
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
            DEALLOCATE(adata,STAT=info)
            NULLIFY(adata)
            IF (info .NE. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_dealloc,'ppm_alloc_2dl',   &
     &             'DATA',__LINE__,info)
            ENDIF
         ENDIF
      ELSE
         !----------------------------------------------------------------------
         !  Unknown iopt
         !----------------------------------------------------------------------
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_alloc_2dl',                      &
     &                  'unknown iopt',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate new memory
      !-------------------------------------------------------------------------
      IF (lalloc) THEN
         ALLOCATE(work(ldl_new(1):ldu_new(1),ldl_new(2):ldu_new(2)),STAT=info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_alloc_2dl',   &
     &           'WORK',__LINE__,info)
             GOTO 9999
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Make a copy of the present contents
      !-------------------------------------------------------------------------
      IF (lcopy) THEN
         ldc(1) = MAX(ldl_new(1),lda(1))
         ldc(2) = MAX(ldl_new(2),lda(2))

         ldd(1) = MIN(ldu_new(1),ldb(1))
         ldd(2) = MIN(ldu_new(2),ldb(2))
         DO j=ldc(2),ldd(2)
            DO i=ldc(1),ldd(1)
               work(i,j) = adata(i,j)
            ENDDO
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  If reallocating, deallocate old data
      !-------------------------------------------------------------------------
      IF (lrealloc) THEN
         DEALLOCATE(adata,STAT=info)
         !NULLIFY(adata)
         IF (info .NE. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_dealloc,'ppm_alloc_2dl',   &
     &           'DATA',__LINE__,info)
         ENDIF
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
      IF (ppm_debug.GE.3) THEN
          CALL substop('ppm_alloc_2dl',t0,info)
      ENDIF
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE alloc_2dl_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE alloc_2dl_d
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE alloc_2dl_sc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE alloc_2dl_dc
#elif __KIND == __INTEGER
      END SUBROUTINE alloc_2dl_i
#elif __KIND == __LONGINT
      END SUBROUTINE alloc_2dl_li
#elif __KIND == __LOGICAL
      END SUBROUTINE alloc_2dl_l
#endif

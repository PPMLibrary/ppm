      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_alloc_4dl
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Lab (ETH Zurich), 
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
      SUBROUTINE ppm_alloc_4dls(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 4D real single arrays
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_alloc_4dld(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 4D real double arrays
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_alloc_4dlsc(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 4D complex single arrays
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_alloc_4dldc(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 4D complex double arrays
#elif __KIND == __INTEGER
      SUBROUTINE ppm_alloc_4dli(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 4D integer arrays
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_alloc_4dll(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 4D logical arrays
#endif
      !!! (pointers) based on absolute lower and upper index bounds.
      !!!
      !!! [NOTE]
      !!! if adata is already associated the memory requirement
      !!! of this routine is 2*size(adata)
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
      REAL(ppm_kind_single)   , DIMENSION(:,:,:,:), POINTER :: adata
#elif __KIND == __DOUBLE_PRECISION
      REAL(ppm_kind_double)   , DIMENSION(:,:,:,:), POINTER :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_single), DIMENSION(:,:,:,:), POINTER :: adata
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_double), DIMENSION(:,:,:,:), POINTER :: adata
#elif __KIND == __INTEGER
      INTEGER                 , DIMENSION(:,:,:,:), POINTER :: adata
#elif __KIND == __LOGICAL
      LOGICAL                 , DIMENSION(:,:,:,:), POINTER :: adata
#endif
      !!! Pointer to array which is to be (re)allocated.
      INTEGER, DIMENSION(:)   , INTENT(IN)    :: ldl
      !!! Lower index limits
      INTEGER, DIMENSION(:)   , INTENT(IN)    :: ldu
      !!! Upper index limits (>ldl(1:4)).
      INTEGER                 , INTENT(IN)    :: iopt
      !!! Allocation mode. One of:
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
      REAL(ppm_kind_single)   , DIMENSION(:,:,:,:), POINTER :: work
#elif __KIND == __DOUBLE_PRECISION
      REAL(ppm_kind_double)   , DIMENSION(:,:,:,:), POINTER :: work
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_single), DIMENSION(:,:,:,:), POINTER :: work
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_double), DIMENSION(:,:,:,:), POINTER :: work
#elif __KIND == __INTEGER
      INTEGER                 , DIMENSION(:,:,:,:), POINTER :: work
#elif __KIND == __LOGICAL
      LOGICAL                 , DIMENSION(:,:,:,:), POINTER :: work
#endif
      INTEGER, DIMENSION(4) :: lda,ldb,ldc,ldd,ldl_new,ldu_new
      INTEGER               :: i,j,k,l
      LOGICAL               :: lcopy,lalloc,lrealloc
      REAL(ppm_kind_double) :: t0
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_alloc_4dl',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         CALL ppm_alloc_argcheck('ppm_alloc_4dl',iopt,ldl,4,info,ldu)
         IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Point to proper work array
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      work => work_4ds
#elif __KIND == __DOUBLE_PRECISION
      work => work_4dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      work => work_4dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      work => work_4ddc
#elif __KIND == __INTEGER
      work => work_4di
#elif __KIND == __LOGICAL
      work => work_4dl
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

            lda(3) = LBOUND(adata,3)
            ldb(3) = UBOUND(adata,3)
            IF (lda(3).NE.ldl(3).OR.ldb(3).NE.ldu(3)) THEN
              lrealloc = .TRUE.
              lalloc   = .TRUE.
              lcopy    = .TRUE.
              ldl_new(3) = ldl(3)
              ldu_new(3) = ldu(3)
            ELSE
              ldl_new(3) = lda(3)
              ldu_new(3) = ldb(3)
            ENDIF

            lda(4) = LBOUND(adata,4)
            ldb(4) = UBOUND(adata,4)
            IF (lda(4).NE.ldl(4).OR.ldb(4).NE.ldu(4)) THEN
              lrealloc = .TRUE.
              lalloc   = .TRUE.
              lcopy    = .TRUE.
              ldl_new(4) = ldl(4)
              ldu_new(4) = ldu(4)
            ELSE
              ldl_new(4) = lda(4)
              ldu_new(4) = ldb(4)
            ENDIF
         ELSE
            lalloc     = .TRUE.
            ldl_new(1) = ldl(1)
            ldl_new(2) = ldl(2)
            ldl_new(3) = ldl(3)
            ldl_new(4) = ldl(4)

            ldu_new(1) = ldu(1)
            ldu_new(2) = ldu(2)
            ldu_new(3) = ldu(3)
            ldu_new(4) = ldu(4)
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

            lda(3) = LBOUND(adata,3)
            ldb(3) = UBOUND(adata,3)
            IF (lda(3).NE.ldl(3).OR.ldb(3).NE.ldu(3)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               ldl_new(3) = ldl(3)
               ldu_new(3) = ldu(3)
            ELSE
               ldl_new(3) = lda(3)
               ldu_new(3) = ldb(3)
            ENDIF

            lda(4) = LBOUND(adata,4)
            ldb(4) = UBOUND(adata,4)
            IF (lda(4).NE.ldl(4).OR.ldb(4).NE.ldu(4)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               ldl_new(4) = ldl(4)
               ldu_new(4) = ldu(4)
            ELSE
               ldl_new(4) = lda(4)
               ldu_new(4) = ldb(4)
            ENDIF
         ELSE
            lalloc     = .TRUE.
            ldl_new(1) = ldl(1)
            ldl_new(2) = ldl(2)
            ldl_new(3) = ldl(3)
            ldl_new(4) = ldl(4)

            ldu_new(1) = ldu(1)
            ldu_new(2) = ldu(2)
            ldu_new(3) = ldu(3)
            ldu_new(4) = ldu(4)
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

            lda(3) = LBOUND(adata,3)
            ldb(3) = UBOUND(adata,3)
            IF (lda(3).GT.ldl(3).OR.ldb(3).LT.ldu(3)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
            ENDIF

            lda(4) = LBOUND(adata,4)
            ldb(4) = UBOUND(adata,4)
            IF (lda(4).GT.ldl(4).OR.ldb(4).LT.ldu(4)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
            ENDIF
            IF (lrealloc) THEN
               ! no dimension must shrink
               ldu_new(1) = MAX(ldb(1),ldu(1))
               ldu_new(2) = MAX(ldb(2),ldu(2))
               ldu_new(3) = MAX(ldb(3),ldu(3))
               ldu_new(4) = MAX(ldb(4),ldu(4))

               ldl_new(1) = MIN(lda(1),ldl(1))
               ldl_new(2) = MIN(lda(2),ldl(2))
               ldl_new(3) = MIN(lda(3),ldl(3))
               ldl_new(4) = MIN(lda(4),ldl(4))
            ELSE
               ldl_new(1) = lda(1)
               ldl_new(2) = lda(2)
               ldl_new(3) = lda(3)
               ldl_new(4) = lda(4)

               ldu_new(1) = ldb(1)
               ldu_new(2) = ldb(2)
               ldu_new(3) = ldb(3)
               ldu_new(4) = ldb(4)
            ENDIF
         ELSE
            lalloc = .TRUE.
            ldl_new(1) = ldl(1)
            ldl_new(2) = ldl(2)
            ldl_new(3) = ldl(3)
            ldl_new(4) = ldl(4)

            ldu_new(1) = ldu(1)
            ldu_new(2) = ldu(2)
            ldu_new(3) = ldu(3)
            ldu_new(4) = ldu(4)
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

            lda(3) = LBOUND(adata,3)
            ldb(3) = UBOUND(adata,3)
            IF (lda(3).GT.ldl(3).OR.ldb(3).LT.ldu(3)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
            ENDIF

            lda(4) = LBOUND(adata,4)
            ldb(4) = UBOUND(adata,4)
            IF (lda(4).GT.ldl(4).OR.ldb(4).LT.ldu(4)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
            ENDIF
            IF (lrealloc) THEN
               ! no dimension must shrink
               ldu_new(1) = MAX(ldb(1),ldu(1))
               ldu_new(2) = MAX(ldb(2),ldu(2))
               ldu_new(3) = MAX(ldb(3),ldu(3))
               ldu_new(4) = MAX(ldb(4),ldu(4))

               ldl_new(1) = MIN(lda(1),ldl(1))
               ldl_new(2) = MIN(lda(2),ldl(2))
               ldl_new(3) = MIN(lda(3),ldl(3))
               ldl_new(4) = MIN(lda(4),ldl(4))
            ELSE
               ldl_new(1) = lda(1)
               ldl_new(2) = lda(2)
               ldl_new(3) = lda(3)
               ldl_new(4) = lda(4)

               ldu_new(1) = ldb(1)
               ldu_new(2) = ldb(2)
               ldu_new(3) = ldb(3)
               ldu_new(4) = ldb(4)
            ENDIF
         ELSE
            lalloc = .TRUE.
            ldl_new(1) = ldl(1)
            ldl_new(2) = ldl(2)
            ldl_new(3) = ldl(3)
            ldl_new(4) = ldl(4)

            ldu_new(1) = ldu(1)
            ldu_new(2) = ldu(2)
            ldu_new(3) = ldu(3)
            ldu_new(4) = ldu(4)
         ENDIF
      ELSEIF (iopt.EQ.ppm_param_dealloc) THEN
         !----------------------------------------------------------------------
         !  deallocate
         !----------------------------------------------------------------------
         IF(ASSOCIATED(adata)) THEN
            DEALLOCATE(adata,STAT=info)
            NULLIFY(adata)
            IF (info .NE. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_dealloc,'ppm_alloc_4dl',   &
     &             'DATA',__LINE__,info)
            ENDIF
         ENDIF
      ELSE
         !----------------------------------------------------------------------
         !  Unknown iopt
         !----------------------------------------------------------------------
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_alloc_4dl',                      &
     &                  'unknown iopt',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate new memory
      !-------------------------------------------------------------------------
      IF (lalloc) THEN
         ALLOCATE(work(ldl_new(1):ldu_new(1),                                  &
     &                 ldl_new(2):ldu_new(2),                                  &
     &                 ldl_new(3):ldu_new(3),                                  &
     &                 ldl_new(4):ldu_new(4)),STAT=info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_alloc_4dl',   &
     &           'WORK',__LINE__,info)
             GOTO 9999
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Make a copy of the present contents
      !-------------------------------------------------------------------------
      IF (lcopy) THEN
         DO k=1,4
            ldc(k) = MAX(ldl_new(k),lda(k))
            ldd(k) = MIN(ldu_new(k),ldb(k))
         ENDDO
         DO l=ldc(4),ldd(4)
            DO k=ldc(3),ldd(3)
               DO j=ldc(2),ldd(2)
                  DO i=ldc(1),ldd(1)
                     work(i,j,k,l) = adata(i,j,k,l)
                  ENDDO
               ENDDO
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
             CALL ppm_error(ppm_err_dealloc,'ppm_alloc_4dl',   &
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
      CALL substop('ppm_alloc_4dl',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_alloc_4dls
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_alloc_4dld
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_alloc_4dlsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_alloc_4dldc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_alloc_4dli
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_alloc_4dll
#endif


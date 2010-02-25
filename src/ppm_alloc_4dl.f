      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_alloc_4dl
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine (re)allocates the memory of
      !                 four-dimensional arrays (pointers) based on absolute
      !                 lower and upper index bounds.
      !
      !  Input        : ldl(1:4) (I) Lower index limits
      !                 ldu(1:4) (I) Upper index limits
      !                              (>ldl(1:4)).
      !                 iopt     (I) allocation mode. One of:
      !                                  ppm_param_alloc_fit
      !                                  ppm_param_alloc_fit_preserve
      !                                  ppm_param_alloc_grow
      !                                  ppm_param_alloc_grow_preserve
      !                                  ppm_param_dealloc
      !
      !  Input/output : adata    (P) Pointer to array which is to be
      !                              (re)allocated.
      !
      !  Output       : info     (I) 0 upon success.
      !
      !  Remarks      : if adata is already associated the memory requirement
      !                 of this routine is 2*size(adata)
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_alloc_4dl.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.8  2006/09/04 18:34:40  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.7  2004/11/11 15:23:01  ivos
      !  Moved allocatable work data to the module.
      !
      !  Revision 1.6  2004/10/01 16:33:31  ivos
      !  cosmetics.
      !
      !  Revision 1.5  2004/10/01 16:08:55  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.4  2004/07/26 07:45:22  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.3  2004/07/19 15:47:02  ivos
      !  bugfix: grow and grow_preserve could shrink non-growing dimensions
      !  since they would only pay attention to the growing one.
      !  Fixed by introducing lda_new,ldl_new,ldu_new.
      !
      !  Revision 1.2  2004/04/06 12:39:14  ivos
      !  Wrapped a line.
      !
      !  Revision 1.1  2004/03/22 12:29:26  kotsalie
      !  I need this routine for multigrid implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_alloc_4dls(adata,ldl,ldu,iopt,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_alloc_4dld(adata,ldl,ldu,iopt,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_alloc_4dlsc(adata,ldl,ldu,iopt,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_alloc_4dldc(adata,ldl,ldu,iopt,info)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_alloc_4dli(adata,ldl,ldu,iopt,info)
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_alloc_4dll(adata,ldl,ldu,iopt,info)
#endif
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
      INTEGER, DIMENSION(:)   , INTENT(IN)    :: ldl,ldu
      INTEGER                 , INTENT(IN)    :: iopt
      INTEGER                 , INTENT(OUT)   :: info
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
         !----------------------------------------------------------------------
         !  Check the iopt
         !----------------------------------------------------------------------
         IF (iopt.NE.ppm_param_alloc_fit          .AND.                        &
     &       iopt.NE.ppm_param_alloc_fit_preserve .AND.                        &
     &       iopt.NE.ppm_param_alloc_grow         .AND.                        &
     &       iopt.NE.ppm_param_alloc_grow_preserve.AND.                        &
     &       iopt.NE.ppm_param_dealloc) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_alloc_4dl',                   &
     &         'unknown iopt',__LINE__,info)
            GOTO 9999
         ENDIF

         !----------------------------------------------------------------------
         !  Check the lda (must be greater than zero)
         !----------------------------------------------------------------------
         IF (iopt.NE.ppm_param_dealloc) THEN
            IF (ldl(1) .GT. ldu(1)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_alloc_4dl',  &
     &              'ldu(1) must be >= ldl(1)',__LINE__,info)
                GOTO 9999
            ENDIF
            IF (ldl(2) .GT. ldu(2)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_alloc_4dl',  &
     &              'ldu(2) must be >= ldl(2)',__LINE__,info)
                GOTO 9999
            ENDIF
            IF (ldl(3) .GT. ldu(3)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_alloc_4dl',  &
     &              'ldu(3) must be >= ldl(3)',__LINE__,info)
                GOTO 9999
            ENDIF
            IF (ldl(4) .GT. ldu(4)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_alloc_4dl',  &
     &              'ldu(4) must be >= ldl(4)',__LINE__,info)
                GOTO 9999
            ENDIF
         ENDIF
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
            DO k=1,4
               lda(k) = LBOUND(adata,k)
               ldb(k) = UBOUND(adata,k)
               IF (lda(k).NE.ldl(k).OR.ldb(k).NE.ldu(k)) THEN
                  lrealloc = .TRUE.
                  lalloc   = .TRUE.
                  lcopy    = .TRUE.
                  ldl_new(k) = ldl(k)
                  ldu_new(k) = ldu(k)
               ELSE
                  ldl_new(k) = lda(k)
                  ldu_new(k) = ldb(k)
               ENDIF            
            ENDDO
         ELSE
            lalloc = .TRUE.
            DO k=1,4
               ldl_new(k) = ldl(k)
               ldu_new(k) = ldu(k)
            ENDDO
         ENDIF
      ELSEIF (iopt.EQ.ppm_param_alloc_fit) THEN
         !----------------------------------------------------------------------
         !  fit memory but skip the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
            DO k=1,4
               lda(k) = LBOUND(adata,k)
               ldb(k) = UBOUND(adata,k)
               IF (lda(k).NE.ldl(k).OR.ldb(k).NE.ldu(k)) THEN
                  lrealloc = .TRUE. 
                  lalloc   = .TRUE. 
                  ldl_new(k) = ldl(k)
                  ldu_new(k) = ldu(k)
               ELSE
                  ldl_new(k) = lda(k)
                  ldu_new(k) = ldb(k)
               ENDIF            
            ENDDO
         ELSE
            lalloc = .TRUE. 
            DO k=1,4
               ldl_new(k) = ldl(k)
               ldu_new(k) = ldu(k)
            ENDDO
         ENDIF
      ELSEIF (iopt.EQ.ppm_param_alloc_grow_preserve) THEN
         !----------------------------------------------------------------------
         !  grow memory and preserve the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
            DO k=1,4
               lda(k) = LBOUND(adata,k)
               ldb(k) = UBOUND(adata,k)
               IF (lda(k).GT.ldl(k).OR.ldb(k).LT.ldu(k)) THEN
                  lrealloc = .TRUE.
                  lalloc   = .TRUE.
                  lcopy    = .TRUE.
               ENDIF            
            ENDDO
            IF (lrealloc) THEN
               ! no dimension must shrink
               DO k=1,4
                  ldu_new(k) = MAX(ldb(k),ldu(k))
                  ldl_new(k) = MIN(lda(k),ldl(k))
               ENDDO
            ELSE
               DO k=1,4
                  ldl_new(k) = lda(k)
                  ldu_new(k) = ldb(k)
               ENDDO
            ENDIF
         ELSE
            lalloc = .TRUE.
            DO k=1,4
               ldl_new(k) = ldl(k)
               ldu_new(k) = ldu(k)
            ENDDO
         ENDIF
      ELSEIF (iopt.EQ.ppm_param_alloc_grow) THEN
         !----------------------------------------------------------------------
         !  grow memory but skip the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
            DO k=1,4
               lda(k) = LBOUND(adata,k)
               ldb(k) = UBOUND(adata,k)
               IF (lda(k).GT.ldl(k).OR.ldb(k).LT.ldu(k)) THEN
                  lrealloc = .TRUE.
                  lalloc   = .TRUE.
               ENDIF            
            ENDDO
            IF (lrealloc) THEN
               ! no dimension must shrink
               DO k=1,4
                  ldu_new(k) = MAX(ldb(k),ldu(k))
                  ldl_new(k) = MIN(lda(k),ldl(k))
               ENDDO
            ELSE
               DO k=1,4
                  ldl_new(k) = lda(k)
                  ldu_new(k) = ldb(k)
               ENDDO
            ENDIF
         ELSE
            lalloc = .TRUE.
            DO k=1,4
               ldl_new(k) = ldl(k)
               ldu_new(k) = ldu(k)
            ENDDO
         ENDIF
      ELSEIF (iopt.EQ.ppm_param_dealloc) THEN
         !----------------------------------------------------------------------
         !  deallocate
         !----------------------------------------------------------------------
         IF(ASSOCIATED(adata)) DEALLOCATE(adata,STAT=info)
         IF (info .NE. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_dealloc,'ppm_alloc_4dl',   &
     &           'DATA',__LINE__,info)
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


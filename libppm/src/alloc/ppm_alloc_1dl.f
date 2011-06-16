      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_alloc_1dl
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_alloc_1dls(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D real single arrays
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_alloc_1dld(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D real double arrays
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_alloc_1dlsc(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D complex single arrays
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_alloc_1dldc(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D complex double arrays
#elif __KIND == __INTEGER
      SUBROUTINE ppm_alloc_1dli(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D integer arrays
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_alloc_1dll(adata,ldl,ldu,iopt,info)
      !!! (Re)allocates the memory of 1D logical arrays
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
      REAL(ppm_kind_single)   , DIMENSION(:), POINTER :: adata
#elif __KIND == __DOUBLE_PRECISION
      REAL(ppm_kind_double)   , DIMENSION(:), POINTER :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_single), DIMENSION(:), POINTER :: adata
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_double), DIMENSION(:), POINTER :: adata
#elif __KIND == __INTEGER
      INTEGER                 , DIMENSION(:), POINTER :: adata
#elif __KIND == __LOGICAL
      LOGICAL                 , DIMENSION(:), POINTER :: adata
#endif
      !!! Pointer to array which is to be (re)allocated.
      INTEGER, DIMENSION(:)   , INTENT(IN)    :: ldl
      !!! Lower index limit in leading dim.
      INTEGER, DIMENSION(:)   , INTENT(IN)    :: ldu
      !!! Upper index limit in leading dim. (>ldl(1)).
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
#elif __KIND == __LOGICAL
      LOGICAL                 , DIMENSION(:), POINTER :: work
#endif
      INTEGER               :: i,k,lda,ldb,ldc,ldd,ldl_new,ldu_new
      LOGICAL               :: lcopy,lalloc,lrealloc
      REAL(ppm_kind_double) :: t0
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_alloc_1dl',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         CALL ppm_alloc_argcheck('ppm_alloc_1dl',iopt,ldl,1,info,ldu)
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
#elif __KIND == __LOGICAL
      work => work_1dl
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
            lda = LBOUND(adata,1)
            ldb = UBOUND(adata,1)
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
      ELSEIF (iopt.EQ.ppm_param_alloc_fit) THEN
         !----------------------------------------------------------------------
         !  fit memory but skip the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
            lda = LBOUND(adata,1)
            ldb = UBOUND(adata,1)
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
      ELSEIF (iopt.EQ.ppm_param_alloc_grow_preserve) THEN
         !----------------------------------------------------------------------
         !  grow memory and preserve the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
            lda = LBOUND(adata,1)
            ldb = UBOUND(adata,1)
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
      ELSEIF (iopt.EQ.ppm_param_alloc_grow) THEN
         !----------------------------------------------------------------------
         !  grow memory but skip the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
            lda = LBOUND(adata,1)
            ldb = UBOUND(adata,1)
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
      ELSEIF (iopt.EQ.ppm_param_dealloc) THEN
         !----------------------------------------------------------------------
         !  deallocate
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
            DEALLOCATE(adata,STAT=info)
            NULLIFY(adata)
            IF (info .NE. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_dealloc,'ppm_alloc_1dl',   &
     &             'DATA',__LINE__,info)
            ENDIF
         ENDIF
      ELSE
         !----------------------------------------------------------------------
         !  Unknown iopt
         !----------------------------------------------------------------------
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_alloc_1dl',                       &
     &                  'unknown iopt',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate new memory
      !-------------------------------------------------------------------------
      IF (lalloc) THEN
         ALLOCATE(work(ldl_new:ldu_new),STAT=info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_alloc_1dl',   &
     &           'WORK',__LINE__,info)
             GOTO 9999
         ENDIF
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
         !NULLIFY(adata)
         IF (info .NE. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_dealloc,'ppm_alloc_1dl',   &
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
      CALL substop('ppm_alloc_1dl',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_alloc_1dls
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_alloc_1dld
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_alloc_1dlsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_alloc_1dldc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_alloc_1dli
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_alloc_1dll
#endif

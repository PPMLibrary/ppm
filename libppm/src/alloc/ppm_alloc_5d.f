      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_alloc_5d
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_alloc_5ds(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 5D real single arrays
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_alloc_5dd(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 5D real double arrays
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_alloc_5dsc(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 5D complex single arrays
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_alloc_5ddc(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 5D complex double arrays
#elif __KIND == __INTEGER
      SUBROUTINE ppm_alloc_5di(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 5D integer arrays
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_alloc_5dl(adata,lda,iopt,info)
      !!! (Re)allocates the memory of 5D logical arrays
#endif
      !!! (pointers) based on the  number of elements.
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
      REAL(ppm_kind_single)   , DIMENSION(:,:,:,:,:), POINTER :: adata
#elif __KIND == __DOUBLE_PRECISION
      REAL(ppm_kind_double)   , DIMENSION(:,:,:,:,:), POINTER :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_single), DIMENSION(:,:,:,:,:), POINTER :: adata
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_double), DIMENSION(:,:,:,:,:), POINTER :: adata
#elif __KIND == __INTEGER
      INTEGER                 , DIMENSION(:,:,:,:,:), POINTER :: adata
#elif __KIND == __LOGICAL
      LOGICAL                 , DIMENSION(:,:,:,:,:), POINTER :: adata
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
      REAL(ppm_kind_single)   , DIMENSION(:,:,:,:,:), POINTER :: work
#elif __KIND == __DOUBLE_PRECISION
      REAL(ppm_kind_double)   , DIMENSION(:,:,:,:,:), POINTER :: work
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_single), DIMENSION(:,:,:,:,:), POINTER :: work
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_double), DIMENSION(:,:,:,:,:), POINTER :: work
#elif __KIND == __INTEGER
      INTEGER                 , DIMENSION(:,:,:,:,:), POINTER :: work
#elif __KIND == __LOGICAL
      LOGICAL                 , DIMENSION(:,:,:,:,:), POINTER :: work
#endif
      INTEGER, DIMENSION(5) :: ldb,ldc,lda_new
      INTEGER               :: i,j,k,l,m
      LOGICAL               :: lcopy,lalloc,lrealloc
      REAL(ppm_kind_double) :: t0
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_alloc_5d',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         CALL ppm_alloc_argcheck('ppm_alloc_5d',iopt,lda,5,info)
         IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Point to proper work array
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      work => work_5ds
#elif __KIND == __DOUBLE_PRECISION
      work => work_5dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      work => work_5dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      work => work_5ddc
#elif __KIND == __INTEGER
      work => work_5di
#elif __KIND == __LOGICAL
      work => work_5dl
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

            ldb(4) = SIZE(adata,4)
            IF (ldb(4).NE.lda(4)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
               lda_new(4) = lda(4)
            ELSE
               lda_new(4) = ldb(4)
            ENDIF

            ldb(5) = SIZE(adata,5)
            IF (ldb(5).NE.lda(5)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
               lda_new(5) = lda(5)
            ELSE
               lda_new(5) = ldb(5)
            ENDIF
         ELSE
            lalloc = .TRUE.
            lda_new(1) = lda(1)
            lda_new(2) = lda(2)
            lda_new(3) = lda(3)
            lda_new(4) = lda(4)
            lda_new(5) = lda(5)
         ENDIF
      ELSEIF (iopt.EQ.ppm_param_alloc_fit) THEN
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

            ldb(4) = SIZE(adata,4)
            IF (ldb(4).NE.lda(4)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lda_new(4) = lda(4)
            ELSE
               lda_new(4) = ldb(4)
            ENDIF

            ldb(5) = SIZE(adata,5)
            IF (ldb(5).NE.lda(5)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lda_new(5) = lda(5)
            ELSE
               lda_new(5) = ldb(5)
            ENDIF
         ELSE
            lalloc = .TRUE.
            lda_new(1) = lda(1)
            lda_new(2) = lda(2)
            lda_new(3) = lda(3)
            lda_new(4) = lda(4)
            lda_new(5) = lda(5)
         ENDIF
      ELSEIF (iopt.EQ.ppm_param_alloc_grow_preserve) THEN
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

            ldb(4) = SIZE(adata,4)
            IF (ldb(4).LT.lda(4)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
            ENDIF

            ldb(5) = SIZE(adata,5)
            IF (ldb(5).LT.lda(5)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
            ENDIF
            IF (lrealloc) THEN
               ! no dimension must shrink
               lda_new(1) = MAX(ldb(1),lda(1))
               lda_new(2) = MAX(ldb(2),lda(2))
               lda_new(3) = MAX(ldb(3),lda(3))
               lda_new(4) = MAX(ldb(4),lda(4))
               lda_new(5) = MAX(ldb(5),lda(5))
            ELSE
               lda_new(1) = ldb(1)
               lda_new(2) = ldb(2)
               lda_new(3) = ldb(3)
               lda_new(4) = ldb(4)
               lda_new(5) = ldb(5)
            ENDIF
         ELSE
            lalloc = .TRUE.
            lda_new(1) = lda(1)
            lda_new(2) = lda(2)
            lda_new(3) = lda(3)
            lda_new(4) = lda(4)
            lda_new(5) = lda(5)
         ENDIF
      ELSEIF (iopt.EQ.ppm_param_alloc_grow) THEN
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

            ldb(4) = SIZE(adata,4)
            IF (ldb(4).LT.lda(4)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
            ENDIF

            ldb(5) = SIZE(adata,5)
            IF (ldb(5).LT.lda(5)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
            ENDIF

            IF (lrealloc) THEN
               ! no dimension must shrink
               lda_new(1) = MAX(ldb(1),lda(1))
               lda_new(2) = MAX(ldb(2),lda(2))
               lda_new(3) = MAX(ldb(3),lda(3))
               lda_new(4) = MAX(ldb(4),lda(4))
               lda_new(5) = MAX(ldb(5),lda(5))
            ELSE
               lda_new(1) = ldb(1)
               lda_new(2) = ldb(2)
               lda_new(3) = ldb(3)
               lda_new(4) = ldb(4)
               lda_new(5) = ldb(5)
            ENDIF
         ELSE
            lalloc = .TRUE.
            lda_new(1) = lda(1)
            lda_new(2) = lda(2)
            lda_new(3) = lda(3)
            lda_new(4) = lda(4)
            lda_new(5) = lda(5)
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
               CALL ppm_error(ppm_err_dealloc,'ppm_alloc_5d',   &
     &             'DATA',__LINE__,info)
            ENDIF
         ENDIF
      ELSE
         !----------------------------------------------------------------------
         !  Unknown iopt
         !----------------------------------------------------------------------
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_alloc_5d',                       &
     &                  'unknown iopt',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate new memory
      !-------------------------------------------------------------------------
      IF (lalloc) THEN
         ALLOCATE(work(lda_new(1),lda_new(2),lda_new(3),lda_new(4),     &
     &       lda_new(5)),STAT=info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_alloc_5d',   &
     &           'WORK',__LINE__,info)
             GOTO 9999
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Make a copy of the old contents
      !-------------------------------------------------------------------------
      IF (lcopy) THEN
         DO k=1,5
            ldc(k) = MIN(lda_new(k),ldb(k))
         ENDDO
         DO m=1,ldc(5)
            DO l=1,ldc(4)
               DO k=1,ldc(3)
                  DO j=1,ldc(2)
                     DO i=1,ldc(1)
                        work(i,j,k,l,m) = adata(i,j,k,l,m)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  If reallocating, deallocate the old data
      !-------------------------------------------------------------------------
      IF (lrealloc) THEN
         DEALLOCATE(adata,STAT=info)
         !NULLIFY(adata)
         IF (info .NE. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_dealloc,'ppm_alloc_5d',   &
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
      CALL substop('ppm_alloc_5d',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_alloc_5ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_alloc_5dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_alloc_5dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_alloc_5ddc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_alloc_5di
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_alloc_5dl
#endif


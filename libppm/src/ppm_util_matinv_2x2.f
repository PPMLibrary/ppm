      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_util_matinv_2x2
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine computes the inverse of a 2x2 matrix
      !                 using explicit formulas.
      !
      !  Input        : Am(2,2)      (F) The 2x2 matrix stored in row-major
      !                                  order (i.e. first index is row
      !                                  number, 2nd is column number).
      !
      !  Input/output :                                            
      !
      !  Output       : Aminv(2,2)   (F) The inverse of Am.
      !                 info         (I) return status
      !
      !  Remarks      : This routine uses the formula inv(A) = adj(a)/det(A).
      !                 The det(A) is computed by expansion of the first
      !                 row. adj(A) is computed from all 1x1
      !                 sub-determinants of A.
      !
      !                 For ppm_debug .GT. 1 the routine checks its own
      !                 result. Maybe - after some time - this can be
      !                 removed.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_matinv_2x2.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.3  2006/09/04 18:34:59  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.2  2004/09/17 14:09:04  ivos
      !  Added result checks at the end if ppm_debug.GT.1.
      !
      !  Revision 1.1  2004/09/15 09:01:41  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_matinv_2x2_s(Am,Aminv,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_matinv_2x2_d(Am,Aminv,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(2,2), INTENT(IN   ) :: Am
      REAL(MK), DIMENSION(2,2), INTENT(  OUT) :: Aminv
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                                :: t0,Adet,Adetinv,lmyeps
      LOGICAL                                 :: correct
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_matinv_2x2',t0,info)
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

      !-------------------------------------------------------------------------
      !  Compute determinant of Am
      !-------------------------------------------------------------------------
      Adet    = (Am(1,1)*Am(2,2)) - (Am(1,2)*Am(2,1))

      !-------------------------------------------------------------------------
      !  Check that Am is invertible
      !-------------------------------------------------------------------------
      IF (Adet .LT. lmyeps) THEN
          info = ppm_error_warning
          CALL ppm_error(ppm_err_mat_singul,'ppm_util_matinv_2x2',     &
     &        'Cannot invert Am!. Returning.',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Invert determinant of Am
      !-------------------------------------------------------------------------
      Adetinv = 1.0_MK/Adet

      !-------------------------------------------------------------------------
      !  Compute inverse of Am using the adjoint formula
      !-------------------------------------------------------------------------
      Aminv(1,1) =  Am(2,2)*Adetinv
      Aminv(1,2) = -Am(1,2)*Adetinv
      Aminv(2,1) = -Am(2,1)*Adetinv
      Aminv(2,2) =  Am(1,1)*Adetinv

      !-------------------------------------------------------------------------
      !  Check result if debug is enabled.
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
         correct = .TRUE.
         IF (ABS(Am(1,1)*Aminv(1,1)+Am(1,2)*Aminv(2,1)-1.0_MK).GT.lmyeps) THEN
            correct = .FALSE.
         ENDIF
         IF (ABS(Am(1,1)*Aminv(1,2)+Am(1,2)*Aminv(2,2)).GT.lmyeps) THEN
            correct = .FALSE.
         ENDIF
         IF (ABS(Am(2,1)*Aminv(1,1)+Am(2,2)*Aminv(2,1)).GT.lmyeps) THEN
            correct = .FALSE.
         ENDIF
         IF (ABS(Am(2,1)*Aminv(1,2)+Am(2,2)*Aminv(2,2)-1.0_MK).GT.lmyeps) THEN
            correct = .FALSE.
         ENDIF
         IF (.NOT. correct) THEN
            info = ppm_error_warning
            CALL ppm_error(ppm_err_test_fail,'ppm_util_matinv_2x2',     &
     &           'Result is probably wrong!',__LINE__,info)
         ELSE
            CALL ppm_write(ppm_rank,'ppm_util_matinv_2x2',     &
     &           'The result is correct up to ppm_eps tolerance',info)
         ENDIF
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_matinv_2x2',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_matinv_2x2_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_matinv_2x2_d
#endif

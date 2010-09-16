      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_util_matinv_3x3
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_matinv_3x3_s(Am,Aminv,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_matinv_3x3_d(Am,Aminv,info)
#endif
      !!! This routine computes the inverse of a 3x3 matrix  using explicit
      !!! formulas.
      !!!
      !!! [NOTE]
      !!! This routine uses the formula `inv(A) = adj(a)/det(A)`.
      !!! The `det(A)` is computed by expansion of the first
      !!! row. `adj(A)` is computed from all 2x2  sub-determinants of A.
      !!!
      !!! [TIP]
      !!! For `ppm_debug > 1` the routine checks its own
      !!! result. Maybe - after some time - this can be removed.
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
      REAL(MK), DIMENSION(3,3), INTENT(IN   ) :: Am
      !!! The 3x3 matrix stored in row-major order (i.e. first index is row
      !!! number, 2nd is column number).
      REAL(MK), DIMENSION(3,3), INTENT(  OUT) :: Aminv
      !!! The inverse of Am
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                                :: t0,Adet,Adetinv,lmyeps
      REAL(MK), DIMENSION(3,3)                :: Udet,res
      LOGICAL                                 :: correct
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_matinv_3x3',t0,info)
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

      !-------------------------------------------------------------------------
      !  Compute all 2x2 subdeterminants of Am. The indices give row and
      !  column number of the supressed element.
      !-------------------------------------------------------------------------
      Udet(1,1) = (Am(2,2)*Am(3,3)) - (Am(3,2)*Am(2,3))
      Udet(1,2) = (Am(2,1)*Am(3,3)) - (Am(3,1)*Am(2,3))
      Udet(1,3) = (Am(2,1)*Am(3,2)) - (Am(3,1)*Am(2,2))

      Udet(2,1) = (Am(1,2)*Am(3,3)) - (Am(3,2)*Am(1,3))
      Udet(2,2) = (Am(1,1)*Am(3,3)) - (Am(3,1)*Am(1,3))
      Udet(2,3) = (Am(1,1)*Am(3,2)) - (Am(3,1)*Am(1,2))

      Udet(3,1) = (Am(1,2)*Am(2,3)) - (Am(2,2)*Am(1,3))
      Udet(3,2) = (Am(1,1)*Am(2,3)) - (Am(2,1)*Am(1,3))
      Udet(3,3) = (Am(1,1)*Am(2,2)) - (Am(2,1)*Am(1,2))

      !-------------------------------------------------------------------------
      !  Compute determinant of Am
      !-------------------------------------------------------------------------
      Adet =         Am(1,1)*Udet(1,1)
      Adet = Adet - (Am(1,2)*Udet(1,2))
      Adet = Adet + (Am(1,3)*Udet(1,3))

      !-------------------------------------------------------------------------
      !  Check that Am is invertible
      !-------------------------------------------------------------------------
      IF (Adet .LT. lmyeps) THEN
          info = ppm_error_warning
          CALL ppm_error(ppm_err_mat_singul,'ppm_util_matinv_3x3',     &
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
      Aminv(1,1) =  Udet(1,1)*Adetinv
      Aminv(1,2) = -Udet(2,1)*Adetinv
      Aminv(1,3) =  Udet(3,1)*Adetinv

      Aminv(2,1) = -Udet(1,2)*Adetinv
      Aminv(2,2) =  Udet(2,2)*Adetinv
      Aminv(2,3) = -Udet(3,2)*Adetinv

      Aminv(3,1) =  Udet(1,3)*Adetinv
      Aminv(3,2) = -Udet(2,3)*Adetinv
      Aminv(3,3) =  Udet(3,3)*Adetinv

      !-------------------------------------------------------------------------
      !  Check result if debug is enabled.
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
         correct = .TRUE.
         res = MATMUL(Am,Aminv)
         IF (ABS(res(1,1)-1.0_MK) .GT. lmyeps) THEN
            correct = .FALSE.
         ENDIF
         IF (ABS(res(1,2)       ) .GT. lmyeps) THEN
            correct = .FALSE.
         ENDIF
         IF (ABS(res(1,3)       ) .GT. lmyeps) THEN
            correct = .FALSE.
         ENDIF
         IF (ABS(res(2,1)       ) .GT. lmyeps) THEN
            correct = .FALSE.
         ENDIF
         IF (ABS(res(2,2)-1.0_MK) .GT. lmyeps) THEN
            correct = .FALSE.
         ENDIF
         IF (ABS(res(2,3)       ) .GT. lmyeps) THEN
            correct = .FALSE.
         ENDIF
         IF (ABS(res(3,1)       ) .GT. lmyeps) THEN
            correct = .FALSE.
         ENDIF
         IF (ABS(res(3,2)       ) .GT. lmyeps) THEN
            correct = .FALSE.
         ENDIF
         IF (ABS(res(3,3)-1.0_MK) .GT. lmyeps) THEN
            correct = .FALSE.
         ENDIF
         IF (.NOT. correct) THEN
            info = ppm_error_warning
            CALL ppm_error(ppm_err_test_fail,'ppm_util_matinv_3x3',     &
     &           'Result is probably wrong!',__LINE__,info)
         ELSE
            CALL ppm_write(ppm_rank,'ppm_util_matinv_3x3',     &
     &           'The result is correct up to ppm_eps tolerance',info)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_matinv_3x3',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_matinv_3x3_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_matinv_3x3_d
#endif

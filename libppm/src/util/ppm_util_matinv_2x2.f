      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_util_matinv_2x2
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
      SUBROUTINE ppm_util_matinv_2x2_s(Am,Aminv,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_matinv_2x2_d(Am,Aminv,info)
#endif
      !!! This routine computes the inverse of a 2x2 matrix using explicit
      !!! formulas.
      !!!
      !!! [NOTE]
      !!! This routine uses the formula `inv(A) = adj(a)/det(A)`.
      !!! The `det(A)` is computed by expansion of the first
      !!! row. `adj(A)` is computed from all 1x1 sub-determinants of A.
      !!!
      !!! [TIP]
      !!! For `ppm_debug > 1` the routine checks its own
      !!! result. Maybe - after some time - this can be
      !!! removed.
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
      !!! The 2x2 matrix stored in row-major order (i.e. first index is row
      !!! number, 2nd is column number).
      REAL(MK), DIMENSION(2,2), INTENT(  OUT) :: Aminv
      !!! The inverse of Am
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
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

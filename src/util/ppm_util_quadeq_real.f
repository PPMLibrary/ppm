      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_util_quadeq_real.f
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

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_quadeq_real_s(coef,roots,qtol,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_quadeq_real_d(coef,roots,qtol,info)
#endif
      !!! This routine solves a quadratic equation using the
      !!! well-known explicit solution formula. The equation
      !!! must only have REAL roots. If complex roots are
      !!! found, the routine returns `info = ppm_error_error` and exits.
      !!!
      !!! [NOTE]
      !!! If the discriminant is `< ppm_eps`, the roots are
      !!! considered equal and returned as exact bit-copy in
      !!! order to allow comparisons outside this routine.
      !!!
      !!! [TIP]
      !!! For `ppm_debug > 1` the routine checks its own  result.
      !!! Maybe - after some time - this can be removed.
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
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(3  ), INTENT(IN   ) :: coef
      !!! Coefficients. The equation being solved is:
      !!! --------------------------
      !!! coef(3)*x^2+coef(2)*x+coef(1) = 0
      !!! --------------------------
      REAL(MK), DIMENSION(2  ), INTENT(  OUT) :: roots
      !!! Contains the two real roots
      REAL(MK)                , INTENT(  OUT) :: qtol
      !!! Numerical tolerance (estimated) of result.
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, ppm_error_error if roots are complex
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                                :: t0,dis,lmyeps,Lf,stol
      REAL(MK), DIMENSION(3  )                :: res
      INTEGER                                 :: i
      LOGICAL                                 :: correct
      CHARACTER(LEN=ppm_char)                 :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_quadeq_real',t0,info)
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif
      roots = 0.0_MK

      !-------------------------------------------------------------------------
      !  Check that the equation is quadratic
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Discriminant of the quadratic equation 
      !-------------------------------------------------------------------------
      res(1) = coef(2)*coef(2)
      res(2) = -4.0_MK*coef(1)*coef(3)
      dis    = res(1) + res(2)

      !-------------------------------------------------------------------------
      !  Estimate numerical tolerance
      !-------------------------------------------------------------------------
      qtol = MAX(lmyeps*MAXVAL(res(1:2)),lmyeps)
      IF (ppm_debug .GT. 0) THEN
         WRITE(mesg,'(A,E12.4)') 'Numerical tolerance: ',qtol
         CALL ppm_write(ppm_rank,'ppm_util_quadeq_real',mesg,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that the equation has only real roots
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (dis .LT. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_quadeq_real',     &
     &            'Equation has complex roots. Exiting.',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  One double root or two different ones?
      !-------------------------------------------------------------------------
      IF (dis .LT. qtol) THEN
         !---------------------------------------------------------------------
         !  Solve for the double root
         !---------------------------------------------------------------------
         IF (ABS(coef(3)) .GT. lmyeps) THEN
            roots(1) = -0.5_MK*coef(2)
            roots(1) = roots(1)/coef(3)
            roots(2) = roots(1)
         ELSEIF (ABS(coef(2)) .GT. lmyeps) THEN
            roots(1) = -coef(1)/coef(2)
            roots(2) = roots(1)
         ENDIF
      ELSE
         !---------------------------------------------------------------------
         !  Solve for the 2 real roots
         !---------------------------------------------------------------------
         dis      = SQRT(dis)
         IF (ABS(coef(3)) .GT. lmyeps) THEN
            IF (coef(2) .LT. 0.0_MK) THEN
               roots(1) = -0.5_MK*(coef(2) - dis)
            ELSEIF (coef(2) .GT. 0.0_MK) THEN
               roots(1) = -0.5_MK*(coef(2) + dis)
            ELSE
               roots(1) = 0.0_MK
            ENDIF
            roots(2) = coef(1)/roots(1)
            roots(1) = roots(1)/coef(3)
         ELSEIF (ABS(coef(2)) .GT. lmyeps) THEN
            roots(1) = 0.0_MK
            roots(2) = -coef(1)/coef(2)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check result if debug is enabled.
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
         correct = .TRUE.
         DO i=1,2
            res(1) = coef(3)*roots(i)*roots(i)
            res(2) = coef(2)*roots(i)
            res(3) = coef(1)
            Lf = ABS(res(1)+res(2)+res(3))
            stol = MAX(MAXVAL(res)*10.0_MK*lmyeps,lmyeps)
            IF (Lf .GT. stol) THEN
               correct = .FALSE.
               info = ppm_error_warning
               WRITE(mesg,'(A,I1,2(A,E12.4))') 'Root ',i,              &
     &              ' is not correct. Error: ',Lf,' Tolerance: ',stol
               CALL ppm_error(ppm_err_test_fail,'ppm_util_quadeq_real',&
     &              mesg,__LINE__,info)
            ENDIF
         ENDDO
         IF (correct .AND. ppm_debug .GT. 0) THEN
            CALL ppm_write(ppm_rank,'ppm_util_quadeq_real',     &
     &           'Roots are correct to tolerance.',info)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_quadeq_real',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (ABS(coef(3)) .LT. lmyeps) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_quadeq_real',     &
     &            'Equation is not quadratic. Exiting.',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_quadeq_real_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_quadeq_real_d
#endif

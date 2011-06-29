      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_util_cubeq_real.f
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
      SUBROUTINE ppm_util_cubeq_real_s(coef,roots,qtol,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_cubeq_real_d(coef,roots,qtol,info)
#endif
      !!! This routine solves a cubic equation using the
      !!! explicit solution formula of Cardano. The equation
      !!! must only have REAL roots. If complex roots are
      !!! found, the routine returns info = ppm_error_error and exits.
      !!!
      !!! [NOTE]
      !!! Two roots that are closer together than ppm_eps are
      !!! considered equal and returned as exact bit-copy in
      !!! order to allow comparisons outside this routine.
      !!!
      !!! [NOTE]
      !!! For ppm_debug .GT. 1 the routine checks its own
      !!! result. Maybe - after some time - this can be removed.
      !!!
      !!! .References
      !!! *********************************************************************
      !!! - http://mathworld.wolfram.com/CubicEquation.html
      !!! *********************************************************************
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
      REAL(MK), DIMENSION(4  ), INTENT(IN   ) :: coef
      !!! Coefficients. The equation being solved is:
      !!! --------------------------
      !!! coef(4)*x^3+coef(3)*x^2+coef(2)*x+coef(1) = 0
      !!! --------------------------
      REAL(MK), DIMENSION(3  ), INTENT(  OUT) :: roots
      !!! Contains the three real roots
      REAL(MK)                , INTENT(  OUT) :: qtol
      !!! Numerical tolerance (estimated)  of result.
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, ppm_error_error if it has complex roots
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                                :: t0,dis,a0,a1,a2,lmyeps,stol
      REAL(MK)                                :: SpT,SpTi2,sqrt3,Q,R
      REAL(MK)                                :: Lf,ota2,rootstmp
      REAL(MK), DIMENSION(4  )                :: res
      COMPLEX(MK)                             :: sqrtD,S,T,SmT
      INTEGER                                 :: i
      LOGICAL                                 :: correct
      REAL(MK), PARAMETER                     :: one_third = 1.0_MK/3.0_MK
      CHARACTER(LEN=ppm_char)                 :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_cubeq_real',t0,info)
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif
      roots = 0.0_MK

      !-------------------------------------------------------------------------
      !  Check that the equation is cubic
      !-------------------------------------------------------------------------
      IF (ABS(coef(4)) .LT. lmyeps) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_util_cubeq_real',     &
     &        'Equation is not cubic. Exiting.',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Normalize the equation
      !-------------------------------------------------------------------------
      IF (coef(4) .NE. 1.0_MK) THEN
         dis = 1.0_MK/coef(4)
         a0  = coef(1)*dis
         a1  = coef(2)*dis
         a2  = coef(3)*dis
      ELSE
         a0  = coef(1)
         a1  = coef(2)
         a2  = coef(3)
      ENDIF

      !-------------------------------------------------------------------------
      !  Solve the cubic equation x^3 + a2*x^2 + a1*x + a0 = 0
      !-------------------------------------------------------------------------
      Q   = 3.0_MK*a1 - (a2*a2)
      Q   = Q/9.0_MK
      R   = 9.0_MK*a2*a1 - 27.0_MK*a0 - 2.0_MK*a2*a2*a2
      R   = R/54.0_MK
      res(1) = Q*Q*Q
      res(2) = R*R

      !-------------------------------------------------------------------------
      !  Estimate numerical tolerance
      !-------------------------------------------------------------------------
      dis = res(1) + res(2)
      qtol = MAX(MAXVAL(res(1:2))*lmyeps,lmyeps)
      IF (ppm_debug .GT. 0) THEN
         WRITE(mesg,'(A,E12.4)') 'Numerical tolerance: ',qtol
         CALL ppm_write(ppm_rank,'ppm_util_cubeq_real',mesg,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that the equation has only real roots
      !-------------------------------------------------------------------------
      IF (dis .GT. qtol) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_util_cubeq_real',     &
     &        'Equation has complex roots. Exiting.',__LINE__,info)
         GOTO 9999
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Solve the cubic equation x^3 + a2*x^2 + a1*x + a0 = 0
      !-------------------------------------------------------------------------
      sqrtD  = SQRT(CMPLX(dis,0.0_MK,MK))
      S      = (R + sqrtD)**one_third
      T      = (R - sqrtD)**one_third

      sqrt3  = SQRT(3.0_MK)
      SpT    = REAL(S+T,MK)
      SpTi2  = 0.5_MK*SpT
      SmT    = S - T
      Lf     = REAL((0.0_MK,0.5_MK)*sqrt3*SmT,MK)
      ota2   = -one_third*a2
      roots(1) = ota2 + SpT
      roots(2) = ota2 - SpTi2 + Lf
      roots(3) = ota2 - SpTi2 - Lf

      !-------------------------------------------------------------------------
      !  Check if we have multiple roots and set them to exact bit-copies
      !  if so (to allow comparisons outtside this routine). This is the
      !  case of dis=0.0. Use the average of the roots to take advantage
      !  of some error cancelling.
      !-------------------------------------------------------------------------
      IF (dis .GT. -qtol) THEN
         res(1) = ABS(roots(1)-roots(2))
         res(2) = ABS(roots(1)-roots(3))
         res(3) = ABS(roots(2)-roots(3))
         IF ((res(1).LT.qtol).AND.(res(2).LT.qtol).AND.(res(3).LT.qtol)) THEN
            ! all 3 are equal
            rootstmp = one_third*(roots(1)+roots(2)+roots(3))
            roots(1) = rootstmp
            roots(2) = roots(1)
            roots(3) = roots(1)
         ELSE
            IF (res(1).LT.qtol) THEN
               ! 1 and 2 are equal
               rootstmp = 0.5_MK*(roots(1)+roots(2))
               roots(1) = rootstmp
               roots(2) = roots(1)
            ELSEIF (res(2) .LT. qtol) THEN
               ! 1 and 3 are equal
               rootstmp = 0.5_MK*(roots(1)+roots(3))
               roots(1) = rootstmp
               roots(3) = roots(1)
            ELSEIF (res(3) .LT. qtol) THEN
               ! 2 and 3 are equal
               rootstmp = 0.5_MK*(roots(2)+roots(3))
               roots(2) = rootstmp
               roots(3) = roots(2)
            ENDIF
         ENDIF
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Check result if debug is enabled.
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
         correct = .TRUE.
         DO i=1,3
            res(3) = roots(i)*roots(i)
            res(1) = coef(4)*res(3)*roots(i)
            res(2) = coef(3)*res(3)
            res(3) = coef(2)*roots(i)
            res(4) = coef(1)
            Lf = ABS(res(1)+res(2)+res(3)+res(4))
            stol = MAX(MAXVAL(res)*10.0_MK*lmyeps,lmyeps)
            IF (Lf .GT. stol) THEN
               correct = .FALSE.
               info = ppm_error_warning
               WRITE(mesg,'(A,I1,2(A,E12.4))') 'Root ',i,              &
     &              ' is not correct. Error: ',Lf,' Tolerance: ',stol
               CALL ppm_error(ppm_err_test_fail,'ppm_util_cubeq_real', &
     &              mesg,__LINE__,info)
            ENDIF
         ENDDO
         IF (correct .AND. ppm_debug .GT. 0) THEN
            CALL ppm_write(ppm_rank,'ppm_util_cubeq_real',     &
     &           'Roots are correct to tolerance.',info)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_cubeq_real',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_cubeq_real_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_cubeq_real_d
#endif

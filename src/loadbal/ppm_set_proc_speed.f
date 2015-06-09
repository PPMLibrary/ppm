      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_set_proc_speed
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
      SUBROUTINE ppm_set_proc_speed_s(proc_speed,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_set_proc_speed_d(proc_speed,info)
#endif
      !!! This routine can be used by the user to set the relative speeds of
      !!! the processors (is used in subs2proc for load balancing).
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(0:), INTENT(IN   ) :: proc_speed
      !!! Relative speeds of all processors from 0 to ppm_nproc-1. The numbers
      !!! must sum up to 1.
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)               :: t0,lmyeps
      REAL(ppm_kind_double)  :: rsum
      INTEGER                :: i
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart('ppm_set_proc_speed',t0,info)

#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that numbers add up to 1
      !-------------------------------------------------------------------------
      rsum = 0.0_ppm_kind_double
      DO i=0,ppm_nproc-1
        rsum = rsum + REAL(proc_speed(i),ppm_kind_double)
      ENDDO
      IF (ABS(rsum - 1.0_MK) .GT. lmyeps) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,'ppm_set_proc_speed', &
     &        'proc_speed must sum up to 1.0',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Set the internal ppm_proc_speed values
      !-------------------------------------------------------------------------
      rsum = 0.0_ppm_kind_double
      DO i=0,ppm_nproc-2
#if   __KIND == __SINGLE_PRECISION
          ppm_proc_speed(i) = REAL(proc_speed(i),ppm_kind_double)
#elif __KIND == __DOUBLE_PRECISION
          ppm_proc_speed(i) = proc_speed(i)
#endif
          rsum = rsum + ppm_proc_speed(i)
      ENDDO
      ppm_proc_speed(ppm_nproc-1) = 1.0_ppm_kind_double - rsum


      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_set_proc_speed',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_set_proc_speed',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (SIZE(proc_speed,1) .LT. ppm_nproc) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_set_proc_speed', &
     &            'proc_speed must be at least of length nproc',__LINE__,info)
             GOTO 8888
          ENDIF
          IF (proc_speed(i) .LT. 0.0_MK) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_set_proc_speed', &
     &           'proc_speed must be >= 0 for all processors',__LINE__,info)
             GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_set_proc_speed_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_set_proc_speed_d
#endif

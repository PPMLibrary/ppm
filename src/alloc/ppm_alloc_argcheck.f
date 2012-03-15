      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_alloc_argcheck
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

      SUBROUTINE ppm_alloc_argcheck(caller,iopt,ldl,dimension,info,ldu)
      !!! Checks the arguments of the ppm_alloc routines

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop

      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CHARACTER(LEN=*)               , INTENT(IN)   :: caller
      !!! name of calling subroutine
      INTEGER, DIMENSION(:)          , INTENT(IN)   :: ldl
      !!! Lower index limit in leading dim.
      INTEGER, OPTIONAL, DIMENSION(:), INTENT(IN)   :: ldu
      !!! Upper index limit in leading dim. (>ldl(1)).
      !!! OPTIONAL in ppm_alloc_*dl.f
      INTEGER                        , INTENT(IN)   :: iopt
      !!! Allocation mode. One of:
      !!!
      !!! * ppm_param_alloc_fit
      !!! * ppm_param_alloc_fit_preserve
      !!! * ppm_param_alloc_grow
      !!! * ppm_param_alloc_grow_preserve
      !!! * ppm_param_dealloc
      INTEGER                        , INTENT(IN)   :: dimension
      !!! Helps to determine the dimension of the caller subroutine (1-5)
      INTEGER                        , INTENT(OUT)  :: info
      !!! Returns status, 0 upon success.
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER               :: i

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      info = 0
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (iopt .NE. ppm_param_alloc_fit           .AND.                       &
     &    iopt .NE. ppm_param_alloc_fit_preserve  .AND.                       &
     &    iopt .NE. ppm_param_alloc_grow          .AND.                       &
     &    iopt .NE. ppm_param_alloc_grow_preserve .AND.                       &
     &    iopt .NE. ppm_param_dealloc) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,'unknown iopt',__LINE__,info)
        GOTO 9999
      ENDIF
      IF (iopt .NE. ppm_param_dealloc) THEN
        IF (PRESENT(ldu)) THEN
            DO i=1,dimension
                IF (ldl(i) .GT. ldu(i)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,caller,  &
     &              'ldu() must be >= ldl()',__LINE__,info)
                GOTO 9999
                ENDIF
            ENDDO
        ELSE
            DO i=1,dimension
                IF (ldl(i) .LT. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,caller,  &
     &              'ldl() must be >= 0',__LINE__,info)
                GOTO 9999
                ENDIF
            ENDDO
        ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_alloc_argcheck

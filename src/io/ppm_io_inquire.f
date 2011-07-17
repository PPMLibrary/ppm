      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_io_inquire
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

      SUBROUTINE ppm_io_inquire(iUnit,isopen,mode,pfmt,stat)
      !!! This routine inquires the state of a given ppm I/O unit.
      !!!
      !!! [NOTE]
      !!! All output arguments are *optional*. It is adivised
      !!! to use explicit argument naming when using this routine. Example:
      !!!
      !!! --------------------------
      !!! CALL ppm_inquire(20,MODE=imode,STAT=info)
      !!! --------------------------

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_io
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER             , INTENT(IN   ) :: iUnit
      !!! IO unit to be inquired
      INTEGER  , OPTIONAL , INTENT(  OUT) :: stat
      !!! Return code. 0 on success.
      INTEGER  , OPTIONAL , INTENT(  OUT) :: pfmt
      !!! I/O format. One of:
      !!!
      !!! * ppm_param_io_ascii
      !!! * ppm_param_io_binary
      !!! * ppm_param_undefined
      INTEGER  , OPTIONAL , INTENT(  OUT) :: mode
      !!! I/O mode. One of:
      !!!
      !!! * ppm_param_io_centralized
      !!! * ppm_param_io_distributed
      !!! * ppm_param_undefined
      LOGICAL  , OPTIONAL , INTENT(  OUT) :: isopen
      !!! Either:
      !!!
      !!! * `.TRUE.` if the unit is open on the local processor.
      !!! * `.FALSE.` otherwise
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)            :: t0
      LOGICAL                          :: lopen
      INTEGER                          :: info
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_io_inquire',t0,info)

      !-------------------------------------------------------------------------
      !  Check argument
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF


      !-------------------------------------------------------------------------
      !  Check if the unit number is open
      !-------------------------------------------------------------------------
      lopen = .FALSE.
      IF(ASSOCIATED(ppm_io_unit)) THEN
          IF (SIZE(ppm_io_unit,1) .GE. iUnit) THEN
              IF (ppm_io_unit(iUnit) .GT. 0) lopen = .TRUE.
          ENDIF
      ENDIF
      IF (PRESENT(isopen)) isopen = lopen

      !-------------------------------------------------------------------------
      !  Get the other information if present
      !-------------------------------------------------------------------------
      IF (lopen) THEN
          IF (PRESENT(pfmt)) pfmt = ppm_io_format(ppm_io_unit(iUnit))
          IF (PRESENT(mode)) mode = ppm_io_mode(ppm_io_unit(iUnit))
      ELSE
          IF (PRESENT(pfmt)) pfmt = ppm_param_undefined
          IF (PRESENT(mode)) mode = ppm_param_undefined
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      IF (PRESENT(stat)) stat = info
      CALL substop('ppm_io_inquire',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (.NOT. ppm_initialized) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_ppm_noinit,'ppm_io_inquire',  &
     &            'Please call ppm_init first!',__LINE__,info)
           GOTO 8888
        ENDIF
        IF (iUnit .LE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_io_inquire',    &
            &        'Unit number needs to be > 0',__LINE__,info)
            IF (PRESENT(isopen)) isopen = .FALSE.
            IF (PRESENT(pfmt)) pfmt = ppm_param_undefined
            IF (PRESENT(mode)) mode = ppm_param_undefined
            GOTO 8888
        ENDIF
 8888   CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_io_inquire

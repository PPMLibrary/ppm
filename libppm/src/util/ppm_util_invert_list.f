      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_util_invert_list
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

      SUBROUTINE ppm_util_invert_list(inlist,outlist,info)
      !!! Creates the inverse list (indiced become entries
      !!! and vice versa) for a 1D integer list.
      !!!
      !!! [NOTE]
      !!! Entries in outlist for which no corresponding entry
      !!! in inlist exists are initialized to
      !!! `-HUGE(outlist(.))` and can be recognized like this.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------

      INTEGER, DIMENSION(:) , POINTER       :: inlist
      !!! Non-negative integer list to be inverted
      INTEGER, DIMENSION(:) , POINTER       :: outlist
      !!! Inverse list
      INTEGER               , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)                 :: t0
      INTEGER, DIMENSION(1)                 :: ldl,ldu
      INTEGER                               :: iopt
      INTEGER                               :: i,j,inmin,inmax,outmin,outmax
      CHARACTER(LEN=ppm_char)               :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_invert_list',t0,info)
      inmin = LBOUND(inlist,1)
      inmax = UBOUND(inlist,1)
      IF (ppm_debug .GE. 2) THEN
          WRITE(mesg,'(A,2I10)') 'Input list bounds ', inmin, inmax
          CALL ppm_write(ppm_rank,'ppm_util_invert_list',  &
    &                    mesg,info)
      ENDIF
      !-------------------------------------------------------------------------
      !  Determine min and max value in inlist
      !-------------------------------------------------------------------------
      outmax = -HUGE(outmax)
      outmin = HUGE(outmin)
      DO i=inmin,inmax
          IF (inlist(i) .GT. outmax) outmax = inlist(i)
          IF (inlist(i) .LT. outmin) outmin = inlist(i)
      ENDDO
      
      !-------------------------------------------------------------------------
      !  Case of an empty input list
      !-------------------------------------------------------------------------
      IF (outmax .LT. outmin) THEN
          ! this is equivalent to testing (inmax .LT. inmin)
          ! We allocate an empty inverse list
          iopt = ppm_param_alloc_fit
          ldu(1) = 0
          CALL ppm_alloc(outlist,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_invert_list',        &
    &              'inverted list OUTLIST',__LINE__,info)
              GOTO 9999
          ENDIF
      !-------------------------------------------------------------------------
      !  General case: non-empty list
      !-------------------------------------------------------------------------
      ELSE
      !-------------------------------------------------------------------------
      !  Allocate outlist. No need to preserve its contents
      !-------------------------------------------------------------------------
          iopt = ppm_param_alloc_fit
          ldl(1) = outmin
          ldu(1) = outmax
          IF (ppm_debug .GE. 2) THEN
              WRITE(mesg,'(A,2I10)') 'Output list bounds ', outmin, outmax
              CALL ppm_write(ppm_rank,'ppm_util_invert_list',  &
    &              mesg,info)
          ENDIF
          CALL ppm_alloc(outlist,ldl,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_invert_list',  &
    &                        'inverted list OUTLIST',__LINE__,info)
              GOTO 9999
          ENDIF

      !-------------------------------------------------------------------------
      !  Initialize outlist
      !-------------------------------------------------------------------------
          DO i=outmin,outmax
              outlist(i) = -HUGE(outlist(i))
          ENDDO

      !-------------------------------------------------------------------------
      !  Build inverse list
      !-------------------------------------------------------------------------
          DO i=inmin,inmax
              outlist(inlist(i)) = i
          ENDDO
      END IF
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_invert_list',t0,info)
      RETURN
      END SUBROUTINE ppm_util_invert_list

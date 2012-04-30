      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_io_delete
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

      SUBROUTINE ppm_io_delete(filename,mode,info)
      !!! This routine deletes a file.

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
      CHARACTER(LEN=*)    , INTENT(IN   ) :: filename
      !!! Name of the file to be deleted
      INTEGER             , INTENT(IN   ) :: mode
      !!! I/O mode. One of:
      !!!
      !!! * ppm_param_io_distributed
      !!! * ppm_param_io_centralized
      !!!
      !!! In the distributed mode, every processor deletes the file
      !!! locally. In the centralized mode, only processor 0 deletes the file.
      INTEGER             , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)            :: t0
      LOGICAL                          :: lexists,lopen
      INTEGER                          :: iUnit
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_io_delete',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_io_delete',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      IF ((mode .NE. ppm_param_io_distributed) .AND.        &
     &    (mode .NE. ppm_param_io_centralized)) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,'ppm_io_delete',     &
     &        'Invalid mode specified.',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that file exists
      !-------------------------------------------------------------------------
      IF (mode .EQ. ppm_param_io_distributed) THEN
          INQUIRE(FILE=filename,EXIST=lexists)
          IF (.NOT. lexists) THEN
             info = ppm_error_notice
             CALL ppm_error(ppm_err_file,'ppm_io_delete',    &
     &           TRIM(filename),__LINE__,info)
             GOTO 9999
          ENDIF
      ELSE
          IF (ppm_rank .EQ. 0) THEN
              INQUIRE(FILE=filename,EXIST=lexists)
              IF (.NOT. lexists) THEN
                 info = ppm_error_notice
                 CALL ppm_error(ppm_err_file,'ppm_io_delete',    &
     &               TRIM(filename),__LINE__,info)
                 GOTO 9999
              ENDIF
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that file is not connected to a unit
      !-------------------------------------------------------------------------
      IF (mode .EQ. ppm_param_io_distributed) THEN
          INQUIRE(FILE=filename,OPENED=lopen)
          IF (lopen) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_delete,'ppm_io_delete',    &
     &           'Cannot delete an open file',__LINE__,info)
             GOTO 9999
          ENDIF
      ELSE
          IF (ppm_rank .EQ. 0) THEN
              INQUIRE(FILE=filename,OPENED=lopen)
              IF (lopen) THEN
                 info = ppm_error_warning
                 CALL ppm_error(ppm_err_delete,'ppm_io_delete',    &
     &               'Cannot delete an open file',__LINE__,info)
                 GOTO 9999
              ENDIF
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Attempt delete
      !-------------------------------------------------------------------------
      IF (mode .EQ. ppm_param_io_distributed) THEN
          !---------------------------------------------------------------------
          !  Find unused unit number
          !---------------------------------------------------------------------
          iUnit = 9
          CALL ppm_io_unused_unit(iUnit,info)
          IF (iUnit .LT. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_open,'ppm_io_delete',    &
     &              'No I/O unit available.',__LINE__,info)
             GOTO 9999
          ENDIF
          !---------------------------------------------------------------------
          !  Delete
          !---------------------------------------------------------------------
          OPEN(iUnit,FILE=filename,STATUS='OLD',ACTION='WRITE',IOSTAT=info)
          CLOSE(iUnit,STATUS='DELETE',IOSTAT=info)
          IF (info .NE. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_delete,'ppm_io_delete',    &
     &           TRIM(filename),__LINE__,info)
             GOTO 9999
          ENDIF
      ELSE
          IF (ppm_rank .EQ. 0) THEN
              !-----------------------------------------------------------------
              !  Find unused unit number
              !-----------------------------------------------------------------
              iUnit = 9
              CALL ppm_io_unused_unit(iUnit,info)
              IF (iUnit .LT. 0) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_open,'ppm_io_delete',    &
     &              'No I/O unit available.',__LINE__,info)
                 GOTO 9999
              ENDIF
              !-----------------------------------------------------------------
              !  Delete
              !-----------------------------------------------------------------
              OPEN(iUnit,FILE=filename,STATUS='OLD',ACTION='WRITE',IOSTAT=info)
              CLOSE(iUnit,STATUS='DELETE',IOSTAT=info)
              IF (info .NE. 0) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_delete,'ppm_io_delete',    &
     &               TRIM(filename),__LINE__,info)
                 GOTO 9999
              ENDIF
          ENDIF
      ENDIF
       
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_io_delete',t0,info)
      RETURN
      END SUBROUTINE ppm_io_delete

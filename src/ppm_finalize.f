      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_finalize
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

      SUBROUTINE ppm_finalize(info)
      !!! This routine terminates the ppm library including
      !!! deallocating memory and closing files.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data, ONLY : ppm_char,ppm_debug,ppm_initialized,          &
      &   ppm_stdout,ppm_stderr,ppm_logfile,ppm_error_error,ppm_param_dealloc, &
      &   ppm_proc_speed,ppm_rank
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_log
      USE ppm_module_mapping_typedef
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(  OUT) :: info
      !!! Returns 0 upon success

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(1) :: lda
      INTEGER               :: iopt
      INTEGER               :: istat

      REAL(ppm_kind_double) :: t0

      LOGICAL :: isopen

      CHARACTER(LEN=ppm_char) :: caller='ppm_finalize'
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      lda(1)=0

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         IF (.NOT.ppm_initialized) THEN
            fail('Please call ppm_init first!',ppm_err_ppm_noinit)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate global arrays (from ppm_module)
      !-------------------------------------------------------------------------
      istat = 0
      iopt = ppm_param_dealloc
      CALL ppm_alloc(ppm_sendbufferd,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_sendbuffers,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_recvbufferd,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_recvbuffers,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_psendbuffer,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_precvbuffer,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_buffer2part,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_isendlist,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_irecvlist,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_buffer_dim,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_buffer_type,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_proc_speed,lda,iopt,info)
      istat = istat + info

      IF (istat.NE.0) THEN
         WRITE(cbuf,'(A,I3,A)') 'for ',istat,' global arrays. Possible memory leak.'
         fail(cbuf,ppm_err_dealloc)
      ENDIF

      !-------------------------------------------------------------------------
      !  Set the global status
      !-------------------------------------------------------------------------
      ppm_initialized = .FALSE.

      IF (ppm_rank.EQ.0) THEN
         stdout_f('(A)',"***********************************************************")
         CALL ppm_log(caller,cbuf,info)
         stdout_f('(A)',"***         Parallel Particle Mesh Library (PPM)        ***")
         CALL ppm_log(caller,cbuf,info)
         stdout_f('(A)',"***                                                     ***")
         CALL ppm_log(caller,cbuf,info)
         stdout_f('(A)',"***  Thank you for using the PPM library in your work.  ***")
         CALL ppm_log(caller,cbuf,info)
         stdout_f('(A)',"***  Please, acknowledge our efforts by including the   ***")
         CALL ppm_log(caller,cbuf,info)
         stdout_f('(A)',"***  following reference when publishing data obtained  ***")
         CALL ppm_log(caller,cbuf,info)
         stdout_f('(A)',"***  using PPM library:                                 ***")
         CALL ppm_log(caller,cbuf,info)
         stdout_f('(A)',"***                                                     ***")
         CALL ppm_log(caller,cbuf,info)
         stdout_f('(A)',"***  I. F. Sbalzarini, J. Walther, M. Bergdorf,         ***")
         CALL ppm_log(caller,cbuf,info)
         stdout_f('(A)',"***  S. E. Hieber, E. M. Kotsalis, P. Koumoutsakos,     ***")
         CALL ppm_log(caller,cbuf,info)
         stdout_f('(A)',"***  `J. Comput. Phys.`, 215, pp. 566-588, (2006)       ***")
         CALL ppm_log(caller,cbuf,info)
         stdout_f('(A)',"***                                                     ***")
         CALL ppm_log(caller,cbuf,info)
         stdout_f('(A)',"***********************************************************")
         CALL ppm_log(caller,cbuf,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Close output units if needed
      !-------------------------------------------------------------------------
      IF (ppm_stdout.GE.0) THEN
         ! do not close it if it is a system standard unit
         IF (ppm_stdout .NE. 6) THEN
            INQUIRE(ppm_stdout,OPENED=isopen)
            IF (isopen) CLOSE(ppm_stdout)
         ENDIF
      ENDIF
      IF (ppm_stderr.GE.0) THEN
         ! do not close it if it is a system standard unit
         IF (ppm_stderr .NE. 0) THEN
            INQUIRE(ppm_stderr,OPENED=isopen)
            IF (isopen) CLOSE(ppm_stderr)
         ENDIF
      ENDIF
      IF (ppm_logfile.GE.0) THEN
         INQUIRE(ppm_logfile,OPENED=isopen)
         IF (isopen) CLOSE(ppm_logfile)
      ENDIF

      RETURN

      END SUBROUTINE ppm_finalize

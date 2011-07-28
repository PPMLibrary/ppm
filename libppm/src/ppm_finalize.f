      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_finalize
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_finalize(info)
      !!! This routine terminates the ppm library including
      !!! deallocating memory and closing files.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_mesh_finalize
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(OUT) :: info
      !!! Returns 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(1) :: lda
      INTEGER               :: iopt
      INTEGER               :: istat
      REAL(ppm_kind_double) :: t0
      LOGICAL               :: isopen
      CHARACTER(LEN=ppm_char) :: mesg
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_finalize',t0,info)
      lda(1) = 0

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_finalize',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
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

      IF (istat .NE. 0) THEN
          WRITE(mesg,'(A,I3,A)') 'for ',istat,' global arrays. Possible memory leak.'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_finalize',mesg,__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate mesh structures (from ppm_module_mesh)
      !-------------------------------------------------------------------------
      CALL ppm_mesh_finalize(info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_finalize',  &
     &        'Mesh deallocation failed',__LINE__,info)
          GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      !  Set the global status
      !-------------------------------------------------------------------------
      ppm_initialized = .FALSE.

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_finalize',t0,info)

      !-------------------------------------------------------------------------
      !  Close output units if needed
      !-------------------------------------------------------------------------
      IF (ppm_stdout .GE. 0) THEN
          ! do not close it if it is a system standard unit
          IF (ppm_stdout .NE. 6) THEN
              INQUIRE(ppm_stdout,opened=isopen)
              IF (isopen) CLOSE(ppm_stdout)
          ENDIF
      ENDIF
      IF (ppm_stderr .GE. 0) THEN
          ! do not close it if it is a system standard unit
          IF (ppm_stderr .NE. 0) THEN
              INQUIRE(ppm_stderr,opened=isopen)
              IF (isopen) CLOSE(ppm_stderr)
          ENDIF
      ENDIF
      IF (ppm_logfile .GE. 0) THEN
          INQUIRE(ppm_logfile,opened=isopen)
          IF (isopen) CLOSE(ppm_logfile)
      ENDIF

      RETURN

      END SUBROUTINE ppm_finalize

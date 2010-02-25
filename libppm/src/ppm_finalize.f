      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_finalize
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine terminates the ppm library including 
      !                 deallocating memory and closing files.
      !
      !  Input        :
      !
      !  Input/output : 
      !
      !  Output       : info    (I) 0 on success.
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_finalize.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.22  2006/09/04 18:34:45  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.21  2004/10/01 16:33:32  ivos
      !  cosmetics.
      !
      !  Revision 1.20  2004/10/01 16:08:59  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.19  2004/07/26 07:46:36  ivos
      !  Changed to use single-interface modules. Updated all USE statements.
      !
      !  Revision 1.18  2004/07/23 12:53:07  ivos
      !  bugfix: IO units are now only closed if they are not std units.
      !  The other confused OS X.
      !
      !  Revision 1.17  2004/07/16 14:46:24  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.16  2004/07/16 14:06:51  ivos
      !  Added CALL to ppm_mesh_finalize to deallocate the PRIVATE variables
      !  of module_mesh.
      !
      !  Revision 1.15  2004/07/16 11:05:12  ivos
      !  ppm_initialized is now set.
      !
      !  Revision 1.14  2004/04/02 15:23:59  ivos
      !  Added DEALLOCs for mesh ghost map lists.
      !
      !  Revision 1.13  2004/03/01 12:02:14  ivos
      !  bogfix: moved closing of output units AFTER call the substop, because
      !  the latter still prints something. Errors occurred on certain systems
      !  if the unit was closed before.
      !
      !  Revision 1.12  2004/02/25 13:55:09  ivos
      !  Added deallocation of ppm_proc_speed.
      !
      !  Revision 1.11  2004/02/24 14:14:11  ivos
      !  Added some deallocs for global arrays and also for all mesh arrays.
      !  added success checks for all deallocs.
      !
      !  Revision 1.10  2004/02/10 08:32:42  hiebers
      !  fixed compilation error when using ifc
      !
      !  Revision 1.9  2004/01/23 17:24:15  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.8  2004/01/23 11:31:21  ivos
      !  Cleanup: (1) updated headers, (2) inserted ppm_error and ppm_write,
      !  (3) added argument checking, (4) added checks after every alloc.
      !
      !  Revision 1.7  2004/01/13 12:31:44  ivos
      !  Added CLOSE for all output units.
      !
      !  Revision 1.6  2003/12/16 12:22:21  ivos
      !  Added DEALLOCs for ppm_psendbuffer, ppm_isendlist, ppm_irecvlist,
      !  ppm_buffer_dim, ppm_buffer_type, ppm_isublist, ppm_nsublist.
      !
      !  Revision 1.5  2003/12/12 17:44:40  ivos
      !  Removed 2 no longer needed DEALLOCs.
      !
      !  Revision 1.4  2003/12/12 16:15:56  ivos
      !  Added DEALLOCs for the global internal topology lists.
      !
      !  Revision 1.3  2003/12/09 09:00:46  ivos
      !  Added DEALLOCATEs for the global variables ppm_isoptimized, 
      !  ppm_icommseq and ppm_ncommseq (declared in ppm_module.f).
      !
      !  Revision 1.2  2003/12/05 14:41:53  ivos
      !  Deallocate the following globals: ppm_isoptimized, ppm_nneighlist,
      !  ppm_ineighlist.
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_finalize(info)
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
      CALL ppm_alloc(ppm_nneighlist,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_ineighlist,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_nsublist,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_isublist,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_nsubs,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_min_subd,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_max_subd,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_min_subs,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_max_subs,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_subs2proc,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_isoptimized,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_icommseq,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_ncommseq,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_user_topoid,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_internal_topoid,lda,iopt,info)
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

      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Initialisation of the ppm library.
      !
      !  Input        : dim     (I) dimension of the problem (2 or 3)
      !                 prec    (I) desired precision 
      !                                 ppm_kind_double or
      !                                 ppm_kind_single
      !                 tolexp  (I) log10 of the numerical tolerance for real
      !                             comparisons (> machine eps). E.g. for a
      !                             tolerance of 10^(-14), tolexp = -14.
      !                             Everything .LT. 10^tolexp will be
      !                             considered zero (0) and its inversion
      !                             is considered dangerous.
      !                 comm    (I) MPI communicator
      !                 debug   (I) The debug level (0,1 or 2)
      !                 logfile (I) OPTIONAL. Unit number where to print
      !                             log messages to.
      !                 stderr  (I) OPTIONAL. Unit number where to print 
      !                             the stderr to
      !                 stdout  (I) OPTIONAL. Unit number where to print 
      !                             the stdout to
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
      !  $Log: ppm_init.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.35  2006/09/04 18:34:47  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.34  2004/10/01 16:33:33  ivos
      !  cosmetics.
      !
      !  Revision 1.33  2004/10/01 16:09:00  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.32  2004/09/22 10:33:17  ivos
      !  Updated comment header for TOLEXP.
      !
      !  Revision 1.31  2004/07/26 11:46:53  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.30  2004/07/26 10:50:19  hiebers
      !  add ppm_pi
      !
      !  Revision 1.29  2004/07/26 07:47:09  ivos
      !  Updated to use single-interface modules. Changed all USE statements.
      !
      !  Revision 1.28  2004/07/16 14:47:41  ivos
      !  Added check for attempted double initialization.
      !
      !  Revision 1.27  2004/07/16 11:05:12  ivos
      !  ppm_initialized is now set.
      !
      !  Revision 1.26  2004/06/10 16:20:00  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.25  2004/06/01 09:26:55  ivos
      !  bigfix: ppm_stdout, ppm_stderr and ppm_logfile are no longer
      !  directly set, but ppm_io_set_unit is used.
      !
      !  Revision 1.24  2004/05/06 07:42:11  ivos
      !  renamed ppm_set_unit to ppm_io_set_unit.
      !
      !  Revision 1.23  2004/04/14 15:45:09  ivos
      !  bugfix: ppm_max_topoid must be initialized to 0 and not -1.
      !  Otherwise the first topology created would be 0 and not 1!
      !
      !  Revision 1.22  2004/04/14 15:02:26  ivos
      !  ppm_topoid, ppm_field_topoid and ppm_max_topoid are now initialized to
      !  -1 because 0 is the ring topology now.
      !
      !  Revision 1.21  2004/04/13 12:38:42  ivos
      !  Now gets and outputs (and logs) host names it is running on.
      !
      !  Revision 1.20  2004/04/05 11:57:52  ivos
      !  Now always setting ppm_myeps in both precisions, irrespective of
      !  ppm_kind.
      !
      !  Revision 1.19  2004/03/01 12:02:30  ivos
      !  Added CALL to ppm_print_defines.
      !
      !  Revision 1.18  2004/02/25 14:44:30  ivos
      !  Changed indices in ppm_proc_speed to 0..nproc-1.
      !
      !  Revision 1.17  2004/02/25 13:54:52  ivos
      !  Added allocation and initialization of ppm_proc_speed.
      !
      !  Revision 1.16  2004/02/20 15:44:38  ivos
      !  bugfix: ppm_topoid and ppm_field_topoid are now both initialized to 0.
      !
      !  Revision 1.15  2004/02/18 17:49:20  walther
      !  Bug fix: ppm_set_unit should be called with ppm_ variables.
      !
      !  Revision 1.14  2004/02/18 17:45:08  walther
      !  Added optional arguments (stdout, stderr, logfile) and ppm_init is now
      !  ppm_set_unit (and using the new ppm_module_io).
      !
      !  Revision 1.13  2004/02/12 17:47:14  ivos
      !  Added OPTIONAL argument to enable log file output of parameters.
      !
      !  Revision 1.12  2004/02/12 17:22:26  ivos
      !  Added debug and log output of parameters.
      !
      !  Revision 1.11  2004/02/12 14:46:01  ivos
      !  Added check and initialization of real tolerance.
      !
      !  Revision 1.10  2004/02/11 14:31:06  ivos
      !  Added initialization of ppm_field_topoid
      !
      !  Revision 1.9  2004/01/23 17:24:15  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.8  2004/01/23 11:31:21  ivos
      !  Cleanup: (1) updated headers, (2) inserted ppm_error and ppm_write,
      !  (3) added argument checking, (4) added checks after every alloc.
      !
      !  Revision 1.7  2004/01/19 11:51:54  ivos
      !  Added error reporting using ppm_error facility.
      !
      !  Revision 1.6  2004/01/13 13:04:09  ivos
      !  Added comment header.
      !
      !  Revision 1.5  2004/01/13 12:18:19  ivos
      !  Added initalization of IO units (call to ppm_set_unit).
      !
      !  Revision 1.4  2003/12/16 13:37:39  ivos
      !  Changed initial value of ppm_topoid from 0 to 1 (otherwise the first
      !  mapping failed because ppm_topoid is only incremented after the send.
      !
      !  Revision 1.3  2003/12/12 16:14:54  ivos
      !  Added initialization for internal topology list counters.
      !
      !  Revision 1.2  2003/12/11 17:28:51  hiebers
      !  Initialized ppm_nproc=1, ppm_rank=0 in serial version
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_init(dim,prec,tolexp,comm,debug,info,logfile,stderr,stdout)
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_log
      USE ppm_module_alloc
      USE ppm_module_print_defines
      USE ppm_module_io_set_unit
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif 
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN)   :: dim,prec,comm,debug,tolexp
      INTEGER, OPTIONAL     :: logfile,stderr,stdout
      INTEGER, INTENT(OUT)  :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER               :: iopt,ilen,istdout,istderr,ilog
      INTEGER, DIMENSION(1) :: ldl,ldu
      REAL(ppm_kind_double) :: t0
      LOGICAL               :: ppm_mpi_init
      CHARACTER(LEN=ppm_char) :: mesg
      CHARACTER(LEN=255)    :: cbuf
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
     
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_init',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_multipleinit,'ppm_init',  &
     &            'Please call ppm_finalize first!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Set unit numbers
      !-------------------------------------------------------------------------
      IF (PRESENT(stdout)) THEN
          istdout = stdout
      ELSE
          istdout = 6
      ENDIF
      IF (PRESENT(stderr)) THEN
          istderr = stderr
      ELSE
          istderr = 0
      ENDIF
      IF (PRESENT(logfile)) THEN
          ilog = logfile
      ELSE
          ! -1 means: do not write log file
          ilog = -1
      ENDIF
 
      !-------------------------------------------------------------------------
      !  Set unit numbers of stdout and stderr and log file
      !-------------------------------------------------------------------------
      CALL ppm_io_set_unit(istdout,istderr,ilog,info)

#if defined __MPI
      !-------------------------------------------------------------------------
      !  Check if MPI has been initialized
      !-------------------------------------------------------------------------
      CALL MPI_Initialized(ppm_mpi_init,info)

      !-------------------------------------------------------------------------
      !  If not exit with error status ppm_error_fatal
      !-------------------------------------------------------------------------
      IF (.NOT.ppm_mpi_init) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_nompi,'ppm_init',  &
     &         'Call MPI_Init before calling ppm_init !',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Get the MPI communicator
      !-------------------------------------------------------------------------
      ppm_comm = comm
      CALL MPI_Comm_Size(ppm_comm,ppm_nproc,info)
      CALL MPI_Comm_Rank(ppm_comm,ppm_rank,info)
      CALL MPI_Get_Processor_Name(cbuf,ilen,info)
      WRITE(mesg,'(2A)') '*** This is the PPM library starting on ',cbuf(1:ilen)
      CALL ppm_log('ppm_init',mesg,info)
      CALL ppm_write(ppm_rank,'ppm_init',mesg,info)
#else
      !SEH: Initialization of serial version
      ppm_nproc = 1
      ppm_rank  = 0
      WRITE(mesg,'(A)') '*** This is the PPM library in single-processor mode'
      CALL ppm_log('ppm_init',mesg,info)
      CALL ppm_write(ppm_rank,'ppm_init',mesg,info)
#endif

      !-------------------------------------------------------------------------
      !  Save the debugging option
      !-------------------------------------------------------------------------
      ppm_debug = debug
      WRITE(mesg,'(A,I2)') 'Debug level set to ',ppm_debug
      CALL ppm_log('ppm_init',mesg,info)
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_init',mesg,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Burp the defines of this library version
      !-------------------------------------------------------------------------
      CALL ppm_print_defines(info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_sub_failed,'ppm_init',     &
     &        'Print defines did not execute correctly',__LINE__,info)
          GOTO 9999
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Check and save the dimensionality of the problem
      !-------------------------------------------------------------------------
      IF     (dim.LT.2) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_wrong_dim,'ppm_init', &
     &         'Space dimension must be greater than 1!',__LINE__,info)
         GOTO 9999
      ELSEIF (dim.GT.3) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_wrong_dim,'ppm_init',  &
     &         'Space dimension must be less than 4!',__LINE__,info)
         GOTO 9999
      ELSE
         ppm_dim = dim
      ENDIF
      WRITE(mesg,'(A,I2)') 'Space dimension set to ',ppm_dim
      CALL ppm_log('ppm_init',mesg,info)
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_init',mesg,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that the user entered a valid precision
      !-------------------------------------------------------------------------
      IF     (prec.NE.ppm_kind_double.AND.prec.NE.ppm_kind_single) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_wrong_prec,'ppm_init', &
     &         'Must be either SINGE or DOUBLE!',__LINE__,info)
         GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  Save the precision
      !-------------------------------------------------------------------------
      ppm_kind = prec
      IF (ppm_kind .EQ. ppm_kind_single) THEN
          WRITE(mesg,'(A)') 'Precision set to SINGLE'
      ELSEIF (ppm_kind .EQ. ppm_kind_double) THEN
          WRITE(mesg,'(A)') 'Precision set to DOUBLE'
      ELSE
          WRITE(mesg,'(A)') 'Precision set to UNKNOWN'
      ENDIF
      CALL ppm_log('ppm_init',mesg,info)
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_init',mesg,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Save the tolerance
      !-------------------------------------------------------------------------
      IF (ppm_kind .EQ. ppm_kind_double) THEN
          IF     (tolexp.LT.-20.OR.tolexp.GT.-3) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_tol_warn,'ppm_init', &
     &           'Usual values are between 10^(-20) and 10^(-3)',__LINE__,info)
          ENDIF 
          IF     (tolexp.LT.INT(LOG10(EPSILON(ppm_myepsd)))) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_tol_warn,'ppm_init', &
     &           'Tolerance must not be smaller than machine epsilon'  &
     &           ,__LINE__,info)
          ENDIF 
          ppm_myepsd = 10.0_ppm_kind_double**REAL(tolexp,ppm_kind_double)
          ppm_myepss = REAL(ppm_myepsd,ppm_kind_single)
          WRITE(mesg,'(A,E15.10)') 'Floating point tolerance set to ',  &
     &        ppm_myepsd
      ELSE
          IF     (tolexp.LT.-12.OR.tolexp.GT.-3) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_tol_warn,'ppm_init', &
     &           'Usual values are between 10^(-12) and 10^(-3)',__LINE__,info)
          ENDIF 
          IF     (tolexp.LT.INT(LOG10(EPSILON(ppm_myepss)))) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_tol_warn,'ppm_init', &
     &           'Tolerance must not be smaller than machine epsilon'  &
     &           ,__LINE__,info)
          ENDIF 
          ppm_myepsd = 10.0_ppm_kind_double**REAL(tolexp,ppm_kind_double)
          ppm_myepss = REAL(ppm_myepsd,ppm_kind_single)
          WRITE(mesg,'(A,E15.10)') 'Floating point tolerance set to ',   &
     &        ppm_myepss
      ENDIF
      CALL ppm_log('ppm_init',mesg,info)
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_init',mesg,info)
      ENDIF

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Save the MPI precision 
      !-------------------------------------------------------------------------
      IF (ppm_kind.EQ.ppm_kind_double) THEN
         ppm_mpi_kind = MPI_DOUBLE_PRECISION
      ELSE
         ppm_mpi_kind = MPI_REAL
      ENDIF 
#endif 

      !-------------------------------------------------------------------------
      !  Definition of PI
      !-------------------------------------------------------------------------
      ppm_pi_d = DACOS(-1.0_ppm_kind_double)
      ppm_pi_s =  ACOS(-1.0_ppm_kind_single)

      !-------------------------------------------------------------------------
      !  Allocate and initialize the processor speed array
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldl(1) = 0
      ldu(1) = ppm_nproc-1
      CALL ppm_alloc(ppm_proc_speed,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_init',  &
     &        'processor speed array PPM_PROC_SPEED',__LINE__,info)
          GOTO 9999
      ENDIF
      ! initially assume all processors to be equally fast
      ppm_proc_speed(0:ppm_nproc-1) =    &
     &   1.0_ppm_kind_double/REAL(ppm_nproc,ppm_kind_double)

      !-------------------------------------------------------------------------
      !  Reset the topology counter
      !-------------------------------------------------------------------------
      ppm_max_topoid   = 0
      ppm_topoid       = -1
      ppm_field_topoid = -1

      !-------------------------------------------------------------------------
      !  Set the global status to initialized
      !-------------------------------------------------------------------------
      ppm_initialized = .TRUE.

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_init',t0,info)
      RETURN
      END SUBROUTINE ppm_init

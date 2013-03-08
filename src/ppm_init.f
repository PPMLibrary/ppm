      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_init
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

      SUBROUTINE ppm_init(dim,prec,tolexp,comm,debug,info,logfile,stderr,stdout)
      !!! Initialisation of the ppm library.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_topo_typedef
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_log
      USE ppm_module_alloc
      USE ppm_module_print_defines
      USE ppm_module_io
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN)   :: dim
      !!! dimension of the problem (2 or 3)
      INTEGER, INTENT(IN)   :: prec
      !!! desired precision  ppm_kind_double or ppm_kind_single
      INTEGER, INTENT(IN)   :: comm
      !!! MPI communicator
      INTEGER, INTENT(IN)   :: debug
      !!! The debug level (0,1 or 2)
      INTEGER, INTENT(IN)   :: tolexp
      !!! log10 of the numerical tolerance for real comparisons (> machine eps).
      !!! E.g. for a tolerance of 10^(-14), tolexp = -14. Everything .LT.
      !!! 10^tolexp will be considered zero (0) and its inversion is considered
      !!! dangerous.
      INTEGER, OPTIONAL     :: logfile
      !!! OPTIONAL. Unit number where to print log messages to.
      INTEGER, OPTIONAL     :: stderr
      !!! OPTIONAL. Unit number where to print the stderr to
      INTEGER, OPTIONAL     :: stdout
      !!! OPTIONAL. Unit number where to print the stdout to
      INTEGER, INTENT(OUT)  :: info
      !!! returns status 0 upon success
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
      !  Save the debugging option
      !-------------------------------------------------------------------------
      ppm_debug = debug

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
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

#if defined __MPI
      !-------------------------------------------------------------------------
      !  Check if MPI has been initialized
      !-------------------------------------------------------------------------
      CALL MPI_Initialized(ppm_mpi_init,info)

      !-------------------------------------------------------------------------
      !  If MPI not initialized exit with error status ppm_error_fatal
      !-------------------------------------------------------------------------
      IF (.NOT.ppm_mpi_init) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_nompi,'ppm_init',  &
         &   'Call MPI_Init before calling ppm_init !',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Get the MPI communicator
      !-------------------------------------------------------------------------
      ppm_comm = comm
      CALL MPI_Comm_Size(ppm_comm,ppm_nproc,info)
      CALL MPI_Comm_Rank(ppm_comm,ppm_rank,info)
      CALL MPI_Get_Processor_Name(cbuf,ilen,info)
#else
      ppm_nproc = 1
      ppm_rank  = 0
#endif

      !-------------------------------------------------------------------------
      !  Set unit numbers of stdout and stderr and log file
      !-------------------------------------------------------------------------
      CALL ppm_io_set_unit(istdout,istderr,ilog,info)

      !-------------------------------------------------------------------------
      !  Now the library can talk...
      !-------------------------------------------------------------------------

      CALL substart('ppm_init',t0,info)

#if defined __MPI
      IF (ppm_debug .GT. 0) THEN
        WRITE(mesg,'(2A)') '*** This is the PPM library starting on ',&
     &                    cbuf(1:ilen)
        CALL ppm_log('ppm_init',mesg,info)
        CALL ppm_write(ppm_rank,'ppm_init',mesg,info)
      ELSE
        IF (ppm_rank .EQ. 0) THEN
            WRITE(mesg,'(2A)') '*** This is the PPM library starting on ',&
     &                    cbuf(1:ilen)
            CALL ppm_log('ppm_init',mesg,info)
            CALL ppm_write(ppm_rank,'ppm_init',mesg,info)
        ENDIF
      ENDIF
#else
      WRITE(mesg,'(A)') '*** This is the PPM library in single-processor mode'
      CALL ppm_log('ppm_init',mesg,info)
      CALL ppm_write(ppm_rank,'ppm_init',mesg,info)
#endif

      !-------------------------------------------------------------------------
      !  Output the debugging option
      !-------------------------------------------------------------------------
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
             GOTO 9999
          ENDIF
          ppm_myepsd = 10.0_ppm_kind_double**REAL(tolexp,ppm_kind_double)
          ppm_myepss = REAL(ppm_myepsd,ppm_kind_single)
          WRITE(mesg,'(A,E17.10)') 'Floating point tolerance set to ',  &
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
             GOTO 9999
          ENDIF
          ppm_myepsd = 10.0_ppm_kind_double**REAL(tolexp,ppm_kind_double)
          ppm_myepss = REAL(ppm_myepsd,ppm_kind_single)
          WRITE(mesg,'(A,E17.10)') 'Floating point tolerance set to ',   &
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
      ppm_next_avail_topo = ppm_param_undefined

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
      CONTAINS
      SUBROUTINE check
        IF (ppm_initialized) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_multipleinit,'ppm_init',  &
     &          'Please call ppm_finalize first!',__LINE__,info)
            GOTO 8888
        ENDIF
 8888   CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_init

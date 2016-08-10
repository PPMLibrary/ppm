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
      USE ppm_module_data, ONLY: ppm_kind,ppm_kind_single,ppm_kind_double,     &
      & ppm_char,ppm_debug,ppm_comm,ppm_nproc,ppm_rank,ppm_dim,ppm_proc_speed, &
      & ppm_initialized,ppm_pi_s,ppm_pi_d,ppm_param_undefined,ppm_mpi_kind,    &
      & ppm_param_alloc_fit,ppm_myepss,ppm_myepsd,                             &
      & ppm_error_fatal,ppm_error_warning,ppm_error_error
      USE ppm_module_topo_typedef, ONLY: ppm_next_avail_topo
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error, ONLY: ppm_error,ppm_err_sub_failed,&
      & ppm_err_wrong_dim,ppm_err_wrong_prec,ppm_err_tol_warn,ppm_err_alloc,&
      & ppm_err_multipleinit,ppm_err_nompi,ppm_err_mpi_fail
      USE ppm_module_write
      USE ppm_module_log
      USE ppm_module_alloc, ONLY: ppm_alloc
      USE ppm_module_mpi
      USE ppm_module_print_defines
      USE ppm_module_io, ONLY: ppm_io_set_unit
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER,           INTENT(IN   ) :: dim
      !!! dimension of the problem (2 or 3)
      INTEGER,           INTENT(IN   ) :: prec
      !!! desired precision  ppm_kind_double or ppm_kind_single
      INTEGER,           INTENT(IN   ) :: comm
      !!! MPI communicator
      INTEGER,           INTENT(IN   ) :: debug
      !!! The debug level (0,1 or 2)
      INTEGER,           INTENT(IN   ) :: tolexp
      !!! log10 of the numerical tolerance for real comparisons (> machine eps).
      !!! E.g. for a tolerance of 10^(-14), tolexp = -14. Everything .LT.
      !!! 10^tolexp will be considered zero (0) and its inversion is considered
      !!! dangerous.
      INTEGER, OPTIONAL, INTENT(IN   ) :: logfile
      !!! OPTIONAL. Unit number where to print log messages to.
      INTEGER, OPTIONAL, INTENT(IN   ) :: stderr
      !!! OPTIONAL. Unit number where to print the stderr to
      INTEGER, OPTIONAL, INTENT(IN   ) :: stdout
      !!! OPTIONAL. Unit number where to print the stdout to
      INTEGER,           INTENT(  OUT) :: info
      !!! returns status 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER               :: iopt,ilen
      INTEGER               :: istdout,istderr,ilog
      INTEGER, DIMENSION(1) :: ldl,ldu

      CHARACTER(LEN=ppm_char) :: mesg
      CHARACTER(LEN=ppm_char) :: caller='ppm_init'

      LOGICAL :: ppm_mpi_init
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
      istdout=MERGE(stdout,6,PRESENT(stdout))
      istderr=MERGE(stderr,0,PRESENT(stderr))
      ! -1 means: do not write log file
      ilog   =MERGE(logfile,-1,PRESENT(logfile))

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         CALL check
         IF (info.NE.0) GOTO 9999
      ENDIF

#if defined __MPI
      !-------------------------------------------------------------------------
      !  Check if MPI has been initialized
      !-------------------------------------------------------------------------
      CALL MPI_Initialized(ppm_mpi_init,info)
      or_fail_MPI("MPI_Initialized")

      !-------------------------------------------------------------------------
      !  If MPI not initialized exit with error status ppm_error_fatal
      !-------------------------------------------------------------------------
      IF (.NOT.ppm_mpi_init) THEN
         fail('Call MPI_Init before calling ppm_init !',ppm_err_nompi,ppm_error=ppm_error_fatal)
      ENDIF

      !-------------------------------------------------------------------------
      !  Get the MPI communicator
      !-------------------------------------------------------------------------
      CALL MPI_Comm_dup(comm,ppm_comm,info)
      or_fail_MPI("MPI_Comm_dup")

      CALL MPI_Comm_size(ppm_comm,ppm_nproc,info)
      or_fail_MPI("MPI_Comm_size")

      CALL MPI_Comm_rank(ppm_comm,ppm_rank,info)
      or_fail_MPI("MPI_Comm_rank")

      CALL MPI_Get_processor_name(mesg,ilen,info)
      or_fail_MPI("MPI_Get_processor_name")
#else
      ppm_nproc = 1
      ppm_rank  = 0
#endif

      !-------------------------------------------------------------------------
      !  Set unit numbers of stdout and stderr and log file
      !-------------------------------------------------------------------------
      CALL ppm_io_set_unit(istdout,istderr,ilog,info)
      or_fail("ppm_io_set_unit")

      !-------------------------------------------------------------------------
      !  Now the library can talk...
      !-------------------------------------------------------------------------

      CALL substart(caller,t0,info)

#if defined __MPI
      IF (ppm_debug.GT.0) THEN
         stdout("*** This is the PPM library starting on ",'mesg(1:ilen)')
         CALL ppm_log(caller,cbuf,info)
      ELSE
        IF (ppm_rank.EQ.0) THEN
           stdout_f('(A)',"**************************************************************")
           CALL ppm_log(caller,cbuf,info)
           stdout_f('(A)',"***       Parallel Particle Mesh Library (PPM)       *********")
           CALL ppm_log(caller,cbuf,info)
           stdout_f('(A)',"***          version:  1.2.2  /  May 2016            ****___**")
           CALL ppm_log(caller,cbuf,info)
           stdout_f('(A,i12,A)',"***       Execution on  / ",ppm_nproc,"    node(s)    ***/ _ \*")
           CALL ppm_log(caller,cbuf,info)
           stdout_f('(A)',"***                                                  **/ .__/*")
           CALL ppm_log(caller,cbuf,info)
           stdout_f('(A)',"*** Contributors                                     */_/___**")
           CALL ppm_log(caller,cbuf,info)
           stdout_f('(A)',"*** ------------                                     ***/ _ \*")
           CALL ppm_log(caller,cbuf,info)
           stdout_f('(A)',"*** CSE-lab at ETH Zurich (group of Professor Petros **/ .__/*")
           CALL ppm_log(caller,cbuf,info)
           stdout_f('(A)',"*** Koumoutsakos), MOSAIC group at MPI-CBG Dresden   */_/___**")
           CALL ppm_log(caller,cbuf,info)
           stdout_f('(A)',"*** (group of Professor Ivo F.Sbalzarini)            **/  ' \*")
           CALL ppm_log(caller,cbuf,info)
           stdout_f('(A)',"*** and Center for Fluid Dynamics at DTU (group of   */_/_/_/*")
           CALL ppm_log(caller,cbuf,info)
           stdout_f('(A)',"*** Professor Jens Walther).                         *********")
           CALL ppm_log(caller,cbuf,info)
           stdout_f('(A)',"**************************************************************")
           CALL ppm_log(caller,cbuf,info)
        ENDIF
      ENDIF
#else
      stdout_f('(A)',"**************************************************************")
      CALL ppm_log(caller,cbuf,info)
      stdout_f('(A)',"***       Parallel Particle Mesh Library (PPM)       *********")
      CALL ppm_log(caller,cbuf,info)
      stdout_f('(A)',"***          version:  1.2.2  /  May 2016            ****___**")
      CALL ppm_log(caller,cbuf,info)
      stdout_f('(A)',"***        Execution in single-processor mode        ***/ _ \*")
      CALL ppm_log(caller,cbuf,info)
      stdout_f('(A)',"***                                                  **/ .__/*")
      CALL ppm_log(caller,cbuf,info)
      stdout_f('(A)',"*** Contributors                                     */_/___**")
      CALL ppm_log(caller,cbuf,info)
      stdout_f('(A)',"*** ------------                                     ***/ _ \*")
      CALL ppm_log(caller,cbuf,info)
      stdout_f('(A)',"*** CSE-lab at ETH Zurich (group of Professor Petros **/ .__/*")
      CALL ppm_log(caller,cbuf,info)
      stdout_f('(A)',"*** Koumoutsakos), MOSAIC group at MPI-CBG Dresden   */_/___**")
      CALL ppm_log(caller,cbuf,info)
      stdout_f('(A)',"*** (group of Professor Ivo F.Sbalzarini)            **/  ' \*")
      CALL ppm_log(caller,cbuf,info)
      stdout_f('(A)',"*** and Center for Fluid Dynamics at DTU (group of   */_/_/_/*")
      CALL ppm_log(caller,cbuf,info)
      stdout_f('(A)',"*** Professor Jens Walther).                         *********")
      CALL ppm_log(caller,cbuf,info)
      stdout_f('(A)',"**************************************************************")
      CALL ppm_log(caller,cbuf,info)
#endif

      !-------------------------------------------------------------------------
      !  Output the debugging option
      !-------------------------------------------------------------------------
      WRITE(cbuf,'(A,I2)') 'Debug level set to ',ppm_debug
      CALL ppm_log(caller,cbuf,info)
      IF (ppm_debug.GT.0) THEN
         CALL ppm_write(ppm_rank,caller,cbuf,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Burp the defines of this library version
      !-------------------------------------------------------------------------
      CALL ppm_print_defines(info)
      or_fail('Print defines did not execute correctly',ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Check and save the dimensionality of the problem
      !-------------------------------------------------------------------------
      IF     (dim.LT.2) THEN
         fail('Space dimension must be greater than 1!',ppm_err_wrong_dim,ppm_error=ppm_error_fatal)
      ELSEIF (dim.GT.3) THEN
         fail('Space dimension must be less than 4!',ppm_err_wrong_dim,ppm_error=ppm_error_fatal)
      ELSE
         ppm_dim = dim
      ENDIF
      WRITE(cbuf,'(A,I2)') 'Space dimension set to ',ppm_dim
      CALL ppm_log(caller,cbuf,info)
      IF (ppm_debug.GT.0) THEN
         CALL ppm_write(ppm_rank,caller,cbuf,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that the user entered a valid precision
      !-------------------------------------------------------------------------
      IF (prec.NE.ppm_kind_double.AND.prec.NE.ppm_kind_single) THEN
         fail('Must be either SINGE or DOUBLE!',ppm_err_wrong_prec,ppm_error=ppm_error_fatal)
      ENDIF

      !-------------------------------------------------------------------------
      !  Save the precision
      !-------------------------------------------------------------------------
      ppm_kind = prec
      SELECT CASE (ppm_kind)
      CASE (ppm_kind_single)
         WRITE(cbuf,'(A)') 'Precision set to SINGLE'

      CASE (ppm_kind_double)
         WRITE(cbuf,'(A)') 'Precision set to DOUBLE'

      CASE DEFAULT
         WRITE(cbuf,'(A)') 'Precision set to UNKNOWN'

      END SELECT
      CALL ppm_log(caller,cbuf,info)
      IF (ppm_debug.GT.0) THEN
         CALL ppm_write(ppm_rank,caller,cbuf,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Save the tolerance
      !-------------------------------------------------------------------------
      SELECT CASE (ppm_kind)
      CASE (ppm_kind_double)
         IF (tolexp.LT.-20.OR.tolexp.GT.-3) THEN
            info = ppm_error_warning
            ppm_fail('Usual values are between 10^(-20) and 10^(-3)',ppm_err_tol_warn,exit_point=no)
         ENDIF
         IF (tolexp.LT.INT(LOG10(EPSILON(ppm_myepsd)))) THEN
            fail('Tolerance must not be smaller than machine epsilon',ppm_err_tol_warn,ppm_error=ppm_error_fatal)
         ENDIF
         ppm_myepsd = 10.0_ppm_kind_double**tolexp
         ppm_myepss = REAL(ppm_myepsd,ppm_kind_single)
         WRITE(cbuf,'(A,E17.10)') 'Floating point tolerance set to ',ppm_myepsd

      CASE DEFAULT
         IF (tolexp.LT.-12.OR.tolexp.GT.-3) THEN
            info = ppm_error_warning
            ppm_fail('Usual values are between 10^(-12) and 10^(-3)',ppm_err_tol_warn,exit_point=no)
         ENDIF
         IF (tolexp.LT.INT(LOG10(EPSILON(ppm_myepss)))) THEN
            fail('Tolerance must not be smaller than machine epsilon',ppm_err_tol_warn,ppm_error=ppm_error_fatal)
         ENDIF
         ppm_myepsd = 10.0_ppm_kind_double**tolexp
         ppm_myepss = REAL(ppm_myepsd,ppm_kind_single)
         WRITE(cbuf,'(A,E17.10)') 'Floating point tolerance set to ',ppm_myepss

      END SELECT
      CALL ppm_log(caller,cbuf,info)
      IF (ppm_debug.GT.0) THEN
         CALL ppm_write(ppm_rank,caller,cbuf,info)
      ENDIF

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Save the MPI precision
      !-------------------------------------------------------------------------
      ppm_mpi_kind = MERGE(MPI_DOUBLE_PRECISION,MPI_REAL,ppm_kind.EQ.ppm_kind_double)
#endif

      !-------------------------------------------------------------------------
      !  Definition of PI
      !-------------------------------------------------------------------------
      ppm_pi_d = DACOS(-1.0_ppm_kind_double)
      ppm_pi_s =  ACOS(-1.0_ppm_kind_single)

      !-------------------------------------------------------------------------
      !  Allocate and initialize the processor speed array
      !-------------------------------------------------------------------------
      iopt  =ppm_param_alloc_fit
      ldl(1)=0
      ldu(1)=ppm_nproc-1
      CALL ppm_alloc(ppm_proc_speed,ldl,ldu,iopt,info)
      or_fail_alloc("processor speed array PPM_PROC_SPEED",ppm_error=ppm_error_fatal)

      ! initially assume all processors to be equally fast
      ppm_proc_speed(0:ppm_nproc-1)=1.0_ppm_kind_double/REAL(ppm_nproc,ppm_kind_double)

      !-------------------------------------------------------------------------
      !  Reset the topology counter
      !-------------------------------------------------------------------------
      ppm_next_avail_topo=ppm_param_undefined

      !-------------------------------------------------------------------------
      !  Set the global status to initialized
      !-------------------------------------------------------------------------
      ppm_initialized=.TRUE.

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
        SUBROUTINE check
          IF (ppm_initialized) THEN
             fail('Please call ppm_finalize first!',ppm_err_multipleinit,exit_point=8888)
          ENDIF
      8888 CONTINUE
        END SUBROUTINE check
      END SUBROUTINE ppm_init

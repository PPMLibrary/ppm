      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_loadbal_inquire
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Lab (ETH Zurich), 
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

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE loadbal_inq_s(ctime,nstep,heuristic,lflush,imbal, &
     &    lredecomp,nredest,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE loadbal_inq_d(ctime,nstep,heuristic,lflush,imbal, &
     &    lredecomp,nredest,info)
#endif
      !!! Inquires about the load balance status and returns advise whether
      !!! redecomposing the problem is recommended (based on some decision
      !!! heuristic the user chooses).
      !!!
      !!! [TIP]
      !!! The user should time (using ppm_time) computations
      !!! for the topology/topologies considered for dyanmic
      !!! remaping. The elapsed time is given to this
      !!! routine.
      !!!
      !!! [NOTE]
      !!! This estimate is currently a scalar, so only one topology/set of
      !!! topologies can be monitored. Introducing a topology set ID and making
      !!! the internal estimates a vector, this can later be extended to
      !!! multiple topology sets if needed. The topology set
      !!! ID for which to evaluate the load balance will then be an additional
      !!! argument to this routine. This routine does two global MPI
      !!! communication operations (`MPI_Allreduce`).
      !!!
      !!! .References
      !!! *************************************************************
      !!! - B. Moon and J. Saltz, Adaptive Runtime
      !!!   Support for Direct Simulation Monte Carlo Methods on
      !!!   Distributed Memory Architectures. Proceedings of
      !!!   the IEEE Scalable High-Performance Computing
      !!!   Conference. 1994. 176--183.
      !!! *************************************************************
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_loadbal
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
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
      REAL(MK)               , INTENT(IN   ) :: ctime
      !!! Elapsed time (as measured by `ppm_time`) for all computation in one
      !!! time step on the local processor.
      INTEGER                , INTENT(IN   ) :: heuristic
      !!! Decision heuristic for redecomposition advise. One of:
      !!! 
      !!! * ppm_param_loadbal_sar (Stop-at-Rise heuristic)
      INTEGER                , INTENT(IN   ) :: nstep
      !!! Number of time steps since last redecomposition. (>0)
      !!! If this routine is not called every time step, linear interpolation of
      !!! the load imbalance will be used to reconstruct missing data points.
      LOGICAL                , INTENT(IN   ) :: lflush
      !!! `TRUE` to flush internal statistics (e.g. the first time this routine
      !!! is called after actually doing a redecomposition of the problem),
      !!! `FALSE` to continue gathering statistics.
      REAL(MK)               , INTENT(  OUT) :: imbal
      !!! Load imbalance defined as the ratio of computation time of the
      !!! bottleneck processor to the average computation time of all processors.
      LOGICAL                , INTENT(  OUT) :: lredecomp
      !!! `TRUE` if the choosen heuristic recommends problem redecomposition.
      !!! Otherwise `FALSE`. Redecomposition means: do the ppm_topo_mktopo again.
      INTEGER                , INTENT(  OUT) :: nredest
      !!! Estimated (linear extrapolation) number of time steps to go until next
      !!! advised redecomposition. -1 is returned if the chosen heuristic
      !!! does not support this kind of information. Be careful with this value!
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)               :: t0,maxtime,avgtime
      REAL(ppm_kind_double)  :: Wn,delta,slp,estsol
      REAL(ppm_kind_double), SAVE :: deltaold = 0.0_ppm_kind_double
      REAL(ppm_kind_double), SAVE :: slpavg = 0.0_ppm_kind_double
      INTEGER, SAVE          :: nold = 0
      INTEGER, SAVE          :: nslp = 0
      CHARACTER(LEN=ppm_char):: mesg
      INTEGER                :: ndiff,i
#ifdef __MPI
      INTEGER                :: MPTYPE
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_loadbal_inquire',t0,info)
      imbal = 0.0_MK
      lredecomp = .FALSE.
      nredest = -1 

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Flush internal information if needed
      !-------------------------------------------------------------------------
      IF (lflush) THEN
          ppm_loadbal_old_sar = -1.0_ppm_kind_double
          ppm_loadbal_runsum  = 0.0_ppm_kind_double
          slpavg  = 0.0_ppm_kind_double
          nold = 0
          nslp = 0
          deltaold = 0.0_ppm_kind_double
          IF (ppm_debug .GT. 1) THEN
              CALL ppm_write(ppm_rank,'ppm_loadbal_inquire',    &
     &            'Internal buffers flushed',info)
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that we have data available
      !-------------------------------------------------------------------------
      IF (ppm_loadbal_dcn .LT. 1) THEN
          info = ppm_error_warning
          CALL ppm_error(ppm_err_no_data,'ppm_loadbal_inquire', &
     &         'Use ppm_set_decomp_cost first to gather statistics',  &
     &         __LINE__,info)
          GOTO 9999
      ENDIF

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Define MPI data type
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      MPTYPE = MPI_REAL
#elif __KIND == __DOUBLE_PRECISION
      MPTYPE = MPI_DOUBLE_PRECISION
#endif
#endif

      !-------------------------------------------------------------------------
      !  Determine the time used by the slowest processor
      !-------------------------------------------------------------------------
#ifdef __MPI
      CALL MPI_Allreduce(ctime,maxtime,1,MPTYPE,MPI_MAX,ppm_comm,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_mpi_fail,'ppm_loadbal_inquire', &
     &         'MPI_ALLREDUCE for max ctime. Nothing computed.',__LINE__,info)
          GOTO 9999
      ENDIF
#else
      maxtime = ctime
#endif

      !-------------------------------------------------------------------------
      !  Determine the average compute time
      !-------------------------------------------------------------------------
#ifdef __MPI
      CALL MPI_Allreduce(ctime,avgtime,1,MPTYPE,MPI_SUM,ppm_comm,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_mpi_fail,'ppm_loadbal_inquire', &
     &         'MPI_ALLREDUCE for avg ctime. Nothing computed.',__LINE__,info)
          GOTO 9999
      ENDIF
      avgtime = avgtime/REAL(ppm_nproc,MK)
#else
      avgtime = ctime
#endif

      !-------------------------------------------------------------------------
      !  Debug output
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(2(A,F12.6))') 'max time = ',maxtime,   &
     &        ' / average time = ',avgtime
          CALL ppm_write(ppm_rank,'ppm_loadbal_inquire',mesg,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute load imbalance
      !-------------------------------------------------------------------------
      imbal = maxtime/avgtime

      !-------------------------------------------------------------------------
      !  Choose decision heuristic
      !-------------------------------------------------------------------------
      IF     (heuristic .EQ. ppm_param_loadbal_sar) THEN
          !---------------------------------------------------------------------
          !  Stop-At-Rise heuristic by Moon:1994
          !---------------------------------------------------------------------
          ndiff = nstep-nold
#if   __KIND == __SINGLE_PRECISION
          delta = REAL(maxtime-avgtime,ppm_kind_double)
#elif __KIND == __DOUBLE_PRECISION
          delta = maxtime-avgtime
#endif
          slp = (delta-deltaold)/REAL(ndiff,ppm_kind_double)
          IF (ndiff .GT. 1) THEN
              !-----------------------------------------------------------------
              !  Use linear interpolation for missing data points
              !-----------------------------------------------------------------
              IF (ppm_debug .GT. 1) THEN
                  WRITE(mesg,'(A,I3)') 'Linear interpolation steps: ',   &
     &                ndiff
                  CALL ppm_write(ppm_rank,'ppm_loadbal_inquire',mesg,info)
              ENDIF
              DO i=1,ndiff
                  ppm_loadbal_runsum = ppm_loadbal_runsum +    &
     &                (deltaold+(REAL(i,ppm_kind_double)*slp))
              ENDDO
          ELSEIF (ndiff .EQ. 1) THEN
              ppm_loadbal_runsum = ppm_loadbal_runsum + delta
          ELSE
              !-----------------------------------------------------------------
              !  Called twice for the same time step or going backward in
              !  time: warn
              !-----------------------------------------------------------------
              info = ppm_error_warning
              CALL ppm_error(ppm_err_rev_time,'ppm_loadbal_inquire',    &
         &       'nstep did not increase since last call!',__LINE__,info)
          ENDIF
          Wn = ppm_loadbal_runsum + ppm_loadbal_decomp_cost
          Wn = Wn/REAL(nstep,ppm_kind_double)
          ! if we have an old value at all...
          IF (ppm_loadbal_old_sar .GT. -0.5_ppm_kind_double) THEN
              IF (ppm_loadbal_old_sar .LT. Wn) THEN
                  lredecomp = .TRUE.
                  nredest   = 0
                  IF (ppm_debug .GT. 1) THEN
                      CALL ppm_write(ppm_rank,'ppm_loadbal_inquire',   &
     &                    'redecomposition advised by SAR heuristic.',info)
                  ENDIF
              ENDIF
          ENDIF
          IF (ppm_debug .GT. 1) THEN
              WRITE(mesg,'(2(A,F12.6))') 'old SAR value: ',   &
     &            ppm_loadbal_old_sar,' / new SAR value: ',Wn
              CALL ppm_write(ppm_rank,'ppm_loadbal_inquire',mesg,info)
          ENDIF
          !---------------------------------------------------------------------
          !  Update running average (with exp. forgetting) of slope
          !---------------------------------------------------------------------
          nslp = nslp + 1
          IF (nslp .EQ. 1) THEN
              slpavg = slp
          ELSE
              slpavg = 0.5_ppm_kind_double*(slpavg+slp)
          ENDIF
          !---------------------------------------------------------------------
          !  Extrapolate expected number of steps until redecomposition
          !---------------------------------------------------------------------
          estsol = 2.0_ppm_kind_double*ppm_loadbal_decomp_cost/slpavg
          IF (estsol .LT. 0.0_ppm_kind_double) THEN
              nredest = -1
              info = ppm_error_notice
              CALL ppm_error(ppm_err_sqrt_neg,'ppm_loadbal_inquire', &
     &            'Balance improved. Cannot extrapolate load imbalance.', &
     &            __LINE__,info)
          ELSE
              nredest = INT(SQRT(estsol)) - nstep
              IF (nredest .LT. 0) nredest = 0
          ENDIF

          !---------------------------------------------------------------------
          !  Update old values
          !---------------------------------------------------------------------
          ppm_loadbal_old_sar = Wn
          deltaold = delta
      ELSE
          !---------------------------------------------------------------------
          !  Unknow heuristic
          !---------------------------------------------------------------------
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,'ppm_loadbal_inquire',    &
     &       'Unknown heuristic specified. Nothing computed.',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Update time step counter
      !-------------------------------------------------------------------------
      nold = nstep

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_loadbal_inquire',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_loadbal_inquire',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
        ENDIF
        IF (ctime .LT. 0.0_MK) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_loadbal_inquire', &
     &          'ctime must be >= 0.0',__LINE__,info)
           GOTO 8888
        ENDIF
        IF (nstep .LT. 1) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_loadbal_inquire', &
     &           'nstep must be > 0',__LINE__,info)
            GOTO 8888
        ENDIF
        IF (nstep .LE. nold) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_loadbal_inquire', &
     &            'nstep did not increase since last call! Use flush?', &
     &            __LINE__,info)
           GOTO 8888
        ENDIF
 8888   CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE loadbal_inq_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE loadbal_inq_d
#endif

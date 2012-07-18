      minclude ppm_header(ppm_module_name="loadbal_inquire_sar")

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE loadbal_inquire_sar_s(ctime,nstep,lflush,lredecomp,nredest, &
     &                                 info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE loadbal_inquire_sar_d(ctime,nstep,lflush,lredecomp,nredest, &
     &                                 info)
#endif
      !!! This routine inquires about the load balance status using the SAR
      !!! decision policy and returns advise whether redecomposing the problem
      !!! domain is recommended.
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
#ifdef __MPI
       INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK)               , INTENT(IN   ) :: ctime
      !!! Elapsed time (as measured by `ppm_time`) for all computations in one
      !!! time step on the local processor.
      INTEGER                , INTENT(IN   ) :: nstep
      !!! Number of time steps since last redecomposition. (>0)
      !!! If this routine is not called every time step, linear interpolation of
      !!! the load imbalance will be used to reconstruct missing data points.
      LOGICAL                , INTENT(IN   ) :: lflush
      !!! `TRUE` to flush internal statistics (e.g. the first time this routine
      !!! is called after actually doing a redecomposition of the problem),
      !!! `FALSE` to continue gathering statistics.
      LOGICAL                , INTENT(INOUT) :: lredecomp
      !!! `TRUE` if the choosen heuristic recommends problem redecomposition.
      !!! Otherwise `FALSE`. Redecomposition means: do the ppm_topo_mktopo again.
      INTEGER                , INTENT(INOUT) :: nredest
      !!! Estimated (linear extrapolation) number of time steps to go until next
      !!! advised redecomposition. -1 is returned if the chosen heuristic
      !!! does not support this kind of information. Be careful with this value!
      !!! DLB will always return 1 (i.e. next step is the DLB step)
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)               :: maxtime,avgtime
      REAL(ppm_kind_double)  :: Wn,delta,slp,estsol
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
      start_subroutine("loadbal_inquire_sar")
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
          ppm_loadbal_slpavg  = 0.0_ppm_kind_double
          ppm_loadbal_nold = 0
          ppm_loadbal_nslp = 0
          ppm_loadbal_deltaold = 0.0_ppm_kind_double
          IF (ppm_debug .GT. 1) THEN
            stdout("Internal buffers flushed")
          ENDIF
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
      !  TODO: Try to avoid global comm
      !-------------------------------------------------------------------------
#ifdef __MPI
      CALL MPI_Allreduce(ctime,maxtime,1,MPTYPE,MPI_MAX,ppm_comm,info)
           or_fail("MPI_ALLREDUCE for max ctime. Nothing computed.")
#else
      maxtime = ctime
#endif

      !-------------------------------------------------------------------------
      !  Determine the average compute time
      !  TODO: Try to avoid global comm
      !-------------------------------------------------------------------------
#ifdef __MPI
      CALL MPI_Allreduce(ctime,avgtime,1,MPTYPE,MPI_SUM,ppm_comm,info)
           or_fail("MPI_ALLREDUCE for avgtime. Nothing computed.")

      avgtime = avgtime/REAL(ppm_nproc,MK)
#else
      avgtime = ctime
#endif

      !-------------------------------------------------------------------------
      !  Debug output
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
         stdout_f('(2(A,F12.6))',"max time=",maxtime,"/average time =",avgtime)
      ENDIF
      !-------------------------------------------------------------------
      !  Stop-At-Rise heuristic by Moon:1994
      !-------------------------------------------------------------------
      ndiff = nstep-ppm_loadbal_nold
#if   __KIND == __SINGLE_PRECISION
      delta = REAL(maxtime-avgtime,ppm_kind_double)
#elif __KIND == __DOUBLE_PRECISION
      delta = maxtime-avgtime
#endif
      slp = (delta-ppm_loadbal_deltaold)/REAL(ndiff,ppm_kind_double)
      IF (ndiff .GT. 1) THEN
        !---------------------------------------------------------------
        !  Use linear interpolation for missing data points
        !---------------------------------------------------------------
        IF (ppm_debug .GT. 1) THEN
!            stdout_f('(A,I3)',"Linear interpolation steps: ",ndiff)
        ENDIF
        DO i=1,ndiff
            ppm_loadbal_runsum = ppm_loadbal_runsum +    &
     &              (ppm_loadbal_deltaold+(REAL(i,ppm_kind_double)*slp))
        ENDDO
      ELSEIF (ndiff .EQ. 1) THEN
        ppm_loadbal_runsum = ppm_loadbal_runsum + delta
      ELSE
        !---------------------------------------------------------------
        !  Called twice for the same time step or going backward in
        !  time: warn
        !---------------------------------------------------------------
        fail("nstep did not increase since last call!")

      ENDIF
      Wn = ppm_loadbal_runsum + ppm_loadbal_decomp_cost
      Wn = Wn/REAL(nstep,ppm_kind_double)
      ! if we have an old value at all...
      IF (ppm_loadbal_old_sar .GT. -0.5_ppm_kind_double) THEN
        IF (ppm_loadbal_old_sar .LT. Wn) THEN
            lredecomp = .TRUE.
            nredest   = 0
            IF (ppm_debug .GT. 1) THEN
                CALL ppm_write(ppm_rank,'loadbal_inquire_sar',   &
     &              'redecomposition advised by SAR heuristic.',info)
            ENDIF
        ENDIF
      ENDIF
      IF (ppm_debug .GT. 1) THEN
        stdout_f('(2(A,F12.6))',"old SAR value: ",  &
     &           ppm_loadbal_old_sar," / new SAR value: ",Wn)
      ENDIF
      !-------------------------------------------------------------------
      !  Update running average (with exp. forgetting) of slope
      !-------------------------------------------------------------------
      ppm_loadbal_nslp = ppm_loadbal_nslp + 1
      IF (ppm_loadbal_nslp .EQ. 1) THEN
          ppm_loadbal_slpavg = slp
      ELSE
          ppm_loadbal_slpavg = 0.5_ppm_kind_double*(ppm_loadbal_slpavg+slp)
      ENDIF
      !-------------------------------------------------------------------
      !  Extrapolate expected number of steps until redecomposition
      !-------------------------------------------------------------------
      estsol=2._ppm_kind_double*ppm_loadbal_decomp_cost/ppm_loadbal_slpavg
      IF (estsol .LT. 0.0_ppm_kind_double) THEN
          nredest = -1
          fail("Balance improved. Cannot extrapolate load imbalance.")

      ELSE
          nredest = INT(SQRT(estsol)) - nstep
          IF (nredest .LT. 0) nredest = 0
      ENDIF

      !-------------------------------------------------------------------
      !  Update old values
      !-------------------------------------------------------------------
      ppm_loadbal_old_sar = Wn
      ppm_loadbal_deltaold = delta

      end_subroutine()
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (.NOT. ppm_initialized) THEN
            fail("Please call ppm_init first!", exit_point=8888)
        ENDIF

        IF (ctime .LT. 0.0_MK) THEN
            fail("ctime must be >= 0.0",exit_point=8888)
        ENDIF
        IF (nstep .LT. 1) THEN
            fail("nstep must be > 0", exit_point=8888)
        ENDIF
        IF (nstep .LE. ppm_loadbal_nold) THEN
            fail("nstep did not increase. Use flush?", exit_point=8888)
        ENDIF

 8888   CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE loadbal_inquire_sar_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE loadbal_inquire_sar_d
#endif

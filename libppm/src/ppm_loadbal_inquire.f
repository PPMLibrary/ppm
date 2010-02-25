      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_loadbal_inquire
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Inquires about the load balance status and returns
      !                 advise whether redecomposing the problem is
      !                 recommended (based on some decision heuristic the
      !                 user chooses).
      !
      !  Input        : ctime      (F) elapsed time (as measured by
      !                                ppm_time) for all computation in one
      !                                time step on the local processor.
      !                 nstep      (I) number of time steps since last 
      !                                redecomposition. (>0)
      !                                If this routine is not called every
      !                                time step, linear interpolation of
      !                                the load imbalance will be used to
      !                                reconstruct missing data points.
      !                 heuristic  (I) Decision heuristic for
      !                                redecomposition advise. One of:
      !                                   ppm_param_loadbal_sar
      !                                       (Stop-at-Rise heuristic by
      !                                       Moon:1994)
      !                 lflush     (L) .TRUE. to flush internal statistics
      !                                (e.g. the first time this routine is
      !                                called after actually doing a
      !                                redecomposition of the problem),
      !                                .FALSE. to continue gathering
      !                                statistics.
      !
      !  Input/output : 
      !
      !  Output       : imbal      (F) load imbalance defined as the ratio
      !                                of computation time of the
      !                                bottleneck processor to the average
      !                                computation time of all processors.
      !                 lredecomp  (L) .TRUE. if the choosen heuristic
      !                                recommends problem redecomposition.
      !                                Otherwise .FALSE.. Redecomposition
      !                                means: do the ppm_topo_mktopo again.
      !                 nredest    (I) estimated (linear extrapolation)
      !                                number of time steps to go until next
      !                                advised redecomposition. -1 is
      !                                returned if the chosen heuristic
      !                                does not support this kind of
      !                                information. Be careful with this 
      !                                value!
      !                 info       (I) 0 on success.
      !
      !  Remarks      : The user should time (using ppm_time) computations
      !                 for the topology/topologies considered for dyanmic
      !                 remaping. The elapsed time is given to this
      !                 routine. This estimate is currently a scalar, so 
      !                 only one topology/set of topologies can be monitored.
      !                 Introducing a topology set ID (and corresponding
      !                 translation lists between internal and external
      !                 numbering of these sets) and making the internal
      !                 estimates a vector, this can later be extended to
      !                 multiple topology sets if needed. The topology set
      !                 ID (in external numbering) for which to evaluate the 
      !                 load balance will then be an additional argument
      !                 to this routine.
      !
      !                 This routine does two global MPI communication
      !                 operations (MPI_Allreduce).
      !
      !  References   : B. Moon and J. Saltz, Adaptive Runtime Support for
      !                 Direct Simulation Monte Carlo Methods on
      !                 Distributed Memory Architectures. Proceedings of
      !                 the IEEE Scalable High-Performance Computing 
      !                 Conference. 1994. 176--183. (bibtex: Moon:1994).
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_loadbal_inquire.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.4  2006/07/21 09:14:13  kotsalie
      !  FRIDAY
      !
      !  Revision 1.3  2004/07/26 07:45:26  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.2  2004/07/23 12:54:47  ivos
      !  Changed a warning to a notice if extrapolation is not possible.
      !
      !  Revision 1.1  2004/07/20 16:48:26  ivos
      !  Initial implementation. Not tested yet.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_loadbal_inquire_s(ctime,nstep,heuristic,lflush,imbal, &
     &    lredecomp,nredest,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_loadbal_inquire_d(ctime,nstep,heuristic,lflush,imbal, &
     &    lredecomp,nredest,info)
#endif
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
      INTEGER                , INTENT(IN   ) :: heuristic,nstep
      LOGICAL                , INTENT(IN   ) :: lflush
      REAL(MK)               , INTENT(  OUT) :: imbal
      LOGICAL                , INTENT(  OUT) :: lredecomp
      INTEGER                , INTENT(  OUT) :: nredest,info
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
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_loadbal_inquire',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (ctime .LT. 0.0_MK) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_loadbal_inquire', &
     &            'ctime must be >= 0.0',__LINE__,info)
             GOTO 9999
          ENDIF
          IF (nstep .LT. 1) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_loadbal_inquire', &
     &            'nstep must be > 0',__LINE__,info)
             GOTO 9999
          ENDIF
          IF (nstep .LE. nold) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_loadbal_inquire', &
     &            'nstep did not increase since last call! Use flush?', &
     &            __LINE__,info)
             GOTO 9999
          ENDIF
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

                  ppm_loadbal_runsum = ndiff*deltaold+REAL(ndiff,ppm_kind_double)*slp
     
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
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_loadbal_inquire_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_loadbal_inquire_d
#endif

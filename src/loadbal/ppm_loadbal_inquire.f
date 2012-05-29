      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_loadbal_inquire
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

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE loadbal_inq_s(ctime,nstep,lflush,lredecomp, &
     &                         nredest,info,heuristic,mov_avg_time,topoid)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE loadbal_inq_d(ctime,nstep,lflush,lredecomp, &
     &                         nredest,info,heuristic,mov_avg_time,topoid)
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
      !!! User must provide the moving average time ('mov_avg_time') and
      !!! topology ID ('topoid') UNLESS the chosen heuristic is SAR
      !!! (stop-at-rise) policy.
      !!! This means that if 'heuristic' is not present OR is equal to
      !!! 'ppm_param_loadbal_dlb', user MUST provide mov_avg_time.

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
      !!! Elapsed time (as measured by `ppm_time`) for all computation in one
      !!! time step on the local processor.
      INTEGER                , INTENT(IN   ) :: nstep
      !!! Number of time steps since last redecomposition. (>0)
      !!! If this routine is not called every time step, linear interpolation of
      !!! the load imbalance will be used to reconstruct missing data points.
      LOGICAL                , INTENT(IN   ) :: lflush
      !!! `TRUE` to flush internal statistics (e.g. the first time this routine
      !!! is called after actually doing a redecomposition of the problem),
      !!! `FALSE` to continue gathering statistics.
      INTEGER,OPTIONAL       , INTENT(IN   ) :: heuristic
      !!! Decision heuristic for redecomposition advise. One of:
      !!!
      !!! * ppm_param_loadbal_sar (Stop-at-Rise heuristic)
      !!! * ppm_param_loadbal_dlb (Dynamic load balancing heuristic)
      !!! If a heuristic is not specified, PPM will first try dynamic load
      !!! balancing. If DLB is bad, it will use SAR heuristic
      INTEGER,OPTIONAL       , INTENT(IN   ) :: topoid
      !!! Topology ID needed for DLB heuristic
      REAL(MK),OPTIONAL      , INTENT(IN   ) :: mov_avg_time
      !!! Moving (running) average time of the simulation.
      LOGICAL                , INTENT(  OUT) :: lredecomp
      !!! `TRUE` if the choosen heuristic recommends problem redecomposition.
      !!! Otherwise `FALSE`. Redecomposition means: do the ppm_topo_mktopo again.
      INTEGER                , INTENT(  OUT) :: nredest
      !!! Estimated (linear extrapolation) number of time steps to go until next
      !!! advised redecomposition. -1 is returned if the chosen heuristic
      !!! does not support this kind of information. Be careful with this value!
      !!! DLB will always return 1 (i.e. next step is the DLB step)
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      start_subroutine("ppm_loadbal_inquire")
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
      !  Does the user specify a heuristic?
      !-------------------------------------------------------------------------
      IF (PRESENT(heuristic)) THEN

        IF (heuristic .EQ. ppm_param_loadbal_sar) THEN
            !-------------------------------------------------------------------------
            !  Check that we have data available
            !-------------------------------------------------------------------------
            IF (ppm_loadbal_dcn .LT. 1) THEN
                fail("Use ppm_set_decomp_cost first to gather statistics")
            ENDIF

            !-------------------------------------------------------------------
            !  Stop-At-Rise heuristic by Moon:1994
            !
            !-------------------------------------------------------------------

            CALL ppm_loadbal_inquire_sar(ctime,nstep,lflush,lredecomp,nredest, &
     &                                   info)
                    or_fail("Something went wrong in ppm_loadbal_inquire_sar")

        ELSEIF (heuristic .EQ. ppm_param_loadbal_dlb) THEN
            !-------------------------------------------------------------------
            !  Using Omer's DLB heuristic
            !-------------------------------------------------------------------
            CALL ppm_loadbal_inquire_dlb(topoid,ctime,mov_avg_time,lflush,&
     &                                   lredecomp,info)
                    or_fail("Something went wrong in ppm_loadbal_inquire_dlb")
            !-------------------------------------------------------------------
            !  DLB heuristic expects load balancing to take place in the next
            !  time step
            !-------------------------------------------------------------------
            nredest=1
        ELSE
            !-------------------------------------------------------------------
            !  Unknow heuristic
            !-------------------------------------------------------------------
            fail("Unknown heuristic specified. Nothing computed.")
        ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Update time step counter
      !-------------------------------------------------------------------------
      ppm_loadbal_nold = nstep

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      end_subroutine()
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (.NOT. ppm_initialized) THEN
            fail("Please call ppm_init first!", exit_point=8888)
        ENDIF
        IF (ctime .LT. 0.0_MK) THEN
            fail("ctime must be >= 0.0", exit_point=8888)
        ENDIF
        IF (mov_avg_time .LT. 0.0_MK) THEN
            fail("moving avg time must be >= 0.0", exit_point=8888)
        ENDIF
        IF (nstep .LT. 1) THEN
            fail("nstep must be > 0",exit_point=8888)
        ENDIF

        IF (nstep .LE. ppm_loadbal_nold) THEN
            fail("nstep did not increase since last step. Use flush?", &
     &                  exit_point=8888)
        ENDIF
 8888   CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE loadbal_inq_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE loadbal_inq_d
#endif

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
      REAL(MK)                               :: imbal_perc,max_ctime,r_neigh
      REAL(MK)                               :: neigh_imbal,recv_ctime
      INTEGER                                :: random_neigh,i,tag1,sendrank
      INTEGER                                :: recvrank
      LOGICAL                                :: l_myneighbor
      TYPE(ppm_t_topo),POINTER               :: topo => NULL()
#ifdef __MPI
      INTEGER                                :: MPTYPE
      INTEGER, DIMENSION(MPI_STATUS_SIZE)    :: status
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      start_subroutine("ppm_loadbal_inquire")
      lredecomp = .FALSE.
      nredest = -1 
      l_myneighbor = .TRUE.
      neigh_imbal = 0._MK
      recv_ctime = 0._MK

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
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
      !  Get the topology first
      !-------------------------------------------------------------------------
      topo    => ppm_topo(topoid)%t

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
            CALL ppm_loadbal_inquire_dlb(topoid,ctime,nstep,mov_avg_time,&
     &                                   max_ctime,imbal_perc,info)
                    or_fail("Something went wrong in ppm_loadbal_inquire_dlb")

            !-------------------------------------------------------------------
            !  If load imbalance is more than 25% (I'm guessing this number),
            !  we need a DLB
            !-------------------------------------------------------------------
            IF (imbal_perc .GT. 0.25_MK) THEN
                lredecomp = .TRUE.
            ELSE
                lredecomp = .FALSE.
            ENDIF
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

      ELSE
            !-------------------------------------------------------------------
            !  No heuristic is provided. PPM will check first if we need LB
            !  by calling Omer's DLB heuristic
            !-------------------------------------------------------------------
            CALL ppm_loadbal_inquire_dlb(topoid,ctime,nstep,mov_avg_time,&
     &                                   max_ctime,imbal_perc,info)
                    or_fail("Something went wrong in ppm_loadbal_inquire_dlb")
            !-------------------------------------------------------------------
            !  We know now the imbalance percentage. We need to check how
            !  the (almost) global load imbalance looks like.
            !  For this every processor chooses a random neighboring
            !  processor and communicates the max_ctimes. Since my neighbor has a
            !  different neighborhood topology, I will (probably) obtain some
            !  an other  If the difference
            !  is big, this will hint a global imbalance. One can later call SAR
            !  or directly do a global remap.
            !  Find a processor first which is not me and not in my neighborhood
            !-------------------------------------------------------------------
            DO WHILE(l_myneighbor)
                CALL RANDOM_NUMBER(r_neigh) ! between 0.0-1.0
                random_neigh = INT(r_neigh*ppm_nproc) ! between 0-ppm_nproc
                IF (random_neigh .NE. ppm_rank) l_myneighbor = .FALSE.
                DO i=1,topo%nneighproc
                    IF (random_neigh .EQ. topo%ineighproc(i)) THEN
                        l_myneighbor = .TRUE.
                        EXIT
                    ENDIF
                ENDDO
            ENDDO

            !-------------------------------------------------------------------
            !  Send & receive the max_ctimes with a neighbor
            !-------------------------------------------------------------------
            tag1 = 678
            sendrank = random_neigh
            recvrank = sendrank
            !-------------------------------------------------------------------
            !  only consider non-negative sendranks
            !-------------------------------------------------------------------
            IF (sendrank.GE.0 .AND. sendrank.NE.ppm_rank) THEN
                CALL MPI_SendRecv(max_ctime ,1,MPTYPE,sendrank,tag1, &
     &                            recv_ctime,1,MPTYPE,recvrank,tag1, &
     &                            ppm_comm,status,info)
                !---------------------------------------------------------------
                !  If a neighborhood's slowest processor is 50% slower than
                !  other neighborhood's slowest processor, do a global mapping
                !---------------------------------------------------------------
                IF (recv_ctime .GT. max_ctime) THEN
                    neigh_imbal = (recv_ctime-max_ctime)/recv_ctime
                ELSE
                    neigh_imbal = (max_ctime-recv_ctime)/max_ctime
                ENDIF

                IF (neigh_imbal .GT. 0.5_MK) THEN
                    ! DO A GLOBAL MAPPING
                ELSE
                    ! DO DLB
                ENDIF
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

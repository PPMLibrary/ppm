      minclude ppm_header(ppm_module_name="loadbal_bc")

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE loadbal_bc_s(topoid,meshid,Pc,my_time,info,&
     &                          vlist,nvlist,pcost,t_comp)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE loadbal_bc_d(topoid,meshid,Pc,my_time,info,&
     &                          vlist,nvlist,pcost,t_comp)
#endif
      !!! This routine executes the "Balancing Circuit" model of dynamic
      !!! load balancing. In this model every processor follows its list
      !!! of communication sequence and equalizes its load with its neighbor.
      !!! Sauerwald:2012 describes this model and they derive the number of
      !!! steps needed for the algorithm to reach a K-close solution.
      !!!
      !!! [NOTE]
      !!! This routine needs to be called after ParticleSet%map_ghost
      !!!
      !!! .References
      !!! *************************************************************
      !!! "Tight Bounds for Randomized Load Balancing on Arbitrary
      !!!  Network Topologies" by Sauerwald et al.
      !!! 2012
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
      USE ppm_module_util_commopt
      USE ppm_module_interfaces
      USE ppm_module_topo_typedef
      USE ppm_module_mapping_typedef
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
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! The topology where the DLB will take place
      INTEGER                 , INTENT(IN   ) :: meshid
      !!! mesh ID for which to compute the cost.
      !!! If -1 is passed, only particles are considered and the mesh is
      !!! ignored.                                                             +
      !!! If there are no particles and mesh_id is -1, the costs are
      !!! computed based on the geometry: cost(sub) = volume(sub).
#if    __KIND == __SINGLE_PRECISION
      TYPE(ppm_t_particles_s) , INTENT(IN   ) :: Pc
#else
      TYPE(ppm_t_particles_d) , INTENT(IN   ) :: Pc
#endif
      !!! Particle data structure (particle set)
      REAL(MK)               , INTENT(IN   )  :: my_time
      !!! The total computation time for this time step
      INTEGER                , INTENT(  OUT)  :: info
      !!! Returns status, 0 upon success
      INTEGER, DIMENSION(:)  , OPTIONAL, INTENT(IN   ) :: nvlist
      !!! Number of particles with which ip has to interact. Index: ip.
      INTEGER, DIMENSION(:,:), OPTIONAL, INTENT(IN   ) :: vlist
      !!! Verlet list. First index: particles with which particle ip interacts.
      !!! Second index: ip. The second index only runs up to the
      !!! largest ip with non-zero nvlist. This is to save memory since the
      !!! last particles are the ghosts and they do not have a verlet list.
      !!! This is only allocated and returned if lstore is .TRUE.
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN) :: pcost
      !!! per-particle costs.
      !!! If it is not present but the nvlist and vlist are given,
      !!! pcost is computed by PPM using the relative communication cost
      !!! If it is not present AND nvlist and vlist are NOT given,
      !!! pcost is 1.0 for each particle
      REAL(MK), OPTIONAL     , INTENT(IN   ) :: t_comp
      !!! Elapsed time (as measured by `ppm_time`) for all computation in one
      !!! time step on the local processor.
!      REAL(MK), OPTIONAL     , INTENT(IN   ) :: t_comm
!      !!! Elapsed time for all COMMUNICATION in one time step on the local proc
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                           :: neighbor_time,cost_ratio,time_ratio
      REAL(MK)                           :: my_cost,n_time,n_cost,balancing_load
      REAL(MK)                           :: ideal_load
      REAL(MK),DIMENSION(1)              :: sendbuf,recvbuf
      ! Neighbor's elapsed time in this iteration
      REAL(MK),DIMENSION(:),POINTER      :: weighted_time=>NULL(),helper=>NULL()
      INTEGER,DIMENSION(:),POINTER       :: n_ineighproc=>NULL()
      INTEGER,DIMENSION(:),POINTER       :: candidate_sublist=>NULL()
      INTEGER,DIMENSION(:),POINTER       :: send_sublist=>NULL()
      INTEGER,DIMENSION(:),POINTER       :: recv_sublist=>NULL()
      INTEGER                            :: balancing_subloc
      INTEGER,DIMENSION(2)               :: temp
      INTEGER                            :: ndiff,i,k,kk,iset,ibuffer,tag1,jj
      INTEGER                            :: nneighs,groupsize,nsteps
      INTEGER                            :: counter,isub,jproc,n_nneighproc,jsub
      INTEGER                            :: balancing_subID,recv_subsize
      TYPE(ppm_t_topo),POINTER           :: topo => NULL()
      LOGICAL                            :: isSender
      LOGICAL,DIMENSION(:),POINTER       :: isNewNeighbor
#ifdef __MPI
      INTEGER                            :: MPTYPE
      INTEGER, DIMENSION(MPI_STATUS_SIZE):: status
#endif
      REAL(MK), DIMENSION(:,:), POINTER  :: xp => NULL()
      ! Particle positions
      INTEGER                            :: Np
      ! Number of particles. Set to <= 0 if mesh-based costs are desired.
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      start_subroutine("loadbal_bc")
!      CALL Pc%get_xp(xp,info)
!       or_fail("Pc%get_xp failed")
      xp => Pc%xp
      Np = Pc%Npart
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
      !  Get the number of my neighbors
      !-------------------------------------------------------------------------
      nneighs =  topo%nneighproc

      !-------------------------------------------------------------------------
      !  Initialize neighbor's time
      !-------------------------------------------------------------------------
      neighbor_time = 0._MK

      !-------------------------------------------------------------------------
      ! First check if the optimal communication protocol is known
      !-------------------------------------------------------------------------
      IF (.NOT. topo%isoptimized) THEN
        ! if not: determine it before running the communication operations
        print*,'COMMUNICATION NOT OPTIMIZED!'
        CALL ppm_util_commopt(topoid,info)
        IF (info.NE.0) GOTO 9999
        IF (ppm_debug .GT. 1) THEN
            DO i=1,nneighs
                stdout_f('(A,I4)',"have neighbor: ",  'topo%ineighproc(i)')
            END DO
            DO i=1,topo%ncommseq
                stdout_f('(A,I4)',"communicate: ", 'topo%icommseq(i)')
            END DO
        ENDIF
      END IF

      !-------------------------------------------------------------------------
      !  Find out the average compute time of a particle
      !  If user provides only the overall time (comp time + comm time),
      !  subtract the ppm_loadbal_comm_time from that.
      !-------------------------------------------------------------------------
      IF (PRESENT(t_comp)) THEN
        ppm_loadbal_comp_part_time = t_comp/Np
      ELSE
        ppm_loadbal_comp_part_time = (my_time-ppm_loadbal_comm_time)/Np
      ENDIF
      stdout("comp_part_time",ppm_loadbal_comp_part_time)
      ppm_loadbal_comp_ratio       = 1._MK / ppm_loadbal_comp_part_time
      !-------------------------------------------------------------------------
      !  Get the estimated cost model
      !-------------------------------------------------------------------------
      CALL ppm_get_cost(topoid,meshid,xp,Np,info,vlist,nvlist,pcost)
       or_fail("Something went wrong in the ppm_get_cost()")
      !-------------------------------------------------------------------------
      !  The optimized communication sequence gives us the pre-described
      !  matching sequence needed by the balancing circuit algorithm.
      !  The number of matching matrices is equal to the number of colors in the
      !  edge coloring algorithm of the communication graph.
      !  This process has to be repeated until we are K-close to an ideal soln.
      !-------------------------------------------------------------------------
!      print*,ppm_rank,': ppm_loadbal_subcostd=',ppm_loadbal_subcostd(1:topo%nsublist)
      nsteps = 1
      DO k=1,1
        groupsize = 0
          DO i=2,ppm_nsendlist
            !-------------------------------------------------------------------
            !  Minimize the load difference as much as possible
            !  TODO: REPLACE THIS W/ A WHILE-LOOP BY USING A THRESHOLD
            !-------------------------------------------------------------------
            DO kk =1,nsteps
                !---------------------------------------------------------------
                !  only consider non-negative send/recv ranks
                !---------------------------------------------------------------
                IF (ppm_isendlist(i).GE.0 .AND. ppm_irecvlist(i).GE.0) THEN
                    tag1 = 123
                    !-----------------------------------------------------------
                    !  Send & receive estimated computational cost
                    !
                    !-----------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
                    sendbuf(1) = ppm_loadbal_proccosts
#else
                    sendbuf(1) = ppm_loadbal_proccostd
#endif
                    CALL MPI_SendRecv(sendbuf,1,MPTYPE,ppm_isendlist(i),tag1, &
     &                                recvbuf,1,MPTYPE,ppm_irecvlist(i),tag1, &
     &                                ppm_comm,status,info)
                    or_fail("Sendrecv problem for time and cost exchange!")

                    stdout("my cost", 'sendbuf(1)',", ",'ppm_irecvlist(i)',&
     &                     " cost:",'recvbuf(1)')
                    !-----------------------------------------------------------
                    !  Now, relate the measured time to the cost model
                    !  Check if the ratios are close to each other.
                    !  If not, adjust the comm_cost/comp_cost ratio
                    !  TODO: Adjust comm_cost/comp_cost ratio
                    !-----------------------------------------------------------
                    n_cost = recvbuf(1)
                    my_cost= sendbuf(1)
                    !-----------------------------------------------------------
                    !  Should I send or receive subdomains?
                    !-----------------------------------------------------------
                    isSender = .FALSE.
                    IF (my_cost .GT. n_cost)  isSender = .TRUE.

                    !-----------------------------------------------------------
                    !  The exact balance in the cost is the arithmetic mean
                    !  and a 'balancing_load' is sent from heavier loaded
                    !  proc to an underloaded one.
                    !-----------------------------------------------------------
                    ideal_load = ABS(my_cost - n_cost)/2._MK
                    !-----------------------------------------------------------
                    !  ...BUT usually I can't find a subdomain to sent whose
                    !  cost would be identical to 'ideal_load'.
                    !  Finding out how many subs to send is an interesting
                    !  problem.
                    !  Criteria on choosing a good subdomain:
                    !  1. No new neighborhoods! Comm graph must remain the SAME
                    !  2. Choose a sub whose cost is closest to 'balancing_load'
                    !  3. You can also 'STAY' (like in a BlackJack game)
                    !     meaning that, if ALL subs you can choose would create
                    !     a greater imbalance, then better do not send ANYTHING!
                    !  4. OR, you may also send multiple subs.
                    !-----------------------------------------------------------
                    IF (isSender) THEN
                        !-------------------------------------------------------
                        !  Check if I have more than one sub, otherwise I can't
                        !  balance the loads
                        !-------------------------------------------------------
                        IF (topo%nsublist.GT.1) THEN
                            !---------------------------------------------------
                            !  Get neighboring processor's list of its neighbors
                            !  First get neighbor's nneighproc, then the actual
                            !  list itself.
                            !---------------------------------------------------
                            tag1 = 124
                            CALL MPI_Recv(n_nneighproc,1,MPI_INTEGER,            &
     &                                    ppm_isendlist(i),tag1,ppm_comm,status, &
     &                                    info)
                                or_fail("n_nneighproc MPI_Recv failed")
                            ALLOCATE(n_ineighproc(1:n_nneighproc),STAT=info)
                                or_fail_alloc("n_ineighproc")
                            tag1 = 125
                            CALL MPI_Recv(n_ineighproc,n_nneighproc,MPI_INTEGER, &
     &                                    ppm_isendlist(i),tag1,ppm_comm,status, &
     &                                    info)
                                or_fail("n_ineighproc MPI_Recv failed")


                            !---------------------------------------------------
                            !  Create a list of possible subdomain candidates
                            !  that are cleared to be sent to the underloaded
                            !  proc.
                            !---------------------------------------------------
                            ALLOCATE(candidate_sublist(1:topo%nsublist),STAT=info)
                                or_fail_alloc("candidate_sublist")
                            candidate_sublist = 0
                            !---------------------------------------------------
                            !  Iterate over my subs to find out which subs are
                            !  candidates ..
                            !---------------------------------------------------
                            counter = 0
                            DO isub=1,topo%nsublist
                                !-----------------------------------------------
                                !  Do I introduce a new neighboring proc to the
                                !  underloaded proc?
                                !  Initially, I assume YES!
                                !-----------------------------------------------
                                ALLOCATE(isNewNeighbor(1:topo%nneighsubs(isub)),&
     &                                   STAT=info)
                                    or_fail_alloc("isNewNeighbor")
                                isNewNeighbor = .TRUE.
                                !-----------------------------------------------
                                !  ... find their neighboring sub
                                !-----------------------------------------------
                                DO jsub=1,topo%nneighsubs(isub)
                                    !-------------------------------------------
                                    !  ... for each neighboring sub check its
                                    !  sub2proc. If even one sub2proc is diff
                                    !  then we can't send it to the underloaded
                                    !  proc as this would violate the initial
                                    !  comm graph. We can't promise to reach
                                    !  a global load balance in that case.
                                    !-------------------------------------------
                                    DO jproc=1,n_nneighproc
                                        IF (topo%ineighsubs(jsub,isub).EQ.       &
     &                                      n_ineighproc(jproc)) THEN
                                            !-----------------------------------
                                            !  This isub can be sent since its
                                            !  neighbors' procs are also
                                            !  underloaded proc's neighbors
                                            !-----------------------------------
                                            isNewNeighbor(jsub) = .FALSE.
                                        ENDIF
                                    ENDDO
                                ENDDO !iterate over neighboring subs

                                !-----------------------------------------------
                                !  Add this sub to the list of possible
                                !  candidates to be sent to the underloaded proc
                                !  candidate_sublist keeps local sub IDs
                                !-----------------------------------------------
                                IF (.NOT. ANY(isNewNeighbor)) THEN
                                    counter = counter + 1
                                    candidate_sublist(counter) = isub
                                ENDIF
                            ENDDO !iterate over all subs
                            !---------------------------------------------------
                            !  Now, I gotta choose which sub to send.
                            !  But I will filter them using a mask.
                            !  The sub with max cost that is lower than
                            !  ideal_load will be sent
                            !---------------------------------------------------
                            stdout("CHOOSING a SUB")
                            CALL ppm_loadbal_choose_sub(ideal_load,candidate_sublist,&
     &                           send_sublist,counter,info)
                                or_fail("ppm_loadbal_choose_sub")
                            !---------------------------------------------------
                            !  Communicate the list of to-be-sent subs with my
                            !  underloaded partner.
                            !  First : Size.......
                            !---------------------------------------------------

                            tag1 = 126
                            print*,'SIZE(send_sublist,1):',SIZE(send_sublist,1)
                            CALL MPI_Send(SIZE(send_sublist,1),1,                 &
     &                              MPI_INTEGER,ppm_isendlist(i),tag1,ppm_comm,   &
     &                              status,info)
                            or_fail("MPI_Send::send_sublist")

                            !---------------------------------------------------
                            !  Second: Actual list
                            !---------------------------------------------------
                            tag1 = 127
                            CALL MPI_Send(send_sublist,SIZE(send_sublist,1),  &
     &                              MPI_INTEGER,ppm_isendlist(i),tag1,ppm_comm,   &
     &                              status,info)
                            or_fail("MPI_Send::send_sublist")
                            !---------------------------------------------------
                            !  If the list has only one subID and that is an
                            !  invalid one (i.e. -1), it means I have
                            !  nothing to send to the underloaded processor
                            !  Thus, I can skip this neighbor in this round
                            !---------------------------------------------------
                            IF (SIZE(send_sublist,1).EQ.1 .AND.                 &
     &                          send_sublist(1).EQ.-1) THEN
                                stdout("Nothing to send....")
                                GOTO 7777
                            !---------------------------------------------------
                            ! ... ELSE I send subdomains one after another
                            !---------------------------------------------------
                            ELSE
                                DO isub=1,SIZE(send_sublist)
                                    IF (send_sublist(isub).NE.-1) THEN
                                        CALL ppm_loadbal_sendsub(topoid,Pc,       &
     &                                        ppm_isendlist(i),send_sublist(isub),&
     &                                        SIZE(topo%ineighsubs,1),topo%prec,  &
     &                                        info)
                                    ENDIF
                                ENDDO
                            ENDIF
                        ELSE
                            fail("I do NOT have enough subs!!")
                        ENDIF

                    ELSE
                        !-------------------------------------------------------
                        !  I'm the underloaded proc.
                        !  I need to send nneighproc and my neighboring processor
                        !  list so that the sendercan select correct sub(s) to
                        !  send me.
                        !-------------------------------------------------------
                        tag1 = 124
                        CALL MPI_Send(n_nneighproc,1,MPI_INTEGER,            &
     &                                ppm_irecvlist(i),tag1,ppm_comm,status, &
     &                                info)
                            or_fail("n_nneighproc MPI_Send failed")
                        ALLOCATE(n_ineighproc(1:n_nneighproc),STAT=info)
                            or_fail_alloc("n_ineighproc")
                        n_ineighproc = 0
                        tag1 = 125
                        CALL MPI_Send(n_ineighproc,n_nneighproc,MPI_INTEGER, &
     &                                ppm_irecvlist(i),tag1,ppm_comm,status, &
     &                                info)
                            or_fail("n_ineighproc MPI_Send failed")
                        !-------------------------------------------------------
                        !  This is the part where the SENDER should be deciding
                        !  which sub to choose. (like waiting for girlfriend to
                        !  choose what to wear)...
                        !  .............
                        !  OK, now I can get the list of the subs I should recv
                        !  First : Size of the list...
                        !------------------------------------------------------
                        tag1 = 126
                        CALL MPI_Recv(recv_subsize,1,MPI_INTEGER,ppm_irecvlist(i),&
     &                                tag1,ppm_comm,status,info)
                        or_fail("MPI_Recv::recv_subsize")
                        !-------------------------------------------------------
                        !  ... allocate the recv_sublist...
                        !-------------------------------------------------------
                        ALLOCATE(recv_sublist(1:recv_subsize),STAT=info)
                            or_fail_alloc("recv_sublist")
                        !-------------------------------------------------------
                        !  ... now the actual list
                        !-------------------------------------------------------
                        tag1 = 127
                        CALL MPI_Recv(recv_sublist,recv_subsize,MPI_INTEGER,  &
     &                                ppm_irecvlist(i),tag1,ppm_comm,status,info)
                        or_fail("MPI_Recv::recv_sublist")
                        !-------------------------------------------------------
                        !  If the list has only one subID and that is an
                        !  invalid one (i.e. -1), it means I have
                        !  nothing to receive from the overloaded processor
                        !  Thus, I can skip this neighbor in this round
                        !-------------------------------------------------------
                        IF (recv_subsize.EQ.1 .AND. recv_sublist(1).EQ.-1) THEN
                            stdout("Nothing to receive....")
                            GOTO 7777
                        !-------------------------------------------------------
                        ! ... ELSE I receive subdomains one after another
                        !-------------------------------------------------------
                        ELSE
                            DO isub=1,recv_subsize
                                IF (recv_sublist(isub).NE.-1) THEN
                                    CALL ppm_loadbal_recvsub(topoid,            &
     &                                      ppm_irecvlist(i),topo%prec,info)
                                    or_fail("ppm_loadbal_recvsub")
                                ENDIF
                            ENDDO
                        ENDIF
                    ENDIF ! isSender?
                ENDIF ! I'm either a sender/receiver in this comm round
            ENDDO !kk =1,nsteps
 7777       CONTINUE
          ENDDO
      ENDDO !threshold to reach k-close solution
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

        IF (my_time .LT. 0.0_MK) THEN
            fail("my_time must be >= 0.0",exit_point=8888)
        ENDIF

 8888   CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE loadbal_bc_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE loadbal_bc_d
#endif

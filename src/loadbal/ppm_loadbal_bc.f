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
      TYPE(ppm_t_particles_s) , INTENT(INOUT) :: Pc
#else
      TYPE(ppm_t_particles_d) , INTENT(INOUT) :: Pc
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
      REAL(MK)                           :: ideal_load,ratio
      REAL(MK),DIMENSION(1)              :: sendbuf,recvbuf
      INTEGER                            :: sendbuf_i,recvbuf_i
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
      INTEGER                            :: send_subsize
      TYPE(ppm_t_topo),POINTER           :: topo => NULL()
      LOGICAL                            :: isSender
      LOGICAL,DIMENSION(:,:),POINTER     :: isNewNeighbor,isEngulfed
#ifdef __MPI
      INTEGER                            :: MPTYPE
      INTEGER, DIMENSION(MPI_STATUS_SIZE):: status
#endif
      REAL(MK), DIMENSION(:,:), POINTER  :: xp => NULL()
      ! Particle positions
      INTEGER                            :: Np,Mp
      ! Number of particles. Set to <= 0 if mesh-based costs are desired.
      INTEGER, DIMENSION(2)              :: ldu
      INTEGER                            :: iopt
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      start_subroutine("loadbal_bc")
      CALL Pc%get_xp(xp,info)
       or_fail("Pc%get_xp failed")
!      xp => Pc%xp
      Np = Pc%Npart
      stdout("========DLB STARTS========")
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
!      stdout("TOTAL TIME:",my_time,"COMM TIME:",ppm_loadbal_comm_time)
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
!      stdout("comp_part_time",ppm_loadbal_comp_part_time)
      ppm_loadbal_comp_ratio       = 1._MK / ppm_loadbal_comp_part_time

      !-------------------------------------------------------------------------
      !  The optimized communication sequence gives us the pre-described
      !  matching sequence needed by the balancing circuit algorithm.
      !  The number of matching matrices is equal to the number of colors in the
      !  edge coloring algorithm of the communication graph.
      !  This process has to be repeated until we are K-close to an ideal soln.
      !-------------------------------------------------------------------------
!      print*,ppm_rank,': ppm_loadbal_subcostd=',ppm_loadbal_subcostd(1:topo%nsublist)
      nsteps = 1
!      stdout("ppm_isendlist",ppm_isendlist,"ppm_nsendlist",ppm_nsendlist)
!      stdout("size(ppm_isendlist)",'size(ppm_isendlist)')
!      stdout("ppm_irecvlist",ppm_irecvlist)
      !-------------------------------------------------------------------------
      !  In the ppm_loadbal_map_subpart routine I'm modifying ppm_isendlist to
      !  be able to use already existing code to send subs during DLB. This
      !  however screws up original ppm_isendlist. So, I'll copy the original
      !  list to recall it later.
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(ppm_loadbal_isendlist,ldu,iopt,info)
       or_fail_alloc("ppm_isendlist")
      CALL ppm_alloc(ppm_loadbal_irecvlist,ldu,iopt,info)
       or_fail_alloc("ppm_irecvlist")
      ppm_loadbal_isendlist = ppm_isendlist
      ppm_loadbal_irecvlist = ppm_irecvlist
      ppm_loadbal_nsendlist = ppm_nsendlist
      DO k=1,1
        groupsize = 0
          DO i=2,ppm_loadbal_nsendlist
            !-------------------------------------------------------------------
            !  Minimize the load difference as much as possible
            !  TODO: REPLACE THIS W/ A WHILE-LOOP BY USING A THRESHOLD
            !
            !-------------------------------------------------------------------

            DO kk =1,nsteps
                !---------------------------------------------------------------
                !  only consider non-negative send/recv ranks
                !---------------------------------------------------------------
!           stdout("modified isend now",'ppm_isendlist(i)',"irecv now",'ppm_irecvlist(i)')
!           stdout("original isend now",'ppm_loadbal_isendlist(i)',"irecv now",'ppm_loadbal_irecvlist(i)')

                IF (ppm_loadbal_isendlist(i).GE.0 .AND. ppm_loadbal_irecvlist(i).GE.0) THEN
                    tag1 = 122
                    !-----------------------------------------------------------
                    !  Send & receive estimated computational cost
                    !
                    !-----------------------------------------------------------
                    sendbuf(1) = my_time
                    CALL MPI_SendRecv(sendbuf,1,MPTYPE,ppm_loadbal_isendlist(i),tag1, &
     &                                recvbuf,1,MPTYPE,ppm_loadbal_irecvlist(i),tag1, &
     &                                ppm_comm,status,info)
                    or_fail("Sendrecv problem for time!")
                    neighbor_time = recvbuf(1)
                    !-----------------------------------------------------------
                    !  Get the estimated cost model
                    !-----------------------------------------------------------
                    ratio = my_time/neighbor_time
                    CALL ppm_get_cost(topoid,meshid,xp,Np,ratio,info,vlist,nvlist,pcost)
                    or_fail("Something went wrong in the ppm_get_cost()")
#if    __KIND == __SINGLE_PRECISION
                    sendbuf(1) = ppm_loadbal_proccosts
#else
                    sendbuf(1) = ppm_loadbal_proccostd
#endif

                    CALL MPI_SendRecv(sendbuf,1,MPTYPE,ppm_loadbal_isendlist(i),tag1, &
     &                                recvbuf,1,MPTYPE,ppm_loadbal_irecvlist(i),tag1, &
     &                                ppm_comm,status,info)
                    or_fail("Sendrecv problem for cost exchange!")


                    stdout("my cost", 'sendbuf(1)',",",'ppm_loadbal_irecvlist(i)',&
     &                     " cost:",'recvbuf(1)')
                    !-----------------------------------------------------------
                    !  Now, relate the measured time to the cost model
                    !  Check if the ratios are close to each other.
                    !  If not, adjust the comm_cost/comp_cost ratio
                    !-----------------------------------------------------------
                    n_cost = recvbuf(1)
                    my_cost= sendbuf(1)
                    !-----------------------------------------------------------
                    !  The exact balance in the cost is the arithmetic mean
                    !  and a 'balancing_load' is sent from heavier loaded
                    !  proc to an underloaded one.
                    !-----------------------------------------------------------
                    ideal_load = ABS(my_cost - n_cost)/2._MK
!                    stdout("ideal load:",ideal_load)
                    !-----------------------------------------------------------
                    !  Should I send or receive subdomains?
                    !-----------------------------------------------------------
                    isSender = .FALSE.
                    IF (my_cost .GT. n_cost) THEN
                        isSender = .TRUE.
                    ELSE IF (my_cost .EQ. n_cost) THEN
                        stdout("Perfectly balanced.. moving on to next round")
                        stdout("Quitting DLB")
                        GOTO 7777
                    ELSE
                        isSender = .FALSE.
                    ENDIF
                    !-----------------------------------------------------------
                    !  If SENDER has only 1 sub, we can't do load balancing
                    !  This info should be told to the RECEIVER
                    !-----------------------------------------------------------
                    stdout("# of subs on this proc:",'topo%nsublist')
                    sendbuf_i  = topo%nsublist
                    tag1 = 123
                    CALL MPI_SendRecv(sendbuf_i,1,MPI_INTEGER,ppm_loadbal_isendlist(i),tag1, &
     &                                recvbuf_i,1,MPI_INTEGER,ppm_loadbal_irecvlist(i),tag1, &
     &                                ppm_comm,status,info)
                    !-----------------------------------------------------------
                    !  Quit DLB if SENDER has only one sub.
                    !-----------------------------------------------------------
                    IF ( (isSender       .AND. sendbuf_i.EQ.1) .OR.               &
     &                   (.NOT. isSender .AND. recvbuf_i.EQ.1)) THEN
                        stdout("Only one sub! Load balancing is not possible!")
                        stdout("Quitting DLB")
                        GOTO 9999
                    ENDIF
                    !-----------------------------------------------------------
                    !  OK, if we made this far, we know that the SENDER has at
                    !  least two subs. So, a DLB is theoretically possible.
                    !  Yet, we need to do more in-depth analysis.
                    !  Finding out how many subs (and which ones) to send is an
                    !  interesting problem.
                    !  Criteria on choosing a good subdomain:
                    !  1. No new neighborhoods! Comm graph must remain the SAME
                    !  2. Choose a sub whose cost is closest to 'balancing_load'
                    !  3. Choose a sub that would not be totally surrounded by
                    !     SENDER. Otherwise, we would create patches.
                    !  4. You can also 'STAY' (like in a BlackJack game)
                    !     meaning that, if ALL subs you can choose would create
                    !     a greater imbalance, then better do not send ANYTHING!
                    !  5. OR, you may also send multiple subs.
                    !-----------------------------------------------------------
                    IF (isSender) THEN
                        !-------------------------------------------------------
                        !  Get RECEIVING processor's list of its neighbors
                        !  First get neighbor's nneighproc, then the actual
                        !  list itself.
                        !-------------------------------------------------------
!                        stdout("MPI_SENDRECV----OK")
                        tag1 = 124
                        CALL MPI_Recv(n_nneighproc,1,MPI_INTEGER,ppm_loadbal_isendlist(i),&
     &                                tag1,ppm_comm,status,info)
                        or_fail("n_nneighproc MPI_Recv failed")
!                        stdout("MPI_RECV----OK")
                        !-------------------------------------------------------
                        !  Once you have the size of the RECEIVING proc's number
                        !  of neighboring processors, allocate the list to be
                        !  ready to receive the actual list
                        !-------------------------------------------------------
                        ldu(1) = n_nneighproc
                        iopt   = ppm_param_alloc_fit
                        CALL ppm_alloc(n_ineighproc,ldu,iopt,info)
                        or_fail_alloc("n_ineighproc")
                        !-------------------------------------------------------
                        !  Get the list of neighboring procs
                        !-------------------------------------------------------
                        tag1 = 125
                        CALL MPI_Recv(n_ineighproc,n_nneighproc,MPI_INTEGER, &
     &                                ppm_loadbal_isendlist(i),tag1,ppm_comm,status,info)
                        or_fail("n_ineighproc MPI_Recv failed")
!                        stdout("my neighbor s neighprocs",n_ineighproc)
!                        stdout("MPI_RECV2----OK")
                        !-------------------------------------------------------
                        !  Create a list of possible subdomain candidates
                        !  that are cleared to be sent to the underloaded proc.
                        !-------------------------------------------------------
                        ldu(1) = topo%nsublist
                        CALL ppm_alloc(candidate_sublist,ldu,iopt,info)
                        or_fail_alloc("candidate_sublist")
                        candidate_sublist = 0
                        !------------------------------------------------------
                        !  I assume initially that every sub of the SENDER
                        !  violates criterion (1). To clear a sub I introduce
                        !  isNewNeighbor logical array.
                        !------------------------------------------------------
                        ldu(1) = topo%nsublist
                        ldu(2) = MAXVAL(topo%nneighsubs)
                        CALL ppm_alloc(isNewNeighbor,ldu,iopt,info)
                        or_fail_alloc("isNewNeighbor")
                        isNewNeighbor = .TRUE.
                        !------------------------------------------------------
                        !  I assume initially that every sub is engulfed totally
                        !  by other SENDER's subs.
                        !------------------------------------------------------
                        ldu(1) = topo%nsublist
                        ldu(2) = MAXVAL(topo%nneighsubs)
                        CALL ppm_alloc(isEngulfed,ldu,iopt,info)
                        or_fail_alloc("isEngilfed")
                        isEngulfed = .TRUE.

                        !-------------------------------------------------------
                        !  Now, iterate over my subs to find out which subs are
                        !  candidates...
                        !  Count the number of possible candidates
                        !-------------------------------------------------------
                        counter = 0
!                        print*,ppm_rank,'topo%sub2proc',topo%sub2proc
!                        print*,ppm_rank,'n_nneighproc:',n_nneighproc
                        DO isub=1,topo%nsublist
                            !---------------------------------------------------
                            !  ... find their neighboring sub
                            !  isub and jsub are both local sub IDs
                            !---------------------------------------------------
!                            stdout("topo%nneighsubs(isub):",'topo%nneighsubs(isub)')
!                            stdout("sizeof ineighsubs",'size(topo%ineighsubs(:,isub),1)'
                            DO jsub=1,topo%nneighsubs(isub)
                                !-----------------------------------------------
                                !  ..for each neighboring sub check its sub2proc
                                !  If even one sub2proc is diff then we can't
                                !  send it to the underloaded proc as this would
                                !  violate the initial comm graph. We can't
                                !  promise to reach a global load balance in that
                                !  case.
                                !  In addition to this, the selected sub should
                                !  not be engulfed by other processor completely.
                                !  So, at least one jsub should belong to Sender
                                !  Iterate over overloaded proc's neighboring
                                !  procs and compare
                                !-----------------------------------------------
                                DO jproc=1,n_nneighproc
!                                    stdout("isub",isub)
!                                    stdout("globalID of isub",'topo%isublist(isub)')
!                                    stdout("globalID of jsub",'topo%ineighsubs(jsub,isub)')

                                    IF (topo%sub2proc(topo%ineighsubs(jsub,isub)).EQ.&
     &                                  n_ineighproc(jproc)       .OR.       &
     &                                  topo%sub2proc(topo%ineighsubs(jsub,isub))   .EQ.       &
     &                                  ppm_loadbal_isendlist(i)) THEN
                                        !---------------------------------------
                                        !  This isub can be sent since its
                                        !  neighbors' procs are also underloaded
                                        !  proc's neighbors OR jsub is already
                                        !  on the underloaded proc
                                        !---------------------------------------
                                        isNewNeighbor(isub,jsub) = .FALSE.
                                    ENDIF
                                ENDDO! jproc -> iterate over RECV's sub's procs
                                IF (.NOT. topo%sub2proc(topo%ineighsubs(jsub,isub)).EQ.&
     &                              ppm_rank) THEN
                                    isEngulfed(isub,jsub) = .FALSE.
                                ENDIF
                            ENDDO! jsub->iterate over neighboring subs
                            !---------------------------------------------------
                            !  Add this sub to the list of possible candidates
                            !  to be sent to the underloaded proc. IFF none of
                            !  its neighboring subs (jsub) have the same proc.
                            !  candidate_sublist keeps local sub IDs
                            !---------------------------------------------------
!                            print*,ppm_rank,isub,isNewNeighbor(isub,:)
                            IF (.NOT. ANY(isNewNeighbor(isub,:)) .AND.   &
     &                          .NOT. ALL(isEngulfed(isub,:))  ) THEN
                                counter = counter + 1
                                candidate_sublist(counter) = isub
                            ENDIF
                        ENDDO ! isub->iterate over all subs
                        !---------------------------------------------------------
                        !  If no candidate subs are eligible then skip DLB
                        !---------------------------------------------------------
                        IF (counter .EQ. 0) GOTO 7777
                        !---------------------------------------------------------
                        !  Now, I gotta choose which sub to send. More on this is
                        !  explained in ppm_loadbal_choose_sub()
                        !  candidate_sublist contains LOCAL IDs of the subs
                        !---------------------------------------------------------

                        ! stdout("CHOOSING a SUB")
!                        stdout("candidate list:",candidate_sublist)
!                        print*,ppm_rank,'candidates:',candidate_sublist
                        CALL ppm_loadbal_choose_sub(topoid,ideal_load,candidate_sublist,&
     &                                            send_sublist,counter,info)
                        or_fail("ppm_loadbal_choose_sub")
                        !---------------------------------------------------------
                        !  Communicate list of to-be-sent subs w/ underloaded proc
                        !  First : Size.......
                        !---------------------------------------------------------
                        send_subsize = SIZE(send_sublist,1)

!                        stdout("SIZE(send_sublist):",send_subsize)
!                        IF (send_sublist(1).NE.-1) THEN
!                            stdout("sublistID to be sent:",'send_sublist')
!                        ELSE
!                            stdout("sending nothing this time")
!                        ENDIF
                        tag1 = 126
                        CALL MPI_Send(send_subsize,1,MPI_INTEGER,ppm_loadbal_isendlist(i),&
     &                              tag1,ppm_comm,status,info)
                        or_fail("MPI_Send::send_sublist")
                        !---------------------------------------------------------
                        !  If no subs are eligible then skip DLB
                        !---------------------------------------------------------
                        IF (send_subsize .EQ. 0) GOTO 7777
                        !---------------------------------------------------------
                        !  Second: Actual list
                        !  send_sublist contains global IDs of subs
                        !---------------------------------------------------------
                        tag1 = 127

                        CALL MPI_Send(send_sublist,send_subsize,MPI_INTEGER,  &
     &                              ppm_loadbal_isendlist(i),tag1,ppm_comm,status,info)
                        or_fail("MPI_Send::send_sublist")
                        !---------------------------------------------------------
                        !  If the list has only one subID and that is an invalid
                        !  one (i.e. -1), it means I have nothing to send to RECV
                        !  Thus, I can skip this neighbor in this round
                        !---------------------------------------------------------
                        IF (send_subsize.LE.1 .AND.send_sublist(1).EQ.-1) THEN
                            stdout("Nothing to send....")
                            GOTO 7777
                        !-------------------------------------------------------
                        !  ... ELSE I send subdomains one after another at most
                        !  one sub is remaining
                        !-------------------------------------------------------
                        ELSE
!                            stdout("send_sublist(isub)",'send_sublist(1)')
                            DO isub=1,send_subsize
                                IF (send_sublist(isub).NE.-1) THEN
                                    CALL ppm_loadbal_map_subpart(topoid,Pc,isSender,&
     &                         send_sublist(isub),ppm_rank,ppm_loadbal_isendlist(i),info)
                                    or_fail("ppm_loadbal_map_subpart--sender")
                                ENDIF
                            ENDDO
                        ENDIF
                    ELSE
                        !-------------------------------------------------------
                        !  I'm the underloaded proc.
                        !  I need to send nneighproc and my neighboring processor
                        !  list so that the sender can select correct sub(s) to
                        !  send me.
                        !  First: size
                        !-------------------------------------------------------
                        tag1 = 124
                        CALL MPI_Send(topo%nneighproc,1,MPI_INTEGER,            &
     &                                ppm_loadbal_irecvlist(i),tag1,ppm_comm,status, &
     &                                info)
                        or_fail("topo%nneighproc MPI_Send failed")
                        stdout("topo%nneighproc",'topo%nneighproc')
!                        stdout("MPI_SEND----OK")
                        !-------------------------------------------------------
                        !  Second: actual list
                        !-------------------------------------------------------
                        tag1 = 125

                        CALL MPI_Send(topo%ineighproc,topo%nneighproc,MPI_INTEGER, &
     &                                ppm_loadbal_irecvlist(i),tag1,ppm_comm,status, &
     &                                info)
                        or_fail("topo%ineighproc MPI_Send failed")
!                        stdout("MPI_SEND2----OK")
!                        stdout("# of neighs",'topo%nneighproc')
!                        stdout("my neighs",'topo%ineighproc(1:topo%nneighproc)')
                        !-------------------------------------------------------
                        !  This is the part where the SENDER should be deciding
                        !  which sub to choose. (like waiting for girlfriend to
                        !  choose what to wear)
                        !  .............
                        !  OK, now I can get the list of the subs I should recv
                        !  First : Size of the list...
                        !------------------------------------------------------
!                        print*,ppm_rank,'topo%sub2proc',topo%sub2proc
                        tag1 = 126
                        CALL MPI_Recv(recv_subsize,1,MPI_INTEGER,ppm_loadbal_irecvlist(i),&
     &                                tag1,ppm_comm,status,info)
                        or_fail("MPI_Recv::recv_subsize")
                        !---------------------------------------------------------
                        !  If no subs are eligible then skip DLB
                        !---------------------------------------------------------
                        IF (recv_subsize .EQ. 0) GOTO 7777
                        !-------------------------------------------------------
                        !  ... allocate recv_sublist...
                        !-------------------------------------------------------
                        ldu(1) = recv_subsize
                        ldu(2) = 0
                        iopt   = ppm_param_alloc_fit
                        CALL ppm_alloc(recv_sublist,ldu,iopt,info)
                        or_fail_alloc("recv_sublist")
                        !-------------------------------------------------------
                        !  ... now the actual list (global sub IDs)
                        !-------------------------------------------------------
                        tag1 = 127
                        CALL MPI_Recv(recv_sublist,recv_subsize,MPI_INTEGER,  &
     &                                ppm_loadbal_irecvlist(i),tag1,ppm_comm,status,info)
                        or_fail("MPI_Recv::recv_sublist")
                        stdout("recv_sublist",recv_sublist)
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
                                    stdout("from",'ppm_loadbal_irecvlist(i)',"to",ppm_rank)
                                    CALL ppm_loadbal_map_subpart(topoid,Pc,isSender,&
     &                          recv_sublist(isub),ppm_loadbal_irecvlist(i),ppm_rank,info)
                                    or_fail("ppm_loadbal_map_subpart")
                                ENDIF
                            ENDDO
                        ENDIF
                    ENDIF ! isSender?
                ENDIF ! I'm either a sender/receiver in this comm round
            ENDDO !kk =1,nsteps

!            CALL Pc%map_ghosts(info)
!            or_fail("ghost mapping failed")

 7777       CONTINUE
!            CALL MPI_Wait(status,ppm_comm,info)
!            or_fail("MPI_Wait failed")
          ENDDO
      ENDDO !threshold to reach k-close solution
      !-------------------------------------------------------------------------
      !  Not only during DLB but also after DLB ppm_isendlist & ppm_irecvlist
      !  need to be converted back to their original settings
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_loadbal_nsendlist
      CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
       or_fail_alloc("ppm_isendlist")
      CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
       or_fail_alloc("ppm_irecvlist")
      ppm_isendlist = ppm_loadbal_isendlist
      ppm_irecvlist = ppm_loadbal_irecvlist
      ppm_nsendlist = ppm_loadbal_nsendlist
      !-------------------------------------------------------------------------
      !  Deallocate arrays
      !-------------------------------------------------------------------------
      IF (isSender) THEN
        DEALLOCATE(n_ineighproc,candidate_sublist,isNewNeighbor,isEngulfed,STAT=info)
        or_fail_dealloc("n_ineighproc,candidate_sublist,isNewNeighbor,isEngulfed")
        IF (ASSOCIATED(send_sublist)) DEALLOCATE(send_sublist,STAT=info)
        or_fail_dealloc("send_sublist")
      ELSE
        IF (ASSOCIATED(recv_sublist)) DEALLOCATE(recv_sublist,STAT=info)
        or_fail_dealloc("recv_sublist")
      ENDIF
      DEALLOCATE(ppm_loadbal_isendlist,ppm_loadbal_irecvlist,STAT=info)
      or_fail_dealloc("ppm_loadbal_irecvlist or ppm_loadbal_isendlist")
      stdout("========DLB ENDED========")
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

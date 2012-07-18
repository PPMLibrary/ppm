      minclude ppm_header(ppm_module_name="loadbal_inquire_dlb")

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE loadbal_inquire_dlb_s(topoid,t_comp,t_comm,my_mat, &
     &                                 max_ctime,imbal_perc,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE loadbal_inquire_dlb_d(topoid,t_comp,t_comm,my_mat, &
     &                                 max_ctime,imbal_perc,info)
#endif
      !!! This routine inquires about the load balance status using XXX dynamic
      !!! load balancing decision policy and returns advise that the total
      !!! workload will be better balanced if ppm_dynamic_loadbal is called.
      !!!
      !!! [NOTE]
      !!! This estimate is currently a scalar, so only one topology/set of
      !!! topologies can be monitored. Introducing a topology set ID and making
      !!! the internal estimates a vector, this can later be extended to
      !!! multiple topology sets if needed. The topology set
      !!! ID for which to evaluate the load balance will then be an additional
      !!! argument to this routine. This routine uses only one-to-one MPI
      !!! operations and thus, it is scalable.
      !!! Also, the timings here are from the previous time step. The change
      !!! in the number of particles from the previous step is incorporated
      !!! in the estimation of the current step's timings. By separating
      !!! computation and communication times from each other, we are able to
      !!! postpone the communication (partial mappings) so that once we decide
      !!! on who sends whom and which subdomain, we can move the subdomain
      !!! and have only one round of partial mapping. Thus, we will be able to
      !!! save time.
      !!! One more thing: Calling DLB in the very first time step would yield
      !!! to severe problems.
      !!!
      !!! .References
      !!! *************************************************************
      !!! "Detecting Application Load Imbalance on High End Massively
      !!! Parallel Systems" by DeRose et al.
      !!! 2007
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
      INTEGER                , INTENT(IN   ) :: topoid
      !!! The topology where the DLB will take place
      REAL(MK)               , INTENT(IN   ) :: t_comp
      !!! Elapsed time (as measured by `ppm_time`) for all computation in one
      !!! time step on the local processor.
      REAL(MK)               , INTENT(IN   ) :: t_comm
      !!! Elapsed time for all COMMUNICATION in one time step on the local proc
      REAL(MK)               , INTENT(IN   ) :: my_mat
      !!! Moving (running) average time of the simulation on this processor.
      !!! This needs to be communicated with the neighboring processors to get
      !!! the max of the my_mat in this local neighborhood
      REAL(MK)               , INTENT(  OUT) :: max_ctime
      !!! The longest time needed by the slowest processor in my neighborhood
      REAL(MK)               , INTENT(  OUT) :: imbal_perc
      !!! The load imbalance percentage of the neighborhood (me + my neighbors)
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                           :: my_imbal,my_imbal_perc
      REAL(MK)                           :: my_weighted_time,max_wtime,min_wtime
      REAL(MK)                           :: avg_mat,x,my_ctime
      REAL(MK)                           :: curr_min,curr_max,npart_ratio
      REAL(MK),DIMENSION(:),POINTER      :: ctime=>NULL()
      REAL(MK),DIMENSION(:),POINTER      :: mat=>NULL()
      ! mat keeps the moving average times of all neighbors
      REAL(MK),DIMENSION(:),POINTER      :: weighted_time=>NULL()
      INTEGER,DIMENSION(1)               :: temp
      INTEGER                            :: ndiff,i,k,iset,ibuffer,tag1
      INTEGER                            :: nneighs,groupsize
      TYPE(ppm_t_topo),POINTER           :: topo => NULL()
#ifdef __MPI
      INTEGER                            :: MPTYPE
      INTEGER, DIMENSION(MPI_STATUS_SIZE):: status
#endif
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      start_subroutine("loadbal_inquire_dlb")

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
      !  Allocate and initialize memory for time arrays
      !-------------------------------------------------------------------------
!      stdout("ppm_nsendlist:",ppm_nsendlist," nneighs:",nneighs)
      ALLOCATE(ctime(1:ppm_nproc),STAT=info)
               or_fail_alloc("ctime")
      ALLOCATE(mat(1:ppm_nproc),STAT=info)
               or_fail_alloc("mat")
      ALLOCATE(weighted_time(1:ppm_nproc),STAT=info)
               or_fail_alloc("weighted_time")

      imbal_perc    =-1._MK
      max_ctime     =-1._MK
      ctime         = 0._MK
      mat           = 0._MK
      weighted_time = 0._MK

      !-------------------------------------------------------------------------
      !  The ratio of the numbers of particles in consequent iterations gives
      !  the approximation for the increase/decrease in the t_comm
      !-------------------------------------------------------------------------
!      stdout(ppm_loadbal_npart_new,ppm_loadbal_npart_old)
      npart_ratio = ppm_loadbal_npart_new/ppm_loadbal_npart_old
!      IF (npart_ratio .GT. 10) THEN
!        stdout("Warning! Suddenly my workload has become 10 times more!! A bug?")
!      ENDIF

      !-------------------------------------------------------------------------
      !  My total iteration time at step i is estimated by time at i-1 *
      !  the change in the number of particles
      !  I cannot update moving average time likewise since I need to have a
      !  'moving average number of particles' info
      !-------------------------------------------------------------------------
      my_ctime = (t_comp + t_comm)*npart_ratio

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
      !  Loop over the neighboring processors according to the optimized
      !  communication sequence.
      !-------------------------------------------------------------------------
      groupsize = 0
      DO i=2,ppm_nsendlist
        !----------------------------------------------------------------------
        !  only consider non-negative send/recv ranks
        !----------------------------------------------------------------------
        IF (ppm_isendlist(i).GE.0 .AND. ppm_irecvlist(i).GE.0) THEN
            tag1 = 123
            !-------------------------------------------------------------------
            !  Send & receive the elasped time of the last time step
            !  +1 is needed because rank 0's ctime will be in ctime(1)
            !-------------------------------------------------------------------
            CALL MPI_SendRecv(my_ctime,1,MPTYPE,ppm_isendlist(i),tag1, &
     &                        ctime(ppm_irecvlist(i)+1),1,MPTYPE,ppm_irecvlist(i), &
     &                        tag1,ppm_comm,status,info)
            or_fail("Sendrecv ctime problem!")
            !-------------------------------------------------------------------
            !  Send & receive the moving average times
            !-------------------------------------------------------------------
            tag1 = 234
            CALL MPI_SendRecv(my_mat,1,MPTYPE,ppm_isendlist(i),tag1, &
     &                        mat(ppm_irecvlist(i)+1),1,MPTYPE,ppm_irecvlist(i), &
     &                        tag1,ppm_comm,status,info)
            or_fail("Sendrecv MAT problem!")

        ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Compute the average of the moving average times of neighborhood
      !  (incl. myself, thus +1)
      !-------------------------------------------------------------------------
      groupsize = nneighs + 1
      avg_mat   = REAL( (SUM(mat)+my_mat) / groupsize)
      !stdout_f('(A,I3,A,F12.4)',"groupsize=",groupsize," avg_mat=",avg_mat)

      !-------------------------------------------------------------------------
      !  Compute max(ctime)
      !-------------------------------------------------------------------------
      IF (MAXVAL(ctime) .GT. my_ctime) THEN
         max_ctime = MAXVAL(ctime)
      ELSE
         max_ctime = my_ctime
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute current load imbalance percentage according to DeRose:2007
      !-------------------------------------------------------------------------
      imbal_perc = (ABS(max_ctime-avg_mat) / max_ctime)*(groupsize/(groupsize-1))

      !-------------------------------------------------------------------------
      !  Get a weighted computation time: (my_ctime + my_mat) / 2
      !  The one with the longest weighted elapsed time should send to the
      !  shortest one.
      !-------------------------------------------------------------------------
      DO i=1,ppm_nproc
         weighted_time(i) = (ctime(i) + mat(i)) / 2._MK
      ENDDO
      !-------------------------------------------------------------------------
      !  Put my weighted_time also in the list
      !-------------------------------------------------------------------------
      my_weighted_time    = (my_ctime + my_mat) / 2._MK
      weighted_time(ppm_rank+1) = my_weighted_time
!      stdout("weighted_times",weighted_time)


      !-------------------------------------------------------------------------
      !  Let's decide now which proc should send to the other one
      !  Zero times are for the processors which are not my neighbors
      !  They must be ignored.
      !-------------------------------------------------------------------------
      curr_min = HUGE(x)
      curr_max = 0._MK
      DO i=1,ppm_nproc
        IF (weighted_time(i) .GT. 0._MK .AND.                         &
     &      weighted_time(i) .LT. curr_min ) THEN

            curr_min = weighted_time(i)
            ppm_loadbal_recvrank = i-1

        ENDIF
        IF (weighted_time(i) .GT. curr_max ) THEN

            curr_max = weighted_time(i)
            ppm_loadbal_sendrank = i-1

        ENDIF
      ENDDO

!      stdout("Recv is:",ppm_loadbal_recvrank," Sender is:",ppm_loadbal_sendrank)

      !-------------------------------------------------------------------------
      !  Make sure that I'm not the guy with the most and least workload at the
      !  same time
      !-------------------------------------------------------------------------
      IF (ppm_loadbal_sendrank.EQ.ppm_loadbal_recvrank .AND. &
     &    ppm_loadbal_sendrank.NE.-1) THEN
        fail("This process has both the most and least workload at the same &
     &        time in its neighborhood. This cannot happen!")
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate local pointers
      !-------------------------------------------------------------------------
      dealloc_pointer(ctime)
      dealloc_pointer(mat)
      dealloc_pointer(weighted_time)

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

        IF (my_ctime .LT. 0.0_MK) THEN
            fail("my_ctime must be >= 0.0",exit_point=8888)
        ENDIF

 8888   CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE loadbal_inquire_dlb_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE loadbal_inquire_dlb_d
#endif

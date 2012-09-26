      minclude ppm_header(ppm_module_name="loadbal_bc")

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE loadbal_bc_s(topoid,my_time,info,t_comp,t_comm)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE loadbal_bc_d(topoid,my_time,info,t_comp,t_comm)
#endif
      !!! This routine executes the "Balancing Circuit" model of dynamic
      !!! load balancing. In this model every processor follows its list
      !!! of communication sequence and equalizes its load with its neighbor.
      !!! Sauerwald:2012 describes this model and they derive the number of
      !!! steps needed for the algorithm to reach a K-close solution.
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
      REAL(MK)               , INTENT(IN   ) :: my_time
      !!! The total computation time for this time step
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      REAL(MK), OPTIONAL     , INTENT(IN   ) :: t_comp
      !!! Elapsed time (as measured by `ppm_time`) for all computation in one
      !!! time step on the local processor.
      REAL(MK), OPTIONAL     , INTENT(IN   ) :: t_comm
      !!! Elapsed time for all COMMUNICATION in one time step on the local proc
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                           :: neighbor_time
      ! Neighbor's elapsed time in this iteration
      REAL(MK),DIMENSION(:),POINTER      :: weighted_time=>NULL()
      INTEGER,DIMENSION(1)               :: temp
      INTEGER                            :: ndiff,i,k,iset,ibuffer,tag1
      INTEGER                            :: nneighs,groupsize,nsteps
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
      start_subroutine("loadbal_bc")

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
      !  The optimized communication sequence gives us the pre-described
      !  matching sequence needed by the balancing circuit algorithm.
      !  The number of matching matrices is equal to the number of colors in the
      !  edge coloring algorithm of the communication graph.
      !  This process has to be repeated until we are K-close to an ideal soln.
      !-------------------------------------------------------------------------
      nsteps = 2
      ! TODO: REPLACE THIS W/ A WHILE-LOOP BY USING A THRESHOLD
      DO k=1,nsteps
      groupsize = 0
          DO i=2,ppm_nsendlist
            !-------------------------------------------------------------------
            !  only consider non-negative send/recv ranks
            !-------------------------------------------------------------------
            IF (ppm_isendlist(i).GE.0 .AND. ppm_irecvlist(i).GE.0) THEN
                tag1 = 123
                !---------------------------------------------------------------
                !  Send & receive the elasped time of the last time step
                !
                !---------------------------------------------------------------
                CALL MPI_SendRecv(my_time,1,MPTYPE,ppm_isendlist(i),tag1, &
         &                        neighbor_time,1,MPTYPE,ppm_irecvlist(i), &
         &                        tag1,ppm_comm,status,info)
                or_fail("Sendrecv ctime problem!")
            ENDIF
            stdout("my time:",my_time,", ",'ppm_irecvlist(i)'," time:",neighbor_time)
            !-------------------------------------------------------------------
            !
            !-------------------------------------------------------------------


          ENDDO
      ENDDO


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

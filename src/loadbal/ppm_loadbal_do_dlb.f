      minclude ppm_header(ppm_module_name="loadbal_do_dlb")

      SUBROUTINE loadbal_do_dlb(topoid,Pc,info)

      !!! Before calling this function ppm_loadbal_inquire with
      !!! ppm_param_loadbal_dlb heuristic has to be called. This ensures that
      !!! ppm_loadbal_sendrank and ppm_loadbal_recvrank are set to correct
      !!! values.
      !!! [TODO]
      !!! If a proc has only one subdomain available, it shouldn't send it
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
      USE ppm_module_particles_typedef
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
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
      !!! Topology ID
      TYPE(ppm_t_particles_d)   , INTENT(IN   ) :: Pc
      !!! Particle set
      INTEGER                 , INTENT(INOUT) :: info
      !!! return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER                            :: i,j,isub,jsub,iproc,nsubs_of_isub
      INTEGER                            :: neigh_index,jsubs_proc,sent_sub
      INTEGER                            :: maxneigh,tag1,tag2,temp_recv
      INTEGER                            :: dlb_s,dlb_r
      LOGICAL                            :: DLBdone=.FALSE.
      TYPE(ppm_t_topo),POINTER           :: topo => NULL()
      INTEGER, DIMENSION(MPI_STATUS_SIZE):: status

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      start_subroutine("loadbal_do_dlb")
      temp_recv = -1
      sent_sub = -1
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Get the topology first
      !-------------------------------------------------------------------------
      topo   => ppm_topo(topoid)%t

      !-------------------------------------------------------------------------
      !  Depending on the outcome from the ppm_loadbal_inquire_dlb, I know
      !  whether I am the processor that needs to send/recv subdomains.
      !  I check first whether I'm the processor that is overloaded and
      !  in need to send subdomains
      !-------------------------------------------------------------------------
      !stdout("ineighsubs:",'SIZE(topo%ineighsubs,1)','SIZE(topo%ineighsubs,2)')
!      stdout("no of subs on this proc",'topo%nsublist')
!      stdout("max no of neighboring subs",'topo%nneighsubs(1)')
      IF (ppm_loadbal_sendrank .EQ. ppm_rank .AND. topo%nsublist .EQ. 1) THEN
        stdout("I have only one sub!")
      ENDIF
      DO i=2,ppm_nsendlist
        !-----------------------------------------------------------------------
        !  For now, pick only one subdomain and send it to my neighbor with
        !  least workload. Follow the optimized comm sequence
        !
        !  only consider non-negative send/recv ranks
        !-----------------------------------------------------------------------
        IF (ppm_isendlist(i).GE.0 .AND. ppm_irecvlist(i).GE.0) THEN
            !-------------------------------------------------------------------
            !  A handshake-like mechanism to ensure that the both the sender &
            !  receiver are on the same page. It means that sender and receiver
            !  should have the same ppm_loadbal_sendrank & ppm_loadbal_recvrank
            !-------------------------------------------------------------------
!            stdout("nsendlist iteration:",i)
!            stdout("->Recv is:",ppm_loadbal_recvrank," Sender is:",ppm_loadbal_sendrank)
            tag1 = 123
            CALL MPI_SendRecv(ppm_loadbal_recvrank,1,MPI_INTEGER,             &
     &            ppm_isendlist(i),tag1,ppm_loadbal_n_recvrank,1,MPI_INTEGER, &
     &            ppm_irecvlist(i),tag1,ppm_comm,status,info)
            or_fail("mpi_sendrecv of ppm_loadbal_n_recvrank has errors")
            tag1 = 321
            CALL MPI_SendRecv(ppm_loadbal_sendrank,1,MPI_INTEGER,             &
     &            ppm_isendlist(i),tag1,ppm_loadbal_n_sendrank,1,MPI_INTEGER, &
     &            ppm_irecvlist(i),tag1,ppm_comm,status,info)
            or_fail("mpi_sendrecv of ppm_loadbal_n_sendrank has errors")
!            stdout("N-Recv:",ppm_loadbal_n_recvrank," N-Sender:",ppm_loadbal_n_sendrank)
            stdout("Recv  :",ppm_loadbal_recvrank,"   Sender:",ppm_loadbal_sendrank)
            !-------------------------------------------------------------------
            !  If I'm not the sender or the receiver in the neighborhood, I pass
            !-------------------------------------------------------------------
            IF (ppm_loadbal_recvrank .NE. ppm_rank .AND.       &
     &          ppm_loadbal_sendrank .NE. ppm_rank ) THEN

                CYCLE
            !-------------------------------------------------------------------
            !  Else if I'm one either receiver or sender, I participate
            !-------------------------------------------------------------------
            ELSE IF (ppm_loadbal_n_recvrank .EQ. ppm_loadbal_recvrank .AND.  &
     &               ppm_loadbal_n_sendrank .EQ. ppm_loadbal_sendrank) THEN
                !---------------------------------------------------------------
                !  If I'm the receiver, I shall receive a subdomain
                !---------------------------------------------------------------
                IF (ppm_loadbal_recvrank .EQ. ppm_rank) THEN
                    stdout("receiver!")
                    CALL ppm_loadbal_recvsub(topoid,topo%prec,info)
                        or_fail("ppm_loadbal_RECVsub has errors")
                !---------------------------------------------------------------
                !  If I'm the sender, I shall send a subdomain
                !---------------------------------------------------------------
                ELSE IF (ppm_loadbal_sendrank .EQ. ppm_rank) THEN
                    !-----------------------------------------------------------
                    !  Iterate over my subs
                    !-----------------------------------------------------------
                    stdout("sender!")
                    DO isub=1,topo%nsublist
                    !-----------------------------------------------------------
                    !  Iterate over the neighboring subs of the isub
                    !-----------------------------------------------------------
                        DO neigh_index=1,topo%nneighsubs(isub)
                        !-------------------------------------------------------
                        !  Check my neighbor sub's (neigh_index) proc assignment.
                        !  If it belongs to this underloaded processor, this
                        !  also means that my sub (isub) is on the
                        !  boundary between me and him. Thus, it can be sent to
                        !  this underloaded processor
                        !-------------------------------------------------------
!                stdout("global IDs of neighboring subs:",'topo%ineighsubs(:,isub)')
!                stdout("isub",isub,"neighbor index",neigh_index)
!                stdout("global ID of this local sub",'topo%isublist(neigh_index)')
!                stdout("proc of jsub",'topo%sub2proc(topo%isublist(neigh_index))')
                            nsubs_of_isub = SIZE(topo%ineighsubs(:,isub))
                            jsubs_proc = topo%sub2proc(            &
     &                                      topo%ineighsubs(neigh_index,isub))
!                            stdout("neighbor sub's process ID:",jsubs_proc)
                            IF (ppm_loadbal_recvrank .EQ. jsubs_proc) THEN
                            !---------------------------------------------------
                            !  I found a neighbor sub that is assigned to another
                            !  processor. So, I found also a sub that is on the
                            !  boundary with that proc
                            !
                            !  Now, send this isub to that neighbor.
                            !  NOTE: Currently, only one sub is sent. So exit
                            !
                            !----------------------------------------------------
                              maxneigh = SIZE(topo%ineighsubs,1)
                              CALL ppm_loadbal_sendsub(topoid,Pc,isub,&
     &                                                 maxneigh,topo%prec,info)
                                  or_fail("ppm_loadbal_SENDsub has errors")
                              sent_sub = topo%isublist(isub)
                              GOTO 7777
                            ENDIF ! ppm_loadbal_recvrank .EQ. jsubs_proc
                        ENDDO
                    ENDDO
                ENDIF ! ppm_loadbal_recvrank .EQ. ppm_rank
            ENDIF ! I don't participate
        ENDIF ! only consider non-negative ranks
      ENDDO
      !-------------------------------------------------------------------------
      !  If DLB is done, update the sub2proc assignment as well
      !-------------------------------------------------------------------------
 7777 CONTINUE
       tag2 = 888
!      DO i=2,ppm_nsendlist
!        !-----------------------------------------------------------------------
!        !  For now, pick only one subdomain and send it to my neighbor with
!        !  least workload. Follow the optimized comm sequence
!        !
!        !  only consider non-negative send/recv ranks
!        !-----------------------------------------------------------------------
!        IF (ppm_isendlist(i).GE.0 .AND. ppm_irecvlist(i).GE.0) THEN
!            IF ( ppm_loadbal_sendrank .EQ. ppm_rank) THEN
!                topo%sub2proc(sent_sub) = ppm_loadbal_recvrank
!                stdout("sending...")
!                DO j=1,topo%nneighproc
!                    CALL MPI_Send(sent_sub,1,MPI_INTEGER,topo%ineighproc(j),tag2,    &
!     &                  ppm_comm,status,info)
!                    or_fail("send of sub2proc has errors")
!                ENDDO
!                stdout("sent!!")
!            ELSE
!                stdout("receiving...")
!                CALL MPI_Recv(sent_sub,1,MPI_INTEGER,ppm_loadbal_sendrank,tag2,      &
!     &              ppm_comm,status,info)
!                or_fail("recv of sub2proc has errors")
!                topo%sub2proc(sent_sub) = ppm_loadbal_recvrank
!                stdout("received!!")
!            ENDIF
!        ENDIF
!      ENDDO
!


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



 8888   CONTINUE
      END SUBROUTINE check

      END SUBROUTINE loadbal_do_dlb


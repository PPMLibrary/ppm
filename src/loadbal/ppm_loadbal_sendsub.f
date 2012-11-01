      minclude ppm_header(ppm_module_name="loadbal_sendsub")

!#if   __KIND == __SINGLE_PRECISION
!      SUBROUTINE loadbal_sendsub_s(topoid,isub,receiver,maxneigh,prec,info)
!#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE loadbal_sendsub(topoid,isub,receiver,maxneigh,prec,info)
!#endif
      !!! This routine is called by the overloaded process. One subdomain is
      !!! sent to the most underloaded process.
      !!! Particles within the subdomain need to be sent as well.
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
      USE ppm_module_interfaces
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
!#if    __KIND == __SINGLE_PRECISION
!      TYPE(ppm_t_particles_s) , INTENT(IN   ) :: Pc
!#else
!      TYPE(ppm_t_particles_d) , INTENT(IN   ) :: Pc
!#endif
      !!! Particle set
      INTEGER                 , INTENT(IN   ) :: receiver
      !!! MPI rank of the receiving underloaded process
      INTEGER                 , INTENT(IN   ) :: isub
      !!! The local ID of the subdomain to be sent
      INTEGER                 , INTENT(IN   ) :: maxneigh
      !!! max number of neighbors of a subdomain on this topology
      INTEGER                 , INTENT(IN   ) :: prec
      !!! Precision. One of:
      !!! * ppm_kind_single
      !!! * ppm_kind_double
      INTEGER                 , INTENT(INOUT) :: info
      !!! return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER                            :: i,k,jsub,iproc,nprocs_of_this_sub
      INTEGER                            :: tag1,itemp,iopt,nsendlist
      INTEGER                            :: new_nsub,old_nsub,nneighsubs
      INTEGER,DIMENSION(2)               :: ldu2
      INTEGER,DIMENSION(1)               :: ldu1
      REAL(ppm_kind_single)              :: subcosts,temp1s
      REAL(ppm_kind_double)              :: subcostd,temp1d
      INTEGER                            :: numelm2=2
      REAL(ppm_kind_single),DIMENSION(2) :: min_subs2,max_subs2,temps2
      REAL(ppm_kind_double),DIMENSION(2) :: min_subd2,max_subd2,tempd2
      REAL(MK),DIMENSION(:,:),POINTER    :: xp => NULL()
      INTEGER,DIMENSION(:),POINTER       :: list_del_parts => NULL()

      INTEGER                            :: numelm3=3
      REAL(ppm_kind_single),DIMENSION(3) :: min_subs3,max_subs3,temps3
      REAL(ppm_kind_double),DIMENSION(3) :: min_subd3,max_subd3,tempd3

      TYPE(ppm_t_topo),POINTER           :: topo => NULL()
      INTEGER,DIMENSION(maxneigh),TARGET :: temp_neigh,ineighsubs
#ifdef __MPI
      INTEGER                            :: MPTYPE
      INTEGER, DIMENSION(MPI_STATUS_SIZE):: status
#endif
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      start_subroutine("loadbal_sendsub")
      topo => ppm_topo(topoid)%t
!      call Pc%get_xp(xp,info)
!       or_fail("Pc%get_xp failed")
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
      IF (prec .EQ. ppm_kind_single) THEN
        MPTYPE = MPI_REAL
      ELSE IF(prec .EQ. ppm_kind_double) THEN
        MPTYPE = MPI_DOUBLE_PRECISION
      ENDIF
#endif
      !-------------------------------------------------------------------------
      !  Current nsublist
      !-------------------------------------------------------------------------
      old_nsub =  topo%nsublist

      !-------------------------------------------------------------------------
      !  We need to have more than 1 processor so that DLB makes sense
      !-------------------------------------------------------------------------
      IF (ppm_nproc.EQ.1) THEN
        stdout("WARNING: you are using only one process. DLB will quit!")
        GOTO 9999
      ENDIF
!      stdout("local subID to be sent",isub)
      !-------------------------------------------------------------------------
      !  Copy subdomain-related information to send them.
      !-------------------------------------------------------------------------
!      IF (ppm_dim .EQ. 2) THEN
!          IF (prec .EQ. ppm_kind_single) THEN
!              min_subs2 = topo%min_subs(:,isub)
!              max_subs2 = topo%max_subs(:,isub)
!              subcosts  = topo%sub_costs( isub)
!          ELSE IF(prec .EQ. ppm_kind_double) THEN
!              min_subd2 = topo%min_subd(:,isub)
!              max_subd2 = topo%max_subd(:,isub)
!              subcostd  = topo%sub_costd( isub)
!          ENDIF
!      ELSE IF (ppm_dim .EQ. 3) THEN
!          IF (prec .EQ. ppm_kind_single) THEN
!              min_subs3 = topo%min_subs(:,isub)
!              max_subs3 = topo%max_subs(:,isub)
!              subcosts  = topo%sub_costs( isub)
!          ELSE IF(prec .EQ. ppm_kind_double) THEN
!              min_subd3 = topo%min_subd(:,isub)
!              max_subd3 = topo%max_subd(:,isub)
!              subcostd  = topo%sub_costd( isub)
!          ENDIF
!      ENDIF
      ineighsubs = 0
      temp_neigh = 0
!      stdout("ineighsubs of isub BEFORE",ineighsubs)
      nneighsubs = topo%nneighsubs(isub)
!      stdout("size of topo%ineighsubs(:,isub)",'size(topo%ineighsubs(:,isub),1)')
      ineighsubs = topo%ineighsubs(:,isub)
!      stdout("ineighsubs of isub AFTER",ineighsubs)
      !-------------------------------------------------------------------------
      !  Send subdomain info to the most underloaded processor
      !-------------------------------------------------------------------------
!      IF (ppm_dim .EQ. 2) THEN
!          IF (prec .EQ. ppm_kind_single) THEN
!            tag1 = 400
!            CALL MPI_Send(min_subs2,numelm2,MPTYPE,receiver,tag1,&
!     &                    ppm_comm,status,info)
!                or_fail("min_sub send failed!")
!            tag1 = 500
!            CALL MPI_Send(max_subs2,numelm2,MPTYPE,receiver,tag1,&
!     &                    ppm_comm,status,info)
!                or_fail("max_sub send failed!")
!            tag1 = 600
!            CALL MPI_Send(subcosts,1,MPTYPE,receiver,tag1,&
!     &                    ppm_comm,status,info)
!                or_fail("subcost send failed!")
!
!          ELSE IF(prec .EQ. ppm_kind_double) THEN
!            tag1 = 400
!            CALL MPI_Send(min_subd2,numelm2,MPTYPE,receiver,tag1,&
!     &                    ppm_comm,status,info)
!                or_fail("min_sub send failed!")
!            tag1 = 500
!            CALL MPI_Send(max_subd2,numelm2,MPTYPE,receiver,tag1,&
!     &                    ppm_comm,status,info)
!                or_fail("max_sub send failed!")
!            tag1 = 600
!            CALL MPI_Send(subcostd,1,MPTYPE,receiver,tag1,&
!     &                    ppm_comm,status,info)
!                or_fail("subcost send failed!")
!          ENDIF
!      ELSE IF (ppm_dim .EQ. 3) THEN
!          IF (prec .EQ. ppm_kind_single) THEN
!            tag1 = 400
!            CALL MPI_Send(min_subs3,numelm3,MPTYPE,receiver,tag1,ppm_comm,&
!     &                    status,info)
!                or_fail("min_sub send failed!")
!            tag1 = 500
!            CALL MPI_Send(max_subs3,numelm3,MPTYPE,receiver,tag1,ppm_comm,&
!     &                    status,info)
!                or_fail("max_sub send failed!")
!            tag1 = 600
!            CALL MPI_Send(subcosts,1,MPTYPE,receiver,tag1,ppm_comm,&
!     &                    status,info)
!                or_fail("subcost send failed!")
!
!          ELSE IF(prec .EQ. ppm_kind_double) THEN
!            tag1 = 400
!            CALL MPI_Send(min_subd3,numelm3,MPTYPE,receiver,tag1,ppm_comm,&
!     &                    status,info)
!                or_fail("min_sub send failed!")
!            tag1 = 500
!            CALL MPI_Send(max_subd3,numelm3,MPTYPE,receiver,tag1,ppm_comm,&
!     &                    status,info)
!                or_fail("max_sub send failed!")
!            tag1 = 600
!            CALL MPI_Send(subcostd,1,MPTYPE,receiver,tag1,ppm_comm,&
!     &                    status,info)
!                or_fail("subcost send failed!")
!          ENDIF
!
!      ENDIF

      tag1 = 700
      CALL MPI_Send(nneighsubs,1,MPI_INTEGER,receiver,tag1, &
     &              ppm_comm,status,info)
        or_fail("nneighsubs send failed!")
      tag1 = 800
      CALL MPI_Send(ineighsubs,nneighsubs,MPI_INTEGER,receiver,&
     &              tag1,ppm_comm,status,info)
        or_fail("nneighsubs send failed!")
!      stdout("subdomain info sent")

      !-----------------------------------------------------------------------------
      !  Copy the sent subdomain to the end in the array so that ppm_alloc
      !  can clear it.
      !-----------------------------------------------------------------------------
!      stdout("isub:",isub," old_nsub:",old_nsub)
!      IF (ppm_dim .EQ. 2) THEN
!          IF (prec .EQ. ppm_kind_single) THEN
!              topo%min_subs(:,isub) = topo%min_subs(:,old_nsub)
!              topo%min_subs(:,old_nsub) = min_subs2
!
!              topo%max_subs(:,isub) = topo%max_subs(:,old_nsub)
!              topo%max_subs(:,old_nsub) = max_subs2
!
!              topo%sub_costs( isub) = topo%sub_costs( old_nsub)
!              topo%sub_costs(old_nsub) = subcosts
!          ELSE IF(prec .EQ. ppm_kind_double) THEN
!
!              topo%min_subd(:,isub) = topo%min_subd(:,old_nsub)
!              topo%min_subd(:,old_nsub) = min_subd2
!
!              topo%max_subd(:,isub) = topo%max_subd(:,old_nsub)
!              topo%max_subd(:,old_nsub) = max_subd2
!
!              topo%sub_costd( isub) = topo%sub_costd( old_nsub)
!              topo%sub_costd(old_nsub) = subcostd
!          ENDIF
!      ELSE IF (ppm_dim .EQ. 3) THEN
!          IF (prec .EQ. ppm_kind_single) THEN
!              topo%min_subs(:,isub) = topo%min_subs(:,old_nsub)
!              topo%min_subs(:,old_nsub) = min_subs3
!
!              topo%max_subs(:,isub) = topo%max_subs(:,old_nsub)
!              topo%max_subs(:,old_nsub) = max_subs3
!
!              topo%sub_costs( isub) = topo%sub_costs( old_nsub)
!              topo%sub_costs(old_nsub) = subcosts
!          ELSE IF(prec .EQ. ppm_kind_double) THEN
!
!              topo%min_subd(:,isub) = topo%min_subd(:,old_nsub)
!              topo%min_subd(:,old_nsub) = min_subd3
!
!              topo%max_subd(:,isub) = topo%max_subd(:,old_nsub)
!              topo%max_subd(:,old_nsub) = max_subd3
!
!              topo%sub_costd( isub) = topo%sub_costd( old_nsub)
!              topo%sub_costd(old_nsub) = subcostd
!          ENDIF
!
!      ENDIF
      !-------------------------------------------------------------------------
      !  New sublist has one LESS sub in it
      !-------------------------------------------------------------------------
      new_nsub =  old_nsub-1
      !-------------------------------------------------------------------------
      !  Continue updating local subdomain info by swaping
      !  isublist, nneighsubs, ineighsubs, nsublist
      !  Only if the local subID is the last sub in the isublist, then
      !  I don't need to do any special thing, ppm_param_alloc_fit_preserve
      !  will handle it for me
      !-------------------------------------------------------------------------
!      stdout("topo%ineighsubs(:,old_nsub) BEFORE",'topo%ineighsubs(:,old_nsub)')
!      stdout("topo%ineighsubs(:,    isub) BEFORE",'topo%ineighsubs(:,isub)')

!      stdout("old_nsub",old_nsub)
      itemp = topo%isublist(isub)
      topo%isublist(    isub) = topo%isublist(old_nsub)
      topo%isublist(old_nsub) = itemp

      itemp = topo%nneighsubs(isub)
      topo%nneighsubs(    isub) = topo%nneighsubs(old_nsub)
      topo%nneighsubs(old_nsub) = itemp


      temp_neigh = topo%ineighsubs(:,isub)
      topo%ineighsubs(:,    isub) = topo%ineighsubs(:,old_nsub)
      topo%ineighsubs(:,old_nsub) = temp_neigh
!
!      stdout("topo%ineighsubs(:,old_nsub) AFTER",'topo%ineighsubs(:,old_nsub)')
!      stdout("topo%ineighsubs(:,    isub) AFTER",'topo%ineighsubs(:,isub)')
!      stdout("topo%ineighsubs(:,new_nsub) AFTER",'topo%ineighsubs(:,new_nsub)')

      topo%nsublist = new_nsub
      stdout("Sub",'topo%isublist(old_nsub)'," sent: total subs:",'topo%nsublist')

      !-------------------------------------------------------------------------
      !  Remove the sent subdomain by shrinking the array by one element
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit_preserve
      ldu2(1) = maxneigh
      ldu2(2) = new_nsub
      ldu1(1) = new_nsub

!      IF(prec .EQ. ppm_kind_single) THEN
!          CALL ppm_alloc(topo%min_subs,ldu2,iopt,info)
!            or_fail_alloc("min_subs re-alloc failed!")
!          CALL ppm_alloc(topo%max_subs,ldu2,iopt,info)
!            or_fail_alloc("max_subs re-alloc failed!")
!          CALL ppm_alloc(topo%sub_costs,ldu1,iopt,info)
!            or_fail_alloc("sub_costs re-alloc failed!")
!
!      ELSE IF(prec .EQ. ppm_kind_double) THEN
!          CALL ppm_alloc(topo%min_subd,ldu2,iopt,info)
!            or_fail_alloc("min_subd re-alloc failed!")
!          CALL ppm_alloc(topo%max_subd,ldu2,iopt,info)
!            or_fail_alloc("max_subd re-alloc failed!")
!          CALL ppm_alloc(topo%sub_costd,ldu1,iopt,info)
!            or_fail_alloc("sub_costd re-alloc failed!")
!      ENDIF
      CALL ppm_alloc(topo%nneighsubs,ldu1,iopt,info)
        or_fail_alloc("nneighsubs re-alloc failed!")
      CALL ppm_alloc(topo%isublist,ldu1,iopt,info)
        or_fail_alloc("isublist re-alloc failed!")
      CALL ppm_alloc(topo%ineighsubs,ldu2,iopt,info)
        or_fail_alloc("ineighsubs re-alloc failed!")

!      stdout("topo%ineighsubs(:,new_nsub) VERY AFTER",'topo%ineighsubs(:,new_nsub)')
!      stdout("size of ineighsubs",'size()'
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
!#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE loadbal_sendsub



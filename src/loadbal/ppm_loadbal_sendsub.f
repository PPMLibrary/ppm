      minclude ppm_header(ppm_module_name="loadbal_sendsub")

      SUBROUTINE loadbal_sendsub(topoid,Pc,isub,maxneigh,prec,info)

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
      TYPE(ppm_t_particles_d) , INTENT(IN   ) :: Pc
      !!! Particle set
      INTEGER                 , INTENT(IN   ) :: isub
      !!! The local ID of the subdomain to be sent
      INTEGER                 , INTENT(IN   ) :: maxneigh
      !!! max number of neighbors of a subdomain
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
      call Pc%get_xp(xp,info)
       or_fail("Pc%get_xp failed")
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

      !-------------------------------------------------------------------------
      !  Copy subdomain-related information to send them.
      !-------------------------------------------------------------------------
      IF (ppm_dim .EQ. 2) THEN
          IF (prec .EQ. ppm_kind_single) THEN
              min_subs2 = topo%min_subs(:,isub)
              max_subs2 = topo%max_subs(:,isub)
              subcosts  = topo%sub_costs( isub)
          ELSE IF(prec .EQ. ppm_kind_double) THEN
              min_subd2 = topo%min_subd(:,isub)
              max_subd2 = topo%max_subd(:,isub)
              subcostd  = topo%sub_costd( isub)
          ENDIF
      ELSE IF (ppm_dim .EQ. 3) THEN
          IF (prec .EQ. ppm_kind_single) THEN
              min_subs3 = topo%min_subs(:,isub)
              max_subs3 = topo%max_subs(:,isub)
              subcosts  = topo%sub_costs( isub)
          ELSE IF(prec .EQ. ppm_kind_double) THEN
              min_subd3 = topo%min_subd(:,isub)
              max_subd3 = topo%max_subd(:,isub)
              subcostd  = topo%sub_costd( isub)
          ENDIF
      ENDIF

      nneighsubs = topo%nneighsubs(isub)
      ineighsubs = topo%ineighsubs(:,isub)
      !-------------------------------------------------------------------------
      !  Send subdomain info to the most underloaded processor
      !-------------------------------------------------------------------------
      IF (ppm_dim .EQ. 2) THEN
          IF (prec .EQ. ppm_kind_single) THEN
            tag1 = 400
            CALL MPI_Send(min_subs2,numelm2,MPTYPE,ppm_loadbal_recvrank,tag1,&
     &                    ppm_comm,status,info)
                or_fail("min_sub send failed!")
            tag1 = 500
            CALL MPI_Send(max_subs2,numelm2,MPTYPE,ppm_loadbal_recvrank,tag1,&
     &                    ppm_comm,status,info)
                or_fail("max_sub send failed!")
            tag1 = 600
            CALL MPI_Send(subcosts,1,MPTYPE,ppm_loadbal_recvrank,tag1,&
     &                    ppm_comm,status,info)
                or_fail("subcost send failed!")

          ELSE IF(prec .EQ. ppm_kind_double) THEN
            tag1 = 400
            CALL MPI_Send(min_subd2,numelm2,MPTYPE,ppm_loadbal_recvrank,tag1,&
     &                    ppm_comm,status,info)
                or_fail("min_sub send failed!")
            tag1 = 500
            CALL MPI_Send(max_subd2,numelm2,MPTYPE,ppm_loadbal_recvrank,tag1,&
     &                    ppm_comm,status,info)
                or_fail("max_sub send failed!")
            tag1 = 600
            CALL MPI_Send(subcostd,1,MPTYPE,ppm_loadbal_recvrank,tag1,&
     &                    ppm_comm,status,info)
                or_fail("subcost send failed!")
          ENDIF
      ELSE IF (ppm_dim .EQ. 3) THEN
          IF (prec .EQ. ppm_kind_single) THEN
            tag1 = 400
            CALL MPI_Send(min_subs3,numelm3,MPTYPE,ppm_loadbal_recvrank,tag1,ppm_comm,&
     &                    status,info)
                or_fail("min_sub send failed!")
            tag1 = 500
            CALL MPI_Send(max_subs3,numelm3,MPTYPE,ppm_loadbal_recvrank,tag1,ppm_comm,&
     &                    status,info)
                or_fail("max_sub send failed!")
            tag1 = 600
            CALL MPI_Send(subcosts,1,MPTYPE,ppm_loadbal_recvrank,tag1,ppm_comm,&
     &                    status,info)
                or_fail("subcost send failed!")

          ELSE IF(prec .EQ. ppm_kind_double) THEN
            tag1 = 400
            CALL MPI_Send(min_subd3,numelm3,MPTYPE,ppm_loadbal_recvrank,tag1,ppm_comm,&
     &                    status,info)
                or_fail("min_sub send failed!")
            tag1 = 500
            CALL MPI_Send(max_subd3,numelm3,MPTYPE,ppm_loadbal_recvrank,tag1,ppm_comm,&
     &                    status,info)
                or_fail("max_sub send failed!")
            tag1 = 600
            CALL MPI_Send(subcostd,1,MPTYPE,ppm_loadbal_recvrank,tag1,ppm_comm,&
     &                    status,info)
                or_fail("subcost send failed!")
          ENDIF

      ENDIF

      tag1 = 700
      CALL MPI_Send(nneighsubs,1,MPI_INTEGER,ppm_loadbal_recvrank,tag1, &
     &              ppm_comm,status,info)
        or_fail("nneighsubs send failed!")
      tag1 = 800
      CALL MPI_Send(ineighsubs,nneighsubs,MPI_INTEGER,ppm_loadbal_recvrank,&
     &              tag1,ppm_comm,status,info)
        or_fail("nneighsubs send failed!")
      stdout("subdomain info sent")
      !-------------------------------------------------------------------------
      !  Send the particles that are in this subdomain.
      !  Allocate the send buffer (ppm_loadbal_xp_) with size of Npart
      !  We will shrink the array once we know the exact number of particles
      !  in the subdomain
      !-------------------------------------------------------------------------
      IF (prec .EQ. ppm_kind_single) THEN
        ALLOCATE(ppm_loadbal_xps(ppm_dim*Pc%Npart),STAT=info)
               or_fail_alloc("ppm_loadbal_xps")
        ppm_loadbal_xps = 0._MK
      ELSE IF (prec .EQ. ppm_kind_double) THEN
        ALLOCATE(ppm_loadbal_xpd(ppm_dim*Pc%Npart),STAT=info)
               or_fail_alloc("ppm_loadbal_xpd")
        ppm_loadbal_xpd = 0._MK
      ENDIF
      ALLOCATE(list_del_parts(Pc%Npart),STAT=info)
            or_fail_alloc("list_del_parts failed")
      !-------------------------------------------------------------------------
      !  First determine the number of to-be-sent particles
      !-------------------------------------------------------------------------
      IF (ppm_dim.EQ.2) THEN
        IF (prec .EQ. ppm_kind_single) THEN
            !-----------------------------------------------------------------------
            !  Loop over the particles in 2D
            !-----------------------------------------------------------------------
            DO k=1,Pc%Npart
                !-------------------------------------------------------------------
                !  and check if they are inside the sub
                !-------------------------------------------------------------------
                IF (xp(1,k).GE.min_subs2(1).AND. &
     &              xp(1,k).LE.max_subs2(1).AND. &
     &              xp(2,k).GE.min_subs2(2).AND. &
     &              xp(2,k).LE.max_subs2(2)) THEN
                    !---------------------------------------------------------------
                    !  In the non-periodic case, allow particles that are
                    !  exactly ON an upper EXTERNAL boundary.
                    !---------------------------------------------------------------
                    IF((xp(1,k).LT.max_subs2(1)                  .OR.  &
     &                 (topo%subs_bc(2,topo%isublist(isub)).EQ.1      .AND. &
     &                  topo%bcdef(2).NE. ppm_param_bcdef_periodic))  .AND. &
     &                 (xp(2,k).LT.max_subs2(2)                  .OR.  &
     &                 (topo%subs_bc(4,topo%isublist(isub)).EQ.1      .AND. &
     &                  topo%bcdef(4).NE. ppm_param_bcdef_periodic))) THEN
                        !-----------------------------------------------------------
                        !  Increase the number of particles to-be-sent and add the
                        !  particle to the buffer
                        !-----------------------------------------------------------
                        ppm_loadbal_subnpart = ppm_loadbal_subnpart + 1
                        list_del_parts(ppm_loadbal_subnpart) = ppm_loadbal_subnpart
                        DO i=0,ppm_dim-1
                            ppm_loadbal_xps(ppm_loadbal_subnpart+i) = &
     &                                     xp(i+1,ppm_loadbal_subnpart)
                        ENDDO
                    ENDIF
                ENDIF
            ENDDO
        ELSE IF (prec .EQ. ppm_kind_double) THEN
            !-----------------------------------------------------------------------
            !  Loop over the particles in 2D
            !-----------------------------------------------------------------------
            DO k=1,Pc%Npart
                !-------------------------------------------------------------------
                !  and check if they are inside the sub
                !-------------------------------------------------------------------
                IF (xp(1,k).GE.min_subd2(1).AND. &
     &              xp(1,k).LE.max_subd2(1).AND. &
     &              xp(2,k).GE.min_subd2(2).AND. &
     &              xp(2,k).LE.max_subd2(2)) THEN
                    !---------------------------------------------------------------
                    !  In the non-periodic case, allow particles that are
                    !  exactly ON an upper EXTERNAL boundary.
                    !---------------------------------------------------------------
                    IF((xp(1,k).LT.max_subd2(1)                  .OR.  &
     &                 (topo%subs_bc(2,topo%isublist(isub)).EQ.1      .AND. &
     &                  topo%bcdef(2).NE. ppm_param_bcdef_periodic))  .AND. &
     &                 (xp(2,k).LT.max_subd2(2)                  .OR.  &
     &                 (topo%subs_bc(4,topo%isublist(isub)).EQ.1      .AND. &
     &                  topo%bcdef(4).NE. ppm_param_bcdef_periodic))) THEN
                        !-----------------------------------------------------------
                        !  Increase the number of particles to-be-sent and add the
                        !  particle to the buffer
                        !-----------------------------------------------------------

                        ppm_loadbal_subnpart = ppm_loadbal_subnpart + 1
                        list_del_parts(ppm_loadbal_subnpart) = ppm_loadbal_subnpart
                        DO i=0,ppm_dim-1
                            ppm_loadbal_xpd(ppm_loadbal_subnpart+i) = &
     &                                     xp(i+1,ppm_loadbal_subnpart)
                        ENDDO
                    ENDIF
                ENDIF
            ENDDO
        ENDIF
      ELSE IF (ppm_dim.EQ.3) THEN
        IF (prec .EQ. ppm_kind_single) THEN
            !-----------------------------------------------------------------------
            !  Loop over the particles in 3D
            !-----------------------------------------------------------------------
            DO k=1,Pc%Npart
                !-------------------------------------------------------------------
                !  and check if they are inside the sub
                !-------------------------------------------------------------------
                IF (xp(1,k).GE.min_subs3(1).AND. &
     &              xp(1,k).LE.max_subs3(1).AND. &
     &              xp(2,k).GE.min_subs3(2).AND. &
     &              xp(2,k).LE.max_subs3(2).AND. &
     &              xp(1,k).GE.min_subs3(3).AND. &
     &              xp(1,k).LE.max_subs3(3)) THEN
                    !---------------------------------------------------------------
                    !  In the non-periodic case, allow particles that are
                    !  exactly ON an upper EXTERNAL boundary.
                    !---------------------------------------------------------------
                    IF((xp(1,k).LT.max_subs3(1)                  .OR.  &
     &                 (topo%subs_bc(2,topo%isublist(isub)).EQ.1      .AND. &
     &                  topo%bcdef(2).NE. ppm_param_bcdef_periodic))  .AND. &
     &                 (xp(2,k).LT.max_subs3(2)                  .OR.  &
     &                 (topo%subs_bc(4,topo%isublist(isub)).EQ.1      .AND. &
     &                  topo%bcdef(4).NE. ppm_param_bcdef_periodic))  .OR.  &
     &                  (xp(3,k).LT.max_subs3(3)                 .OR.  &
     &                 (topo%subs_bc(6,topo%isublist(isub)).EQ.1      .AND. &
     &                  topo%bcdef(6).NE. ppm_param_bcdef_periodic))) THEN
                        !-----------------------------------------------------------
                        !  Increase the number of particles to-be-sent and add the
                        !  particle to the buffer
                        !-----------------------------------------------------------
                        ppm_loadbal_subnpart = ppm_loadbal_subnpart + 1
                        list_del_parts(ppm_loadbal_subnpart) = ppm_loadbal_subnpart
                        DO i=0,ppm_dim-1
                            ppm_loadbal_xps(ppm_loadbal_subnpart+i) = &
     &                                     xp(i+1,ppm_loadbal_subnpart)
                        ENDDO
                    ENDIF
                ENDIF
            ENDDO
        ELSE IF (prec .EQ. ppm_kind_double) THEN
            !-----------------------------------------------------------------------
            !  Loop over the particles in 3D
            !-----------------------------------------------------------------------
            DO k=1,Pc%Npart
                !-------------------------------------------------------------------
                !  and check if they are inside the sub
                !-------------------------------------------------------------------
                IF (xp(1,k).GE.min_subd3(1).AND. &
     &              xp(1,k).LE.max_subd3(1).AND. &
     &              xp(2,k).GE.min_subd3(2).AND. &
     &              xp(2,k).LE.max_subd3(2).AND. &
     &              xp(1,k).GE.min_subd3(3).AND. &
     &              xp(1,k).LE.max_subd3(3)) THEN
                    !---------------------------------------------------------------
                    !  In the non-periodic case, allow particles that are
                    !  exactly ON an upper EXTERNAL boundary.
                    !---------------------------------------------------------------
                    IF((xp(1,k).LT.max_subd3(1)                  .OR.  &
     &                 (topo%subs_bc(2,topo%isublist(isub)).EQ.1      .AND. &
     &                  topo%bcdef(2).NE. ppm_param_bcdef_periodic))  .AND. &
     &                 (xp(2,k).LT.max_subd3(2)                  .OR.  &
     &                 (topo%subs_bc(4,topo%isublist(isub)).EQ.1      .AND. &
     &                  topo%bcdef(4).NE. ppm_param_bcdef_periodic))  .OR.  &
     &                  (xp(3,k).LT.max_subd3(3)                 .OR.  &
     &                 (topo%subs_bc(6,topo%isublist(isub)).EQ.1      .AND. &
     &                  topo%bcdef(6).NE. ppm_param_bcdef_periodic))) THEN
                        !-----------------------------------------------------------
                        !  Increase the number of particles to-be-sent and add the
                        !  particle to the buffer
                        !-----------------------------------------------------------
                        ppm_loadbal_subnpart = ppm_loadbal_subnpart + 1
                        list_del_parts(ppm_loadbal_subnpart) = ppm_loadbal_subnpart
                        DO i=0,ppm_dim-1
                            ppm_loadbal_xpd(ppm_loadbal_subnpart+i) = &
     &                                     xp(i+1,ppm_loadbal_subnpart)
                        ENDDO
                    ENDIF
                ENDIF
            ENDDO
        ENDIF
      ENDIF
      !-----------------------------------------------------------------------------
      !  Shorten the list of to-be-deleted particles
      !-----------------------------------------------------------------------------
      iopt = ppm_param_alloc_grow_preserve
      ldu1(1) = ppm_loadbal_subnpart
      CALL ppm_alloc(list_del_parts,ldu1,iopt,info)
        or_fail_alloc("list_del_parts re-alloc failed!")
      !-----------------------------------------------------------------------------
      !  Delete the particles
      !-----------------------------------------------------------------------------
      CALL Pc%del_parts(list_del_parts,ppm_loadbal_subnpart,info)
        or_fail("Error while at deleting particles")
      !-----------------------------------------------------------------------------
      !  Before we send the buffer, we should also trim it so that the size will be
      !  shorter (we are allowed to send only an 1D array)
      !-----------------------------------------------------------------------------
      ldu1(1) = ppm_dim*ppm_loadbal_subnpart
      IF (prec .EQ. ppm_kind_single) THEN
        CALL ppm_alloc(ppm_loadbal_xps,ldu1,iopt,info)
            or_fail_alloc("ppm_loadbal_xps re-alloc failed!")
      ELSE IF (prec .EQ. ppm_kind_double) THEN
        CALL ppm_alloc(ppm_loadbal_xpd,ldu1,iopt,info)
            or_fail_alloc("ppm_loadbal_xpd re-alloc failed!")
      ENDIF

      !-----------------------------------------------------------------------------
      !  Send the buffer, first the size then the data
      !-----------------------------------------------------------------------------
      tag1 = 900
      CALL MPI_Send(ldu1(1),1,MPI_INTEGER,ppm_loadbal_recvrank,tag1, &
     &              ppm_comm,status,info)
        or_fail("# of particle coordinates send failed!")
      tag1 = 1000
      IF (prec .EQ. ppm_kind_single) THEN
        CALL MPI_Send(ppm_loadbal_xps,ldu1(1),MPTYPE,ppm_loadbal_recvrank,&
     &              tag1,ppm_comm,status,info)
            or_fail("ppm_loadbal_xps send failed!")
      ELSE IF (prec .EQ. ppm_kind_double) THEN
        CALL MPI_Send(ppm_loadbal_xpd,ldu1(1),MPTYPE,ppm_loadbal_recvrank,&
     &              tag1,ppm_comm,status,info)
            or_fail("ppm_loadbal_xpd send failed!")
      ENDIF
      !-----------------------------------------------------------------------------
      !  Copy the sent subdomain to the last in the array so that ppm_alloc
      !  can clear it.
      !-----------------------------------------------------------------------------
      stdout("isub:",isub," old_nsub:",old_nsub)
      IF (ppm_dim .EQ. 2) THEN
          IF (prec .EQ. ppm_kind_single) THEN
              topo%min_subs(:,isub) = topo%min_subs(:,old_nsub)
              topo%min_subs(:,old_nsub) = min_subs2

              topo%max_subs(:,isub) = topo%max_subs(:,old_nsub)
              topo%max_subs(:,old_nsub) = max_subs2

              topo%sub_costs( isub) = topo%sub_costs( old_nsub)
              topo%sub_costs(old_nsub) = subcosts
          ELSE IF(prec .EQ. ppm_kind_double) THEN

              topo%min_subd(:,isub) = topo%min_subd(:,old_nsub)
              topo%min_subd(:,old_nsub) = min_subd2

              topo%max_subd(:,isub) = topo%max_subd(:,old_nsub)
              topo%max_subd(:,old_nsub) = max_subd2

              topo%sub_costd( isub) = topo%sub_costd( old_nsub)
              topo%sub_costd(old_nsub) = subcostd
          ENDIF
      ELSE IF (ppm_dim .EQ. 3) THEN
          IF (prec .EQ. ppm_kind_single) THEN
              topo%min_subs(:,isub) = topo%min_subs(:,old_nsub)
              topo%min_subs(:,old_nsub) = min_subs3

              topo%max_subs(:,isub) = topo%max_subs(:,old_nsub)
              topo%max_subs(:,old_nsub) = max_subs3

              topo%sub_costs( isub) = topo%sub_costs( old_nsub)
              topo%sub_costs(old_nsub) = subcosts
          ELSE IF(prec .EQ. ppm_kind_double) THEN

              topo%min_subd(:,isub) = topo%min_subd(:,old_nsub)
              topo%min_subd(:,old_nsub) = min_subd3

              topo%max_subd(:,isub) = topo%max_subd(:,old_nsub)
              topo%max_subd(:,old_nsub) = max_subd3

              topo%sub_costd( isub) = topo%sub_costd( old_nsub)
              topo%sub_costd(old_nsub) = subcostd
          ENDIF

      ENDIF
      !-------------------------------------------------------------------------
      !  New sublist has one LESS sub in it
      !-------------------------------------------------------------------------
      new_nsub =  old_nsub-1
      !-------------------------------------------------------------------------
      !  Continue updating local subdomain info by swaping
      !  isublist, nneighsubs, ineighsubs, nsublist
      !-------------------------------------------------------------------------
      itemp = topo%isublist(isub)
      topo%isublist(    isub) = topo%isublist(old_nsub)
      topo%isublist(old_nsub) = itemp

      itemp = topo%nneighsubs(isub)
      topo%nneighsubs(    isub) = topo%nneighsubs(old_nsub)
      topo%nneighsubs(old_nsub) = itemp


      temp_neigh = topo%ineighsubs(:,isub)
      topo%ineighsubs(:,    isub) = topo%ineighsubs(:,old_nsub)
      topo%ineighsubs(:,old_nsub) = temp_neigh

      topo%nsublist = new_nsub

      !-------------------------------------------------------------------------
      !  Remove the sent subdomain by shrinking the array by one element
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_grow_preserve
      ldu2(1) = ppm_dim
      ldu2(2) = new_nsub
      ldu1(1) = new_nsub

      IF(prec .EQ. ppm_kind_single) THEN
          CALL ppm_alloc(topo%min_subs,ldu2,iopt,info)
            or_fail_alloc("min_subs re-alloc failed!")
          CALL ppm_alloc(topo%max_subs,ldu2,iopt,info)
            or_fail_alloc("max_subs re-alloc failed!")
          CALL ppm_alloc(topo%sub_costs,ldu1,iopt,info)
            or_fail_alloc("sub_costs re-alloc failed!")

      ELSE IF(prec .EQ. ppm_kind_double) THEN
          CALL ppm_alloc(topo%min_subd,ldu2,iopt,info)
            or_fail_alloc("min_subd re-alloc failed!")
          CALL ppm_alloc(topo%max_subd,ldu2,iopt,info)
            or_fail_alloc("max_subd re-alloc failed!")
          CALL ppm_alloc(topo%sub_costd,ldu1,iopt,info)
            or_fail_alloc("sub_costd re-alloc failed!")
      ENDIF
      CALL ppm_alloc(topo%nneighsubs,ldu1,iopt,info)
        or_fail_alloc("nneighsubs re-alloc failed!")
      CALL ppm_alloc(topo%isublist,ldu1,iopt,info)
        or_fail_alloc("isublist re-alloc failed!")
      CALL ppm_alloc(topo%ineighsubs,ldu2,iopt,info)
        or_fail_alloc("ineighsubs re-alloc failed!")

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
      END SUBROUTINE loadbal_sendsub

